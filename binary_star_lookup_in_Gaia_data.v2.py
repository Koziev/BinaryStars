import time

import requests
import pandas as pd
import lxml
from astropy.io.votable import parse
from astroquery.vizier import Vizier
from astropy.table import Table
from sklearn.neighbors import KDTree
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
import sqlitedict
import tqdm
import numpy as np
import terminaltables


def pc_to_au(distance_pc):
    """
    Convert distance from parsecs to astronomical units.

    Parameters:
        distance_pc (float): Distance in parsecs.

    Returns:
        float: Distance in astronomical units.
    """
    return distance_pc * 206265


# Function to convert RA and Dec from degrees to radians
def deg_to_rad(deg):
    return deg * (np.pi / 180)


# Function to calculate Cartesian coordinates
def spherical_to_cartesian(ra, dec, distance):
    ra_rad = deg_to_rad(ra)
    dec_rad = deg_to_rad(dec)
    x = distance * np.cos(dec_rad) * np.cos(ra_rad)
    y = distance * np.cos(dec_rad) * np.sin(ra_rad)
    z = distance * np.sin(dec_rad)
    return x, y, z


# Function to calculate distance between two stars
def distance_between_stars(ra1, dec1, parallax1, ra2, dec2, parallax2):
    # Convert parallax to distance
    distance1 = 1000 / parallax1
    distance2 = 1000 / parallax2

    # Convert RA and Dec to Cartesian coordinates
    x1, y1, z1 = spherical_to_cartesian(ra1, dec1, distance1)
    x2, y2, z2 = spherical_to_cartesian(ra2, dec2, distance2)

    # Calculate Euclidean distance
    distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    return distance


def calc_visible_separation(star1, star2) -> float:
    return angular_separation(star1['ra'], star1['dec'], star2['ra'], star2['dec'])


def angular_separation(ra1, dec1, ra2, dec2):
    """
    Calculate the angular separation between two celestial objects.

    Parameters:
        ra1 (float): Right ascension of the first object in degrees.
        dec1 (float): Declination of the first object in degrees.
        ra2 (float): Right ascension of the second object in degrees.
        dec2 (float): Declination of the second object in degrees.

    Returns:
        float: Angular separation in degrees.
    """
    # Convert degrees to radians
    ra1_rad = np.radians(ra1)
    dec1_rad = np.radians(dec1)
    ra2_rad = np.radians(ra2)
    dec2_rad = np.radians(dec2)

    # Calculate the angular separation using the haversine formula
    delta_ra = ra2_rad - ra1_rad
    delta_dec = dec2_rad - dec1_rad

    a = np.sin(delta_dec / 2)**2 + np.cos(dec1_rad) * np.cos(dec2_rad) * np.sin(delta_ra / 2)**2
    c = 2 * np.arcsin(np.sqrt(a))

    # Convert the result back to degrees
    angular_sep_deg = np.degrees(c)
    return angular_sep_deg


def proper_motion_divergence(pm1, pm2):
    m = (pm1 + pm2)*0.5 + 1e-8
    div1 = abs(pm1 - m) / m
    div2 = abs(pm2 - m) / m
    return (div1 + div2) * 0.5


def pow2(x):
    return x*x

def query_star_name(ra, dec):
    coords = f"{ra} {dec}"
    result_table = Simbad.query_region(coords, radius="0d0m10s")
    if result_table:
        min_dist = 1e8
        best_star_name = None
        for row in result_table:
            dist = pow2(ra-row['ra']) + pow2(dec-row['dec'])
            if dist < min_dist:
                min_dist = dist
                best_star_name = row['main_id']

        return best_star_name
    else:
        return 'n/a'



# for row in table1:
#     ra = float(row['_RA'])
#     dec = float(row['_DE'])
#
#     # These values represent the smallest angular separations that Gaia can resolve in its measurements of right ascension and declination.
#     # Resolution in RA: ~2.78 × 10⁻⁵ degrees (at the equator).
#     # Resolution in Dec: ~2.78 × 10⁻⁵ degrees.
#     patch_size = 1e-3
#     ra_min = ra - patch_size
#     ra_max = ra + patch_size
#     dec_min = dec - patch_size
#     dec_max = dec + patch_size
#
#     # Define the query to retrieve star data
#     query = f"""
#     SELECT source_id, non_single_star, ra, dec, parallax, parallax_error, phot_g_mean_mag, pmra, pmdec --, phot_bp_mean_mag, phot_rp_mean_mag
#     FROM gaiadr3.gaia_source
#     WHERE ra >= {ra_min} AND ra <= {ra_max}
#           AND dec >= {dec_min} AND dec <= {dec_max}
#     """
#
#     # Execute the query
#     job = Gaia.launch_job(query)
#
#     # Retrieve the results as an Astropy Table
#     patch_stars = job.get_results()
#     if len(patch_stars) > 0:
#         print('-'*10 + ' HIT ' + '-'*10)
#
#
#         print('DEBUG@82')
#
# exit(0)


stars_db_fp = 'gaia_data.sqlite'
patches_db = sqlitedict.SqliteDict(stars_db_fp, tablename='PATCHES')
stars_db = sqlitedict.SqliteDict(stars_db_fp, tablename='GAIA')
star_names_db = sqlitedict.SqliteDict(stars_db_fp, tablename='STAR_NAMES')
star_names_db.clear()

distances = []

ra_subranges = 20
dec_subranges = 20
patch_coun = 0
total_stars_count = 0
for ira in range(ra_subranges):
    ra_subrange = 24.0/ra_subranges
    ra_min = ira * ra_subrange
    ra_max = ra_min + ra_subrange

    for idec in range(dec_subranges):
        dec_subrange = 180.0 / dec_subranges
        dec_min = -90 + idec*dec_subrange
        dec_max = dec_min + dec_subrange

        patch_coun += 1

        patch_key = f'{ra_min} {ra_max} {dec_min} {dec_max}'
        if patch_key not in patches_db:
            print(f'Start processing patch #{patch_coun}/{ra_subranges * dec_subranges}')

            # Define the query to retrieve star data
            query = f"""
            SELECT source_id, non_single_star, ra, dec, parallax, parallax_error, phot_g_mean_mag, pmra, pmdec --, phot_bp_mean_mag, phot_rp_mean_mag
            FROM gaiadr3.gaia_source
            WHERE parallax > 1 and phot_g_mean_mag < 12
                  AND ra >= {ra_min} AND ra < {ra_max}
                  AND dec >= {dec_min} AND dec < {dec_max}
            """

            # Execute the query
            time.sleep(1)
            job = Gaia.launch_job(query)

            # Retrieve the results as an Astropy Table
            patch_stars = job.get_results()
            patches_db[patch_key] = patch_stars
            patches_db.commit()

            for star in patch_stars:
                star_id = star['source_id']
                stars_db[star_id] = star
            stars_db.commit()

            print('Number of stars loaded so far: {}'.format(len(stars_db)))

print('Total number of stars loaded: {}'.format(len(stars_db)))

# -----------------------------------------------------------------------------------------------

# Initialize a KD-tree structure for fast search for neighbouring stars

star_coords = []
star_ids = []
for gaia_id, star in tqdm.tqdm(stars_db.items(), desc='stars coordinates', total=len(stars_db)):
    star_coords.append((float(star['ra']), float(star['dec'])))
    star_ids.append(gaia_id)

star_kdtree = KDTree(star_coords, leaf_size=2)  # leaf_size is a tuning parameter

binary_systems = []

found_pairs = set()

for gaia_id, star1 in tqdm.tqdm(stars_db.items(), desc='Look for binaries', total=len(stars_db)):
    query_point = (float(star1['ra']), float(star1['dec']))

    # Find the nearest neighbors
    distances, indices = star_kdtree.query([query_point], k=10)  # k=1 means find the 1 nearest neighbor

    # We run 1 point per query, so
    distances = distances[0]
    indices = indices[0]

    ra1 = float(star1['ra'])
    dec1 = float(star1['dec'])
    par1 = float(star1['parallax'])
    par_err1 = float(star1['parallax_error'])

    for index2 in indices:
        star2 = stars_db[star_ids[index2]]
        if star1["source_id"] != star2["source_id"] and (star1["source_id"], star2["source_id"]) not in found_pairs:
            ra2 = float(star2['ra'])
            dec2 = float(star2['dec'])
            par2 = float(star2['parallax'])
            par_err2 = float(star2['parallax_error'])

            a_s = angular_separation(ra1, dec1, ra2, dec2)
            # Gaia's effective angular resolution is 100–200 mas.
            # For bright stars, positional precision is sub-milliarcsecond.
            # Two objects are considered physically separated if their angular separation is greater than Gaia's resolution and their physical separation (considering distance) is significant.
            if a_s > 0.00006:
                dist = distance_between_stars(ra1, dec1, par1, ra2, dec2, par2)
                # A reasonable upper limit for the separation of a stable binary system is 10,000 AU.
                if pc_to_au(dist) < 10_000:

                    # # Filter for stars with similar proper motion
                    if proper_motion_divergence(star1['pmra'], star2['pmra']) < 1.0:
                        if proper_motion_divergence(star1['pmdec'], star2['pmdec']) < 1.0:
                            found_pairs.add((star1["source_id"], star2["source_id"]))
                            found_pairs.add((star2["source_id"], star1["source_id"]))

                            binary_systems.append((dist, star1, star2))

                            # -----------------------------------------------------------------
                            table = [['Orbital separ., au', 'Visible separ., mas', 'Star_1 RA', 'Star_1 DEC', 'Star_1 plx', 'Star_1 name', 'Star_2 RA', 'Star_2 DEC', 'Star_2 plx', 'Star_2 name']]
                            table2 = [table[0] + ['Link']]
                            system_count = 0
                            for dist, row1, row2 in sorted(binary_systems, key=lambda z: z[0])[:100]:
                                system_count += 1

                                # milli-arcseconds (mas)=degrees×60×60×1000
                                vis_separ = calc_visible_separation(row1, row2) * 60 * 60 * 1000
                                #print('Visible separation: {} mas'.format(vis_separ))
                                #print('Orbital separation: {} au ({} pc)'.format(round(pc_to_au(dist), 3), dist))

                                star_key = "{} {}".format(float(row1['ra']), float(row1['dec']))
                                if star_key in star_names_db:
                                    star1_name = star_names_db[star_key]
                                else:
                                    star1_name = query_star_name(row1['ra'], row1['dec'])
                                    star_names_db[star_key] = star1_name
                                    star_names_db.commit()

                                star_key = "{} {}".format(float(row2['ra']), float(row2['dec']))
                                if star_key in star_names_db:
                                    star2_name = star_names_db[star_key]
                                else:
                                    star2_name = query_star_name(row2['ra'], row2['dec'])
                                    star_names_db[star_key] = star2_name
                                    star_names_db.commit()

                                #print('star#1: id={} name=«{}»  ra={} dec={} par={}'.format(row1['source_id'], star1_name, round(row1['ra'], 2), round(row1['dec'], 2), round(row1['parallax'], 2)))
                                #print('star#2: id={} name=«{}»  ra={} dec={} par={}'.format(row2['source_id'], star2_name, round(row2['ra'], 2), round(row2['dec'], 2), round(row2['parallax'], 2)))
                                row = (round(pc_to_au(dist), 1), round(vis_separ, 2),
                                       round(row1['ra'], 4), round(row1['dec'], 4), round(row1['parallax'], 3), star1_name,
                                       round(row2['ra'], 4), round(row2['dec'], 4), round(row2['parallax'], 3), star2_name,
                                       )
                                table.append(row)

                                ra = 0.5*(row1['ra'] + row2['ra'])
                                dec = 0.5*(row1['dec'] + row2['dec'])

                                # SCALE: Arcseconds per pixel (e.g., 0.4 for SDSS).
                                image_size = 256

                                region = vis_separ * 0.02
                                scale = region / image_size
                                link = f"https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl?Survey=DSS&Position={round(ra,4)},{round(dec,4)}&Size={round(scale,4)}&Pixels={image_size}&Return=JPEG"

                                # if -5 < dec < 60:
                                #     link = f"https://skyserver.sdss.org/dr16/SkyServerWS/ImgCutout/getjpeg?ra={round(ra,4)}&dec={round(dec,4)}&scale={round(region,4)}&width={image_size}&height={image_size}"
                                # else:
                                #     link = ''
                                table2.append(list(row) + [f'[image]({link})'])

                            #print('-'*80)
                            #print('Most close stars: {} systems'.format(len(binary_systems)))
                            # Show top-10 systems (table header goes first)
                            print(terminaltables.AsciiTable(table[:11]).table)

                            with open('binary_systems.md', 'w') as f:
                                f.write(terminaltables.GithubFlavoredMarkdownTable(table2).table+'\n')



