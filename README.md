# BinaryStars

## Searching for Visual Binaries Using Gaia Data

This repository contains the [code](binary_star_lookup_in_Gaia_data.v2.py) for an experiment aimed at identifying binary star systems using data from the [Gaia orbital observatory](https://en.wikipedia.org/wiki/Gaia_(spacecraft)), specifically the [Data Release 3 (DR3)](https://www.cosmos.esa.int/web/gaia/dr3#). The Gaia DR3 dataset provides detailed information about stars, including their celestial coordinates (right ascension and declination) and trigonometric parallax, which can be used to calculate the actual distance to each star.

By converting these celestial coordinates and parallax measurements into **three-dimensional Euclidean coordinates**, we can compute the physical distance between pairs of stars. If the calculated distance between two stars falls below a reasonable threshold (approximately **10,000 astronomical units (AU)** — though this value is subject to further discussion and refinement), the stars are likely gravitationally bound and may form part of a [binary or multiple star system](https://en.wikipedia.org/wiki/Visual_binary).

## Results

The output of this project is a [table](binary_systems.md) listing the identified **visual binary star systems**. For each pair of stars, the table includes the following details:

1. **Visible Separation**: The angular separation between the stars, measured in **milli-arcseconds (mas)**.
2. **Physical Separation**: The calculated physical distance between the stars, expressed in **astronomical units (AU)**.
3. **Celestial Coordinates and Parallax**:
   - Right ascension (RA) and declination (Dec) for each star.
   - Trigonometric parallax (`plx`) for each star, used to determine their distances.
4. **Star Identifiers**: The names of the stars, as resolved using the [Simbad catalog](https://aas.aanda.org/articles/aas/pdf/2000/07/ds1821.pdf).
5. **Sky Image**: A link to an image of the sky region centered approximately on the binary system, providing a visual representation of the pair. [NASA SkyView](https://skyview.gsfc.nasa.gov) is used to create these links.

## Example images

**Eta Cassiopeiae** denoted in the tables as `* eta Cas`  
![Eta Cassiopea](https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl?Survey=DSS&Position=12.2832,57.8142&Size=1.0411&Pixels=256&Return=JPEG)

HD 113A and HD 113B  
![HD 113A, HD 113B  ](https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl?Survey=DSS&Position=1.4808,18.0749&Size=0.2661&Pixels=256&Return=JPEG)


## Key Features of the Code

#### 1. Efficient Data Retrieval from Gaia DR3
The **Gaia DR3 database** contains information for nearly **1.5 billion celestial objects**. Due to timeout restrictions on SQL queries and the physical limitations of local hardware, it is impractical to retrieve and analyze the entire dataset locally. To address this, the following constraints were applied:
- **Parallax Threshold**: Only stars with a parallax of at least **1 milli-arcsecond (mas)** were included, ensuring a focus on relatively nearby stars.
- **Brightness Threshold**: Stars dimmer than a **photometric magnitude (phot_g_mean_mag)** of 12 were excluded, reducing the sample to approximately **100,000 stars**.

To make the data retrieval process scalable, the celestial sphere was divided into 400 patches by splitting the right ascension and declination coordinates into 20 subranges each. A separate SQL query is executed for each patch, and the results are aggregated into a single local list for analysis.

#### 2. Data Persistence and Debugging
- **Loading Time**: Retrieving the filtered star data, with a 1-second timeout between queries to avoid overloading the server, takes approximately 30 minutes.
- **Local Storage**: To avoid redundant downloads and facilitate debugging, the retrieved data is stored in a local SQLite3 database. Subsequent analyses extract data directly from this database, ensuring efficiency and reproducibility. Results of Simbad queries are also persisted in this local DB.  

#### 3. Nearest Neighbor Search Optimization
- **Brute-Force Complexity**: A naive nearest-neighbor search has a complexity of **O(n<sup>2</sup>)**, which is computationally infeasible for 100,000 objects.
- **Acceleration with KD-Tree**: To achieve efficient nearest-neighbor searches, a [KD-Tree](https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KDTree.html) data structure is employed. This reduces the search complexity significantly, enabling rapid identification of potential binary star systems.



   
