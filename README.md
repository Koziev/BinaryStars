# BinaryStars

### Searching for Visual Binaries Using Gaia Data

This repository contains the [code](binary_star_lookup_in_Gaia_data.v2.py) for an experiment aimed at identifying binary star systems using data from the [Gaia orbital observatory](https://en.wikipedia.org/wiki/Gaia_(spacecraft)), specifically the [Data Release 3 (DR3)](https://www.cosmos.esa.int/web/gaia/dr3#). The Gaia DR3 dataset provides detailed information about stars, including their celestial coordinates (right ascension and declination) and trigonometric parallax, which can be used to calculate the actual distance to each star.

By converting these celestial coordinates and parallax measurements into **three-dimensional Euclidean coordinates**, we can compute the physical distance between pairs of stars. If the calculated distance between two stars falls below a reasonable threshold (approximately **10,000 astronomical units (AU)** — though this value is subject to further discussion and refinement), the stars are likely gravitationally bound and may form part of a [binary or multiple star system](https://en.wikipedia.org/wiki/Visual_binary).

### Results

The output of this project is a [table](binary_systems.md) listing the identified **visual binary star systems**. For each pair of stars, the table includes the following details:

1. **Visible Separation**: The angular separation between the stars, measured in **milli-arcseconds (mas)**.
2. **Physical Separation**: The calculated physical distance between the stars, expressed in **astronomical units (AU)**.
3. **Celestial Coordinates and Parallax**:
   - Right ascension (RA) and declination (Dec) for each star.
   - Trigonometric parallax (`plx`) for each star, used to determine their distances.
4. **Star Identifiers**: The names of the stars, as resolved using the [Simbad catalog](https://aas.aanda.org/articles/aas/pdf/2000/07/ds1821.pdf).
5. **Sky Image**: A link to an image of the sky region centered approximately on the binary system, providing a visual representation of the pair. [NASA SkyView](https://skyview.gsfc.nasa.gov) is used to create these links.

   
