# BinaryStars

### Searching for Visual Binaries Using Gaia Data

This repository contains the code for an experiment aimed at identifying binary star systems using data from the [Gaia orbital observatory](https://en.wikipedia.org/wiki/Gaia_(spacecraft)), specifically the [Data Release 3 (DR3)](https://www.cosmos.esa.int/web/gaia/dr3#). The Gaia DR3 dataset provides detailed information about stars, including their celestial coordinates (right ascension and declination) and trigonometric parallax, which can be used to calculate the actual distance to each star.

By converting these celestial coordinates and parallax measurements into **three-dimensional Euclidean coordinates**, we can compute the physical distance between pairs of stars. If the calculated distance between two stars falls below a reasonable threshold (approximately **10,000 astronomical units (AU)** — though this value is subject to further discussion and refinement), the stars are likely gravitationally bound and may form part of a [binary or multiple star system](https://en.wikipedia.org/wiki/Visual_binary).


