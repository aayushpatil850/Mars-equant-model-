# Mars-Equant-Model
# Mars Orbit Prediction Using Equant Model

This repository contains the implementation of an Equant-based model for predicting Mars' orbit. The model estimates key orbital parameters, calculates errors, and visualizes the orbit using heliocentric longitudes of Mars.

## Project Overview

In ancient astronomy, an equant was a point used to explain observed variations in planetary motion. This project applies the equant concept to predict Mars' orbit, iterating over various parameters to find the best-fitting orbit that minimizes the error between predicted and observed heliocentric longitudes.

### Data

The data used for this project is stored in `01_data_mars_opposition_updated.csv`. The dataset contains observations of Mars' heliocentric longitudes in degrees, minutes, and seconds, along with dates of opposition events. The heliocentric longitudes are converted into degrees and then radians for computation purposes.Link for the dataset is https://courses.iisc.ac.in/pluginfile.php/83489/mod_assign/introattachment/0/01_data_mars_opposition_updated.csv?forcedownload=1 

## Model

The Mars orbit model estimates five key parameters:

1. **C (Central Longitude)**: The heliocentric longitude of the central point.
2. **R (Orbit Radius)**: The radius of Mars' orbit.
3. **E1, E2 (Equant parameters)**: Parameters that describe the position of the equant.
4. **Z (Initial Angle Offset)**: The angular displacement from the equant to Mars.
5. **S (Angular Speed)**: The angular speed of Mars' movement.

### Key Functions

- **getIntersectionPoint**: Computes the intersection points based on the equant model.
- **MarsEquantModel**: Implements the Mars orbit prediction model using the equant.
- **bestOrbitInnerParams**: Finds the optimal parameters (`C`, `E1`, `E2`, `Z`) for a given radius `R` and angular speed `S`.
- **bestS**: Iterates over possible values of `S` to find the best angular speed.
- **bestR**: Iterates over possible values of `R` to find the optimal radius of Mars' orbit.
- **bestMarsOrbitParams**: Combines all parameter tuning to find the optimal values for `R`, `S`, `C`, `E1`, `E2`, and `Z` by minimizing the maximum error between predicted and observed heliocentric longitudes.

### Visualization

- **plotMarsOrbit**: Visualizes the Mars orbit using the calculated parameters, highlighting key elements like the sun, center, and equant.

