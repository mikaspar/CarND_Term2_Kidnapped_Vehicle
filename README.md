# Kidnapped Vehicle Project

## Project Introduction

A robot has been kidnapped and transported to a new location! Luckily it has a map of this location, a (noisy) GPS estimate of its initial location, and lots of (noisy) sensor and control data.

In this project I implemented a 2 dimensional particle filter in C++. The particle filter ist given a map and some initial localization information (analogous to what a GPS would provide). At each time step filter gets observation and control data.

## Particle Filter Algorithm

1. A Particle structure and Particle Filter class are initialized with random values and sensor standard deviations.
2. Prediction step for all particles is performed based on measurement data.
3. Transformation of observation from vehicle to map coordinate system ist performed and only landmarks in the sensor range are taken into account.
4. The particle weights are calculated based on multivariate (x,y) Gaussian probability function.
5. The particles are resampled proportional to their weights by means of the wheel algorithm.

### Result






