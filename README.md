# Overview

This code performs radiative transfer simulations using libradtran, a well-established open source tool for atmospheric radiative transfer simulations. The code can be used to compute the brightness temperatures at different microwave frequencies for a set of positions and times.
Requirements

- Python 3
- libradtran
- numpy
- xarray
- tqdm
- concurrent.futures

Inputs

- A list of timestamps
- A list of latitudes
- A list of longitudes
- A list of altitudes
- A list of atmospheric files

Outputs

The code outputs a netCDF file containing the brightness temperatures at different frequencies for the specified locations and times.
Usage

Performance

The code is parallelized using the concurrent.futures library, which allows the libradtran simulations to run in parallel. The number of workers can be specified using the max_workers parameter in the ProcessPoolExecutor.
Future work


