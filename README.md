# SAR-Profile-Data

This repository applies a series of quality control measures to SAR Tropical Cyclone Winds data. The main script is process_tc_case2, which utilizes functions from the 'Functions' folder.

The key quality control steps include:

Coastal Buffer Removal: Data within 1 km and 3 km buffer zones from the coastline are removed to eliminate coastal effects.

Azimuthal Statistical Filter: A 3-sigma filter is applied azimuthally. For each radial distance from the storm center, data points with values exceeding three standard deviations from the mean at that radius are identified and removed.
