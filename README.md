# StabilityMetrics

Calculates resistance, recovery time, recovery rate and variability at multiple spatial scales using MODIS version 6 EVI data. This data has been preprocessed by Jon Yearsley. These measures are calculated in Resilience_Xkm_method2standarddev.R where X is the pixel resolution.

Longest recovery time and slowest recovery rate are calculated in Longest_slowest_1km.R.

These measures are calculated using EVI anomalies, which are generated in Multiscale_anomalies.R


Graphics to investigate these metrics are in StabilityXMetrics.R where X is the pixel resolution of investigation and EventsInTime.R.
