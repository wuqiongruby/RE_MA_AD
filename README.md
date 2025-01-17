# Random-effects meta-analysis of aggregate data on medical events

This repository contains the example code for the manuscript "Random-effects meta-analysis of aggregate data on medical events". 

The repository is divided by the two data set used for analysis:
* Data set with type 2 diabetes as indication: this is the same data set used in Holzhauer B, Meta-analysis of aggregate data on medical events. Statistics in Medicine. 2017;36(5):723-737. Original data set is available on this paper´s DOI: https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.7181
* Data set with oncology indication: the data is not available due to data protection regulations of data source provider.

Each directory contains the code for simulation study and real data application. We used Rstan to perform Bayesian analysis, therefore you will find R.files paired with stan.files for different parameter settings.
