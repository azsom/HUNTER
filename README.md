HUNTER
======

HUNTER - Habitable zone UNcerTainty EstimatoR

This code performs an MC sampling to estimate what fraction of rocky planets could harbor liquid water on their surfaces. HUNTER is described in detail in Zsom, 2014 (under review). The code uses precalculated all-convective T-P profiles from a 1D climate code described in Zsom et al., 2013. Currently HUNTER comes with tables for H2, N2, and CO2 dominated atmospheres.

The code was developed with python 2.7, numpy 1.8.1, scipy 0.12.0 ,matplotlib 1.3.1, and pickle. HUNTER.py is the main file of the code. Currently it is set up to run one simulation that generates 10^5 scenarios. The outcome is saved in /results in human-readable form. plot_data.py contains some routines that will help you to visualize the data. Modify this file to your needs.

HUNTER.py - main file to run the simulation
MC_functions.py - contains all subroutines called by HUNTER.py. The subroutine 'occurrence_rate' performs the gaussian kernel density estimation on the data set of Dressing&Charbonneau 2013. If you want to use the code on a different data set, include a similar subroutine. This file contains all PDFs as well. If you want to experiment with new PDF forms, add it here. The heart of the code is MC_sample. The code flow is illustrated in Fig. 3 of Zsom, 2014 (under review).
data.py - contains some constants and material properties (not everything here is used by HUNTER)
numerics.py - contains numerical methods that are not readily available in numpy or scipy to my knowledge
