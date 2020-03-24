# Causal-mechanism-of-extreme-river-discharges

This is the repository for the code related to the paper 

Mhalla, L., Chavez-Demoulin, V., Dupuis, D.J. (2020) "Causal mechanism of extreme river discharges in the upper Danube basin network", Journal of the Royal Statistical Society Series C, to appear


Data 
----

Data set in the paper consists of daily water discharge measurements (in m3/s) at 31 stations in the upper Danube basin, from 1960 to 2010. The dataset can be downloaded from the supplementary material of Asadi, P., Davison, A. C. and Engelke, S. (2015) Extremes on river networks. The Annals of Applied Statistics, 9, 2023–2050. 

[https://projecteuclid.org/euclid.aoas/1453994189#supplemental]

 - The declustered data used in the paper can be loaded in R: load("Data/declustered_data.Rdata").
 - The causal scores of all pairs of stations and their bootstrap replicates can be loaded in R: load("Data/score_edges_all_pairs.Rdata").

Functions
----

 - The CausEV method is implemented in the R file Functions/JRSSC_functions. The functions therein allow to compute the causal score of Eq.(15) in the paper.
 - The data analysis leading to the causal graph of Figure 7 (right) is implemented in the R file JRSSC_Danube_causality.R.


Third-party code 
----

We follow

Asadi, P., Davison, A. C. and Engelke, S. (2015) Extremes on river networks. The Annals of Applied Statistics, 9, 2023–2050.

in order to decluster the daily observations and end-up with independent extreme events over the catchment.

We used selected parts of their code from https://projecteuclid.org/euclid.aoas/1453994189#supplemental.

The declustering procedure using a 9-day window is described in the R file Functions/Declustering_Asadi_et_al_2015.R.

