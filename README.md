# ehte


The Python package `ehte` is designed to assess the heterogeneity of treatment effects between study arms. To evaluate the presence and quantify the degree of treatment response variability, it calculates the variance in the differences in cumulative responses between the placebo and active treatment arms across various percentiles. It tests null hypotheses ($H_0$): suggesting homogeneity in treatment effect.  An alternative hypothesis indicating heterogeneity in treatment effect.

eHTE_Estimator estimates the standard deviation of individual treatment effects (ITE) in two ways. 1. Uses actual patient data to determine percentiles. 2. It calculates percentiles at intervals from 3 to 97 by 2 (a total of 48 percentiles).

See [tutorial](https://htmlpreview.github.io/?https://github.com/stomioka/ehte/blob/main/tutorial.html).

**SAS ehte macro and example SAS code** are in [sas-ehte-macros](https://github.com/stomioka/ehte/tree/main/sas-ehte-macros)

**See our paper** 

[Siegel JS](https://scholar.google.com/citations?hl=en&user=DrHXk5wAAAAJ), Zhong J, [Tomioka S](https://scholar.google.com/citations?user=830xBqsAAAAJ), Ogirala A, Faraone SV, Szabo ST, Koblan KS, [Hopkins SC.](https://scholar.google.com/citations?hl=en&user=oBdcDJ0AAAAJ) Estimating heterogeneity of treatment effect in psychiatric clinical trials. medRxiv [Preprint]. 2024 Apr 23:2024.04.23.24306211. doi: 10.1101/2024.04.23.24306211. PMID: 38712180; PMCID: PMC11071592.
