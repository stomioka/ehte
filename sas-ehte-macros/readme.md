# Readme SAS eHTE

1. eHTE-macro.sas: essential. It contains all the necessary macros the example file will call.
2. eHTE_example: essential. It demonstrates how to read in a dataset which is already in the right format
3. psil201_chg.sas7bdat, a sample dataset if we can share with you.

When you run the SAS, you should get the output like below. The first plot is a histogram plot with kernel density curves overlayed. The 2<sup>nd</sup> plot is the cumulative distribution curves for each arm including the placebo arms. Then you have a small table which contains the active treatment arm which is coded trt01pn=2, and the stdev of ITE (4.15006), the eHTE (0.33547) and p-value (0.0615) for the testing of the H0: no treatment heterogeneity. All the three plots should also be created automatically to the directory you specified in eHTE_example.sas.

![overlapping-hist-score](https://github.com/user-attachments/assets/e070cd4f-24de-4ef3-9903-4786dde99663)
![cum-scores](https://github.com/user-attachments/assets/a91b9ac5-617f-4bfa-949f-4b690f42f9a0)


**HTE and p-value**

| **Obs** | **trt01pn** | **sigma** | **eHTE** | **pvalue** |
| --- | --- | --- | --- | --- |
| **1** | 2   | 4.15006 | 0.33547 | 0.0615 |

![scores-distr](https://github.com/user-attachments/assets/b51791d2-e334-4d98-8eeb-054829bdf7f4)
