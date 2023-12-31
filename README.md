# Web appendix

## Introduction
This repository contains code and output to reproduce results in the paper [A multistate approach to disability insurance
reserving with information delays](https://arxiv.org/abs/2312.14324) by Oliver Lunding Sandqvist. 

## Code structure

* [reserving.R](<reserving.R>): Loads the LEC-DK19 data and the estimated parameters, sample insured for the different categories, calculate reserves using the proposed and naive procedure, save the results and a barchart of the reserves decomposed by categories. This code is used to generate results in Section 5 of the paper. The output is the following.
    * [Results/dataResult.Rda](<Results/dataResult.Rda>): Calculated reserves for the sampled insured.
    * [Figures/Barchart.png](<Figures/Barchart.png>): Barchart of the calculated reserves.
