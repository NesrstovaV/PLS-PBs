# PLS-based-Principal-Balances
The method to construct principal balances based on PLS is introduced in the paper "Nesrstová, V., Wilms, I., Palarea-Albaladejo, J., Filzmoser, P., Martín-Fernández, J.A., Friedecký, D., & Hron, K. (2023). Principal Balances of Compositional Data for Regression
and Classification using Partial Least Squares.", which is currently (July 2023) under revision.

## Requirements
Calculations were performed in following software, using the listed libraries:
- R version: 4.1.3 (2022-03-10) -- "One Push-Up"
- Essential libraries:
    - pls: for PLS; version 2.8-0
    - compositions: to handle compositional data (CoDa); version 2.0-4
    - selbal: version 0.1.0

## Available scripts
1) Scripts with functions:
   
a) Functions for PLS PBs:

fBalChip_PLS.R            -> Main function to obtain PLS PBs -> calls the following 3 functions within its code

fBPMaxOrthNewChip_PLS.R

fBalChipman_PLS.R  

fBPUpChi_PLS.r

Functions to calculate the original PCA PBs introduced in "J. A. Martín-Fernández, V. Pawlowsky-Glahn, J. J. Egozcue, and R. Tolosona-Delgado. Advances in
principal balances for compositional data. Mathematical Geosciences, 50(3):273–298, 2018." are also available at http://www.compositionaldata.com

fBalChip.R            -> Main function to obtain PCA PBs -> calls the following 3 functions within its code

fBPMaxOrthNewChip.R

fBalChipman.R  

fBPUpChi.r


b) Auxiliary functions:
Auxiliary_functions.R

2) Simulation scripts
   
Each simulation scenario is in a separate script. Each script covers CV and the comparison of first balances in for PCA PBs, PLS PBs and Selbal.
- Sim_1.R: one block of markers
- Sim_2.R: 4 same-sized blocks of markers
- Sim_3.R: different-sized blocks of markers

3) An easy-to-use toy example
   
Model_example.R

An easy-to-use script, which shows how to use the function fBalChip_PLS to obtain PLS PBs. CoDa data set X is generated via ilr coordinates.

   

