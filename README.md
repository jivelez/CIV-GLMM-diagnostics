# CIV-GLMM-diagnostics

### What's this?

This repository contains the R code and .qmd files used to generate the plots and numerical results, via simulation, presented in the article:

Ángel, J. A.; Vélez, J.I. (2026). Finite-Sample Diagnostics for Random-Effects Misspecification in Poisson Generalized Linear Mixed Models. Axioms 2026.

---


### About the files

#### Numerical experiments 

The files:

* `simulations_CIV_M0.R`
* `simulations_CIV_M1.R`
* `simulations_CIV_M2.R`

correspond to the instructions required to simulate the scenarios and estimate the empirical probability of Type I error and Power for the CIV test.

To run them, use the following instructions from the terminal:

```bash
cd ~/path/to/file/
Rscript simulations_CIV_M0.R
Rscript simulations_CIV_M1.R
Rscript simulations_CIV_M2.R
```

The file

* `civ_white_simulation_src.R`

contains the `R` code to compare the CIV test and White's version. The results are presented in **Table 1** of the article.

To run it, use the following from the terminal:


```bash
cd ~/path/to/file/
Rscript civ_white_simulation_src.R
```


#### Figures

The figures reported in the article can be reproduced through the following files:

* `figures_1_2.qmd`
* `figures_3_4_5.qmd`

Simply download them, open them in `RStudio` or `Positron`, and render the document.



### Help

Should you require assistance, please contact the co-author, [Jorge I. Vélez](https://jorgeivanvelez.netlify.app/).

