
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ivdoctr

<!-- badges: start -->

<!-- badges: end -->

Instrumental variables are a useful tool in causal inference. In order
to be valid, researchers impose both formal beliefs (instrumental
relevance and the exclusion restriction) and informal beliefs (e.g.,
correlation between endogenous treatment and the error term). The goal
of `ivdoctr` is to formalize those beliefs and quantifies how sensitive
a researcher’s instrumental variables are to measurement error and
instrument endogeneity. Using data and researcher’s beliefs on
measurement error and instrument endogeneity, this package generates the
space of consistent beliefs across measurement error, instrument
endogeneity, and instrumental relevance for IV regressions.

## Installation

You can install the released version of ivdoctr from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("ivdoctr")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fditraglia/ivdoctr")
```

## Example

This section illustrates how to use the `ivdoctr` package in practice.
This example comes from “The colonial origins of comparative
development: An empirical investigation” by Acemoglu, Johnson, and
Robinson (2001). The authors study the effect of institutions on GDP per
capita across 64 countries. Since institutional quality is endogenous,
they use differences in mortality rates of early western settlers across
colonies as an instrumental variable. The regression specification is as
follows:

![](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign*%7D%20%5Clog\(GDP%20%2F%20capita\)%20%26%3D%20%5Calpha_0%20%2B%20%5Cbeta%20Institutions%20%2B%20u%20%20%5C%5C%20Institutions%20%26%3D%20%5Calpha_1%20%2B%20%5Cpi%20%5Clog\(Settler%20Mortality\)%20%2B%20v%20%5Cend%7Balign*%7D)

The authors state that there is likely a positive correlation between
institutional quality and the error term, which could come from reverse
causality (e.g., wealthier societies can afford better institutions) or
omitted variables (e.g., rule of law or British culture are positively
correlated with present-day institutional quality). This positive
correlation is a researcher belief that can be input into `ivdoctr`
using the `r_TstarU_restriction` argument that accepts a 2-column matrix
of bounds. For the exercise, we use 0.9 as the conservative upper bound
on the extent of the endogeneity.

The authors also state that up to 40% of the measure of “institutions”
is noise. Measurement error is
![1-](https://render.githubusercontent.com/render/math?math=1-%5Ckappa),
so
![\[0.6, 1\]](https://render.githubusercontent.com/render/math?math=%5Ckappa%20%5Cin%20%5B0.6%2C%201%5D)
is the translation of this belief. The code below implements these
beliefs and runs the estimation and saves the TeX table to
“colonial.tex”:

``` r
library(ivdoctr)
endog <- matrix(c(0, 0.9), nrow = 1)
meas <- matrix(c(0.6, 1), nrow = 1)

colonial_example <- makeExample(y_name = "logpgp95", T_name = "avexpr", 
                                z_name = "logem4", data = colonial,
                                controls = NULL, robust = FALSE,
                                r_TstarU_restriction = endog,
                                k_restriction = meas,
                                example_name = "Colonial Origins")
                                
makeTable("colonial.tex", binary = FALSE, colonial_example)
```

To explore the surface of estimates consistent with the researcher’s
beliefs, `ivdoctr` also generates an interactive 3D plot of the surface:

``` r
library(ivdoctr)
endog <- matrix(c(0, 0.9), nrow = 1)
meas <- matrix(c(0.6, 1), nrow = 1)

plot_3d_beta(y_name = "logpgp95", T_name = "avexpr", 
             z_name = "logem4", data = colonial, 
             r_TstarU_restriction = endog,
             k_restriction = meas)
```

## Usage

This package exports three main functions:

  - `makeExample()`: Generates the TeX code that can be used to build a
    regression summary table.
  - `makeTable()`: Generates the TeX code for a stand-alone regression
    table and saves it to the specified file. Uses the output of
    `makeExample()` as its input.
  - `plot_3d_beta()`: Generates an interactive 3D plot illustrating the
    relationship between the causal estimates, instrument endogeneity,
    instrument invalidity, and measurement error.

Both `makeExample` and `plot_3d_beta` use the same primary inputs. Users
input the name of the dataset (`data`), the name of the dependent
variable (`y_name`), the name of the treatment variable(s) (`T_name`),
the name(s) of the instrument(s) (`z_name`), and the names of the
control variables (`controls`). Without any additional arguments, the
functions will output the identified set. If users have beliefs over
measurement error and/or instrument endogeneity, they can specify those
using `k_restriction` and `r_TstarU_restriction`, respectively.
