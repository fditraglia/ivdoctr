#' Acemoglu, Johnson & Robinson (2001) Dataset
#'
#' @description Cross-country dataset used to construct Table 4 of Acemoglu, Johnson & Robinson (2001).
#' @format A data frame with 64 rows and 9 variables:
#' \describe{
#'    \item{shortnam}{three letter country abbreviation, e.g. AUS for Australia}
#'    \item{africa}{dummy variable =1 if country is in Africa}
#'    \item{lat_abst}{absolute distance to equator (scaled between 0 and 1)}
#'    \item{rich4}{dummy variable, =1 for "Neo-Europes" (AUS, CAN, NZL, USA)}
#'    \item{avexpr}{Average protection against expropriation risk.
#'    Measures risk of government appropriation of foreign private investment
#'    on a scale from 0 (least risk) to 10 (most risk). Averaged over all years
#'    from 1985-1995.}
#'    \item{logpgp95}{Natural logarithm of per capita GDP in 1995 at purchasing
#'    power parity}
#'    \item{logem4}{Natural logarithm of European settler mortality}
#'    \item{asia}{dummy variable, =1 if country is in Asia}
#'    \item{loghjypl}{Natural logarithm of output per worker in 1988}
#'  }
#' @source \url{http://economics.mit.edu/faculty/acemoglu/data/ajr2001}
#' @references \url{https://www.aeaweb.org/articles.php?doi=10.1257/aer.91.5.1369}
"colonial"

#' Becker & Woessmann (2009) Dataset
#'
#' @description This is some data!
#' @source \url{https://www.cesifo-group.de/ifoHome/facts/iPEHD-Ifo-Prussian-Economic-History-Database/publications.html}
#' @references \url{https://www.cesifo-group.de/ifoHome/facts/iPEHD-Ifo-Prussian-Economic-History-Database.html}
#' \url{http://qje.oxfordjournals.org/content/124/2/531.short}
"weber"

#' Blackburn and Neumark (1992) Wage dataset
#'
#' @description A subset of Blackburn and Neumark's (1992) dataset provided by
#' Jeffrey Wooldridge as an example in his introductory econometrics textbook.
#' The subset contains only observations from 1980.
#' @format A data frame with 935 rows and 17 variables:
#' \describe{
#'    \item{wage}{monthly earnings}
#'    \item{hours}{average weekly hours}
#'    \item{IQ}{IQ score}
#'    \item{KWW}{knowledge of world work score}
#'    \item{educ}{years of education}
#'    \item{exper}{years of work experience}
#'    \item{tenure}{years with current employer}
#'    \item{age}{age in years}
#'    \item{married}{=1 if married}
#'    \item{black}{=1 if black}
#'    \item{south}{=1 if live in south}
#'    \item{urban}{=1 if live in SMSA}
#'    \item{sibs}{number of siblings}
#'    \item{brthord}{birth order}
#'    \item{meduc}{mother's education}
#'    \item{feduc}{father's education}
#'    \item{lwage}{natural log of wage}
#' }
#' @source \url{http://www.cengage.com/aise/economics/wooldridge_3e_datasets/}
#' @references \url{http://qje.oxfordjournals.org/content/107/4/1421.full.pdf+html}
"wage2"

