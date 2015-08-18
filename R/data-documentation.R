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
#' @description Data on Prussian counties in 1871 from Becker & Woessmann's (2009) paper "Was Weber Wrong? A Human Capital Theory of Protestant Economic History."
#' @format A data frame with 452 rows and 44 variables:
#' \describe{
#'    \item{kreiskey1871}{kreiskey1871}
#'    \item{county1871}{County name in 1871}
#'    \item{rbkey}{District key}
#'    \item{lat_rad}{Latitude (in rad)}
#'    \item{lon_rad}{Longitude (in rad)}
#'    \item{kmwittenberg}{Distance to Wittenberg (in km)}
#'    \item{zupreussen}{Year in which county was annexed by Prussia}
#'    \item{hhsize}{Average household size}
#'    \item{gpop}{Population growth from 1867-1871 in percentage points}
#'    \item{f_prot}{Percent Protestants}
#'    \item{f_jew}{Percent Jews}
#'    \item{f_rw}{Percent literate}
#'    \item{f_miss}{Percent missing education information}
#'    \item{f_young}{Percent below the age of 10}
#'    \item{f_fem}{Percent female}
#'    \item{f_ortsgeb}{Percent born in municipality}
#'    \item{f_pruss}{Percent of Prussian origin}
#'    \item{f_blind}{Percent blind}
#'    \item{f_deaf}{Percent deaf-mute}
#'    \item{f_dumb}{Percent insane}
#'    \item{f_urban}{Percent of county population in urban areas}
#'    \item{lnpop}{Natural logarithm of total population size}
#'    \item{lnkmb}{Natural logarithm of distance to Berlin (km)}
#'    \item{poland}{Dummy variable, =1 if county is Polish-speaking}
#'    \item{latlon}{Latitude * Longitude * 100}
#'    \item{f_over_3km}{Percent of pupils farther than 3km from school}
#'    \item{f_mine}{Percent of labor force employed in mining}
#'    \item{inctaxpc}{Income tax revenue per capita in 1877}
#'    \item{perc_secB}{Percentage of labor force employed in manufacturing in 1882}
#'    \item{perc_secC}{Percentage of labor force employed in services in 1882}
#'    \item{perc_secCnB}{Percentage of labor force employed in manufacturing and
#'    services in 1882}
#'    \item{lnyteacher}{100 * Natural logarithm of male elementary school
#'    teachers in 1886}
#'    \item{rhs}{Dummy variable, =1 if Imperial of Hanseatic city in 1517}
#'    \item{yteacher}{Income of male elementary school teachers in 1886}
#'    \item{yteacher}{Income of male elementary school teachers in 1886}
#'    \item{pop}{Total population size}
#'    \item{kmb}{Distance to Berlin (km)}
#'    \item{uni1517}{Dummy variable, =1 if University in 1517}
#'    \item{reichsstadt}{Dummy variable, =1 if Imperial city in 1517}
#'    \item{hansestadt}{Dummy variable, =1 if Hanseatic city in 1517}
#'    \item{f_cath}{Percentage of Catholics}
#'    \item{sh_al_in_tot}{Share of municipalities beginning with letter A to L}
#'    \item{ncloisters1517_pkm2}{Monasteries per square kilometer in 1517}
#'    \item{school1517}{Dummy variable, =1 if school in 1517}
#'    \item{dnpop1500}{City population in 1500}
#'  }
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
#'    \item{age_sq}{age in years squared}
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

