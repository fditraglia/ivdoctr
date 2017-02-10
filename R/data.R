#' @importFrom haven read_dta
#' @importFrom utils read.table
NULL

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

if ("colonial.rda" %in% list.files("./data")) {
} else {
  # Our example is based on Table 4 of the paper
  download.file("http://economics.mit.edu/files/5136", "./data/colonial.zip")
  unzip("./data/colonial.zip", files = "maketable4.dta", exdir = "./data")
  colonial <- haven::read_dta("./data/maketable4.dta")
  # Use the base sample of countries
  colonial <- subset(colonial, baseco == 1)
  colonial$baseco <- NULL
  # Saving and cleaning up workspace
  save(colonial, file = "./data/colonial.rda", compress = TRUE)
  rm(colonial)
  system("rm ./data/colonial.zip")
  system("rm ./data/maketable4.dta")
}

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
#'    \item{f_over3km}{Percent of pupils farther than 3km from school}
#'    \item{f_mine}{Percent of labor force employed in mining}
#'    \item{inctaxpc}{Income tax revenue per capita in 1877}
#'    \item{perc_secB}{Percentage of labor force employed in manufacturing in 1882}
#'    \item{perc_secC}{Percentage of labor force employed in services in 1882}
#'    \item{perc_secBnC}{Percentage of labor force employed in manufacturing and
#'    services in 1882}
#'    \item{lnyteacher}{100 * Natural logarithm of male elementary school
#'    teachers in 1886}
#'    \item{rhs}{Dummy variable, =1 if Imperial of Hanseatic city in 1517}
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

if ("weber.rda" %in% list.files("./data")) {
} else {
  download.file("https://www.cesifo-group.de/dms/ifodoc/iPEHD/Datasets/ipehd_qje2009_data_tables.zip", "./data/weber.zip", method = "curl") # need curl since https
  unzip("./data/weber.zip", files = "ipehd_qje2009_master.dta",
        exdir = "./data")
  weber <- haven::read_dta("./data/ipehd_qje2009_master.dta")
  # Converting non-ASCII characters to UTF-8 encoding
  Encoding(weber$county1871) <- "latin1"
  weber$county1871 <- iconv(weber$county1871, "latin1", "UTF-8")
  attr(weber, "var.labels")[9] <- "Popul. growth 1867-1871 (in %)"
  # Saving and cleaning up workspace
  save(weber, file = "./data/weber.rda", compress = TRUE)
  rm(weber)
  system("rm ./data/weber.zip")
  system("rm ./data/ipehd_qje2009_master.dta")
}

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

if ("wage2.rda" %in% list.files("./data")) {
} else {
  download.file("http://www.cengage.com/aise/economics/wooldridge_3e_datasets/textfiles.ZIP",
                "./data/wooldridge.zip")
  unzip("./data/wooldridge.zip", files = c("WAGE2.raw", "WAGE2.DES"),
        exdir = "./data")
  wage2 <- utils::read.table("./data/WAGE2.raw")
  # The columns names are explained in WAGE2.DES
  names(wage2) <- c("wage",    #monthly earnings
                    "hours",   #average weekly hours
                    "IQ",      #IQ score
                    "KWW",     #knowledge of world work score
                    "educ",    #years of education
                    "exper",   #years of work experience
                    "tenure",  #years with current employer
                    "age",     #age in years
                    "married", #=1 if married
                    "black",   #=1 if black
                    "south",   #=1 if live in south
                    "urban",   #=1 if live in SMSA
                    "sibs",    #number of siblings
                    "brthord", #birth order
                    "meduc",   #mother's education
                    "feduc",   #father's education
                    "lwage")   #natural log of wage
  # Create age-squared for non-linear earnings profile
  wage2$age_sq <- wage2$age ^ 2
  # Savingand cleaning up workspace
  save(wage2, file = "./data/wage2.rda", compress = TRUE)
  rm(wage2)
  system("rm ./data/wooldridge.zip")
  system("rm ./data/WAGE2.DES")
  system("rm ./data/WAGE2.raw")
}
