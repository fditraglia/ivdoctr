#' @import data.table
NULL
#' Acemoglu, Johnson, and Robinson (2001) Dataset
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
  print(paste0("Please download data from http://economics.mit.edu/files/5136 ",
               "and file an issue at https://github.com/fditraglia/ivdoctr/issues",
               "so that we can fix this problem"))
  # # Our example is based on Table 4 of the paper
  # download.file("http://economics.mit.edu/files/5136", "./data/colonial.zip")
  # unzip("./data/colonial.zip", files = "maketable4.dta", exdir = "./data")
  # colonial <- setDT(haven::read_dta("./data/maketable4.dta"))
  # # Use the base sample of countries
  # colonial <- colonial[baseco == 1][, "baseco" := NULL]
  # # Saving and cleaning up workspace
  # save(colonial, file = "./data/colonial.rda", compress = TRUE)
  # rm(colonial)
  # system("rm ./data/colonial.zip ./data/maketable4.dta")
}

#' Becker and Woessmann (2009) Dataset
#'
#' @description Data on Prussian counties in 1871 from Becker and Woessmann's (2009) paper "Was Weber Wrong? A Human Capital Theory of Protestant Economic History."
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
  print(paste0("Please download data from ",
               "https://www.cesifo-group.de/dms/ifodoc/iPEHD/Datasets/ipehd_qje2009_data_tables.zip",
               "and file an issue at https://github.com/fditraglia/ivdoctr/issues",
               "so that we can fix this problem"))
  # download.file("https://www.cesifo-group.de/dms/ifodoc/iPEHD/Datasets/ipehd_qje2009_data_tables.zip",
  #               "./data/weber.zip")
  # unzip("./data/weber.zip", files = "ipehd_qje2009_master.dta", exdir = "./data")
  # weber <- setDT(haven::read_dta("./data/ipehd_qje2009_master.dta"))
  # # Converting non-ASCII characters to UTF-8 encoding
  # enc2utf8(weber$county1871)
  # # Saving and cleaning up workspace
  # save(weber, file = "./data/weber.rda", compress = TRUE)
  # rm(weber)
  # system("rm ./data/weber.zip ./data/ipehd_qje2009_master.dta")
}

#' Burde and Linden (2013, AEJ Applied) Dataset
#'
#' @description Replicates IV using controls from Table 2
#' @format A data frame with 687 rows and 17 variables:
#' \describe{
#'    \item{enrolled}{Indicator if child is enrolled in formal school. Outcome.}
#'    \item{testscore}{Normalized test score}
#'    \item{buildschool}{Indicator if village is treated. Instrument.}
#'    \item{headchild}{Indicator if child is child of head of household}
#'    \item{nhh}{Number of household members}
#'    \item{female}{Female indicator}
#'    \item{age}{Child's age}
#'    \item{yrsvill}{Time family has lived in village}
#'    \item{farsi}{Indicator for speaking Farsi}
#'    \item{tajik}{Indicator for speaking Tajik}
#'    \item{farmers}{Indicator for if head of household is a farmer}
#'    \item{land}{Number of jeribs of land owned}
#'    \item{agehead}{Head of household age}
#'    \item{educhead}{Years of education for head of household}
#'    \item{sheep}{Number of sheep and goats owned}
#'    \item{chagcharan}{Indicator if village is in Chagcharan district}
#'    \item{distschool}{Distance to nearest non-community based school}
#' }
#' @source Provided by author.
#' @references \url{http://www.jstor.org/stable/3083335}
"afghan"

if ("afghan.rda" %in% list.files("./data")) {
} else {
  print(paste0("Please download data from ",
               "https://www.aeaweb.org/aej/app/data/2012-0252_data.zip",
               "and file an issue at https://github.com/fditraglia/ivdoctr/issues",
               "so that we can fix this problem"))
  # download.file("https://www.aeaweb.org/aej/app/data/2012-0252_data.zip",
  #               "./data/afghan.zip")
  # unzip("./data/afghan.zip", exdir = "./data", junkpaths = TRUE,
  #       files = "Data_20120252 2015-06-15/afghanistan_anonymized_data.dta")
  # afghan <- setDT(haven::read_dta("./data/afghanistan_anonymized_data.dta"))
  #
  # setnames(afghan, c("hhid07", "headchild", "female", "age", "yrsvill", "agehead",
  #                    "educhead", "nhh", "land", "sheep", "farsi", "tajik",
  #                    "farmers", "test_ind", "headchild07", "female07", "age07",
  #                    "agehead07", "educhead07", "land07", "sheep07", "yrsvill07",
  #                    "farsi07", "tajik07", "farmers07", "nhh07", "test_ind07",
  #                    "obs07", "obs", "buildschool", "c", "chagcharan",
  #                    "enrolled07", "enrolled", "distschool07", "distschool",
  #                    "testscore07", "testscore", "hhid", "childid"))
  #
  # # Remove outliers following the authors' STATA code
  # outlier <- with(afghan, (nhh07 > 20 & obs07 == 1) |
  #                   (land07 > 10 & obs07 == 1) |
  #                   (sheep07 > 50 & obs07 == 1) |
  #                   (nhh > 20 & obs == 1) |
  #                   (land > 10 & obs == 1) |
  #                   (sheep > 50 & obs == 1))
  # afghan <- afghan[!outlier, .(enrolled, testscore, buildschool, headchild, nhh,
  #                              female, age, yrsvill, farsi, tajik, farmers, land,
  #                              agehead, educhead, sheep, distschool, chagcharan)]
  # # Remove missing and subset to girls
  # afghan <- na.omit(afghan)
  # afghan <- subset(afghan, female == 1)
  # save(afghan, file = "./data/afghan.rda", compress = TRUE)
  # system("rm ./data/afghan.zip ./data/*.dta")
  # rm(afghan)
}
