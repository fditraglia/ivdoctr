#' @import data.table
#' @import utils
NULL

# Colonial Origins Example ===========================================
if ("colonial.rda" %in% list.files("./data")) {
} else {
  # Our example is based on Table 4 of the paper
  download.file("http://economics.mit.edu/files/5136", "./data/colonial.zip")
  unzip("./data/colonial.zip", files = "maketable4.dta", exdir = "./data")
  colonial <- data.table::setDT(foreign::read.dta("./data/maketable4.dta"))
  # Use the base sample of countries
  colonial <- colonial[baseco == 1]
  colonial[, "baseco" := NULL]
  # Saving and cleaning up workspace
  devtools::use_data(colonial, overwrite = TRUE, compress = "xz")
  rm(colonial)
  system("rm ./data/colonial.zip")
  system("rm ./data/maketable4.dta")
}

# Was Weber Wrong Example ===========================================
if ("weber.rda" %in% list.files("./data")) {
} else {
  download.file("https://www.cesifo-group.de/dms/ifodoc/iPEHD/Datasets/ipehd_qje2009_data_tables.zip", "./data/weber.zip", method = "curl") # need curl since https
  unzip("./data/weber.zip", files = "ipehd_qje2009_master.dta",
        exdir = "./data")
  weber <- data.table::setDT(foreign::read.dta("./data/ipehd_qje2009_master.dta"))
  # Converting non-ASCII characters to UTF-8 encoding
  Encoding(weber$county1871) <- "latin1"
  weber$county1871 <- iconv(weber$county1871, "latin1", "UTF-8")
  attr(weber, "var.labels")[9] <- "Popul. growth 1867-1871 (in %)"
  # Saving and cleaning up workspace
  devtools::use_data(weber, overwrite = TRUE, compress = "xz")
  rm(weber)
  system("rm ./data/weber.zip")
  system("rm ./data/ipehd_qje2009_master.dta")
}

# Wage2 Example from Wooldridge =====================================
if ("wage2.rda" %in% list.files("./data")) {
} else {
  download.file("http://www.cengage.com/aise/economics/wooldridge_3e_datasets/textfiles.ZIP",
                "./data/wooldridge.zip")
  unzip("./data/wooldridge.zip", files = c("WAGE2.raw", "WAGE2.DES"),
        exdir = "./data")
  wage2 <- data.table::fread("./data/WAGE2.raw")
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
  wage2[, "age_sq" := age ^ 2]
  # Savingand cleaning up workspace
  devtools::use_data(wage2, overwrite = TRUE, compress = "xz")
  rm(wage2)
  system("rm ./data/wooldridge.zip")
  system("rm ./data/WAGE2.DES")
  system("rm ./data/WAGE2.raw")
}
