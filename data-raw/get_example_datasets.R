# Colonial Origins Example ===========================================

# Our example is based on Table 4 of the paper
download.file("http://economics.mit.edu/files/5136",
              "./data-raw/colonial.zip")
unzip("./data-raw/colonial.zip", files = "maketable4.dta",
      exdir = "./data-raw")
colonial <- foreign::read.dta("./data-raw/maketable4.dta")
# Use the base sample of countries
colonial <- subset(colonial, baseco == 1)
colonial$baseco <- NULL # drop this column since it isn't used in Table 4
devtools::use_data(colonial)
rm(colonial)
system("rm ./data-raw/colonial.zip")



# Was Weber Wrong Example ===========================================
download.file("https://www.cesifo-group.de/dms/ifodoc/iPEHD/Datasets/ipehd_qje2009_data_tables.zip", "./data-raw/weber.zip",
              method = "libcurl") # need libcurl since https
unzip("./data-raw/weber.zip", files = "ipehd_qje2009_master.dta",
      exdir = "./data-raw")
weber <- foreign::read.dta("./data-raw/ipehd_qje2009_master.dta")
devtools::use_data(weber)
rm(weber)
system("rm ./data-raw/weber.zip")

# Wage2 Example from Wooldridge =====================================

download.file("http://www.cengage.com/aise/economics/wooldridge_3e_datasets/textfiles.ZIP", "./data-raw/wooldridge.zip")
unzip("./data-raw/wooldridge.zip", files = c("WAGE2.raw", "WAGE2.DES"), exdir = "./data-raw")

wage2 <- read.table("./data-raw/WAGE2.raw", na.strings = ".")
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
devtools::use_data(wage2)
rm(wage2)
system("rm ./data-raw/wooldridge.zip")
