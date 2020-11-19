devtools::load_all(here::here())

# We load the databases
# https://unstats.un.org/unsd/methodology/m49/overview/
sdg_ <- read.csv("./data-raw/UNSD — Methodology.csv", sep = ",", header = TRUE)
sdg_$cnt <- sdg_$Country.or.Area
sdg_$iso3 <- sdg_$ISO.alpha3.Code
sdg_$least_developped <- ifelse(sdg_$Least.Developed.Countries..LDC. == "x", 1, 0)
sdg_$region <- sdg_$Region.Name
sdg_$subregion <- sdg_$Sub.region.Name
sdg_$intregion <- sdg_$Intermediate.Region.Name

sdg0 <- subset(sdg_, select = c("cnt", "iso3", "region", "subregion", "intregion", "least_developped")); nrow(sdg0)

data(gbd)
sdg1 <- sdg0[sdg0$iso3 %in% gbd$iso, ]; nrow(sdg1)
# we add taiwan
twn <- data.frame(cnt = "Taiwan", iso3 = "TWN", region = "Asia", subregion = "Eastern Asia", intregion = NA, least_developped = 0)
sdg <- rbind(sdg1, twn)
usethis::use_data(sdg, overwrite = TRUE)

#  ---- We save sdg ----
sdg_ <- read.csv("./data-raw/SDGCountry.csv", sep = ",", header = TRUE)
sdg_$iso3 <- sdg_$Country.Code
sdg_$cnt <- sdg_$Short.Name
sdg_$region <- sdg_$Region
sdg <- subset(sdg_, select = c("country" , "iso3", "region"))
twn <- data.frame(cnt = "Taiwan", iso3 = "TWN", region = "Asia", subregion = "Eastern Asia", intregion = NA, least_developped = 0)
sdg <- rbind(sdg, twn)

sdg_v2 <- sdg[as.character(sdg$iso3) %in% as.character(sdg_$iso), ]
usethis::use_data(sdg_v2, overwrite = TRUE)

# WHO Regions
# database sent by Lynnmarie
who_ <- readxl::read_excel("./data-raw/who_regions.xlsx")
who_$country <- who_$Country
who_$iso3 <- who_$iso
who_$region <- who_$WHORegion_abb
who_$region_name <- who_$who
who_ <- subset(who_, select = c("country", "iso3", "region", "region_name"))
who <- na.omit(who_); nrow(who)

# high-income stratificaiton from World Bank
# https://datahelpdesk.worldbank.org/knowledgebase/articles/906519-world-bank-country-and-lending-groups
wb_ <- read.csv("./data-raw/world_bank_high-income.csv")
wb <- subset(wb_, GroupName == 'High income')

who$region_name_hi <- as.factor(ifelse(who$iso3 %in% wb$CountryCode, 
                                       as.character("High Income Region"), 
                                       as.character(who$region_name)))

usethis::use_data(who, overwrite = TRUE)












