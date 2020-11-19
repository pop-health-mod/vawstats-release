devtools::load_all(here::here())

# We load the databases (women, 1990 to 2020, 15+ by 5yr age group)
# https://population.un.org/wpp/DataQuery/
WStd <- read.csv('./data-raw/WorldDenom.csv', sep = ',', header = TRUE)
WStd$Pop <- as.numeric(gsub(" ", "", WStd$Pop, fixed = TRUE))
WStd <- subset(WStd, select = c("Location", "Year", "Age", "LAge", "UAge", "Pop"))
write.csv(WStd, "./data-raw/WorldDenom.csv")

WStd$AgeM <- WStd$LAge + (WStd$UAge + 1 - WStd$LAge) / 2
WStd20 <- subset(WStd, Year == 2020 & LAge < 70)
WStd20$Pop[WStd20$LAge == 65] <- sum(WStd$Pop[WStd$Year == 2020 & WStd$LAge >= 70])
WStd20$W <- WStd20$Pop / sum(WStd20$Pop)

usethis::use_data(WStd20, overwrite = TRUE)

# Pre-processing WPP 2019 to obtain age_stratum
# download this from WPP 219, the link for
# https://population.un.org/wpp/Download/Standard/CSV/
wpp19_ <- read.csv('./data-raw/WPP2019_PopulationByAgeSex_Medium.csv', sep = ',', header = TRUE)
wpp19 <- wpp19_[wpp19_$Time %in% c(seq(1990, 2020, 5), 2018) & 
                  wpp19_$AgeGrpStart >= 15, ]; nrow(wpp19)
wpp19$iso3 <- countrycode::countrycode(wpp19$Location, origin = "country.name", destination = "iso3c")
wpp19_na <- wpp19[is.na(wpp19$iso3) & wpp19$Time == 2018 & wpp19$AgeGrpStart == 15, ]

#unique(wpp19_cnt$Location)
loc_exclude <- c("United States of America (and dependencies)", "United Kingdom (and dependencies)", 
             "New Zealand (and dependencies)", "Netherlands (and dependencies)", 
             "Denmark (and dependencies)", "China (and dependencies)",
             "Less developed regions, excluding China", 
             "SIDS Atlantic, and Indian Ocean, Mediterranean and South China Sea (AIMS)")
wpp19_cnt <- wpp19[!is.na(wpp19$iso3) & !(wpp19$Location %in% loc_exclude), ]; nrow(wpp19_cnt)
wpp19_cnt$loage <- wpp19_cnt$AgeGrpStart
wpp19_cnt$hiage <- wpp19_cnt$AgeGrpStart + 4
wpp19_cnt$women <- wpp19_cnt$PopFemale
wpp19_cnt$year <- wpp19_cnt$Time

wpp19_cnt <- wpp19_cnt[, names(wpp19_cnt) %in% c("iso3", "loage", "hiage", "women", "year")]
write.csv(wpp19_cnt, "./data-raw/age_stratums_by_iso3.csv")

# For age-standardization by region
Pop <- read.csv('./data-raw/age_stratums_by_iso3.csv', sep = ',', header = TRUE)
Pop <- na.omit(Pop) # NAs are for the year priors to 1995 (no women categoreis beyond 80+)
Pop$Age <- Pop$loage + (Pop$hiage + 1 - Pop$loage) / 2
Pop10_ <- subset(Pop, year == 2010)
Pop10 <- subset(Pop10_, loage < 70)

Pop_now <- subset(Pop, year == 2018)

# We sum all women in the 65+ categories
C_in_Pop <- as.character(unique(Pop$iso3))
for (i in 1:length(C_in_Pop)) {
  Pop10$women[Pop10$loage == 65 & as.character(Pop10$iso3) == C_in_Pop[i]] <-
        sum(Pop10_$women[Pop10_$loage >= 65 & as.character(Pop10_$iso3) == C_in_Pop[i]])
 }

# We load gbd to relate super region with region and countries
gbd_ <- read.csv('./data-raw/IHME_GBD_2015_LOCATION_HIERARCHIES_Y2016M12D01.CSV', header=TRUE)
gbd_lr <- subset(gbd_, level == 2 & subnational_location == 0)
gbd_lc <- subset(gbd_, level == 3 & subnational_location == 0)
gbd_lc$location_code <- droplevels(gbd_lc$location_code)

gbd <- NULL
for (i in 1:length(unique(gbd_lc$location_code))){
  gbd <- rbind(gbd, data.frame(iso = gbd_lc$location_code[i], NameCnt = gbd_lc$location_name[i],
                    gbd = gbd_lr$location_name[gbd_lr$location_id == gbd_lc$parent_id[i]]))
}

library(plyr)
gbd$Region <- revalue(gbd$gbd, c(
  'High-income Asia Pacific' = "Asia Pacific, High Income",
  'Central Asia' = "Asia, Central",
  "East Asia" = "Asia, East",
  'South Asia' = "Asia, South",
  "Southeast Asia" = "Asia, Southeast",
  'Australasia' = "Australasia",
  'Caribbean' = "Caribbean",
  'Central Europe' = "Europe, Central",
  'Eastern Europe' = "Europe, Eastern",
  'Western Europe' = "Europe, Western",
  'Andean Latin America' = "Latin America, Andean",
  'Central Latin America' = "Latin America, Central",
  'Southern Latin America' = "Latin America, Southern",
  'Tropical Latin America' = "Latin America, Tropical",
  'North Africa and Middle East' = "North Africa/Middle East",
  'High-income North America' = "North America, High Income",
  'Oceania' = "Oceania",
  'Central Sub-Saharan Africa' = "Sub-Saharan Africa, Central",
  'Eastern Sub-Saharan Africa' = "Sub-Saharan Africa, East",
  'Southern Sub-Saharan Africa' = "Sub-Saharan Africa, Southern",
  'Western Sub-Saharan Africa' = "Sub-Saharan Africa, West"))

gbd$SuperRegion <- revalue(gbd$Region, c(
  "Asia Pacific, High Income" = "High Income",
  "Asia, Central" = "Central Europe, Eastern Europe & Central Asia",
  "Asia, East" = "South-East Asia, East Asia & Oceania",
  "Asia, South" = "South Asia",
  "Asia, Southeast" = "South-East Asia, East Asia & Oceania",
  "Australasia" = "High Income",
  "Caribbean" = "Latin America & Caribbean",
  "Europe, Central" = "Central Europe, Eastern Europe & Central Asia",
  "Europe, Eastern" = "Central Europe, Eastern Europe & Central Asia",
  "Europe, Western" = "High Income",
  "Latin America, Andean" = "Latin America & Caribbean",
  "Latin America, Central" = "Latin America & Caribbean",
  "Latin America, Southern" = "High Income",
  "Latin America, Tropical" = "Latin America & Caribbean",
  "North Africa/Middle East" = "North Africa & Middle East",
  "North America, High Income" = "High Income",
  "Oceania" = "South-East Asia, East Asia & Oceania",
  "Sub-Saharan Africa, Central" = "Sub-Saharan Africa",
  "Sub-Saharan Africa, East" = "Sub-Saharan Africa",
  "Sub-Saharan Africa, Southern" = "Sub-Saharan Africa",
  "Sub-Saharan Africa, West" = "Sub-Saharan Africa"))

#  ---- We save gbd ----
usethis::use_data(gbd, overwrite = TRUE)

data(gbd)
# Formatting WPP
C_in_Pop <- as.character(unique(Pop$iso3))
C_in_GBD <- as.character(unique(gbd$iso))
# 12 countries in World Pop Prospect not in GBD
# Aruba, Curacao, Western Sahara, Guadeloupe, French Guiana, Hong Kong, Macao, Martinique, Mayotte, New Caledonia, French Polynesia, Reunion
C_in_Pop[!(C_in_Pop %in% C_in_GBD)]
# 5 countries in GBD are not in World Pop Prospect (b/c pop is smaller than 90K)
to_add <- C_in_GBD[!(C_in_GBD %in% C_in_Pop)]; to_add

# We add the 5 countries in GBD that are not in WPP (only for 2015)
mat_to_add <- expand.grid(X = NA, iso3 = to_add,
                          loage = seq(15, 104, by = 5),
                          hiage = NA, women = NA,
                          year = c(2010, 2018))
mat_to_add$hiage <- mat_to_add$loage + 4
mat_to_add$Age <- mat_to_add$loage + 2.5
mat_to_add <- mat_to_add[order(mat_to_add$iso3, mat_to_add$year, mat_to_add$loage), ]

mat_to_add_2010 <- subset(mat_to_add, year == 2010 & loage <= 65)
mat_to_add_now <- subset(mat_to_add, year == 2018)

Pop10 <- rbind(Pop10, mat_to_add_2010)
Pop_now <- rbind(Pop_now, mat_to_add_now)

to_assign <- data.frame(iso = c("ABW", "CUW", "ESH", 
                                 "GLP", "GUF", "HKG",
                                 "MAC", "MTQ", "MYT", 
                                 "NCL", "PYF", "REU"),
                        NameCnt = NA,
                        gbd = NA,
                        Region = c('Caribbean', 'Caribbean', 'North Africa/Middle East',
                                   'Caribbean', 'Caribbean', 'Asia, East',
                                   'Asia, East', 'Caribbean',"Sub-Saharan Africa, East",
                                   "Oceania","Oceania", "Asia, Southeast"),
                        SuperRegion = c('Latin America & Caribbean', 'Latin America & Caribbean', 'North Africa & Middle East',
                                        'Latin America & Caribbean', 'Latin America & Caribbean', 'South-East Asia, East Asia & Oceania',
                                        'South-East Asia, East Asia & Oceania','Latin America & Caribbean',"Sub-Saharan Africa",
                                        'South-East Asia, East Asia & Oceania','South-East Asia, East Asia & Oceania', 'South-East Asia, East Asia & Oceania'))
gbd <- rbind(gbd, to_assign)

# We assgined GBD region and super region to WPP
Pop10$GBD_R <- NA
Pop10$GBD_SR <- NA
C_in_Pop10 <- as.character(unique(Pop10$iso3))
for (i in 1:length(C_in_Pop10)) {
  R.i <- gbd$Region[gbd$iso == C_in_Pop10[i]]
  SR.i <- gbd$SuperRegion[gbd$iso == C_in_Pop10[i]]
  if (length(R.i) > 0) {
    Pop10$GBD_R[Pop10$iso3 == C_in_Pop10[i]] <- as.character(R.i) }
  if (length(SR.i) > 0) {
    Pop10$GBD_SR[Pop10$iso3 == C_in_Pop10[i]] <- as.character(SR.i)  }
}

Pop_now$GBD_R <- NA
Pop_now$GBD_SR <- NA
C_in_Pop_now <- as.character(unique(Pop_now$iso3))
for (i in 1:length(C_in_Pop_now)) {
  R.i <- gbd$Region[gbd$iso == C_in_Pop_now[i]]
  SR.i <- gbd$SuperRegion[gbd$iso == C_in_Pop_now[i]]
  if (length(R.i) > 0) {
    Pop_now$GBD_R[Pop_now$iso3 == C_in_Pop_now[i]] <- as.character(R.i) }
  if (length(SR.i) > 0) {
    Pop_now$GBD_SR[Pop_now$iso3 == C_in_Pop_now[i]] <- as.character(SR.i) }
}


# ---- Population weight by region ----
EverSex_ <- read.csv("./data-raw/ever_had_sex.csv", header = TRUE)
EverSex <- subset(EverSex_, year == 2010)

WPopRegion10 <- NULL
WPopRegion_now <- NULL
UniqueR <- as.character(unique(gbd$Region))
for (i in 1:length(UniqueR)) {
  region10.i <- subset(Pop10, GBD_R == UniqueR[i])
  region10.i$WomenSex <- 0
  region_now.i <- subset(Pop_now, GBD_R == UniqueR[i])
  region_now.i$WomenSex <- 0
  for (j in 1:dim(region10.i)[1]) {
    HadSex10.j <- EverSex$ever_had_sex[EverSex$iso3 == as.character(region10.i$iso3)[j] &
                        EverSex$age_group == ifelse(region10.i$loage[j] > 45, 45, region10.i$loage[j])]
    if (length(HadSex10.j) > 0) {
      region10.i$WomenSex[j] <- region10.i$women[j] * HadSex10.j }
    }
  for (j in 1:dim(region_now.i)[1]) {
    HadSex_now.j <- EverSex$ever_had_sex[EverSex$iso3 == as.character(region_now.i$iso3)[j] &
                        EverSex$age_group == ifelse(region_now.i$loage[j] > 45, 45, region_now.i$loage[j])]
    if (length(HadSex_now.j) > 0) {
      region_now.i$WomenSex[j] <- region_now.i$women[j] * HadSex_now.j }
    }
  #' <<<<<<<<<<<<<<<<<<<< *** >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # If you want to take into account <ever had sex>, change '$women' for '$WomenSex'
  #' <<<<<<<<<<<<<<<<<<<< *** >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  agg10.i <- aggregate(region10.i$women, by = list(region10.i$loage), FUN = sum, na.rm = TRUE)
  agg10.i$region <- UniqueR[i]
  agg10.i$loage <- agg10.i$Group.1 + 2.5
  agg10.i$Denom <- agg10.i$x + 0.5 # we add 500 to each count to avoid zero values.
  WPopRegion10 <- rbind(WPopRegion10, agg10.i)
  agg_now.i <- aggregate(region_now.i$women, by = list(region_now.i$loage), FUN = sum, na.rm = TRUE)
  agg_now.i$region <- UniqueR[i]
  agg_now.i$loage <- agg_now.i$Group.1 + 2.5
  agg_now.i$Denom <- agg_now.i$x + 0.5 # we add 500 to each count to avoid zero values.
  WPopRegion_now <- rbind(WPopRegion_now, agg_now.i)
}

# We reformat in long form
wom_by_age_and_region_2010 <- NULL
for (i in 1:length(UniqueR)){
  region10.i <- subset(WPopRegion10, region == UniqueR[i])
  wom_by_age_and_region_2010 <- cbind(wom_by_age_and_region_2010, region10.i$Denom)
}

#  ---- We save denom by region for 2010 ----
usethis::use_data(wom_by_age_and_region_2010, overwrite = TRUE)


# We inpute country with population size smaller than 90K (we assume 45K)
# and the age distribution of their region. (For 2020)
miss_iso <- unique(Pop_now$iso3[which(is.na(Pop_now$women))])
for (i in 1:length(miss_iso)) {
  region_i <- gbd$Region[gbd$iso == as.character(miss_iso[i])]
  denom_r <- WPopRegion_now[WPopRegion_now$region == region_i, ]
  denom_r$w <- denom_r$Denom / sum(denom_r$Denom)
  Pop_now$women[Pop_now$iso3 == miss_iso[i]] <- 45 * denom_r$w
}
denom_cnt_age_now <- Pop_now

#  ---- We save denom by country for 2020 ----
usethis::use_data(denom_cnt_age_now, overwrite = TRUE)

# World standard population now
wpp_std_now <- aggregate(Pop_now$women, by = list(Pop_now$loage), FUN = sum, na.rm = TRUE)
wpp_std_now$AgeM <- wpp_std_now$Group.1 + 2.5
wpp_std_now$Pop <- wpp_std_now$x
wpp_std_now <- wpp_std_now[, -c(1, 2)]

#  ---- We save population standard in NOW ----
usethis::use_data(wpp_std_now, overwrite = TRUE)

# Ever Sex Weighting
EverSex_ <- read.csv("./data-raw/ever_had_sex.csv", header = TRUE)
denom_ever_sex_2010 <- subset(EverSex_, year == 2010)

#  ---- We save ever had sex for 2010 ----
usethis::use_data(denom_ever_sex_2010, overwrite = TRUE)

