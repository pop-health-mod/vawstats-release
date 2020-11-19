library(rdhs)
#  devtools::install_github("mrc-ide/demogsurv")

 rdhs::set_rdhs_config(data_frame = "data.table::as.data.table",
                        email = toString("xxx@xxx"),
                        project = "Global estimates of violence against women",
                        config_path = "rdhs.json",
                        global = FALSE)


library(haven)
library(survey)
library(demogsurv)
library(data.table)
library(sf)
library(tidyverse)

#' Identify HIV behaviour characteristic (domestic violence)
survchar <- dhs_survey_characteristics()
survchar[grepl("Domestic violence", SurveyCharacteristicName)][order(SurveyCharacteristicID)]

#' Before downloading the individual level data, we can query the API
data_api <- dhs_data(indicatorIds = c("DV_FSVL_W_POS",  # Percentage of ever married women whose current or most recent husband or partner committed any form of physical and/or sexual violence
                                      "DV_FSVL_W_PPS",  # Percentage of ever married women for whom any husband or partner committed any physical and/or sexual violence
                                      "DV_SPVL_W_POS",  # Percentage of ever married women who have ever experienced physical or sexual violence committed by their husband or partner
                                      "DV_SPV1_W_POS"), # Percentage of ever married women who have experienced physical or sexual violence committed by their husband/partner in the 12 months preceding the survey
                   breakdown = "national", surveyYearStart = 2000)
head(data_api); dim(data_api)
length(unique(data_api$SurveyId))
# write.csv(data_api, file = "~/Desktop/data_vaw_api.csv")

# ---- downloaing microdata ----
#' Identify surveys with HIV testing from desired countries
surveys <- dhs_surveys(surveyCharacteristicIds = 6, surveyYearStart = 2000)

#' IPV
#' View(dhs_indicators(returnFields=c("IndicatorId","Label","Definition")))
dat <- dhs_data(indicatorIds = c("DV_FSVL_W_POS",  # Percentage of ever married women whose current or most recent husband or partner committed any form of physical and/or sexual violence
                                 "DV_FSVL_W_PPS",  # Percentage of ever married women for whom any husband or partner committed any physical and/or sexual violence
                                 "DV_SPVL_W_POS",  # Percentage of ever married women who have ever experienced physical or sexual violence committed by their husband or partner
                                 "DV_SPV1_W_POS",  # Percentage of ever married women who have experienced physical or sexual violence committed by their husband/partner in the 12 months preceding the survey
                                 "DV_EXSV_W_EVR",  # Women who ever experienced sexual violence
                                 "DV_EXSV_W_12M",  # Women who experienced sexual violence in past 12 months
                                 "DV_PCSV_W_CHP"),  # Sexual violence committed by current husband/partner
                  surveyYearStart = 1990)
setdiff(surveys$SurveyId, dat$SurveyId)

survids <- union(dat$SurveyId, surveys$SurveyId)



#' Identify datasets for these surveys
#' Remove India temporarily (too big)
survids <- survids[survids != "IA2015DHS"]
ird <- dhs_datasets(fileType = "IR", fileFormat = "flat")[SurveyId %in% survids]

#' Get local path to dataset (download if needed)
ird$path <- unlist(get_datasets(ird))


#' Variables to extract from individual / male recode datasets
ir <- readRDS(tail(ird$path, 1))

irvars <- c("v001",  # Cluster number.
            "v002",  # Household number.
            "v012",  # Current age of respondent.
            "v022",  # Sample strata for SE.
            "v044",  # Selected for domestic violence module.
            "v102",  # rural/urban.
            "v005",  # Individual weight.
            "v501",  # Current marital status of respondent.
            "v502",  # Women currently in union.
            "v503",  # Whether the respondent has been married or lived with a man once or more than once.
            "v525",  # Ever had sex.
            "d005",  # Weight for domestic violence.
            "d105a", # Spouse ever pushed, shook or threw something.
            "d105b", # Spouse ever slapped.
            "d105c", # Spouse ever punched with fist or something harmful. (Severe)
            "d105d", # Spouse ever kicked or dragged. (Severe)
            "d105e", # Spouse ever tried to strangle or burn. (Severe)
            "d105f", # Spouse ever threatened with knife/gun or other weapon. (Severe)
            "d105g", # Spouse ever attacked with knife/gun or other weapon. (Severe)
            "d105h", # Spouse ever physically forced sex when not wanted. (Severe)
            "d105i", # Spouse ever forced other sexual acts when not wanted. (Severe)
            "d105j", # Spouse ever twisted her arm or pulled her hair.
            "d106",  # Experienced any less severe violence.
            "d107",  # Experienced any severe violence.
            "d108",  # Experienced any sexual violence.
            "d130a", # Previous husband: ever hit, slap, kick or physcially hurt repsondent.
            "d130b", # Previous husband: physically forced to have sex or to perform sexual acts.
            "d124",  # Anyone (besides partner) forced respondent to have sexual intercourse in past 12 months
            "d125",  # Ever forced to perform unwanted sexual acts
            "d127")  # Person who forced respondent into first sexual act

allvars <- c("SurveyId", "CountryName", "SurveyYear", irvars)


#' Load and merge datasets
datlst <- list()
for(survid in ird$SurveyId){

  print(survid)

  ir <- readRDS(ird[SurveyId == survid]$path)
  # We only keep women that are selected for domestic violence module
  if (!is.null(ir$v044)) {
      ir <- ir[ir$v044 == 1, ]
  }

  ir <- ir[, intersect(irvars, names(ir))]
  ir[setdiff(irvars, names(ir))] <- NA

  dat <- ir
  dat[setdiff(allvars, names(dat))] <- NA
  dat$SurveyId <- survid
  dat$CountryName <- ird[SurveyId == survid]$CountryName
  dat$SurveyYear <- ird[SurveyId == survid]$SurveyYear

  datlst[[survid]] <- dat[allvars]
}

set_na <- function(x, na_codes = 9){ x[x %in% na_codes] <- NA; x }


#' Create VAW indicators
#' - Ever IPV
recode_dhs <- function(dat) {

  name <- dat$SurveyId[1]
  print(name)

  dat$country <- dat$CountryName
  dat$surveyid <- dat$SurveyId
  dat$survyear <- dat$SurveyYear
  dat$psu <- as.numeric(as.character(dat$v001))
  dat$stratum <- 1 # dat$v022
  dat$age <- dat$v012
  dat$selected_domvio <- dat$v044
  dat$restype <- as_factor(dat$v102)
  dat$indweight <- dat$v005 / 1e6
  dat$vawweight <- dat$d005 / 1e6

  if (all(is.na(dat$vawweight))) { return(NULL) }

  #' ## Sexual behaviour outcomes
  dat$cmu <- ifelse(dat$v502 == 1, 1, 0)       # currently married or union
  dat$evermu <- ifelse(dat$v501 != 0, 1, 0)    # ever married or union
  dat$eversex <- ifelse(dat$v525 != 0, 1, 0)   # Need to verify handling of missing / unknown

  #' VAW
  dat$d106c <- ifelse(is.na(set_na(dat$d106, 9)), -999, dat$d106)
  dat$d107c <- ifelse(is.na(set_na(dat$d107, 9)), -999, dat$d107)
  dat$d108c <- ifelse(is.na(set_na(dat$d108, 9)), -999, dat$d108)
  dat$cmr_ipv <- ifelse(dat$d106c > 0 | dat$d107c > 0 | dat$d108c > 0, 1,
                 ifelse(dat$d106c == -999 | dat$d107c == -999 | dat$d108c == -999, NA, 0))
  if(dat$SurveyId[1] == 'JO2012DHS') {
    # sexual violence coded only as 1 or NA. Seems that they don't have 0.
      dat$cmr_ipv <- ifelse(dat$d106c > 0 | dat$d107c > 0 | dat$d108c > 0, 1,
                 ifelse(dat$d106c == -999 | dat$d107c == -999, NA, 0))
  }
  dat$d130ac <- ifelse(is.na(set_na(dat$d130a, 9)), -999, dat$d130a)
  dat$d130bc <- ifelse(is.na(set_na(dat$d130b, 9)), -999, dat$d130b)
  dat$anypartner <- ifelse(dat$d130ac > 0 | dat$d130bc > 0, 1,
                    ifelse(dat$d130ac == -999 | dat$d130bc == -999, NA, 0))
  dat$ever_ipv <- dat$cmr_ipv
  dat$ever_ipv[dat$anypartner == 1] <- 1

  # If sexual violence not recorded, we set to NA
  if (all(is.na(dat$d108))) { dat$ever_ipv <- NA }
  
  # Any severe violence
  dat$ever_severe <- ifelse(dat$d107c == -999, NA,
                            ifelse(dat$d107c == 1, 1, 0))
  
  # Last 12 months IPV (Physical and/or sexual)
  dat$d105ac <- ifelse(dat$d105a %in% c(4, 9), -999, dat$d105a)
  dat$d105bc <- ifelse(dat$d105b %in% c(4, 9), -999, dat$d105b)
  dat$d105cc <- ifelse(dat$d105c %in% c(4, 9), -999, dat$d105c)
  dat$d105dc <- ifelse(dat$d105d %in% c(4, 9), -999, dat$d105d)
  dat$d105ec <- ifelse(dat$d105e %in% c(4, 9), -999, dat$d105e)
  dat$d105fc <- ifelse(dat$d105f %in% c(4, 9), -999, dat$d105f)
  dat$d105gc <- ifelse(dat$d105g %in% c(4, 9), -999, dat$d105g)
  dat$d105hc <- ifelse(dat$d105h %in% c(4, 9), -999, dat$d105h)
  dat$d105ic <- ifelse(dat$d105i %in% c(4, 9), -999, dat$d105i)
  dat$d105jc <- ifelse(dat$d105j %in% c(4, 9), -999, dat$d105j)
  
  dat$p12m_ipv <- ifelse(dat$d105ac %in% c(1, 2) | 
                         dat$d105bc %in% c(1, 2) | 
                         dat$d105cc %in% c(1, 2) | 
                         dat$d105dc %in% c(1, 2) | 
                         dat$d105ec %in% c(1, 2) | 
                         dat$d105fc %in% c(1, 2) | 
                         dat$d105gc %in% c(1, 2) | 
                         dat$d105hc %in% c(1, 2) | 
                         dat$d105ic %in% c(1, 2) | 
                         dat$d105jc %in% c(1, 2), 1,
                  ifelse(dat$d105ac %in% c(-999) & 
                         dat$d105bc %in% c(-999) & 
                         dat$d105cc %in% c(-999) & 
                         dat$d105dc %in% c(-999) & 
                         dat$d105ec %in% c(-999) & 
                         dat$d105fc %in% c(-999) & 
                         dat$d105gc %in% c(-999) &
                         dat$d105hc %in% c(-999) & 
                         dat$d105ic %in% c(-999) &
                         dat$d105jc %in% c(-999), NA, 0))
  
  dat$p12m_severe <- ifelse(dat$d105cc %in% c(1, 2) | 
                         dat$d105dc %in% c(1, 2) | 
                         dat$d105ec %in% c(1, 2) | 
                         dat$d105fc %in% c(1, 2) | 
                         dat$d105gc %in% c(1, 2) |
                         dat$d105hc %in% c(1, 2) | 
                         dat$d105ic %in% c(1, 2), 1,
                  ifelse(dat$d105cc %in% c(-999) & 
                         dat$d105dc %in% c(-999) & 
                         dat$d105ec %in% c(-999) &
                         dat$d105fc %in% c(-999) &
                         dat$d105gc %in% c(-999) &
                         dat$d105hc %in% c(-999) &
                         dat$d105ic %in% c(-999), NA, 0))
  
  # NPSV
  dat$d124_ <- set_na(dat$d124, na_codes = 6)
  dat$d125_ <- set_na(dat$d125, na_codes = 6)
  dat$d127_ <- ifelse(dat$d127 >= 95, NA, dat$d127)
  dat$p12m_npsv <- ifelse(dat$d125_ == 0, 0,
                   ifelse(dat$d124_ == 1, 1, 0))
  dat$ever_npsv <- ifelse(dat$d125_ == 1 & dat$d127_ !=1 & dat$d127_ != 2, 1, 0)
  
  dat <- dat[c("country", "surveyid", "survyear", "psu", "stratum", "age", "restype",
        "indweight", "vawweight", "selected_domvio",
        "cmu", "evermu", "eversex",
        "ever_ipv", "cmr_ipv", "p12m_ipv", 
        "p12m_severe", "ever_severe",
        "p12m_npsv", "ever_npsv")]
  return(dat)
}

datrc <- lapply(datlst, recode_dhs)

calc_cor <- function(dat,
                     agegr = c(15, 50),
                     weights = ~ vawweight) {

  if (is.null(dat)) {
    return(NULL) }
  dat <- dat[!is.na(dat[[all.vars(weights)]]) &
             !is.na(dat[[all.vars(~ ever_ipv)]]) &
             !is.na(dat[[all.vars(~ ever_npsv)]]), ]
  if (!nrow(dat)) {
        return(NULL) }

  # If all reports of are NA, it probably means that the question was not implemented correctly or that was not collected at all.
  if(sum(dat$ever_ipv, na.rm = TRUE) == 0) {
    return(NULL) }
  if(sum(dat$ever_npsv, na.rm = TRUE) == 0) {
    return(NULL) }
  
  dat$agegr <- cut(dat$age, agegr, demogsurv:::.epis_labels(agegr), TRUE, FALSE)
  dat$ipv <- ifelse(dat$evermu == 0, 0, dat$ever_ipv)
  dat$npsv <- dat$ever_npsv
  res <- weights::wtd.cor(dat$ipv, dat$npsv, weight = dat$indweight)
  
  if (length(unique(dat$psu)) == 1) { 
    design <- svydesign(ids = ~ 0, data = dat, strata = ~ stratum, weights = ~ vawweight, nest = TRUE) 
  } else {
    design <- svydesign(ids = ~ psu, data = dat, strata = ~ stratum, weights = ~ vawweight, nest = TRUE)    
  }
  
  val <- svyby(~ ipv, ~ surveyid + country + survyear + agegr,
               design, unwtd.count, na.rm = TRUE)
  val$cor <- as.data.frame(res)$correlation
  val$se <- as.data.frame(res)$std.err

  return(val)
}


calc_bin <- function(dat,
                     agegr = c(15, 50),
                     weights = ~ vawweight,
                     formula = ~ ever_ipv) {

  if(is.null(dat)) {
    return(NULL) }
  dat <- dat[!is.na(dat[[all.vars(weights)]]) &
             !is.na(dat[[all.vars(formula)]]), ]
  if(!nrow(dat)) {
        return(NULL) }

  # If all reports of are NA, it probably means that the question was not implemented correctly or that was not collected at all.
  if(sum(dat$ever_ipv, na.rm = TRUE) == 0) {
    return(NULL) }

  dat$agegr <- cut(dat$age, agegr, demogsurv:::.epis_labels(agegr), TRUE, FALSE)

  if(length(unique(dat$psu)) == 1) {
    des_cmu <- svydesign(ids = ~ 0, data = dat[dat$cmu == 1, ], strata = ~ stratum, weights = ~ vawweight, nest = TRUE)
    des_evermu <- svydesign(ids = ~ 0, data = dat[dat$evermu == 1, ], strata = ~ stratum, weights = ~ vawweight, nest = TRUE)
  } else {
    des_cmu <- svydesign(ids = ~ psu, data = dat[dat$cmu == 1, ], strata = ~ stratum, weights = ~ vawweight, nest = TRUE)
    des_evermu <- svydesign(ids = ~ psu, data = dat[dat$evermu == 1, ], strata = ~ stratum, weights = ~ vawweight, nest = TRUE)
  }

  val_cmu <- svyby(formula, ~ surveyid + country + survyear + agegr,
               des_cmu, svyciprop, vartype = c("se", "ci"), method = "logit", na.rm = TRUE)
  val_cmu$denom <- 'currently married or union'
  names(val_cmu) <- sub("^se[^x].*", "se", names(val_cmu))
  rownames(val_cmu) <- NULL

  cnt <- svyby(formula, ~ surveyid + country + survyear + agegr,
               des_cmu, unwtd.count, na.rm=TRUE)
  cnt$se <- NULL
  val_cmu <- merge(val_cmu, cnt)

  val_evermu <- svyby(formula, ~ surveyid + country + survyear + agegr,
               des_evermu, svyciprop, vartype = c("se", "ci"), method = "logit", na.rm = TRUE)
  val_evermu$denom <- 'ever married or union'
  names(val_evermu) <- sub("^se[^x].*", "se", names(val_evermu))
  rownames(val_evermu) <- NULL
  cnt <- svyby(formula, ~ surveyid + country + survyear + agegr,
               des_evermu, unwtd.count, na.rm = TRUE)
  cnt$se <- NULL
  val_evermu <- merge(val_evermu, cnt)

  val <- rbind(val_cmu, val_evermu)
  return(val)
}

# for severe in last 12 months
calc_bin_svr <- function(dat,
                     agegr = c(15, 50),
                     weights = ~ vawweight,
                     formula = ~ p12m_severe) {

  if(is.null(dat)) {
    return(NULL) }
  dat <- dat[!is.na(dat[[all.vars(weights)]]) &
             !is.na(dat[[all.vars(formula)]]), ]
  if(!nrow(dat)) {
        return(NULL) }

  # If all reports of are NA, it probably means that the question was not implemented correctly or that was not collected at all.
  if(sum(dat$p12m_severe, na.rm = TRUE) == 0) {
    return(NULL) }

  dat$agegr <- cut(dat$age, agegr, demogsurv:::.epis_labels(agegr), TRUE, FALSE)

  if(length(unique(dat$psu)) == 1) {
    des <- svydesign(ids = ~ 0, data = dat, strata = ~ stratum, weights = ~ vawweight, nest = TRUE)
  } else {
    des<- svydesign(ids = ~ psu, data = dat, strata = ~ stratum, weights = ~ vawweight, nest = TRUE)
  }

  val <- svyby(formula, ~ surveyid + country + survyear + agegr,
               des, svyciprop, vartype = c("se", "ci"), method = "logit", na.rm = TRUE)
  val$outcome <- as.character(formula)[2]
  names(val) <- sub("^se[^x].*", "se", names(val))
  names(val) <- sub("^p12m[^x].*", "est", names(val))
  rownames(val) <- NULL

  cnt <- svyby(formula, ~ surveyid + country + survyear + agegr,
               des, unwtd.count, na.rm = TRUE)
  cnt$se <- NULL
  val <- merge(val, cnt)

  return(val)
}

# for npsv
calc_bin_npsv <- function(dat,
                     agegr = c(15, 50),
                     weights = ~ vawweight,
                     formula = ~ p12m_npsv) {

  if(is.null(dat)) {
    return(NULL) }
  dat <- dat[!is.na(dat[[all.vars(weights)]]) &
             !is.na(dat[[all.vars(formula)]]), ]
  if(!nrow(dat)) {
        return(NULL) }

  # If all reports of are NA, it probably means that the question was not implemented correctly or that was not collected at all.
  if(sum(dat$p12m_severe, na.rm = TRUE) == 0) {
    return(NULL) }

  dat$agegr <- cut(dat$age, agegr, demogsurv:::.epis_labels(agegr), TRUE, FALSE)

  if(length(unique(dat$psu)) == 1) {
    des <- svydesign(ids = ~ 0, data = dat, strata = ~ stratum, weights = ~ vawweight, nest = TRUE)
  } else {
    des <- svydesign(ids = ~ psu, data = dat, strata = ~ stratum, weights = ~ vawweight, nest = TRUE)
  }

  val <- svyby(formula, ~ surveyid + country + survyear + agegr,
               des, svyciprop, vartype = c("se", "ci"), method = "logit", na.rm = TRUE)
  val$outcome <- as.character(formula)[2]
  names(val) <- sub("^se[^x].*", "se", names(val))
  names(val) <- sub("^npsv[^x].*", "est", names(val))
  rownames(val) <- NULL

  cnt <- svyby(formula, ~ surveyid + country + survyear + agegr,
               des, unwtd.count, na.rm = TRUE)
  cnt$se <- NULL
  val <- merge(val, cnt)

  return(val)
}


# correlation ipv & npsv
devtools::load_all(here::here())
library(parallel)
options(mc.cores = 4)
ipv_npsv <- mclapply(datrc, calc_cor)
cor_ipv_npsv <- do.call(rbind, ipv_npsv)
rownames(cor_ipv_npsv) <- NULL

cor_ipv_npsv$iso3 <- countrycode::countrycode(sourcevar = cor_ipv_npsv$country, origin = "country.name",
                                            destination = "iso3c")
usethis::use_data(cor_ipv_npsv, overwrite = TRUE)

# We pool the survey
library(parallel)
options(mc.cores = 12)

everipv <- mclapply(datrc, calc_bin, formula = ~ ever_ipv)
everipv <- do.call(rbind, everipv)
rownames(everipv) <- NULL

p12mipv <- mclapply(datrc, calc_bin, formula = ~ p12m_ipv)
p12mipv <- do.call(rbind, p12mipv)
rownames(p12mipv) <- NULL

# Severe vs. not svr ipv in p12m
p12m_svr_all <- mclapply(datrc, calc_bin_svr, formula = ~ p12m_ipv)
p12m_svr_all <- do.call(rbind, p12m_svr_all)
rownames(p12m_svr_all) <- NULL

p12m_svr_severe <- mclapply(datrc, calc_bin_svr, formula = ~ p12m_severe)
p12m_svr_severe <- do.call(rbind, p12m_svr_severe)
rownames(p12m_svr_severe) <- NULL

p12m_svr <- rbind(p12m_svr_all, p12m_svr_severe)

# Among 15-19
everipv1519 <- mclapply(datrc, calc_bin, formula = ~ ever_ipv, agegr = c(15, 20))
everipv1519 <- do.call(rbind, everipv1519)
rownames(everipv1519) <- NULL

p12mipv1519 <- mclapply(datrc, calc_bin, formula = ~ p12m_ipv, agegr = c(15, 20))
p12mipv1519 <- do.call(rbind, p12mipv1519)
rownames(p12mipv1519) <- NULL

everipv1519$prv <- everipv1519$ever_ipv
everipv1519$viotime <- "ever"
p12mipv1519$prv <- p12mipv1519$p12m_ipv
p12mipv1519$viotime <- "p12m"

ipv1519 <- rbind(everipv1519[everipv1519$denom == "ever married or union", !(names(everipv1519) %in% c("ever_ipv"))],
                 p12mipv1519[p12mipv1519$denom == "ever married or union", !(names(p12mipv1519) %in% c("p12m_ipv"))])

# npsv
ever_npsv <- mclapply(datrc, calc_bin_npsv, formula = ~ ever_npsv, agegr = c(15, 50))
ever_npsv <- do.call(rbind, ever_npsv)
rownames(ever_npsv) <- NULL

p12m_npsv <- mclapply(datrc, calc_bin_npsv, formula = ~ p12m_npsv, agegr = c(15, 50))
p12m_npsv <- do.call(rbind, p12m_npsv)
rownames(p12m_npsv) <- NULL

# We save in data()
xwalk_ever <- everipv
xwalk_p12m <- p12mipv
xwalk_p12msvr <- p12m_svr
xwalk_1519 <- ipv1519
xwalk_ever_npsv <- ever_npsv
xwalk_p12m_npsv <- p12m_npsv

# Inputting iso3 based on country name
devtools::load_all(here::here())

xwalk_ever$iso3 <- countrycode::countrycode(sourcevar = xwalk_ever$country, origin = "country.name",
                                            destination = "iso3c")
xwalk_p12m$iso3 <- countrycode::countrycode(sourcevar = xwalk_p12m$country, origin = "country.name",
                                            destination = "iso3c")
xwalk_p12msvr$iso3 <- countrycode::countrycode(sourcevar = xwalk_p12msvr$country, origin = "country.name",
                                            destination = "iso3c")
xwalk_1519$iso3 <- countrycode::countrycode(sourcevar = xwalk_1519$country, origin = "country.name",
                                            destination = "iso3c")
xwalk_ever_npsv$iso3 <- countrycode::countrycode(sourcevar = xwalk_ever_npsv$country, origin = "country.name",
                                            destination = "iso3c")
xwalk_p12m_npsv$iso3 <- countrycode::countrycode(sourcevar = xwalk_p12m_npsv$country, origin = "country.name",
                                            destination = "iso3c")
iso_gbd <- read.csv("./data-raw/table_gbd_iso3.csv", header = TRUE)

xwalk_ever$region <- NA
for (i in 1:nrow(xwalk_ever)) {
  xwalk_ever$region[i] <- as.character(iso_gbd$region[which(as.character(xwalk_ever$iso3[i]) == as.character(iso_gbd$iso3))])
}
xwalk_ever$SuperRegion <- class_reg_to_superreg(xwalk_ever$region)

xwalk_p12m$region <- NA
for (i in 1:nrow(xwalk_p12m)) {
  xwalk_p12m$region[i] <- as.character(iso_gbd$region[which(as.character(xwalk_p12m$iso3[i]) == as.character(iso_gbd$iso3))])
}
xwalk_p12m$SuperRegion <- class_reg_to_superreg(xwalk_p12m$region)

xwalk_p12msvr$region <- NA
for (i in 1:nrow(xwalk_p12msvr)) {
  xwalk_p12msvr$region[i] <- as.character(iso_gbd$region[which(as.character(xwalk_p12msvr$iso3[i]) == as.character(iso_gbd$iso3))])
}
xwalk_p12msvr$SuperRegion <- class_reg_to_superreg(xwalk_p12msvr$region)

xwalk_1519$region <- NA
for (i in 1:nrow(xwalk_1519)) {
  xwalk_1519$region[i] <- as.character(iso_gbd$region[which(as.character(xwalk_1519$iso3[i]) == as.character(iso_gbd$iso3))])
}
xwalk_1519$SuperRegion <- class_reg_to_superreg(xwalk_1519$region)

xwalk_ever_npsv$region <- NA
for (i in 1:nrow(xwalk_ever_npsv)) {
  xwalk_ever_npsv$region[i] <- as.character(iso_gbd$region[which(as.character(xwalk_ever_npsv$iso3[i]) == as.character(iso_gbd$iso3))])
}
xwalk_ever_npsv$SuperRegion <- class_reg_to_superreg(xwalk_ever_npsv$region)

xwalk_p12m_npsv$region <- NA
for (i in 1:nrow(xwalk_p12m_npsv)) {
  xwalk_p12m_npsv$region[i] <- as.character(iso_gbd$region[which(as.character(xwalk_p12m_npsv$iso3[i]) == as.character(iso_gbd$iso3))])
}
xwalk_p12m_npsv$SuperRegion <- class_reg_to_superreg(xwalk_p12m_npsv$region)


usethis::use_data(xwalk_ever, overwrite = TRUE)
usethis::use_data(xwalk_p12m, overwrite = TRUE)
usethis::use_data(xwalk_p12msvr, overwrite = TRUE)
usethis::use_data(xwalk_1519, overwrite = TRUE)
usethis::use_data(xwalk_ever_npsv, overwrite = TRUE)
usethis::use_data(xwalk_p12m_npsv, overwrite = TRUE)

# calculating design effect
length(unique(everipv$surveyid))

def <- everipv[everipv$denom == 'ever married or union', ]
def$s_srs <- def$ever_ipv * (1 - def$ever_ipv) / def$counts
def$def <- def$se^2 / def$s_srs

# not sure why but "AF2015DHS" has a really high DEFT not in line with that in its DHS report (p. 322).
def$def[def$surveyid == "AF2015DHS"] <- 3.010
# not sure why but "NG2013DHS" has a really high DEFT not in line with that in its DHS report (p. 388).
def$def[def$surveyid == "NG2013DHS"] <- 2.438

summary(def$def)
mean(def$def)
median(def$def)

length(unique(def$surveyid))


# For last year ipv
def <- p12mipv[p12mipv$denom == 'ever married or union', ]
def$s_srs <- def$p12m_ipv * (1 - def$p12m_ipv) / def$counts
def$def <- def$se^2 / def$s_srs

summary(def$def)
mean(def$def)
median(def$def)

length(unique(def$surveyid))

# for ever npsv
def <- ever_npsv
def$s_srs <- def$ever_npsv * (1 - def$ever_npsv) / def$counts
def$def <- def$se^2 / def$s_srs
def$def <- ifelse(def$def < 1, 1, def$def)

summary(def$def)
mean(def$def)
median(def$def)

length(unique(def$surveyid))

# for past year npsv
def <- p12m_npsv
def$s_srs <- def$p12m_npsv * (1 - def$p12m_npsv) / def$counts
def$def <- def$se^2 / def$s_srs
def$def <- ifelse(def$def < 1, 1, def$def)

summary(def$def)
mean(def$def)
median(def$def)

length(unique(def$surveyid))
