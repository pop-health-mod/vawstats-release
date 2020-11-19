devtools::load_all(here::here())

# --- Reading Data Here ----
# data not included... not yet publicly released
Data <- read.csv('ipv_update_final_names_2020_11_12.csv', sep = ',', header = TRUE)
Data <- subset(Data, startyr >= 2000 &
                 genacis == "Not GENACIS study data" &
                 loqual == "Specific acts")
dim(Data)
Data$Age65 <- ifelse(Data$AgeM > 65, 65, Data$AgeM)

# Reclassify GBD regions in GBD super regions.
Data$SuperRegion <- class_reg_to_superreg(Data$region)


# ---- X-Walk for Severe Violence ----
IPV <- subset(Data, viotime == "Past yr" & # violence == "Physical and/or sexual IPV" &
                pstat == "Ever-partnered"); dim(IPV)

library(MatchIt)
IPVnona <- ipv_no_na(IPV); dim(IPVnona)
IPVout <- MatchIt::matchit(I(severe == 'Only severe violence') ~ 
                      ID + violence + pstat + spouseonly + currpart + geo + 
                      viotime + startyr + loage + hiageI, 
                  method = 'exact', data = IPVnona)
summary(IPVout)$nn

C_Severe <- MatchIt::match.data(IPVout)
C_Severe$Study <- C_Severe$subclass
length(unique(C_Severe$iso3))
sum(C_Severe$denom_imp)

# Standard meta-analysis
C_meta_ <- meta_proc(C_Severe, var = "severe",
                    event.name = "Only severe violence"); nrow(C_meta_)
C_meta <- sel_meta_unique(C_meta_, var = "severe"); nrow(C_meta)

MA <- meta::metabin(event.e = C_meta$ess.event.e, n.e = C_meta$ess.n.e, 
              event.c = C_meta$ess.event.c, n.c = C_meta$ess.n.c, 
              studlab = C_meta$author_year, sm = 'OR', hakn = FALSE,
              byvar = C_meta$SuperRegion)
summary(MA)

pdf("Meta Past Yr - Severe Violence ERROR.pdf", width = 10, height = 15)
  meta::forest(MA, fontsize = 6, sortvar = C_meta$SuperRegion, comb.random = FALSE, subgroup = TRUE, lab.e = "Severe", lab.c = "All")
dev.off()


# We redo with DHS data only
data("xwalk_p12msvr")
dhs <- dplyr::distinct(Data[Data$author == "DHS", c("iso3", "startyr", "author_year")])

# Fix the years and citations
colnames(xwalk_p12msvr)[3] <- "startyr"

xwalk_p12msvr <- merge(xwalk_p12msvr, dhs, by = c("iso3", "startyr"), all.x = T)
xwalk_p12msvr$author_year <- as.character(xwalk_p12msvr$author_year)
xwalk_p12msvr$author_year <- ifelse(is.na(xwalk_p12msvr$author_year), 
                                 sprintf("%s - (DHS, %s)",xwalk_p12msvr$country, xwalk_p12msvr$startyr),
                                 xwalk_p12msvr$author_year)
xwalk_p12msvr$ess <- neff(xwalk_p12msvr$est, xwalk_ever$ci_l, xwalk_p12msvr$ci_u, wilson = TRUE)
xwalk_p12msvr$Num <- xwalk_p12msvr$num_ess <- round(xwalk_p12msvr$est * xwalk_p12msvr$counts, 0 )
xwalk_p12msvr$denom_imp <- xwalk_p12msvr$den_ess <- xwalk_p12msvr$counts

library(MatchIt)
IPVout <- matchit(I(outcome == 'p12m_severe') ~ surveyid, 
                    method = 'exact', data = xwalk_p12msvr)
summary(IPVout)$nn

C_Severe <- MatchIt::match.data(IPVout)
C_Severe$Study <- C_Severe$subclass
length(unique(C_Severe$iso3))
sum(C_Severe$denom_imp)

# Standard meta-analysis
C_meta <- meta_proc(C_Severe, var = "outcome",
                     event.name = "p12m_severe",
                     is.dhs = TRUE); nrow(C_meta)

MAsvr <- meta::metabin(event.e = C_meta$ess.event.e, n.e = C_meta$ess.n.e, 
              event.c = C_meta$ess.event.c, n.c = C_meta$ess.n.c, 
              studlab = C_meta$author_year, sm = 'OR', hakn = FALSE,
              byvar = C_meta$SuperRegion)
summary(MAsvr)

pdf("Meta Past Yr - Severe Violence DHS.pdf", width = 10, height = 25)
  meta::forest(MAsvr, fontsize = 6, sortvar = C_meta$SuperRegion,
               comb.random = TRUE, comb.fixed = FALSE, subgroup = TRUE,
               layout = "JAMA", lab.e = "Severe", lab.c = "All")
dev.off()


LogOR_Severe <- out_meta(MAsvr)
LogOR_Severe


# ---- X-Walk for Past 2 Year ----
IPV <- subset(Data, parametertype == "Past year IPV" & violence == "Physical and/or sexual IPV" &
              viotime != "Past 5 yrs"); dim(IPV)

# Not required... all in optimal set.
LogOR_Past2Yr <- NA


# ---- X-walk for Physical ----
IPV <- subset(Data, viotime == "Past yr" & 
             (violence == 'Physical IPV only' | violence == "Physical and/or sexual IPV") &
                # We only perform the X-Walk on optimal set (reference case)
                IPV_exw == "All" &
                severe == 'Not only severe violence'); dim(IPV)
dim(IPV)

# Phyisical only first
library(MatchIt)
IPVnona <- ipv_no_na(IPV); dim(IPVnona)
IPVout <- MatchIt::matchit(I(violence == "Physical IPV only") ~ 
                    ID + pstat + spouseonly + currpart + geo +
                    viotime + startyr + loage + hiageI, 
                method = "exact", data = IPVnona)
summary(IPVout)$nn

C_Violence_p <- match.data(IPVout)
C_Violence_p$Study <- C_Violence_p$subclass
sum(C_Violence_p$denom_imp)

# Standard meta-analysis
C_meta_p <- meta_proc(C_Violence_p, var = "violence",
                    event.name = "Physical IPV only"); nrow(C_meta_p)
C_meta <- sel_meta_unique(C_meta_p, var = "violence"); nrow(C_meta)

MAp <- meta::metabin(event.e = C_meta$ess.event.e, n.e = C_meta$ess.n.e, 
                     event.c = C_meta$ess.event.c, n.c = C_meta$ess.n.c, 
               studlab = C_meta$author_year, sm = 'OR', hakn = TRUE, byvar = C_meta$SuperRegion)
summary(MAp)

pdf('Meta Past Yr - Physical Violence.pdf', width = 9, height = 52)
  meta::forest(MAp, fontsize = 6, comb.fixed = FALSE, comb.random = TRUE, 
               layout = "JAMA", lab.e = "Physical", lab.c = "Physical/Sexual", subgroup = TRUE)
dev.off()

LogOR_Physical <- out_meta(MAp)
LogOR_Physical



# ---- X-walk for Sexual ----
IPV <- subset(Data, viotime == "Past yr" & 
             (violence == 'Sexual IPV only' | violence == "Physical and/or sexual IPV") &
                # We only perform the X-Walk on optimal set (reference case)
                IPV_exw == "All" &
                severe == 'Not only severe violence'); dim(IPV)
dim(IPV)

# Sexual only 
library(MatchIt)
IPVnona <- ipv_no_na(IPV); dim(IPVnona)
IPVout <- MatchIt::matchit(I(violence == "Sexual IPV only") ~ 
                    ID + pstat + spouseonly + currpart + geo +
                    viotime + startyr + loage + hiageI, 
                method = "exact", data = IPVnona)
summary(IPVout)$nn

C_Violence_s <- match.data(IPVout)
C_Violence_s$Study <- C_Violence_s$subclass
sum(C_Violence_s$denom_imp)

# Standard meta-analysis
C_meta_s <- meta_proc(C_Violence_s, var = "violence",
                    event.name = "Sexual IPV only"); nrow(C_meta_s)
C_meta <- sel_meta_unique(C_meta_s, var = "violence"); nrow(C_meta)

MAs <- meta::metabin(event.e = C_meta$ess.event.e, n.e = C_meta$ess.n.e, 
                     event.c = C_meta$ess.event.c, n.c = C_meta$ess.n.c, 
               studlab = C_meta$author_year, sm = 'OR', hakn = TRUE, byvar = C_meta$SuperRegion)
summary(MAs)

pdf('Meta Past Yr - Sexual Violence.pdf', width = 10, height = 52)
  meta::forest(MAs, fontsize = 6, comb.fixed = FALSE, comb.random = TRUE, 
               layout = "JAMA", lab.e = "Sexual", lab.c = "Physical/Sexual", subgroup = TRUE)
dev.off()

LogOR_Sexual <- out_meta(MAs)
LogOR_Sexual


# ---- X-Walk for All Women (Denom) ----
IPV <- subset(Data, viotime == "Past yr" & parametertype == "Past year IPV" &
                violence != "Emotional abuse only" &
                pstat != "Currently partnered only"); dim(IPV)
# Here, we need to include all response type (including severe), otherwhise, not enough variation.

summary(IPV$pstat)

# All Women vs. Ever-Partnered
library(MatchIt)
IPVnona <- ipv_no_na(IPV); dim(IPVnona)
IPVout_aw <- matchit(I(pstat == 'All women') ~ 
                    ID + spouseonly + geo + severe + violence +
                    startyr + loage + hiageI + currpart, 
                  method = 'exact', data = IPVnona)
summary(IPVout_aw)$nn

C_AllWom <- match.data(IPVout_aw)
C_AllWom$Study <- C_AllWom$subclass
length(unique(C_AllWom$iso3))
sum(C_AllWom$denom_imp)

# Standard meta-analysis
C_meta_ <- meta_proc(C_AllWom, var = "pstat",
                     event.name = "All women"); nrow(C_meta_)
C_meta <- sel_meta_unique(C_meta_, var = "pstat"); nrow(C_meta)

MA_aw <- meta::metabin(event.e = C_meta$ess.event.e, n.e = C_meta$ess.n.e, 
                       event.c = C_meta$ess.event.c, n.c = C_meta$ess.n.c, 
                       studlab = C_meta$ID, sm = 'OR', hakn = TRUE, byvar = C_meta$SuperRegion)
summary(MA_aw)

pdf('Meta Past Yr - Denom All Women - ERROR.pdf', width = 10, height = 10)
  meta::forest(MA_aw, fontsize = 6, comb.fixed = FALSE, comb.random = TRUE, 
               lab.e = "All women", lab.c = "Ever-partnered", subgroup = TRUE)
dev.off()


# ### === ### === ### === ### === ### === ### === ### === ### ===
# We will use the OR from Ever IPV to adjust the denominator.
# ### === ### === ### === ### === ### === ### === ### === ### ===
IPV <- subset(Data, viotime == 'Ever' & 
                violence == "Physical and/or sexual IPV" &
                pstat != "Currently partnered only"); dim(IPV)
# Here, we need to include all response type (including severe), otherwhise, not enough variation.

# All Women vs. Ever-Partnered
library(MatchIt)
IPVnona <- ipv_no_na(IPV); dim(IPVnona)
IPVout <- matchit(I(pstat == 'All women') ~ 
                    ID + spouseonly + currpart + geo + severe +
                    viotime + startyr + loage + hiageI, 
                  method = 'exact', data = IPVnona)
summary(IPVout)$nn

C_AllWom <- match.data(IPVout)
C_AllWom$Study <- C_AllWom$subclass
length(unique(C_AllWom$iso3))
sum(C_AllWom$denom_imp)

# Standard meta-analysis
C_meta_ <- meta_proc(C_AllWom, var = "pstat",
                     event.name = "All women"); nrow(C_meta_)
C_meta <- sel_meta_unique(C_meta_, var = "pstat"); nrow(C_meta)

MA_aw <- meta::metabin(event.e = C_meta$ess.event.e, n.e = C_meta$ess.n.e, 
                       event.c = C_meta$ess.event.c, n.c = C_meta$ess.n.c, 
                       studlab = C_meta$ID, sm = 'OR', hakn = TRUE, byvar = C_meta$SuperRegion)
meta::forest(MA_aw)
summary(MA_aw)

# We don't plot because we are using the estimate from Ever IPV.

LogOR_AllW <- out_meta(MA_aw)
LogOR_AllW


# ---- X-Walk for Currently Partnered (Denom) ----
IPV <- subset(Data, viotime == "Past yr" & violence == "Physical and/or sexual IPV" &
                pstat != "All women"); dim(IPV) 
# Here, we need to include all response type (including severe), otherwhise, not enough variation.

summary(IPV$pstat)

library(MatchIt)
IPVnona <- ipv_no_na(IPV); dim(IPVnona)
IPVout <- matchit(I(pstat == 'Currently partnered only') ~ 
                    ID + currpart + spouseonly + geo + severe +
                    viotime + startyr + loage + hiageI, 
                  method = 'exact', data = IPVnona)
summary(IPVout)$nn

C_CurPar <- match.data(IPVout)
C_CurPar$Study <- C_CurPar$subclass
length(unique(C_CurPar$iso3))
sum(C_CurPar$denom_imp)

# Not systematically abstracted from DHS study so we use the DHS microdata
data("xwalk_p12m")
dhs <- dplyr::distinct(Data[Data$author == "DHS", c("iso3", "startyr", "author_year")])

# Fix the years and citations
colnames(xwalk_p12m)[3] <- "startyr"

xwalk_p12m <- merge(xwalk_p12m, dhs, by = c("iso3", "startyr"), all.x = T)
xwalk_p12m$author_year <- as.character(xwalk_p12m$author_year)
xwalk_p12m$author_year <- ifelse(is.na(xwalk_p12m$author_year), 
                                 sprintf("%s - (DHS, %s)",xwalk_p12m$country, xwalk_p12m$startyr),
                                 xwalk_p12m$author_year)
xwalk_p12m$ess <- neff(xwalk_p12m$p12m, xwalk_p12m$ci_l, xwalk_p12m$ci_u, wilson = TRUE)
xwalk_p12m$Num <- xwalk_p12m$num_ess <- round(xwalk_p12m$p12m * xwalk_p12m$counts, 0 )
xwalk_p12m$denom_imp <- xwalk_p12m$den_ess <- xwalk_p12m$counts

library(MatchIt)
IPVout <- MatchIt::matchit(I(denom == 'currently married or union') ~ surveyid, 
                           method = 'exact', data = xwalk_p12m)
summary(IPVout)$nn

C_CurPar <- match.data(IPVout)
C_CurPar$Study <- C_CurPar$subclass

# Standard meta-analysis
C_meta <- meta_proc(C_CurPar, var = "denom",
                    event.name = "currently married or union",
                    is.dhs = TRUE); nrow(C_meta)

MA_cp <- meta::metabin(event.e = C_meta$ess.event.e, n.e = C_meta$ess.n.e, 
                       event.c = C_meta$ess.event.c, n.c = C_meta$ess.n.c, 
                       studlab = C_meta$author_year, sm = 'OR', hakn = TRUE, byvar = C_meta$SuperRegion)
summary(MA_cp)

pdf('Meta Past Yr - Denom Currently Partnered.pdf', width = 9, height = 25)
meta::forest(MA_cp, fontsize = 6, sortvar = C_meta$SuperRegion, comb.fixed = FALSE, comb.random = TRUE,
             layout = "JAMA", subgroup = TRUE, lab.e = "Currently Partnered", lab.c = "Ever-Partnered")
dev.off()

LogOR_CurP <- out_meta(MA_cp)
LogOR_CurP



# ---- X-Walk for Currpart ----
IPV <- subset(Data, viotime == "Past yr" & violence == "Physical and/or sexual IPV"); dim(IPV)

library(MatchIt)
IPVnona <- ipv_no_na(IPV); dim(IPVnona)
IPVnona$currpart <- ifelse(IPVnona$currpart != "Not only asked violence from current partner", 
                           "Asked about violence from current or most recent partner only", "Not only asked violence from current partner")
IPVout <- MatchIt::matchit(I(currpart == 'Asked about violence from current or most recent partner only') ~ 
                             ID + spouseonly + pstat + geo + severe +
                             viotime + startyr + loage + hiageI, 
                           method = 'exact', data = IPVnona)
summary(IPVout)$nn

C_cmr <- match.data(IPVout)
C_cmr$Study <- C_cmr$subclass
length(unique(C_cmr$iso3)); length(unique(C_cmr$ID))
sum(C_cmr$denom_imp)

# Standard meta-analysis
C_meta_ <- meta_proc(C_cmr, var = "currpart",
                     event.name = "Asked about violence from current or most recent partner only"); nrow(C_meta_)
C_meta <- sel_meta_unique(C_meta_, var = "currpart"); nrow(C_meta)

MA_cmr <- meta::metabin(event.e = C_meta$ess.event.e, n.e = C_meta$ess.n.e, 
                        event.c = C_meta$ess.event.c, n.c = C_meta$ess.n.c, 
                        studlab = C_meta$author_year, sm = 'OR', hakn = TRUE, byvar = C_meta$SuperRegion)
summary(MA_cmr)

pdf('Meta Past Yr - Violence from CMR partner only.pdf', width = 9, height = 15)
meta::forest(MA_cmr, fontsize = 6, sortvar = C_meta$SuperRegion, comb.fixed = FALSE, comb.random = TRUE,
             subgroup = TRUE, layout = "JAMA",
             lab.e = "Current Partner Only", lab.c = "Any Partner")
dev.off()

LogOR_CMRpartner <- out_meta(MA_cmr)
LogOR_CMRpartner



# ---- X-Walk for Urban/Rural ----
# For Geo, we need to bring back the original database where we don't remove rural/urban if national is available.
#devtools::load_all("~/Google Drive/McGill/Research/VAW Statistics/vawstats/R")
#setwd("~/Google Drive/McGill/Research/VAW Statistics/Data")

Data_UR <- read.csv('ipv_urban_rural_update_final_names_2020_11_12.csv', sep  = ',', header = TRUE); length(unique(Data_UR$ID)); sum(Data_UR$denom_imp)
Data_UR <- subset(Data_UR, startyr >= 2000 & 
                 genacis == "Not GENACIS study data" &
                 loqual == "Specific acts")
dim(Data_UR)
Data_UR$Age65 <- ifelse(Data_UR$AgeM > 65, 65, Data_UR$AgeM)
Data$geo <- relevel(Data$geo, ref = 'National')

Data_UR$SuperRegion <- class_reg_to_superreg(Data_UR$region)



#' ----  Regression ----
IPV <- subset(Data_UR, viotime == 'Past yr' & 
                violence == "Physical and/or sexual IPV" &
                geo != 'Mixed'); dim(IPV)

IPVnona <- ipv_no_na(IPV); dim(IPVnona)

# both together
# we match only if we have all three.
my_match <- match_geo(IPVnona)
dim(my_match); length(unique(my_match$subclass))

C_meta <- sel_meta_unique(my_match, var = "geo"); dim(C_meta)

# Random effect
C_meta$rural <- ifelse(C_meta$geo == "Rural", 1 , 0)
C_meta$urban <- ifelse(C_meta$geo == "Urban", 1 , 0)

# Regression based random effect meta-analysis
fit <- lme4::glmer(cbind(Num, denom_imp - Num) ~ (1 | ID) + 
                     rural + urban +
                    (rural + urban - 1 | SuperRegion),
             family = 'binomial', data = C_meta)
coef(fit)$SuperRegion[, c("urban","rural")]
coef(summary(fit))[c("urban", "rural"), ]

LogOR_Urban <- out_meta_glmer_geo(fit, geo = "urban")
LogOR_Urban

LogOR_Rural <- out_meta_glmer_geo(fit, geo = "rural") 
LogOR_Rural


# ---- Urban ----
C_data_urb <- sel_meta_unique(my_match[my_match$geo != "Rural", ], var = "geo")
C_meta_urb <- meta_proc_geo(C_data_urb, var = "geo", event.name = "Urban")
MA_urb <- meta::metabin(event.e = C_meta_urb$ess.event.e, n.e = C_meta_urb$ess.n.e, 
                        event.c = C_meta_urb$ess.event.c, n.c = C_meta_urb$ess.n.c, 
                        studlab = C_meta_urb$author_year, sm = 'OR', hakn = TRUE,
                        byvar = C_meta_urb$SuperRegion)
summary(MA_urb)

pdf('Meta Past Yr - Urban.pdf', width = 10, height = 26)
meta::forest(MA_urb, fontsize = 6, sortvar = C_meta_urb$SuperRegion, comb.fixed = FALSE,
             comb.random = FALSE, layout = "JAMA",
             subgroup = TRUE, lab.e = "Urban", lab.c = "National")
dev.off()


# ---- Rural ----
C_data_rur <- sel_meta_unique(my_match[my_match$geo != "Urban", ], var = "geo")
C_meta_rur <- meta_proc_geo(C_data_rur, var = "geo", event.name = "Rural")
MA_rur <- meta::metabin(event.e = C_meta_rur$ess.event.e, n.e = C_meta_rur$ess.n.e, 
                        event.c = C_meta_rur$ess.event.c, n.c = C_meta_rur$ess.n.c, 
                        studlab = C_meta_rur$author_year, sm = 'OR', hakn = TRUE,
                        byvar = C_meta_rur$SuperRegion)
summary(MA_rur)

pdf('Meta Past Yr - Rural.pdf', width = 10, height = 26)
meta::forest(MA_rur, fontsize = 6, sortvar = C_meta_rur$SuperRegion, comb.fixed = FALSE,
             comb.random = FALSE, layout = "JAMA",
             subgroup = TRUE, lab.e = "Rural", lab.c = "National")
dev.off()


# Checking if all consistent
MA_rur$data$.studlab[which(MA_rur$TE > 0 & MA_urb$TE > 0)]
  exp(MA_rur$TE)[which(MA_rur$TE > 0 & MA_urb$TE > 0)]
  exp(MA_urb$TE)[which(MA_rur$TE > 0 & MA_urb$TE > 0)]
MA_rur$data$.studlab[which(MA_rur$TE < 0 & MA_urb$TE < 0)]
  exp(MA_urb$TE)[which(MA_rur$TE < 0 & MA_urb$TE < 0)]
  exp(MA_rur$TE)[which(MA_rur$TE < 0 & MA_urb$TE < 0)]
  

# ---- Saving Outputs! ----
save(LogOR_Severe, LogOR_Physical, LogOR_Sexual, LogOR_AllW, LogOR_CMRpartner, LogOR_CurP,
     LogOR_Rural, LogOR_Urban,
    file = "LogOR-PastYr 2020-10-30.RData")







