devtools::load_all(here::here())

# --- Reading Data Here ----
# data not included... not yet publicly released
Data <- read.csv('ipv_update_final_names_2020_11_12.csv', sep = ',', header = TRUE, encoding = "UTF-8")

Data <- subset(Data, startyr >= 2000 & 
                 genacis == "Not GENACIS study data" &
                 loqual == "Specific acts")
dim(Data)
Data$Age65 <- ifelse(Data$AgeM > 65, 65, Data$AgeM)

# Reclassify GBD regions in GBD super regions.
Data$SuperRegion <- class_reg_to_superreg(Data$region)


# ---- X-Walk for Severe Violence ----
IPV <- subset(Data, viotime == "Ever" & violence == "Physical and/or sexual IPV")
dim(IPV)


library(MatchIt)
IPVnona <- ipv_no_na(IPV); dim(IPVnona)
IPVout <- MatchIt::matchit(I(severe == 'Not only severe violence') ~ 
                    ID + pstat + spouseonly + currpart + geo + 
                      viotime + startyr + loage + hiageI, 
                  method = 'exact', data = IPVnona)
summary(IPVout)$n

C_Severe <- MatchIt::match.data(IPVout)
C_Severe$Study <- C_Severe$subclass

length(unique(C_Severe$iso3))
sum(C_Severe$denom_imp)


# Regression-based meta-analsyis
C_Severe$X <- ifelse(C_Severe$severe == 'Only severe violence', 1, 0)
fit <- lme4::glmer(cbind(Num, floor(denom_imp) - Num) ~  -1 + (1 | Study) + X + (-1 + X | Study),
           family = 'binomial', data = C_Severe)
summary(fit)
coef(summary(fit))

match_or <- ranef_sr(fit, C_Severe)
supreg_severe <- meta_sr(fit, match_or)
supreg_severe

# Standard meta-analysis
C_meta_ <- meta_proc(C_Severe, var = "severe",
                    event.name = "Only severe violence"); nrow(C_meta_)
C_meta <- sel_meta_unique(C_meta_, var = "severe"); nrow(C_meta)

MA <- meta::metabin(event.e = C_meta$ess.event.e, n.e = C_meta$ess.n.e, 
              event.c = C_meta$ess.event.c, n.c = C_meta$ess.n.c, 
              studlab = C_meta$author_year, sm = 'OR', hakn = FALSE,
              byvar = C_meta$SuperRegion)
summary(MA)

pdf('Meta Ever - Severe Violence.pdf', width = 9, height = 12)
  meta::forest(x = MA, fontsize = 6, sortvar = C_meta$SuperRegion,
               comb.random = TRUE, comb.fixed = FALSE,
              layout = "JAMA",
              subgroup = TRUE, lab.e = "Severe", lab.c = "All")
dev.off()

LogOR_Severe <- out_meta(MA)
LogOR_Severe



# ---- X-Walk for Physical ----
IPV <- subset(Data, viotime == 'Ever' & 
                (violence == 'Physical IPV only' | violence == "Physical and/or sexual IPV") &
                # We only perform the X-Walk on optimal set (reference case)
                IPV_exw == "All" &
                severe == 'Not only severe violence'); dim(IPV)

# Physical only first
IPVnona <- ipv_no_na(IPV); dim(IPVnona)
IPVout <- MatchIt::matchit(I(violence == 'Physical IPV only') ~ 
                    ID + pstat + spouseonly + currpart + geo + 
                    viotime + startyr + loage + hiageI, 
                  method = 'exact', data = IPVnona)
summary(IPVout)$nn

C_Violence_p <- MatchIt::match.data(IPVout)
C_Violence_p$Study <- C_Violence_p$subclass
sum(C_Violence_p$denom_imp)


# Standard meta-analysis
C_meta_ <- meta_proc(C_Violence_p, var = "violence",
                    event.name = "Physical IPV only"); nrow(C_meta_)
C_meta <- sel_meta_unique(C_meta_, var = "violence"); nrow(C_meta)

MAp <- meta::metabin(event.e = C_meta$ess.event.e, n.e = C_meta$ess.n.e, 
                     event.c = C_meta$ess.event.c, n.c = C_meta$ess.n.c, 
               studlab = C_meta$author_year, sm = 'OR', hakn = TRUE, byvar = C_meta$SuperRegion)
summary(MAp)

pdf('Meta Ever - Physical Violence.pdf', width = 9, height = 54)
  meta::forest(MAp, fontsize = 6, comb.fixed = FALSE, comb.random = TRUE, 
               layout = "JAMA",
               lab.e = "Physical", lab.c = "Physical/Sexual", subgroup = TRUE)
dev.off()

LogOR_Physical <- out_meta(MAp)
LogOR_Physical


# ---- X-Walk for Sexual ----
IPV <- subset(Data, viotime == 'Ever' & 
                (violence == 'Sexual IPV only' | violence == "Physical and/or sexual IPV") &
                # We only perform the X-Walk on optimal set (reference case)
                IPV_exw == "All" &
                severe == 'Not only severe violence'); dim(IPV)

library(MatchIt)
IPVnona <- ipv_no_na(IPV); dim(IPVnona)
IPVout <- MatchIt::matchit(I(violence == 'Sexual IPV only') ~ 
                             ID + pstat + spouseonly + currpart + geo + 
                             viotime + startyr + loage + hiageI, 
                  method = 'exact', data = IPVnona)
summary(IPVout)$nn

C_Violence_s <- MatchIt::match.data(IPVout)
C_Violence_s$Study <- C_Violence_s$subclass
sum(C_Violence_s$denom_imp)

# Standard meta-analysis
C_meta_ <- meta_proc(C_Violence_s, var = "violence",
                    event.name = "Sexual IPV only"); nrow(C_meta_)
C_meta <- sel_meta_unique(C_meta_, var = "violence"); nrow(C_meta)

MAs <- meta::metabin(event.e = C_meta$ess.event.e, n.e = C_meta$ess.n.e, 
                     event.c = C_meta$ess.event.c, n.c = C_meta$ess.n.c, 
               studlab = C_meta$author_year, sm = 'OR', hakn = TRUE, byvar = C_meta$SuperRegion)
summary(MAs)

pdf('Meta Ever - Sexual Violence.pdf', width = 9, height = 54)
  meta::forest(MAs, fontsize = 6, comb.fixed = FALSE, comb.random = TRUE,
               layout = "JAMA", lab.e = "Sexual", lab.c = "Physical/Sexual",
               subgroup = TRUE)
dev.off()

LogOR_Sexual <- out_meta(MAs)
LogOR_Sexual


# ---- X-Walk for All Women (Denom) ----
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
                 studlab = C_meta$author_year, sm = 'OR', hakn = TRUE, byvar = C_meta$SuperRegion)
summary(MA_aw)

pdf('Meta Ever - Denom All Women.pdf', width = 9, height = 8)
  meta::forest(MA_aw, fontsize = 6, sortvar = C_meta$SuperRegion, comb.fixed = FALSE, comb.random = TRUE, 
               layout = "JAMA", subgroup = TRUE, lab.e = "All women", lab.c = "Ever-Partnered")
dev.off()

LogOR_AllW <- out_meta(MA_aw)
LogOR_AllW


# ---- X-Walk for Currently Partnered (Denom) ----
IPV <- subset(Data, viotime == 'Ever' & violence == "Physical and/or sexual IPV" &
                pstat != "All women"); dim(IPV) 

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

# We get the citations for the DHS studies

# Much more power if we use DHS. (only 32 studies otherwise)
data("xwalk_ever")
dhs <- dplyr::distinct(Data[Data$author == "DHS", c("iso3", "startyr", "author_year")])

# Fix the years and citations
colnames(xwalk_ever)[3] <- "startyr"

xwalk_ever <- merge(xwalk_ever, dhs, by = c("iso3", "startyr"), all.x = T)
xwalk_ever$author_year <- as.character(xwalk_ever$author_year)
xwalk_ever$author_year <- ifelse(is.na(xwalk_ever$author_year), 
                                 sprintf("%s - (DHS, %s)",xwalk_ever$country, xwalk_ever$startyr),
                                 xwalk_ever$author_year)
xwalk_ever$ess <- neff(xwalk_ever$ever_ipv, xwalk_ever$ci_l, xwalk_ever$ci_u, wilson = TRUE)
xwalk_ever$Num <- xwalk_ever$num_ess <- round(xwalk_ever$ever_ipv * xwalk_ever$counts, 0 )
xwalk_ever$denom_imp <- xwalk_ever$den_ess <- xwalk_ever$counts

library(MatchIt)
IPVout <- MatchIt::matchit(I(denom == 'currently married or union') ~ surveyid, 
                    method = 'exact', data = xwalk_ever)
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

pdf('Meta Ever - Denom Currently Partnered.pdf', width = 9, height = 25)
  meta::forest(MA_cp, fontsize = 6, sortvar = C_meta$SuperRegion, comb.fixed = FALSE, comb.random = TRUE,
               layout = "JAMA", subgroup = TRUE, lab.e = "Currently Partnered", lab.c = "Ever-Partnered")
dev.off()

LogOR_CurP <- out_meta(MA_cp)
LogOR_CurP

# ---- X-Walk for Currpart ----
IPV <- subset(Data, viotime == 'Ever' & 
                violence == "Physical and/or sexual IPV"); dim(IPV)

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

pdf('Meta Ever - Violence from CMR partner only.pdf', width = 9, height = 15)
  meta::forest(MA_cmr, fontsize = 6, sortvar = C_meta$SuperRegion, comb.fixed = FALSE, comb.random = TRUE, subgroup = TRUE,
               layout = "JAMA", lab.e = "Current Partner Only", lab.c = "Any Partner")
dev.off()

LogOR_CMRpartner <- out_meta(MA_cmr)
LogOR_CMRpartner






# ---- X-Walk for Urban/Rural ----
# For Geo, we need to bring back the original database where we don't remove rural/urban if national is available.

Data_UR <- read.csv('ipv_urban_rural_update_final_names_2020_11_12.csv', sep  = ',', header = TRUE); length(unique(Data_UR$ID)); sum(Data_UR$denom_imp)
Data_UR <- subset(Data_UR, startyr >= 2000 & 
                 genacis == "Not GENACIS study data" &
                 loqual == "Specific acts")
dim(Data_UR)
Data_UR$Age65 <- ifelse(Data_UR$AgeM > 65, 65, Data_UR$AgeM)
Data_UR$geo <- relevel(Data_UR$geo, ref = 'National')

Data_UR$SuperRegion <- class_reg_to_superreg(Data_UR$region)

#' ----  Regression ----
IPV <- subset(Data_UR, viotime == 'Ever' & 
                violence == "Physical and/or sexual IPV" &
                geo != 'Mixed'); dim(IPV)

IPVnona <- ipv_no_na(IPV); dim(IPVnona)

# both together
# we match only if we have all three.
my_match <- match_geo(IPVnona)
dim(my_match); length(unique(my_match$subclass))

C_meta <- sel_meta_unique(my_match, var = "geo"); dim(C_meta)

fit_geo <- glm(cbind(Num, denom_imp - Num) ~  I(geo == 'Rural') + I(geo == 'Urban') + as.factor(ID), 
                   family = 'binomial', data = C_meta)
exp(coef(fit_geo)[2:3])

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


# Standard meta-analysis
# ----    Urban ----
C_data_urb <- sel_meta_unique(my_match[my_match$geo != "Rural", ], var = "geo")
C_meta_urb <- meta_proc_geo(C_data_urb, var = "geo", event.name = "Urban")
MA_urb <- meta::metabin(event.e = C_meta_urb$event.e, n.e = C_meta_urb$n.e, 
                        event.c = C_meta_urb$event.c, n.c = C_meta_urb$n.c, 
                 studlab = C_meta_urb$author_year, sm = 'OR', hakn = TRUE, byvar = C_meta_urb$SuperRegion)
summary(MA_urb)

pdf('Meta Ever - Urban.pdf', width = 10, height = 30)
  meta::forest(MA_urb, fontsize = 6, sortvar = C_meta_urb$SuperRegion, comb.fixed = FALSE, comb.random = FALSE,
               layout = "JAMA", subgroup = TRUE, lab.e = "Urban", lab.c = "National")
dev.off()


# ---- Rural ----
C_data_rur <- sel_meta_unique(my_match[my_match$geo != "Urban", ], var = "geo")
C_meta_rur <- meta_proc_geo(C_data_rur, var = "geo", event.name = "Rural")
MA_rur <- meta::metabin(event.e = C_meta_rur$event.e, n.e = C_meta_rur$n.e, 
                        event.c = C_meta_rur$event.c, n.c = C_meta_rur$n.c, 
                 studlab = C_meta_rur$author_year, sm = 'OR', hakn = TRUE, byvar = C_meta_rur$SuperRegion)
summary(MA_rur)

pdf('Meta Ever - Rural.pdf', width = 10, height = 30)
  meta::forest(MA_rur, fontsize = 6, sortvar = C_meta_rur$SuperRegion, comb.fixed = FALSE, comb.random = FALSE,
               layout = "JAMA",
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
save(LogOR_Severe, LogOR_Physical, LogOR_Sexual, LogOR_AllW, LogOR_CurP,
    LogOR_CMRpartner, LogOR_Rural, LogOR_Urban,
     file = "LogOR-Ever 2020-10-30.RData")


  
  