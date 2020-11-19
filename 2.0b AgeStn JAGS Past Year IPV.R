devtools::load_all(here::here())

# ---- Load the data ----
# data not included... not yet publicly released
data(wom_by_age_and_region_2010)
data(vwa_data) 

# Rename.
Data <- vwa_data; dim(Data)

# ---- Past Yr IPV ----
IPV_ <- subset(Data, parametertype == "Past year IPV" &
                      violence != "Emotional abuse only" &
                      violence != "Any sexual violence only" &
                      violence != "Nonpartner sexual violence only" &
                      genacis == "Not GENACIS study data" &
                      loqual == "Specific acts" &
                      startyr >= 2000 & endyr <= 2018); dim(IPV_)
length(unique(IPV_$ID)); length(unique(IPV_$iso3))

# We remove NA
any(is.na(IPV_))
IPV_ <- na.omit(subset(IPV_, select = c(ID, country, region, SuperRegion, iso3, startyr, agerange, Time,
                                        num_ess, den_ess, LCI, UCI, Prv, loage, hiageI, Age65, AgeW, AgeM, denom_imp, 
                                        parametertype, pstat, currpart, viotime, violence, severe, 
                                        loqual, notviostudy1, nointrain,
                                        spouseonly, currpart, geo, imputation_denom))); dim(IPV_)          
any(is.na(IPV_))

# Population covered
data(denom_cnt_age_now)
pop_cov <- subset(denom_cnt_age_now, iso3 %in% unique(IPV_$iso3))
round(sum(pop_cov$women) / sum(denom_cnt_age_now$women) * 100, 0)

# ---- Selecting Optimal Set ----
IPV_$IDage <- paste(IPV_$ID, IPV_$agerange, sep = '_')
nrow(IPV_); length(unique(IPV_$ID)); length(unique(IPV_$IDage))
sum(IPV_$denom_imp)

IPV <- optimal_set(IPV_, outcome = "Past Year IPV")
nrow(IPV); 
length(unique(IPV$ID)); 
length(unique(IPV$IDage))
sum(IPV$denom_imp)
IPV1 <- IPV

# We remove larger age-group if age-disaggregated is available
IPV <- remove_larger_agegrp(IPV)
nrow(IPV); length(unique(IPV$ID))
length(unique(IPV$iso3))
sum(IPV$denom_imp)

# What we removed
ID1 <- as.character(unique(IPV1$ID))
ID2 <- as.character(unique(IPV$ID))
ID1[!(ID1 %in% ID2)] # Should be zero; we don't exclude any study.

# We examine all variables (watchout for NAs)
summary(IPV$violence)
summary(IPV$viotime)
summary(IPV$severe)
summary(IPV$pstat)
summary(IPV$spouseonly)
summary(IPV$currpart)
summary(IPV$geo)
summary(IPV$nointrain)
sum(IPV$nointrain  == "Interviewers not trained") / nrow(IPV)
length(unique(IPV$ID[which(IPV$nointrain  == "Interviewers not trained")])) / length(unique(IPV$ID))
summary(IPV$loqual)
summary(IPV$notviostudy1)
summary(IPV$notviostudy1[!duplicated(IPV$ID)])

# We order
IPV <- IPV[order(IPV$region, IPV$iso3, IPV$startyr, IPV$ID, IPV$violence, 
                 IPV$severe, IPV$pstat, IPV$geo, IPV$viotime, IPV$AgeM), ]
IPV$RowID <- seq(1, dim(IPV)[1], 1)

I <- length(unique(IPV$ID)); I
C <- length(unique(IPV$iso3)); C
R <- length(unique(IPV$region)); R
SR <- length(unique(IPV$SuperRegion)); SR
sum(IPV$denom_imp)
nrow(IPV)


# Countries with 2 estimates
study_cnt <- data.frame(iso3 = unique(IPV$iso3))
for (i in 1:nrow(study_cnt)) {
  cnt.i <- IPV[IPV$iso3 == study_cnt$iso3[i], ] 
  study_cnt$count[i] <- length(unique(cnt.i$ID))
}
sum(ifelse(study_cnt$count == 1, 1, 0), na.rm = TRUE)
sum(ifelse(study_cnt$count == 2, 1, 0), na.rm = TRUE)
sum(ifelse(study_cnt$count == 3, 1, 0), na.rm = TRUE)
sum(ifelse(study_cnt$count >= 4, 1, 0), na.rm = TRUE)

round(sum(ifelse(study_cnt$count == 1, 1, 0), na.rm = TRUE) / nrow(study_cnt) * 100, 0)
round(sum(ifelse(study_cnt$count == 2, 1, 0), na.rm = TRUE) / nrow(study_cnt) * 100, 0)
round(sum(ifelse(study_cnt$count == 3, 1, 0), na.rm = TRUE) / nrow(study_cnt) * 100, 0)
round(sum(ifelse(study_cnt$count >= 4, 1, 0), na.rm = TRUE) / nrow(study_cnt) * 100, 0)

# number of observations in each age category (by WHO region)
age_dist <- age_dist_who(IPV)
write.csv(age_dist, file = "VAW age distribution - Past Year IPV.csv")

# ---- Map of data availability ----
world_map_data(study_id = IPV$ID, iso3 = IPV$iso3, 
               outcome = "Past Yr IPV", suffix = "2000-2018")

# ---- Data formatting for JAGS ----
IPV$sup <- IPV$SuperRegion
IPV$reg <- IPV$region
IPV$cnt <- IPV$iso3
IPV$Super <- recode.cluster(IPV$SuperRegion)
IPV$Region <- recode.cluster(IPV$region)
IPV$Country <- recode.cluster(IPV$iso3)
IPV$Study <- recode.cluster(as.factor(IPV$ID))
SRLookUp_Cnt <- unique(IPV[c("Country", "Super")])[, "Super"]
SRLookUp <- unique(IPV[c("Region", "Super")])[, "Super"]
RLookUp <- unique(IPV[c("Country", "Region")])[, "Region"]
CLookUp <- unique(IPV[c("Study", "Country")])[, "Country"]
IPV$National <- ifelse(IPV$geo == "National", 1, 2)
NatLookUp <- unique(IPV[c("Study", "National")])[, "National"]

library(splines)
Center_a <- 30
Knots_a <- c(20) - Center_a
BoundKnots_a <- c(15, 50) - Center_a
Spl_a <- ns(x = (IPV$Age65 - Center_a), knots = Knots_a, Boundary.knots = BoundKnots_a)
Spl_W <- ns(seq(17.5, 67.5, by = 5) - Center_a, knots = Knots_a, Boundary.knots = BoundKnots_a)
Ndof <- dim(Spl_a)[2]

IPV$LI <- indexing_age_li(IPV$loage)
IPV$UI <- indexing_age_ui(IPV$hiageI)
IPV$UI <- ifelse(IPV$UI < IPV$LI, IPV$LI, IPV$UI)
IPV$ColIndex <- indexing_region(IPV$region)

Center_t <- 2018
Knots_t <- c(2011) - Center_t
BoundKnots_t <- c(2000, 2018) - Center_t
Spl_t <- ns(IPV$Time - Center_t, knots = Knots_t, Boundary.knots = BoundKnots_t)

time_dist <- aggregate(IPV$Time, by = list(IPV$ID), FUN = mean)
quantile(time_dist$x, probs = c(0.1, 0.9))

# ---- We Add Fixed Effects ----
load('LogOR-PastYr 2020-10-30.RData')

# Watchout the minus sign in front to get the correct reference category
XDat <- data.frame(Severe = ifelse(IPV$severe == "Only severe violence",
                      OR_to_use(IPV$SuperRegion, LogOR_Severe), 0),
                    SexVioOnly = ifelse(IPV$violence == "Sexual IPV only",
                      OR_to_use(IPV$SuperRegion, LogOR_Sexual), 0),
                    PhyVioOnly = ifelse(IPV$violence == "Physical IPV only",
                      OR_to_use(IPV$SuperRegion, LogOR_Physical), 0),
                    CurPartner = ifelse(IPV$pstat == "Currently partnered only",
                      OR_to_use(IPV$SuperRegion, LogOR_CurP), 0),
                    AllWom = ifelse(IPV$pstat == "All women",
                      OR_to_use(IPV$SuperRegion, LogOR_AllW), 0),
                    CMRpartner = ifelse(IPV$currpart == "Asked about violence from current or most recent partner only", 
                      OR_to_use(IPV$SuperRegion, LogOR_CMRpartner), 0),
                    # PastTwo = ifelse(IPV$viotime == "Past yr", 0, LogOR_Past2Yr$LogOR[LogOR_Past2Yr$SuperR == "Overall"]),
                    Rural = ifelse(IPV$geo == "Rural",
                      OR_to_use(IPV$SuperRegion, LogOR_Rural), 0),
                    Urban = ifelse(IPV$geo == "Urban",
                      OR_to_use(IPV$SuperRegion, LogOR_Urban), 0))

# Summary of adjustment factors
summary(XDat); dim(XDat)
# Observation level
apply(XDat, MARGIN = 2, FUN = function(x) sum(ifelse(x != 0, 1, 0), na.rm = TRUE))
sum(ifelse(rowSums(XDat) == 0, 1, 0), na.rm = TRUE)
not_optim <- sum(ifelse(rowSums(XDat) != 0, 1, 0), na.rm = TRUE) / nrow(XDat)
print(paste(round(not_optim * 100, 0), "% are not in optimal set"))

# Study level
xdat_id <- xdat_time <- xdat_geo <- XDat
xdat_id$ID <- IPV$ID
xdat_id$Age <- IPV$agerange
xdat_id$X <- rowSums(XDat)
xdat_time$ID <- IPV$ID
xdat_time$Time <- IPV$Time
xdat_geo$ID <- IPV$ID
xdat_geo$geo <- IPV$geo
X <- aggregate(xdat_id$X, by = list(xdat_id$ID), FUN = sum)
round(sum(ifelse(X$x != 0, 1, 0)) / nrow(X) * 100, 0)

UID <- X$Group.1
xstu <- xstu_t <- xstu_g <- NULL
for (i in 1:nrow(X)) {
  dat_i <- xdat_id[xdat_id$ID == UID[i], names(xdat_id) %in% names(XDat)]
  dat_t <- xdat_time[xdat_time$ID == UID[i], ]
  dat_g <- xdat_geo[xdat_geo$ID == UID[i], ]
  if (!all(apply(dat_i, 2, function(x) length(unique(x)) == 1)[!(names(dat_i) %in% c("Rural", "Urban"))])) {
    print(paste("Study with different covariates:", "(", i, ")", xdat_id$ID[xdat_id$ID == UID[i]][1]))
  }
  xstu <- rbind(xstu, dat_i[1, ])
  xstu_t <- rbind(xstu_t, dat_t[1, ])
  xstu_g <- c(xstu_g, any(dat_g$geo == "National"))
}
apply(xstu, MARGIN = 2, FUN = function(x) sum(ifelse(x != 0, 1, 0), na.rm = TRUE))
round(apply(xstu, MARGIN = 2, FUN = function(x) sum(ifelse(x != 0, 1, 0), na.rm = TRUE)) / nrow(xstu) * 100, 1)
not_optim <- sum(ifelse(rowSums(xstu) != 0, 1, 0), na.rm = TRUE) / nrow(xstu)
print(paste(round((1 - not_optim) * 100), "% of studies are in optimal set"))

summary(xstu_t$Time)
median(xstu_t$Time)
sum(ifelse(xstu_t$Time < 2005, 1, 0))
sum(ifelse(xstu_t$Time >= 2005 & xstu_t$Time < 2010, 1, 0))
sum(ifelse(xstu_t$Time >= 2010 & xstu_t$Time < 2015, 1, 0))
sum(ifelse(xstu_t$Time >= 2015, 1, 0))

round(sum(ifelse(xstu_t$Time < 2005, 1, 0)) / nrow(xstu_t) * 100, 0)
round(sum(ifelse(xstu_t$Time >= 2005 &  xstu_t$Time < 2010, 1, 0)) / nrow(xstu_t) * 100, 0)
round(sum(ifelse(xstu_t$Time >= 2010 & xstu_t$Time < 2015, 1, 0)) / nrow(xstu_t) * 100, 0)
round(sum(ifelse(xstu_t$Time >= 2015, 1, 0)) / nrow(xstu_t) * 100, 0)

# country-years
length(unique(paste(IPV$iso3, round(IPV$Time))))
# nationally representative
sum(ifelse(xstu_g, 1, 0))
round(mean(xstu_g) * 100, 0)

SDat_original <- SDat <- list(N = nrow(IPV), SR = SR, R = R, C = C, I = I, Ndof = Ndof,
             Num = floor(IPV$num_ess), Denom = floor(IPV$den_ess), 
             Nt = ncol(Spl_t),
             Spl_W = Spl_W, W = wom_by_age_and_region_2010, li = IPV$LI, ui = IPV$UI, ColIndex = IPV$ColIndex,
             Spl_t = Spl_t,
             XDat = rowSums(XDat),
             SRLookUp = SRLookUp, RLookUp = RLookUp, CLookUp = CLookUp, NatLookUp = NatLookUp,
             Country = IPV$Country,
             Study = IPV$Study)

sdat_past <- SDat_original
sdat_past$sup <- IPV$SuperRegion
sdat_past$cnt <- IPV$iso3
sdat_past$reg <- IPV$region
sdat_past$nat <- IPV$National
sdat_past$IPV <- IPV
saveRDS(sdat_past, file = "sdat_past.rds")
