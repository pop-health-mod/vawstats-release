devtools::load_all(here::here())

# loading database
# data not included... not yet publicly released
data(gbd)
sdat_ever <- readRDS("sdat_ever.rds")
sdat_past <- readRDS("sdat_past.rds")

#' Prediction for joint modeling
#' We need to create a list of all countries in any of the 2 datasets.
#' and then add one "typical" country by region (for prediction).
Center_a <- 30
uic_ <- unique(c(as.character(sdat_ever$IPV$iso3), as.character(sdat_past$IPV$iso3)))
age <- c(17.5, 22.5, 27.5, 42.5) 
la <- length(age)
liprd <- indexing_age_li(age - 2.5)
uiprd <- indexing_age_ui(age + 1.5)

Center_t <- 2018
time <- c(2000L, 2006L, 2012L, 2018L)
n_per <- length(time)

# Adding one "typical country" per region    
xprd <- NULL
uir <- as.character(unique(gbd$Region))
for (i in 1:length(uir)) {
  xprd.i <- data.frame(num = NA,
                      denom = 1,
                      cnt = rep(paste("CNT - ", uir[i]), n_per * la), 
                      reg = rep(uir[i], n_per * la),
                      sup = rep(gbd$SuperRegion[gbd$Region == uir[i]][1], n_per * la),
                      age = rep(age, n_per),
                      study = rep(paste("New", uir[i]), n_per * la),
                      li = rep(liprd, n_per),
                      ui = rep(uiprd, n_per),
                      ColIndex = rep(indexing_region(gbd$Region[gbd$Region == uir[i]][1]), n_per * la),
                      nat = 1,
                      time = sort(rep(time, la)))
  xprd <- rbind(xprd, xprd.i)  
}

unique(sdat_ever$IPV$iso3[which(!(as.character(sdat_ever$IPV$iso3) %in% as.character(sdat_past$IPV$iso3)))])
unique(sdat_past$IPV$iso3[which(!(as.character(sdat_past$IPV$iso3) %in% as.character(sdat_ever$IPV$iso3)))])

country_to_add <- c("ISR", "MAR", "SDN", "CAN", "SWZ", "NER", # countries with past year available but not ever
                    "MYS", "IRQ")                      # countries with ever available but not past year
country_to_add <- c(
  as.character(unique(sdat_ever$IPV$iso3[which(!(as.character(sdat_ever$IPV$iso3) %in% as.character(sdat_past$IPV$iso3)))])),
  as.character(unique(sdat_past$IPV$iso3[which(!(as.character(sdat_past$IPV$iso3) %in% as.character(sdat_ever$IPV$iso3)))])))
country_to_add

xprd_new <- NULL
for (i in 1:length(country_to_add)) {
xprd_new_i <- data.frame(num = rep(NA, n_per * la),
                      denom = rep(1, n_per * la),
                      cnt = rep(country_to_add[i], n_per * la), 
                      reg = rep(gbd$Region[gbd$iso == country_to_add[i]][1], n_per * la),
                      sup = rep(gbd$SuperRegion[gbd$iso == country_to_add[i]][1], n_per * la),
                      age = rep(age, n_per),
                      study = rep(paste("New", country_to_add[i]), n_per * la),
                      li = rep(liprd, n_per),
                      ui = rep(uiprd, n_per),
                      ColIndex = rep(indexing_region(gbd$Region[gbd$iso == country_to_add[i]][1]), n_per * la),
                      nat = 1,
                      time = sort(rep(time, la)))
xprd_new <- rbind(xprd_new, xprd_new_i)
}
xprd <- rbind(xprd, xprd_new)

xprd$study <- paste(xprd$study, xprd$time)
xprd <- xprd[order(xprd$age), ]
ind_1519 <- max(which(xprd$age == 17.5))


# loading Ever IPV data
ipv_e <- data.frame(num = sdat_ever$Num,
                    denom = sdat_ever$Denom,
                    cnt = sdat_ever$cnt,
                    reg = sdat_ever$reg,
                    sup = sdat_ever$sup,
                    study = sdat_ever$Study,
                    li = sdat_ever$li,
                    ui = sdat_ever$ui,
                    ColIndex = sdat_ever$ColIndex,
                    nat = sdat_ever$nat,
                    time = sdat_ever$IPV$Time)  
  
  ipv_e <- rbind(ipv_e, xprd[, !(names(xprd) %in% c("age"))])

  N_e <- nrow(ipv_e)
  SR_e <- length(unique(ipv_e$sup))
  R_e <- length(unique(ipv_e$reg))
  C_e <- length(unique(ipv_e$cnt))
  I_e <- length(unique(ipv_e$study))
  Num_e <- ipv_e$num
  Denom_e <- ipv_e$denom
  li_e <- ipv_e$li
  ui_e <- ipv_e$ui
  ColIndex_e <- ipv_e$ColIndex
  Spl_W_e <- sdat_ever$Spl_W
  W_e <- sdat_ever$W
  Ndof_e <- sdat_ever$Ndof
  XDat_e <- c(sdat_ever$XDat, rep(0, nrow(xprd)))
  Spl_t_e <- splines::ns(ipv_e$time - Center_t,
                       knots = attr(sdat_ever$Spl_t, "knots"),
                       Boundary.knots = attr(sdat_ever$Spl_t, "Boundary.knots"))
  Nt_e <- ncol(Spl_t_e)
    
  ipv_e$Super <- recode.cluster(ipv_e$sup)
  ipv_e$Region <- recode.cluster(ipv_e$reg)
  ipv_e$Country <- recode.cluster(ipv_e$cnt)
  ipv_e$Study <- recode.cluster(ipv_e$study)

  SRLookUp_e <- unique(ipv_e[c("Region", "Super")])[, "Super"] 
  RLookUp_e <- unique(ipv_e[c("Country", "Region")])[, "Region"]
  CLookUp_e <- unique(ipv_e[c("Study", "Country")])[, "Country"]
  NatLookUp_e <- unique(ipv_e[c("Study", "nat")])[, "nat"]


# Past Year IPV data
ipv_p <- data.frame(num = sdat_past$Num,
                    denom = sdat_past$Denom,
                    cnt = sdat_past$cnt,
                    reg = sdat_past$reg,
                    sup = sdat_past$sup,
                    study = sdat_past$Study,
                    li = sdat_past$li,
                    ui = sdat_past$ui,
                    ColIndex = sdat_past$ColIndex,
                    nat = sdat_past$nat,
                    time = sdat_past$IPV$Time)  
  
  ipv_p <- rbind(ipv_p, xprd[, !(names(xprd) %in% c("age"))])

  N_p <- nrow(ipv_p)
  SR_p <- length(unique(ipv_p$sup)) 
  R_p <- length(unique(ipv_p$reg))
  C_p <- length(unique(ipv_p$cnt))
  I_p <- length(unique(ipv_p$study))
  Num_p <- ipv_p$num
  Denom_p <- ipv_p$denom
  li_p <- ipv_p$li
  ui_p <- ipv_p$ui
  ColIndex_p <- ipv_p$ColIndex
  Spl_W_p <- sdat_past$Spl_W
  W_p <- sdat_past$W
  Ndof_p <- sdat_past$Ndof
  Nt_p <- sdat_past$Nt
  XDat_p <- c(sdat_past$XDat, rep(0, nrow(xprd)))
  Spl_t_p <- splines::ns(ipv_p$time - Center_t,
                       knots = attr(sdat_past$Spl_t, "knots"),
                       Boundary.knots = attr(sdat_past$Spl_t, "Boundary.knots"))
  Nt_p <- ncol(Spl_t_p)
    
  ipv_p$Super <- recode.cluster(ipv_p$sup)  
  ipv_p$Region <- recode.cluster(ipv_p$reg)
  ipv_p$Country <- recode.cluster(ipv_p$cnt)
  ipv_p$Study <- recode.cluster(ipv_p$study)

  SRLookUp_p <- unique(ipv_p[c("Region", "Super")])[, "Super"]
  RLookUp_p <- unique(ipv_p[c("Country", "Region")])[, "Region"]
  CLookUp_p <- unique(ipv_p[c("Study", "Country")])[, "Country"]
  NatLookUp_p <- unique(ipv_p[c("Study", "nat")])[, "nat"]
  

# Putting all in a list
SDat <- SDat_original <- list(W = wom_by_age_and_region_2010, 
             # Ever
             N_e = N_e, SR_e = SR_e, R_e = R_e, C_e = C_e, I_e = I_e,
             Ndof_e = Ndof_e, Nt_e = Nt_e,
             Num_e = Num_e, Denom_e = Denom_e,
             Spl_W_e = Spl_W_e, li_e = li_e, ui_e = ui_e, ColIndex_e = ColIndex_e,
             XDat_e = XDat_e,
             SRLookUp_e = SRLookUp_e, RLookUp_e = RLookUp_e, CLookUp_e = CLookUp_e, NatLookUp_e = NatLookUp_e,
             Study_e = ipv_e$Study,
             Country_e = ipv_e$Country,
             Spl_t_e = Spl_t_e,
             # Past year
             N_p = N_p, SR_p = SR_p, R_p = R_p, C_p = C_p, I_p = I_p, 
             Ndof_p = Ndof_p, Nt_p = Nt_p,
             Num_p = Num_p, Denom_p = Denom_p,
             Spl_W_p = Spl_W_p, li_p = li_p, ui_p = ui_p, ColIndex_p = ColIndex_p,
             XDat_p = XDat_p,
             SRLookUp_p = SRLookUp_p, RLookUp_p = RLookUp_p, CLookUp_p = CLookUp_p, NatLookUp_p = NatLookUp_p,
             Study_p = ipv_p$Study,
             Country_p = ipv_p$Country,
             Spl_t_p = Spl_t_p,
             # for "past <= ever" constraint
             N_const = nrow(xprd), start_e = sdat_ever$N, start_p = sdat_past$N,
             N_age = ind_1519,
             ones = rep(1, nrow(xprd)),
             ones_age = rep(1, ind_1519),
             n_age_e = nrow(Spl_W_e),
             n_age_p = nrow(Spl_W_p))


# ### ### ### ### ### ### ### ### ###
# ---- We Fit Here ----
# ### ### ### ### ### ### ### ### ###
library(rjags)
modelstring = "
model {
  # --- Ever IPV ----

  # Likelihood
  for (i in 1:N_e) {
    Num_e[i] ~ dbin(linpred_e[i], Denom_e[i])
    # Splines times
    Spl_t_e_coeff[i] <- inprod(Spl_t_e[i, ], b_t_c_e[Country_e[i], ])
    # Age-Standardization
    total_weight_e[i] <- sum(W[li_e[i]:ui_e[i], ColIndex_e[i]])
      for (j in 1:n_age_e) {
        logit(age_prev_e[i, j]) <- a_i_e[Study_e[i]] + XDat_e[i] + Spl_t_e_coeff[i] + age_spl_country_e[Country_e[i], j]
        logit(prv_age_constraint_e[i, j]) <- a_c_e[Country_e[i]] + Spl_t_e_coeff[i] + age_spl_country_e[Country_e[i], j]
      }
      linpred_e[i] <- inprod(age_prev_e[i, li_e[i]:ui_e[i]], 
                             W[li_e[i]:ui_e[i], ColIndex_e[i]] / total_weight_e[i])
      # For constraint
      prv_e[i] <- inprod(prv_age_constraint_e[i, li_e[i]:ui_e[i]], 
                         W[li_e[i]:ui_e[i], ColIndex_e[i]] / total_weight_e[i])
    }

  # Priors
    a_g_e ~ dnorm(0, 0.001)
  # Time Trends
    for (t in 1:Nt_e) { b_t_e[t] ~ dnorm(0, 0.005)  }
  # Age Spline
    for (j in 1:Ndof_e) { b_spl_e[j] ~ dnorm(0, 0.005) }

  # Random Intercepts & Slope
  # Study
  for (i in 1:I_e) {
    a_i_e[i] ~ dnorm(a_c_e[CLookUp_e[i]], tau_ai_e[NatLookUp_e[i]])
  }
  # Country
  for (c in 1:C_e) {
    a_c_e[c] ~ dnorm(a_r_e[RLookUp_e[c]], tau_ac_e)
    # indexing to get matrix of contry-specific and age-specific spline coefficients
    for (k in 1:n_age_e) {
      age_spl_country_e[c, k] <- inprod(Spl_W_e[k, ],  b_spl_c_e[c, ]) }
    for (t in 1:Nt_e) {
      b_t_c_e[c, t] ~ dnorm(b_t_r_e[RLookUp_e[c], t], tau_t_c_e[t]) }
    for (j in 1:Ndof_e) {
      b_spl_c_e[c, j] ~ dnorm(b_spl_r_e[RLookUp_e[c], j], tau_spl_c_e[j]) }
    }
  # Region
  for (r in 1:R_e) {
    a_r_e[r] ~ dnorm(a_sr_e[SRLookUp_e[r]], tau_ar_e)
    for (t in 1:Nt_e) {
      b_t_r_e[r, t] ~ dnorm(b_t_sr_e[SRLookUp_e[r], t], tau_t_r_e[t]) }
    for (j in 1:Ndof_e) {
      b_spl_r_e[r, j] ~ dnorm(b_spl_sr_e[SRLookUp_e[r], j], tau_spl_r_e[j]) }
  }
  # Super Region
  for (s in 1:SR_e) {
    a_sr_e[s] ~ dnorm(a_g_e, tau_sr_e)
    for (t in 1:Nt_e) {
      b_t_sr_e[s, t] ~ dnorm(b_t_e[t], tau_t_sr_e[t]) }
    for (j in 1:Ndof_e) {
      b_spl_sr_e[s, j] ~ dnorm(b_spl_e[j], tau_spl_sr_e[j]) }
  }

  # Priors for precision of Random Effects
  # Variance study
    tau_ai_e[1] <- pow(sd_ai_e[1], -2)
    tau_ai_e[2] <- pow(sd_ai_e[2], -2)
      sd_ai_e[1] ~ dt(0, 1/25, 1)I(0, )
      sd_ai_e[2] <- sd_ai_e[1] + sd_ai_NotNat_e
      sd_ai_NotNat_e ~ dt(0, 1/25, 1)I(0, )
  # Variance country
    tau_ac_e <- pow(sd_ac_e, -2)
      sd_ac_e ~ dt(0, 1/25, 1)I(0, )
  # Variance region
    tau_ar_e <- pow(sd_ar_e, -2)
      sd_ar_e ~ dt(0, 1/25, 1)I(0, )
  # Variance super region
    tau_sr_e <- pow(sd_sr_e, -2)
      sd_sr_e ~ dt(0, 1/25, 1)I(0, )
  # Variance age spline at region and super regions
  for (j in 1:Ndof_e) {
    tau_spl_sr_e[j] <- pow(sd_spl_sr_e[j], -2)
      sd_spl_sr_e[j] ~ dt(0, 1/25, 1)I(0, )
    tau_spl_r_e[j] <- pow(sd_spl_r_e[j], -2)
      sd_spl_r_e[j] ~ dt(0, 1/25, 1)I(0, )
    tau_spl_c_e[j] <- pow(sd_spl_c_e[j], -2)
      sd_spl_c_e[j] ~ dt(0, 1/25, 1)I(0, )
  }
  # Variance time trend at country and super region
  for (t in 1:Nt_e) {
    tau_t_sr_e[t] <- pow(sd_t_sr_e[t], -2)
      sd_t_sr_e[t] ~ dt(0, 1/25, 1)I(0, )
    tau_t_r_e[t] <- pow(sd_t_r_e[t], -2)
      sd_t_r_e[t] ~ dt(0, 1/25, 1)I(0, )
    tau_t_c_e[t] <- pow(sd_t_c_e[t], -2)
      sd_t_c_e[t] ~ dt(0, 1/25, 1)I(0, )
  }

  # --- Past Year IPV ----
  # Likelihood
  for (i in 1:N_p) {
    Num_p[i] ~ dbin(linpred_p[i], Denom_p[i])
    # Splines times
    Spl_t_p_coeff[i] <- inprod(Spl_t_p[i, ], b_t_c_p[Country_p[i], ])
    # Age-Standardization
    total_weight_p[i] <- sum(W[li_p[i]:ui_p[i], ColIndex_p[i]])
      for (j in 1:n_age_p) {
        logit(age_prev_p[i, j]) <- a_i_p[Study_p[i]] + XDat_p[i] + Spl_t_p_coeff[i] + age_spl_country_p[Country_p[i], j]
        logit(prv_age_constraint_p[i, j]) <- a_c_p[Country_p[i]] + Spl_t_p_coeff[i] + age_spl_country_p[Country_p[i], j]
      }
      linpred_p[i] <- inprod(age_prev_p[i, li_p[i]:ui_p[i]], 
                           W[li_p[i]:ui_p[i], ColIndex_p[i]] / total_weight_p[i])
      # For constraint
      prv_p[i] <- inprod(prv_age_constraint_p[i, li_p[i]:ui_p[i]], 
                           W[li_p[i]:ui_p[i], ColIndex_p[i]] / total_weight_p[i])
  }

  # Priors
    a_g_p ~ dnorm(0, 0.001)
  # Time Trends
    for (t in 1:Nt_p) { b_t_p[t] ~ dnorm(0, 0.001)  }
  # Age Spline
    for (j in 1:Ndof_p) { b_spl_p[j] ~ dnorm(0, 0.001) }

  # Random Intercepts & Slope
  # Study
  for (i in 1:I_p) {
    a_i_p[i] ~ dnorm(a_c_p[CLookUp_p[i]], tau_ai_p[NatLookUp_p[i]])
  }
  # Country
  for (c in 1:C_p) {
    a_c_p[c] ~ dnorm(a_r_p[RLookUp_p[c]], tau_ac_p)
    # indexing to get matrix of contry-specific and age-specific spline coefficients
    for (k in 1:n_age_e) {
      age_spl_country_p[c, k] <- inprod(Spl_W_p[k, ],  b_spl_c_p[c, ]) }
    for (t in 1:Nt_p) {
      b_t_c_p[c, t] ~ dnorm(b_t_r_p[RLookUp_p[c], t], tau_t_c_p[t]) }
    for (j in 1:Ndof_p) {
      b_spl_c_p[c, j] ~ dnorm(b_spl_r_p[RLookUp_p[c], j], tau_spl_c_p[j]) }
    }
  # Region
  for (r in 1:R_p) {
    a_r_p[r] ~ dnorm(a_sr_p[SRLookUp_p[r]], tau_ar_p)
    for (t in 1:Nt_p) {
      b_t_r_p[r, t] ~ dnorm(b_t_sr_p[SRLookUp_p[r], t], tau_t_r_p[t]) }
    for (j in 1:Ndof_p) {
      b_spl_r_p[r, j] ~ dnorm(b_spl_sr_p[SRLookUp_p[r], j], tau_spl_r_p[j]) }
  }
  # Super Region
  for (s in 1:SR_p) {
    a_sr_p[s] ~ dnorm(a_g_p, tau_sr_p)
    for (t in 1:Nt_p) {
      b_t_sr_p[s, t] ~ dnorm(b_t_p[t], tau_t_sr_p[t]) }
    for (j in 1:Ndof_p) {
      b_spl_sr_p[s, j] ~ dnorm(b_spl_p[j], tau_spl_sr_p[j]) }   
  }

  # Priors for precision of Random Effects
  # Variance study
    tau_ai_p[1] <- pow(sd_ai_p[1], -2)
    tau_ai_p[2] <- pow(sd_ai_p[2], -2)
      sd_ai_p[1] ~ dt(0, 1/25, 1)I(0, )
      sd_ai_p[2] <- sd_ai_p[1] + sd_ai_NotNat_p
      sd_ai_NotNat_p ~ dt(0, 1/25, 1)I(0, )
  # Variance country  
    tau_ac_p <- pow(sd_ac_p, -2)
      sd_ac_p ~ dt(0, 1/25, 1)I(0, )
  # Variance region
    tau_ar_p <- pow(sd_ar_p, -2)
      sd_ar_p ~ dt(0, 1/25, 1)I(0, )
  # Variance super region
    tau_sr_p <- pow(sd_sr_p, -2)
      sd_sr_p ~ dt(0, 1/25, 1)I(0, )
  # Variance age spline at region and super regions
  for (j in 1:Ndof_p) {
    tau_spl_sr_p[j] <- pow(sd_spl_sr_p[j], -2)
      sd_spl_sr_p[j] ~ dt(0, 1/25, 1)I(0, )
    tau_spl_r_p[j] <- pow(sd_spl_r_p[j], -2)
      sd_spl_r_p[j] ~ dt(0, 1/25, 1)I(0, )
    tau_spl_c_p[j] <- pow(sd_spl_c_p[j], -2)
      sd_spl_c_p[j] ~ dt(0, 1/25, 1)I(0, )
  }
  # Variance time trend at country and super region
  for (t in 1:Nt_p) {
    tau_t_sr_p[t] <- pow(sd_t_sr_p[t], -2)
      sd_t_sr_p[t] ~ dt(0, 1/25, 1)I(0, )
    tau_t_r_p[t] <- pow(sd_t_r_p[t], -2)
      sd_t_r_p[t] ~ dt(0, 1/25, 1)I(0, )
    tau_t_c_p[t] <- pow(sd_t_c_p[t], -2)
      sd_t_c_p[t] ~ dt(0, 1/25, 1)I(0, )
  }
  
  # --- The constraint to have past year <= ever prevalence ---- 
  for (z in 1:N_const) {
    ones[z] ~ dbern(const[z])
    const[z] <- step(prv_e[start_e + z] - prv_p[start_p + z])
  }
  # --- The constraint to have the ratio of lifetime / past year IPV less than 3 for the 15-19 years old ---
  for (a in 1:N_age) {
    ones_age[a] ~ dbern(const_age[a])
    const_age[a] <- step(3 - rr_age[a])
    rr_age[a] <- prv_e[start_e + a] / prv_p[start_p + a]
  }
}"

# We sample the model
load.module("glm")
n_chain <- 4
n_adapt <- 10000
n_burn <- 5000
n_iter <- 50000
n_thin <- 20

param <- c(# Ever
          'a_g_e', 'a_i_e', 'a_c_e', 'a_r_e', 'a_sr_e',
          'b_spl_e', 'b_spl_sr_e', 'b_spl_r_e', 'b_spl_c_e',      
          'b_t_e', 'b_t_sr_e', 'b_t_r_e', 'b_t_c_e',
          'sd_ai_e', 'sd_ac_e', 'sd_ar_e', 'sd_sr_e', 'sd_ai_NotNat_e',
          'sd_spl_sr_e', 'sd_spl_r_e', 'sd_spl_c_e',
          'sd_t_sr_e', 'sd_t_r_e', 'sd_t_c_e',
          'linpred_e',
          # Past year
          'a_g_p', 'a_i_p', 'a_c_p', 'a_r_p', 'a_sr_p',
          'b_spl_p', 'b_spl_sr_p', 'b_spl_r_p', 'b_spl_c_p',       
          'b_t_p', 'b_t_sr_p', 'b_t_r_p', 'b_t_c_p',
          'sd_ai_p', 'sd_ac_p', 'sd_ar_p', 'sd_sr_p', 'sd_ai_NotNat_p',
          'sd_spl_sr_p', 'sd_spl_r_p', 'sd_spl_c_p',
          'sd_t_sr_p', 'sd_t_r_p', 'sd_t_c_p',
          'linpred_p',
          'rr_age')


# ---- We Add Fixed Effects ----
load('LogOR-Ever 2020-10-30.RData')
LogOR_Severe_e <- LogOR_Severe
LogOR_Sexual_e <- LogOR_Sexual
LogOR_Physical_e <- LogOR_Physical
LogOR_CurP_e <- LogOR_CurP
LogOR_AllW_e <- LogOR_AllW
LogOR_CMRpartner_e <- LogOR_CMRpartner
LogOR_Rural_e <- LogOR_Rural
LogOR_Urban_e <- LogOR_Urban

load('LogOR-PastYr 2020-10-30.RData')
LogOR_Severe_p <- LogOR_Severe
LogOR_Sexual_p <- LogOR_Sexual
LogOR_Physical_p <- LogOR_Physical
LogOR_CurP_p <- LogOR_CurP
LogOR_AllW_p <- LogOR_AllW
LogOR_CMRpartner_p <- LogOR_CMRpartner
LogOR_Rural_p <- LogOR_Rural
LogOR_Urban_p <- LogOR_Urban


# ---- LHS Sample for X-Walk ----
#' --- Set to 1 if you only want to fit without uncertainty in the x-walk
n_mi <- 10
#' ---
lhs_severe_e <- OR_to_use_lhs(LogOR_Severe_e, n = n_mi, truncate = TRUE)
lhs_sexual_e <- OR_to_use_lhs(LogOR_Sexual_e, n = n_mi, truncate = TRUE)
lhs_physical_e <- OR_to_use_lhs(LogOR_Physical_e, n = n_mi, truncate = TRUE)
lhs_curpart_e <- OR_to_use_lhs(LogOR_CurP_e, n = n_mi, truncate = TRUE)
lhs_allwom_e <- OR_to_use_lhs(LogOR_AllW_e, n = n_mi, truncate = FALSE)
lhs_cmr_e <- OR_to_use_lhs(LogOR_CMRpartner_e, n = n_mi, truncate = TRUE)
lhs_rural_e <- OR_to_use_lhs(LogOR_Rural_e, n = n_mi, truncate = FALSE)
lhs_urban_e <- OR_to_use_lhs(LogOR_Urban_e, n = n_mi, truncate = FALSE)


lhs_severe_p <- OR_to_use_lhs(LogOR_Severe_p, n = n_mi, truncate = TRUE)
lhs_sexual_p <- OR_to_use_lhs(LogOR_Sexual_p, n = n_mi, truncate = TRUE)
lhs_physical_p <- OR_to_use_lhs(LogOR_Physical_p, n = n_mi, truncate = TRUE)
lhs_curpart_p <- OR_to_use_lhs(LogOR_CurP_p, n = n_mi, truncate = TRUE)
lhs_allwom_p <- OR_to_use_lhs(LogOR_AllW_p, n = n_mi, truncate = FALSE)
lhs_cmr_p <- OR_to_use_lhs(LogOR_CMRpartner_p, n = n_mi, truncate = TRUE)
lhs_rural_p <- OR_to_use_lhs(LogOR_Rural_p, n = n_mi, truncate = FALSE)
lhs_urban_p <- OR_to_use_lhs(LogOR_Urban_p, n = n_mi, truncate = FALSE)

# ---- Multiple Imputation Here ----
res_lhs_e <- list()
res_lhs_p <- list()
XDat_e <- list()
XDat_p <- list()
for (i in 1:n_mi) {
  XDat_e[[i]] <- data.frame(Severe = ifelse(sdat_ever$IPV$severe == "Only severe violence",
      OR_to_use(sdat_ever$IPV$SuperRegion, lhs_severe_e, lhs = TRUE, index = i), 0),
    SexVioOnly = ifelse(sdat_ever$IPV$violence == "Sexual IPV only",
      OR_to_use(sdat_ever$IPV$SuperRegion, lhs_sexual_e, lhs = TRUE, index = i), 0),
    PhyVioOnly = ifelse(sdat_ever$IPV$violence == "Physical IPV only",
      OR_to_use(sdat_ever$IPV$SuperRegion, lhs_physical_e, lhs = TRUE, index = i), 0),
    CurPartner = ifelse(sdat_ever$IPV$pstat == "Currently partnered only",
      OR_to_use(sdat_ever$IPV$SuperRegion, lhs_curpart_e, lhs = TRUE, index = i), 0),
    AllWom = ifelse(sdat_ever$IPV$pstat == "All women",
      OR_to_use(sdat_ever$IPV$SuperRegion, lhs_allwom_e, lhs = TRUE, index = i), 0),
    AnyHusband = ifelse(sdat_ever$IPV$currpart == "Asked about violence from current or most recent partner only",
      OR_to_use(sdat_ever$IPV$SuperRegion, lhs_cmr_e, lhs = TRUE, index = i), 0),
    Rural = ifelse(sdat_ever$IPV$geo == "Rural",
      OR_to_use(sdat_ever$IPV$SuperRegion, lhs_rural_e, lhs = TRUE, index = i), 0),
    Urban = ifelse(sdat_ever$IPV$geo == "Urban",
      OR_to_use(sdat_ever$IPV$SuperRegion, lhs_urban_e, lhs = TRUE, index = i), 0))

  XDat_p[[i]] <- data.frame(Severe = ifelse(sdat_past$IPV$severe == "Only severe violence",
      OR_to_use(sdat_past$IPV$SuperRegion, lhs_severe_p, lhs = TRUE, index = i), 0),
    SexVioOnly = ifelse(sdat_past$IPV$violence == "Sexual IPV only",
      OR_to_use(sdat_past$IPV$SuperRegion, lhs_sexual_p, lhs = TRUE, index = i), 0),
    PhyVioOnly = ifelse(sdat_past$IPV$violence == "Physical IPV only",
      OR_to_use(sdat_past$IPV$SuperRegion, lhs_physical_p, lhs = TRUE, index = i), 0),
    CurPartner = ifelse(sdat_past$IPV$pstat == "Currently partnered only",
      OR_to_use(sdat_past$IPV$SuperRegion, lhs_curpart_p, lhs = TRUE, index = i), 0),
    AllWom = ifelse(sdat_past$IPV$pstat == "All women",
      OR_to_use(sdat_past$IPV$SuperRegion, lhs_allwom_p, lhs = TRUE, index = i), 0),
    AnyHusband = ifelse(sdat_past$IPV$currpart == "Asked about violence from current or most recent partner only",
      OR_to_use(sdat_past$IPV$SuperRegion, lhs_cmr_p, lhs = TRUE, index = i), 0),
    Rural = ifelse(sdat_past$IPV$geo == "Rural",
      OR_to_use(sdat_past$IPV$SuperRegion, lhs_rural_p, lhs = TRUE, index = i), 0),
    Urban = ifelse(sdat_past$IPV$geo == "Urban",
      OR_to_use(sdat_past$IPV$SuperRegion, lhs_urban_p, lhs = TRUE, index = i), 0))
  
  SDat$XDat_e <- c(rowSums(XDat_e[[i]]), rep(0, SDat$N_e - SDat$start_e))
  SDat$XDat_p <- c(rowSums(XDat_p[[i]]), rep(0, SDat$N_p - SDat$start_p))
  
  if (n_mi == 1) {
  jags.fit <- jags.model(textConnection(modelstring), n.adapt = n_adapt, 
                         data = SDat_original, n.chains = n_chain)
  } else {
  jags.fit <- jags.model(textConnection(modelstring), n.adapt = n_adapt, 
                         data = SDat, n.chains = n_chain)
  }
                         
  # We update the model using 5,000 iterations as burn-in
  update(jags.fit, n_burn)
  
  # To extract random samples from the posterior distribution of the parameters of a jags model
  MCMC <- jags.samples(jags.fit, variable.names = param, n.iter = n_iter, thin = n_thin)
 
  res_lhs_e[[i]] <- process_mcmc(MCMC, report = "ever") 
  res_lhs_p[[i]] <- process_mcmc(MCMC, report = "past")
}

  res_ever <- process_lhs(res_lhs_e)
  res_past <- process_lhs(res_lhs_p)


# Diagnostics
DoDiagnostic <- FALSE
if (DoDiagnostic == TRUE) {
  # --- EVER 
  library(rjags)
  # using the last imputation - ever
  # multiply by n_imp to get the real number of samples
  MCMC_e <- process_mcmc(MCMC, report = "ever")
  effectiveSize(MCMC_e[["a_g"]])
  summary(effectiveSize(MCMC_e[["a_i"]]))
  summary(effectiveSize(MCMC_e[["a_c"]]))
  summary(effectiveSize(MCMC_e[["a_r"]]))
  summary(effectiveSize(MCMC_e[["a_sr"]]))
  
  effectiveSize(MCMC_e[["b_spl"]])
  summary(effectiveSize(MCMC_e[["b_spl_sr"]]))
  summary(effectiveSize(MCMC_e[["b_spl_r"]]))
  summary(effectiveSize(MCMC_e[["b_spl_c"]]))
  
  effectiveSize(MCMC_e[["b_t"]])
  summary(effectiveSize(MCMC_e[["b_t_sr"]]))
  summary(effectiveSize(MCMC_e[["b_t_r"]]))    
  summary(effectiveSize(MCMC_e[["b_t_c"]]))   

  effectiveSize(MCMC_e[["sd_ai"]])
  effectiveSize(MCMC_e[["sd_ac"]])
  effectiveSize(MCMC_e[["sd_ar"]])
  effectiveSize(MCMC_e[["sd_sr"]])
  effectiveSize(MCMC_e[["sd_spl_sr"]])  
  effectiveSize(MCMC_e[["sd_spl_r"]])  
  effectiveSize(MCMC_e[["sd_spl_c"]])  
  effectiveSize(MCMC_e[["sd_t_sr"]])  
  effectiveSize(MCMC_e[["sd_t_r"]])  
  effectiveSize(MCMC_e[["sd_t_c"]])
 
  gelman.diag(MCMC_e[["a_i"]])
  gelman.diag(MCMC_e[["a_c"]])    
  gelman.diag(MCMC_e[["a_r"]])  
  
  # using the last imputation - past
  MCMC_p <- process_mcmc(MCMC, report = "past")
  effectiveSize(MCMC_p[["a_g"]])
  summary(effectiveSize(MCMC_p[["a_i"]]))
  summary(effectiveSize(MCMC_p[["a_c"]]))
  summary(effectiveSize(MCMC_p[["a_r"]]))
  summary(effectiveSize(MCMC_p[["a_sr"]]))
  
  effectiveSize(MCMC_p[["b_spl"]])
  summary(effectiveSize(MCMC_p[["b_spl_sr"]]))
  summary(effectiveSize(MCMC_p[["b_spl_r"]]))
  summary(effectiveSize(MCMC_p[["b_spl_c"]]))
  
  effectiveSize(MCMC_p[["b_t"]])
  summary(effectiveSize(MCMC_p[["b_t_sr"]]))
  summary(effectiveSize(MCMC_p[["b_t_r"]]))    
  summary(effectiveSize(MCMC_p[["b_t_c"]]))   

  effectiveSize(MCMC_p[["sd_ai"]])
  effectiveSize(MCMC_p[["sd_ac"]])
  effectiveSize(MCMC_p[["sd_ar"]])
  effectiveSize(MCMC_p[["sd_sr"]])
  
  effectiveSize(MCMC_p[["sd_spl_sr"]])
  effectiveSize(MCMC_p[["sd_spl_r"]])  
  effectiveSize(MCMC_p[["sd_spl_c"]])  
  
  effectiveSize(MCMC_p[["sd_t_sr"]]) 
  effectiveSize(MCMC_p[["sd_t_r"]])  
  effectiveSize(MCMC_p[["sd_t_c"]])  
  
  gelman.diag(MCMC_p[["a_i"]])
  gelman.diag(MCMC_p[["a_c"]])    
  gelman.diag(MCMC_p[["a_r"]])  
  
  # plotting traceplots (takes a while to execute)
  Col <- c(rgb(20, 180, 235, 255, max = 255), rgb(134, 179, 0, 255, max = 255), rgb(204, 0, 67, 255, max = 255),
           rgb(255, 173, 51, 255, max = 255), rgb(255, 102, 53, 255, max = 255), rgb(135, 132, 148, 255, max = 255),
           rgb(128, 0, 128, 150, max = 255), rgb(102, 0, 53, 255, max = 255),
           rgb(89, 165, 165, max = 255), rgb(178, 102, 102, 255, max = 255))[1:n_chain]
  SDat <- SDat_original
  jpeg(file = "traceplots ever jnt.jpeg", width = 21, height = 345, units = "in", res = 100)
    par(mar = c(2, 1, 2, 1), oma = c(0, 1, 0, 0))
  n_parm <- 1 + SDat$SR_e + SDat$R_e + SDat$C_e + SDat$I_e + 
      SDat$Ndof_e + SDat$Ndof_e * (SDat$SR_e + SDat$R_e + SDat$C_e) +
      SDat$Nt_e + SDat$Nt_e * (SDat$SR_e + SDat$R_e + SDat$C_e) +    
      5  + SDat$Ndof_e * 3 + SDat$Nt_e * 3 + 1
    nf <- layout(matrix(c(1:(10 * ceiling(n_parm / 10))), ncol = 10, nrow = ceiling(n_parm / 10), byrow = TRUE),
            widths = rep(lcm(5), 10), heights = rep(lcm(5), ceiling(n_parm / 10)), respect = TRUE)
    traceplot(res_ever[["a_g"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["a_sr"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["a_r"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["a_c"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["a_i"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["b_spl"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["b_spl_sr"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["b_spl_r"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["b_spl_c"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["b_t"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["b_t_sr"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["b_t_r"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["b_t_c"]], smooth = TRUE, col = Col) 
    traceplot(res_ever[["sd_sr"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["sd_ar"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["sd_ac"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["sd_ai"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["sd_ai2"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["sd_spl_sr"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["sd_spl_r"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["sd_spl_c"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["sd_t_sr"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["sd_t_r"]], smooth = TRUE, col = Col)
    traceplot(res_ever[["sd_t_c"]], smooth = TRUE, col = Col)
  dev.off()
  
  jpeg(file = "traceplots past jnt.jpeg", width = 21, height = 340, units = "in", res = 150)
  par(mar = c(2, 1, 2, 1), oma = c(0, 1, 0, 0))
  n_parm <- 1 + SDat$SR_p + SDat$R_p + SDat$C_p + SDat$I_p + 
      SDat$Ndof_p + SDat$Ndof_p * (SDat$SR_p + SDat$R_p + SDat$C_p) +
      SDat$Nt_p + SDat$Nt_p * (SDat$SR_p + SDat$R_p + SDat$C_p) +  
      5  + SDat$Ndof_p * 3 + SDat$Nt_p * 3 + 1
    nf <- layout(matrix(c(1:(10 * ceiling(n_parm / 10))), ncol = 10, nrow = ceiling(n_parm / 10), byrow = TRUE),
            widths = rep(lcm(5), 10), heights = rep(lcm(5), ceiling(n_parm / 10)), respect = TRUE)
    traceplot(res_past[["a_g"]], smooth = TRUE, col = Col)
    traceplot(res_past[["a_sr"]], smooth = TRUE, col = Col)
    traceplot(res_past[["a_r"]], smooth = TRUE, col = Col)
    traceplot(res_past[["a_c"]], smooth = TRUE, col = Col)
    traceplot(res_past[["a_i"]], smooth = TRUE, col = Col)
    traceplot(res_past[["b_spl"]], smooth = TRUE, col = Col)
    traceplot(res_past[["b_spl_sr"]], smooth = TRUE, col = Col)
    traceplot(res_past[["b_spl_r"]], smooth = TRUE, col = Col)
    traceplot(res_past[["b_spl_c"]], smooth = TRUE, col = Col)
    traceplot(res_past[["b_t"]], smooth = TRUE, col = Col)
    traceplot(res_past[["b_t_sr"]], smooth = TRUE, col = Col)
    traceplot(res_past[["b_t_r"]], smooth = TRUE, col = Col)
    traceplot(res_past[["b_t_c"]], smooth = TRUE, col = Col) 
    traceplot(res_past[["sd_sr"]], smooth = TRUE, col = Col)
    traceplot(res_past[["sd_ar"]], smooth = TRUE, col = Col)
    traceplot(res_past[["sd_ac"]], smooth = TRUE, col = Col)
    traceplot(res_past[["sd_ai"]], smooth = TRUE, col = Col)
    traceplot(res_past[["sd_ai2"]], smooth = TRUE, col = Col)
    traceplot(res_past[["sd_spl_sr"]], smooth = TRUE, col = Col)
    traceplot(res_past[["sd_spl_r"]], smooth = TRUE, col = Col)
    traceplot(res_past[["sd_spl_c"]], smooth = TRUE, col = Col)
    traceplot(res_past[["sd_t_sr"]], smooth = TRUE, col = Col)
    traceplot(res_past[["sd_t_r"]], smooth = TRUE, col = Col)
    traceplot(res_past[["sd_t_c"]], smooth = TRUE, col = Col)
    dev.off()
}

#' ======================================
# ---- Graphs Ever ----
#' ======================================
post_pred_chk(res_ever, sdat_ever, IPV = sdat_ever$IPV, outcome = "Ever IPV", suffix = "2020-10-30_JNT") 

plot_age_pattern(res_ever, sdat_ever, IPV = ipv_e, outcome = "Ever IPV", x_max = 65,
                 Center_a = Center_a, Center_t = Center_t, suffix = "2020-10-30_JNT")

pred_ever <- pool_pred(res_ever, sdat_ever, IPV = ipv_e, Center_a = Center_a, Center_t = Center_t,
                       year_pred = 2018, outcome = "Ever IPV") 
agg_cnt_ever <- pool_cnt(pred_ever, IPV = sdat_ever$IPV, outcome = "Ever IPV", suffix = "2020-10-30_JNT") 
agg_reg_ever <- pool_reg(pred_ever, ll_age = 15, ul_age = 49, outcome = "Ever IPV", suffix = "2020-10-30_JNT"); agg_reg_ever 
agg_sdg_ever <- pool_reg_sdg(pred_ever, ll_age = 15, ul_age = 49, outcome = "Ever IPV", suffix = "2020-10-30_JNT"); agg_sdg_ever
agg_who_ever <- pool_reg_who(pred_ever, ll_age = 15, ul_age = 49, outcome = "Ever IPV", suffix = "2020-10-30_JNT"); agg_who_ever
agg_who_ever <- pool_reg_who(pred_ever, ll_age = 15, ul_age = 104, outcome = "Ever IPV", suffix = "2020-10-30_JNT"); agg_who_ever
agg_glb_ever <- pool_glb(pred_ever, ll_age = 15, ul_age = 49, outcome = "Ever IPV", suffix = "2020-10-30_JNT"); agg_glb_ever

agg_glb_ever_exc <- pool_glb(pred_ever, ll_age = 15, ul_age = 49, 
                             outcome = "Ever IPV", cnt_exclude = c("CHN"), suffix = "excl. CHN 2020-10-30_JNT", save_results = FALSE)

#' ---- New Checks - 15-49 by Cnt*Time ----
plot_time_trend(res_ever, sdat_ever, IPV = ipv_e, Center_t = Center_t,
                outcome = "Ever IPV", suffix = "2020-10-30_JNT")
agg_cnt_time_ever <- pool_pred_time(res_ever, sdat_ever, IPV = ipv_e, Center_a = Center_a, Center_t = Center_t,
                outcome = "Ever IPV", suffix = "2020-10-30_JNT", max_iter = 25000, age_grp = "15-49")
plot_agg_cnt_time(agg_cnt_time_ever, sdat_ever$IPV, outcome = "Ever IPV", suffix = "2020-10-30_JNT")


# ---- Ever - Heatmap and Map all countries  ----
agg_cnt_ever <- read.csv("Estimates Age by Cnt - Ever IPV 2020-10-30_JNT.csv")
agg_cnt_ever$cnt <- countrycode::countrycode(sourcevar = agg_cnt_ever$Country, origin = "iso3c", destination = "un.name.en")
agg_cnt_ever$cnt[agg_cnt_ever$cnt == "Lao People’s Democratic Republic"] <- "Lao People's Democratic Republic"

heatmap_vaw(agg_cnt_ever, cnt_var = "cnt", outcome = "Ever IPV", suffix = "2020-10-30_JNT", 
                        width = 8.5, height = 17)
world_map_prv(agg_cnt_ever, ipv = sdat_ever$IPV, age_grp = "15-49", outcome = "Ever IPV", width = 11, height = 8.5,
                           suffix = "2020-10-30_JNT")
prv_plot_cnt(agg_cnt_ever, sdat_ever$IPV, age_grp = "15-49", outcome = "Ever IPV", width = 8.5, height = 19,
                           suffix = "2020-10-30_JNT", xlim = c(0, 70))
  
agg_reg_ever <- read.csv("Estimates by GBD Region (15-49) - Ever IPV 2020-10-30_JNT.csv")
prv_plot_reg(agg_reg_ever, sdat_ever$IPV, age_grp = "15-49", outcome = "Ever IPV", width = 8.5, height = 6,
                           suffix = "2020-10-30_JNT", xlim = c(0, 70))

#' ======================================
# ---- Graphs Past ----
#' ======================================
post_pred_chk(res_past, sdat_past, IPV = sdat_past$IPV, outcome = "Past Year IPV", suffix = "2020-10-30_JNT")  

plot_age_pattern(res_past, sdat_past, IPV = sdat_past$IPV, outcome = "Past Year IPV", x_max = 65, 
                 Center_a = Center_a, Center_t = Center_t, suffix = "2020-10-30_JNT")

pred_past <- pool_pred(res_past, sdat_past, IPV = ipv_p, Center_a = Center_a, Center_t = Center_t, 
                  year_pred = 2018, outcome = "Past Year IPV") 
agg_cnt_past <- pool_cnt(pred_past, IPV = sdat_past$IPV, outcome = "Past Year IPV", suffix = "2020-10-30_JNT") 
agg_reg_past <- pool_reg(pred_past, ll_age = 15, ul_age = 49, outcome = "Past Year IPV", suffix = "2020-10-30_JNT"); agg_reg_past
agg_sdg_past <- pool_reg_sdg(pred_past, ll_age = 15, ul_age = 49, outcome = "Past Year IPV", suffix = "2020-10-30_JNT"); agg_sdg_past
agg_who_past <- pool_reg_who(pred_past, ll_age = 15, ul_age = 49, outcome = "Past Year IPV", suffix = "2020-10-30_JNT"); agg_who_past
agg_who_past <- pool_reg_who(pred_past, ll_age = 15, ul_age = 104, outcome = "Past Year IPV", suffix = "2020-10-30_JNT"); agg_who_past
agg_glb_past <- pool_glb(pred_past, ll_age = 15, ul_age = 49, outcome = "Past Year IPV", suffix = "2020-10-30_JNT"); agg_glb_past

agg_glb_past_exc <- pool_glb(pred_past, ll_age = 15, ul_age = 49, 
                             outcome = "Past Year IPV", cnt_exclude = c("CHN"), suffix = "excl. CHN 2020-10-30_JNT", save_results = FALSE)


#' ---- New Checks - 15-49 by Cnt*Time ----
plot_time_trend(res_past, sdat_past, IPV = ipv_p, Center_t = Center_t,
                outcome = "Past Year IPV", suffix = "2020-10-30_JNT")
agg_cnt_time_past <- pool_pred_time(res_past, sdat_past, IPV = ipv_p, Center_a = Center_a, Center_t = Center_t, 
                outcome = "Past Year IPV", suffix = "2020-10-30_JNT", max_iter = 25000, age_grp = "15-49")
plot_agg_cnt_time(agg_cnt_time_past, sdat_past$IPV, outcome = "Past Year IPV", suffix = "2020-10-30_JNT")


# ---- Past - Heatmap and Map all countries  ----
agg_cnt_past <- read.csv("Estimates Age by Cnt - Past Year IPV 2020-10-30_JNT.csv")
agg_cnt_past$cnt <- countrycode::countrycode(sourcevar = agg_cnt_past$Country, origin = "iso3c", destination = "un.name.en")
agg_cnt_past$cnt[agg_cnt_past$cnt == "Lao People’s Democratic Republic"] <- "Lao People's Democratic Republic"

heatmap_vaw(agg_cnt_past, cnt_var = "cnt", outcome = "Past Year IPV", suffix = "2020-10-30_JNT", 
                        width = 8.5, height = 17)
world_map_prv(agg_cnt_past, ipv = sdat_past$IPV, age_grp = "15-49", outcome = "Past Year IPV", width = 11, height = 8.5,
                           suffix = "2020-10-30_JNT")
prv_plot_cnt(agg_cnt_past, sdat_past$IPV, age_grp = "15-49", outcome = "Past Year IPV", width = 8.5, height = 19,
                           suffix = "2020-10-30_JNT", xlim = c(0, 70))

agg_reg_past <- read.csv("Estimates by GBD Region (15-49) - Past Year IPV 2020-10-30_JNT.csv")
prv_plot_reg(agg_reg_past, sdat_past$IPV, age_grp = "15-49", outcome = "Past Year IPV", width = 8.5, height = 6,
                           suffix = "2020-10-30_JNT", xlim = c(0, 70))


# Checking that ever is always greater than past year.
all_agg <- data.frame(agg_cnt_ever, agg_cnt_past)
ifelse(sum(ifelse(all_agg$Median < all_agg$Median.1, 1, 0)) > 0, "error ever < past", "all is good!")
ifelse(sum(ifelse(all_agg$Median < all_agg$Median.1, 1, 0)) > 0, "error ever < past", "all is good!")


#' ============================================
# ---- Out-of-sample validation - COUNTRY ----
#' ============================================
library(rjags)
seeds <- seq(2001, 2020, 1)
ever_oos <- NULL
past_oos <- NULL

for (i in 1:length(seeds)) { 
set.seed(seeds[i]) 
SDat_oos <- SDat_original
unique_cnt_e <- unique(ipv_e$Country[1:SDat_oos$start_e])
cnt_to_remove_e <- sample(unique_cnt_e, size = floor(0.2 * length(unique_cnt_e)))
SDat_oos$Num_e <- ifelse(ipv_e$Country %in% cnt_to_remove_e, NA, SDat_oos$Num_e)

unique_cnt_p <- unique(ipv_p$Country[1:SDat_oos$start_p])
cnt_to_remove_p <- sample(unique_cnt_p, size = floor(0.2 * length(unique_cnt_p)))
SDat_oos$Num_p <- ifelse(ipv_p$Country %in% cnt_to_remove_p, NA, SDat_oos$Num_p)

load.module("glm")
n_chain <- 4
n_adapt <- 10000
n_burn <- 5000
n_iter <- 50000
n_thin <- 20
param_oos <- c("linpred_e", "linpred_p")

jags.fit.oos <- jags.model(textConnection(modelstring), n.adapt = n_adapt, data = SDat_oos, n.chains = n_chain)
update(jags.fit.oos, n_burn)
MCMC.oos <- jags.samples(jags.fit.oos, variable.names = param_oos, n.iter = n_iter, thin = n_thin)
linpred_oos_e <- coda::as.mcmc.list(MCMC.oos$linpred_e)
lprd_oos_e <- do.call(rbind, linpred_oos_e)[, 1:SDat_oos$start_e]
linpred_oos_p <- coda::as.mcmc.list(MCMC.oos$linpred_p)
lprd_oos_p <- do.call(rbind, linpred_oos_p)[, 1:SDat_oos$start_p]

ever_oos_i <- out_of_sample(lprd_oos_e, sdat_ever, sdat_ever$IPV, cnt_to_remove_e, cnt = TRUE, 
              outcome = "Ever IPV", suffix = "CNT 2020-10-30_JNT", width = 17) 
past_oos_i <- out_of_sample(lprd_oos_p, sdat_past, sdat_past$IPV, cnt_to_remove_p, cnt = TRUE, 
              outcome = "Past Year IPV", suffix = "CNT 2020-10-30_JNT", width = 17) 

ever_oos <- rbind(ever_oos, ever_oos_i)
past_oos <- rbind(past_oos, past_oos_i)
}

ever_oos_cnt <- ever_oos
past_oos_cnt <- past_oos

write.csv(ever_oos_cnt, file = "Out-of-sample - Ever IPV COUNTRIES 2020-10-30_JNT.csv")
write.csv(past_oos_cnt, file = "Out-of-sample - Past IPV COUNTRIES 2020-10-30_JNT.csv")

ever_oos_cnt
apply(ever_oos_cnt[, names(ever_oos_cnt) %in% c("me", "mae", "lower", "upper", "cov")], 2, FUN = median)
past_oos_cnt
apply(past_oos_cnt[, names(past_oos_cnt) %in% c("me", "mae", "lower", "upper", "cov")], 2, FUN = median)


#' ============================================
# ---- Out-of-sample validation - STUDIES ----
#' ============================================
library(rjags)
seeds <- seq(2001, 2020, 1)
ever_oos <- NULL
past_oos <- NULL

for (i in 1:length(seeds)) { 
set.seed(seeds[i]) 
SDat_oos <- SDat_original
unique_stu_e <- unique(ipv_e$Study[1:SDat_oos$start_e])
stu_to_remove_e <- sample(unique_stu_e, size = floor(0.2 * length(unique_stu_e)))
SDat_oos$Num_e <- ifelse(ipv_e$Study %in% stu_to_remove_e, NA, SDat_oos$Num_e)

unique_stu_p <- unique(ipv_p$Study[1:SDat_oos$start_p])
stu_to_remove_p <- sample(unique_stu_p, size = floor(0.2 * length(unique_stu_p)))
SDat_oos$Num_p <- ifelse(ipv_p$Study %in% stu_to_remove_p, NA, SDat_oos$Num_p)

load.module("glm")
n_chain <- 4
n_adapt <- 10000
n_burn <- 5000
n_iter <- 50000
n_thin <- 20
param_oos <- c("linpred_e", "linpred_p")
jags.fit.oos <- jags.model(textConnection(modelstring), n.adapt = n_adapt, data = SDat_oos, n.chains = n_chain)
update(jags.fit.oos, n_burn)
MCMC.oos <- jags.samples(jags.fit.oos, variable.names = param_oos, n.iter = n_iter, thin = n_thin)
linpred_oos_e <- coda::as.mcmc.list(MCMC.oos$linpred_e)
lprd_oos_e <- do.call(rbind, linpred_oos_e)[, 1:SDat_oos$start_e]
linpred_oos_p <- coda::as.mcmc.list(MCMC.oos$linpred_p)
lprd_oos_p <- do.call(rbind, linpred_oos_p)[, 1:SDat_oos$start_p]

ever_oos_i <- out_of_sample(lprd_oos_e, sdat_ever, sdat_ever$IPV, stu_to_remove_e, cnt = FALSE, 
              outcome = "Ever IPV", suffix = "STUDIES 2020-10-30_JNT", width = 17) 
past_oos_i <- out_of_sample(lprd_oos_p, sdat_past, sdat_past$IPV, stu_to_remove_p, cnt = FALSE,
              outcome = "Past Year IPV", suffix = "STUDIES 2020-10-30_JNT", width = 17) 

ever_oos <- rbind(ever_oos, ever_oos_i)
past_oos <- rbind(past_oos, past_oos_i)
}

ever_oos_stu <- ever_oos
past_oos_stu <- past_oos
write.csv(ever_oos_stu, file = "Out-of-sample - Ever IPV STUDIES 2020-10-30_JNT.csv")
write.csv(past_oos_stu, file = "Out-of-sample - Past IPV STUDIES 2020-10-30_JNT.csv")

ever_oos_stu
past_oos_stu
apply(ever_oos_stu[, names(ever_oos_stu) %in% c("me", "mae", "lower", "upper", "cov")], 2, FUN = median)
apply(past_oos_stu[, names(past_oos_stu) %in% c("me", "mae", "lower", "upper", "cov")], 2, FUN = median)


