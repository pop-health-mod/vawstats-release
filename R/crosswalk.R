ipv_no_na <- function(IPV, add = "") {
  colsel <- c("Num", "denom_imp", 
              "num_ess", "den_ess",
              "ID", "parametertype", 
              "iso3", "region", "SuperRegion", "startyr",
              "loage", "hiageI",
              "violence", "viotime", "severe", "geo", "pstat",
              "currpart", "spouseonly", "IPV_exw", "author_year", add)
  val <- na.omit(IPV[, names(IPV) %in% colsel])
  return(val)
}

npsv_no_na <- function(IPV, add = "") {
  colsel <- c("Num", "denom_imp", 
              "num_ess", "den_ess",
              "ID", "parametertype", 
              "iso3", "region", "SuperRegion", "startyr",
              "loage", "hiageI",
              "sexviotime15", "npsvever_bf", "npsvsince15_bf",
              "violence", "viotime", "severe", "geo", "pstat", "author_year", add)
  val <- na.omit(IPV[, names(IPV) %in% colsel])
  return(val)
}

meta_proc <- function(C_data, var, event.name, is.dhs = FALSE, npsv = FALSE) { 
  C_meta <- NULL
  UID <- unique(C_data$subclass)
  for (i in 1:length(unique(C_data$subclass))) {
    Dat.i <- subset(C_data, subclass == UID[i])
    event.e <- Dat.i$Num[Dat.i[, names(Dat.i) %in% var] == event.name]
    n.e <- Dat.i$denom_imp[Dat.i[, names(Dat.i) %in% var] == event.name]
    event.c <- Dat.i$Num[Dat.i[, names(Dat.i) %in% var] != event.name]
    n.c <- Dat.i$denom_imp[Dat.i[, names(Dat.i) %in% var] != event.name]
    
    ess.event.e <- Dat.i$num_ess[Dat.i[, names(Dat.i) %in% var] == event.name]
    ess.n.e <- Dat.i$den_ess[Dat.i[, names(Dat.i) %in% var] == event.name]
    ess.event.c <- Dat.i$num_ess[Dat.i[, names(Dat.i) %in% var] != event.name]
    ess.n.c <- Dat.i$den_ess[Dat.i[, names(Dat.i) %in% var] != event.name]
    
    if (is.dhs == TRUE) {
      C_meta.i <- data.frame(Study = Dat.i$Study[1], ID = Dat.i$surveyid[1], author_year = Dat.i$author_year[1],
                             iso3 = Dat.i$iso3[1], event.e, n.e, event.c, n.c,
                             ess.event.e, ess.n.e, ess.event.c, ess.n.c, 
                             SuperRegion = Dat.i$SuperRegion[1])      
    } else {
      if (npsv == TRUE) {
      C_meta.i <- data.frame(Study = Dat.i$Study[1], ID = Dat.i$ID[1], author_year = Dat.i$author_year[1],
                             iso3 = Dat.i$iso3[1], 
                             event.e, n.e, event.c, n.c, 
                             ess.event.e, ess.n.e, ess.event.c, ess.n.c, 
                             SuperRegion = Dat.i$SuperRegion[1],
                             severe = Dat.i$severe[1], geo = Dat.i$geo[1], 
                             pstat = Dat.i$pstat[1], viotime = Dat.i$viotime[1], violence = Dat.i$violence[1],
                             loage = Dat.i$loage[1], hiage = Dat.i$hiageI[1], startyr = Dat.i$startyr[1])        
      } 
      if (npsv == FALSE) {
        C_meta.i <- data.frame(Study = Dat.i$Study[1], ID = Dat.i$ID[1], author_year = Dat.i$author_year[1],
                               iso3 = Dat.i$iso3[1], 
                               event.e, n.e, event.c, n.c, 
                               ess.event.e, ess.n.e, ess.event.c, ess.n.c, 
                               SuperRegion = Dat.i$SuperRegion[1],
                               currpart = Dat.i$currpart[1], spouseonly = Dat.i$spouseonly[1],
                               severe = Dat.i$severe[1], geo = Dat.i$geo[1], IPV_exw = Dat.i$IPV_exw[1],
                               pstat = Dat.i$pstat[1], viotime = Dat.i$viotime[1], violence = Dat.i$violence[1],
                               loage = Dat.i$loage[1], hiage = Dat.i$hiageI[1], startyr = Dat.i$startyr[1])           
      }
    }
    C_meta <- rbind(C_meta, C_meta.i)
  }
  C_meta$SuperRegion <- droplevels(C_meta$SuperRegion)
  return(C_meta)
}

meta_proc_geo <- function(C_data, var, event.name) { 
  C_meta <- NULL
  UID <- unique(C_data$ID)
  for (i in 1:length(unique(C_data$ID))) {
    Dat.i <- subset(C_data, ID == UID[i])
    event.e <- Dat.i$Num[Dat.i[, names(Dat.i) %in% var] == event.name]
    n.e <- Dat.i$denom_imp[Dat.i[, names(Dat.i) %in% var] == event.name]
    event.c <- Dat.i$Num[Dat.i[, names(Dat.i) %in% var] != event.name]
    n.c <- Dat.i$denom_imp[Dat.i[, names(Dat.i) %in% var] != event.name]
    
    ess.event.e <- Dat.i$num_ess[Dat.i[, names(Dat.i) %in% var] == event.name]
    ess.n.e <- Dat.i$den_ess[Dat.i[, names(Dat.i) %in% var] == event.name]
    ess.event.c <- Dat.i$num_ess[Dat.i[, names(Dat.i) %in% var] != event.name]
    ess.n.c <- Dat.i$den_ess[Dat.i[, names(Dat.i) %in% var] != event.name]
    
      C_meta.i <- data.frame(Study = Dat.i$Study[1], ID = Dat.i$ID[1], author_year = Dat.i$author_year[1],
                             iso3 = Dat.i$iso3[1], 
                           event.e, n.e, event.c, n.c, 
                           ess.event.e, ess.n.e, ess.event.c, ess.n.c, 
                           SuperRegion = Dat.i$SuperRegion[1],
                           currpart = Dat.i$currpart[1], spouseonly = Dat.i$spouseonly[1],
                           severe = Dat.i$severe[1], geo = Dat.i$geo[1], 
                           pstat = Dat.i$pstat[1], violence = Dat.i$violence[1], IPV_exw = Dat.i$IPV_exw[1],
                           loage = Dat.i$loage[1], hiage = Dat.i$hiageI[1], startyr = Dat.i$startyr[1])
    C_meta <- rbind(C_meta, C_meta.i)
  }
  C_meta$SuperRegion <- droplevels(C_meta$SuperRegion)
  return(C_meta)
}

meta_proc_geo_npsv <- function(C_data, var, event.name) { 
  C_meta <- NULL
  UID <- unique(C_data$ID)
  for (i in 1:length(unique(C_data$ID))) {
    Dat.i <- subset(C_data, ID == UID[i])
    event.e <- Dat.i$Num[Dat.i[, names(Dat.i) %in% var] == event.name]
    n.e <- Dat.i$denom_imp[Dat.i[, names(Dat.i) %in% var] == event.name]
    event.c <- Dat.i$Num[Dat.i[, names(Dat.i) %in% var] != event.name]
    n.c <- Dat.i$denom_imp[Dat.i[, names(Dat.i) %in% var] != event.name]
    
    ess.event.e <- Dat.i$num_ess[Dat.i[, names(Dat.i) %in% var] == event.name]
    ess.n.e <- Dat.i$den_ess[Dat.i[, names(Dat.i) %in% var] == event.name]
    ess.event.c <- Dat.i$num_ess[Dat.i[, names(Dat.i) %in% var] != event.name]
    ess.n.c <- Dat.i$den_ess[Dat.i[, names(Dat.i) %in% var] != event.name]
    
    C_meta.i <- data.frame(Study = Dat.i$Study[1], ID = Dat.i$ID[1], author_year = Dat.i$author_year[1],
                           iso3 = Dat.i$iso3[1], 
                           event.e, n.e, event.c, n.c, 
                           ess.event.e, ess.n.e, ess.event.c, ess.n.c, 
                           SuperRegion = Dat.i$SuperRegion[1],
                           severe = Dat.i$severe[1], geo = Dat.i$geo[1], 
                           pstat = Dat.i$pstat[1], violence = Dat.i$violence[1],
                           loage = Dat.i$loage[1], hiage = Dat.i$hiageI[1], startyr = Dat.i$startyr[1])
    C_meta <- rbind(C_meta, C_meta.i)
  }
  C_meta$SuperRegion <- droplevels(C_meta$SuperRegion)
  return(C_meta)
}

sel_meta_unique <- function(C_meta, var, npsv = FALSE) {
  unid <- unique(C_meta$ID)
  new_data <- NULL
  
  if (npsv == FALSE) {
    for (id in unid) {
    dat.i <- C_meta[C_meta$ID == id, ]
    if (var != "pstat" & nrow(dat.i) > 1 & any(dat.i$pstat == "Ever-partnered")) {
      dat.i <- dat.i[dat.i$pstat == "Ever-partnered", ]
    }
    if (var != "violence" & nrow(dat.i) > 1 & any(dat.i$viotime == "Physical and/or sexual IPV")) {
      dat.i <- dat.i[dat.i$viotime == "Physical and/or sexual IPV", ]
    }
    if (var != "currpart" & nrow(dat.i) > 1 & any(dat.i$currpart == "Not only asked violence from current partner")) {
      dat.i <- dat.i[dat.i$currpart == "Not only asked violence from current partner", ]
    }
    if (var != "spouseonly" & nrow(dat.i) > 1 & any(dat.i$spouseonly == "Not only asked violence from spouse")) {
      dat.i <- dat.i[dat.i$spouseonly == "Not only asked violence from spouse", ]
    }
    if (var != "geo" & nrow(dat.i) > 1 & any(dat.i$geo == "National")) {
      dat.i <- dat.i[dat.i$geo == "National", ]
    }
    if (var != "severe" & nrow(dat.i) > 1 & any(dat.i$severe == "Not only severe violence")) {
      dat.i <- dat.i[dat.i$severe == "Not only severe violence", ]
    }
    if (var != "viotime" & nrow(dat.i) > 1 & any(dat.i$viotime == "Ever")) {
      dat.i <- dat.i[dat.i$viotime == "Ever", ]
    }
    if (var != "IPV_exw" & nrow(dat.i) > 1 & any(dat.i$IPV_exw == "All")) {
      dat.i <- dat.i[dat.i$IPV_exw == "All", ]
    }
    new_data <- rbind(new_data, dat.i)
  }
  new_data <- new_data[order(new_data$iso3, new_data$startyr), ]
  }
  
  if (npsv == TRUE) {
  for (id in unid) {
    dat.i <- C_meta[C_meta$ID == id, ]
    if (var != "pstat" & nrow(dat.i) > 1 & any(dat.i$pstat == "All women")) {
      dat.i <- dat.i[dat.i$pstat == "All women", ]
    }
    if (var != "geo" & nrow(dat.i) > 1 & any(dat.i$geo == "National")) {
      dat.i <- dat.i[dat.i$geo == "National", ]
    }
    if (var != "severe" & nrow(dat.i) > 1 & any(dat.i$severe == "Not only severe violence")) {
      dat.i <- dat.i[dat.i$severe == "Not only severe violence", ]
    }
    if (var != "viotime" & nrow(dat.i) > 1 & any(dat.i$viotime == "Ever")) {
      dat.i <- dat.i[dat.i$viotime == "Ever", ]
    }
    if (nrow(dat.i) > 1 & any(dat.i$npsvsince15_bf == "Since 15_bf")) {
      dat.i <- dat.i[dat.i$npsvsince15_bf == "Since 15_bf", ]
    }
    if (nrow(dat.i) > 1 & any(dat.i$sexviotime15 == "Since 15")) {
      dat.i <- dat.i[dat.i$sexviotime15 == "Since 15", ]
    } 
    if (nrow(dat.i) > 1 & any(dat.i$npsvever_bf == "Ever_bf")) {
      dat.i <- dat.i[dat.i$npsvever_bf == "Ever_bf", ]
    }
    new_data <- rbind(new_data, dat.i)
  }
  new_data <- new_data[order(new_data$iso3, new_data$startyr), ]
  }
  return(new_data)
  }
  
aggregate_age_group_for_meta <- function(cmeta) {
  uid <- unique(cmeta$ID)
  val <- NULL
  for (i in uid) {
    cmeta_i <- subset(cmeta, ID == i)
    if (nrow(cmeta_i > 1)) { 
    cmeta_agg <- cmeta_i[1, ]  
    cmeta_agg$event.e <- sum(cmeta_i$event.e) 
    cmeta_agg$n.e <- sum(cmeta_i$n.e) 
    cmeta_agg$event.c <- sum(cmeta_i$event.c) 
    cmeta_agg$n.c <- sum(cmeta_i$n.c) 
    cmeta_agg$ess.event.e <- sum(cmeta_i$ess.event.e) 
    cmeta_agg$ess.n.e <- sum(cmeta_i$ess.n.e) 
    cmeta_agg$ess.event.c <- sum(cmeta_i$ess.event.c) 
    cmeta_agg$ess.n.c <- sum(cmeta_i$ess.n.c) 
    cmeta_agg$loage <- min(cmeta_i$loage)
    cmeta_agg$hiage <- max(cmeta_i$hiage)
    cmeta_i <- cmeta_agg
    }
    val <- rbind(val, cmeta_i) 
  }
  return(val)
}

match_geo <- function(IPV, npsv = FALSE) {
  IPV$id_match <- paste(IPV$ID, IPV$severe, IPV$viotime, IPV$startyr, IPV$pstat, IPV$spouseonly, IPV$currpart, IPV$loage, IPV$hiageI, IPV$IPV_exw)
  if (npsv == TRUE) {
    IPV$id_match <- paste(IPV$ID, IPV$severe, IPV$viotime, IPV$startyr, IPV$pstat, IPV$loage, IPV$hiageI,
                          IPV$sexviotime15, IPV$npsvever_bf, IPV$npsvsince15_bf)
  }
  
  UID <- unique(IPV$id_match)
  my_match <- NULL
  for (i in 1:length(UID)) {
    study.i <- subset(IPV, id_match == UID[i])
    study.i$subclass <- study.i$Study <- i
    if (any(is.na(study.i$geo))) { next }
    if (any(study.i$geo == 'National' | study.i$geo == 'Mixed') &
       any(study.i$geo == 'Urban') & any(study.i$geo == 'Rural')) {
      my_match <- rbind(my_match, study.i)
    }
  }
  return(my_match)
  }

out_meta <- function(MA) {
  log_or_sr <- data.frame(SuperR = MA$bylevs, LogOR = MA$TE.fixed.w, 
                             OR = exp(MA$TE.fixed.w), SE = MA$seTE.fixed.w, 
                             LogORr = MA$TE.random.w, SEr = MA$seTE.random.w)
  # we remove regions where there are fewer than 3 studies
  select <- NULL
  for (i in MA$bylevs) {
    if (length(unique(MA$data$.studlab[MA$data$.byvar %in% i])) > 3) {
    select <- c(select, i) }
  }
  log_or_sr <- log_or_sr[log_or_sr$SuperR %in% select, ] 
  log_or_all <- data.frame(SuperR = "Overall", LogOR = MA$TE.fixed, 
                        OR = exp(MA$TE.fixed), SE = MA$seTE.fixed, 
                        LogORr = MA$TE.random, SEr = MA$seTE.random)
  log_or <- rbind(log_or_sr, log_or_all)
  return(log_or)
}

out_meta_glmer_geo <- function(fit, geo = "rural") {

   if (geo == "urban") {
      log_or_urb_o <- coef(summary(fit))["urban", "Estimate"]
      se_urb_o <- coef(summary(fit))["urban", "Std. Error"]
      log_or_urb_r <- coef(fit)$SuperRegion[, c("urban")]
      varfix_urb <- vcov(fit)["urban", "urban"]
      varcm_urb <- attr(lme4::ranef(fit, condVar = TRUE)$SuperRegion, "postVar")[2, 2, ]

      se_urb <- sqrt(varfix_urb + varcm_urb)
      log_or_urban_r <- data.frame(SuperR = rownames(coef(fit)$SuperRegion), LogOR = log_or_urb_r, 
                             OR = exp(log_or_urb_r), SE = se_urb, 
                             LogORr = log_or_urb_r, SEr = se_urb)
      log_or_urban_r <- log_or_urban_r[order(log_or_urban_r$SuperR), ]
      log_or_urban_o <- data.frame(SuperR = "Overall", LogOR = log_or_urb_o, 
                        OR = exp(log_or_urb_o), SE = se_urb_o, 
                        LogORr = log_or_urb_o, SEr = se_urb_o)
      log_or_urban <- rbind(log_or_urban_r, log_or_urban_o)
      return(log_or_urban)
   }
  if (geo == "rural") {
      log_or_rur_o <- coef(summary(fit))["rural", "Estimate"]
      se_rur_o <- coef(summary(fit))["rural", "Std. Error"]   
      log_or_rur_r <- coef(fit)$SuperRegion[, c("rural")]  
      varfix_rur <- vcov(fit)["rural", "rural"]
      varcm_rur <- attr(lme4::ranef(fit, condVar = TRUE)$SuperRegion, "postVar")[1, 1, ]
      se_rur <- sqrt(varfix_rur + varcm_rur)
      
      log_or_rural_r <- data.frame(SuperR = rownames(coef(fit)$SuperRegion), LogOR = log_or_rur_r, 
                             OR = exp(log_or_rur_r), SE = se_rur, 
                             LogORr = log_or_rur_r, SEr = se_rur)
      log_or_rural_r <- log_or_rural_r[order(log_or_rural_r$SuperR), ]
      log_or_rural_o <- data.frame(SuperR = "Overall", LogOR = log_or_rur_o, 
                        OR = exp(log_or_rur_o), SE = se_rur_o, 
                        LogORr = log_or_rur_o, SEr = se_rur_o)
      log_or_rural <- rbind(log_or_rural_r, log_or_rural_o)
      return(log_or_rural) }
}

ranef_sr <- function(fit, matched_data) {
  val <- matched_data[matched_data$X == 1, c("Study", "SuperRegion", "iso3", "ID")]
  val <- val[order(val$Study), ]
  val$re <- coef(fit)$Study[, "X"]
  val$se <- sqrt(as.vector(attr(lme4::ranef(fit, condVar = TRUE)[[1]], "postVar")$X))
  val$or <- exp(val$re)
  val$lci <- exp(val$re - qnorm(0.975) * val$se)
  val$uci <- exp(val$re + qnorm(0.975) * val$se)
  val <- val[order(val$SuperRegion, val$iso3, val$Study), ]
}

meta_sr <- function(fit, val) {
  unique_sr <- unique(val$SuperRegion)
  val$w <- 1 / (val$se^2)
  tau <- (fit@pp$theta[2])^2
  valsr <- NULL
  for (i in unique_sr) {
    val.i <- val[val$SuperRegion %in% i, ]
    if (nrow(val.i) > 1) {
    val.i$w_re <- 1 / (val.i$se^2 + tau)
    mean.i <- sum(val.i$w_re * val.i$re) / sum(val.i$w_re)
    valsr.i <- data.frame(superregion = i, 
                          logor = mean.i, se = sqrt(1 / sum(val.i$w_re)),
                          or = exp(mean.i), 
                          lci = exp(mean.i - sqrt(1 / sum(val.i$w_re)) * qnorm(0.975)),
                          uci = exp(mean.i + sqrt(1 / sum(val.i$w_re)) * qnorm(0.975)))
    } else {
    valsr.i <- data.frame(superregion = i, 
                            logor = val.i$re, se = val.i$se,
                            or = exp(val.i$re), 
                            lci = val.i$lci,
                            uci = val.i$uci)
      
    }
    valsr <- rbind(valsr, valsr.i)
  }
  val.o <- data.frame(superregion = "Overall", 
             logor = coef(summary(fit))[, "Estimate"], se = coef(summary(fit))[, "Std. Error"],
             or = exp(coef(summary(fit))[, "Estimate"]), 
             lci = exp(coef(summary(fit))[, "Estimate"] - qnorm(0.975) * coef(summary(fit))[, "Std. Error"]),
             uci = exp(coef(summary(fit))[, "Estimate"] + qnorm(0.975) * coef(summary(fit))[, "Std. Error"]))
  valsr <- rbind(valsr, val.o)
  return(valsr)
}


# plot_meta <- function(val) { 
#   dotchart(val$or, labels = val$ID, xlim = c(0, 1))
#   for (i in 1:nrow(val)){
#     lines(x = c(val$lci[i], val$uci[i]), y = c(i, i))
#   }
#   }