
#' Utility functions

#' Function to recode clustering variable into unique numeric id from 1-to-N
recode.cluster <- function(Var) {
  id <- unique(Var)
  n <- length(id)
  clust.id <- rep('NA', length(Var))
  for (i in 1:n) {
    id.i <- id[i]
    clust.id[which(Var == id.i)] <- i
  }
  as.numeric(clust.id)	}

# inverse logit
inv.logit <- function(x) {
  val <- 1 / (1 + exp(-x))
  return(val)
}

#' Function to calculate the effective sample size from a survey from its CI
neff <- function(p_hat, lci, uci, wilson = FALSE) {
  if (any(is.na(c(p_hat, lci, uci)))) { return(NA) }
  if ( any(p_hat < 0) | any(p_hat > 1) | any(lci < 0) | any(lci > 1) | any(uci < 0) | any(uci > 1) ) {
    stop('Error, must be between 0 and 1') }
  
  n_ess <- NULL
  
  if (wilson == FALSE) {
  for(i in 1:length(p_hat)) {
    if (p_hat[i] == 0) { p_hat[i] <- 1e-3 }
    nlci <- (p_hat[i] * (1 - p_hat[i])) / ((p_hat[i] - lci[i]) / qnorm(0.975))^2
    nuci <- (p_hat[i] * (1 - p_hat[i])) / ((uci[i] - p_hat[i]) / qnorm(0.975))^2
    n.i <- mean(c(nlci, nuci))
    n_ess <- c(n_ess, ceiling(n.i)) }}
  
  if (wilson == TRUE) {
  for(i in 1:length(p_hat)) {
    n.i <- - ((uci[i] - 1) * uci[i] * qnorm(0.975)^2) / (uci[i] - p_hat[i])^2
    n_ess <- c(n_ess, ceiling(n.i)) }}
  
  return(n_ess)
}

  abstract_pd <- function(x) {
    # Look at https://rdrr.io/cran/rjags/src/R/dic.R for details
    dev <- sum(x$deviance)
    psum <- sum(x[[2]])
    pD <- as.numeric(dev) + as.numeric(psum)
    return(round(pD, 2))
  }
  
#' Function to attribute particular OR to adjustment factors by super region or overall
OR_to_use <- function(GBD_SR, LogOR, lhs = FALSE, index) {
  match_id <- match(x = GBD_SR, table = LogOR$SuperR)
  match_id <- ifelse(is.na(match_id), which(LogOR$SuperR == 'Overall'), match_id)
  if (lhs == FALSE) {
    val <- LogOR$LogORr[match_id]
  } else {
    val <- LogOR[match_id, index + 1]
  }
  return(val)
}

#' Function to attribute particular OR to adjustment factors by super region or overall
OR_to_use_lhs <- function(LogOR, n = 10, truncate = TRUE) {
  or_mn <- LogOR$LogORr
  or_sd <- LogOR$SEr
  lhs_rnd <- lhs::randomLHS(n = n, k = length(or_mn))
  or_rnd <- NULL
  for (i in 1:length(or_mn)) {
    or_rnd.i <- qnorm(lhs_rnd[, i], mean = or_mn[i], sd = or_sd[i])
    if (truncate == TRUE) {
      or_rnd.i <- truncnorm::qtruncnorm(p = lhs_rnd[, i], a = - Inf, b = 0, 
                             mean = or_mn[i], sd = or_sd[i])
    }
    or_rnd <- cbind(or_rnd, or_rnd.i)
  }

 val <- data.frame(SuperR = as.character(LogOR$SuperR), t(as.matrix(or_rnd)))
 return(val)
}
 


#' Indexing Age for use in spline (lower bound).
indexing_age_li <- function(loage) {
  li <- ifelse(loage <= 18, 1,
        ifelse(loage > 18 & loage <= 22, 2, # 20
        ifelse(loage > 22 & loage <= 27, 3, # 25
        ifelse(loage > 27 & loage <= 32, 4,
        ifelse(loage > 32 & loage <= 37, 5,
        ifelse(loage > 37 & loage <= 42, 6,
        ifelse(loage > 42 & loage <= 47, 7,
        ifelse(loage > 47 & loage <= 52, 8,
        ifelse(loage > 52 & loage <= 57, 9,
        ifelse(loage > 57 & loage <= 62, 10,
        ifelse(loage > 62, 11, NA) ))))))))))
  return(li) }


#' Indexing Age for use in spline (upper bound).
indexing_age_ui <- function(hiageI) {
  ui <- ifelse(hiageI < 22, 1, #19
        ifelse(hiageI >= 22 & hiageI < 27, 2, #24
        ifelse(hiageI >= 27 & hiageI < 32, 3,
        ifelse(hiageI >= 32 & hiageI < 37, 4,
        ifelse(hiageI >= 37 & hiageI < 42, 5,
        ifelse(hiageI >= 42 & hiageI < 47, 6,
        ifelse(hiageI >= 47 & hiageI < 52, 7,
        ifelse(hiageI >= 52 & hiageI < 57, 8,
        ifelse(hiageI >= 57 & hiageI < 62, 9,
        ifelse(hiageI >= 62 & hiageI < 67, 10,
        ifelse(hiageI >= 67, 11, NA)))))))))))
  return(ui) }

#' Indexing GBD region.
indexing_region <- function(region) {
  UniqueR <- c("Asia, East", "Asia, Southeast", "Oceania", "Asia, Central",
              "Europe, Central", "Europe, Eastern", "Asia Pacific, High Income",
              "Australasia", "Europe, Western", "Latin America, Southern",
              "North America, High Income", "Caribbean", "Latin America, Andean",
              "Latin America, Central", "Latin America, Tropical", "North Africa/Middle East",
              "Asia, South", "Sub-Saharan Africa, Central", "Sub-Saharan Africa, East",
              "Sub-Saharan Africa, Southern", "Sub-Saharan Africa, West")
  col_index <- ifelse(as.character(region) == UniqueR[1], 1,
                ifelse(as.character(region) == UniqueR[2], 2,
                ifelse(as.character(region) == UniqueR[3], 3,
                ifelse(as.character(region) == UniqueR[4], 4,
                ifelse(as.character(region) == UniqueR[5], 5,
                ifelse(as.character(region) == UniqueR[6], 6,
                ifelse(as.character(region) == UniqueR[7], 7,
                ifelse(as.character(region) == UniqueR[8], 8,
                ifelse(as.character(region) == UniqueR[9], 9,
                ifelse(as.character(region) == UniqueR[10], 10,
                ifelse(as.character(region) == UniqueR[11], 11,
                ifelse(as.character(region) == UniqueR[12], 12,
                ifelse(as.character(region) == UniqueR[13], 13,
                ifelse(as.character(region) == UniqueR[14], 14,
                ifelse(as.character(region) == UniqueR[15], 15,
                ifelse(as.character(region) == UniqueR[16], 16,
                ifelse(as.character(region) == UniqueR[17], 17,
                ifelse(as.character(region) == UniqueR[18], 18,
                ifelse(as.character(region) == UniqueR[19], 19,
                ifelse(as.character(region) == UniqueR[20], 10, 21))))))))))))))))))))
  return(col_index) }


#' Indexing GBD region.
class_reg_to_superreg <- function(region) {
super_region <- plyr::revalue(region, c(
  "Asia Pacific, High Income" = 'High Income',
  "Asia, Central" = 'Central Europe, Eastern Europe & Central Asia',
  "Asia, East" = 'South-East Asia, East Asia & Oceania',
  "Asia, South" = 'South Asia',
  "Asia, Southeast" = 'South-East Asia, East Asia & Oceania',
  "Australasia" = 'High Income',
  "Caribbean" = 'Latin America & Caribbean',
  "Europe, Central" = 'Central Europe, Eastern Europe & Central Asia',
  "Europe, Eastern" = 'Central Europe, Eastern Europe & Central Asia',
  "Europe, Western" = 'High Income',
  "Latin America, Andean" = 'Latin America & Caribbean',
  "Latin America, Central" = 'Latin America & Caribbean',
  "Latin America, Southern" = 'High Income',
  "Latin America, Tropical" = 'Latin America & Caribbean',
  "North Africa/Middle East" = 'North Africa & Middle East',
  "North America, High Income" = 'High Income',
  "Oceania" = 'South-East Asia, East Asia & Oceania',
  "Sub-Saharan Africa, Central" = 'Sub-Saharan Africa',
  "Sub-Saharan Africa, East" = 'Sub-Saharan Africa',
  "Sub-Saharan Africa, Southern" = 'Sub-Saharan Africa',
  "Sub-Saharan Africa, West" = 'Sub-Saharan Africa'))
  return(super_region) }

recode_gbd <- function(region) {
  reclassed_region <- plyr::revalue(region, c(  
        "Andean Latin America" = "Latin America, Andean", 
        "Australasia" = "Australasia", 
        "Caribbean" = "Caribbean", 
        "Central Asia" = "Asia, Central", 
        "Central Europe" = "Europe, Central", 
        "Central Latin America" = "Latin America, Central", 
        "Central Sub-Saharan Africa" = "Sub-Saharan Africa, Central", 
        "East Asia" = "Asia, East", 
        "Eastern Europe" = "Europe, Eastern", 
        "Eastern Sub-Saharan Africa" = "Sub-Saharan Africa, East", 
        "High-income Asia Pacific" = "Asia Pacific, High Income", 
        "High-income North America" = "North America, High Income", 
        "North Africa and Middle East" = "North Africa/Middle East", 
        "Oceania" = "Oceania", 
        "South Asia" = "Asia, South",
        "Southeast Asia" = "Asia, Southeast", 
        "Southern Latin America" = "Latin America, Southern", 
        "Southern Sub-Saharan Africa" = "Sub-Saharan Africa, Southern", 
        "Tropical Latin America" = "Latin America, Tropical", 
        "Western Europe" = "Europe, Western", 
        "Western Sub-Saharan Africa" = "Sub-Saharan Africa, West"))
  return(reclassed_region)
}



#' Wolrd Map data availability.
world_map_data <- function(study_id, iso3, outcome = "Ever IPV", width = 11, height = 8.5,
                           to_plot = TRUE, suffix = "") {
    if (!require(gpclib)) install.packages("gpclib", type="source")
      #library(gpclib)
      #gpclibPermit()
      
    data(world_map)
    data("who_lines")
    data("who_poly")
    world_map$ISO_A3 <- world_map$iso3
    map <- world_map

    availability <- aggregate(study_id, by = list(iso3), FUN = function(x) (length(unique(x))))
    availability$iso3 <- as.character(availability$Group.1)
    c_no_data <- world_map$ISO_A3 %in% availability$iso3
    no_data <- data.frame(Group.1 = world_map$ISO_A3[!c_no_data], x = 0,
                          iso3 = world_map$ISO_A3[!c_no_data])
    #no_data <- no_data[-which(no_data$iso3 == "-99")]
    availability <- rbind(availability, no_data)
    availability$iso3 <- as.character(availability$iso3)

    availability$cat <- ifelse(availability$x >= 4, "4 or +", availability$x)
    #col <- wesanderson::wes_palette(name = "FantasticFox1", n = 5, type = 'discrete')
    #col <- nord::nord(palette = "lumina", n = 4)
    #col <- nationalparkcolors::park_palette("Acadia", n = 4)
    #col <- NineteenEightyR::sunset1(n = 4, alpha = 0.8)
    col <- wesanderson::wes_palette("Zissou1", 5, type = "continuous")

    cols <- c("0" = "grey85", "1" = col[1], "2" = col[2],
              "3" = col[3], "4 or +" = col[4])
    
    gg <- ggplot2::ggplot()
    gg <- gg + ggplot2::geom_map(data = map, map = map,
                        fill = "#ffffff", color = NA,
                        ggplot2::aes(x = map$long, y = map$lat, map_id = map$id, group = map$group))
    gg <- gg + ggplot2::geom_map(data = availability, map = map, color = "white", size = 0.15,
                        ggplot2::aes(fill = cat, group = iso3, map_id = iso3))
    gg <- gg + ggplot2::geom_map(data = who_poly, map = who_poly, fill = "white", color = NA,
                        ggplot2::aes(x = long, y = lat, map_id = id))
    gg <- gg + ggplot2::geom_path(data = who_lines, color = "grey85", size = 0.35,
                        ggplot2::aes(x = long, y = lat, group = as.numeric(group)))
    gg <- gg + ggplot2::scale_fill_manual(values = cols,
                name = 'Number \n of Studies')
    gg <- gg + ggplot2::labs(x = '', y = '',
                    title = paste("Data Availability ", outcome))
    gg <- gg + ggplot2::coord_equal(ratio = 1)
    gg <- gg + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                  axis.text.y = ggplot2::element_blank(),
                  axis.ticks = ggplot2::element_blank(),
                  rect = ggplot2::element_blank()) +
               ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    gg <- gg + ggplot2::theme(legend.position = "right")
    gg <- gg + ggplot2::theme(plot.title = ggplot2::element_text(size = 18))
    
    ggplot2::ggsave(filename = paste("Map. Data Availability - ", outcome, " ", suffix, ".pdf", sep = ""),
           device = "pdf",  width = width, height = height)

    if (to_plot == TRUE) { print(gg) }

}


#' Selecting Optimal Set
optimal_set <- function(data, outcome = "Ever IPV") {
  if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV",
                       "Past Year NPSV"))) { stop("Incorrect outcome name") }


  if (outcome == "Ever IPV") {
    # special case of Egypt DHS with weird age group.
    idx <- which(data$iso3 == "EGY" & data$startyr == 2005 & data$severe == "Only severe violence")
    data <- data[-idx, ]
    idx <- which(data$iso3 == "EGY" & data$startyr == 2014 & data$violence != "Physical and/or sexual IPV")
    data <- data[-idx, ]  
    # special case of Egypt ESP with weird age group... there are weird non-overlapping age groups. we have 60-99 but also 55-64, 65-74, and 75-99
    idx <- which(data$iso3 == "ESP" & data$startyr == 2014 &
                 (data$agerange %in% c("55-64", "65-74", "75-NA")))
    data <- data[-idx, ]  
    # Bolivia 2003 - special case
    idx <- which(data$iso3 == "BOL" & data$startyr == 2003 & data$violence != "Physical and/or sexual IPV")
    data <- data[-idx, ]     
    # DRC 2013 - special case
    idx <- which(data$iso3 == "COD" & data$startyr == 2013 & data$violence != "Physical and/or sexual IPV")
    data <- data[-idx, ]   
    # PRY 2008 - special case
    idx <- which(data$iso3 == "PRY" & data$startyr == 2008 & data$violence != "Physical and/or sexual IPV")
    data <- data[-idx, ]  
  }
  
  if (outcome == "Past Year IPV") {
    # PRY 2008 - special case
    idx <- which(data$iso3 == "PRY" & data$startyr == 2008 & data$violence != "Physical and/or sexual IPV")
    data <- data[-idx, ]  
    # CIV 2011 - special case
    idx1 <- which(data$iso3 == "CIV" & data$startyr == 2011 & data$violence != "Physical and/or sexual IPV")
    data <- data[-idx1, ]
    idx2 <- which(data$iso3 == "CIV" & data$startyr == 2011 & data$currpart != "Asked about violence from current or most recent partner only")
    data <- data[-idx2, ]
  }  
  
  if (outcome == "NPSV") {
    # GRD 2018 - special case
    idx <- which(data$iso3 == "GRD" & data$startyr == 2018 & data$parametertype != "SexVio Ever")
    data <- data[-idx, ]  
    print('assuming that NPSV is for Ever NPSV, correct special case of GRD 2018 otherwise')
  }    
  
  # We don't put <hiageI> in the definition on purpose (see below)
  data$IDage <- paste(data$ID, data$loage, sep = '_')
  data$RowID <- seq(1:nrow(data))
  UID <- unique(data$IDage)
  vaw <- NULL


  if (outcome == "Ever IPV") {
    for (i in 1:length(UID)) {
      ipv.i <- data[data$IDage == UID[i], ]
      
      if (nrow(ipv.i) > 1) {
          # If a broad age group and a small one, we prefer age-specific information (delete maximum)
          if (length(unique(ipv.i$agerange)) > 1) {
              min_age <- ipv.i$AgeM %in% min(ipv.i$AgeM)
              ipv.i <- ipv.i[min_age, ] } 
          # We select non-severe
          if (nrow(ipv.i) > 1 & any(ipv.i$severe == "Not only severe violence")) {
            ipv.i <- ipv.i[ipv.i$severe == "Not only severe violence", ] }
          # We select "Physical and/or sexual IPV"
          if (nrow(ipv.i) > 1 & any(ipv.i$violence == "Physical and/or sexual IPV")) {
            ipv.i <- ipv.i[ipv.i$violence == "Physical and/or sexual IPV", ] }
          # --- If no "Physical and/or sexual IPV", we select "Physical IPV"
          if (nrow(ipv.i) > 1 & any(ipv.i$violence == "Physical IPV only")) {
            ipv.i <- ipv.i[ipv.i$violence == "Physical IPV only", ] }
          # We select "Ever-Partnered"
          if (nrow(ipv.i) > 1 & any(ipv.i$pstat == "Ever-partnered")) {
            ipv.i <- ipv.i[ipv.i$pstat == "Ever-partnered", ] }
          # --- If no "Ever-Partnered", we select "All-Women"
          if (nrow(ipv.i) > 1 & any(ipv.i$pstat == "All women")) {
            ipv.i <- ipv.i[ipv.i$pstat == "All women", ] }
          # We select "Not only asked violence from current partner"
          if (nrow(ipv.i) > 1 & any(ipv.i$currpart == "Not only asked violence from current partner")) {
            ipv.i <- ipv.i[ipv.i$currpart == "Not only asked violence from current partner", ] }
          # We select "Not only asked violence from spouse"
          if (nrow(ipv.i) > 1 & any(ipv.i$spouseonly == "Not only asked violence from spouse")) {
            ipv.i <- ipv.i[ipv.i$spouseonly == "Not only asked violence from spouse", ] }
      }
      vaw <- rbind(vaw, ipv.i) }}

  if (outcome == "Past Year IPV") {
    for (i in 1:length(UID)) {
      ipv.i <- data[data$IDage == UID[i], ]
      if (nrow(ipv.i) > 1) {
          # If a broad age group and a small one, we prefer age-specific information
          if (length(unique(ipv.i$agerange)) > 1) {
            min_age <- ipv.i$AgeM %in% min(ipv.i$AgeM)
            ipv.i <- ipv.i[min_age, ] } 
        # --- We select non-severe
        if (nrow(ipv.i) > 1 & any(ipv.i$severe == "Not only severe violence")) {
          ipv.i <- ipv.i[ipv.i$severe == "Not only severe violence", ] }
        # --- We select "Physical and/or sexual IPV"
        if (nrow(ipv.i) > 1 & any(ipv.i$violence == "Physical and/or sexual IPV")) {
          ipv.i <- ipv.i[ipv.i$violence == "Physical and/or sexual IPV", ] }
        # If no "Physical and/or sexual IPV", we select "Physical IPV"
        if (nrow(ipv.i) > 1 & any(ipv.i$violence == "Physical IPV only")) {
          ipv.i <- ipv.i[ipv.i$violence == "Physical IPV only", ] }
        # --- We select "Ever-Partnered"
        if (nrow(ipv.i) > 1 & any(ipv.i$pstat == "Ever-partnered")) {
          ipv.i <- ipv.i[ipv.i$pstat == "Ever-partnered", ] }
        # If no "Ever-Partnered", we select "All-Women"
        if (nrow(ipv.i) > 1 & any(ipv.i$pstat == "All women")) {
          ipv.i <- ipv.i[ipv.i$pstat == "All women", ] }
        # --- We select "Not only asked violence from current partner"
        if (nrow(ipv.i) > 1 & any(ipv.i$currpart == "Not only asked violence from current partner")) {
          ipv.i <- ipv.i[ipv.i$currpart == "Not only asked violence from current partner", ] }
        # --- We select "Not only asked violence from spouse"
        if (nrow(ipv.i) > 1 & any(ipv.i$spouseonly == "Not only asked violence from spouse")) {
          ipv.i <- ipv.i[ipv.i$spouseonly == "Not only asked violence from spouse", ] }
        # ---We select "Past Year" over "Past 2 Yr"
        if (nrow(ipv.i) > 1 & any(ipv.i$viotime == "Past yr")) {
          ipv.i <- ipv.i[ipv.i$viotime == "Past yr", ] }
        }
      vaw <- rbind(vaw, ipv.i) }}

  if (outcome == "NPSV") {
    for (i in 1:length(UID)) {
      ipv.i <- data[data$IDage == UID[i], ]
      if (nrow(ipv.i) > 1) {
          # If a broad age group and a small one, we prefer age-specific information
          if (length(unique(ipv.i$agerange)) > 1) {
            min_age <- ipv.i$AgeM %in% min(ipv.i$AgeM)
            ipv.i <- ipv.i[min_age, ] } 
          # --- We select "non-severe" of "severe"
          if (nrow(ipv.i) > 1 & any(ipv.i$severe == "Not only severe violence")) {
            ipv.i <- ipv.i[ipv.i$severe == "Not only severe violence", ] }
          # --- We select "Ever" over "Past Yr"
          if (nrow(ipv.i) > 1 & any(ipv.i$viotime == "Ever")) {
            ipv.i <- ipv.i[ipv.i$viotime == "Ever", ] }
          # --- We select "boyfriend since aged 15" over other
          if (nrow(ipv.i) > 1 & any(ipv.i$npsvsince15_bf == "Since 15_bf")) {
            ipv.i <- ipv.i[ipv.i$npsvsince15_bf == "Since 15_bf", ] } 
          if (nrow(ipv.i) > 1 & any(ipv.i$sexviotime15 == "Since 15")) {
            ipv.i <- ipv.i[ipv.i$sexviotime15 == "Since 15", ] }
          if (nrow(ipv.i) > 1 & any(ipv.i$npsvever_bf == "Ever_bf")) {
            ipv.i <- ipv.i[ipv.i$npsvever_bf == "Ever_bf", ]  } 
      }
      vaw <- rbind(vaw, ipv.i) }}

  if (outcome == "Past Year NPSV") {
    for (i in 1:length(UID)) {
      ipv.i <- data[data$IDage == UID[i], ]
      if (nrow(ipv.i) > 1) {
      # --- We select "Ever" over "Past Yr"
        if (nrow(ipv.i) > 1 & any(ipv.i$viotime == "Past yr")) {
            ipv.i <- ipv.i[ipv.i$viotime == "Past yr", ] }
      # If a broad age group and a small one, we prefer age-specific information
        if (length(unique(ipv.i$agerange)) > 1) {
            min_age <- ipv.i$AgeM %in% min(ipv.i$AgeM)
            ipv.i <- ipv.i[min_age, ] } 
      # --- We select "non-severe" of "severe"
        if (nrow(ipv.i) > 1 & any(ipv.i$severe == "Not only severe violence")) {
            ipv.i <- ipv.i[ipv.i$severe == "Not only severe violence", ] }
      # --- We select "boyfriend since aged 15" over other
        if (nrow(ipv.i) > 1 & any(ipv.i$npsvsince15_bf == "Since 15_bf")) {
            ipv.i <- ipv.i[ipv.i$npsvsince15_bf == "Since 15_bf", ] } 
        if (nrow(ipv.i) > 1 & any(ipv.i$sexviotime15 == "Since 15")) {
            ipv.i <- ipv.i[ipv.i$sexviotime15 == "Since 15", ] }
        if (nrow(ipv.i) > 1 & any(ipv.i$npsvever_bf == "Ever_bf")) {
            ipv.i <- ipv.i[ipv.i$npsvever_bf == "Ever_bf", ]  } 
      }
      vaw <- rbind(vaw, ipv.i) }}  
  
      # In all cases, we remove urban/rural if national is present.
      study_id <- unique(data$ID)
      vaw_final <- NULL
      for (i in 1:length(study_id)) {
        ipv.i <- vaw[vaw$ID == study_id[i], ]
        geo_rep <- unique(ipv.i$geo)
        if (length(geo_rep) > 1 & "National" %in% geo_rep) {
          ipv.i <- ipv.i[ipv.i$geo == "National", ] }
        geo_rep <- unique(ipv.i$geo)
        # If Rural/Urban/Mixed are available, we only select Mixed
        if (length(geo_rep) > 1 & "Mixed" %in% geo_rep) {
          ipv.i <- ipv.i[ipv.i$geo == "Mixed", ] }
        vaw_final <- rbind(vaw_final, ipv.i)
        }

  return(vaw_final)
}


#' Removing larger age group if age-disagregated data available
remove_larger_agegrp <- function(data) {
    ToRemove <- NULL
    data$RowID <- seq(1, dim(data)[1], 1)
    data$NewID <- with(data, paste(ID, geo))
    UID <- unique(data$NewID)
    for (i in 1:length(UID)) {
      Dat.l <- subset(data, NewID == UID[i])
      Min <- min(Dat.l$loage)
      Max <- max(Dat.l$hiageI)
      LargerAgeClass <- which(Dat.l$loage == Min & abs(Dat.l$hiageI - Max) <= 1)
      ID <- Dat.l$RowID[LargerAgeClass]
      if (nrow(Dat.l) > 1 & length(LargerAgeClass) > 0) { 
        ToRemove <- c(ToRemove, ID) 
      # If they are remaining wider age groups, we remove them too (i.e., 15-49 and 15-59 with 5-year age group in between)
       if (length(LargerAgeClass) > 0) {
         Dat.lj <- Dat.l[-LargerAgeClass, ] } 
       if (length(LargerAgeClass) == 0) { 
        Dat.lj <- Dat.l }
        if (nrow(Dat.lj) > 1 & (length(unique(Dat.lj$loage)) != nrow(Dat.lj))) {
          dup_idx <- which(Dat.lj$loage %in% Dat.lj$loage[duplicated(Dat.lj$loage)])
          Dat.lji <- Dat.lj[dup_idx, ]
          dup_idz <- unique(Dat.lji$loage)
          for (z in 1:length(dup_idz)) {
            Dat.ljz <- Dat.lji[Dat.lji$loage == dup_idz[z], ]
            dup_z <- Dat.ljz$RowID[which.max(Dat.ljz$AgeW)]
            if (length(dup_z) > 0) {
              ToRemove <- c(ToRemove, dup_z)
            }}}}
      # if we have 35-39 and 45-49 but have other age groups that encompasses them.
        if (nrow(Dat.l) > 1 & any(Dat.l$agerange == "35-39") & any(Dat.l$agerange == "30-39")) {
          dup_rem1 <- Dat.l$RowID[Dat.l$agerange == "35-39"]
          print(paste("Removing 35-39 and keeping 30-39 in ", Dat.l$iso3[1], " (", Dat.l$Time[1], ")", sep = ""))
          ToRemove <- c(ToRemove, dup_rem1)
        }
       if (nrow(Dat.l) > 1 & any(Dat.l$agerange == "45-49") & any(Dat.l$agerange == "40-49")) {
          dup_rem2 <- Dat.l$RowID[Dat.l$agerange == "45-49"]
          print(paste("Removing 45-49 and keeping 40-49 in ", Dat.l$iso3[1], " (", Dat.l$Time[1], ")", sep = ""))
          ToRemove <- c(ToRemove, dup_rem2)
       }
      # Sometime the 15-19 was not abstracted but the 15-49 and the 20-24 and others were. We remove 15-49
         if (nrow(Dat.l) > 1 & any(Dat.l$agerange == "15-49") & any(Dat.l$agerange == "20-24")) {
          dup_rem3 <- Dat.l$RowID[Dat.l$agerange == "15-49"]
          print(paste("Removing 15-49 and keeping 20-24 and up in ", Dat.l$iso3[1], " (", Dat.l$Time[1], ")", sep = ""))
          ToRemove <- c(ToRemove, dup_rem3)
         }
      }

    if(length(ToRemove) > 1) {
      data <- data[-ToRemove, ] }

    return(data)
}

#' output age group by who region and high-income
age_dist_who <- function(IPV) {
    who <- addition_to_cnt("who") 
    
    IPV$age_agg <- ifelse(IPV$hiageI < 25, "<25",
               ifelse(IPV$loage >= 25 & IPV$hiageI < 35, "25-34",
               ifelse(IPV$loage >= 35 & IPV$hiageI < 50, "35-49",
               ifelse(IPV$loage >= 50 & IPV$hiageI < 65, "50-64",
               ifelse(IPV$loage >= 65, "65+", 
               ifelse(IPV$loage < 35 & IPV$hiageI < 50, "other (less 50)",
               ifelse(IPV$loage >= 50 & IPV$hiageI >= 65, "other (above 50)",
                "any other")))))))
    
    IPV$age_agg <- factor(IPV$age_agg, 
                          levels = c("<25", "25-34", "35-49" ,"50-64", "65+",
                                     "other (less 50)", "other (above 50)", "any other"))
    
    IPV$WHO <- who$region_name[match(IPV$iso3, who$iso3)]
    IPV$WHO_HI <- who$region_name_hi[match(IPV$iso3, who$iso3)]
    
    age_dist <- NULL
    ureg <- unique(IPV$WHO_HI)
    uage <- levels(IPV$age_agg)
    for (i in 1:length(ureg)) {
      ipv_i <- IPV[IPV$WHO_HI == ureg[i], ]
      for (j in 1:length(uage)) {
        ipv_ij <- ipv_i[ipv_i$age_agg == uage[j], ]
        age_dist_ij <- data.frame(region = ureg[i], age_grp = uage[j], n = nrow(ipv_ij), proportion = round(nrow(ipv_ij) / nrow(ipv_i) * 100))
        age_dist <- rbind(age_dist, age_dist_ij)
        }
       #overall
       if (i == length(ureg)) {
       for (j in 1:length(uage)) {
          ipv_ij <- IPV[IPV$age_agg == uage[j], ]
          age_dist_ij <- data.frame(region = "Overall", age_grp = uage[j], n = nrow(ipv_ij), proportion = round(nrow(ipv_ij) / nrow(IPV) * 100))
          age_dist <- rbind(age_dist, age_dist_ij)
        }}
    }
    return(age_dist)
}

n_studies_who <- function(IPV) {
  who <- addition_to_cnt("who") 
  
  IPV$WHO <- who$region_name[match(IPV$iso3, who$iso3)]
  IPV$WHO_HI <- who$region_name_hi[match(IPV$iso3, who$iso3)]
  
  n_dist <- data.frame(region = "Overall", n = length(unique(IPV$ID)), proportion = 100)
  ureg <- unique(IPV$WHO_HI)
  for (i in 1:length(ureg)) {
    ipv_i <- IPV[IPV$WHO_HI == ureg[i], ]
    n_dist_i <- data.frame(region = ureg[i], n = length(unique(ipv_i$ID)), proportion = length(unique(ipv_i$ID)) / length(unique(IPV$ID)) * 100)
    n_dist <- rbind(n_dist, n_dist_i)
  }
  return(n_dist)
}

#' Process MCMC
process_mcmc <- function(MCMC, report = "ever") {
  
  if (!any(names(MCMC) == "a_g_e")) {
  # converting to mcmc list
  a_g <- coda::as.mcmc.list(MCMC$a_g)
  a_i <- coda::as.mcmc.list(MCMC$a_i)
  a_c <- coda::as.mcmc.list(MCMC$a_c)
  a_r <- coda::as.mcmc.list(MCMC$a_r)
  a_sr <- coda::as.mcmc.list(MCMC$a_sr) 
  b_spl <- coda::as.mcmc.list(MCMC$b_spl)
  b_spl_sr <- coda::as.mcmc.list(MCMC$b_spl_sr)
  b_spl_r <- coda::as.mcmc.list(MCMC$b_spl_r)
  b_spl_c <- coda::as.mcmc.list(MCMC$b_spl_c)
  b_t <- coda::as.mcmc.list(MCMC$b_t) 
  b_t_sr <- coda::as.mcmc.list(MCMC$b_t_sr)
  b_t_r <- coda::as.mcmc.list(MCMC$b_t_r) 
  b_t_c <- coda::as.mcmc.list(MCMC$b_t_c) 

  sd_ai <- coda::as.mcmc.list(MCMC$sd_ai)
  sd_ai2 <- coda::as.mcmc.list(MCMC$sd_ai_NotNat)
  sd_ac <- coda::as.mcmc.list(MCMC$sd_ac)
  sd_ar <- coda::as.mcmc.list(MCMC$sd_ar)
  sd_sr <- coda::as.mcmc.list(MCMC$sd_sr)
  sd_spl_sr <- coda::as.mcmc.list(MCMC$sd_spl_sr)
  sd_spl_r <- coda::as.mcmc.list(MCMC$sd_spl_r)
  sd_spl_c <- coda::as.mcmc.list(MCMC$sd_spl_c)
  sd_t_sr <- coda::as.mcmc.list(MCMC$sd_t_sr) 
  sd_t_r <- coda::as.mcmc.list(MCMC$sd_t_r) 
  sd_t_c <- coda::as.mcmc.list(MCMC$sd_t_c)
 
  linpred <- coda::as.mcmc.list(MCMC$linpred)
 }

  if (any(names(MCMC) == "a_g_e") & report == "ever") {
  # converting to mcmc list
  a_g <- coda::as.mcmc.list(MCMC$a_g_e)
  a_i <- coda::as.mcmc.list(MCMC$a_i_e)
  a_c <- coda::as.mcmc.list(MCMC$a_c_e)
  a_r <- coda::as.mcmc.list(MCMC$a_r_e)
  a_sr <- coda::as.mcmc.list(MCMC$a_sr_e) 
  b_spl <- coda::as.mcmc.list(MCMC$b_spl_e)
  b_spl_sr <- coda::as.mcmc.list(MCMC$b_spl_sr_e)
  b_spl_r <- coda::as.mcmc.list(MCMC$b_spl_r_e)
  b_spl_c <- coda::as.mcmc.list(MCMC$b_spl_c_e)
  b_t <- coda::as.mcmc.list(MCMC$b_t_e) 
  b_t_sr <- coda::as.mcmc.list(MCMC$b_t_sr_e)
  b_t_r <- coda::as.mcmc.list(MCMC$b_t_r_e)
  b_t_c <- coda::as.mcmc.list(MCMC$b_t_c_e) 

  sd_ai <- coda::as.mcmc.list(MCMC$sd_ai_e)
  sd_ai2 <- coda::as.mcmc.list(MCMC$sd_ai_NotNat_e)
  sd_ac <- coda::as.mcmc.list(MCMC$sd_ac_e)
  sd_ar <- coda::as.mcmc.list(MCMC$sd_ar_e)
  sd_sr <- coda::as.mcmc.list(MCMC$sd_sr_e)
  sd_spl_sr <- coda::as.mcmc.list(MCMC$sd_spl_sr_e)
  sd_spl_r <- coda::as.mcmc.list(MCMC$sd_spl_r_e)
  sd_spl_c <- coda::as.mcmc.list(MCMC$sd_spl_c_e)
  sd_t_sr <- coda::as.mcmc.list(MCMC$sd_t_sr_e)  
  sd_t_r <- coda::as.mcmc.list(MCMC$sd_t_r_e)  
  sd_t_c <- coda::as.mcmc.list(MCMC$sd_t_c_e)

  linpred <- coda::as.mcmc.list(MCMC$linpred_e)
  }
  
  if (any(names(MCMC) == "a_g_p") & report == "past") {
  # Converting to mcmc list
  a_g <- coda::as.mcmc.list(MCMC$a_g_p)
  a_i <- coda::as.mcmc.list(MCMC$a_i_p)
  a_c <- coda::as.mcmc.list(MCMC$a_c_p)
  a_r <- coda::as.mcmc.list(MCMC$a_r_p)
  a_sr <- coda::as.mcmc.list(MCMC$a_sr_p) 
  b_spl <- coda::as.mcmc.list(MCMC$b_spl_p)
  b_spl_sr <- coda::as.mcmc.list(MCMC$b_spl_sr_p)
  b_spl_r <- coda::as.mcmc.list(MCMC$b_spl_r_p)
  b_spl_c <- coda::as.mcmc.list(MCMC$b_spl_c_p)
  b_t <- coda::as.mcmc.list(MCMC$b_t_p) 
  b_t_sr <- coda::as.mcmc.list(MCMC$b_t_sr_p)
  b_t_r <- coda::as.mcmc.list(MCMC$b_t_r_p)
  b_t_c  <- coda::as.mcmc.list(MCMC$b_t_c_p) 

  sd_ai <- coda::as.mcmc.list(MCMC$sd_ai_p)
  sd_ai2 <- coda::as.mcmc.list(MCMC$sd_ai_NotNat_p)
  sd_ac <- coda::as.mcmc.list(MCMC$sd_ac_p)
  sd_ar <- coda::as.mcmc.list(MCMC$sd_ar_p)
  sd_sr <- coda::as.mcmc.list(MCMC$sd_sr_p)
  sd_spl_sr <- coda::as.mcmc.list(MCMC$sd_spl_sr_p)
  sd_spl_r <- coda::as.mcmc.list(MCMC$sd_spl_r_p)
  sd_spl_c <- coda::as.mcmc.list(MCMC$sd_spl_c_p)
  sd_t_sr <- coda::as.mcmc.list(MCMC$sd_t_sr_p)  
  sd_t_r <- coda::as.mcmc.list(MCMC$sd_t_r_p)  
  sd_t_c <- coda::as.mcmc.list(MCMC$sd_t_c_p)

  linpred <- coda::as.mcmc.list(MCMC$linpred_p)
  }
  
  # Combining all chains into one vector
  glb <- do.call(rbind, a_g)
  stu <- do.call(rbind, a_i)
  cnt <- do.call(rbind, a_c)
  reg <- do.call(rbind, a_r)
  sup <- do.call(rbind, a_sr)
  
  b_spline <- do.call(rbind, b_spl)
  b_spline_sup <- do.call(rbind, b_spl_sr)
  b_spline_reg <- do.call(rbind, b_spl_r)
  b_spline_cnt <- do.call(rbind, b_spl_c)
  b_time <- do.call(rbind, b_t)
  b_time_sup <- do.call(rbind, b_t_sr)
  b_time_reg <- do.call(rbind, b_t_r)
  b_time_cnt <- do.call(rbind, b_t_c)  
  
  sd_sup <- do.call(rbind, sd_sr)
  sd_reg <- do.call(rbind, sd_ar)
  sd_cnt <- do.call(rbind, sd_ac)
  sd_spline_sup <- do.call(rbind, sd_spl_sr) 
  sd_spline_reg <- do.call(rbind, sd_spl_r)
  sd_spline_cnt <- do.call(rbind, sd_spl_c)
  sd_time_sup <- do.call(rbind, sd_t_sr)
  sd_time_reg <- do.call(rbind, sd_t_r)
  sd_time_cnt <- do.call(rbind, sd_t_c)

  lprd <- do.call(rbind, linpred)

  val <- list(
              # original parameters
              a_g = a_g, a_i = a_i, a_c = a_c, a_r = a_r, a_sr = a_sr,
              b_spl = b_spl, b_spl_sr = b_spl_sr, b_spl_r = b_spl_r,  b_spl_c = b_spl_c,
              b_t = b_t, b_t_sr = b_t_sr, b_t_r = b_t_r, b_t_c = b_t_c,
              sd_ai = sd_ai, sd_ai2 = sd_ai2, sd_ac = sd_ac, sd_ar = sd_ar, sd_sr = sd_sr, 
              sd_spl_sr = sd_spl_sr, sd_spl_r = sd_spl_r, sd_spl_c = sd_spl_c,
              sd_t_sr = sd_t_sr, sd_t_r = sd_t_r, sd_t_c = sd_t_c, 
              linpred = linpred,
              # formatted parameters
              glb = glb, stu = stu, cnt = cnt, reg = reg, sup = sup,
              b_spline = b_spline, b_spline_sup = b_spline_sup, b_spline_reg = b_spline_reg, b_spline_cnt = b_spline_cnt,
              b_time = b_time, b_time_sup = b_time_sup, b_time_reg = b_time_reg, b_time_cnt = b_time_cnt,
              sd_spline_sup = sd_spline_sup, sd_spline_reg = sd_spline_reg, sd_spline_cnt = sd_spline_cnt,
              sd_time_sup = sd_time_sup, sd_time_reg = sd_time_reg, sd_time_cnt = sd_time_cnt, 
              sd_sup = sd_sup, sd_reg = sd_reg, sd_cnt = sd_cnt, 
              lprd = lprd) 
  return(val)
}


#' Process MCMC
process_mcmc_npsv <- function(MCMC) {
  
  if (any(names(MCMC) == "b_t_r")) {
  # converting to mcmc list
  a_g <- coda::as.mcmc.list(MCMC$a_g)
  a_i <- coda::as.mcmc.list(MCMC$a_i)
  a_c <- coda::as.mcmc.list(MCMC$a_c)
  a_r <- coda::as.mcmc.list(MCMC$a_r)
  a_sr <- coda::as.mcmc.list(MCMC$a_sr) 
  b_spl <- coda::as.mcmc.list(MCMC$b_spl)
  b_spl_sr <- coda::as.mcmc.list(MCMC$b_spl_sr)
  b_spl_r <- coda::as.mcmc.list(MCMC$b_spl_r)
  b_spl_c <- coda::as.mcmc.list(MCMC$b_spl_c)
  b_t <- coda::as.mcmc.list(MCMC$b_t) 
  b_t_sr <- coda::as.mcmc.list(MCMC$b_t_sr)
  b_t_r <- coda::as.mcmc.list(MCMC$b_t_r) 
  b_t_c <- coda::as.mcmc.list(MCMC$b_t_c) 

  sd_ai <- coda::as.mcmc.list(MCMC$sd_ai)
  sd_ai2 <- coda::as.mcmc.list(MCMC$sd_ai_NotNat)
  sd_ac <- coda::as.mcmc.list(MCMC$sd_ac)
  sd_ar <- coda::as.mcmc.list(MCMC$sd_ar)
  sd_sr <- coda::as.mcmc.list(MCMC$sd_sr)
  sd_spl_sr <- coda::as.mcmc.list(MCMC$sd_spl_sr)
  sd_spl_r <- coda::as.mcmc.list(MCMC$sd_spl_r)
  sd_spl_c <- coda::as.mcmc.list(MCMC$sd_spl_c)
  sd_t_sr <- coda::as.mcmc.list(MCMC$sd_t_sr) 
  sd_t_r <- coda::as.mcmc.list(MCMC$sd_t_r) 
  sd_t_c <- coda::as.mcmc.list(MCMC$sd_t_c)
 
  linpred <- coda::as.mcmc.list(MCMC$linpred)
  
  # Combining all chains into one vector
  glb <- do.call(rbind, a_g)
  stu <- do.call(rbind, a_i)
  cnt <- do.call(rbind, a_c)
  reg <- do.call(rbind, a_r)
  sup <- do.call(rbind, a_sr)
  
  b_spline <- do.call(rbind, b_spl)
  b_spline_sup <- do.call(rbind, b_spl_sr)
  b_spline_reg <- do.call(rbind, b_spl_r)
  b_spline_cnt <- do.call(rbind, b_spl_c)
  b_time <- do.call(rbind, b_t)
  b_time_sup <- do.call(rbind, b_t_sr)
  b_time_reg <- do.call(rbind, b_t_r)
  b_time_cnt <- do.call(rbind, b_t_c)  
  
  sd_sup <- do.call(rbind, sd_sr)
  sd_reg <- do.call(rbind, sd_ar)
  sd_cnt <- do.call(rbind, sd_ac)
  sd_spline_sup <- do.call(rbind, sd_spl_sr) 
  sd_spline_reg <- do.call(rbind, sd_spl_r)
  sd_spline_cnt <- do.call(rbind, sd_spl_c)
  sd_time_sup <- do.call(rbind, sd_t_sr)
  sd_time_reg <- do.call(rbind, sd_t_r)
  sd_time_cnt <- do.call(rbind, sd_t_c)

  lprd <- do.call(rbind, linpred)

  val <- list(
              # original parameters
              a_g = a_g, a_i = a_i, a_c = a_c, a_r = a_r, a_sr = a_sr,
              b_spl = b_spl, b_spl_sr = b_spl_sr, b_spl_r = b_spl_r,  b_spl_c = b_spl_c,
              b_t = b_t, b_t_sr = b_t_sr, b_t_r = b_t_r, b_t_c = b_t_c,
              sd_ai = sd_ai, sd_ai2 = sd_ai2, sd_ac = sd_ac, sd_ar = sd_ar, sd_sr = sd_sr, 
              sd_spl_sr = sd_spl_sr, sd_spl_r = sd_spl_r, sd_spl_c = sd_spl_c,
              sd_t_sr = sd_t_sr, sd_t_r = sd_t_r, sd_t_c = sd_t_c, 
              linpred = linpred,
              # formatted parameters
              glb = glb, stu = stu, cnt = cnt, reg = reg, sup = sup,
              b_spline = b_spline, b_spline_sup = b_spline_sup, b_spline_reg = b_spline_reg, b_spline_cnt = b_spline_cnt,
              b_time = b_time, b_time_sup = b_time_sup, b_time_reg = b_time_reg, b_time_cnt = b_time_cnt,
              sd_spline_sup = sd_spline_sup, sd_spline_reg = sd_spline_reg, sd_spline_cnt = sd_spline_cnt,
              sd_time_sup = sd_time_sup, sd_time_reg = sd_time_reg, sd_time_cnt = sd_time_cnt, 
              sd_sup = sd_sup, sd_reg = sd_reg, sd_cnt = sd_cnt, 
              lprd = lprd) 
 }

  if (all(names(MCMC) != "b_t_r")) {
  # converting to mcmc list
  a_g <- coda::as.mcmc.list(MCMC$a_g)
  a_i <- coda::as.mcmc.list(MCMC$a_i)
  a_c <- coda::as.mcmc.list(MCMC$a_c)
  a_r <- coda::as.mcmc.list(MCMC$a_r)
  a_sr <- coda::as.mcmc.list(MCMC$a_sr) 
  b_spl <- coda::as.mcmc.list(MCMC$b_spl)
  b_spl_sr <- coda::as.mcmc.list(MCMC$b_spl_sr)
  b_spl_r <- coda::as.mcmc.list(MCMC$b_spl_r)
  b_t <- coda::as.mcmc.list(MCMC$b_t) 
  b_t_sr <- coda::as.mcmc.list(MCMC$b_t_sr)


  sd_ai <- coda::as.mcmc.list(MCMC$sd_ai)
  sd_ai2 <- coda::as.mcmc.list(MCMC$sd_ai_NotNat)
  sd_ac <- coda::as.mcmc.list(MCMC$sd_ac)
  sd_ar <- coda::as.mcmc.list(MCMC$sd_ar)
  sd_sr <- coda::as.mcmc.list(MCMC$sd_sr)
  sd_spl_sr <- coda::as.mcmc.list(MCMC$sd_spl_sr)
  sd_spl_r <- coda::as.mcmc.list(MCMC$sd_spl_r)
  sd_t_sr <- coda::as.mcmc.list(MCMC$sd_t_sr) 

  linpred <- coda::as.mcmc.list(MCMC$linpred)
  
  # Combining all chains into one vector
  glb <- do.call(rbind, a_g)
  stu <- do.call(rbind, a_i)
  cnt <- do.call(rbind, a_c)
  reg <- do.call(rbind, a_r)
  sup <- do.call(rbind, a_sr)
  
  b_spline <- do.call(rbind, b_spl)
  b_spline_sup <- do.call(rbind, b_spl_sr)
  b_spline_reg <- do.call(rbind, b_spl_r)
  b_time <- do.call(rbind, b_t)
  b_time_sup <- do.call(rbind, b_t_sr)

  sd_sup <- do.call(rbind, sd_sr)
  sd_reg <- do.call(rbind, sd_ar)
  sd_cnt <- do.call(rbind, sd_ac)
  sd_spline_sup <- do.call(rbind, sd_spl_sr) 
  sd_spline_reg <- do.call(rbind, sd_spl_r)
  sd_time_sup <- do.call(rbind, sd_t_sr)

  lprd <- do.call(rbind, linpred)

  val <- list(
              # original parameters
              a_g = a_g, a_i = a_i, a_c = a_c, a_r = a_r, a_sr = a_sr,
              b_spl = b_spl, b_spl_sr = b_spl_sr, b_spl_r = b_spl_r,
              b_t = b_t, b_t_sr = b_t_sr, 
              sd_ai = sd_ai, sd_ai2 = sd_ai2, sd_ac = sd_ac, sd_ar = sd_ar, sd_sr = sd_sr, 
              sd_spl_sr = sd_spl_sr, sd_spl_r = sd_spl_r, 
              sd_t_sr = sd_t_sr, 
              linpred = linpred,
              # formatted parameters
              glb = glb, stu = stu, cnt = cnt, reg = reg, sup = sup,
              b_spline = b_spline, b_spline_sup = b_spline_sup, b_spline_reg = b_spline_reg, 
              b_time = b_time, b_time_sup = b_time_sup,
              sd_spline_sup = sd_spline_sup, sd_spline_reg = sd_spline_reg, 
              sd_time_sup = sd_time_sup, 
              sd_sup = sd_sup, sd_reg = sd_reg, sd_cnt = sd_cnt, 
              lprd = lprd) 
    }

 

  return(val)
}

#' Process lhs. Combining all samples from the MI
process_lhs <- function(res_lhs) {
  n_mi <- length(res_lhs)
  res <- res_lhs[[1]]
  if (n_mi >= 2) {
  for (i in 2:n_mi) {
    res$a_g <- as.mcmc(c(res$a_g, res_lhs[[i]]$a_g))
    res$a_i <- as.mcmc(c(res$a_i, res_lhs[[i]]$a_i))
    res$a_c <- as.mcmc(c(res$a_c, res_lhs[[i]]$a_c))
    res$a_r <- as.mcmc(c(res$a_r, res_lhs[[i]]$a_r))
    res$a_sr <- as.mcmc(c(res$a_sr, res_lhs[[i]]$a_sr))
    res$b_spl <- as.mcmc(c(res$b_spl, res_lhs[[i]]$b_spl))
    res$b_spl_sr <- as.mcmc(c(res$b_spl_sr, res_lhs[[i]]$b_spl_sr))
    res$b_spl_r <- as.mcmc(c(res$b_spl_r, res_lhs[[i]]$b_spl_r))
    res$b_spl_c <- as.mcmc(c(res$b_spl_c, res_lhs[[i]]$b_spl_c))
    res$b_t <- as.mcmc(c(res$b_t, res_lhs[[i]]$b_t))    
    res$b_t_sr <- as.mcmc(c(res$b_t_sr, res_lhs[[i]]$b_t_sr))    
    res$b_t_r <- as.mcmc(c(res$b_t_r, res_lhs[[i]]$b_t_r))    
    res$b_t_c <- as.mcmc(c(res$b_t_c, res_lhs[[i]]$b_t_c))  
    
    res$sd_ai <- as.mcmc(c(res$sd_ai, res_lhs[[i]]$sd_ai))
    res$sd_ai2 <- as.mcmc(c(res$sd_ai2, res_lhs[[i]]$sd_ai2))
    res$sd_ac <- as.mcmc(c(res$sd_ac, res_lhs[[i]]$sd_ac))
    res$sd_ar <- as.mcmc(c(res$sd_ar, res_lhs[[i]]$sd_ar))
    res$sd_sr <- as.mcmc(c(res$sd_sr, res_lhs[[i]]$sd_sr))
    res$sd_spl_sr <- as.mcmc(c(res$sd_spl_sr, res_lhs[[i]]$sd_spl_sr))
    res$sd_spl_r <- as.mcmc(c(res$sd_spl_r, res_lhs[[i]]$sd_spl_r))
    res$sd_spl_c <- as.mcmc(c(res$sd_spl_c, res_lhs[[i]]$sd_spl_c))
    res$sd_t_sr <- as.mcmc(c(res$sd_t_sr, res_lhs[[i]]$sd_t_sr))   
    res$sd_t_r <- as.mcmc(c(res$sd_t_r, res_lhs[[i]]$sd_t_r))   
    res$sd_t_c <- as.mcmc(c(res$sd_t_c, res_lhs[[i]]$sd_t_c))
 
    res$linpred <- as.mcmc(c(res$linpred, res_lhs[[i]]$linpred))
    
    res$glb <- rbind(res$glb, res_lhs[[i]]$glb)
    res$stu <- rbind(res$stu, res_lhs[[i]]$stu)
    res$cnt <- rbind(res$cnt, res_lhs[[i]]$cnt)
    res$reg <- rbind(res$reg, res_lhs[[i]]$reg)
    res$sup <- rbind(res$sup, res_lhs[[i]]$sup) 
    
    res$b_spline <- rbind(res$b_spline, res_lhs[[i]]$b_spline)
    res$b_spline_sup <- rbind(res$b_spline_sup, res_lhs[[i]]$b_spline_sup)
    res$b_spline_reg <- rbind(res$b_spline_reg, res_lhs[[i]]$b_spline_reg)
    res$b_spline_cnt <- rbind(res$b_spline_cnt, res_lhs[[i]]$b_spline_cnt)
    res$b_time <- rbind(res$b_time, res_lhs[[i]]$b_time)
    res$b_time_sup <- rbind(res$b_time_sup, res_lhs[[i]]$b_time_sup)
    res$b_time_reg <- rbind(res$b_time_reg, res_lhs[[i]]$b_time_reg)
    res$b_time_cnt <- rbind(res$b_time_cnt, res_lhs[[i]]$b_time_cnt)  
    
    res$sd_sup <- rbind(res$sd_sup, res_lhs[[i]]$sd_sup)
    res$sd_reg <- rbind(res$sd_reg, res_lhs[[i]]$sd_reg)
    res$sd_cnt <- rbind(res$sd_cnt, res_lhs[[i]]$sd_cnt)
    res$sd_spline_sup <- rbind(res$sd_spline_sup, res_lhs[[i]]$sd_spline_sup)
    res$sd_spline_reg <- rbind(res$sd_spline_reg, res_lhs[[i]]$sd_spline_reg)
    res$sd_spline_cnt <- rbind(res$sd_spline_cnt, res_lhs[[i]]$sd_spline_cnt)
    res$sd_time_sup <- rbind(res$sd_time_sup, res_lhs[[i]]$sd_time_sup)
    res$sd_time_reg <- rbind(res$sd_time_reg, res_lhs[[i]]$sd_time_reg)
    res$sd_time_cnt <- rbind(res$sd_time_cnt, res_lhs[[i]]$sd_time_cnt)

    res$lprd <- rbind(res$lprd, res_lhs[[i]]$lprd) 
  }}
  return(res)
} 


#' Posterior predictive checks
post_pred_chk <- function(res, SDat, IPV, outcome, suffix = "", covstats = TRUE) {
  if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
    stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }
 
  # Pre-processing step
  Nb.iter <- dim(res[["glb"]])[1]
  y_ppc <- NULL
  pb <- txtProgressBar(1, SDat$N, style = 3)
  for (i in 1:SDat$N) {
    p.i <- rbinom(Nb.iter, SDat$Denom[i], res[["lprd"]][, i]) / SDat$Denom[i]
    y_ppc.i <- quantile(p.i, probs = c(0.5, 0.025, 0.975))
    y_ppc <- rbind(y_ppc, y_ppc.i)
    setTxtProgressBar(pb, i)
  }
  close(pb)

  Prd <- data.frame(Median = y_ppc[, 1], LCI = y_ppc[, 2], UCI = y_ppc[, 3])
  Prd$RowID <- seq(1, nrow(Prd), 1)

  # Color Specification
  if (outcome == "Ever IPV") { a <- 6; b <- 4 }
  if (outcome == "Past Year IPV") { a <- 6; b <- 3 }
  if (outcome == "NPSV") { a <- 6; b <- 2 }
  Col <- c(rgb(20, 180, 235, 255, max = 255), rgb(134, 179, 0, 255, max = 255), rgb(204, 0, 67, 255, max = 255),
          rgb(255, 173, 51, 255, max = 255), rgb(255, 102, 53, 255, max = 255), rgb(135, 132, 148, 255, max = 255))
  ColT <- c(rgb(20, 180, 235, 150, max = 255), rgb(134, 179, 0, 150, max = 255), rgb(204, 0, 67, 150, max = 255),
          rgb(255, 173, 51, 150, max = 255), rgb(255, 102, 53, 150, max = 255), rgb(135, 132, 148, 150, max = 255))

  if (!is.null(dev.list())) { dev.off() }
  pdf(file = paste("Posterior Checks - ", outcome, " ", suffix, ".pdf", sep = ""),
      width = 10, height = 7.5)
  UR <- unique(IPV$reg)
  if (!is.null(IPV$severe)) {
  IPV$Severe <- as.character(ifelse(IPV$severe == "Only severe violence",
                                    "Severe", "All Severity"))
  } else { IPV$Severe <- NULL }
  for (i in 1:SDat$R) {
    Sel <- which(IPV$reg == UR[i])
    IPV.i <- IPV[Sel, ]
    Prd.i <- Prd[Sel, ]
    XLab.i <- IPV.i[!duplicated(IPV.i$iso3), ]
    xlim.i <- c(min(Prd.i$RowID), 
                ifelse((max(Prd.i$RowID) - min(Prd.i$RowID)) < 10, 
                       min(Prd.i$RowID) + 15.5, max(Prd.i$RowID) + 0.5))
    par(mar = c(3, 4, 2, 0))
    plot(Prd.i$Median ~ Prd.i$RowID,
      type = "n", ylim = c(-0.30, 1.25), xlim = xlim.i,
      xlab = "", ylab = "Prevalence", main = paste(outcome, "-", IPV.i$region[1], sep = " "),
      axes = FALSE)
    axis(1, at = c(min(Prd.i$RowID), max(Prd.i$RowID) + 1), labels = FALSE)
    axis(1, at = XLab.i$RowID, labels = NA, las = 2, cex.axis = 0.3)
    text(x = XLab.i$RowID, par("usr")[3] - 0.05, labels = XLab.i$iso3, 
         srt = 45, cex = 0.6, pos = 1, xpd = TRUE)
    axis(2, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2))
    points(IPV.i$Prv / 100 ~ IPV.i$RowID, pch = 17, cex = 0.5, col = Col[a])
    points(I(Prd.i$Median) ~ I(Prd.i$RowID + 0.25), pch = 16, cex = 0.5, col = Col[b])
    for (j in 1:dim(Prd.i)[1]) {
      segments(x0 = IPV.i$RowID[j], y0 = IPV.i$LCI[j], x1 = IPV.i$RowID[j], y1 = IPV.i$UCI[j], col = ColT[a])
      segments(x0 = Prd.i$RowID[j] + 0.25, y0 = Prd.i$LCI[j], x1 = Prd.i$RowID[j] + 0.25, y1 = Prd.i$UCI[j], col = ColT[b])
    }
    legend(x = min(IPV.i$RowID), y = 1.25, pch = c(17, 16), col = c(Col[a], Col[b]), lwd = 2,
      legend = c("Data", "Modeled Estimates"), bty = "n")
   # above graph 
   IPV.i$VioType <- paste(IPV.i$ID, IPV.i$startyr, IPV.i$violence)
      IPVTxtVioT <- IPV.i[!duplicated(IPV.i$VioType), ]
      text(IPVTxtVioT$RowID, y = 0.95, labels = as.character(IPVTxtVioT$violence), srt = 90, cex = 0.35)
  IPV.i$StudyPstat <- paste(IPV.i$ID, IPV.i$startyr, IPV.i$violence, IPV.i$Severe, IPV.i$pstat, sep = "-")
      IPV_pstat <- IPV.i[!duplicated(IPV.i$StudyPstat), ]
      text(IPV_pstat$RowID, y = 0.75, labels = IPV_pstat$pstat, srt = 90, cex = 0.35)
  # below graph    
  text(IPV.i$RowID, y = -0.05, labels = paste(IPV.i$loage, IPV.i$hiage, sep = "-"), srt = 90, cex = 0.3)      
  IPV.i$StudyYr <- paste(IPV.i$ID, IPV.i$startyr, sep = "-")
      IPV_yr <- IPV.i[!duplicated(IPV.i$StudyYr), ]
      text(IPV_yr$RowID, y = -0.10, labels = IPV_yr$startyr, srt = 45, cex = 0.35)
  IPV.i$GEO <- paste(IPV.i$ID, IPV.i$startyr, IPV.i$violence, IPV.i$Severe, IPV.i$pstat, IPV.i$geo, sep = "-")
      IPV_GEO <- IPV.i[!duplicated(IPV.i$GEO), ]
      text(IPV_GEO$RowID, y = -0.16, labels = IPV_GEO$geo, srt = 90, cex = 0.35)
  IPV.i$StudySvr <- paste(IPV.i$ID, IPV.i$startyr, IPV.i$violence, IPV.i$Severe, sep = "-")
      IPV_Severe <- IPV.i[!duplicated(IPV.i$StudySvr), ]
      text(IPV_Severe$RowID, y = -0.25, labels = IPV_Severe$Severe, srt = 90, cex = 0.35)
  # This line is for diplaying the time
  if (outcome == "NPSV") {
    IPV.i$PastYr <- paste(IPV.i$ID, IPV.i$startyr, IPV.i$violence, IPV.i$Severe, IPV.i$pstat, IPV.i$geo, IPV.i$viotime, sep = "-")
      IPV_PastYr <- IPV.i[!duplicated(IPV.i$PastYr), ]
      text(IPV_PastYr$RowID, y = -0.20, labels = IPV_PastYr$viotime, srt = 45, cex = 0.35)  
  } }
dev.off()
    
    if (covstats == TRUE) {
      # overall
      me <- median(IPV$Prv/100 - Prd$Median)
      mae <- median(abs(IPV$Prv/100 - Prd$Median ))
      mre <- median((IPV$Prv/100 - Prd$Median) / (IPV$Prv/100 + 1e-100))
      mare <- median(abs((IPV$Prv/100 - Prd$Median)) / (IPV$Prv/100 + 1e-100))
      stats <- ifelse(IPV$Prv/100 >= Prd$LCI & IPV$Prv/100 <= Prd$UCI, 1, 0)
      li <- ifelse(IPV$Prv/100 < Prd$LCI, 1, 0)
      ui <- ifelse(IPV$Prv/100 > Prd$UCI, 1, 0)
      stats_overall <- data.frame(region = "Overall", n = nrow(IPV), 
                                  me = me, mae = mae, mre = mre, mare = mare,
                                  cov = mean(stats), lower = mean(li), upper = mean(ui))
      # by regions
      UR <- unique(IPV$region)
      stats_r <- NULL
      for (i in 1:length(UR)) {
            Sel <- which(IPV$region == UR[i])
            IPV.i <- IPV[Sel, ]
            Prd.i <- Prd[Sel, ]
            me.i <- median(IPV.i$Prv / 100 - Prd.i$Median)
            mae.i <- median(abs(IPV.i$Prv / 100 - Prd.i$Median ))
            mre.i <- median((IPV.i$Prv / 100 - Prd.i$Median) / (IPV.i$Prv / 100), na.rm = TRUE)
            mare.i <- median(abs((IPV.i$Prv / 100 - Prd.i$Median)) / (IPV.i$Prv / 100), na.rm = TRUE)
            stats.i <- ifelse(IPV.i$Prv / 100 >= Prd.i$LCI & IPV.i$Prv / 100 <= Prd.i$UCI, 1, 0)
            li.i <- ifelse(IPV.i$Prv / 100 < Prd.i$LCI, 1, 0)
            ui.i <- ifelse(IPV.i$Prv / 100 > Prd.i$UCI, 1, 0)
            stats_region <- data.frame(region = UR[i], n = nrow(IPV.i), 
                                       me = me.i, mae = mae.i, mre = mre.i, mare = mare.i,
                                       lower = mean(li.i), upper = mean(ui.i), cov = mean(stats.i))
            stats_r <- rbind(stats_r, stats_region)
      }
      
     # by time period
      time_per <- c(2005, 2010, 2015, 2020)
      stats_t <- NULL
      for (i in 1:length(time_per)) {
            name_time_per <- paste(time_per[i] - 5, "-", time_per[i] - 1, sep = "")
            Sel <- which(IPV$Time >= (time_per[i] - 5) & IPV$Time < time_per[i])
            IPV.i <- IPV[Sel, ]
            Prd.i <- Prd[Sel, ]
            me.i <- median(IPV.i$Prv / 100 - Prd.i$Median)
            mae.i <- median(abs(IPV.i$Prv / 100 - Prd.i$Median ))
            mre.i <- median((IPV.i$Prv / 100 - Prd.i$Median) / (IPV.i$Prv / 100), na.rm = TRUE)
            mare.i <- median(abs((IPV.i$Prv / 100 - Prd.i$Median)) / (IPV.i$Prv / 100), na.rm = TRUE)
            stats.i <- ifelse(IPV.i$Prv / 100 >= Prd.i$LCI & IPV.i$Prv / 100 <= Prd.i$UCI, 1, 0)
            li.i <- ifelse(IPV.i$Prv / 100 < Prd.i$LCI, 1, 0)
            ui.i <- ifelse(IPV.i$Prv / 100 > Prd.i$UCI, 1, 0)
            stats_time <- data.frame(region = name_time_per, n = nrow(IPV.i), 
                                       me = me.i, mae = mae.i, mre = mre.i, mare = mare.i,
                                       lower = mean(li.i), upper = mean(ui.i), cov = mean(stats.i))
            stats_t <- rbind(stats_t, stats_time)
      }

      stats <- rbind(stats_r, stats_t, stats_overall)
      stats$me <- round(stats$me * 100, 1)
      stats$mae <- round(stats$mae * 100, 1)
      stats$mre <- round(stats$mre * 100, 1)                     
      stats$mare <- round(stats$mare * 100, 1)                       
      stats$cov <- round(stats$cov * 100, 1)
      stats$lower <- round(stats$lower * 100, 1)
      stats$upper <- round(stats$upper * 100, 1)
      write.csv(stats, paste("In-samples validation - ", outcome, " ", suffix, ".csv", sep = ""))
      return(stats)
    }
}

# widely applicable information criterion
waic <- function(pred, y, n) { 
    ll <- dbinom(x = y, size = n, p = pred, log = TRUE) 
    lik <- exp(ll)
    waic1 <- log(apply(lik, 1, mean))
    waic2 <- apply(ll, 1, sd)
    waic_pw <- -2 * (waic1 - waic2) # case-wise waic
    elpd_waic <- sum(waic1) - sum(waic2) # total waic
    val <- -2 * elpd_waic
    return(val)
}
  
# DIC for splines
splines_dic <- function(SDat_original, Center_a = 30, 
                        BoundKnots = c(15, 50) - Center_a,
                        modelstring, n_iter = 150000, n_adapt = 5000,
                        n_burn = 5000, n_thin = 10, n_chain = 3) {
  
  SDat_dic <- SDat_original
  param <- c('linpred')
  
  Knots <- 20 - Center_a 
  Spl <- ns(IPV$Age65 - Center_a, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Spl_W <- ns(seq(15, 65, by = 5)  - Center_a, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Ndof <- ncol(Spl)
  jags.fit <- rjags::jags.model(textConnection(modelstring), n.adapt = n_adapt, data = SDat_dic, n.chains = n_chain)
  update(jags.fit, n_burn)
  dic_at20 <- rjags::dic.samples(jags.fit, n_iter, nthin = n_thin, "popt") # Penalized expected deviance
  mcmc <- rjags::jags.samples(jags.fit, variable.names = param, n.iter = n_iter, thin = n_thin)
  lprd <- do.call(rbind, coda::as.mcmc.list(mcmc$linpred))
  waic_at20 <- waic(pred = lprd, y = SDat_original$Num, n = SDat_original$Denom)

  Knots <- 25 - Center_a 
  Spl <- ns(IPV$Age65 - Center_a, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Spl_W <- ns(seq(15, 65, by = 5)  - Center_a, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Ndof <- ncol(Spl)
  jags.fit <- rjags::jags.model(textConnection(modelstring), n.adapt = n_adapt, data = SDat_dic, n.chains = n_chain)
  update(jags.fit, n_burn)
  dic_at25 <- rjags::dic.samples(jags.fit, n_iter, nthin = n_thin, "popt") # Penalized expected deviance
  mcmc <- rjags::jags.samples(jags.fit, variable.names = param, n.iter = n_iter, thin = n_thin)
  lprd <- do.call(rbind, coda::as.mcmc.list(mcmc$linpred))
  waic_at25 <- waic(pred = lprd, y = SDat_original$Num, n = SDat_original$Denom)
  
  Knots <- 30 - Center_a 
  Spl <- ns(IPV$Age65 - Center_a, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Spl_W <- ns(seq(15, 65, by = 5)  - Center_a, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Ndof <- ncol(Spl)
  jags.fit <- rjags::jags.model(textConnection(modelstring), n.adapt = n_adapt, data = SDat_dic, n.chains = n_chain)
  update(jags.fit, n_burn)
  dic_at30 <- rjags::dic.samples(jags.fit, n_iter, nthin = n_thin, "popt") # Penalized expected deviance
  mcmc <- rjags::jags.samples(jags.fit, variable.names = param, n.iter = n_iter, thin = n_thin)
  lprd <- do.call(rbind, coda::as.mcmc.list(mcmc$linpred))
  waic_at30 <- waic(pred = lprd, y = SDat_original$Num, n = SDat_original$Denom)
  
  Knots <- 35 - Center_a 
  Spl <- ns(IPV$Age65 - Center_a, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Spl_W <- ns(seq(15, 65, by = 5) - Center_a, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Ndof <- ncol(Spl)
  jags.fit <- rjags::jags.model(textConnection(modelstring), n.adapt = n_adapt, data = SDat_dic, n.chains = n_chain)
  update(jags.fit, n_burn)
  dic_at35 <- rjags::dic.samples(jags.fit, n_iter, nthin = n_thin, "popt") # Penalized expected deviance
  mcmc <- rjags::jags.samples(jags.fit, variable.names = param, n.iter = n_iter, thin = n_thin)
  lprd <- do.call(rbind, coda::as.mcmc.list(mcmc$linpred))
  waic_at35 <- waic(pred = lprd, y = SDat_original$Num, n = SDat_original$Denom)
  
  Knots <- c(20, 35) - Center_a 
  Spl <- ns(IPV$Age65 - Center_a, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Spl_W <- ns(seq(15, 65, by = 5) - Center_a, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Ndof <- ncol(Spl)
  jags.fit <- rjags::jags.model(textConnection(modelstring), n.adapt = n_adapt, data = SDat_dic, n.chains = n_chain)
  update(jags.fit, n_burn)
  dic_at2035 <- rjags::dic.samples(jags.fit, n_iter, nthin = n_thin, "popt") # Penalized expected deviance
  mcmc <- rjags::jags.samples(jags.fit, variable.names = param, n.iter = n_iter, thin = n_thin)
  lprd <- do.call(rbind, coda::as.mcmc.list(mcmc$linpred))
  waic_at2035 <- waic(pred = lprd, y = SDat_original$Num, n = SDat_original$Denom)
  
  Knots <- c(20, 40) - Center_a 
  Spl <- ns(IPV$Age65 - Center_a, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Spl_W <- ns(seq(15, 65, by = 5) - Center_a, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Ndof <- ncol(Spl)
  jags.fit <- rjags::jags.model(textConnection(modelstring), n.adapt = n_adapt, data = SDat_dic, n.chains = n_chain)
  update(jags.fit, n_burn)
  dic_at2040 <- rjags::dic.samples(jags.fit, n_iter, nthin = n_thin, "popt") # Penalized expected deviance
  mcmc <- rjags::jags.samples(jags.fit, variable.names = param, n.iter = n_iter, thin = n_thin)
  lprd <- do.call(rbind, coda::as.mcmc.list(mcmc$linpred))
  waic_at2040 <- waic(pred = lprd, y = SDat_original$Num, n = SDat_original$Denom)
  
  Knots <- c(25, 35) - Center_a 
  Spl <- ns(IPV$Age65 - Center_a, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Spl_W <- ns(seq(15, 65, by = 5)  - Center_a, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Ndof <- ncol(Spl)
  jags.fit <- rjags::jags.model(textConnection(modelstring), n.adapt = n_adapt, data = SDat_dic, n.chains = n_chain)
  update(jags.fit, n_burn)
  dic_at2535 <- rjags::dic.samples(jags.fit, n_iter, nthin = n_thin, "popt") # Penalized expected deviance
  mcmc <- rjags::jags.samples(jags.fit, variable.names = param, n.iter = n_iter, thin = n_thin)
  lprd <- do.call(rbind, coda::as.mcmc.list(mcmc$linpred))
  waic_at2535 <- waic(pred = lprd, y = SDat_original$Num, n = SDat_original$Denom)
  
  Knots <- c(25, 40) - Center_a 
  Spl <- ns(IPV$Age65 - Center_a, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Spl_W <- ns(seq(15, 65, by = 5)  - Center_a, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Ndof <- ncol(Spl)
  jags.fit <- rjags::jags.model(textConnection(modelstring), n.adapt = n_adapt, data = SDat_dic, n.chains = n_chain)
  update(jags.fit, n_burn)
  dic_at2540 <- rjags::dic.samples(jags.fit, n_iter, nthin = n_thin, "popt") # Penalized expected deviance
  mcmc <- rjags::jags.samples(jags.fit, variable.names = param, n.iter = n_iter, thin = n_thin)
  lprd <- do.call(rbind, coda::as.mcmc.list(mcmc$linpred))
  waic_at2540 <- waic(pred = lprd, y = SDat_original$Num, n = SDat_original$Denom)
  
  diclst <- list(dic_at20 = dic_at20,
                 dic_at25 = dic_at25, 
                 dic_at30 = dic_at30,
                 dic_at35 = dic_at35,
                 dic_at2035 = dic_at2035,
                 dic_at2040 = dic_at2040,
                 dic_at2535 = dic_at2535, 
                 dic_at2540 = dic_at2540)
  
  val <- data.frame(knots = c("one at 20 years",
                              "one at 25 years",
                               "one at 30 years",
                               "one at 35 years",
                               "two at 20 & 35 years",
                               "two at 20 & 40 years",
                               "two at 25 & 35 years",
                               "two at 25 & 40 years"),
                     waic = c(waic_at20,
                              waic_at25,
                              waic_at30,
                              waic_at35,
                              waic_at2035,
                              waic_at2040,
                              waic_at2535,
                              waic_at2540),
                     penalized_dev = c(abstract_pd(dic_at20),
                                       abstract_pd(dic_at25),
                                       abstract_pd(dic_at30),
                                       abstract_pd(dic_at35),
                                       abstract_pd(dic_at2035),  
                                       abstract_pd(dic_at2040),                                       
                                       abstract_pd(dic_at2535),
                                       abstract_pd(dic_at2540)))

    mindic <- which.min(val$penalized_dev)
    val$dif <- NA
    val$dif_se <- NA
    for (i in 1:(length(diclst))) {
      if (i != mindic) {
        val$dif[i]
        n1 <- names(diclst[[mindic]]$deviance)
        n2 <- names(diclst[[i]]$deviance)
        if (!identical(n1, n2)) {
        ord1 <- order(diclst[[mindic]])
        ord2 <- order(diclst[[i]])
        diclst[[i]]$deviance[ord1] <- diclst[[i]]$deviance[ord2]
        diclst[[i]]$penalty[ord1] <- diclst[[i]]$penalty[ord2]
        }
        delta <-  (sapply(diclst[[i]]$deviance, mean) + sapply(diclst[[i]]$penalty, mean)) -
                  (sapply(diclst[[mindic]]$deviance, mean) + sapply(diclst[[mindic]]$penalty, mean))
        
        val$dif[i] <- round(sum(delta), 2)
        val$dif_se[i] <- round(sqrt(length(delta)) * sd(delta), 2)
      }}
                     
  return(val)
}

# DIC for splines
splines_dic_time <- function(SDat_original, Center_t = 2018, 
                        BoundKnots = c(2000, 2018) - Center_t,
                        modelstring, n_iter = 75000, n_adapt = 5000,
                        n_burn = 5000, n_thin = 10, n_chain = 4) {
  
  SDat_dic <- SDat_original
  
  Knots <- 2008 - Center_t 
  SDat_dic$Spl_t <- ns(IPV$Time - Center_t, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Nt <- ncol(SDat_dic$Spl_t)
  jags.fit <- rjags::jags.model(textConnection(modelstring), n.adapt = n_adapt, data = SDat_dic, n.chains = n_chain)
  update(jags.fit, n_burn)
  dic_at08 <- rjags::dic.samples(jags.fit, n_iter, nthin = n_thin, "popt") # Penalized expected deviance

  Knots <- 2009 - Center_t 
  SDat_dic$Spl_t <- ns(IPV$Time - Center_t, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Nt <- ncol(SDat_dic$Spl_t)
  jags.fit <- rjags::jags.model(textConnection(modelstring), n.adapt = n_adapt, data = SDat_dic, n.chains = n_chain)
  update(jags.fit, n_burn)
  dic_at09 <- rjags::dic.samples(jags.fit, n_iter, nthin = n_thin, "popt") # Penalized expected deviance
  
  Knots <- 2010 - Center_t 
  SDat_dic$Spl_t <- ns(IPV$Time - Center_t, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Nt <- ncol(SDat_dic$Spl_t)
  jags.fit <- rjags::jags.model(textConnection(modelstring), n.adapt = n_adapt, data = SDat_dic, n.chains = n_chain)
  update(jags.fit, n_burn)
  dic_at10 <- rjags::dic.samples(jags.fit, n_iter, nthin = n_thin, "popt") # Penalized expected deviance
  
  Knots <- 2011 - Center_t 
  SDat_dic$Spl_t <- ns(IPV$Time - Center_t, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Nt <- ncol(SDat_dic$Spl_t)
  jags.fit <- rjags::jags.model(textConnection(modelstring), n.adapt = n_adapt, data = SDat_dic, n.chains = n_chain)
  update(jags.fit, n_burn)
  dic_at11 <- rjags::dic.samples(jags.fit, n_iter, nthin = n_thin, "popt") # Penalized expected deviance
  
  Knots <- 2012 - Center_t 
  SDat_dic$Spl_t <- ns(IPV$Time - Center_t, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Nt <- ncol(SDat_dic$Spl_t)
  jags.fit <- rjags::jags.model(textConnection(modelstring), n.adapt = n_adapt, data = SDat_dic, n.chains = n_chain)
  update(jags.fit, n_burn)
  dic_at12 <- rjags::dic.samples(jags.fit, n_iter, nthin = n_thin, "popt") # Penalized expected deviance
  
  Knots <- 2013 - Center_t 
  SDat_dic$Spl_t <- ns(IPV$Time - Center_t, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Nt <- ncol(SDat_dic$Spl_t)
  jags.fit <- rjags::jags.model(textConnection(modelstring), n.adapt = n_adapt, data = SDat_dic, n.chains = n_chain)
  update(jags.fit, n_burn)
  dic_at15 <- rjags::dic.samples(jags.fit, n_iter, nthin = n_thin, "popt") # Penalized expected deviance
  
  Knots <- c(2007.5, 2012.5) - Center_t 
  SDat_dic$Spl_t <- ns(IPV$Time - Center_t, knots = Knots, Boundary.knots = BoundKnots)
  SDat_dic$Nt <- ncol(SDat_dic$Spl_t)
  jags.fit <- rjags::jags.model(textConnection(modelstring), n.adapt = n_adapt, data = SDat_dic, n.chains = n_chain)
  update(jags.fit, n_burn)
  dic_at0712 <- rjags::dic.samples(jags.fit, n_iter, nthin = n_thin, "popt") # Penalized expected deviance
  
  diclst <- list(dic_at08 = dic_at08,
                 dic_at09 = dic_at09, 
                 dic_at10 = dic_at10,
                 dic_at11 = dic_at11,
                 dic_at12 = dic_at12,
                 dic_at13 = dic_at13,
                 dic_at0712 = dic_at0712)
  
  val <- data.frame(knots = c("one in 2005",
                              "one in 2007.5",
                               "one in 2010",
                               "one in 2012.5",
                               "one in 2015",
                               "two in 2007.5 and 2012.5"),
                    penalized_dev = c(abstract_pd(dic_at05),
                                      abstract_pd(dic_at07),
                                      abstract_pd(dic_at10),
                                      abstract_pd(dic_at12),
                                      abstract_pd(dic_at15),  
                                      abstract_pd(dic_at0712)))

    mindic <- which.min(val$penalized_dev)
    val$dif <- NA
    val$dif_se <- NA
    for (i in 1:(length(diclst))) {
      if (i != mindic) {
        val$dif[i]
        n1 <- names(diclst[[mindic]]$deviance)
        n2 <- names(diclst[[i]]$deviance)
        if (!identical(n1, n2)) {
        ord1 <- order(diclst[[mindic]])
        ord2 <- order(diclst[[i]])
        diclst[[i]]$deviance[ord1] <- diclst[[i]]$deviance[ord2]
        diclst[[i]]$penalty[ord1] <- diclst[[i]]$penalty[ord2]
        }
        delta <-  (sapply(diclst[[i]]$deviance, mean) + sapply(diclst[[i]]$penalty, mean)) -
                  (sapply(diclst[[mindic]]$deviance, mean) + sapply(diclst[[mindic]]$penalty, mean))
        
        val$dif[i] <- round(sum(delta), 2)
        val$dif_se[i] <- round(sqrt(length(delta)) * sd(delta), 2)
      }}
                     
  return(val)
}


#' Plotting the age pattern by regions
plot_age_pattern <- function(res, SDat, IPV, outcome, x_max = 65,
                             Center_a = 30, Center_t = 2018,
                             year_pr = 2018,
                             suffix = "") {
  if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
    stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }

    # We create the database to predict
  age <- seq(15, 85, 1)
  age65 <- ifelse(age >= 65, 65, age)
  Xage <- splines::ns(age65 - Center_a,
                      knots = attr(SDat$Spl_W, "knots"),
                      Boundary.knots = attr(SDat$Spl_W, "Boundary.knots"))
  
  # manipulating time trends 
  is_period <- length(attr(SDat$Spl_t, "period")) > 0 
  if (is_period) {
    XTime <- as.matrix(ifelse(year_pr < attr(SDat$Spl_t, "cutoff"), 0, 1))
  }
  if (!is_period & length(attr(SDat$Spl_t, "knots")) > 0) {
  XTime <- as.matrix(splines::ns(year_pr - Center_t,
                      knots = attr(SDat$Spl_t, "knots"),
                      Boundary.knots = attr(SDat$Spl_t, "Boundary.knots")))    
  } 
  if (!is_period & length(attr(SDat$Spl_t, "knots")) == 0) {
  XTime <- as.matrix(splines::ns(year_pr - Center_t,
                      df = ncol(SDat$Spl_t),
                      Boundary.knots = attr(SDat$Spl_t, "Boundary.knots")))      
  }
    
  # Regions
  XPred <- NULL
  UIR <- unique(IPV$reg)
  XPred <- NULL
  for (i in 1:SDat$R) {
    XPrd.i <- data.frame(Age = age, Spl_a = Xage, Spl_t = XTime)
    XPrd.i$region <- as.character(UIR[i])
    XPrd.i$Region <- IPV$Region[which(IPV$reg == XPrd.i$region[1])][1]
    XPred <- rbind(XPred, XPrd.i)
  }
  if (any(names(XPred) == "X1")) { colnames(XPred)[which(names(XPred) == "X1")] <- "Spl_t.1" }
  if (any(names(XPred) == "Spl_t")) { colnames(XPred)[which(names(XPred) == "Spl_t")] <- "Spl_t.1" }

  n_cnt <- length(unique(IPV$cnt))
  n_reg <- length(unique(IPV$reg))
  n_sup <- length(unique(IPV$sup))
      
  # We predict
  y.i <- inv.logit(res[["glb"]])
  y_pred <- NULL
  pb <- txtProgressBar(1, dim(XPred)[1], style = 3)
  ind_spl_a <- names(XPred) %in% paste("Spl_a", rep(1:SDat$Ndof), sep = ".")
  ind_spl_t <- names(XPred) %in% paste("Spl_t", rep(1:SDat$Nt), sep = ".")
  for (i in 1:dim(XPred)[1]) {
    ind_a <- seq(XPred$Region[i], SDat$Ndof * n_reg - n_reg + XPred$Region[i], n_reg)
    ind_t <- seq(XPred$Region[i], SDat$Nt * n_reg - n_reg + XPred$Region[i], n_reg)
    y.i <- inv.logit(res[["reg"]][, XPred$Region[i]] +
                    as.vector(res[["b_time_reg"]][, ind_t] %*% t(XPred[i, ind_spl_t])) +
                    as.vector((res[["b_spline_reg"]][, ind_a]) %*% t(XPred[i, ind_spl_a])))
    y_pred.i <- quantile(y.i, probs = c(0.5, 0.025, 0.975))
    y_pred <- rbind(y_pred, y_pred.i)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  row.names(y_pred) <- NULL
  colnames(y_pred) <- c("Median", "LCI", "UCI")
  XPAge <- cbind(XPred, y_pred)

  # Super Regions
  UIR <- unique(IPV$sup)
  XPred <- NULL
  for (i in 1:SDat$SR) {
    XPrd.i <- data.frame(Age = age, Spl_a = Xage, Spl_t = XTime)
    XPrd.i$SuperRegion <- as.character(UIR[i])
    XPrd.i$Super <- IPV$Super[which(IPV$sup == XPrd.i$SuperRegion[1])][1]
    XPred <- rbind(XPred, XPrd.i)
  }
  if (any(names(XPred) == "X1")) { colnames(XPred)[which(names(XPred) == "X1")] <- "Spl_t.1" }
  if (any(names(XPred) == "Spl_t")) { colnames(XPred)[which(names(XPred) == "Spl_t")] <- "Spl_t.1" }
  
  # We predict
  y.i <- inv.logit(res[["glb"]])
  y_pred <- NULL
  pb <- txtProgressBar(1, dim(XPred)[1], style = 3)
  ind_spl_a <- names(XPred) %in% paste("Spl_a", rep(1:SDat$Ndof), sep = ".")
  ind_spl_t <- names(XPred) %in% paste("Spl_t", rep(1:SDat$Nt), sep = ".")
  for (i in 1:dim(XPred)[1]) {
    ind_a <- seq(XPred$Super[i], SDat$Ndof * n_sup - n_sup + XPred$Super[i], n_sup)
    ind_t <- seq(XPred$Super[i], SDat$Nt * n_sup - n_sup + XPred$Super[i], n_sup)
    y.i <- inv.logit(res[["sup"]][, XPred$Super[i]] +
                    as.vector(res[["b_time_sup"]][, ind_t] %*% t(XPred[i, ind_spl_t])) +
                    as.vector((res[["b_spline_sup"]][, ind_a]) %*% t(XPred[i, ind_spl_a])))
    y_pred.i <- quantile(y.i, probs = c(0.5, 0.025, 0.975))
    y_pred <- rbind(y_pred, y_pred.i)
  }
  row.names(y_pred) <- NULL
  colnames(y_pred) <- c("Median", "LCI", "UCI")
  XSRAge <- cbind(XPred, y_pred)
  
  # overall
  if (!is.matrix(XTime)) { XTime <- as.matrix(XTime) }
  OAge <- NULL
  for (i in 1:dim(Xage)[1]) {
    y.i <- inv.logit(res[["glb"]] + 
                       as.vector(res[["b_time"]] %*% XTime[1, ]) +
                       as.vector(res[["b_spline"]] %*% Xage[i, ]))
    OAge.i <- quantile(y.i, probs = c(0.5, 0.025, 0.975))
    OAge <- rbind(OAge, OAge.i)
  }
  OAge <- as.data.frame(OAge)
  row.names(OAge) <- NULL
  colnames(OAge) <- c("Median", "LCI", "UCI")
  OPAge <- cbind(age, OAge)

  # We Graph
  if (outcome == "Ever IPV") { a <- 3 }
  if (outcome == "Past Year IPV") { a <- 2 }
  if (outcome == "NPSV") { a <- 1 }
  Col <- c(rgb(134, 179, 0, 255, max = 255), rgb(204, 0, 67, 255, max = 255),  rgb(255, 173, 51, 255, max = 255))
  ColT <- c(rgb(134, 179, 0, 125, max = 255), rgb(204, 0, 67, 125, max = 255), rgb(255, 173, 51, 125, max = 255))

  # First graph - overall only
  if (!is.null(dev.list())) { dev.off() }
  plot(OPAge$Median ~ OPAge$age,
      type = "n", ylim = c(0, 1), xlim = c(15, x_max), axes = FALSE,
      xlab = "Age", ylab = "Prevalence", main = paste(outcome, " - Overall Age Trend", sep = " "))
    lines(OPAge$Median ~ OPAge$age, col = Col[a], lwd = 2)
    polygon(x = c(OPAge$age, rev(OPAge$age)),
      y = c(OPAge$LCI, rev(OPAge$UCI)), col = ColT[a], border = NA)
    axis(1, seq(15, x_max, by = 10), labels = seq(15, x_max, by = 10))
    axis(2, seq(0, 1, by = 0.25), labels = seq(0, 100, by = 25)) 
  dev.copy(pdf, paste("Age Pattern Overall - ", outcome, " ", suffix, ".pdf", sep = ""), width = 6, height = 6)
  dev.off()

  # Age pattern by super region
  dev.off()
  pdf(paste("Age Pattern by Super Regions - ", outcome, " ", suffix, ".pdf", sep = ""), width = 10, height = 7.5)
  matrix_layout <- matrix(c(0, 0, 0, 1, 1, 0, 0, 0,
                            0, 2, 2, 3, 3, 4, 4, 0,
                            5, 5, 6, 6, 7, 7, 8, 8),
                          ncol = 8, nrow = 3, byrow = TRUE)    
  nf <- layout(matrix_layout, widths = rep(lcm(5 / 2), 8), heights = rep(lcm(5), 3), 
               respect = TRUE)
  
  par(mar = c(2, 1.5, 2, 1), oma = c(0, 3, 0, 0))
  # The overall graph
  OPAge <- subset(OPAge, age <= x_max)
  plot(OPAge$Median ~ OPAge$age,
       type = "n", ylim = c(-0, 1), xlim = c(14, x_max + 1),
       xlab = "Age", ylab = "Prevalence", main = "Overall Age Trend",
       axes = FALSE)
  lines(OPAge$Median ~ OPAge$age, col = Col[a], lwd = 2)
  polygon(x = c(OPAge$age, rev(OPAge$age)),
          y = c(OPAge$LCI, rev(OPAge$UCI)), col = ColT[a], border = NA)
  axis(1, seq(15, x_max, by = 10), labels = seq(15, x_max, by = 10))
  axis(2, seq(0, 1, by = 0.25), labels = seq(0, 100, by = 25))
  # The regional graphs
  UR <- unique(XSRAge$SuperRegion)
  for (i in 1:SDat$SR) {
    Xi <- subset(XSRAge, SuperRegion == UR[i] & Age <= x_max)
    plot(Xi$Median ~ Xi$Age,
         type = "n", ylim = c(-0, 1), xlim = c(14, x_max + 1),
         xlab = "", ylab = "", main = Xi$SuperRegion[1], axes = FALSE, cex.main = 1.1) 
    lines(Xi$Median ~ Xi$Age, col = Col[a], lwd = 2)
    polygon(x = c(Xi$Age, rev(Xi$Age)),
            y = c(Xi$LCI, rev(Xi$UCI)), col = ColT[a], border = NA)
    axis(1, seq(15, x_max, by = 10), labels = seq(15, x_max, by = 10))
    axis(2, seq(0, 1, by = 0.25), labels = seq(0, 100, by = 25))
  }
  mtext("Age", side = 1, line = -1.1, outer = TRUE)
  mtext(paste("Prevalence (%) ", outcome, sep = ""), side = 2, line = 0, outer = TRUE)
  dev.off()
  
  # Third - overal and region-specific patterns
  pdf(paste("Age Pattern by Regions - ", outcome, " ", suffix, ".pdf", sep = ""), width = 14.5, height = 8.5)
  matrix_layout <- matrix(c(0, 0, 0, 1, 0, 0, 0,
                            2, 3, 4, 5, 6, 7, 8,
                            9, 10, 11, 12, 13, 14, 15,
                            16, 17, 18, 19, 20, 21, 22),
                  ncol = 7, nrow = 4, byrow = TRUE)    
  nf <- layout(matrix_layout, widths = rep(lcm(5), 7), heights = rep(lcm(5), 4), 
               respect = TRUE)

    par(mar = c(2, 1.5, 2, 1), oma = c(0, 3, 0, 0))
    # The overall graph
    OPAge <- subset(OPAge, age <= x_max)
    plot(OPAge$Median ~ OPAge$age,
      type = "n", ylim = c(-0, 1), xlim = c(14, x_max + 1),
      xlab = "Age", ylab = "Prevalence", main = "Overall Age Trend",
      axes = FALSE)
    lines(OPAge$Median ~ OPAge$age, col = Col[a], lwd = 2)
    polygon(x = c(OPAge$age, rev(OPAge$age)),
      y = c(OPAge$LCI, rev(OPAge$UCI)), col = ColT[a], border = NA)
    axis(1, seq(15, x_max, by = 10), labels = seq(15, x_max, by = 10))
    axis(2, seq(0, 1, by = 0.25), labels = seq(0, 100, by = 25))
    # The regional graphs
    UR <- unique(XPAge$region)
    for (i in 1:SDat$R) {
      Xi <- subset(XPAge, region == UR[i] & Age <= x_max)
      plot(Xi$Median ~ Xi$Age,
        type = "n", ylim = c(-0, 1), xlim = c(14, x_max + 1),
        xlab = "", ylab = "", main = Xi$region[1], axes = FALSE, cex.main = 1.1) 
      lines(Xi$Median ~ Xi$Age, col = Col[a], lwd = 2)
      polygon(x = c(Xi$Age, rev(Xi$Age)),
        y = c(Xi$LCI, rev(Xi$UCI)), col = ColT[a], border = NA)
    axis(1, seq(15, x_max, by = 10), labels = seq(15, x_max, by = 10))
    axis(2, seq(0, 1, by = 0.25), labels = seq(0, 100, by = 25))
    }
    mtext("Age", side = 1, line = -1.1, outer = TRUE)
    mtext(paste("Prevalence (%) ", outcome, sep = ""), side = 2, line = 0, outer = TRUE)
  dev.off()
}

#' Plotting the age pattern by regions
plot_time_trend <- function(res, SDat, IPV, Center_t = 2018, 
                            outcome, suffix = "") {
  if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
    stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }
  # We create the database to predict
  time <- seq(2000, 2019, by = 1)

  # manipulating time trends   
  is_period <- length(attr(SDat$Spl_t, "period")) > 0  
  if (is_period) {
  XTime <- as.matrix(ifelse(time < attr(SDat$Spl_t, "cutoff"), 0, 1))
  }
  if (!is_period & length(attr(SDat$Spl_t, "knots")) > 0) {
  XTime <- as.matrix(splines::ns(time - Center_t,
                      knots = attr(SDat$Spl_t, "knots"),
                      Boundary.knots = attr(SDat$Spl_t, "Boundary.knots")))    
  } 
  if (!is_period & length(attr(SDat$Spl_t, "knots")) == 0) {
  XTime <- as.matrix(splines::ns(time - Center_t,
                      df = ncol(SDat$Spl_t),
                      Boundary.knots = attr(SDat$Spl_t, "Boundary.knots")))     
  }
  
  n_cnt <- length(unique(IPV$cnt))
  n_reg <- length(unique(IPV$reg))
  n_sup <- length(unique(IPV$sup))
  
  # Regions
  XPred <- NULL
  UIR <- unique(IPV$reg)
  XPred <- NULL
  for (i in 1:SDat$R) {
    XPrd.i <- data.frame(Time = time, Spl_t = XTime)
    XPrd.i$region <- as.character(UIR[i])
    XPrd.i$Region <- IPV$Region[which(IPV$reg == XPrd.i$region[1])][1]
    XPred <- rbind(XPred, XPrd.i)
  }
  if (any(names(XPred) == "X1")) { colnames(XPred)[which(names(XPred) == "X1")] <- "Spl_t.1" }
  if (any(names(XPred) == "Spl_t")) { colnames(XPred)[which(names(XPred) == "Spl_t")] <- "Spl_t.1" }
  
  # We predict
  y.i <- inv.logit(res[["glb"]])
  y_pred <- NULL
  pb <- txtProgressBar(1, dim(XPred)[1], style = 3)
  ind_spl_t <- names(XPred) %in% paste("Spl_t", rep(1:SDat$Nt), sep = ".")
  for (i in 1:dim(XPred)[1]) {
    ind_t <- seq(XPred$Region[i], SDat$Nt * n_reg - n_reg + XPred$Region[i], n_reg)
    y.i <- inv.logit(res[["reg"]][, XPred$Region[i]] +
                    as.vector(res[["b_time_reg"]][, ind_t] %*% t(XPred[i, ind_spl_t])))
    y_pred.i <- quantile(y.i, probs = c(0.5, 0.025, 0.975))
    y_pred <- rbind(y_pred, y_pred.i)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  row.names(y_pred) <- NULL
  colnames(y_pred) <- c("Median", "LCI", "UCI")
  XRTime <- cbind(XPred, y_pred)
  
  # Super Region
  UIR <- unique(IPV$sup) 
  NR <- SDat$SR
  XPred <- NULL
  for (i in 1:NR) {
    XPrd.i <- data.frame(Time = time, Spl_t = XTime)
    XPrd.i$SuperRegion <- as.character(UIR[i])
    XPrd.i$Super <- IPV$Super[which(IPV$sup == XPrd.i$SuperRegion[1])][1]
    XPred <- rbind(XPred, XPrd.i)
  }
  if (any(names(XPred) == "X1")) { colnames(XPred)[which(names(XPred) == "X1")] <- "Spl_t.1" }
  if (any(names(XPred) == "Spl_t")) { colnames(XPred)[which(names(XPred) == "Spl_t")] <- "Spl_t.1" }
  
  # We predict
  y_pred <- NULL
  pb <- txtProgressBar(1, nrow(XPred), style = 3)
  ind_spl_t <- names(XPred) %in% paste("Spl_t", rep(1:SDat$Nt), sep = ".")
  for (i in 1:nrow(XPred)) {
    ind_t <- seq(XPred$Super[i], SDat$Nt * n_sup - n_sup + XPred$Super[i], n_sup)
    y.i <- inv.logit(res[["sup"]][, XPred$Super[i]] + 
                    as.vector(res[["b_time_sup"]][, ind_t] %*% t(XPred[i, ind_spl_t]))) 
    y_pred.i <- quantile(y.i, probs = c(0.5, 0.025, 0.975))
    y_pred <- rbind(y_pred, y_pred.i)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  row.names(y_pred) <- NULL
  colnames(y_pred) <- c("Median", "LCI", "UCI")
  XPTime <- cbind(XPred, y_pred)
  
  OTime <- NULL
  XPred <- data.frame(Time = time, Spl_t = XTime)
  if (any(names(XPred) == "X1")) { colnames(XPred)[which(names(XPred) == "X1")] <- "Spl_t.1" }
  if (any(names(XPred) == "Spl_t")) { colnames(XPred)[which(names(XPred) == "Spl_t")] <- "Spl_t.1" }
  
  ind_spl_t <- names(XPred) %in% paste("Spl_t", rep(1:SDat$Nt), sep = ".")
  for (i in 1:nrow(XPred)) {
    y.i <- inv.logit(res[["glb"]] + 
                       as.vector(res[["b_time"]] %*% t(XPred[i, ind_spl_t])))
    OTime.i <- quantile(y.i, probs = c(0.5, 0.025, 0.975))
    OTime <- rbind(OTime, OTime.i)
  }
  OTime <- as.data.frame(OTime)
  row.names(OTime) <- NULL
  colnames(OTime) <- c("Median", "LCI", "UCI")
  OPTime <- cbind(time, OTime)
  
  # We Graph
  if (outcome == "Ever IPV") { a <- 3 }
  if (outcome == "Past Year IPV") { a <- 2 }
  if (outcome == "NPSV") { a <- 1 }
  Col <- c(rgb(134, 179, 0, 255, max = 255), rgb(204, 0, 67, 255, max = 255),  rgb(255, 173, 51, 255, max = 255))
  ColT <- c(rgb(134, 179, 0, 125, max = 255), rgb(204, 0, 67, 125, max = 255), rgb(255, 173, 51, 125, max = 255))
  
  # First graph - overall only
  if (!is.null(dev.list())) { dev.off() }
  plot(OPTime$Median ~ OPTime$time,
       type = "n", ylim = c(0, 1), xlim = c(2000, 2020), axes = FALSE,
       xlab = "Year", ylab = "Prevalence", main = paste(outcome, " - Overall Time Trend", sep = " "))
  lines(OPTime$Median ~ OPTime$time, col = Col[a], lwd = 2)
  polygon(x = c(OPTime$time, rev(OPTime$time)),
          y = c(OPTime$LCI, rev(OPTime$UCI)), col = ColT[a], border = NA)
  axis(1, seq(min(time), 2020, by = 5), labels = seq(min(time), 2020, by = 5))
  axis(2, seq(0, 1, by = 0.25), labels = seq(0, 100, by = 25)) 
  dev.copy(pdf, paste("Time Trend Overall - ", outcome, " ", suffix, ".pdf", sep = ""), width = 6, height = 6)
  dev.off()
  
  # Second graph - overal and region-specific patterns
  dev.off()
  pdf(paste("Time Trend by Super Regions - ", outcome, " ", suffix, ".pdf", sep = ""), width = 10, height = 7)
  matrix_layout <- matrix(c(0, 0, 0, 1, 1, 0, 0, 0,
                            0, 2, 2, 3, 3, 4, 4, 0,
                            5, 5, 6, 6, 7, 7, 8, 8),
                          ncol = 8, nrow = 3, byrow = TRUE)    
  nf <- layout(matrix_layout, widths = rep(lcm(5/2), 8), heights = rep(lcm(5), 4), 
               respect = TRUE)
  
  par(mar = c(2, 1.5, 2, 1), oma = c(0, 3, 0, 0))
  # The overall graph
  plot(OPTime$Median ~ OPTime$time,
       type = "n", ylim = c(0, 1), xlim = c(2000, 2019), axes = FALSE,
       xlab = "Year", ylab = "Prevalence", main = paste(outcome, " - Overall Time Trend", sep = " "))
  lines(OPTime$Median ~ OPTime$time, col = Col[a], lwd = 2)
  polygon(x = c(OPTime$time, rev(OPTime$time)),
          y = c(OPTime$LCI, rev(OPTime$UCI)), col = ColT[a], border = NA)
  axis(1, seq(min(time), 2020, by = 5), labels = seq(min(time), 2020, by = 5))
  axis(2, seq(0, 1, by = 0.25), labels = seq(0, 100, by = 25)) 
  # The regional graphs
  UIR <- unique(XPTime$SuperRegion)
  for (i in 1:SDat$SR) {
    Xi <- subset(XPTime, SuperRegion == UIR[i])
    plot(Xi$Median ~ Xi$Time,
         type = "n", ylim = c(-0, 1), xlim = c(2000, 2019),
         xlab = "", ylab = "", main = Xi$SuperRegion[1], axes = FALSE, cex.main = 1.1) 
    lines(Xi$Median ~ Xi$Time, col = Col[a], lwd = 2)
    polygon(x = c(Xi$Time, rev(Xi$Time)),
            y = c(Xi$LCI, rev(Xi$UCI)), col = ColT[a], border = NA)
    axis(1, seq(min(time), 2020, by = 5), labels = seq(min(time), 2020, by = 5))
    axis(2, seq(0, 1, by = 0.25), labels = seq(0, 100, by = 25))
  }
  mtext("Year", side = 1, line = -1.1, outer = TRUE)
  mtext(paste("Prevalence (%) ", outcome, sep = ""), side = 2, line = 0, outer = TRUE)
  dev.off()
  
  # Third - overal and region-specific patterns
  pdf(paste("Time Trend by Regions - ", outcome, " ", suffix, ".pdf", sep = ""), width = 14.5, height = 8.5)
  matrix_layout <- matrix(c(0, 0, 0, 1, 0, 0, 0,
                            2, 3, 4, 5, 6, 7, 8,
                            9, 10, 11, 12, 13, 14, 15,
                            16, 17, 18, 19, 20, 21, 22),
                          ncol = 7, nrow = 4, byrow = TRUE)    
  nf <- layout(matrix_layout, widths = rep(lcm(5), 7), heights = rep(lcm(5), 4), 
               respect = TRUE)
  
  par(mar = c(2, 1.5, 2, 1), oma = c(0, 3, 0, 0))
  # The overall graph
  plot(OPTime$Median ~ OPTime$time,
       type = "n", ylim = c(0, 1), xlim = c(2000, 2019), axes = FALSE,
       xlab = "Year", ylab = "Prevalence", main = paste(outcome, " - Overall Time Trend", sep = " "))
  lines(OPTime$Median ~ OPTime$time, col = Col[a], lwd = 2)
  polygon(x = c(OPTime$time, rev(OPTime$time)),
          y = c(OPTime$LCI, rev(OPTime$UCI)), col = ColT[a], border = NA)
  axis(1, seq(min(time), 2020, by = 5), labels = seq(min(time), 2020, by = 5))
  axis(2, seq(0, 1, by = 0.25), labels = seq(0, 100, by = 25)) 
  # The regional graphs
  UR <- unique(XRTime$region)
  for (i in 1:SDat$R) {
    Xi <- subset(XRTime, region == UR[i])
    plot(Xi$Median ~ Xi$Time,
         type = "n", ylim = c(-0, 1), xlim = c(2000, 2019),
         xlab = "", ylab = "", main = Xi$region[1], axes = FALSE, cex.main = 1.1) 
    lines(Xi$Median ~ Xi$Time, col = Col[a], lwd = 2)
    polygon(x = c(Xi$Time, rev(Xi$Time)),
            y = c(Xi$LCI, rev(Xi$UCI)), col = ColT[a], border = NA)
    axis(1, seq(min(time), 2020, by = 5), labels = seq(min(time), 2020, by = 5))
    axis(2, seq(0, 1, by = 0.25), labels = seq(0, 100, by = 25))
  }
  mtext("Year", side = 1, line = -1.1, outer = TRUE)
  mtext(paste("Prevalence (%) ", outcome, sep = ""), side = 2, line = 0, outer = TRUE)
  dev.off()
}

addition_to_cnt <- function(aggregate = "gbd") {
  if (!(aggregate %in% c("gbd", "who", "sdg", "unicef", "unfpa"))) {
    stop("aggregate must be gbd, who, or sdg")
  }
  
  if (aggregate == "gbd") {
      data(gbd)
      # We add a few countries for plotting results
      addition <- data.frame(iso = c("COK", "PLW", "NRU", "TUV", "NIU",
                                 "KNA", 
                                 "MCO", "SMR",
                                 "HKG"),
                              NameCnt = c("Cook Islands", "Palau", "Nauru", "Tuvalu", "Niue",
                                         "Saint Kitts and Nevis",  
                                         "Monaco", "San Marino",
                                         "Hong Kong"),
                              gbd = c(rep("Oceania", 5),
                                      rep("Caribbean", 1),
                                      rep("Western Europe", 2),
                                      rep("East Asia", 1)),
                              Region = c(rep("Oceania", 5),
                                      rep("Caribbean", 1),
                                      rep("Europe, Western", 2),
                                      rep("Asia, East", 1)),
                              SuperRegion = c(rep("South-East Asia, East Asia & Oceania", 5),
                                      rep("Latin America & Caribbean", 1),
                                      rep("Central Europe, Eastern Europe & Central Asia", 2),
                                      rep("South-East Asia, East Asia & Oceania", 1)))
  val <- rbind(gbd, addition)   
  }
  
  if (aggregate == "who") {
      data(who)
      #already included: "Cook Islands", "Palau", "Nauru", "Tuvalu", "Niue", "Saint Kitts and Nevis", "Monaco", "San Marino",
      addition <- data.frame(country = c("Hong Kong"),
                         iso3 = c("HKG"),
                         region = c(rep("WPRO", 1)),
                         region_name = c(rep("Western Pacific Region", 1)),
                         region_name_hi = c(rep("High Income Region", 1)))
      val <- rbind(who, addition)
  }
  
  if (aggregate == "sdg") {
      data(sdg)
      sdg <- sdg[, !(colnames(sdg) %in% c("intregion"))]
      addition <- data.frame(cnt = c("Cook Islands", "Palau", "Nauru", "Tuvalu", "Niue",
                                         "Saint Kitts and Nevis",  
                                         "Monaco", "San Marino",
                                         "Hong Kong"),
                         iso3 = c("COK", "PLW", "NRU", "TUV", "NIU",
                                 "KNA", 
                                 "MCO", "SMR",
                                 "HKG"),
                         region = c(rep("Oceania", 5),
                                      rep("Americas", 1),
                                      rep("Europe", 2),
                                      rep("Asia", 1)),
                         subregion = c(rep("Polynesia", 5),
                                      rep("Latin America and the Caribbean", 1),
                                      rep("Southern Europe", 2),
                                      rep("Eastern Asia", 1)),
                         least_developped = c(rep(1, 5),
                                       rep(0, 4)))
      val <- rbind(sdg, addition)
  }

  return(val)
}

#' Aggregate by country
pool_pred <- function(res, SDat, IPV, Center_a = 30, Center_t = 2018, 
                     outcome, year_pred = 2018,
                     all_women = FALSE, 
                     extroplated_ever_sex_wgt_2018 = TRUE) {
  
    if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
      stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }

    data("denom_cnt_age_now")
    if (extroplated_ever_sex_wgt_2018 == FALSE) { 
      data("denom_ever_sex_2010")
      print("you are using 2010 ever had sex weights")
      denom_ever_sex <- denom_ever_sex_2010 }
    if (extroplated_ever_sex_wgt_2018 == TRUE) {  
      data("denom_ever_sex_2018") 
      print("you are using 2018 ever had sex weights")
      denom_ever_sex <- denom_ever_sex_2018 }
    data("wpp_std_now")

    gbd <- addition_to_cnt("gbd")
    who <- addition_to_cnt("who") 
    sdg <- addition_to_cnt("sdg")

    n_cnt <- length(unique(IPV$cnt))
    n_reg <- length(unique(IPV$reg))
    n_sup <- length(unique(IPV$sup))
  
  # manipulating time trends    
  is_period <- length(attr(SDat$Spl_t, "period")) > 0 
  if (is_period) {
    XTime <- as.matrix(ifelse(year_pred < attr(SDat$Spl_t, "cutoff"), 0, 1))
  }
  if (!is_period & length(attr(SDat$Spl_t, "knots")) > 0) {
  XTime <- as.matrix(splines::ns(year_pred - Center_t,
                      knots = attr(SDat$Spl_t, "knots"),
                      Boundary.knots = attr(SDat$Spl_t, "Boundary.knots")))    
  } 
  if (!is_period & length(attr(SDat$Spl_t, "knots")) == 0) {
  XTime <- as.matrix(splines::ns(year_pred - Center_t,
                      df = ncol(SDat$Spl_t),
                      Boundary.knots = attr(SDat$Spl_t, "Boundary.knots")))       
  }
  
    # Some countries don't have data on ever had sex, we take the world average
    EverSexAvg <- aggregate(denom_ever_sex$ever_had_sex,
                            by = list(denom_ever_sex$age_group), 
                            FUN = mean, na.rm = TRUE)
    colnames(EverSexAvg) <- c('age_group','ever_had_sex')

    age <- seq(17.5, 104, 5)
    age65 <- ifelse(age > 65, 65, age)
    Xage <- splines::ns(age65 - Center_a,
                        knots = attr(SDat$Spl_W, "knots"),
                        Boundary.knots = attr(SDat$Spl_W, "Boundary.knots"))

    SRLkUp <- as.data.frame(unique(IPV[c("sup", "Super")]))
      colnames(SRLkUp) <- c("sup", "Super")
    RLkUp <- as.data.frame(unique(IPV[c("reg", "Region")]))
      colnames(RLkUp) <- c("reg", "Region")
    CLkUp <- as.data.frame(unique(IPV[c("cnt", "Country")]))
      colnames(CLkUp) <- c("cnt", "Country")
    
    UIC <- unique(gbd$iso)
    LA <- length(age)
    XPred <- NULL
    # We create the prediction matrix of Age-by-Country
    for (i in 1:length(UIC)) {
      XPrd.i <- data.frame(
            Age = age, Spl_a = Xage,
            ISO = rep(UIC[i], LA),
            GBD = rep(gbd$gbd[gbd$iso == UIC[i]], LA))
      XPrd.i$SuperRegion <- SRLkUp$Super[as.character(SRLkUp$sup) == as.character(gbd$SuperRegion[gbd$iso == UIC[i]])]
      MissRegion <- RLkUp$Region[as.character(RLkUp$reg) == as.character(gbd$Region[gbd$iso == UIC[i]])]
        XPrd.i$Region <- ifelse(length(MissRegion) == 0, NA, MissRegion)
      MissCountry <- CLkUp$Country[as.character(CLkUp$cnt) == as.character(UIC[i])]
        XPrd.i$Country <- ifelse(length(MissCountry) == 0, NA, MissCountry)
      XPred <- rbind(XPred, XPrd.i)
     }
     XPred$GBD <- droplevels(XPred$GBD)
     
     # We predict using the individual iterations
     y_pred <- array(NA, dim = c(length(unique(XPred$ISO)),
        length(unique(XPred$Age)),
        length(res[["glb"]])))
     y_denom <- array(NA, dim = c(length(unique(XPred$ISO)),
        length(unique(XPred$Age))))
      XPred$ISO_New <- recode.cluster(XPred$ISO)
      YP_Index <- NULL
        
      # if we are not modeling time and age trends at country levels
      if (all(names(res) != "sd_time_cnt")) {
         res[["b_time_cnt"]]  <- matrix(0, nrow = nrow(res[["b_time_reg"]]), ncol = SDat$Nt * n_cnt)
         res[["sd_time_cnt"]] <- matrix(0, nrow = nrow(res[["sd_time_reg"]]), ncol = SDat$Nt) 
         for (x in 1:n_cnt) {
           ind_t_cnt <- seq(x, SDat$Nt * n_cnt - n_cnt + x, n_cnt)
           ind_t_reg <- seq(SDat$RLookUp[x], SDat$Nt * n_reg - n_reg + SDat$RLookUp[x], n_reg)
           for (j in 1:SDat$Nt) {
            res[["b_time_cnt"]][, ind_t_cnt[j]] <- res[["b_time_reg"]][, ind_t_reg[j]]
           }
         }
      }
      if (all(names(res) != "sd_spline_cnt")) {
         res[["b_spline_cnt"]]  <- matrix(0, nrow = nrow(res[["b_spline_reg"]]), ncol = SDat$Ndof * n_cnt)
         res[["sd_spline_cnt"]] <- matrix(0, nrow = nrow(res[["sd_spline_reg"]]), ncol = SDat$Ndof) 
       for (x in 1:n_cnt) {
           ind_a_cnt <- seq(x, SDat$Ndof * n_cnt - n_cnt + x, n_cnt)
           ind_a_reg <- seq(SDat$RLookUp[x], SDat$Ndof * n_reg - n_reg + SDat$RLookUp[x], n_reg)
           for (j in 1:SDat$Nt) {
            res[["b_spline_cnt"]][, ind_a_cnt[j]] <- res[["b_spline_reg"]][, ind_a_reg[j]]
           }
        }
      }

      pb <- txtProgressBar(1, nrow(XPred), style = 3)
      ind_spl_a <- names(XPred) %in% paste("Spl_a", rep(1:SDat$Ndof), sep = ".")
      for (i in 1:nrow(XPred)) {
        Region_NA <- XPred$Region[i]
        Country_NA <- XPred$Country[i]
        YP_Index.i <- data.frame(It = i, Country = XPred$ISO_New[i], GBD = XPred$GBD[i],
          Region = XPred$Region[i], Super = XPred$SuperRegion[i], D1 = XPred$ISO_New[i], ISO = XPred$ISO[i])
        YP_Index <- rbind(YP_Index, YP_Index.i)
        # If we are predicting a country with data 
        if (!is.na(Country_NA)) {
          ind_a <- seq(XPred$Country[i], SDat$Ndof * n_cnt - n_cnt + XPred$Country[i], n_cnt)
          ind_t <- seq(XPred$Country[i], SDat$Nt * n_cnt - n_cnt + XPred$Country[i], n_cnt)
          y.i <- inv.logit(res[["cnt"]][, XPred$Country[i]] +
            as.vector(res[["b_time_cnt"]][, ind_t] %*% t(XTime)) +                
            as.vector(res[["b_spline_cnt"]][, ind_a] %*% t(XPred[i, ind_spl_a])))           
        }
        # If we are predicting a new country or region
        if (is.na(XPred$Country[i])) {
          # New country in a new Region
          if (is.na(Region_NA)) {
            var_new_cnt <- rnorm(length(res[["sd_cnt"]]), mean = 0, sd = res[["sd_cnt"]])
            var_new_reg <- rnorm(length(res[["sd_reg"]]), mean = 0, sd = res[["sd_reg"]])
            var_new_spline_cnt <- matrix(NA, nrow = nrow(res[["sd_spline_cnt"]]), ncol = SDat$Ndof) 
            var_new_spline_reg <- matrix(NA, nrow = nrow(res[["sd_spline_reg"]]), ncol = SDat$Ndof)
            for (z in 1:SDat$Ndof) {
              var_new_spline_cnt[, z] <- rnorm(nrow(res[["sd_spline_cnt"]]), mean = 0, sd = res[["sd_spline_cnt"]][, z]) 
              var_new_spline_reg[, z] <- rnorm(nrow(res[["sd_spline_reg"]]), mean = 0, sd = res[["sd_spline_reg"]][, z])
            }
            ind_new_a <- seq(XPred$SuperRegion[i], SDat$Ndof * n_sup - n_sup + XPred$SuperRegion[i], n_sup)
            
            spl_new_reg <- res[["b_spline_sup"]][, ind_new_a] + var_new_spline_cnt + var_new_spline_reg
           # time
           ind_new_t <- seq(XPred$SuperRegion[i], SDat$Nt * n_sup - n_sup + XPred$SuperRegion[i], n_sup)
           var_new_time_cnt <- matrix(0, nrow = nrow(res[["sd_time_cnt"]]), ncol = SDat$Nt) 
           var_new_time_reg <- matrix(0, nrow = nrow(res[["sd_time_reg"]]), ncol = SDat$Nt)
           for (z in 1:SDat$Nt) {
              var_new_time_reg[, z] <- rnorm(nrow(res[["sd_time_reg"]]), mean = 0, sd = res[["sd_time_reg"]][, z])
              var_new_time_cnt[, z] <- rnorm(nrow(res[["sd_time_cnt"]]), mean = 0, sd = res[["sd_time_cnt"]][, z]) 
            }
          time_new_trend <- res[["b_time_sup"]][, ind_new_t] + var_new_time_reg + var_new_time_cnt
          y.i <- inv.logit(res[["sup"]][, XPred$SuperRegion[i]] + var_new_reg + var_new_cnt +
                            as.vector(time_new_trend %*% t(XTime)) + 
                            as.vector(spl_new_reg %*% t(XPred[i, ind_spl_a])))
          } 
          # New Country in a known region
          if (!is.na(Region_NA)) {
            var_new_cnt <- rnorm(length(res[["sd_cnt"]]), mean = 0, sd = res[["sd_cnt"]])
            var_new_spline_cnt <- matrix(NA, nrow = nrow(res[["sd_spline_cnt"]]), ncol = SDat$Ndof)
            for (z in 1:SDat$Ndof) {
              var_new_spline_cnt[, z] <- rnorm(nrow(res[["sd_spline_cnt"]]), mean = 0, sd = res[["sd_spline_cnt"]][, z])
            }
            ind_new_a <- seq(XPred$Region[i], SDat$Ndof * n_reg - n_reg + XPred$Region[i], n_reg)
            spl_new_cnt <- res[["b_spline_reg"]][, ind_new_a] + var_new_spline_cnt
          # time
          ind_new_t <- seq(XPred$Region[i], SDat$Nt * n_reg - n_reg + XPred$Region[i], n_reg)
          var_new_time_cnt <- matrix(NA, nrow = nrow(res[["sd_time_cnt"]]), ncol = SDat$Nt)
          for (z in 1:SDat$Nt) {
            var_new_time_cnt[, z] <- rnorm(nrow(res[["sd_time_cnt"]]), mean = 0, sd = res[["sd_time_cnt"]][, z])
           }
          time_new_trend <- res[["b_time_reg"]][, ind_new_t] + var_new_time_cnt
            y.i <- inv.logit(res[["reg"]][, XPred$Region[i]] + var_new_cnt + 
                   as.vector(time_new_trend %*% t(XTime)) + 
                   as.vector(spl_new_cnt %*% t(XPred[i, ind_spl_a])))
          }
        }
        # Adding the Population Weights to create a new numerator/denominator.
        # For IPV
        if (outcome != "NPSV") {
          W.i <- subset(denom_cnt_age_now, as.character(iso3) == as.character(XPred$ISO[i]))
          EverSex.i <- subset(denom_ever_sex, as.character(iso3) == as.character(XPred$ISO[i]))
          if (nrow(EverSex.i) == 0) { EverSex.i <- EverSexAvg } # There are some countries without data for Ever_Had_Sex.
          if (nrow(W.i) > 0) {
            closest_age <- which.min(abs(EverSex.i$age_group - XPred$Age[i])) # b/c ever sex is only up to age 45-49.
            w.i <- (W.i$women[W.i$Age == XPred$Age[i]] + 0.5) * EverSex.i$ever_had_sex[closest_age]
            y_pred[XPred$ISO_New[i], which(age == XPred$Age[i]), ] <- y.i * w.i
            y_denom[XPred$ISO_New[i], which(age == XPred$Age[i])] <- w.i
            
            # if we shall express IPV among all women
            if (all_women == TRUE) {
              y_pred[XPred$ISO_New[i], which(age == XPred$Age[i]), ] <- (y.i * EverSex.i$ever_had_sex[closest_age]) * (W.i$women[W.i$Age == XPred$Age[i]] + 0.5)
              y_denom[XPred$ISO_New[i], which(age == XPred$Age[i])] <- (W.i$women[W.i$Age == XPred$Age[i]] + 0.5)
            }
            
          } else {  # Countries with a population less than 90K are not takent into account.
            # we take the population of the world and assume that the country has a 45K / 2 = 22.5 pop size
            W.i <- wpp_std_now$Pop[wpp_std_now$AgeM == XPred$Age[i]] / sum(wpp_std_now$Pop) * 22.5 
            closest_age <- which.min(abs(EverSex.i$age_group - XPred$Age[i])) # b/c ever sex is only up to age 45-49.
            w.i <- W.i * EverSex.i$ever_had_sex[closest_age]
            y_pred[XPred$ISO_New[i], 
                    which(age == XPred$Age[i]), ] <- y.i * w.i
            y_denom[XPred$ISO_New[i], 
                     which(age == XPred$Age[i])] <- w.i
            
            # if we shall express IPV among all women
            if (all_women == TRUE) {
              y_pred[XPred$ISO_New[i], which(age == XPred$Age[i]), ] <- (y.i * EverSex.i$ever_had_sex[closest_age]) * W.i
              y_denom[XPred$ISO_New[i], which(age == XPred$Age[i])] <- W.i
            }
            
          } }
        # For NPSV
        if (outcome == "NPSV") {
          W.i <- subset(denom_cnt_age_now, as.character(iso3) == as.character(XPred$ISO[i]))
          if (nrow(W.i) > 0) {
            w.i <- W.i$women[W.i$Age == XPred$Age[i]] + 0.5
            y_pred[XPred$ISO_New[i], which(age == XPred$Age[i]), ] <- y.i * w.i
            y_denom[XPred$ISO_New[i], which(age == XPred$Age[i])] <- w.i
          } else { # Countries with a population less than 90K are not takent into account.
            w.i <- wpp_std_now$Pop[wpp_std_now$AgeM == XPred$Age[i]] / sum(wpp_std_now$Pop) * 22.5
            y_pred[XPred$ISO_New[i], 
                    which(age == XPred$Age[i]), ] <- y.i * w.i
            y_denom[XPred$ISO_New[i], 
                     which(age == XPred$Age[i])] <- w.i
          }}
        setTxtProgressBar(pb, i)
      }
      close(pb)

      # we add the index to find WHO and SDG regions
      YP_Index$WHO <- who$region_name[match(YP_Index$ISO, who$iso3)]
      YP_Index$WHO_HI <- who$region_name_hi[match(YP_Index$ISO, who$iso3)]
      YP_Index$SDG <- sdg$subregion[match(YP_Index$ISO, sdg$iso3)]      
      YP_Index$SDGsup <- sdg$region[match(YP_Index$ISO, sdg$iso3)]   
      YP_Index$SDGdev <- sdg$least_developped[match(YP_Index$ISO, sdg$iso3)]   
      
    data(unicef)
    data(unfpa)
     YP_Index$unfpa <- unfpa$regions[match(YP_Index$ISO, unfpa$iso3)]
     YP_Index$unicef_prg <- unicef$program_region[match(YP_Index$ISO, unicef$iso3)]
     YP_Index$unicef_reg1 <- unicef$reporting_region1[match(YP_Index$ISO, unicef$iso3)]
     YP_Index$unicef_reg2 <- unicef$reporting_region2[match(YP_Index$ISO, unicef$iso3)]
     YP_Index$unfpa <- ifelse(is.na(YP_Index$unfpa), "not included", YP_Index$unfpa)
     YP_Index$unicef_prg <- ifelse(is.na(YP_Index$unicef_prg), "not included", YP_Index$unicef_prg)     
     YP_Index$unicef_reg1 <- ifelse(is.na(YP_Index$unicef_reg1), "not included", YP_Index$unicef_reg1) 
     YP_Index$unicef_reg2 <- ifelse(is.na(YP_Index$unicef_reg2), "not included", YP_Index$unicef_reg2)
     
      YPInd <- YP_Index[!duplicated(YP_Index[, -1]), ]
      
      val <- list(YPInd = YPInd, y_denom = y_denom, y_pred = y_pred)
      return(val)
}


pool_cnt <- function(pred, IPV, outcome,
                     plot_curves = TRUE, save_results = TRUE, suffix = "") {
  if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
    stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }

    YPInd <- pred[["YPInd"]]
    y_denom <- pred[["y_denom"]]
    y_pred <- pred[["y_pred"]]

    age <- seq(17.5, 104, 5)

    # We do the aggregation by Country-AgeGrp here
      AgeGrp <- data.frame(
        loage = c(15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 15, 15, 15, 25, 35, 50),
        hiage = c(19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 99, 49, 99, 24, 34, 49, 99))
      rownames(AgeGrp) <- NULL
      AggCnt <- NULL

    pb <- txtProgressBar(1, nrow(y_pred), style = 3)
    for (i in 1:nrow(y_pred)) {
      Prev <- NULL
      for (j in 1:nrow(AgeGrp)) {
        Pop.ij <- y_denom[YPInd$D1[i], which(age > AgeGrp$loage[j] & age < AgeGrp$hiage[j])]
        if (any(is.na(Pop.ij))) {
          Prev <- rbind(Prev, c(NA, NA, NA))
        } else {
          Age.ij <- y_pred[YPInd$D1[i], which(age > AgeGrp$loage[j] & age < AgeGrp$hiage[j]), ]
          if (dim(as.matrix(Age.ij))[2] > 1) {
            Age.ij <- colSums(y_pred[YPInd$D1[i], which(age > AgeGrp$loage[j] & age < AgeGrp$hiage[j]), ])
          }
          W.ij <- Age.ij / sum(Pop.ij)
          Prev.ij <- quantile(W.ij, probs = c(0.5, 0.025, 0.975))
          Prev <- rbind(Prev, Prev.ij)
        }}
      rownames(Prev) <- NULL
      colnames(Prev) <- c("Median", "LCI", "UCI")
      Agg.i <- cbind(AgeGrp, Prev)
      Agg.i$Country <- as.character(YPInd$ISO[i])
      Agg.i$GBD <- YPInd$GBD[i]
      AggCnt <- rbind(AggCnt, Agg.i)
      setTxtProgressBar(pb, i)
    }
    close(pb)

    AggCnt$name <- countrycode::countrycode(sourcevar = AggCnt$Country, 
                   origin = "iso3c", destination = "un.name.en", warn = TRUE,
                   custom_match = c(
                     "ASM" = "American Samoa",
                     "BMU" = "Bermuda",
                     "COK" = "Cook Islands",
                     "GRL" = "Greenland",
                     "GUM" = "Guam",
                     "HKG" = "Hong Kong (S.A.R. China)",
                     "MNP" = "Northern Mariana Islands",
                     "NIU" = "Niue",
                     "PRI" = "Puerto Rico",
                     "PSE" = "Palestinian Territory, Occupied",
                     "TWN" = "Taiwan (Republic of China)",
                     "VIR" = "Virgin Islands (USA)"))
    AggCnt$with_data <- ifelse(as.character(AggCnt$Country) %in% unique(IPV$iso3), "-", "No data")

  if (save_results == TRUE) {
    write.csv(AggCnt, file = paste("Estimates Age by Cnt - ", outcome, " ", suffix, ".csv", sep = ""))
  }

  if (plot_curves == TRUE) {
    #' Plot age by country.
    AggCntToPlot <- subset(AggCnt, (hiage - loage) < 5)
    AggCntToPlot$Age <- AggCntToPlot$loage + (AggCntToPlot$hiage - AggCntToPlot$loage) / 2
    plot(AggCntToPlot$Median ~ AggCntToPlot$Age,
         type = "n", ylim = c(0, 1),
         xlab = "Age Groups", ylab = outcome, axes = TRUE)
    UID <- unique(AggCntToPlot$Country)
    for (i in 1:length(UID)) {
      X.i <- subset(AggCntToPlot, Country == UID[i])
      lines(X.i$Median ~ X.i$Age, col = rgb(50, 50, 50, 100, max = 255), lwd = 2)
    } }

  return(AggCnt)
}


pool_reg <- function(pred, ll_age = 15, ul_age = 49,
                     outcome, save_results = TRUE, suffix = "") {

  if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
    stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }

  YPInd <- pred[["YPInd"]]
  y_denom <- pred[["y_denom"]]
  y_pred <- pred[["y_pred"]]
  age <- seq(17.5, 104, by = 5)

  RPrv <- NULL
  R <- unique(YPInd$GBD)
  pb <- txtProgressBar(1, length(R), style = 3)
  for (i in R) {
    D1 <- YPInd$D1[YPInd$GBD == i]
    # Prevalence in age group
    sum.ij <- colSums(colSums(y_pred[D1, which(age >= ll_age & age <= ul_age), ]))
    wom.ij <- sum.ij * 1000
    denom.ij <- sum(colSums(y_denom[D1, which(age >= ll_age & age <= ul_age)]))
    RPrv.ij <- sum.ij / denom.ij
    RPrv.i <- data.frame(GBD = i,
                        Median = quantile(RPrv.ij, probs = 0.5),
                        LCI = quantile(RPrv.ij, probs = 0.025),
                        UCI = quantile(RPrv.ij, probs = 0.975),
                        women = quantile(wom.ij, probs = 0.5),
                        women_lci = quantile(wom.ij, probs = 0.025),
                        women_uci = quantile(wom.ij, probs = 0.975))
  RPrv <- rbind(RPrv, RPrv.i)
  setTxtProgressBar(pb, which(R == i))
  }
  close(pb)
  rownames(RPrv) <- NULL
  RPrv_AgeStd <- RPrv

  RPrv_AgeStd$women <- round(RPrv_AgeStd$women / 1000, 0) * 1000
  RPrv_AgeStd$women_lci <- round(RPrv_AgeStd$women_lci / 1000, 0) * 1000
  RPrv_AgeStd$women_uci <- round(RPrv_AgeStd$women_uci / 1000, 0) * 1000
    
  if (save_results == TRUE) {
    write.csv(RPrv_AgeStd , file = paste("Estimates by GBD Region (",
                                   ll_age, "-", ul_age, ") - ", outcome,
                                   " ", suffix, ".csv", sep = "")) }
  return(RPrv_AgeStd)
}


pool_reg_who <- function(pred, ll_age = 15, ul_age = 49,
                     outcome, save_results = TRUE, suffix = "") {

  if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
    stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }

  YPInd <- pred[["YPInd"]]
  y_denom <- pred[["y_denom"]]
  y_pred <- pred[["y_pred"]]
  age <- seq(17.5, 104, by = 5)

  # not including high-income as a seperate category
  RPrv <- NULL
  R <- unique(as.character(YPInd$WHO))
  pb <- txtProgressBar(1, length(R), style = 3)
  for (i in R) {
    D1 <- YPInd$D1[YPInd$WHO == i]
    # Prevalence in age group
    sum.ij <- colSums(colSums(y_pred[D1, which(age >= ll_age & age <= ul_age), ]))
    wom.ij <- sum.ij * 1000    
    denom.ij <- sum(colSums(y_denom[D1, which(age >= ll_age & age <= ul_age)]))
    RPrv.ij <- sum.ij / denom.ij
    RPrv.i <- data.frame(WHO = i,
                        Median = quantile(RPrv.ij, probs = 0.5),
                        LCI = quantile(RPrv.ij, probs = 0.025),
                        UCI = quantile(RPrv.ij, probs = 0.975),
                        women = quantile(wom.ij, probs = 0.5),
                        women_lci = quantile(wom.ij, probs = 0.025),
                        women_uci = quantile(wom.ij, probs = 0.975))
  RPrv <- rbind(RPrv, RPrv.i)
  setTxtProgressBar(pb, which(R == i))
  }
  close(pb)
  rownames(RPrv) <- NULL
  RPrv_AgeStd <- RPrv
  
  RPrv_AgeStd$women <- round(RPrv_AgeStd$women / 1000, 0) * 1000
  RPrv_AgeStd$women_lci <- round(RPrv_AgeStd$women_lci / 1000, 0) * 1000
  RPrv_AgeStd$women_uci <- round(RPrv_AgeStd$women_uci / 1000, 0) * 1000
  
  if (save_results == TRUE) {
    write.csv(RPrv_AgeStd , file = paste("Estimates by WHO Regions (",
                                   ll_age, "-", ul_age, ") - ", outcome,
                                   " ", suffix, ".csv", sep = "")) }


  # including high-income as a separate category
  RPrv <- NULL
  R <- unique(as.character(YPInd$WHO_HI))
  pb <- txtProgressBar(1, length(R), style = 3)
  for (i in R) {
    D1 <- YPInd$D1[YPInd$WHO_HI == i]
    # Prevalence in age group
    sum.ij <- colSums(colSums(y_pred[D1, which(age >= ll_age & age <= ul_age), ]))
    wom.ij <- sum.ij * 1000
    denom.ij <- sum(colSums(y_denom[D1, which(age >= ll_age & age <= ul_age)]))
    RPrv.ij <- sum.ij / denom.ij
    RPrv.i <- data.frame(WHO_HI = i,
                        Median = quantile(RPrv.ij, probs = 0.5),
                        LCI = quantile(RPrv.ij, probs = 0.025),
                        UCI = quantile(RPrv.ij, probs = 0.975),
                        women = quantile(wom.ij, probs = 0.5),
                        women_lci = quantile(wom.ij, probs = 0.025),
                        women_uci = quantile(wom.ij, probs = 0.975))
  RPrv <- rbind(RPrv, RPrv.i)
  setTxtProgressBar(pb, which(R == i))
  }
  close(pb)
  rownames(RPrv) <- NULL
  RPrv_AgeStd <- RPrv

  RPrv_AgeStd$women <- round(RPrv_AgeStd$women / 1000, 0) * 1000
  RPrv_AgeStd$women_lci <- round(RPrv_AgeStd$women_lci / 1000, 0) * 1000
  RPrv_AgeStd$women_uci <- round(RPrv_AgeStd$women_uci / 1000, 0) * 1000
  
  if (save_results == TRUE) {
    write.csv(RPrv_AgeStd , file = paste("Estimates by WHO & Regions High Income (",
                                   ll_age, "-", ul_age, ") - ", outcome,
                                   " ", suffix, ".csv", sep = "")) }
  return(RPrv_AgeStd)
}


pool_reg_sdg <- function(pred, ll_age = 15, ul_age = 49,
                     outcome, save_results = TRUE, suffix = "") {

  if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
    stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }

  YPInd <- pred[["YPInd"]]
  y_denom <- pred[["y_denom"]]
  y_pred <- pred[["y_pred"]]
  age <- seq(17.5, 104, by = 5)

  RPrv <- NULL
  R <- unique(YPInd$SDG)
  pb <- txtProgressBar(1, length(R), style = 3)
  for (i in R) {
    D1 <- YPInd$D1[YPInd$SDG == i]
    # Prevalence in age group
    sum.ij <- colSums(colSums(y_pred[D1, which(age >= ll_age & age <= ul_age), ]))
    wom.ij <- sum.ij * 1000
    denom.ij <- sum(colSums(y_denom[D1, which(age >= ll_age & age <= ul_age)]))
    RPrv.ij <- sum.ij / denom.ij
  RPrv.i <- data.frame(SDG = i,
                        Median = quantile(RPrv.ij, probs = 0.5),
                        LCI = quantile(RPrv.ij, probs = 0.025),
                        UCI = quantile(RPrv.ij, probs = 0.975),
                        women = quantile(wom.ij, probs = 0.5),
                        women_lci = quantile(wom.ij, probs = 0.025),
                        women_uci = quantile(wom.ij, probs = 0.975))
  RPrv <- rbind(RPrv, RPrv.i)
  setTxtProgressBar(pb, which(R == i))
  }
  close(pb)
  rownames(RPrv) <- NULL
  RPrv_AgeStd <- RPrv

  # super region SDG
  RPrv <- NULL
  R <- unique(YPInd$SDGsup)
  pb <- txtProgressBar(1, length(R), style = 3)
  for (i in R) {
    D1 <- YPInd$D1[YPInd$SDGsup == i]
    # Prevalence in age group
    sum.ij <- colSums(colSums(y_pred[D1, which(age >= ll_age & age <= ul_age), ]))
    wom.ij <- sum.ij * 1000
    denom.ij <- sum(colSums(y_denom[D1, which(age >= ll_age & age <= ul_age)]))
    RPrv.ij <- sum.ij / denom.ij
    RPrv.i <- data.frame(SDG_SuperRegion = i,
                        Median = quantile(RPrv.ij, probs = 0.5),
                        LCI = quantile(RPrv.ij, probs = 0.025),
                        UCI = quantile(RPrv.ij, probs = 0.975),
                        women = quantile(wom.ij, probs = 0.5),
                        women_lci = quantile(wom.ij, probs = 0.025),
                        women_uci = quantile(wom.ij, probs = 0.975))
  RPrv <- rbind(RPrv, RPrv.i)
  setTxtProgressBar(pb, which(R == i))
  }
  close(pb)
  rownames(RPrv) <- NULL
  SDGsup_AgeStd <- RPrv
  
  RPrv <- NULL
  R <- unique(YPInd$SDGdev)
  pb <- txtProgressBar(1, length(R), style = 3)
  for (i in R) {
    D1 <- YPInd$D1[YPInd$SDGdev == i]
    # Prevalence in age group
    sum.ij <- colSums(colSums(y_pred[D1, which(age >= ll_age & age <= ul_age), ]))
    wom.ij <- sum.ij * 1000
    denom.ij <- sum(colSums(y_denom[D1, which(age >= ll_age & age <= ul_age)]))
    RPrv.ij <- sum.ij / denom.ij
    RPrv.i <- data.frame(Least_Developped = i,
                        Median = quantile(RPrv.ij, probs = 0.5),
                        LCI = quantile(RPrv.ij, probs = 0.025),
                        UCI = quantile(RPrv.ij, probs = 0.975),
                        women = quantile(wom.ij, probs = 0.5),
                        women_lci = quantile(wom.ij, probs = 0.025),
                        women_uci = quantile(wom.ij, probs = 0.975))
  RPrv <- rbind(RPrv, RPrv.i)
  setTxtProgressBar(pb, which(R == i))
  }
  close(pb)
  rownames(RPrv) <- NULL
  SDGdev_AgeStd <- RPrv
  
  RPrv_AgeStd$women <- round(RPrv_AgeStd$women / 1000, 0) * 1000
  RPrv_AgeStd$women_lci <- round(RPrv_AgeStd$women_lci / 1000, 0) * 1000
  RPrv_AgeStd$women_uci <- round(RPrv_AgeStd$women_uci / 1000, 0) * 1000
  
  SDGsup_AgeStd$women <- round(SDGsup_AgeStd$women / 1000, 0) * 1000
  SDGsup_AgeStd$women_lci <- round(SDGsup_AgeStd$women_lci / 1000, 0) * 1000
  SDGsup_AgeStd$women_uci <- round(SDGsup_AgeStd$women_uci / 1000, 0) * 1000
  
  SDGdev_AgeStd$women <- round(SDGdev_AgeStd$women / 1000, 0) * 1000
  SDGdev_AgeStd$women_lci <- round(SDGdev_AgeStd$women_lci / 1000, 0) * 1000
  SDGdev_AgeStd$women_uci <- round(SDGdev_AgeStd$women_uci / 1000, 0) * 1000
  
  if (save_results == TRUE) {
    write.csv(RPrv_AgeStd, file = paste("Estimates by SDG Regions (",
                                   ll_age, "-", ul_age, ") - ", outcome,
                                   " ", suffix, ".csv", sep = ""))
    write.csv(SDGsup_AgeStd, file = paste("Estimates by SDG Super Regions (",
                                   ll_age, "-", ul_age, ") - ", outcome,
                                   " ", suffix, ".csv", sep = ""))
    write.csv(SDGdev_AgeStd, file = paste("Estimates by SDG Least-Developped (",
                                   ll_age, "-", ul_age, ") - ", outcome,
                                   " ", suffix, ".csv", sep = ""))
  }
  
  return(RPrv_AgeStd)
}

pool_reg_unicef <- function(pred, ll_age = 15, ul_age = 49,
                            outcome, save_results = TRUE, suffix = "") {

  if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
    stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }

  YPInd <- pred[["YPInd"]]
  y_denom <- pred[["y_denom"]]
  y_pred <- pred[["y_pred"]]
  age <- seq(17.5, 104, by = 5)

  RPrv <- NULL
  R <- unique(YPInd$unicef_prg)
  R <- R[-which(R == "not included")]
  pb <- txtProgressBar(1, length(R), style = 3)
  for (i in R) {
    D1 <- YPInd$D1[YPInd$unicef_prg == i]
    # Prevalence in age group
    sum.ij <- colSums(colSums(y_pred[D1, which(age >= ll_age & age <= ul_age), ]))
    wom.ij <- sum.ij * 1000
    denom.ij <- sum(colSums(y_denom[D1, which(age >= ll_age & age <= ul_age)]))
    RPrv.ij <- sum.ij / denom.ij
    RPrv.i <- data.frame(program_region = i,
                        Median = quantile(RPrv.ij, probs = 0.5),
                        LCI = quantile(RPrv.ij, probs = 0.025),
                        UCI = quantile(RPrv.ij, probs = 0.975),
                        women = quantile(wom.ij, probs = 0.5),
                        women_lci = quantile(wom.ij, probs = 0.025),
                        women_uci = quantile(wom.ij, probs = 0.975))
    RPrv <- rbind(RPrv, RPrv.i)
    setTxtProgressBar(pb, which(R == i))
  }
  close(pb)
  rownames(RPrv) <- NULL
  RPrv_program_region <- RPrv

  RPrv_program_region$women <- round(RPrv_program_region$women / 1000, 0) * 1000
  RPrv_program_region$women_lci <- round(RPrv_program_region$women_lci / 1000, 0) * 1000
  RPrv_program_region$women_uci <- round(RPrv_program_region$women_uci / 1000, 0) * 1000
  
  # super region SDG
  RPrv <- NULL
  R <- unique(YPInd$unicef_reg1)
  R <- R[-which(R == "not included")]
  pb <- txtProgressBar(1, length(R), style = 3)
  for (i in R) {
    D1 <- YPInd$D1[YPInd$unicef_reg1 == i]
    # Prevalence in age group
    sum.ij <- colSums(colSums(y_pred[D1, which(age >= ll_age & age <= ul_age), ]))
    wom.ij <- sum.ij * 1000
    denom.ij <- sum(colSums(y_denom[D1, which(age >= ll_age & age <= ul_age)]))
    RPrv.ij <- sum.ij / denom.ij
    RPrv.i <- data.frame(reporting_region1 = i,
                        Median = quantile(RPrv.ij, probs = 0.5),
                        LCI = quantile(RPrv.ij, probs = 0.025),
                        UCI = quantile(RPrv.ij, probs = 0.975),
                        women = quantile(wom.ij, probs = 0.5),
                        women_lci = quantile(wom.ij, probs = 0.025),
                        women_uci = quantile(wom.ij, probs = 0.975))
  RPrv <- rbind(RPrv, RPrv.i)
  setTxtProgressBar(pb, which(R == i))
  }
  close(pb)
  rownames(RPrv) <- NULL
  RPrv_reporting_region1 <- RPrv

  RPrv_reporting_region1$women <- round(RPrv_reporting_region1$women / 1000, 0) * 1000
  RPrv_reporting_region1$women_lci <- round(RPrv_reporting_region1$women_lci / 1000, 0) * 1000
  RPrv_reporting_region1$women_uci <- round(RPrv_reporting_region1$women_uci / 1000, 0) * 1000
  
  RPrv <- NULL
  R <- unique(YPInd$unicef_reg2)
  R <- R[-which(R == "not included")]
  pb <- txtProgressBar(1, length(R), style = 3)
  for (i in R) {
    D1 <- YPInd$D1[YPInd$unicef_reg2 == i]
    # Prevalence in age group
    sum.ij <- colSums(colSums(y_pred[D1, which(age >= ll_age & age <= ul_age), ]))
    wom.ij <- round(sum.ij) * 1000
    denom.ij <- sum(colSums(y_denom[D1, which(age >= ll_age & age <= ul_age)]))
    RPrv.ij <- sum.ij / denom.ij
    RPrv.i <- data.frame(reporting_region2 = i,
                        Median = quantile(RPrv.ij, probs = 0.5),
                        LCI = quantile(RPrv.ij, probs = 0.025),
                        UCI = quantile(RPrv.ij, probs = 0.975),
                        women = quantile(wom.ij, probs = 0.5),
                        women_lci = quantile(wom.ij, probs = 0.025),
                        women_uci = quantile(wom.ij, probs = 0.975))
  RPrv <- rbind(RPrv, RPrv.i)
  setTxtProgressBar(pb, which(R == i))
  }
  close(pb)
  rownames(RPrv) <- NULL
  RPrv_reporting_region2 <- RPrv

  RPrv_reporting_region2$women <- round(RPrv_reporting_region2$women / 1000, 0) * 1000
  RPrv_reporting_region2$women_lci <- round(RPrv_reporting_region2$women_lci / 1000, 0) * 1000
  RPrv_reporting_region2$women_uci <- round(RPrv_reporting_region2$women_uci / 1000, 0) * 1000
  
  if (save_results == TRUE) {
    write.csv(RPrv_program_region , file = paste("Estimates by UNICEF Program Regions (",
                                   ll_age, "-", ul_age, ") - ", outcome,
                                   " ", suffix, ".csv", sep = ""))
    write.csv(RPrv_reporting_region1 , file = paste("Estimates by UNICEF Reporting 1 Regions (",
                                   ll_age, "-", ul_age, ") - ", outcome,
                                   " ", suffix, ".csv", sep = ""))
    write.csv(RPrv_reporting_region2 , file = paste("Estimates by UNICEF Reporting 2 Regions (",
                                   ll_age, "-", ul_age, ") - ", outcome,
                                   " ", suffix, ".csv", sep = ""))
  }
  
  return(RPrv_program_region)
}

pool_reg_unfpa <- function(pred, ll_age = 15, ul_age = 49,
                           outcome, save_results = TRUE, suffix = "") {

  if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
    stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }

  YPInd <- pred[["YPInd"]]
  y_denom <- pred[["y_denom"]]
  y_pred <- pred[["y_pred"]]
  age <- seq(17.5, 104, by = 5)

  RPrv <- NULL
  R <- unique(YPInd$unfpa)
  R <- R[-which(R == "not included")]
  pb <- txtProgressBar(1, length(R), style = 3)
  for (i in R) {
    D1 <- YPInd$D1[YPInd$unfpa == i]
    # Prevalence in age group
    sum.ij <- colSums(colSums(y_pred[D1, which(age >= ll_age & age <= ul_age), ]))
    wom.ij <- sum.ij * 1000
    denom.ij <- sum(colSums(y_denom[D1, which(age >= ll_age & age <= ul_age)]))
    RPrv.ij <- sum.ij / denom.ij
    RPrv.i <- data.frame(unfpa = i,
                        Median = quantile(RPrv.ij, probs = 0.5),
                        LCI = quantile(RPrv.ij, probs = 0.025),
                        UCI = quantile(RPrv.ij, probs = 0.975),
                        women = quantile(wom.ij, probs = 0.5),
                        women_lci = quantile(wom.ij, probs = 0.025),
                        women_uci = quantile(wom.ij, probs = 0.975))
    RPrv <- rbind(RPrv, RPrv.i)
    setTxtProgressBar(pb, which(R == i))
  }
  close(pb)
  rownames(RPrv) <- NULL
  RPrv_unfpa <- RPrv

  RPrv_unfpa$women <- round(RPrv_unfpa$women / 1000, 0) * 1000
  RPrv_unfpa$women_lci <- round(RPrv_unfpa$women_lci / 1000, 0) * 1000
  RPrv_unfpa$women_uci <- round(RPrv_unfpa$women_uci / 1000, 0) * 1000
  
  if (save_results == TRUE) {
    write.csv(RPrv_unfpa, file = paste("Estimates by UNFPA Regions (",
                                   ll_age, "-", ul_age, ") - ", outcome,
                                   " ", suffix, ".csv", sep = ""))
  }
  
  return(RPrv_unfpa)
}

pool_glb <- function(pred, ll_age = 15, ul_age = 49, outcome,
                     cnt_exclude = NULL, suffix, save_results = TRUE) {

  if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
    stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }

  if (!is.null(cnt_exclude)) {
      to_select <- !(pred$YPInd$ISO %in% cnt_exclude) 
      pred$YPInd <- pred$YPInd[to_select, ]
      pred$y_denom <- pred$y_denom[to_select, ]
      pred$y_pred <- pred$y_pred[to_select, , ]
  }
  
  YPInd <- pred[["YPInd"]]
  y_denom <- pred[["y_denom"]]
  y_pred <- pred[["y_pred"]]
  age <- seq(17.5, 104, by = 5)

  # Prevalence 15+
  sum.i <- colSums(colSums(y_pred[, , ]))
  wom_all <- sum.i
  denom.i <- sum(colSums(y_denom[, ]))
  Prv_all <- sum.i / denom.i
  # Prevalence 15-49
  sum.i <- colSums(colSums(y_pred[, which(age >= ll_age & age <= ul_age), ]))
  wom_1549 <- sum.i
  denom.i <- sum(colSums(y_denom[, which(age >= ll_age & age <= ul_age)]))
  Prv_1549 <- sum.i / denom.i      
  
  res_all <- quantile(Prv_all, probs = c(0.5, 0.025, 0.975))
  res_1549 <- quantile(Prv_1549, probs = c(0.5, 0.025, 0.975))
    broad_prev <- rbind(res_1549, res_all)
  women_all <- quantile(wom_all, probs = c(0.5, 0.025, 0.975))
  women_1549 <- quantile(wom_1549, probs = c(0.5, 0.025, 0.975))
    broad_women <- rbind(women_1549, women_all)
  broad_prev <- data.frame(age_grp = c(paste(ll_age, "-", ul_age, sep = ""), "15+"),
                      median = broad_prev[, 1], lci = broad_prev[, 2], uci = broad_prev[, 3],
                      women = broad_women[, 1], women_lci = broad_women[, 2], women_uci = broad_women[, 3])
       
  age_grp_lb <- seq(15, 65, by = 5) 
  age_grp_ub <- c(seq(19, 64, by = 5), 104)
  
  prv_age_stat <- NULL
  pb <- txtProgressBar(1, length(age_grp_lb), style = 3)
  for (i in 1:length(age_grp_lb)) {
    sum_age.i <- colSums(y_pred[ , which(age >= age_grp_lb[i] & age <= age_grp_ub[i]), ])
    wom_age.i <- sum_age.i
    denom_age.i <- sum(y_denom[ , which(age >= age_grp_lb[i] & age <= age_grp_ub[i])])
    if (!is.null(nrow(sum_age.i))) {  sum_age.i <- colSums(sum_age.i) 
                                      wom_age.i <- sum_age.i }
    if (!is.null(nrow(denom_age.i))) { denom_age.i <- sum(denom_age.i) }
    prev <- quantile(sum_age.i / denom_age.i, probs = c(0.5, 0.025, 0.975))
    women <- quantile(wom_age.i,  probs = c(0.5, 0.025, 0.975))
    prev_age.i <- data.frame(age_grp = paste(age_grp_lb[i], "-", age_grp_ub[i], sep = ""),
                             median = prev[1], lci = prev[2], uci = prev[3],
                             women = women[1], women_lci = women[2], women_uci = women[3])
    prv_age_stat <- rbind(prv_age_stat, prev_age.i)
    setTxtProgressBar(pb, i)
  }
  close(pb)

  age_grp_lb <- seq(15, 55, by = 10) 
  age_grp_ub <- seq(24, 64, by = 10)

  prv_age_stat_10yr <- NULL
  pb <- txtProgressBar(1, length(age_grp_lb), style = 3)
  for (i in 1:length(age_grp_lb)) {
    sum_age.i <- colSums(y_pred[ , which(age >= age_grp_lb[i] & age <= age_grp_ub[i]), ])
    wom_age.i <- sum_age.i
    denom_age.i <- sum(y_denom[ , which(age >= age_grp_lb[i] & age <= age_grp_ub[i])])
    if (!is.null(nrow(sum_age.i))) {  sum_age.i <- colSums(sum_age.i) 
                                      wom_age.i <- sum_age.i }
    if (!is.null(nrow(denom_age.i))) { denom_age.i <- sum(denom_age.i) }
    prev <- quantile(sum_age.i / denom_age.i, probs = c(0.5, 0.025, 0.975))
    women <- quantile(wom_age.i,  probs = c(0.5, 0.025, 0.975))
    prev_age.i <- data.frame(age_grp = paste(age_grp_lb[i], "-", age_grp_ub[i], sep = ""),
                             median = prev[1], lci = prev[2], uci = prev[3],
                             women = women[1], women_lci = women[2], women_uci = women[3])
    prv_age_stat_10yr <- rbind(prv_age_stat_10yr, prev_age.i)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  val <- rbind(broad_prev, prv_age_stat, prv_age_stat_10yr)
  val$women <- round(val$women) * 1000
  val$women_lci <- round(val$women_lci) * 1000
  val$women_uci <- round(val$women_uci) * 1000
  rownames(val) <- NULL
    
  if (save_results == TRUE) {
    if (is.null(cnt_exclude)) {
        write.csv(val, file = paste("Estimates Global - ", outcome,
                                   " ", suffix, ".csv", sep = ""))
    } else {
    write.csv(val, file = paste("Estimates Global - excl ", cnt_exclude, " - ", outcome,
                                   " ", suffix, ".csv", sep = ""))      
    }}
   return(val)
}


#' Aggregate by country AND time
pool_pred_time <- function(res, SDat, IPV, Center_a = 30, Center_t = 2018, outcome,
                           save_results = TRUE, suffix = "", max_iter = 1e5, age_grp = "15-49",
                           extroplated_ever_sex_wgt_2018 = TRUE) {

  if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
    stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }

  data("denom_cnt_age_now")
    if (extroplated_ever_sex_wgt_2018 == FALSE) { 
      data("denom_ever_sex_2010")
      print("you are using 2010 ever had sex weights")
      denom_ever_sex <- denom_ever_sex_2010 }
    if (extroplated_ever_sex_wgt_2018 == TRUE) {  
      data("denom_ever_sex_2018") 
      print("you are using 2018 ever had sex weights")
      denom_ever_sex <- denom_ever_sex_2018 }
  data("wpp_std_now")

  n_cnt <- length(unique(IPV$cnt))
  n_reg <- length(unique(IPV$reg))
  n_sup <- length(unique(IPV$sup))
  
  gbd <- addition_to_cnt("gbd")
  who <- addition_to_cnt("who") 
  sdg <- addition_to_cnt("sdg")

  # Some countries don't have data on ever had sex, we take the world average
  EverSexAvg <- aggregate(denom_ever_sex$ever_had_sex,
                          by = list(denom_ever_sex$age_group), FUN = mean, na.rm = TRUE)
  colnames(EverSexAvg) <- c('age_group','ever_had_sex')

  if (!(age_grp %in% c("15+", "15-49"))) { 
      stop("specify correct age group") }
  if (age_grp == "15+") { 
      age <- seq(17.5, 104.5, 5) }
  if (age_grp == "15-49") { 
      age <- seq(17.5, 47.5, 5) }
  age65 <- ifelse(age > 65, 65, age)
  Xage <- splines::ns(age65 - Center_a,
                      knots = attr(SDat$Spl_W, "knots"),
                      Boundary.knots = attr(SDat$Spl_W, "Boundary.knots"))
  XPred <- NULL
  UIC <- unique(gbd$iso)
  LA <- length(age)

  SRLkUp <- as.data.frame(unique(IPV[c("sup", "Super")]))
      colnames(SRLkUp) <- c("SuperRegion", "Super")
  RLkUp <- as.data.frame(unique(IPV[c("reg", "Region")]))
      colnames(RLkUp) <- c("region", "Region")
  CLkUp <- as.data.frame(unique(IPV[c("cnt", "Country")]))
      colnames(CLkUp) <- c("iso3", "Country")
      
  pop_yr <- seq(2000, 2019, 1)
  # manipulating time trends   
  is_period <- length(attr(SDat$Spl_t, "period")) > 0 
  if (is_period) {
  XTime <- as.matrix(ifelse(pop_yr < attr(SDat$Spl_t, "cutoff"), 0, 1))
  }
  if (!is_period & length(attr(SDat$Spl_t, "knots")) > 0) {
  XTime <- as.matrix(splines::ns(pop_yr - Center_t,
                      knots = attr(SDat$Spl_t, "knots"),
                      Boundary.knots = attr(SDat$Spl_t, "Boundary.knots")))    
  } 
  if (!is_period & length(attr(SDat$Spl_t, "knots")) == 0) {
  XTime <- as.matrix(splines::ns(pop_yr - Center_t,
                      df = ncol(SDat$Spl_t),
                      Boundary.knots = attr(SDat$Spl_t, "Boundary.knots")))       
  }
  
  unique_yr <- unique(pop_yr)
  n_yr <- length(unique_yr)
  
  # We create the prediction matrix of Age * Country * Time
  pb <- txtProgressBar(1, length(UIC), style = 3)
  XPred <- NULL
    for (i in 1:length(UIC)) {
      for (j in 1:n_yr) {
      id_gbd_i <- which(as.character(gbd$iso) == as.character(UIC[i]))
      XPrd.ij <- data.frame(Age = age, Spl_a = Xage,
        Spl_t = do.call("rbind", replicate(LA, XTime[j, ], simplify = FALSE)),
        Time = pop_yr[j],
        ISO = rep(UIC[i], LA),
        GBD = rep(gbd$gbd[id_gbd_i], LA))
      XPrd.ij$SuperRegion <- SRLkUp$Super[as.character(SRLkUp$SuperRegion) == as.character(gbd$SuperRegion[id_gbd_i])]
      MissRegion <- RLkUp$Region[as.character(RLkUp$region) == as.character(gbd$Region[id_gbd_i])]
      XPrd.ij$Region <- ifelse(length(MissRegion) == 0, NA, MissRegion)
      MissCountry <- CLkUp$Country[as.character(CLkUp$iso3) == as.character(UIC[i])]
      XPrd.ij$Country <- ifelse(length(MissCountry) == 0, NA, MissCountry)
      XPred <- rbind(XPred, XPrd.ij)
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)
    XPred$GBD <- droplevels(XPred$GBD)

    if (any(names(XPred) == "X1")) { colnames(XPred)[which(names(XPred) == "X1")] <- "Spl_t.1" }
    if (any(names(XPred) == "Spl_t")) { colnames(XPred)[which(names(XPred) == "Spl_t")] <- "Spl_t.1" }

    # We predict using the individual iterations
    # To speed up, we only keep a max of max_iter = 25,000 iterations
    if (length(res[["glb"]]) > max_iter) { 
      n_iter <- max_iter
      print("Caution, you are not using all iterations - increase max_iter")
    } else {
      n_iter <- length(res[["glb"]])
    }
    iter_sel <- ceiling(seq(1, length(res[["glb"]]), length.out = n_iter))
    
    y_predt <- array(NA, dim = c(length(unique(XPred$ISO)),
                              length(unique(XPred$Age)), 
                              n_yr, 
                              n_iter))
    y_denomt <- array(NA, dim = c(length(unique(XPred$ISO)), 
                                  length(unique(XPred$Age))))
    XPred$ISO_New <- recode.cluster(XPred$ISO)
    YP_Index <- NULL
      
    
    pb <- txtProgressBar(1, nrow(XPred), style = 3)
    ind_spl_a <- names(XPred) %in% paste("Spl_a", rep(1:SDat$Ndof), sep = ".")
    ind_spl_t <- names(XPred) %in% paste("Spl_t", rep(1:SDat$Nt), sep = ".")
    for (i in 1:nrow(XPred)) {
        Region_NA <- XPred$Region[i]
        Country_NA <- XPred$Country[i]
        YP_Index.i <- data.frame(It = i, 
                                 Age = XPred$Age[i],
                                 GBD = XPred$GBD[i], 
                                 Country = XPred$ISO[i],
                                 Region = XPred$Region[i], 
                                 D1 = XPred$ISO_New[i],
                                 Time = XPred$Time[i])
        YP_Index <- rbind(YP_Index, YP_Index.i)
        # If we are predicting a country with data
        if (!is.na(XPred$Country[i])) {
          ind_a <- seq(XPred$Country[i], SDat$Ndof * n_cnt - n_cnt + XPred$Country[i], n_cnt)
          ind_t <- seq(XPred$Country[i], SDat$Nt * n_cnt - n_cnt + XPred$Country[i], n_cnt)
          y.i <- inv.logit(res[["cnt"]][iter_sel, XPred$Country[i]] + 
                  as.vector(res[["b_time_cnt"]][iter_sel, ind_t] %*% t(XPred[i, ind_spl_t])) +
                  as.vector((res[["b_spline_cnt"]][iter_sel, ind_a]) %*% t(XPred[i, ind_spl_a])))
        }
        # If we are predicting a new country or region
        if (is.na(Country_NA)) {
          # New country in a new region
          if (is.na(Region_NA)) {
            # intercepts
            var_new_cnt <- rnorm(n_iter, mean = 0, sd = res[["sd_cnt"]][iter_sel])
            var_new_reg <- rnorm(n_iter, mean = 0, sd = res[["sd_reg"]][iter_sel])
            # age
            ind_new_a <- seq(XPred$SuperRegion[i], SDat$Ndof * n_sup - n_sup + XPred$SuperRegion[i], n_sup)            
            var_new_spline_cnt <- matrix(NA, nrow = n_iter, ncol = SDat$Ndof)
            var_new_spline_reg <- matrix(NA, nrow = n_iter, ncol = SDat$Ndof)
            for (z in 1:SDat$Ndof) {
              var_new_spline_reg[, z] <- rnorm(n_iter, mean = 0, sd = res[["sd_spline_reg"]][iter_sel, z])
              var_new_spline_cnt[, z] <- rnorm(n_iter, mean = 0, sd = res[["sd_spline_cnt"]][iter_sel, z])
            }
            spl_new_reg <- res[["b_spline_sup"]][iter_sel, ind_new_a] + var_new_spline_reg + var_new_spline_cnt
           # time
           ind_new_t <- seq(XPred$SuperRegion[i], SDat$Nt * n_sup - n_sup + XPred$SuperRegion[i], n_sup)
           var_new_time_cnt <- matrix(0, nrow = n_iter, ncol = SDat$Nt)
           var_new_time_reg <- matrix(0, nrow = n_iter, ncol = SDat$Nt)
           for (z in 1:SDat$Nt) {
              var_new_time_reg[, z] <- rnorm(n_iter, mean = 0, sd = res[["sd_time_reg"]][iter_sel, z])
              var_new_time_cnt[, z] <- rnorm(n_iter, mean = 0, sd = res[["sd_time_cnt"]][iter_sel, z])
            }
          time_new_trend <- res[["b_time_sup"]][iter_sel, ind_new_t] + var_new_time_reg + var_new_time_cnt
          # outcome
          y.i <- inv.logit(res[["sup"]][iter_sel, XPred$SuperRegion[i]] + var_new_reg + var_new_cnt +
                   as.vector(time_new_trend %*% t(XPred[i, ind_spl_t])) + 
                   as.vector(spl_new_reg %*% t(XPred[i, ind_spl_a])))
          } 
          # New Country in a known region
          if (!is.na(Region_NA)) {
          # intercept  
          var_new_cnt <- rnorm(n_iter, mean = 0, sd = res[["sd_cnt"]][iter_sel])
          # age
          ind_new_a <- seq(XPred$Region[i], SDat$Ndof * n_reg - n_reg + XPred$Region[i], n_reg)
          var_new_spline_cnt <- matrix(NA, nrow = n_iter, ncol = SDat$Ndof)
          for (z in 1:SDat$Ndof) {
            var_new_spline_cnt[, z] <- rnorm(n_iter, mean = 0, sd = res[["sd_spline_cnt"]][iter_sel, z])
          }
          spl_new_reg <- res[["b_spline_reg"]][iter_sel, ind_new_a] + var_new_spline_cnt
          # time
          ind_new_t <- seq(XPred$Region[i], SDat$Nt * n_reg - n_reg + XPred$Region[i], n_reg)
          var_new_time_cnt <- matrix(NA, nrow = n_iter, ncol = SDat$Nt)
          for (z in 1:SDat$Nt) {
            var_new_time_cnt[, z] <- rnorm(n_iter, mean = 0, sd = res[["sd_time_cnt"]][iter_sel, z])
           }
          time_new_trend <- res[["b_time_reg"]][iter_sel, ind_new_t] + var_new_time_cnt
          # outcome
          y.i <- inv.logit(res[["reg"]][iter_sel, XPred$Region[i]] + var_new_cnt +
                     as.vector(time_new_trend %*% t(XPred[i, ind_spl_t])) + 
                     as.vector(res[["b_spline_reg"]][iter_sel, ind_new_a] %*% t(XPred[i, ind_spl_a])))
          }
        }
        # Adding the Population Weights to create a new numerator/denominator.
        # For IPV
        if (outcome != "NPSV") {
          W.i <- subset(denom_cnt_age_now, as.character(iso3) == as.character(XPred$ISO[i]))
          EverSex.i <- subset(denom_ever_sex, as.character(iso3) == as.character(XPred$ISO[i]))
          if (nrow(EverSex.i) == 0) { EverSex.i <- EverSexAvg } # There are some countries without data for Ever_Had_Sex.
          if (nrow(W.i) > 0) {
            closest_age <- which.min(abs(EverSex.i$age_group - XPred$Age[i])) # b/c ever sex is only up to age 45-49.
            w.i <- (W.i$women[W.i$Age == XPred$Age[i]] + 0.5) * EverSex.i$ever_had_sex[closest_age]
            y_predt[XPred$ISO_New[i], 
                    which(age == XPred$Age[i]), 
                    which(unique_yr == XPred$Time[i]), ] <- y.i * w.i
            y_denomt[XPred$ISO_New[i], 
                     which(age == XPred$Age[i])] <- w.i
          } else {  # Countries with a population less than 90K are not takent into account.
            W.i <- wpp_std_now$Pop[wpp_std_now$AgeM == XPred$Age[i]] / sum(wpp_std_now$Pop) * 22.5 
            closest_age <- which.min(abs(EverSex.i$age_group - XPred$Age[i])) # b/c ever sex is only up to age 45-49.
            w.i <- W.i * EverSex.i$ever_had_sex[closest_age]
            y_predt[XPred$ISO_New[i], 
                    which(age == XPred$Age[i]), 
                    which(unique_yr == XPred$Time[i]), ] <- y.i * w.i
            y_denomt[XPred$ISO_New[i], 
                     which(age == XPred$Age[i])] <- w.i
          } }
        # For NPSV
        if (outcome == "NPSV") {
          W.i <- subset(denom_cnt_age_now, as.character(iso3) == as.character(XPred$ISO[i]))
          if (dim(W.i)[1] > 0) {
            w.i <- W.i$women[W.i$Age == XPred$Age[i]] + 0.5
            y_predt[XPred$ISO_New[i], 
                    which(age == XPred$Age[i]), 
                    which(unique_yr == XPred$Time[i]), ] <- y.i * w.i
            y_denomt[XPred$ISO_New[i], 
                     which(age == XPred$Age[i])] <- w.i
          } else { # Countries with a population less than 90K are assumed to have the world standard pop.
            W.i <- wpp_std_now$Pop[wpp_std_now$AgeM == XPred$Age[i]] / sum(wpp_std_now$Pop) * 22.5 
            y_predt[XPred$ISO_New[i], 
                    which(age == XPred$Age[i]), 
                    which(unique_yr == XPred$Time[i]), ] <- y.i * w.i
            y_denomt[XPred$ISO_New[i], 
                     which(age == XPred$Age[i])] <- w.i
          } } 
        setTxtProgressBar(pb, i)
        }
      close(pb)
      
      # we add the index to find WHO and SDG regions
      YP_Index$WHO <- who$region_name[match(YP_Index$Country, who$iso3)]
      YP_Index$SDG <- sdg$subregion[match(YP_Index$Country, sdg$iso3)]      
      YP_Index$SDGsup <- sdg$region[match(YP_Index$Country, sdg$iso3)]   
      YP_Index$SDGdev <- sdg$least_developped[match(YP_Index$Country, sdg$iso3)]   
      
      # ---- Country Year ----
      print("processing country-year")
      # We only keep an index and aggregate over the specified age group
      YPInd <- YP_Index[!duplicated(YP_Index[, !(names(YP_Index) %in% c("It", "Age"))]), , ]
      #  Aggregate by Country-AgeGrp-Year
      AggCntTime <- NULL
      pb <- txtProgressBar(1, nrow(YPInd), style = 3)
      for (i in 1:nrow(YPInd)) {
        Prev <- NULL
        Pop.ij <- y_denomt[YPInd$D1[i], ]
        if (any(is.na(Pop.ij))) {
          next
        }
        Age.ij <- colSums(y_predt[YPInd$D1[i], , 
                                  which(unique_yr == YPInd$Time[i]), ])
        W.ij <- Age.ij / sum(Pop.ij)
        Prev.ij <- quantile(W.ij, probs = c(0.5, 0.025, 0.975))
        Prev <- rbind(Prev, Prev.ij)
        Agg.i <- data.frame(Median = Prev[, 1], LCI = Prev[, 2], UCI = Prev[, 3])
        Agg.i$Country <- YPInd$Country[i]
        Agg.i$GBD <- YPInd$GBD[i]
        Agg.i$Time <- YPInd$Time[i]
        if (is.na(YPInd$Country[i])) {
          Agg.i$ISO <- NA
        } else {
          Agg.i$ISO <- YPInd$Country[i]
        }
        AggCntTime <- rbind(AggCntTime, Agg.i)
        setTxtProgressBar(pb, i)
      }
      close(pb)

      AggCntTime$name <- countrycode::countrycode(sourcevar = AggCntTime$Country, 
                   origin = "iso3c", destination = "un.name.en", warn = TRUE,
                   custom_match = c(
                     "ASM" = "American Samoa",
                     "BMU" = "Bermuda",
                     "COK" = "Cook Islands",
                     "GRL" = "Greenland",
                     "GUM" = "Guam",
                     "HKG" = "Hong Kong (S.A.R. China)",
                     "MNP" = "Northern Mariana Islands",
                     "NIU" = "Niue",
                     "PRI" = "Puerto Rico",
                     "PSE" = "Palestinian Territory, Occupied",
                     "TWN" = "Taiwan (Republic of China)",
                     "VIR" = "Virgin Islands (USA)"))
          
    if (save_results == TRUE) {
      write.csv(AggCntTime, file = paste("Estimates Age-Cnt-Time - ", outcome, " ", suffix, ".csv", sep = ""))
    }

      # ---- GBD ----
      print("processing gbd")
      # we aggregate by GBD regions.
      YPInd <- YP_Index[!duplicated(YP_Index[, !(names(YP_Index) %in% c("It", "Age"))]), , ]
      YPgbd <- YP_Index[!duplicated(YP_Index[, !(names(YP_Index) %in% c("It", "Age", "Country", "Region", "D1", "WHO", "SDG", "SDGsup", "SDGdev"))]), , ]
      AggGBDTime <- NULL
      pb <- txtProgressBar(1, nrow(YPgbd), style = 3)
      for (i in 1:nrow(YPgbd)) {
        cnt_in_gbd <- unique(YPInd$D1[YPInd$GBD == YPgbd$GBD[i]])
        Pop.ij <- y_denomt[cnt_in_gbd, ]
        Age.ij <- colSums(colSums(y_predt[cnt_in_gbd, , 
                                  which(unique_yr == YPgbd$Time[i]), ]))
        W.ij <- Age.ij / sum(Pop.ij)
        Prev <- quantile(W.ij, probs = c(0.5, 0.025, 0.975))
        Agg.i <- data.frame(Median = Prev[1], LCI = Prev[2], UCI = Prev[3])
        Agg.i$GBD <- YPgbd$GBD[i]
        Agg.i$Time <- YPgbd$Time[i]
        AggGBDTime <- rbind(AggGBDTime, Agg.i)
        setTxtProgressBar(pb, i)
      }
      close(pb)
      rownames(AggGBDTime) <- NULL
      
     # change between 2000 and 2018 by GBD
      YPInd <- YP_Index[!duplicated(YP_Index[, !(names(YP_Index) %in% c("It", "Age"))]), , ]
      YPgbd <- YP_Index[!duplicated(YP_Index[, !(names(YP_Index) %in% c("It", "Age", "Country", "Region", "D1", "WHO", "SDG", "SDGsup", "SDGdev", "Time"))]), , ]
      gbd_time_change <- NULL
      for (i in 1:nrow(YPgbd)) {
        cnt_in_gbd <- unique(YPInd$D1[YPInd$GBD == YPgbd$GBD[i]])
        Pop.ij <- y_denomt[cnt_in_gbd, ]
        Prv1.ij <- colSums(colSums(y_predt[cnt_in_gbd, , 
                                  which(unique_yr == 2000), ])) / sum(Pop.ij)
        Prv2.ij <- colSums(colSums(y_predt[cnt_in_gbd, , 
                                  which(unique_yr == 2018), ])) / sum(Pop.ij)  
        W.ij <- Prv2.ij - Prv1.ij
        Change <- quantile(W.ij, probs = c(0.5, 0.025, 0.975))
        Agg.i <- data.frame(Median = Change[1], LCI = Change[2], UCI = Change[3])
        Agg.i$GBD <- YPgbd$GBD[i]
        gbd_time_change <- rbind(gbd_time_change, Agg.i)
      }
      rownames(gbd_time_change) <- NULL
      
      # ---- WHO ----
      # we aggregate by who regions.
      YPInd <- YP_Index[!duplicated(YP_Index[, !(names(YP_Index) %in% c("It", "Age"))]), , ]
      YPwho <- YP_Index[!duplicated(YP_Index[, !(names(YP_Index) %in% c("It", "Age", "Country", "Region", "D1", "GBD", "SDG", "SDGsup", "SDGdev"))]), , ]
      AggWHOTime <- NULL
      pb <- txtProgressBar(1, nrow(YPwho), style = 3)
      for (i in 1:nrow(YPwho)) {
        cnt_in_who <- unique(YPInd$D1[YPInd$WHO == YPwho$WHO[i]])
        Pop.ij <- y_denomt[cnt_in_who, ]
        Age.ij <- colSums(colSums(y_predt[cnt_in_who, , 
                                  which(unique_yr == YPwho$Time[i]), ]))
        W.ij <- Age.ij / sum(Pop.ij)
        Prev <- quantile(W.ij, probs = c(0.5, 0.025, 0.975))
        Agg.i <- data.frame(Median = Prev[1], LCI = Prev[2], UCI = Prev[3])
        Agg.i$WHO <- YPwho$WHO[i]
        Agg.i$Time <- YPwho$Time[i]
        AggWHOTime <- rbind(AggWHOTime, Agg.i)
        setTxtProgressBar(pb, i)
      }
      close(pb)      
      rownames(AggWHOTime) <- NULL
      
      # ---- SDG ----
      print("processing SDG")
      # we aggregate by sdg regions.
      YPInd <- YP_Index[!duplicated(YP_Index[, !(names(YP_Index) %in% c("It", "Age"))]), , ]
      YPsdg <- YP_Index[!duplicated(YP_Index[, !(names(YP_Index) %in% c("It", "Age", "Country", "Region", "D1", "GBD", "WHO", "SDGsup", "SDGdev"))]), , ]
      AggSDGTime <- NULL
      pb <- txtProgressBar(1, nrow(YPsdg), style = 3)
      for (i in 1:nrow(YPsdg)) {
        cnt_in_sdg <- unique(YPInd$D1[YPInd$SDG == YPsdg$SDG[i]])
        Pop.ij <- y_denomt[cnt_in_sdg, ]
        Age.ij <- colSums(colSums(y_predt[cnt_in_sdg, , 
                                  which(unique_yr == YPsdg$Time[i]), ]))
        W.ij <- Age.ij / sum(Pop.ij)
        Prev <- quantile(W.ij, probs = c(0.5, 0.025, 0.975))
        Agg.i <- data.frame(Median = Prev[1], LCI = Prev[2], UCI = Prev[3])
        Agg.i$WHO <- YPsdg$WHO[i]
        Agg.i$Time <- YPsdg$Time[i]
        AggSDGTime <- rbind(AggSDGTime, Agg.i)
        setTxtProgressBar(pb, i)
      }
      close(pb)       
      rownames(AggSDGTime) <- NULL
      
      # we aggregate by SDG SUPER regions.
      YPInd <- YP_Index[!duplicated(YP_Index[, !(names(YP_Index) %in% c("It", "Age"))]), , ]
      YPsdg <- YP_Index[!duplicated(YP_Index[, !(names(YP_Index) %in% c("It", "Age", "Country", "Region", "D1", "GBD", "WHO", "SDG", "SDGdev"))]), , ]
      AggSDGsupTime <- NULL
      pb <- txtProgressBar(1, nrow(YPsdg), style = 3)
      for (i in 1:nrow(YPsdg)) {
        cnt_in_sdg <- unique(YPInd$D1[YPInd$SDGsup == YPsdg$SDGsup[i]])
        Pop.ij <- y_denomt[cnt_in_sdg, ]
        Age.ij <- colSums(colSums(y_predt[cnt_in_sdg, , 
                                  which(unique_yr == YPsdg$Time[i]), ]))
        W.ij <- Age.ij / sum(Pop.ij)
        Prev <- quantile(W.ij, probs = c(0.5, 0.025, 0.975))
        Agg.i <- data.frame(Median = Prev[1], LCI = Prev[2], UCI = Prev[3])
        Agg.i$SDGsup <- YPsdg$SDGsup[i]
        Agg.i$Time <- YPsdg$Time[i]
        AggSDGsupTime <- rbind(AggSDGsupTime, Agg.i)
        setTxtProgressBar(pb, i)
      }
      close(pb)       
      rownames(AggSDGsupTime) <- NULL
      
      # we aggregate by least developped
      YPInd <- YP_Index[!duplicated(YP_Index[, !(names(YP_Index) %in% c("It", "Age"))]), , ]
      YPsdg <- YP_Index[!duplicated(YP_Index[, !(names(YP_Index) %in% c("It", "Age", "Country", "Region", "D1", "GBD", "WHO", "SDG", "SDGsup"))]), , ]
      AggSDGdevTime <- NULL
      pb <- txtProgressBar(1, nrow(YPsdg), style = 3)
      for (i in 1:nrow(YPsdg)) {
        cnt_in_sdg <- unique(YPInd$D1[YPInd$SDGdev == YPsdg$SDGdev[i]])
        Pop.ij <- y_denomt[cnt_in_sdg, ]
        Age.ij <- colSums(colSums(y_predt[cnt_in_sdg, , 
                                  which(unique_yr == YPsdg$Time[i]), ]))
        W.ij <- Age.ij / sum(Pop.ij)
        Prev <- quantile(W.ij, probs = c(0.5, 0.025, 0.975))
        Agg.i <- data.frame(Median = Prev[1], LCI = Prev[2], UCI = Prev[3])
        Agg.i$SDGdev <- YPsdg$SDGdev[i]
        Agg.i$Time <- YPsdg$Time[i]
        AggSDGdevTime <- rbind(AggSDGdevTime, Agg.i)
        setTxtProgressBar(pb, i)
      }
      close(pb)   
      rownames(AggSDGdevTime) <- NULL
      
      # ---- global time trend ----
        global_time <- NULL
        for (i in unique_yr) { 
        Pop.ij <- y_denomt[unique(YPInd$D1), ]
        Age.ij <- colSums(colSums(y_predt[unique(YPInd$D1), , 
                                  which(unique_yr == i), ]))
        W.ij <- Age.ij / sum(Pop.ij)
        Prev <- quantile(W.ij, probs = c(0.5, 0.025, 0.975))
        Agg.i <- data.frame(Median = Prev[1], LCI = Prev[2], UCI = Prev[3])
        Agg.i$Time <- i
        global_time <- rbind(global_time, Agg.i)
        }
        vaw_2000 <- colSums(colSums(y_predt[ , , which(unique_yr == 2000), ])) / sum(y_denomt)
        vaw_2018 <- colSums(colSums(y_predt[ , , which(unique_yr == 2018), ])) / sum(y_denomt) 
        change <-  quantile(vaw_2018 - vaw_2000, probs = c(0.5, 0.025, 0.975))
        change_2000_18 <- data.frame(Median = change[1], LCI = change[2], UCI = change[3])
        change_2000_18$Time <- as.character("2000 to 2018")
        global_time <- rbind(global_time, change_2000_18)
        rownames(global_time) <- NULL
        print(global_time)
        
    if (save_results == TRUE) {
      write.csv(AggGBDTime, file = paste("Estimates GBD-Time - ", outcome, " ", suffix, ".csv", sep = ""))
      write.csv(AggWHOTime, file = paste("Estimates WHO-Time - ", outcome, " ", suffix, ".csv", sep = ""))
      write.csv(AggSDGTime, file = paste("Estimates SDG-Time - ", outcome, " ", suffix, ".csv", sep = ""))
      write.csv(AggSDGTime, file = paste("Estimates SDG-Time - ", outcome, " ", suffix, ".csv", sep = ""))
      write.csv(AggSDGsupTime, file = paste("Estimates SDG Regions-Time - ", outcome, " ", suffix, ".csv", sep = ""))
      write.csv(AggSDGdevTime, file = paste("Estimates SDG Least-Developped-Time - ", outcome, " ", suffix, ".csv", sep = ""))
      write.csv(gbd_time_change, file = paste("Estimates Changes 2000-2018 by GBD - ", outcome, " ", suffix, ".csv", sep = ""))
      write.csv(global_time, file = paste("Estimates VAW 2000-2018 Global - ", outcome, " ", suffix, ".csv", sep = ""))
    }
      
      return(AggCntTime)
}


plot_agg_cnt_time <- function(agg_cnt_time, IPV, outcome, suffix = "") {

  if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
    stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }

  ToPlot <- NULL
  if (!("spouseonly" %in% names(IPV))) {
    IPV$spouseonly <- NA
  }
  IPV_sel <- subset(IPV, select = c("loage", "hiageI", "iso3", "country", "num_ess", "den_ess", "Time",
      "severe", "violence", "pstat", "spouseonly", "loqual", "viotime", "geo"))
  C_wdata <- unique(IPV_sel$iso3)
    
  for (i in 1:length(C_wdata)) {
    Cnt.i <- IPV_sel[IPV_sel$iso3 == C_wdata[i], ]
    Cnt.i$Duplicat <- paste(Cnt.i$Time, Cnt.i$severe, Cnt.i$violence, Cnt.i$pstat, Cnt.i$nointrain,
      Cnt.i$spouseonly, Cnt.i$loqual, Cnt.i$viotime, Cnt.i$geo)
    Unique.i <- unique(Cnt.i$Duplicat)
    for (j in 1:length(Unique.i)) {
      Cnt.ij <- subset(Cnt.i, Duplicat == Unique.i[j])
      Dim.ij <- nrow(Cnt.ij)
      Cnt1549 <- Cnt.ij[1, ]
      if (Dim.ij == 1) {
        CI.ij <- binom.test(round(Cnt1549$num_ess, 0), floor(Cnt1549$den_ess))$conf.int
        Cnt1549$Num <- Cnt1549$num_ess
        Cnt1549$Den <- Cnt1549$den_ess
        Cnt1549$Prv <- Cnt1549$Num / Cnt1549$Den
        Cnt1549$LCI <- CI.ij[1]
        Cnt1549$UCI <- CI.ij[2]
      }
      if (Dim.ij > 1) {
        Num.ij <- sum(Cnt.ij$num_ess[Cnt.ij$loage >= 15 & Cnt.ij$hiageI < 50])
        Den.ij <- sum(Cnt.ij$den_ess[Cnt.ij$loage >= 15 & Cnt.ij$hiageI < 50])
        if (Den.ij == 0) {
          Num.ij <- sum(Cnt.ij$num_ess[Cnt.ij$loage >= 15])
          Den.ij <- sum(Cnt.ij$den_ess[Cnt.ij$loage >= 15])
          Cnt1549$loage <- min(Cnt.ij$loage[Cnt.ij$loage >= 15])
          Cnt1549$hiageI <- max(Cnt.ij$hiageI)
          Cnt1549$Num <- Num.ij
          Cnt1549$Den <- Den.ij
          Cnt1549$Prv <- Cnt1549$Num / Cnt1549$Den
          CI.ij <- binom.test(round(Cnt1549$Num, 0), floor(Cnt1549$Den))$conf.int
          Cnt1549$LCI <- CI.ij[1]
          Cnt1549$UCI <- CI.ij[2]
        }
        if (Den.ij > 0) {
          Cnt1549$loage <- min(Cnt.ij$loage[Cnt.ij$loage >= 15])
          Cnt1549$hiageI <- max(Cnt.ij$hiageI[Cnt.ij$hiageI < 50])
          Cnt1549$Num <- Num.ij
          Cnt1549$Den <- Den.ij
          Cnt1549$Prv <- Cnt1549$Num / Cnt1549$Den
          CI.ij <- binom.test(round(Cnt1549$Num, 0), floor(Cnt1549$Den))$conf.int
          Cnt1549$LCI <- CI.ij[1]
          Cnt1549$UCI <- CI.ij[2]
        }}
      ToPlot <- rbind(ToPlot, Cnt1549)
    }}

  # We Graph
  Col <- c(rgb(20, 180, 235, 255, max = 255), rgb(134, 179, 0, 255, max = 255), rgb(204, 0, 67, 255, max = 255),
    rgb(255, 173, 51, 255, max = 255), rgb(255, 102, 53, 255, max = 255), rgb(135, 132, 148, 255, max = 255))
  ColT <- c(rgb(20, 180, 235, 50, max = 255), rgb(134, 179, 0, 50, max = 255), rgb(204, 0, 67, 50, max = 255),
    rgb(255, 173, 51, 50, max = 255), rgb(255, 102, 53, 50, max = 255), rgb(135, 132, 148, 50, max = 255))

  if (!is.null(dev.list())) { dev.off() }
  # All on the same page
  n_cnt <- length(C_wdata)
  if (outcome == "Ever IPV") { 
    a <- 4
      pdf(paste("Fit All - Ever IPV ", suffix, ".pdf", sep = ""),
          width = 14.5, height = 46)
      C_wdata <- unique(IPV_sel$iso3[order(IPV_sel$country)])
      length(C_wdata) + 1
      nf <- layout(matrix(c(1:(7 * ceiling((n_cnt + 1) / 7))), ncol = 7, nrow = ceiling((n_cnt + 1) / 7), byrow = TRUE),
                   widths = rep(lcm(5), 7), heights = rep(lcm(5), ceiling((n_cnt + 1) / 7)), respect = TRUE)
      par(mar = c(2, 1, 2, 1), oma = c(0, 1, 0, 0))
      plot(c(1) ~ 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Legend", xlim = c(0, 1), ylim = c(0, 1))
      legend(0, 0.9, legend = c("Ever-partnered", "All women", "Currently partnered",
                                "Nationally representative", "Sub-national", "All other types of estimate"),
             pch = c(1, 0, 2, NA, NA, 3), lty = c(NA, NA, NA, 1, 1, NA), lwd = c(NA, NA, NA, 4, 4, NA),
             col = c(rep("grey35", 3), Col[1], "darkseagreen4", "grey35"), bty = "n")
      text(0, 0.1, "Only physical and/or sexual \n IPV estimates are plotted \n (others are indicated by +)", pos = 4, cex = 0.9)

      for (i in 1:length(C_wdata)) {
        Data.all.i <- subset(ToPlot, iso3 == C_wdata[i])
        Data.i <- subset(Data.all.i, severe == "Not only severe violence" &
                           violence == "Physical and/or sexual IPV" & viotime == "Ever")
        Data_not_gs <- subset(Data.all.i, severe != "Not only severe violence" |
                           violence != "Physical and/or sexual IPV" | viotime != "Ever")
        Pred.i <- subset(agg_cnt_time, as.character(ISO) == as.character(C_wdata[i]))
        time_lim <- c(2000, 2019.25)
        country_title <- Data.all.i$country[1]
        if (country_title == "DEMOCRATIC REPUBLIC OF THE CONGO") { country_title <- "DRC" }
        if (country_title == "COTE D'IVOIRE") { country_title <- "CÔTE D'IVOIRE" }
        if (country_title == "CENTRAL AFRICAN REPUBLIC") { country_title <- "CAR" }
        plot(Pred.i$Median ~ Pred.i$Time,
             lwd = 2, type = "l", col = Col[a], xlim = time_lim, ylim = c(0, 1),
             xlab = "", ylab = "Prevalence Ever IPV", main = country_title)
        polygon( x = c(Pred.i$Time, rev(Pred.i$Time)),
          y = c(Pred.i$LCI, rev(Pred.i$UCI)), col = ColT[a], border = NA)

        if (nrow(Data_not_gs) > 0) { 
          for (j in 1:nrow(Data_not_gs)) { 
          col_sub <- ifelse(Data_not_gs$geo == "National", Col[1], "darkseagreen4") 
          points(Data_not_gs$Prv[j] ~ Data_not_gs$Time[j], pch = 3, col = col_sub)
          segments(x0 = Data_not_gs$Time[j], y0 = Data_not_gs$LCI[j], 
                   x1 = Data_not_gs$Time[j], y1 = Data_not_gs$UCI[j], 
                   col = col_sub, lwd = 1)
        }}
        
        if (dim(Data.i)[1] < 1) { next }
        for (j in 1:dim(Data.i)[1]) {
          pch.j <- ifelse(Data.i$pstat[j] == "Ever-partnered", 16,
                          ifelse(Data.i$pstat[j] == "All women", 15, 17))
          col.j <- ifelse(Data.i$geo[j] == "National", Col[1], "darkseagreen4")
          points(Data.i$Prv[j] ~ Data.i$Time[j], pch = pch.j, col = col.j, cex = 0.8)
          segments(x0 = Data.i$Time[j], y0 = Data.i$LCI[j], x1 = Data.i$Time[j], y1 = Data.i$UCI[j], col = col.j, lwd = 1)
          text(x = Data.i$Time[j], y = 0.01, paste(Data.i$loage[j], Data.i$hiageI[j], sep = "-"), cex = 0.5)
        } }
      dev.off() }

  if (outcome == "Past Year IPV") { 
    a <- 3
    pdf(paste("Fit All - Past Year IPV ", suffix, ".pdf", sep = ""),
        width = 14.5, height = 46)
    C_wdata <- unique(IPV_sel$iso3[order(IPV_sel$country)])
    length(C_wdata) + 1
    nf <- layout(matrix(c(1:(7 * ceiling((n_cnt + 1) / 7))), ncol = 7, 
                        nrow = ceiling((n_cnt + 1) / 7), byrow = TRUE),
      widths = rep(lcm(5), 7), heights = rep(lcm(5), ceiling((n_cnt + 1) / 7)), respect = TRUE)
    par(mar = c(2, 1, 2, 1), oma = c(0, 1, 0, 0))
    plot(c(1) ~ 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Legend", xlim = c(0, 1), ylim = c(0, 1))
    legend(0, 0.9,legend = c("Ever-partnered", "All women", "Currently partnered",
        "Nationally representative", "Sub-national", "All other types of estimate"),
      pch = c(1, 0, 2, NA, NA, 3), lty = c(NA, NA, NA, 1, 1, NA), lwd = c(NA, NA, NA, 4, 4, NA),
      col = c(rep("grey35", 3), Col[1], "darkseagreen4", "grey35"), bty = "n")
    text(0, 0.1, "Only physical and/or sexual \n IPV estimates are plotted \n (others are indicated by +)", pos = 4, cex = 0.9)

     for (i in 1:length(C_wdata)) {
        Data.all.i <- subset(ToPlot, iso3 == C_wdata[i])
        Data.i <- subset(Data.all.i, severe == "Not only severe violence" &
                          violence == "Physical and/or sexual IPV" & viotime == "Past yr")
        Data_not_gs <- subset(Data.all.i, severe != "Not only severe violence" |
                           violence != "Physical and/or sexual IPV" | viotime != "Past yr") 
      Pred.i <- subset(agg_cnt_time, as.character(ISO) == as.character(C_wdata[i]))
      time_lim <- c(2000, 2019.25)
      country_title <- Data.all.i$country[1]
      if (country_title == "DEMOCRATIC REPUBLIC OF THE CONGO") { country_title <- "DR CONGO" }
      if (country_title == "COTE D'IVOIRE") { country_title <- "CÔTE D'IVOIRE" }
      if (country_title == "CENTRAL AFRICAN REPUBLIC") { country_title <- "CAR" }
      plot(Pred.i$Median ~ Pred.i$Time,
        lwd = 2, type = "l", col = Col[a], xlim = time_lim, ylim = c(0, 1),
        xlab = "", ylab = "Prevalence Past Year IPV", main = country_title)
      polygon(x = c(Pred.i$Time, rev(Pred.i$Time)),
        y = c(Pred.i$LCI, rev(Pred.i$UCI)), col = ColT[a], border = NA)

      if (nrow(Data_not_gs) > 0) { 
          for (j in 1:nrow(Data_not_gs)) { 
          col_sub <- ifelse(Data_not_gs$geo == "National", Col[1], "darkseagreen4") 
          points(Data_not_gs$Prv[j] ~ Data_not_gs$Time[j], pch = 3, col = col_sub)
          segments(x0 = Data_not_gs$Time[j], y0 = Data_not_gs$LCI[j], 
                   x1 = Data_not_gs$Time[j], y1 = Data_not_gs$UCI[j], 
                   col = col_sub, lwd = 1)
          }}
 
        if (nrow(Data.i) < 1) { next }     
        for (j in 1:dim(Data.i)[1]) {
          pch.j <- ifelse(Data.i$pstat[j] == "Ever-partnered", 16,
                          ifelse(Data.i$pstat[j] == "All women", 15, 17))
          col.j <- ifelse(Data.i$geo[j] == "National", Col[1], "darkseagreen4")
          points(Data.i$Prv[j] ~ Data.i$Time[j], pch = pch.j, col = col.j, cex = 0.8)
          segments(x0 = Data.i$Time[j], y0 = Data.i$LCI[j], x1 = Data.i$Time[j], y1 = Data.i$UCI[j], col = col.j, lwd = 1)
          text(x = Data.i$Time[j], y = 0.01, paste(Data.i$loage[j], Data.i$hiageI[j], sep = "-"), cex = 0.5)
        } }
    dev.off() }

  if (outcome == "NPSV") { 
    a <- 2
    pdf(paste("Fit All - NPSV ", suffix, ".pdf", sep = ""),
      width = 14.5, height = 40)
    C_wdata <- unique(IPV_sel$iso3[order(IPV_sel$country)])
    length(C_wdata) + 1
    nf <- layout(matrix(c(1:(7 * ceiling((n_cnt + 1) / 7))), ncol = 7, nrow = ceiling((n_cnt + 1) / 7), byrow = TRUE),
                 widths = rep(lcm(5), 7), heights = rep(lcm(5), ceiling((n_cnt + 1) / 7)), respect = TRUE)
    par(mar = c(2, 1, 2, 1), oma = c(0, 1, 0, 0))
    plot(c(1) ~ 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Legend", xlim = c(0, 1), ylim = c(0, 1))
    legend(0, 0.9,
           legend = c("All severity (ever)", "All severity (past year)", "Severe (ever)", "Severe (past year)",
                      "Nationally representative", "Sub-national"),
           pch = c(16, 1, 17, 2, NA, NA), lty = c(NA, NA, NA, NA, 1, 1), lwd = c(NA, NA, NA, NA, 4, 4),
           col = c(rep("grey35", 4), Col[1], "darkorange2"), bty = "n")

  for (i in 1:length(C_wdata)) {
    Data.all.i <- subset(ToPlot, iso3 == C_wdata[i])
    Data.i <- subset(ToPlot, iso3 == C_wdata[i])
    Pred.i <- subset(agg_cnt_time, as.character(ISO) == as.character(C_wdata[i]))
    time_lim <- c(2000, 2019.25)
    country_title <- Data.all.i$country[1]
    if (country_title == "DEMOCRATIC REPUBLIC OF THE CONGO") { country_title <- "DRC" }
    if (country_title == "COTE D'IVOIRE") { country_title <- "CÔTE D'IVOIRE" }
    if (country_title == "CENTRAL AFRICAN REPUBLIC") { country_title <- "CAR" }
    plot(Pred.i$Median ~ Pred.i$Time,
      lwd = 2, type = "l", col = Col[a], xlim = time_lim, ylim = c(0, 1),
      xlab = "", ylab = "Prevalence NPSV", main = country_title)
    polygon(x = c(Pred.i$Time, rev(Pred.i$Time)),
        y = c(Pred.i$LCI, rev(Pred.i$UCI)), col = ColT[a], border = NA)
      if (nrow(Data.i) < 1) { next }
      for (j in 1:nrow(Data.i)) {
        pch.j <- ifelse(Data.i$severe[j] == "Not only severe violence" & Data.i$viotime[j] == "Ever", 16,
          ifelse(Data.i$severe[j] == "Not only severe violence" & Data.i$viotime[j] != "Ever", 1,
            ifelse(Data.i$severe[j] == "Only severe violence" & Data.i$viotime[j] == "Ever", 17, 2)))
        col.j <- ifelse(Data.i$geo[j] == "National", Col[1], "darkorange2")
        points(Data.i$Prv[j] ~ Data.i$Time[j], pch = pch.j, col = col.j, cex = 0.8)
        segments(x0 = Data.i$Time[j], y0 = Data.i$LCI[j], x1 = Data.i$Time[j], y1 = Data.i$UCI[j], col = col.j, lwd = 1)
        text(x = Data.i$Time[j], y = 0.01, paste(Data.i$loage[j], Data.i$hiageI[j], sep = "-"), cex = 0.5)
      }}

    dev.off() }
}

simpleCap <- function(x) {
      cap <- function(s) { 
            x <- tolower(s)
            x <- strsplit(x, " ")[[1]]
            val <- paste(toupper(substring(x, 1, 1)), substring(x, 2), sep = "", collapse = " ")
            return(val)
      }
      sapply(x, FUN = cap)
}

out_of_sample <- function(lprd_oos, SDat_oos, IPV, id_to_remove, cnt = FALSE, outcome = "Ever IPV", suffix = "", width = 21) {
  if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
    stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }
 
  # Pre-processing step
  Nb.iter <- dim(lprd_oos)[1]
  y_ppc <- NULL
  pb <- txtProgressBar(1, SDat_oos$N, style = 3)
  for (i in 1:SDat_oos$N) {
    p.i <- rbinom(Nb.iter, SDat_oos$Denom[i], lprd_oos[, i]) / SDat_oos$Denom[i]
    y_ppc.i <- quantile(p.i, probs = c(0.5, 0.025, 0.975))
    y_ppc <- rbind(y_ppc, y_ppc.i)
    setTxtProgressBar(pb, i)
  }
  close(pb)

  Prd <- data.frame(Median = y_ppc[, 1], LCI = y_ppc[, 2], UCI = y_ppc[, 3])
  if (cnt == TRUE) {
      rows_to_keep <- which(IPV$Country %in% id_to_remove) }
  if (cnt != TRUE) {
      rows_to_keep <- which(IPV$Study %in% id_to_remove) }
  Prd_oos <- Prd[rows_to_keep, ]
  IPV_oos <- IPV[rows_to_keep, ]
  Prd_oos$RowID <- seq(1, nrow(Prd_oos), 1)
  IPV_oos$RowID <- seq(1, nrow(IPV_oos), 1)  
  
  # Color Specification
  if (outcome == "Ever IPV") { a <- 6; b <- 4 }
  if (outcome == "Past Year IPV") { a <- 6; b <- 3 }
  if (outcome == "NPSV") { a <- 6; b <- 2 }
  Col <- c(rgb(20, 180, 235, 255, max = 255), rgb(134, 179, 0, 255, max = 255), rgb(204, 0, 67, 255, max = 255),
          rgb(255, 173, 51, 255, max = 255), rgb(255, 102, 53, 255, max = 255), rgb(135, 132, 148, 255, max = 255))
  ColT <- c(rgb(20, 180, 235, 150, max = 255), rgb(134, 179, 0, 150, max = 255), rgb(204, 0, 67, 150, max = 255),
          rgb(255, 173, 51, 150, max = 255), rgb(255, 102, 53, 150, max = 255), rgb(135, 132, 148, 150, max = 255))
  
  pdf(file = paste("Out-of-Sample - ", outcome, " ", suffix, ".pdf", sep = ""),
      width = width, height = 7.5)
    XLab <- IPV_oos[!duplicated(IPV_oos$iso3), ]
    xlim <- c(min(Prd_oos$RowID), 
              ifelse((max(Prd_oos$RowID) - min(Prd_oos$RowID)) < 10, 
                       min(Prd_oos$RowID) + 15.5, max(Prd_oos$RowID) + 0.5))
    plot(Prd_oos$Median ~ Prd_oos$RowID,
      type = "n", ylim = c(-0.25, 1.25), xlim = xlim,
      xlab = "", ylab = "Prevalence", main = paste(outcome, "- Out-of-sample", sep = " "),axes = FALSE)
    axis(1, at = c(min(Prd_oos$RowID), max(Prd_oos$RowID) + 1), labels = FALSE)
    axis(1, at = XLab$RowID, labels = XLab$iso3, las = 2, cex = 0.8)
    axis(2, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2))
    points(IPV_oos$Prv / 100 ~ IPV_oos$RowID, pch = 17, cex = 0.5, col = Col[a])
    points(I(Prd_oos$Median) ~ I(Prd_oos$RowID + 0.2), pch = 16, cex = 0.5, col = Col[b])
    for (j in 1:dim(Prd_oos)[1]) {
      segments(x0 = IPV_oos$RowID[j], y0 = IPV_oos$LCI[j], x1 = IPV_oos$RowID[j], y1 = IPV_oos$UCI[j], col = ColT[a])
      segments(x0 = Prd_oos$RowID[j] + 0.2, y0 = Prd_oos$LCI[j], x1 = Prd_oos$RowID[j] + 0.2, y1 = Prd_oos$UCI[j], col = ColT[b])
    }
    legend(x = min(IPV_oos$RowID), y = 1.25, pch = c(17, 16), col = c(Col[a], Col[b]), lwd = 2,
      legend = c("Data", "Out-of-Sample Predictions"), bty = "n")
    text(IPV_oos$RowID, y = -0.05, labels = paste(IPV_oos$loage, IPV_oos$hiage, sep = "-"), srt = 90, cex = 0.3)
    IPV_oos$StudyYr <- paste(IPV_oos$ID, IPV_oos$startyr, sep = "-")
      IPV_yr <- IPV_oos[!duplicated(IPV_oos$StudyYr), ]
      text(IPV_yr$RowID, y = -0.1, labels = IPV_yr$startyr, srt = 0, cex = 0.35)
    IPV_oos$VioType <- paste(IPV_oos$ID, IPV_oos$startyr, IPV_oos$violence)
      IPVTxtVioT <- IPV_oos[!duplicated(IPV_oos$VioType), ]
      text(IPVTxtVioT$RowID, y = 0.95, labels = as.character(IPVTxtVioT$violence), srt = 90, cex = 0.35)
    IPV_oos$StudySvr <- paste(IPV_oos$ID, IPV_oos$startyr, IPV_oos$violence, IPV_oos$Severe, sep = "-")
      IPV_Severe <- IPV_oos[!duplicated(IPV_oos$StudySvr), ]
      text(IPV_Severe$RowID, y = -0.25, labels = IPV_Severe$Severe, srt = 90, cex = 0.35)
    IPV_oos$StudyPstat <- paste(IPV_oos$ID, IPV_oos$startyr, IPV_oos$violence, IPV_oos$Severe, IPV_oos$pstat, sep = "-")
      IPV_pstat <- IPV_oos[!duplicated(IPV_oos$StudyPstat), ]
      text(IPV_pstat$RowID, y = 0.75, labels = IPV_pstat$pstat, srt = 90, cex = 0.35)
    IPV_oos$GEO <- paste(IPV_oos$ID, IPV_oos$startyr, IPV_oos$violence, IPV_oos$Severe, IPV_oos$pstat, IPV_oos$geo, sep = "-")
      IPV_GEO <- IPV_oos[!duplicated(IPV_oos$GEO), ]
      text(IPV_GEO$RowID, y = -0.15, labels = IPV_GEO$geo, srt = 90, cex = 0.35)
    # This line is for past year IPV only
    IPV_oos$PastYr <- paste(IPV_oos$ID, IPV_oos$startyr, IPV_oos$violence, IPV_oos$Severe, IPV_oos$pstat, IPV_oos$geo, IPV_oos$viotime, sep = "-")
      IPV_PastYr <- IPV_oos[!duplicated(IPV_oos$PastYr), ]
      text(IPV_PastYr$RowID, y = -0.2, labels = IPV_PastYr$viotime, srt = 0, cex = 0.35)
    dev.off()
  
  me <- median(IPV_oos$Prv/100 - Prd_oos$Median)
  mae <- median(abs(IPV_oos$Prv/100 - Prd_oos$Median))
  mre <- median((IPV_oos$Prv/100 - Prd_oos$Median) / (IPV_oos$Prv/100), na.rm = TRUE)
  mare <- median(abs((IPV_oos$Prv/100 - Prd_oos$Median)) / (IPV_oos$Prv/100) , na.rm = TRUE)
  cov <- ifelse(IPV_oos$Prv/100 >= Prd_oos$LCI & IPV_oos$Prv/100 <= Prd_oos$UCI, 1, 0)
  li <- ifelse(IPV_oos$Prv/100 < Prd_oos$LCI, 1, 0)
  ui <- ifelse(IPV_oos$Prv/100 > Prd_oos$UCI, 1, 0)
  stats <- data.frame(region = "Overall", n_excl = length(id_to_remove), n_obs = nrow(IPV_oos), 
                                  me = me, mae = mae, mre = mre, mare = mare,
                                  lower = mean(li), upper = mean(ui), cov = mean(cov))
  stats$me <- stats$me * 100
  stats$mae <- stats$mae * 100
  stats$mre <- stats$mre * 100                    
  stats$mare <- stats$mare * 100                      
  stats$lower <- stats$lower * 100
  stats$upper <- stats$upper * 100
  stats$cov <- stats$cov * 100
  write.csv(stats, paste("Out-of-samples validation - ", outcome, " ", suffix, ".csv", sep = ""))
  return(stats)
}


heatmap_vaw <- function(agg_cnt, cnt_var = "cnt", outcome = "", suffix = "", 
                        width = 8.5, height = 11.5) {
    require(magrittr)
    if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
    stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }
  
  # if some countries have NA name, we remove
  agg_cnt <- na.omit(agg_cnt)
  agg_cnt$cnt <- agg_cnt[, names(agg_cnt) %in% cnt_var]
  agg_cnt$cnt[agg_cnt$cnt == "Côte D’ivoire"] <- "Cote d'Ivoire"
  agg_cnt$median <- agg_cnt$Median * 100
  agg_cnt$agerange <- paste(agg_cnt$loage, agg_cnt$hiage, sep = "-")
  agg_cnt$agerange[agg_cnt$agerange == "75-99"] <- "75+"
  agg <- subset(agg_cnt, !(agerange %in% c("15-49","15-24","25-34", "35-49","50-99")))
  
  if (outcome == "Ever IPV") {
      agg <- agg %>%
      # convert state to factor and reverse order of levels
      plyr::mutate(cnt = factor(cnt, levels = rev(sort(unique(cnt))))) %>%
      # create a new variable from count
      plyr::mutate(median_cat = cut(median, 
                   breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100),
                   labels = c("0-4%", "5-9%", "10-14%", "15-19%", "20-24%", "25-29%", "30-34%", "35-39%", "40-44%", "45-49%",
                              "50-59%", "60-69%", "70-79%", "80-100%"))) %>%
      # change level order
      plyr::mutate(median_cat = factor(as.character(median_cat), levels = rev(levels(median_cat))))
  } else {
      agg <- agg %>%
      # convert state to factor and reverse order of levels
      plyr::mutate(cnt = factor(cnt, levels = rev(sort(unique(cnt))))) %>%
      # create a new variable from count
      plyr::mutate(median_cat = cut(median, 
                   breaks = c(0, 2.5, 5, 7.5, 10, 15, 20, 25, 100),
                   labels = c("<2.5%", "2.5-4.9%", "5.0-7.4%", "7.5-9.9%", "10-14.9%", "15.0-19.9%", "20.0-24.9%", "25.0-100%"))) %>%
      # change level order
      plyr::mutate(median_cat = factor(as.character(median_cat), levels = rev(levels(median_cat)))) 
  }
  
   pal <- rev(wesanderson::wes_palette("Zissou1", length(unique(agg$median_cat)), type = "continuous"))
  # pal <- wesanderson::wes_palette("Zissou1", 5, type = "continuous")
  
  pp <- ggplot2::ggplot(agg, ggplot2::aes(x = agerange, y = cnt, fill = median_cat)) +
    # add border white colour of line thickness 0.25
       ggplot2::geom_tile(colour = "white", size = 0.5) +
       ggplot2::facet_grid(GBD ~ ., scales = "free_y", space = "free_y", switch = "x") +
       ggplot2::guides(fill = ggplot2::guide_legend(title = "Prevalence (%)")) +
    # remove x and y axis labels
       ggplot2::labs(x = "", y = "", title = outcome) +
    # remove extra space
       ggplot2::scale_y_discrete(expand = c(0, 0)) +
    # set a base size for all fonts
      ggplot2::theme_grey(base_size = 8) +
      ggplot2::scale_fill_manual(values = pal, na.value = "grey90") +
      # ggplot2::scale_fill_gradientn(colours = pal, na.value = "grey90") +
    # theme options
      ggplot2::theme(
        strip.text.y = ggplot2::element_text(angle = 0),
        # bold font for legend text
        legend.text = ggplot2::element_text(face = "bold"),
        # set thickness of axis ticks
        axis.ticks = ggplot2::element_line(size = 0.4),
        # remove plot background
        plot.background = ggplot2::element_blank(),
        # remove plot border
        panel.border = ggplot2::element_blank(),
        # strip facet backgroud
        strip.background = ggplot2::element_rect(fill = "grey95"))

    # Watermak
    pp <- pp + ggplot2::annotate("text", x = Inf, y = -Inf, label = "Preliminary",
            hjust = 1, vjust = -1, angle = 0, col = rgb(1, 0, 0, 0.2), cex = 9,
            fontface = "bold", alpha = 0.25)
    
   ggplot2::ggsave(pp, filename = paste("Heatmap - ", outcome, " ", suffix, ".pdf", sep = ""),
           device = "pdf",  width = width, height = height)

  print(pp)
}



#' Wolrd Map Results.
world_map_prv <- function(agg_cnt, ipv, age_grp = "15-49", outcome = "Ever IPV", 
                          width = 11, height = 8.5, suffix = "", watermark = FALSE,
                          new_proj = TRUE) {
    require(magrittr)
    if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
    stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }
  
    if (new_proj == FALSE) {
      data("world_map")
      data("who_lines")
      data("who_poly") 
    }
    if (new_proj == TRUE) {
      data("world_map_proj")
      data("who_lines_proj")
      data("who_poly_proj") 
      world_map <- world_map_proj
      who_lines <- who_lines_proj
      who_poly <- who_poly_proj      
    }  

    world_map$ISO_A3 <- world_map$iso3
    map <- world_map
    
    agg_cnt$median <- agg_cnt$Median * 100
    agg_cnt$agerange <- paste(agg_cnt$loage, agg_cnt$hiage, sep = "-")
    agg <- subset(agg_cnt, agerange == age_grp)
    agg$cnt <- agg$Country
    
      #agg_qtl <- quantile(agg$median, probs = seq(0.1, 1, by = 0.1))
      #agg_qtl_name <- round(agg_qtl, 0)
      #agg_qtl_name <- c(paste(0, "-", agg_qtl_name[1] - 0.1, "%", sep = ""),
      #                  paste(agg_qtl_name[1:(length(agg_qtl_name) - 1)], "-", agg_qtl_name[-1] - 0.1, "%", sep = ""))
      agg_qtl <- unique(c(seq(10, 40, by = 5), 100))
      agg_qtl_name <- agg_qtl
      agg_qtl_name <- c(paste(0, "-", agg_qtl_name[1] - 1, "%", sep = ""),
                        paste(agg_qtl_name[1:(length(agg_qtl_name) - 2)], "-", agg_qtl_name[-c(1, length(agg_qtl_name))] - 1, "%", sep = ""),
                        paste("\u2265", agg_qtl_name[(length(agg_qtl_name) - 1)], "%", sep =""))
      agg <- agg %>%
      # convert state to factor and reverse order of levels
      plyr::mutate(cnt = factor(cnt, levels = rev(sort(unique(cnt))))) %>%
      # create a new variable from count
      plyr::mutate(median_cat = cut(median, 
                   breaks = c(0, agg_qtl),
                   labels = agg_qtl_name)) %>%
      # change level order
      plyr::mutate(median_cat = factor(as.character(median_cat), levels = rev(levels(median_cat)))) 

  if (outcome == "Past Year IPV") {
      agg_qtl <- c(seq(5, 25, by = 5), 100)
      agg_qtl_name <- agg_qtl
      agg_qtl_name <- c(paste(0, "-", agg_qtl_name[1] - 1, "%", sep = ""),
                        paste(agg_qtl_name[1:(length(agg_qtl_name) - 2)], "-", agg_qtl_name[-c(1, length(agg_qtl_name))] - 1, "%", sep = ""),
                        paste("\u2265", agg_qtl_name[(length(agg_qtl_name) - 1)], "%", sep =""))
      agg <- agg %>%
      # convert state to factor and reverse order of levels
      plyr::mutate(cnt = factor(cnt, levels = rev(sort(unique(cnt))))) %>%
      # create a new variable from count
      plyr::mutate(median_cat = cut(median, 
                   breaks = c(0, agg_qtl),
                   labels = agg_qtl_name)) %>%
      # change level order
      plyr::mutate(median_cat = factor(as.character(median_cat), levels = rev(levels(median_cat)))) 
    
  }      
      
  if (outcome == "NPSV") {
      agg <- agg %>%
      # convert state to factor and reverse order of levels
      plyr::mutate(cnt = factor(cnt, levels = rev(sort(unique(cnt))))) %>%
      # create a new variable from count
      plyr::mutate(median_cat = cut(median, 
                   breaks = c(0, 2, 4, 6, 9, 12, 100),
                   labels = c("<2%", "2-3%", "4-5%", "6-8%", "9-11%",
                              "\u2265 12%"))) %>%
      # change level order (\u2265 not working for bigger equal)
      plyr::mutate(median_cat = factor(as.character(median_cat), levels = rev(levels(median_cat)))) 
  }
  
    # Removing countries with no data
    agg$with_data <- ifelse(as.character(agg$cnt) %in% unique(ipv$iso3), 1, 0)
    agg$median_cat <- ordered(agg$median_cat, levels = rev(levels(agg$median_cat)))
    levels(agg$median_cat) <- c(levels(agg$median_cat), "No data")
    agg$median_cat[agg$with_data == 0] <- "No data"

    # We plot
    pal <- (wesanderson::wes_palette("Zissou1", length(unique(agg$median_cat)) - 1, type = "continuous"))
    pal <- c(pal, "grey85")
    
    gg <- ggplot2::ggplot()
    gg <- gg + ggplot2::geom_map(data = map, map = map,
                        fill = "#ffffff", color = NA,
                        ggplot2::aes(x = long, y = lat, map_id = id, group = group))
    gg <- gg + ggplot2::geom_map(data = agg, map = map, color = "white", size = 0.15,
                        ggplot2::aes(fill = median_cat, group = Country, map_id = Country))
    gg <- gg + ggplot2::scale_fill_manual(values = pal, name = "Median prevalence (%)")
    gg <- gg + ggplot2::geom_map(data = who_poly, map = who_poly, fill = "white", color = NA,
                        ggplot2::aes(x = long, y = lat, map_id = id))
    gg <- gg + ggplot2::geom_path(data = who_lines, color = "grey85", size = 0.35, 
                        ggplot2::aes(x = long, y = lat, group = as.numeric(group)))
    gg <- gg + ggplot2::labs(x = '', y = '', title = paste(outcome, " (", age_grp, ")", sep = ""))
    gg <- gg + ggplot2::coord_equal(ratio = 1)
    gg <- gg + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                  axis.text.y = ggplot2::element_blank(),
                  axis.ticks = ggplot2::element_blank(),
                  rect = ggplot2::element_blank()) +
               ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    gg <- gg + ggplot2::theme(legend.position = "right")
    gg <- gg + ggplot2::theme(plot.title = ggplot2::element_text(size = 18))
    
    # Watermak
    if (watermark == TRUE) {
    gg <- gg + ggplot2::annotate("text", x = median(map$long), y = median(map$lat), label = "Preliminary",
            angle = 45, col = rgb(1, 0, 0, 0.2), cex = 9,
            fontface = "bold", alpha = 0.25)
    }
    
    #require(Cairo)
    ggplot2::ggsave(filename = paste("Map. Prevalence - ", outcome, " ", suffix, ".png", sep = ""),
           plot = gg, device = "png", dpi = 500, units = "in",
           width = width, height = height)

    print(gg) 
}


prv_plot_cnt <- function(agg_cnt, IPV, age_grp = "15-49", outcome = "Ever IPV", width = 8.5, height = 19,
                        suffix = "", xlim = c(0, 60)) {
    if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
    stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }
  
    agg_cnt$median <- agg_cnt$Median * 100
    agg_cnt$lci <- agg_cnt$LCI * 100
    agg_cnt$uci <- agg_cnt$UCI *100
    agg_cnt$agerange <- paste(agg_cnt$loage, agg_cnt$hiage, sep = "-")
    agg <- subset(agg_cnt, agerange == age_grp)
    agg <- na.omit(agg)
    agg$region <- recode_gbd(agg$GBD)
    agg$sr <- class_reg_to_superreg(agg$region)

    to_insert <- agg[!duplicated(agg$region), ]
    row_insert <- to_insert
    row_insert[,] <- NA
    row_insert$sr <- NA
    row_insert$region <- to_insert$region
    row_insert$sr <- class_reg_to_superreg(row_insert$region)
    row_insert$order <- 0
    row_insert$gbd_label <- to_insert$region
    row_insert$cnt <- "A"
    
    agg$order <- 1
    agg$gbd_label <- NA
    agg <- rbind(agg, row_insert)
    agg$cnt[!is.na(agg$gbd_label)] <- NA
    agg <- agg[ order(as.character(agg$region), as.character(agg$order), as.character(agg$cnt)), ]
    agg$row_id <- rev(seq(1:nrow(agg)))
    
    agg$with_data <- ifelse(as.character(agg$Country) %in% unique(IPV$iso3), 1, 0)
    
    col <- wesanderson::wes_palette(name = "FantasticFox1", n = 5, type = 'discrete')

    if (!is.null(dev.list())) { dev.off() }
    par(mar = c(2, 7, 2, 2), oma = c(2, 5, 2, 2))
    plot(agg$row_id ~ agg$median, xlim = xlim, xlab = "", ylab = "", 
         axes = FALSE, pch = "")
    # Watermark
      text(x = grconvertX(0.5, from = "npc"),  # align to center of plot X axis
       y = grconvertY(0.5, from = "npc"),      # align to center of plot Y axis
        labels = "Preliminary",                # our watermark
        cex = 3, font = 2,                     # large, bold font - hard to miss
        col = rgb(1, 0, 0, 0.2),               # translucent (0.2 = 20%) red color
        srt = 45) 
      
    mtext(text = paste(outcome, " (", age_grp, ")", sep = ""), side = 3 ,line = -2, font = 2)
    mtext(text = "Prevalence (%)", side = 1, line = 0, cex = 0.75)
      axis(1, at = seq(0, xlim[2], 5), labels = seq(0, xlim[2], 5), cex.axis = 0.7, line = -2, tck = -0.01)
      axis(2, at = agg$row_id, labels = agg$cnt, lwd = 0, las = 1, cex.axis = 0.5, line = 0)
      axis(2, at = agg$row_id, labels = agg$gbd_label, lwd = 0, las = 1, cex.axis = 0.5, line = 0, font = 2, adj = 1) 
    # countries with data
    agg1 <- agg[agg$with_data == 1, ]
    segments(x0 = agg1$lci, x1 = agg1$uci, y0 = agg1$row_id, y1 = agg1$row_id, lwd = 1.5, col = col[3])
    points(agg1$row_id ~ agg1$median, pch = 15, col = col[3], cex = 0.5)
    # countries without data
    agg0 <- agg[agg$with_data == 0, ]
    segments(x0 = agg0$lci, x1 = agg0$uci, y0 = agg0$row_id, y1 = agg0$row_id, lwd = 1.5, col = col[2])
    points(agg0$row_id ~ agg0$median, pch = 16, col = col[1], cex = 0.5)
    # legend
    legend(x = xlim[2] * 0.6, y = max(agg$row_id), legend = c("Countries with data", "Countries without data"),
           pch = c(15, 16), lwd = 1, col = c(col[3], col[2]),
           bty = "n", cex = 0.75)
    dev.copy(pdf, paste("Prevalence CNT - ", outcome, " ", suffix, ".pdf", sep = ""), width = width, height = height)
    dev.off()
}


prv_plot_reg <- function(agg_reg, IPV, age_grp = "15-49", outcome = "Ever IPV", width = 8.5, height = 6,
                         suffix = "", xlim = c(0, 50)) {
    if (!(outcome %in% c("Ever IPV", "Past Year IPV", "NPSV"))) {
    stop("Outcome must be 'Ever IPV', 'Past Year IPV', or 'NPSV'") }
  
    agg_reg$median <- agg_reg$Median * 100
    agg_reg$lci <- agg_reg$LCI * 100
    agg_reg$uci <- agg_reg$UCI *100
    agg <- agg_reg
    agg <- na.omit(agg)
    agg$region <- recode_gbd(agg$GBD)
    agg$sr <- class_reg_to_superreg(agg$region)
      
    to_insert <- agg[!duplicated(agg$sr), ]
    row_insert <- to_insert
    row_insert[,] <- NA
    row_insert$sr <- to_insert$sr
    row_insert$region <- NA
    row_insert$order <- 0
    row_insert$gbd_label <- to_insert$sr
    agg$order <- 1
    agg$gbd_label <- NA
    agg <- rbind(agg, row_insert)
    agg$cnt[!is.na(agg$gbd_label)] <- NA
    
    agg <- agg[order(as.character(agg$sr), agg$order, as.character(agg$region)), ]
    agg$row_id <- rev(seq(1:nrow(agg)))
    
    agg$with_data <- ifelse(as.character(agg$region) %in% unique(IPV$region), 1, 0)
    
    col <- wesanderson::wes_palette(name = "FantasticFox1", n = 5, type = 'discrete')
    
    if (!is.null(dev.list())) { dev.off() }
    par(mar = c(2, 7, 2, 2), oma = c(2, 5, 2, 2))
    plot(agg$row_id ~ agg$median, xlim = xlim, xlab = "", ylab = "", 
         axes = FALSE, pch = "")
        # Watermark
      text(x = grconvertX(0.5, from = "npc"),  # align to center of plot X axis
       y = grconvertY(0.5, from = "npc"),      # align to center of plot Y axis
        labels = "Preliminary",                # our watermark
        cex = 3, font = 2,                     # large, bold font - hard to miss
        col = rgb(1, 0, 0, 0.2),               # translucent (0.2 = 20%) red color
        srt = 45) 
      
    mtext(text = paste(outcome, " (", age_grp, ")", sep = ""), side = 3 ,line = 0, font = 2)
    mtext(text = "Prevalence (%)", side = 1, line = 2, cex = 0.75)
      axis(1, at = seq(0, xlim[2], 5), labels = seq(0, xlim[2], 5), cex.axis = 0.7, line = 0, tck = -0.02)
      axis(2, at = agg$row_id, labels = agg$region, lwd = 0, las = 1, cex.axis = 0.5, line = 0)
      axis(2, at = agg$row_id, labels = agg$gbd_label, lwd = 0, las = 1, cex.axis = 0.5, line = 0, font = 2, adj = 1) 
    # countries with data
    agg1 <- agg[agg$with_data == 1, ]
    segments(x0 = agg1$lci, x1 = agg1$uci, y0 = agg1$row_id, y1 = agg1$row_id, lwd = 1.5, col = col[3])
    points(agg1$row_id ~ agg1$median, pch = 15, col = col[3], cex = 0.5)
    # countries without data
    agg0 <- agg[agg$with_data == 0, ]
    segments(x0 = agg0$lci, x1 = agg0$uci, y0 = agg0$row_id, y1 = agg0$row_id, lwd = 1.5, col = col[2])
    points(agg0$row_id ~ agg0$median, pch = 16, col = col[1], cex = 0.5)
    # legend
    legend(x = xlim[2] * 0.65, y = max(agg$row_id) - 0.5, legend = c("Regions with data", "Regions without data"),
           pch = c(15, 16), lwd = 1, col = c(col[3], col[2]),
           bty = "n", cex = 0.7)
    dev.copy(pdf, paste("Prevalence REG - ", outcome, " ", suffix, ".pdf", sep = ""), width = width, height = height)
    dev.off()
}

plot_global_trend <- function(ever_trend, past_trend, width = 6, height = 5, ylim = c(0, 50)) {

  ever_trend$outcome <- "ever ipv"
  past_trend$outcome <- "past ipv"
  
  trend_ <- rbind(ever_trend, past_trend)
  slope_id <- which(trend_$Time == "2000 to 2018")
  slope <- trend_[slope_id, ]
  print(slope)
  est_2000_2018 <- which(trend_$Time != "2000 to 2018" & trend_$Time != "2019")
  trend <- trend_[est_2000_2018, ]
  trend$year <- as.numeric(as.character(trend$Time))
  trend$median <- trend$Median * 100
  trend$lci <- trend$LCI * 100
  trend$uci <- trend$UCI * 100
  trend_ever <- subset(trend, outcome == "ever ipv")
  trend_past <- subset(trend, outcome == "past ipv")
  
  png("global_time_trend.png", width, height, unit = "in", res = 400)
  par(mfrow = c(1, 1), family = "sans", yaxs = "i", yaxs = "i", oma = c(0, 0, 0, 0), mar = c(4, 4, 2, 0))
  
  plot(trend_past$median ~ trend_past$year, type = 'n', axes = FALSE, yaxs = "i", xaxs = "i",
       xlim = c(2000, 2019), ylim = ylim, xlab = "Year", ylab = "Standardized prevalence (%)")
  axis(1, at = seq(2000, 2020, by = 5), labels = seq(2000, 2020, by = 5), las = 1)
  axis(2, at = seq(0, ylim[2], by = 10), labels = seq(0, ylim[2], by = 10), las = 1)
  #abline(h = seq(0, 100, by = 10), col = "lightgray", lwd = 0.5)
      # past
      polygon(x = c(trend_past$year, rev(trend_past$year)),
        y = c(trend_past$lci, rev(trend_past$uci)),
        col = scales::alpha("#E44253", 0.5), border = NA)
      lines(trend_past$median ~ trend_past$year, lwd = 1, col = scales::alpha("#E44253", 1))
      # ever
      polygon(x = c(trend_ever$year, rev(trend_ever$year)),
        y = c(trend_ever$lci, rev(trend_ever$uci)),
        col = scales::alpha("#90C4AA", 0.5), border = NA)
      lines(trend_ever$median ~ trend_ever$year, lwd = 1, col = scales::alpha("#90C4AA", 1))
     
      # segments(x0 = trend_past$year + 0.1, y0 = trend_past$lci,
      #          x1 = trend_past$year + 0.1, y1 = trend_past$uci,
      #          col = scales::alpha("#E44253", 0.75), lwd = 2)
      # points(trend_past$median ~ I(trend_past$year + 0.1), pch = 19, cex = 0.6,
      #        col = scales::alpha("#E44253", 1))
      # # ever
      # segments(x0 = trend_ever$year - 0.1, y0 = trend_ever$lci,
      #          x1 = trend_ever$year - 0.1, y1 = trend_ever$uci,
      #          col = scales::alpha("#90C4AA", 0.75), lwd = 2)
      # points(trend_ever$median ~ I(trend_ever$year - 0.1), pch = 19, cex = 0.6,
      #       col = scales::alpha("#90C4AA", 1))
      
  legend("topright", lwd = 3, col =  c(scales::alpha("#90C4AA", 1), scales::alpha("#E44253", 1)),
         legend = c("Lifetime", "Past year"), bty = "n")
  dev.off()
  
}

plot_global_age <- function(ever_age, past_age, width = 3, height = 4, ylim = c(0, 50)) {
  ever_age$outcome <- "ever ipv"
  past_age$outcome <- "past ipv"
  
  trend_ <- rbind(ever_age, past_age)
  trend_$age <- ifelse(trend_$age_grp == "15-19", 17.5,
                ifelse(trend_$age_grp == "20-24", 22.5,
                ifelse(trend_$age_grp == "25-29", 27.5,
                ifelse(trend_$age_grp == "30-34", 32.5,
                ifelse(trend_$age_grp == "35-39", 37.5,
                ifelse(trend_$age_grp == "40-44", 42.5,
                ifelse(trend_$age_grp == "45-49", 47.5,
                ifelse(trend_$age_grp == "50-54", 52.5,
                ifelse(trend_$age_grp == "55-59", 57.5,
                ifelse(trend_$age_grp == "60-64", 62.5, 
                ifelse(trend_$age_grp == "15+", 70, NA)))))))))))
  trend <- trend_[!(trend_$age_grp %in% c("15-49", "15+", "65-104")), ]
  trend <- trend[order(trend$age), ]
  trend$median <- trend$median * 100
  trend$lci <- trend$lci * 100
  trend$uci <- trend$uci * 100
  trend_ever <- subset(trend, outcome == "ever ipv")
  trend_past <- subset(trend, outcome == "past ipv")
  
  png("global_age_trend.png", width, height, unit = "in", res = 600)
  par(mfrow = c(1, 1), family = "sans", yaxs = "i", yaxs = "i", oma = c(0, 0, 0, 0), mar = c(4, 4, 2, 0))
  
  plot(trend_past$median ~ trend_past$age, type = 'n', axes = FALSE, yaxs = "i", xaxs = "i",
      xlim = c(17.5, 65), ylim = ylim, xlab = "Age groups", ylab = "Prevalence in 2018 (%)")
  axis(1, at = trend_past$age[order(trend_past$age)], 
       labels = trend_past$age_grp[order(trend_past$age)], las = 2, cex.axis = 0.8)
  axis(2, at = seq(0, ylim[2], by = 10), labels = seq(0, ylim[2], by = 10), las = 1)
  #abline(h = seq(0, 100, by = 10), col = "lightgray", lwd = 0.5)
      # past
      polygon(x = c(trend_past$age, rev(trend_past$age)),
        y = c(trend_past$lci, rev(trend_past$uci)),
        col = scales::alpha("#E44253", 0.5), border = NA)
      lines(trend_past$median ~ trend_past$age, lwd = 1, col = scales::alpha("#E44253", 1))
      # ever
      polygon(x = c(trend_ever$age, rev(trend_ever$age)),
        y = c(trend_ever$lci, rev(trend_ever$uci)),
        col = scales::alpha("#90C4AA", 0.5), border = NA)
      lines(trend_ever$median ~ trend_ever$age, lwd = 1, col = scales::alpha("#90C4AA", 1))

  legend("topright", lwd = 3, col = c(scales::alpha("#90C4AA", 1), scales::alpha("#E44253", 1)),
         legend = c("Lifetime", "Past year"), bty = "n")
  dev.off()
  
  png("global_age_trend_points.png", width, height, unit = "in", res = 600)
  par(mfrow = c(1, 1), family = "sans", yaxs = "i", yaxs = "i", oma = c(0, 0, 0, 0), mar = c(4, 4, 2, 0))
  
  plot(trend_past$median ~ trend_past$age, type = 'n', axes = FALSE, yaxs = "i", xaxs = "i",
      xlim = c(15, 65), ylim = ylim, xlab = "Age groups", ylab = "Prevalence (%)")
  axis(1, at = c(15, trend_past$age[order(trend_past$age)]), 
       labels = c(NA, trend_past$age_grp[order(trend_past$age)]), las = 2, cex.axis = 0.8)
  axis(2, at = seq(0, ylim[2], by = 10), labels = seq(0, ylim[2], by = 10), las = 1)
  #abline(h = seq(0, 100, by = 10), col = "lightgray", lwd = 0.5)
      # past
      segments(x0 = trend_past$age + 0.3, y0 = trend_past$lci,
               x1 = trend_past$age + 0.3, y1 = trend_past$uci,
               col = scales::alpha("#E44253", 0.75), lwd = 2)
      points(trend_past$median ~ I(trend_past$age + 0.3), pch = 19, cex = 0.6,
             col = scales::alpha("#E44253", 1))
      # ever
      segments(x0 = trend_ever$age - 0.3, y0 = trend_ever$lci,
               x1 = trend_ever$age - 0.3, y1 = trend_ever$uci,
               col = scales::alpha("#90C4AA", 0.75), lwd = 2)
      points(trend_ever$median ~ I(trend_ever$age - 0.3), pch = 19, cex = 0.6,
            col = scales::alpha("#90C4AA", 1))
  
  legend("topright", lwd = 3, col = c(scales::alpha("#90C4AA", 1), scales::alpha("#E44253", 1)),
         legend = c("Lifetime", "Past year"), bty = "n")
  dev.off()
}

pool_glb_ipv_nspv <- function(pred_ever, pred_npsv, cor, ll_age = 15, ul_age = 49, 
                              suffix, save_results = TRUE) {

  YPInd_ever <- pred_ever[["YPInd"]]
  y_denom_ever <- pred_ever[["y_denom"]]
  y_pred_ever <- pred_ever[["y_pred"]]
  YPInd_npsv <- pred_npsv[["YPInd"]]
  y_denom_npsv <- pred_npsv[["y_denom"]]
  y_pred_npsv <- pred_npsv[["y_pred"]]  
  
  age <- seq(17.5, 104, by = 5)

  # === ==== === === ===
  # prevalence 15-49
  y_pred_both_1549 <- array(NA, dim = c(dim(y_pred_ever)[1], dim(y_pred_ever)[3]))
  denom_both_1549 <- array(NA, dim = dim(y_pred_ever)[1])
  
  pb <- txtProgressBar(1, dim(y_pred_ever)[1], style = 3)
  for (i in 1:dim(y_pred_ever)[1]) {
    ipv_ever <- colSums(y_pred_ever[i, which(age >= ll_age & age <= ul_age), ]) / sum(y_denom_ever[i, which(age >= ll_age & age <= ul_age)])
    npsv_ever <- colSums(y_pred_npsv[i, which(age >= ll_age & age <= ul_age), ]) / sum(y_denom_npsv[i, which(age >= ll_age & age <= ul_age)])

    comb <- ipv_ever + npsv_ever - 
      (ipv_ever * npsv_ever + cor * sqrt(ipv_ever * (1 - ipv_ever) * npsv_ever * (1 - npsv_ever)))
    y_pred_both_1549[i, ] <- comb * sum(y_denom_ever[i, which(age >= ll_age & age <= ul_age)])
    denom_both_1549[i] <- sum(y_denom_ever[i, which(age >= ll_age & age <= ul_age)])
  setTxtProgressBar(pb, i)
  }
  close(pb)

  # Prevalence 1549
  tot_num_1549 <- colSums(y_pred_both_1549[, ])
  tot_den_1549 <- sum(denom_both_1549)
  prv_ <- tot_num_1549 / tot_den_1549
  prv_1549 <- quantile(prv_, probs = c(0.5, 0.025, 0.975))  
  
  # === ==== === === ===
  # 15+
  y_pred_both_15p <- array(NA, dim = c(dim(y_pred_ever)[1], dim(y_pred_ever)[3]))
  denom_both_15p <- array(NA, dim = dim(y_pred_ever)[1])
  
  pb <- txtProgressBar(1, dim(y_pred_ever)[1], style = 3)
  for (i in 1:dim(y_pred_ever)[1]) {
    ipv_ever <- colSums(y_pred_ever[i, , ]) / sum(y_denom_ever[i, ])
    npsv_ever <- colSums(y_pred_npsv[i, , ]) / sum(y_denom_npsv[i, ])

    comb <- ipv_ever + npsv_ever - 
      (ipv_ever * npsv_ever + cor * sqrt(ipv_ever * (1 - ipv_ever) * npsv_ever * (1 - npsv_ever)))
    y_pred_both_15p[i, ] <- comb * sum(y_denom_ever[i, ])
    denom_both_15p[i] <- sum(y_denom_ever[i, ])
  setTxtProgressBar(pb, i)
  }
  close(pb)

  # Prevalence 15+  
  tot_num_15p <- colSums(y_pred_both_15p[, ])
  tot_den_15p <- sum(denom_both_15p)
  prv_ <- tot_num_15p / tot_den_15p
  prv_1599 <- quantile(prv_, probs = c(0.5, 0.025, 0.975))
  
  global_res <- data.frame(age_grp = c(paste(ll_age, "-", ul_age, sep = ""), "15+"),
                    median = c(prv_1549[1], prv_1599[1]), 
                    lci = c(prv_1549[2], prv_1599[2]), 
                    uci = c(prv_1549[3], prv_1599[3]),
                    women = c(quantile(tot_num_1549, probs = 0.5), 
                              quantile(tot_num_15p, probs = 0.5)),
                    women_lci = c(quantile(tot_num_1549, probs = 0.025),
                                  quantile(tot_num_15p, probs = 0.025)),
                    women_uci = c(quantile(tot_num_1549, probs = 0.975),
                                  quantile(tot_num_15p, probs = 0.975)))
       
  # by 10-year age groups
  age_grp_lb <- seq(15, 65, by = 10) 
  age_grp_ub <- c(seq(24, 64, by = 10), 104)
  
  prv_age_stat <- NULL
  pb <- txtProgressBar(1, length(age_grp_lb), style = 3)
  for (i in 1:length(age_grp_lb)) {
    y_pred_both <- array(NA, dim = c(dim(y_pred_ever)[1], dim(y_pred_ever)[3]))
    denom_both <- array(NA, dim = dim(y_pred_ever)[1])
    for (j in 1:dim(y_pred_ever)[1]) {
      ipv_ever <- colSums(y_pred_ever[j, which(age >= age_grp_lb[i] & age <= age_grp_ub[i]), ]) / 
                      sum(y_denom_ever[j, which(age >= age_grp_lb[i] & age <= age_grp_ub[i])])
      npsv_ever <- colSums(y_pred_npsv[j, which(age >= age_grp_lb[i] & age <= age_grp_ub[i]), ]) / 
                      sum(y_denom_npsv[j, which(age >= age_grp_lb[i] & age <= age_grp_ub[i])])

      comb <- ipv_ever + npsv_ever - 
          (ipv_ever * npsv_ever + cor * sqrt(ipv_ever * (1 - ipv_ever) * npsv_ever * (1 - npsv_ever)))
      y_pred_both[j, ] <- comb * sum(y_denom_ever[j, which(age >= age_grp_lb[i] & age <= age_grp_ub[i])])
      denom_both[j] <- sum(y_denom_ever[j, which(age >= age_grp_lb[i] & age <= age_grp_ub[i])])
    }
    prv_age_stat_i <- data.frame(age_grp = paste(age_grp_lb[i], "-", age_grp_ub[i], sep = ""),
                                  median = quantile(colSums(y_pred_both[, ] / sum(denom_both)), probs = 0.5),
                                  lci = quantile(colSums(y_pred_both[, ] / sum(denom_both)), probs = 0.025),
                                  uci = quantile(colSums(y_pred_both[, ] / sum(denom_both)), probs = 0.975),
                                  women = quantile(colSums(y_pred_both[, ]), probs = 0.5),
                                  women_lci = quantile(colSums(y_pred_both[, ]), probs = 0.025),
                                  women_uci = quantile(colSums(y_pred_both[, ]), probs = 0.975))
    prv_age_stat <- rbind(prv_age_stat, prv_age_stat_i)
    setTxtProgressBar(pb, i)
  }
  close(pb)  
  
  global_res <- rbind(global_res, prv_age_stat)
  global_res$women <- round(global_res$women) * 1000
  global_res$women_lci <- round(global_res$women_lci) * 1000
  global_res$women_uci <- round(global_res$women_uci) * 1000 
  
  if (save_results == TRUE) {
    write.csv(global_res, file = paste("Estimates Global IPV & NPSV ", suffix, ".csv", sep = ""))      
  }
  
  # Among 15+: By WHO region (including high-income as a separate category)
  comb_who <- NULL
  R <- unique(as.character(YPInd_ever$WHO_HI))
  pb <- txtProgressBar(1, length(R), style = 3)
  for (i in R) {
    D1_i <- YPInd_ever$D1[YPInd_ever$WHO_HI == i]
    tot_num_i <- colSums(y_pred_both_15p[D1_i, ])
    tot_den_i <- sum(denom_both_15p[D1_i])
    prv_i <- tot_num_i / tot_den_i
    comb_who_i <- data.frame(WHO_HI = i, age_grp = "15+",
                        median = quantile(prv_i, probs = 0.5),
                        lci = quantile(prv_i, probs = 0.025),
                        uci = quantile(prv_i, probs = 0.975),
                        women = round(quantile(tot_num_i, probs = 0.5)) * 1000,
                        women_lci = round(quantile(tot_num_i, probs = 0.025)) * 1000,
                        women_uci = round(quantile(tot_num_i, probs = 0.975)) * 1000)
  comb_who <- rbind(comb_who, comb_who_i)
  setTxtProgressBar(pb, which(R == i))
  }
  close(pb)
  rownames(comb_who) <- NULL
  comb_who_15p <- comb_who

  if (save_results == TRUE) {
    write.csv(comb_who_15p, file = paste("Estimates Global IPV & NPSV by WHO Regions (15+)",
              suffix, ".csv", sep = "")) }
  
  # Among 1549: By WHO region (including high-income as a separate category)
  comb_who <- NULL
  R <- unique(as.character(YPInd_ever$WHO_HI))
  pb <- txtProgressBar(1, length(R), style = 3)
  for (i in R) {
    D1_i <- YPInd_ever$D1[YPInd_ever$WHO_HI == i]
    tot_num_i <- colSums(y_pred_both_1549[D1_i, ])
    tot_den_i <- sum(denom_both_1549[D1_i])
    prv_i <- tot_num_i / tot_den_i
    comb_who_i <- data.frame(WHO_HI = i, age_grp = paste(ll_age, "-", ul_age, sep = ""),
                        median = quantile(prv_i, probs = 0.5),
                        lci = quantile(prv_i, probs = 0.025),
                        uci = quantile(prv_i, probs = 0.975),
                        women = round(quantile(tot_num_i, probs = 0.5)) * 1000,
                        women_lci = round(quantile(tot_num_i, probs = 0.025)) * 1000,
                        women_uci = round(quantile(tot_num_i, probs = 0.975)) * 1000)
  comb_who <- rbind(comb_who, comb_who_i)
  setTxtProgressBar(pb, which(R == i))
  }
  close(pb)
  rownames(comb_who) <- NULL
  comb_who_1549 <- comb_who
  
  if (save_results == TRUE) {
    write.csv(comb_who_1549, file = paste("Estimates Global IPV & NPSV by WHO Regions (",
              paste(ll_age, "-", ul_age, sep = ""), ")",
              suffix, ".csv", sep = "")) } 
  
  return(global_res)
}