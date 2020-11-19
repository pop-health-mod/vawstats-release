devtools::load_all(here::here())

who_general <- rgdal::readOGR("./data-raw/Shapefiles/general_2013.shp")
#who_general <- sp::spTransform(who_general, sp::CRS("+proj=wintri"))
who_id <- data.frame(id = who_general$ISO_3_CODE, iso3 = who_general$ISO_3_CODE, nom = who_general$CNTRY_TERR)
who_general <- broom::tidy(who_general, region = "ISO_3_CODE")
world_map <- who_general
world_map <- merge(who_general, who_id, by = "id")

mask_line <- rgdal::readOGR("./data-raw/Shapefiles/maskline_general_2013.shp")
#mask_line <- sp::spTransform(mask_line, sp::CRS("+proj=wintri"))
lines_id <- data.frame(id = as.numeric(mask_line$ID), cnt = mask_line$COUNTRY)
who_lines <- broom::tidy(mask_line, region = "ID")
who_lines <- merge(who_lines, lines_id, by = "id")

mask_poly <- rgdal::readOGR("./data-raw/Shapefiles/maskpoly_general_2013.shp")
#mask_poly <- sp::spTransform(mask_poly, sp::CRS("+proj=wintri"))
who_poly <- broom::tidy(mask_poly, region = "AREA")

    gg <- ggplot2::ggplot()
    gg <- gg + ggplot2::geom_map(data = world_map, map = world_map,
                        fill = "#ffffff", color = NA,
                        ggplot2::aes(x = world_map$long, y = world_map$lat, map_id = world_map$id, group = world_map$group))
    gg <- gg + ggplot2::geom_map(data = who_poly, map = who_poly, fill = "white", color = NA,
                        ggplot2::aes(x = long, y = lat, map_id = id))
    gg <- gg + ggplot2::geom_path(data = who_lines, color = "grey50",
                        ggplot2::aes(x = long, y = lat, group = as.numeric(group)))
    gg

    
usethis::use_data(world_map, overwrite = TRUE)
usethis::use_data(who_lines, overwrite = TRUE)
usethis::use_data(who_poly, overwrite = TRUE)

# URL <- "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_admin_0_countries.geojson"
# fil <- basename(URL)
# if (!file.exists(fil)) download.file(URL, fil)

# world <- rgdal::readOGR(dsn = "ne_50m_admin_0_countries.geojson", layer = "OGRGeoJSON")
# 
# # remove antarctica
# world <- world[!world$ISO_A3 %in% c("ATA"),]
# world$ISO_A3 <-  factor(world$ISO_A3, levels = c(levels(world$ISO_A3), "FRA", "NOR"))
# world$ISO_A3[world$ADMIN == "France"] <- "FRA"
# world$ISO_A3[world$ADMIN == "Norway"] <- "NOR"
# world <- sp::spTransform(world, sp::CRS("+proj=wintri"))
# world$ISO_A3 <- as.character(world$ISO_A3)
# world_map_old <- world
# 
# usethis::use_data(world_map_old, overwrite = TRUE)



