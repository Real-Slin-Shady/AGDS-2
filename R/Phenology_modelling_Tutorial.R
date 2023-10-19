# I will use the phenocamr package which
# interfaces with the phenocam network API
# to download time series of vegetation
# greenness and derived phenology metrics
library(phenocamr)
library(dplyr)
library(ggplot2)
# download greenness time series,
# calculate phenology (phenophases),
# amend with DAYMET data
phenocamr::download_phenocam(
  site = "harvard$",
  veg_type = "DB",
  roi_id = "1000",
  daymet = TRUE,
  phenophase = TRUE,
  trim = 2022,
  out_dir = tempdir()
)

harvard_phenocam_data <- readr::read_csv(
  file.path(tempdir(), "harvard_DB_1000_3day.csv"),
  comment = "#"
)

# reading in harvard phenology only retaining
# spring (rising) phenology for the GCC 90th
# percentile time series (the default)
harvard_phenology <- readr::read_csv(
  file.path(
    tempdir(),
    "harvard_DB_1000_3day_transition_dates.csv"
  ),
  comment = "#"
) |>
  dplyr::filter(
    direction == "rising",
    gcc_value == "gcc_90"
  )










# return mean daily temperature as well
# as formal dates (for plotting)
harvard_temp <- harvard_phenocam_data |>
  group_by(year) |>
  dplyr::mutate(
    tmean = (tmax..deg.c. + tmin..deg.c.)/2
  ) |>
  dplyr::mutate(
    date = as.Date(date),
    gdd = cumsum(ifelse(tmean >= 5, tmean - 5, 0))
  ) |>
  dplyr::select(
    date,
    year,
    tmean,
    gdd
  ) |>
  dplyr::ungroup()

# convert the harvard phenology data and only
# retain required data
harvard_phenology <- harvard_phenology |>
  mutate(
    doy = as.numeric(format(as.Date(transition_25),"%j")),
    year = as.numeric(format(as.Date(transition_25),"%Y"))
  ) |>
  select(
    year,
    doy,
    transition_25,
    threshold_25
  )







gdd_model <- function(temp, par) {
  # split out parameters from a simple
  # vector of parameter values
  temp_threshold <- par[1]
  gdd_crit <- par[2]

  # accumulate growing degree days for
  # temperature data
  gdd <- cumsum(ifelse(temp > temp_threshold, temp - temp_threshold, 0))

  # figure out when the number of growing
  # degree days exceeds the minimum value
  # required for leaf development, only
  # return the first value
  doy <- unlist(which(gdd >= gdd_crit)[1])

  return(doy)
}








# run model and compare to true values
# returns the RMSE
rmse_gdd <- function(par, data) {

  # split out data
  drivers <- data$drivers
  validation <- data$validation

  # calculate phenology predictions
  # and put in a data frame
  predictions <- drivers |>
    group_by(year) |>
    summarise(
      predictions = gdd_model(
        temp = tmean,
        par = par
      )
    )

  predictions <- left_join(predictions, validation, by = "year")

  rmse <- predictions |>
    summarise(
      rmse = sqrt(mean((predictions - doy)^2, na.rm = TRUE))
    ) |>
    pull(rmse)

  # return rmse value
  return(rmse)
}






# starting model parameters
par = c(0, 130)

# limits to the parameter space
lower <- c(-10,0)
upper <- c(45,500)

# data needs to be provided in a consistent
# single data file, a nested data structure
# will therefore accept non standard data formats
data <- list(
  drivers = harvard_temp,
  validation = harvard_phenology
)

# optimize the model parameters
optim_par = GenSA::GenSA(
  par = par,
  fn = rmse_gdd,
  lower = lower,
  upper = upper,
  control = list(
    max.call = 4000
  ),
  data = data
)$par



predictions <- harvard_temp |>
  group_by(year) |>
  summarize(
    prediction = gdd_model(
      temp = tmean,
      par = optim_par
    )
  )


# join predicted with observed data
validation <- left_join(predictions, harvard_phenology)

ggplot2::ggplot(validation) +
  geom_smooth(
    aes(
      doy,
      prediction
    ),
    colour = "grey25",
    method = "lm"
  ) +
  geom_point(
    aes(
      doy,
      prediction
    )
  ) +
  geom_abline(
    intercept=0,
    slope=1,
    linetype="dotted"
  ) +
  labs(
    x = "Observed leaf-out date (DOY)",
    y = "Predicted leaf-out date (DOY)"
  ) +
  theme_bw()  +
  theme(
    legend.position = "none"
  )





library(daymetr)

# Download daily data
daymetr::download_daymet_tiles(
  tiles = 11935,
  start = 2012,
  end = 2012,
  param = c("tmin","tmax"),
  path = paste0(here::here(), "/data-raw/"),
  silent = TRUE
)

# calculate the daily mean values
r <- daymetr::daymet_grid_tmean(
  path = paste0(here::here(), "/data-raw/"),
  product = 11935,
  year = 2012,
  internal = TRUE
)


# reproject to lat lon
r <- terra::project(
  r,
  "+init=epsg:4326"
)

# subset to first 180 days
ma_nh_temp <- terra::subset(
  r,
  1:180
)

predicted_phenology <- terra::app(
  ma_nh_temp,
  fun = gdd_model,
  par = optim_par
)


library(leaflet)
library(terra)

# set te colour scale manually
pal <- colorNumeric(
  "magma",
  values(predicted_phenology),
  na.color = "transparent"
)

# build the leaflet map
# using ESRI tile servers
# and the loaded demo raster
leaflet() |>
  addProviderTiles(providers$Esri.WorldImagery, group = "World Imagery") |>
  addProviderTiles(providers$Esri.WorldTopoMap, group = "World Topo") |>
  addRasterImage(
    predicted_phenology,
    colors = pal,
    opacity = 0.8,
    group = "Phenology model results"
  ) |>
  addLayersControl(
    baseGroups = c("World Imagery","World Topo"),
    position = "topleft",
    options = layersControlOptions(collapsed = FALSE),
    overlayGroups = c("Phenology model results")
  ) |>
  addLegend(
    pal = pal,
    values = terra::values(predicted_phenology),
    title = "DOY")
