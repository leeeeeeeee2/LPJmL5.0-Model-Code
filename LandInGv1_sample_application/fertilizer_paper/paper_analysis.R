################################################################################
## Copyright (C) 2022 Potsdam Institute for Climate Impact Research (PIK),    ##
## see COPYRIGHT file.                                                        ##
##                                                                            ##
## This file is part of LandInG and licensed under GNU AGPL Version 3 or      ##
## later. See LICENSE file or go to http://www.gnu.org/licenses/              ##
## Contact: https://github.com/PIK-LPJmL/LandInG/                             ##
################################################################################

################################################################################
## Script to do analysis of generated fertilizer and manure inputs for        ##
## documentation paper                                                        ##
################################################################################

rm(list = ls())

library(udunits2)
library(RColorBrewer)
library(abind)
library(fields)
library(spatstat)

################################################################################
## Unit used in generated file. By default, LPJmL uses "g/m2"                 ##
cft_output_unit <- "g/m2"
## Grid file name                                                             ##
cft_gridname <- c(
  "5arcmin" = "../gadm_paper/grid_gadm_5arcmin.bin",
  "30arcmin" = "../gadm_paper/grid_gadm_30arcmin_predefined_grid.bin"
)
## Fertilizer input file name                                                 ##
cft_fertname <- c(
  "5arcmin" = "fert_N_default_cft_aggregation_5min_1860-2017.clm",
  "30arcmin" = "fert_N_default_cft_aggregation_30min_1860-2017.clm"
)
## Optional file accounting for non-allocated fertilizer (e.g. above maximum  ##
## application threshold). Set to NULL if not applicable.                     ##
cft_fert_noalloc_name <- c(
  "5arcmin" = "fert_noalloc_N_default_cft_aggregation_5min_1860-2017.clm",
  "30arcmin" = "fert_noalloc_N_default_cft_aggregation_30min_1860-2017.clm"
)
cft_fert_noalloc_unit <- "g/m2"
## Manure input file name                                                     ##
cft_manuname <- c(
  "5arcmin" = "manure_N_default_cft_aggregation_5min_1860-2017.clm",
  "30arcmin" = "manure_N_default_cft_aggregation_30min_1860-2017.clm"
)
## Optional file accounting for non-allocated manure (e.g. missing cropland   ##
## or above maximum application threshold). Set to NULL if not applicable.    ##
cft_manu_noalloc_name <- c(
  "5arcmin" = "manure_noalloc_N_default_cft_aggregation_5min_1860-2017.clm",
  "30arcmin" = "manure_noalloc_N_default_cft_aggregation_30min_1860-2017.clm"
)
cft_manu_noalloc_global_name <- c(
  "5arcmin" = "manure_noalloc_N_default_cft_aggregation_5min_1860-2017.csv",
  "30arcmin" = "manure_noalloc_N_default_cft_aggregation_30min_1860-2017.csv"
)
cft_manu_noalloc_unit <- "g"
## Corresponding land use input                                               ##
cft_landuse_name <- c(
  "5arcmin" = file.path(
    "..", "landuse_paper",
    "cft_default_cft_aggregation_20200417_20200127_5min_1500-2017.bin"
  ),
  "30arcmin" = file.path(
    "..", "landuse_paper",
    "cft_default_cft_aggregation_20200417_20200127_30min_1500-2017.bin"
  )
)
cft_landuse_unit <- ""
## Corresponding fallow land input. Set to NULL if fallow land is included in ##
## cft_bands and therefore part of cft_landuse_name                           ##
cft_fallow_name <- c(
  "5arcmin" = file.path(
    "..", "landuse_paper",
    "cft_fallow_default_cft_aggregation_20200417_20200127_5min_1500-2017.bin"
  ),
  "30arcmin" = file.path(
    "..", "landuse_paper",
    "cft_fallow_default_cft_aggregation_20200417_20200127_30min_1500-2017.bin"
  )
)
cft_fallow_unit <- ""
## Crop band names
cft_bands <- c(
  "rainfed temperate cereals",
  "rainfed rice",
  "rainfed maize",
  "rainfed tropical cereals",
  "rainfed pulses",
  "rainfed temperate roots",
  "rainfed tropical roots",
  "rainfed oil crops sunflower",
  "rainfed oil crops soybean",
  "rainfed oil crops groundnut",
  "rainfed oil crops rapeseed",
  "rainfed sugar cane",
  "rainfed other crops",
  "rainfed pasture/managed grass",
  "rainfed bio-energy grass",
  "rainfed bio-energy tree",
  "irrigated temperate cereals",
  "irrigated rice",
  "irrigated maize",
  "irrigated tropical cereals",
  "irrigated pulses",
  "irrigated temperate roots",
  "irrigated tropical roots",
  "irrigated oil crops sunflower",
  "irrigated oil crops soybean",
  "irrigated oil crops groundnut",
  "irrigated oil crops rapeseed",
  "irrigated sugar cane",
  "irrigated other crops",
  "irrigated pasture/managed grass",
  "irrigated bio-energy grass",
  "irrigated bio-energy tree"
)
# Do not use separate fallow input if fallow is included in cft_bands
if (length(grep("fallow", cft_bands)) == 2)
  cft_fallow_name <- NULL
################################################################################

if (file.exists("../lpjml_format_helper_functions.R")) {
  source("../lpjml_format_helper_functions.R")
} else {
  stop("Please update path to script with LPJmL input format helper function")
}

################################################################################
## Paper analysis                                                             ##
cft_gridheader <- cft_griddata <- cft_gridarea <- list()
for (resol in names(cft_gridname)) {
  cft_gridheader[[resol]] <- tmpheader <- read_header(cft_gridname[resol])
  zz <- file(cft_gridname[resol], "rb")
  seek(zz, get_headersize(tmpheader))
  cft_griddata[[resol]] <- matrix(
    readBin(
      zz,
      what = get_datatype(tmpheader)$type,
      size = get_datatype(tmpheader)$size,
      n = prod(tmpheader$header[c("ncell", "nbands")]),
      endian = tmpheader$endian
    ) * tmpheader$header["scalar"],
    ncol = tmpheader$header["nbands"],
    byrow = tmpheader$header["order"] != 4
  )
  close(zz)
  cft_gridarea[[resol]] <- cellarea(
    cft_griddata[[resol]][, 2],
    tmpheader$header["cellsize_lon"],
    tmpheader$header["cellsize_lat"]
  )
  # Add unit conversion here if cft_output_unit is not per m2
}
## Calculate global application sums                                          ##
first_pass_RData <- "paper_analysis_1st_pass.RData"
reload <- TRUE
if (file.exists(first_pass_RData)) {
  cat("Checking if data loaded by previous script run can be reused from",
      sQuote(first_pass_RData), "\n")
  reload <- FALSE
  tmp_env <- new.env()
  load(first_pass_RData, envir = tmp_env)
  check_vars <-  c("cft_gridname", "cft_fertname", "cft_fert_noalloc_name",
                   "cft_manuname", "cft_manu_noalloc_name",
                   "cft_manu_noalloc_global_name", "cft_landuse_name",
                   "cft_fallow_name")
  for (check in check_vars) {
    reload <- reload || !identical(get(check), tmp_env[[check]])
  }
  rm (tmp_env)
  if (reload) {
    cat("Data cannot be reused because input files have changed\n")
  } else {
    load(first_pass_RData)
  }
}
if (reload) {
  cft_fert_header <- cft_fert_noalloc_header <- cft_manu_header <-
    cft_manu_noalloc_header <- cft_landuse_header <- cft_fallow_header <- list()
  cft_fert_sum <- cft_fert_noalloc_sum <- cft_manu_sum <-
    cft_manu_noalloc_sum <- cft_landuse_sum <- cft_manu_fallow_sum <- list()
  cft_fert_max <- cft_manu_max <- list()
  for (resol in names(cft_gridname)) {
    cat("** Resolution", sQuote(resol), "**\n")
    # Fertilizer data
    cft_fert_header[[resol]] <- tmpheader <- read_header(cft_fertname[resol])
    # Fertilizer not allocated
    if (resol %in% names(cft_fert_noalloc_name)) {
      cft_fert_noalloc_header[[resol]] <- tmpheader <-
        read_header(cft_fert_noalloc_name[resol])
      if (any(tmpheader$header != cft_fert_header[[resol]]$header)) {
        stop(
          "Mismatch between ", sQuote(cft_fertname[resol]), " and ",
          sQuote(cft_fert_noalloc_name[resol])
        )
      }
    }
    # Land use data
    cft_landuse_header[[resol]] <- tmpheader <-
      read_header(cft_landuse_name[resol])
    header_par <- c("ncell", "nbands", "cellsize_lon", "cellsize_lat",
                    "firstcell")
    if (
      any(tmpheader$header[header_par] !=
        cft_fert_header[[resol]]$header[header_par]
      )
    ) {
      stop(
        "Mismatch between ", sQuote(cft_fertname[resol]), " and ",
        sQuote(cft_landuse_name[resol])
      )
    }
    # Manure data
    cft_manu_header[[resol]] <- tmpheader <- read_header(cft_manuname[resol])
    if (
      any(tmpheader$header[header_par] !=
        cft_fert_header[[resol]]$header[header_par])
    ) {
      stop(
        "Mismatch between ", sQuote(cft_manuname[resol]), " and ",
        sQuote(cft_landuse_name[resol])
      )
    }
    # Manure not allocated
    if (resol %in% names(cft_manu_noalloc_name)) {
      cft_manu_noalloc_header[[resol]] <- tmpheader <-
        read_header(cft_manu_noalloc_name[resol])
      header_par <- c("ncell", "cellsize_lon", "cellsize_lat", "firstcell")
      if (
        any(tmpheader$header[header_par] !=
          cft_manu_header[[resol]]$header[header_par])
      ) {
        stop(
          "Mismatch between ", sQuote(cft_manu_noalloc_name[resol]), " and ",
          sQuote(cft_manu_name[resol])
        )
      }
    }
    # Fallow data
    if (resol %in% names(cft_fallow_name)) {
      cft_fallow_header[[resol]] <- tmpheader <-
        read_header(cft_fallow_name[resol])
      header_par <- c("ncell", "cellsize_lon", "cellsize_lat", "firstcell")
      if (
        any(tmpheader$header[header_par] !=
          cft_manu_header[[resol]]$header[header_par])
      ) {
        stop(
          "Mismatch between ", sQuote(cft_fallow_name[resol]), " and ",
          sQuote(cft_manu_name[resol])
        )
      }
    }
    firstyear <- max(
      cft_fert_header[[resol]]$header["firstyear"],
      ifelse(
        is.null(cft_fert_noalloc_header[[resol]]),
        -Inf,
        cft_fert_noalloc_header[[resol]]$header["firstyear"]
      ),
      cft_landuse_header[[resol]]$header["firstyear"],
      cft_manu_header[[resol]]$header["firstyear"],
      ifelse(
        is.null(cft_manu_noalloc_header[[resol]]),
        -Inf,
        cft_manu_noalloc_header[[resol]]$header["firstyear"]
      ),
      ifelse(
        is.null(cft_fallow_header[[resol]]),
        -Inf,
        cft_fallow_header[[resol]]$header["firstyear"]
      )
    )
    lastyear <- min(
      sum(cft_fert_header[[resol]]$header[c("firstyear", "nyear")]),
      ifelse(
        is.null(cft_fert_noalloc_header[[resol]]),
        Inf,
        sum(cft_fert_noalloc_header[[resol]]$header[c("firstyear", "nyear")])
      ),
      sum(cft_landuse_header[[resol]]$header[c("firstyear", "nyear")]),
      sum(cft_manu_header[[resol]]$header[c("firstyear", "nyear")]),
      ifelse(
        is.null(cft_manu_noalloc_header[[resol]]),
        Inf,
        sum(cft_manu_noalloc_header[[resol]]$header[c("firstyear", "nyear")])
      ),
      ifelse(
        is.null(cft_fallow_header[[resol]]),
        Inf,
        sum(cft_fallow_header[[resol]]$header[c("firstyear", "nyear")])
      )
    ) - 1
    nyear <- lastyear - firstyear + 1
    cft_fert_sum[[resol]] <- array(
      dim = c(cft_fert_header[[resol]]$header["nbands"], nyear)
    )
    if (cft_fert_header[[resol]]$header["nbands"] == length(cft_bands)) {
      dimnames(cft_fert_sum[[resol]]) <- list(cft_bands,
                                              seq(firstyear, lastyear))
    } else {
      dimnames(cft_fert_sum[[resol]]) <- list(NULL, seq(firstyear, lastyear))
    }
    cft_fert_max[[resol]] <- cft_manu_sum[[resol]] <- cft_manu_max[[resol]] <-
      cft_landuse_sum[[resol]] <- cft_fert_sum[[resol]]
    if (!is.null(cft_fert_noalloc_header[[resol]]))
      cft_fert_noalloc_sum[[resol]] <- cft_fert_sum[[resol]]
    if (!is.null(cft_manu_noalloc_header[[resol]])) {
      cft_manu_noalloc_sum[[resol]] <- array(
        dim = c(cft_manu_noalloc_header[[resol]]$header["nbands"], nyear)
      )
      if (cft_manu_noalloc_header[[resol]]$header["nbands"] ==
        length(cft_bands)
      ) {
        dimnames(cft_manu_noalloc_sum[[resol]]) <-
          list(cft_bands, seq(firstyear, lastyear))
      } else {
        dimnames(cft_manu_noalloc_sum[[resol]]) <-
          list(NULL, seq(firstyear, lastyear))
      }
    }
    if (!is.null(cft_fallow_header[[resol]])) {
      cft_manu_fallow_sum[[resol]] <- array(
        dim = c(cft_fallow_header[[resol]]$header["nbands"], nyear)
      )
      if (cft_fallow_header[[resol]]$header["nbands"] == length(cft_bands)) {
        dimnames(cft_manu_fallow_sum[[resol]]) <-
          list(cft_bands, seq(firstyear, lastyear))
      } else {
        dimnames(cft_manu_fallow_sum[[resol]]) <-
          list(NULL, seq(firstyear, lastyear))
      }
    }
    zz <- list()
    for (v in c("landuse", "fert", "manu", "fert_noalloc", "manu_noalloc",
                "fallow")) {
      tmpheader <- switch(
        v,
        landuse = cft_landuse_header[[resol]],
        fert = cft_fert_header[[resol]],
        manu = cft_manu_header[[resol]],
        fert_noalloc = cft_fert_noalloc_header[[resol]],
        manu_noalloc = cft_manu_noalloc_header[[resol]],
        fallow = cft_fallow_header[[resol]]
      )
      if (is.null(tmpheader))
        next
      zz[[v]] <- file(
        description = switch(
          v,
          landuse = cft_landuse_name[resol],
          fert = cft_fertname[resol],
          manu = cft_manuname[resol],
          fert_noalloc = cft_fert_noalloc_name[resol],
          manu_noalloc = cft_manu_noalloc_name[resol],
          fallow = cft_fallow_name[resol]
        ),
        "rb"
      )
      pos <- prod(tmpheader$header[c("ncell", "nbands")]) *
        (firstyear - tmpheader$header["firstyear"]) *
        get_datatype(tmpheader)$size + get_headersize(tmpheader)
      seek(zz[[v]], pos)
    }
    for (y in seq_len(nyear)) {
      cat(firstyear + y - 1, "\n")
      for (v in c("landuse", "fert", "manu", "fert_noalloc", "manu_noalloc",
                  "fallow")) {
        # Make sure "landuse" is loaded before other variables that require it.
        # Also load fallow land after manure rate.
        tmpheader <- switch(
          v,
          landuse = cft_landuse_header[[resol]],
          fert = cft_fert_header[[resol]],
          manu = cft_manu_header[[resol]],
          fert_noalloc = cft_fert_noalloc_header[[resol]],
          manu_noalloc = cft_manu_noalloc_header[[resol]],
          fallow = cft_fallow_header[[resol]]
        )
        if (is.null(tmpheader))
          next
        var_sum <- switch(
          v,
          fert = cft_fert_sum,
          manu = cft_manu_sum,
          fert_noalloc = cft_fert_noalloc_sum,
          manu_noalloc = cft_manu_noalloc_sum,
          landuse = cft_landuse_sum,
          fallow = cft_manu_fallow_sum
        )
        var_max <- switch(
          v,
          fert = cft_fert_max,
          manu = cft_manu_max,
          NULL
        )
        var_unit <- switch(
          v,
          fert = cft_output_unit,
          manu = cft_output_unit,
          fert_noalloc = cft_fert_noalloc_unit,
          manu_noalloc = cft_manu_noalloc_unit,
          landuse = cft_landuse_unit,
          fallow = cft_fallow_unit
        )
        tmpdata <- matrix(
          readBin(
            zz[[v]],
            what = get_datatype(tmpheader)$type,
            size = get_datatype(tmpheader)$size,
            n = prod(tmpheader$header[c("ncell", "nbands")]),
            endian = tmpheader$endian
          ) * tmpheader$header["scalar"],
          ncol = tmpheader$header["nbands"],
          byrow = tmpheader$header["nbands"]
        )
        crop_max <- apply(tmpdata, 2, max)
        if (!is.null(var_max)) {
          var_max[[resol]][, y] <- crop_max
          assign(
            switch(
              v,
              fert = "cft_fert_max",
              manu = "cft_manu_max",
              "trash"
            ),
            var_max
          )
        }
        if (ud.are.convertible(var_unit, "g/m2")) {
          if (v == "manu") {
            # Save manure rate for use with fallow land. Using max over all
            # bands assumes that all crops receive the same manure rate.
            if (ncol(tmpdata) == length(cft_bands)) {
              manurate_data <- cbind(
                apply(tmpdata[, grep("rainfed", cft_bands)], 1, max),
                apply(tmpdata[, grep("irrigated", cft_bands)], 1, max)
              )
            } else {
              manurate_data <- apply(tmpdata, 1, max)
            }
          }
          # Multiply by area to derive absolute amount
          tmpdata <- tmpdata * area_data
        } else if (v == "manu" && ud.are.convertible(var_unit, "kg")) {
          manurate_data <- tmpdata / area_data
          manurate_data[which(area_data == 0)] <- 0
          if (ncol(manurate_data) == length(cft_bands)) {
            manurate_data <- cbind(
              apply(manurate_data[, grep("rainfed", cft_bands)], 1, max),
              apply(manurate_data[, grep("irrigated", cft_bands)], 1, max)
            )
          } else {
            manurate_data <- apply(manurate_data, 1, max)
          }
          gc(full = FALSE)
        }
        if (v == "landuse" && ud.are.convertible(var_unit, "1")) {
          area_data <- tmpdata <- ud.convert(tmpdata, var_unit, "1") *
            cft_gridarea[[resol]]
          if (any(rowSums(area_data) > (cft_gridarea[[resol]] * 1.0001))) {
            warning(
              "Land use areas exceed cell area in ",
              length(
                which(rowSums(area_data) > (cft_gridarea[[resol]] * 1.0001))
              ),
              " cells by up to ",
              max(rowSums(area_data) - cft_gridarea[[resol]]), " m2",
              call. = FALSE, immediate. = TRUE
            )
          }
        } else if (v == "landuse") {
          area_data <- tmpdata
          if (any(rowSums(area_data) > (cft_gridarea[[resol]] * 1.0001))) {
            warning(
              "Land use areas exceed cell area in ",
              length(
                which(rowSums(area_data) > (cft_gridarea[[resol]] * 1.0001))
              ),
              " cells by up to ",
              max(rowSums(area_data) - cft_gridarea[[resol]]), " m2",
              call. = FALSE, immediate. = TRUE
            )
          }
        }
        if (v == "fallow" && ud.are.convertible(var_unit, "1")) {
          # Convert to fallow land areas
          tmpdata <- ud.convert(tmpdata, var_unit, "1") * cft_gridarea[[resol]]
          if (any(rowSums(tmpdata) > (cft_gridarea[[resol]] * 1.0001))) {
            warning(
              "Fallow areas exceed cell area in ",
              length(which(rowSums(tmpdata) > (cft_gridarea[[resol]] * 1.0001))),
              " cells by up to ",
              max(rowSums(tmpdata) - cft_gridarea[[resol]]), " m2",
              call. = FALSE, immediate. = TRUE
            )
          }
        }
        if (v == "fallow") {
          # Assume same manure application on fallow as on crops.
          if (length(manurate_data) > length(tmpdata))
            stop("Mismatch betwwen manurate_data and tmpdata")
          tmpdata <- tmpdata * manurate_data
        }
        var_sum[[resol]][, y] <- colSums(tmpdata)
        assign(
          switch(
            v,
            landuse = "cft_landuse_sum",
            fert = "cft_fert_sum",
            manu = "cft_manu_sum",
            fert_noalloc = "cft_fert_noalloc_sum",
            manu_noalloc = "cft_manu_noalloc_sum",
            fallow = "cft_manu_fallow_sum"
          ),
          var_sum
        )
      }
      rm(area_data, manurate_data)
    }
    sapply(zz, close)
    rm(zz)
  }
  ## Global application sums                                                  ##
  cft_manu_noalloc_global_data <- list()
  for (resol in names(cft_manu_noalloc_global_name)) {
    if (file.exists(cft_manu_noalloc_global_name[resol])) {
      cft_manu_noalloc_global_data[[resol]] <-
        read.csv(cft_manu_noalloc_global_name[resol])
      # Check that years from CSV cover all years from gridded inputs
      if (!is.null(cft_manu_noalloc_sum[[resol]]) &&
        !all(dimnames(cft_manu_noalloc_sum[[resol]])[[2]] %in%
          cft_manu_noalloc_global_data[[resol]][, "Year"])
      ) {
        warning(
          "Global sums read from ", sQuote(cft_manu_noalloc_global_name[resol]),
          " do not match gridded data read from ",
          sQuote(cft_manu_noalloc_name[resol]),
          ". Missing year(s): ",
          toString(
            setdiff(
              dimnames(cft_manu_noalloc_sum[[resol]])[[2]],
              cft_manu_noalloc_global_data[[resol]][, "Year"]
            )
          ),
          call. = FALSE, immediate. = TRUE
        )
        # Do not use incomplete data
        cft_manu_noalloc_global_data[[resol]] <- NULL
      }
    } else {
      warning(
        "Global file ", sQuote(cft_manu_noalloc_global_name[resol]),
        " does not exist. Check cft_manu_noalloc_global_name[",
        dQuote(resol), "]",
        call. = FALSE, immediate. = TRUE
      )
    }
  }
  savevars <- c(
    "cft_gridname", "cft_fertname", "cft_fert_noalloc_name",
    "cft_manuname", "cft_manu_noalloc_name", "cft_manu_noalloc_global_name",
    "cft_landuse_name", "cft_fallow_name",
    "cft_fert_header", "cft_fert_noalloc_header",
    "cft_manu_header", "cft_manu_noalloc_header",
    "cft_landuse_header", "cft_fallow_header",
    "cft_fert_sum", "cft_fert_noalloc_sum",
    "cft_manu_sum", "cft_manu_noalloc_sum",
    "cft_landuse_sum", "cft_manu_fallow_sum",
    "cft_fert_max", "cft_manu_max", "cft_manu_noalloc_global_data",
    "firstyear", "lastyear"
  )
  save(list = savevars, file = first_pass_RData)
}

app_col <- brewer.pal(9, "Greens")[9]
noalloc_col <- rev(brewer.pal(9, "YlOrBr"))[6:8]
plot_seq <- c("fert", "manu")
cat("Global N application sums saved to paper_global_N_application_sums.pdf\n")
pdf("paper_global_N_application_sums.pdf",
    width = 6.5, height = 2 * length(plot_seq), paper = "special",
    pointsize = 8
)
par(
  mfrow = c(2, length(cft_gridname)),
  mar = c(0, 0, 1.5, 0) + 0.1,
  oma = c(3, 3, 0, 0) + 0.1
)
p <- 1
for (v in plot_seq) {
  app_sums <- sapply(get(paste0("cft_", v, "_sum")), colSums)
  noalloc_sums <- array(
    dim = c(dim(app_sums), 3),
    dimnames = c(
      dimnames(app_sums),
      list(c("preprocessing", "threshold", "fallow"))
    )
  )
  for (resol in names(cft_gridname)) {
    if (v == "manu" && !is.null(cft_manu_noalloc_global_data[[resol]])) {
      r <- match(
        rownames(app_sums),
        cft_manu_noalloc_global_data[[resol]][, "Year"]
      )
      cols <- c(
        grep("preprocessing", colnames(cft_manu_noalloc_global_data[[resol]])),
        grep("threshold", colnames(cft_manu_noalloc_global_data[[resol]]))
      )
      noalloc_sums[, resol, c("preprocessing", "threshold")] <-
        as.matrix(cft_manu_noalloc_global_data[[resol]][r, cols])
    } else if (!is.null(get(paste0("cft_", v, "_noalloc_sum"))[[resol]])) {
      noalloc_sums[, resol, "threshold"] <- 
        colSums(get(paste0("cft_", v, "_noalloc_sum"))[[resol]])
    }
    if (v == "manu" && !is.null(cft_manu_fallow_sum[[resol]])) {
      noalloc_sums[, resol, "fallow"] <- colSums(cft_manu_fallow_sum[[resol]])
    }
  }
  ylim <- c(
    0,
    max(app_sums + apply(noalloc_sums, c(1, 2), sum, na.rm = TRUE), na.rm = TRUE)
  )
  if (v == plot_seq[1]) {
    global_application_table <- abind(app_sums, noalloc_sums, along = 3)
  } else {
    global_application_table <- abind(
      global_application_table,
      abind(app_sums, noalloc_sums, along = 3),
      along = 4
    )
    dimnames(global_application_table)[[3]][1] <- "applied"
  }
  xlim <- range(as.integer(rownames(app_sums)))
  for (resol in names(cft_gridname)) {
    plot(1, 1, xlim = xlim, ylim = ylim, type = "n", xlab = "", ylab = "",
         axes = FALSE)
    bottom <- rep(0, nrow(app_sums))
    top <- app_sums[, resol]
    x <- c(as.integer(rownames(app_sums)), rev(as.integer(rownames(app_sums))))
    y <- c(bottom, rev(top))
    polygon(x, y, density = NA, col = app_col, border = NA)
    for (t in seq_len(dim(noalloc_sums)[3])) {
      if (anyNA(noalloc_sums[, resol, t]))
        next
      bottom <- top
      top <- bottom + noalloc_sums[, resol, t]
      y <- c(bottom, rev(top))
      polygon(x, y, density = NA, col = noalloc_col[t], border = NA)
    }
    if (resol == names(cft_gridname)[1]) {
      axis(2, at = pretty(ylim), labels = pretty(ylim) * 1e-12,
           mgp = c(2, 0.5, 0))
      title(
        ylab = paste0(
          "Global ",
          switch(v, manu = "manure", fert = "fertilizer"),
          " application [Tg]"
        ),
        mgp = c(2, 0.5, 0), xpd = NA
      )
      if (v == "manu") {
        leg_col <- app_col
        leg_text <- c(applied = "Applied nitrogen")
        noalloc_text <- c(
          preprocessing = "Removed by pre-processing (manure)",
          threshold = "Removed by application threshold",
          fallow = "Cropland as fallow land"
        )
        noalloc_present <- apply(
          noalloc_sums,
          3,
          function(indata) !all(is.na(indata))
        )
        if (!noalloc_present["preprocessing"]) {
          noalloc_text["threshold"] <-
            "Removed by pre-processing or application threshold"
        }
        if (any(noalloc_present)) {
          leg_col <- c(leg_col, noalloc_col[noalloc_present])
          leg_text <- c(leg_text, noalloc_text[noalloc_present])
        }
        legend("topleft", legend = rev(leg_text), fill = rev(leg_col), bty = "n")
      }
    }
    if (v == plot_seq[length(plot_seq)]) {
      axis(1, mgp = c(2, 0.5, 0))
      title(xlab = "Years", mgp = c(2, 0.5, 0), xpd = NA)
    }
    if (v == plot_seq[1]) {
      title(
        main = paste(sub("arcmin", "'", sub("arcsec", "''", resol)), "dataset"),
        line = 0.25
      )
    }
    box()
  }
}
dimnames(global_application_table)[[4]] <- plot_seq
dev.off()
global_application_table <- aperm(global_application_table, c(1, 2, 4, 3))
global_application_rel_table <- global_application_table / c(apply(
  global_application_table, c(1, 2, 3), sum, na.rm = TRUE
))

## Global application sum
cat("Total N application in Tg:\n")
first_fert_year <- apply(
  sapply(cft_fert_sum, colSums),
  2,
  function(indata) min(which(indata > 0))
) + firstyear - 1

stat_years <- unique(c(firstyear, first_fert_year, 1949, 1985, 1995, lastyear))
print(global_application_table[as.character(stat_years), , , "applied"] * 1e-12)

## Difference fertilizer application between datasets in %
cat("Difference between fertilizer datasets (30arcmin compared to...) in %:\n")
print(
  apply(
    global_application_table[, "30arcmin", "fert", "applied"] /
      global_application_table[, , "fert", "applied"] - 1,
    2,
    range,
    na.rm = TRUE
  ) * 100
)

## Maximum loss to fertilizer threshold in Tg and %
cat("Maximum fertilizer lost to application threshold:\n")
max_loss_year_abs <- apply(
  global_application_table[, , "fert", "threshold"],
  2,
  which.max
) + firstyear - 1
max_loss_year_rel <- apply(
  global_application_rel_table[, , "fert", "threshold"],
  2,
  which.max
) + firstyear - 1
max_loss_year <- union(max_loss_year_abs, max_loss_year_rel)
print(
  cbind(
    Tg = global_application_table[as.character(max_loss_year), , "fert", "threshold"] * 1e-12,
    "%" = global_application_rel_table[as.character(max_loss_year), , "fert", "threshold"] * 100
  )
)

## Difference manure application between datasets
cat("Difference between manure datasets (30arcmin compared to...) in %:\n")
print(
  apply(
    global_application_table[, "30arcmin", "manu", "applied"] /
      global_application_table[, , "manu", "applied"] - 1,
    2,
    range,
    na.rm = TRUE
  ) * 100
)

## Range of manure loss to preprocessing in %
cat("Range of manure loss to preprocessing in %\n")
print(
  apply(global_application_rel_table[, , "manu", "preprocessing"], 2, range) * 100
)
## Range of manure loss due to application threshold in %
cat("Range of manure loss due to application threshold in %\n")
print(
  apply(global_application_rel_table[, , "manu", "threshold"], 2, range) * 100
)
## Range of manure loss on fallow land in %
cat("Range of manure loss on fallow land in %\n")
print(
  apply(global_application_rel_table[, , "manu", "fallow"], 2, range) * 100
)


## Bin fertilizer/manure application rates and calculate growing area in each ##
## bin.                                                                       ##
second_pass_RData <- "paper_analysis_2nd_pass.RData"
reload <- TRUE
first_bin_year <- max(firstyear, 1950)
last_bin_year <- min(lastyear, 2017)
nyear_bin <- last_bin_year - first_bin_year + 1
bin_n <- 100
bin_seq <- seq(0, 1, length.out = bin_n + 1)
if (file.exists(second_pass_RData)) {
  cat("Checking if data loaded by previous script run can be reused from",
      sQuote(second_pass_RData), "\n")
  reload <- FALSE
  tmp_env <- new.env()
  load(second_pass_RData, envir = tmp_env)
  check_vars <-  c("cft_gridname", "cft_fertname", "cft_fert_noalloc_name",
                   "cft_manuname", "cft_manu_noalloc_name",
                   "cft_manu_noalloc_global_name", "cft_landuse_name",
                   "cft_fallow_name", "first_bin_year", "last_bin_year", "bin_n")
  for (check in check_vars) {
    reload <- reload || !identical(get(check), tmp_env[[check]])
  }
  rm (tmp_env)
  if (reload) {
    cat("Data cannot be reused because input files have changed\n")
  } else {
    load(second_pass_RData)
  }
}
if (reload) {
  cft_fert_bin <- cft_manu_bin <- cft_fert_bin_thr <- cft_manu_bin_thr <- list()
  for (resol in names(cft_gridname)) {
    cat("** Resolution", sQuote(resol), "**\n")
    zz <- list()
    for (v in c("landuse", "fert", "manu")) {
      tmpheader <- switch(
        v,
        landuse = cft_landuse_header[[resol]],
        fert = cft_fert_header[[resol]],
        manu = cft_manu_header[[resol]],
        fert_noalloc = cft_fert_noalloc_header[[resol]],
        manu_noalloc = cft_manu_noalloc_header[[resol]],
        fallow = cft_fallow_header[[resol]]
      )
      zz[[v]] <- file(
        description = switch(
          v,
          landuse = cft_landuse_name[resol],
          fert = cft_fertname[resol],
          manu = cft_manuname[resol],
        ),
        "rb"
      )
      pos <- prod(tmpheader$header[c("ncell", "nbands")]) *
        (first_bin_year - tmpheader$header["firstyear"]) *
        get_datatype(tmpheader)$size + get_headersize(tmpheader)
      seek(zz[[v]], pos)
    }
    cft_fert_bin[[resol]] <- array(
      dim = c(bin_n, cft_fert_header[[resol]]$header["nbands"], nyear_bin)
    )
    if (cft_fert_header[[resol]]$header["nbands"] == length(cft_bands)) {
      dimnames(cft_fert_bin[[resol]]) <- list(
        NULL, cft_bands, seq(first_bin_year, last_bin_year)
      )
    } else {
      dimnames(cft_fert_bin[[resol]]) <- list(
        NULL, NULL, seq(first_bin_year, last_bin_year)
      )
    }
    cft_manu_bin[[resol]] <- cft_fert_bin_thr[[resol]] <-
      cft_manu_bin_thr[[resol]] <- array(
      dim = dim(cft_fert_bin[[resol]]),
      dimnames = dimnames(cft_fert_bin[[resol]])
    )
    for (y in seq_len(nyear_bin)) {
      cat(first_bin_year + y - 1, "\n")
      for (v in c("landuse", "fert", "manu")) {
        tmpheader <- switch(
          v,
          landuse = cft_landuse_header[[resol]],
          fert = cft_fert_header[[resol]],
          manu = cft_manu_header[[resol]],
          fert_noalloc = cft_fert_noalloc_header[[resol]],
          manu_noalloc = cft_manu_noalloc_header[[resol]],
          fallow = cft_fallow_header[[resol]]
        )
        var_max <- switch(
          v,
          fert = cft_fert_max,
          manu = cft_manu_max,
          NULL
        )
        var_unit <- switch(
          v,
          fert = cft_output_unit,
          manu = cft_output_unit,
          fert_noalloc = cft_fert_noalloc_unit,
          manu_noalloc = cft_manu_noalloc_unit,
          landuse = cft_landuse_unit,
          fallow = cft_fallow_unit
        )
        var_bin <- switch(
          v,
          fert = cft_fert_bin,
          manu = cft_manu_bin,
          NULL
        )
        var_bin_thr <- switch(
          v,
          fert = cft_fert_bin_thr,
          manu = cft_manu_bin_thr,
          NULL
        )
        tmpdata <- matrix(
          readBin(
            zz[[v]],
            what = get_datatype(tmpheader)$type,
            size = get_datatype(tmpheader)$size,
            n = prod(tmpheader$header[c("ncell", "nbands")]),
            endian = tmpheader$endian
          ) * tmpheader$header["scalar"],
          ncol = tmpheader$header["nbands"],
          byrow = tmpheader$header["nbands"]
        )
        if (v != "landuse") {
          # Normalize to maximum value during during binning period. Use same
          # normalization for all resolutions.
          yseq <- match(
            seq(first_bin_year, last_bin_year),
            dimnames(var_max[[resol]])[[2]]
          )
          crop_max <- apply(var_max[[resol]][, yseq, drop = FALSE], 1, max)
          for (r2 in setdiff(names(var_max), resol)) {
            yseq <- match(
              seq(first_bin_year, last_bin_year),
              dimnames(var_max[[r2]])[[2]]
            )
            crop_max <- pmax(
              crop_max,
              apply(var_max[[r2]][, yseq, drop = FALSE], 1, max)
            )
          }
          var_bin_thr[[resol]][, , y] <- rep(crop_max, each = bin_n) * bin_seq[-1]
          for (b in seq_along(crop_max)) {
            if (crop_max[b] > 0) {
              tmpdata[, b] <- tmpdata[, b] / crop_max[b]
            }
            tmpdata[, b] <- findInterval(tmpdata[, b], bin_seq,
                                         all.inside = TRUE)
            for (bin in unique(tmpdata[, b])) {
              bincells <- which(tmpdata[, b] == bin)
              var_bin[[resol]][bin, b, y] <- sum(area_data[bincells, b])
              rm(bincells)
            }
          }
          assign(
            switch(
              v,
              fert = "cft_fert_bin",
              manu = "cft_manu_bin",
            ),
            var_bin
          )
          assign(
            switch(
              v,
              fert = "cft_fert_bin_thr",
              manu = "cft_manu_bin_thr",
            ),
            var_bin_thr
          )
        } else {
          if (ud.are.convertible(var_unit, "1")) {
            tmpdata <- ud.convert(tmpdata, var_unit, "1") * cft_gridarea[[resol]]
          }
          area_data <- tmpdata
        }
      }
      rm(area_data)
    }
    sapply(zz, close)
    rm(zz)
  }
  savevars <- c("cft_gridname", "cft_fertname", "cft_fert_noalloc_name",
                "cft_manuname", "cft_manu_noalloc_name",
                "cft_manu_noalloc_global_name", "cft_landuse_name",
                "cft_fallow_name",
                "first_bin_year", "last_bin_year", "bin_n",
                "cft_fert_bin", "cft_manu_bin",
                "cft_fert_bin_thr", "cft_manu_bin_thr")
  save(list = savevars, file = second_pass_RData)
}

plot_crops <- grep(
  "rainfed",
  grep("temperate cereals|maize|soybean", cft_bands, value = TRUE),
  value = TRUE
)
plot_crop_names <- sub("rainfed.+soybean.*", "rainfed soybean", plot_crops)
names(plot_crop_names) <- plot_crops
area_pal <- colorRampPalette(brewer.pal(9, "YlGn"))
cat("Areas per N rate saved to paper_nitrogen_rates_crop.pdf\n")
pdf(
  "paper_nitrogen_rates_crop.pdf",
  width = 6.5, height = 2 * length(plot_crops), paper = "special",
  pointsize = 8
)
par(
  mfrow = c(length(plot_crops), length(cft_fert_bin)),
  mar = c(0, 0, 1.25, 0.75) + 0.1,
  oma = c(3.5, 3.6, 0, 8) + 0.1
)
for (cft in plot_crops) {
  crop_max <- 0
  area_max <- 0
  area_min <- Inf
  for (resol in names(cft_fert_bin)) {
    yseq <- match(
      seq(first_bin_year, last_bin_year),
      dimnames(cft_fert_max[[resol]])[[2]]
    )
    crop_max <- max(crop_max, max(cft_fert_max[[resol]][cft, yseq]))
    area_max <- max(area_max, max(cft_fert_bin[[resol]][, cft, ], na.rm = TRUE))
    tmparea <- cft_fert_bin[[resol]][, cft, ]
    tmparea[which(tmparea == 0)] <- NA
    area_min <- min(area_min, min(tmparea, na.rm = TRUE))
  }
  for (resol in names(cft_fert_bin)) {
    image(
      x = seq(first_bin_year - 0.5, length.out = nyear_bin + 1),
      y = seq(0, crop_max * 10, length.out = bin_n + 1),
      z = t(cft_fert_bin[[resol]][, cft, ]),
      zlim = c(0, area_max),
      col = area_pal(255),
      axes = FALSE
    )
    box()
    if (resol == names(cft_fert_bin)[1]) {
      axis(2, mgp = c(2, 0.5, 0), cex.axis = 1.5)
      title(
        ylab = expression(paste("Fertilizer application rate [", kg/m^2, "]")),
        mgp = c(2, 0.5, 0), xpd = NA, cex.lab = 1.5
      )
      legend("topleft", legend = plot_crop_names[cft], cex = 1.5, seg.len = 0)
    }
    if (cft == plot_crops[length(plot_crops)]) {
      axis(1,  mgp = c(2, 1, 0), cex.axis = 1.5)
      title(xlab = "Years", mgp = c(2.5, 1, 0), cex.lab = 1.5, xpd = NA)
    }
    if (cft == plot_crops[1]) {
      title(
        main = sub("arcsec", "'' dataset", sub("arcmin", "' dataset", resol)),
        line = 0.1, cex.main = 1.5
      )
    }
  }
  image.plot(zlim = c(0, area_max * 1e-6), col = area_pal(255), add = TRUE,
             breaks = seq(0, area_max * 1e-6, length.out = 256),
             legend.mar = 0,
             axis.args = list(
               cex.axis = 1.5,
               at = pretty(c(0, area_max * 1e-6)),
               labels = format(pretty(c(0, area_max * 1e-6)), scientific = FALSE)
             ),
             #smallplot = c(1, 1.05, 0.2, 0.8),
             legend.only = TRUE, horizontal = FALSE)
  mtext(expression(paste("Growing area [", km^2, "]")), side = 4, line = 7.5,
        xpd = NA)
}
dev.off()


cat("Areas per N rate (logarithmic scale) saved to paper_nitrogen_rates_crop_log.pdf\n")
area_pal <- colorRampPalette(brewer.pal(9, "YlGn"))
pdf(
  "paper_nitrogen_rates_crop_log.pdf",
  width = 6.5, height = 2 * length(plot_crops) + 0.5, paper = "special",
  pointsize = 8
)
par(
  mfrow = c(length(plot_crops), length(cft_fert_bin)),
  mar = c(0, 0, 1.25, 0.75) + 0.1,
  oma = c(3.5, 3.6, 0, 8) + 0.1
)
for (cft in plot_crops) {
  crop_max <- 0
  area_max <- 0
  area_min <- Inf
  for (resol in names(cft_fert_bin)) {
    yseq <- match(
      seq(first_bin_year, last_bin_year),
      dimnames(cft_fert_max[[resol]])[[2]]
    )
    crop_max <- max(crop_max, max(cft_fert_max[[resol]][cft, yseq]))
    area_max <- max(area_max, max(cft_fert_bin[[resol]][, cft, ], na.rm = TRUE))
    tmparea <- cft_fert_bin[[resol]][, cft, ]
    tmparea[which(tmparea == 0)] <- NA
    area_min <- min(area_min, min(tmparea, na.rm = TRUE))
  }
  area_min <- max(area_min, 1e6)
  area_max <- area_max * 1e-6
  area_min <- area_min * 1e-6
  for (resol in names(cft_fert_bin)) {
    plotdata <- cft_fert_bin[[resol]][, cft, ] * 1e-6
    plotdata[which(plotdata == 0)] <- NA
    plotdata <- log10(plotdata)
    image(
      x = seq(first_bin_year - 0.5, length.out = nyear_bin + 1),
      y = seq(0, crop_max * 10, length.out = bin_n + 1),
      z = t(plotdata),
      zlim = c(log10(area_min), log10(area_max)),
      col = area_pal(255),
      axes = FALSE, useRaster = TRUE
    )
    box()
    if (resol == names(cft_fert_bin)[1]) {
      axis(2, mgp = c(2, 0.5, 0), cex.axis = 1.5)
      title(
        ylab = expression(paste("Fertilizer application rate [", kg/ha, "]")),
        mgp = c(2, 0.5, 0), xpd = NA, cex.lab = 1.5
      )
      legend("topleft", legend = plot_crop_names[cft], cex = 1.5, seg.len = 0)
    }
    if (cft == plot_crops[length(plot_crops)]) {
      axis(1,  mgp = c(2, 1, 0), cex.axis = 1.5)
      title(xlab = "Years", mgp = c(2.5, 1, 0), cex.lab = 1.5, xpd = NA)
    }
    if (cft == plot_crops[1]) {
      title(
        main = sub("arcsec", "'' dataset", sub("arcmin", "' dataset", resol)),
        line = 0.1, cex.main = 1.5
      )
    }
  }
  image.plot(zlim = c(log10(area_min), log10(area_max)), col = area_pal(255),
             breaks = seq(log10(area_min), log10(area_max), length.out = 256),
             legend.mar = 0, add = TRUE,
             axis.args = list(
               cex.axis = 1.5,
               labels = axisTicks(log10(c(area_min, area_max)), log = TRUE, nint = 10),
               at = log10(axisTicks(log10(c(area_min, area_max)), log = TRUE, nint = 10))
               #at = pretty(c(0, area_max * 1e-6)),
               #labels = format(pretty(c(0, area_max * 1e-6)), scientific = FALSE)
             ),
             #smallplot = c(1, 1.05, 0.2, 0.8),
             legend.only = TRUE, horizontal = FALSE)
  mtext(expression(paste("Growing area [", km^2, "]")), side = 4, line = 7.5,
        xpd = NA)
}
dev.off()

## Plot cumulative distributions for areas with certain N application rates   ##
## for time slices                                                            ##
ts_width <- 5
ts_firstyear <- first_bin_year
ts_lastyear <- last_bin_year
# Only include full time slices
ts_start <- seq(ts_firstyear, ts_lastyear - ts_width + 1, by = ts_width)
ts_pal <- colorRampPalette(brewer.pal(11, "Spectral")[-1])
ts_lastyear <- min(ts_lastyear, max(ts_start + ts_width - 1))
probs <- c(0.5, 0.75, 0.9)
legend_space <- 0.6
ts_quantiles <- array(
  dim = c(2, length(probs), length(ts_start), length(plot_crops),
          length(cft_fert_bin)),
  dimnames = list(c("rate", "area"), probs, ts_start, plot_crops,
                  names(cft_fert_bin))
)
pdf(
  "paper_nitrogen_rates_crop_cdf.pdf",
  width = 6.5, height = 2.1 * length(plot_crops), paper = "special",
  pointsize = 8
)
par(
  mfrow = c(length(plot_crops), length(cft_fert_bin)),
  mar = c(0.5, 0, 1.25, 0.75) + 0.1,
  oma = c(3 + legend_space / par("cin")[2] * 1.5, 3.6, 0, 0) + 0.1
)
for (cft in plot_crops) {
  area_sum <- 0
  fert_max <- 0
  for (resol in names(cft_fert_sum)) {
    yseq <- match(
      seq(ts_firstyear, ts_lastyear),
      dimnames(cft_landuse_sum[[resol]])[[2]]
    )
    area_sum <- max(area_sum, cft_landuse_sum[[resol]][cft, yseq] * 1e-6,
                    na.rm = TRUE)
    yseq <- match(
      seq(ts_firstyear, ts_lastyear),
      dimnames(cft_fert_max[[resol]])[[2]]
    )
    fert_max <- max(fert_max, cft_fert_max[[resol]][cft, yseq],
                    na.rm = TRUE)
  }
  if (area_sum > 0) {
    # Always scale to million km2
    fact <- 1e-6
    ylab <- expression(paste("Cumulative growing area [", 10^6, " ", km^2, "]"))
  } else {
    fact <- 1
    ylab <- expression(paste("Cumulative growing area [", km^2, "]"))
  }
  for (resol in names(cft_fert_bin)) {
    plot(1, 1, type = "n", xlim = c(0, fert_max * 10), ylim = c(0, area_sum * fact),
        xlab = "", ylab = "", axes = FALSE)
    if (resol == names(cft_fert_bin)[1]) {
      axis(2, mgp = c(2, 0.5, 0), cex.axis = 1.5)
      title(ylab = ylab,
            mgp = c(2, 0.5, 0), xpd = NA, cex.lab = 1.5)
    }
    if (resol == names(cft_fert_bin)[length(cft_fert_bin)]) {
      legend("bottomright", legend = plot_crop_names[cft], bty = "n", cex = 1.5)
    }
    if (cft == plot_crops[1]) {
      title(
        main = sub("arcsec", "'' dataset", sub("arcmin", "' dataset", resol)),
        line = 0.1, cex.main = 1.5
      )
    }
    # Calculate mean areas across time slice years. Replace NA with 0.
    for (ts in seq_along(ts_start)) {
      yseq <- match(
        seq(ts_start[ts], length.out = ts_width),
        dimnames(cft_fert_bin[[resol]])[[3]]
      )
      ts_data <- cft_fert_bin[[resol]][, cft, yseq]
      ts_data <- pmax(ts_data, 0, na.rm = TRUE)
      ts_data <- apply(ts_data, 1, mean)
      # Set all values after last non-zero value to NA
      if (ts_data[bin_n] == 0) {
        ts_data[seq(max(which(ts_data > 0)) + 1, bin_n)] <- NA
      }
      if (length(unique(cft_fert_bin_thr[[resol]][bin_n, cft, yseq])) != 1) {
        stop("Changing bin rates")
      }
      x <- unique(cft_fert_bin_thr[[resol]][bin_n, cft, yseq]) * bin_seq[-1] * 10
      lines(x = x, y = cumsum(ts_data) * 1e-6 * fact,
            col = ts_pal(length(ts_start))[ts])
      ts_quantiles[, , ts, cft, resol] <- rbind(
        weighted.quantile(x, ts_data, probs = probs),
        sum(ts_data, na.rm = TRUE) * 1e-6 * fact * probs
      )
      points(x = ts_quantiles["rate", , ts, cft, resol],
             y = ts_quantiles["area", , ts, cft, resol],
             pch = seq_along(probs), col = ts_pal(length(ts_start))[ts])
    }
    box()
    axis(1, mgp = c(2, 1, 0), cex.axis = 1.5)
    if (cft == plot_crops[length(plot_crops)]) {
      title(
        xlab = expression(paste("Fertilizer application rate [", kg/ha, "]")),
        mgp = c(2.75, 1, 0), cex.lab = 1.5, xpd = NA
      )
    }
  }
}
par(oma = c(0, 0, 0, 0) + 0.1, mar = c(0, 0, 0, 0) + 0.1,
    fig = c(0, 1, 0 , grconvertY(legend_space, "inches", "ndc")),
    new = TRUE)
plot.new()
legend(
  "topleft",
  legend = paste(ts_start, ts_start + ts_width -1, sep = "-"),
  lty = 1, col = ts_pal(length(ts_start)), cex = 1.5, bty = "n", ncol = 6,
  title = "Time slices", title.adj = 0
)
legend(
  "topright",
  legend = paste(probs * 100, "% quantile"),
  pch = seq_along(probs), bty = "n", cex = 1.5, pt.cex = 1,
  title = "Area-weighted quantiles", title.adj = 0
)
abline(h = 1, xpd = NA)
dev.off()

cat("Application rate [kg/ha] quantiles over time:\n")
print(aperm(ts_quantiles["rate",,,,,drop = FALSE], c(2,3,1,5,4)))

cat("Growing areas [km2] corresponding to rate quantiles:\n")
print(aperm(ts_quantiles["area",,,,,drop = FALSE] / fact, c(2,3,1,5,4)))

## Version with area per bin, not cumulative area                             ##
pdf(
  "paper_nitrogen_rates_crop_ts_sums.pdf",
  width = 6.5, height = 2 * length(plot_crops), paper = "special",
  pointsize = 8
)
par(
  mfrow = c(length(plot_crops), length(cft_fert_bin)),
  mar = c(0.5, 0, 1.25, 0.75) + 0.1,
  oma = c(3, 3.6, 0, 0) + 0.1
)

for (cft in plot_crops) {
  area_sum <- 0
  area_max <- 0
  fert_max <- 0
  for (resol in names(cft_fert_sum)) {
    yseq <- match(
      seq(ts_firstyear, ts_lastyear),
      dimnames(cft_landuse_sum[[resol]])[[2]]
    )
    area_sum <- max(area_sum, cft_landuse_sum[[resol]][cft, yseq] * 1e-6,
                    na.rm = TRUE)
    yseq <- match(
      seq(ts_firstyear, ts_lastyear),
      dimnames(cft_fert_max[[resol]])[[2]]
    )
    fert_max <- max(fert_max, cft_fert_max[[resol]][cft, yseq],
                    na.rm = TRUE)
    yseq <- match(
      seq(ts_firstyear, ts_lastyear),
      dimnames(cft_fert_bin[[resol]])[[3]]
    )
    area_max <- max(
      area_max,
      max(
        apply(
          pmax(
            array(
              cft_fert_bin[[resol]][, cft, yseq],
              dim = c(bin_n, ts_width, length(yseq) / ts_width)
            ),
            0,
            na.rm = TRUE
          ) * 1e-6,
          c(1, 3), mean
        )
      )
    )
  }
  for (resol in names(cft_fert_bin)) {
    plot(1, 1, type = "n", xlim = c(0, fert_max * 10), ylim = c(0, area_max),
        xlab = "", ylab = "", axes = FALSE)
    if (resol == names(cft_fert_bin)[1]) {
      axis(2, mgp = c(2, 0.5, 0), cex.axis = 1.5)
      title(ylab = expression(paste("Growing area [", km^2, "]")),
            mgp = c(2, 0.5, 0), xpd = NA, cex.lab = 1.5)
      legend("topright", legend = plot_crop_names[cft], bty = "n", cex = 1.5)
    }
    if (cft == plot_crops[1]) {
      title(
        main = sub("arcsec", "'' dataset", sub("arcmin", "' dataset", resol)),
        line = 0.1, cex.main = 1.5
      )
    }
    # Calculate mean areas across time slice years. Replace NA with 0.
    for (ts in seq_along(ts_start)) {
      yseq <- match(
        seq(ts_start[ts], length.out = ts_width),
        dimnames(cft_fert_bin[[resol]])[[3]]
      )
      ts_data <- cft_fert_bin[[resol]][, cft, yseq]
      ts_data <- pmax(ts_data, 0, na.rm = TRUE)
      ts_data <- apply(ts_data, 1, mean)
      # Set all values after last non-zero value to NA
      if (ts_data[bin_n] == 0) {
        ts_data[seq(max(which(ts_data > 0)) + 1, bin_n)] <- NA
      }
      if (length(unique(cft_fert_bin_thr[[resol]][bin_n, cft, yseq])) != 1) {
        stop("Changing bin rates")
      }
      x <- unique(cft_fert_bin_thr[[resol]][bin_n, cft, yseq]) * bin_seq[-1] * 10
      lines(x = x, y = ts_data * 1e-6,
            col = ts_pal(length(ts_start))[ts])
    }
    box()
    axis(1, mgp = c(2, 1, 0), cex.axis = 1.5)
    if (cft == plot_crops[length(plot_crops)]) {
      title(
        xlab = expression(paste("Fertilizer application rate [", kg/ha, "]")),
        mgp = c(2.75, 1, 0), cex.lab = 1.5, xpd = NA
      )
    }
  }
}
legend(
  "topright",
  legend = paste(ts_start, ts_start + ts_width -1, sep = "-"),
  lty = 1, col = ts_pal(length(ts_start)), cex = 1.25, bty = "n",
  title = "Time slice"
)
dev.off()
