################################################################################
## Copyright (C) 2022 Potsdam Institute for Climate Impact Research (PIK),    ##
## see COPYRIGHT file.                                                        ##
##                                                                            ##
## This file is part of LandInG and licensed under GNU AGPL Version 3 or      ##
## later. See LICENSE file or go to http://www.gnu.org/licenses/              ##
## Contact: https://github.com/PIK-LPJmL/LandInG/                             ##
################################################################################

################################################################################
## Script to do analysis of generated grid and administrative masks for       ##
## documentation paper                                                        ##
################################################################################

rm(list = ls())

################################################################################
## Set up file names                                                          ##
gridname <- c(
  "5arcmin" = "../gadm_paper/grid_gadm_5arcmin.bin",
  "30arcmin" = "../gadm_paper/grid_gadm_30arcmin_predefined_grid.bin",
  old_version = file.path(
    "/p", "projects", "lpjml", "input", "historical", "input_VERSION2",
    "grid.bin"
  )
)
soilname <- c(
  "5arcmin" = "../soil_paper/soil_5arcmin_13_types.bin",
  "30arcmin" = "../soil_paper/soil_30arcmin_13_types.bin",
  old_version = file.path(
    "/p", "projects", "lpjml", "input", "historical", "input_VERSION2",
    "soil_new_67420.bin"
  )
)
landfrac_name <- c(
  "5arcmin" = "../gadm_paper/landfrac_gadm_5arcmin.bin",
  "30arcmin" = "../gadm_paper/landfrac_gadm_30arcmin_predefined_grid.bin",
  # Not available for Schaphoff18 dataset, use 30arcmin
  old_version = "../gadm_paper/landfrac_gadm_30arcmin_predefined_grid.bin"
)
soil_aggregation_name <- c(
  "5arcmin" = "../soil_paper/soil_aggregation_5arcmin_13_types.RData",
  "30arcmin" = "../soil_paper/soil_aggregation_30arcmin_13_types.RData",
  old_version = NULL
)
## Soil texture types in LPJmL                                                ##
lpjml_soiltypes <- c(
  "clay", "silty clay", "sandy clay", "clay loam", "silty clay loam",
  "sandy clay loam", "loam", "silt loam", "sandy loam", "silt", "loamy sand",
  "sand", "rock and ice"
)
################################################################################


################################################################################
## Helper functions for LPJmL input format                                    ##
## The script lpjml_format_helper_functions.R is saved in the parent          ##
## directory by default.                                                      ##
if (file.exists("../lpjml_format_helper_functions.R")) {
  source("../lpjml_format_helper_functions.R")
} else {
  stop("Please update path to script with LPJmL input format helper function")
}
################################################################################


################################################################################
## Load data from files                                                       ##
griddata <- gridarea <- soildata <- landfrac_data <- soil_aggregation_data <-
  list()

for (resol in names(gridname)) {
  tmpheader <- read_header(gridname[resol])
  zz <- file(gridname[resol], "rb")
  seek(zz, get_headersize(tmpheader))
  griddata[[resol]] <- matrix(
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
  gridarea[[resol]] <- cellarea(
    griddata[[resol]][, 2],
    tmpheader$header["cellsize_lon"],
    tmpheader$header["cellsize_lat"]
  )
  if (resol %in% names(landfrac_name)) {
    tmpheader <- read_header(landfrac_name[resol])
    zz <- file(landfrac_name[resol], "rb")
    seek(zz, get_headersize(tmpheader))
    landfrac_data[[resol]] <- matrix(
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
  }
  # Soil input for Schaphoff18 is without header
  if (resol != "old_version" &&  resol %in% names(soilname)) {
    tmpheader <- read_header(soilname[resol])
    zz <- file(soilname[resol], "rb")
    seek(zz, get_headersize(tmpheader))
  } else if (resol %in% names(soilname)) {
    tmpheader$header["datatype"] <- 0
    tmpheader$header["nbands"] <- 1
    zz <- file(soilname[resol], "rb")
  }
  if (typeof(get_datatype(tmpheader)$type) == "raw") {
    # Need to convert raw into integer
    soildata[[resol]] <- matrix(
      as.integer(
        readBin(
          zz,
          what = get_datatype(tmpheader)$type,
          size = get_datatype(tmpheader)$size,
          n = prod(tmpheader$header[c("ncell", "nbands")]),
          endian = tmpheader$endian
        )
      ) * tmpheader$header["scalar"],
      ncol = tmpheader$header["nbands"],
      byrow = tmpheader$header["order"] != 4
    )
  } else {
    soildata[[resol]] <- matrix(
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
  }
  close(zz)
  if (resol %in% names(soil_aggregation_name)) {
    tmpenv <- new.env()
    load(soil_aggregation_name[resol], envir = tmpenv)
    soil_aggregation_data[[resol]] <- tmpenv
  }
}
################################################################################


################################################################################
## Calculate global area share of each soil type in each dataset              ##
soil_texture_share <- data.frame(
  parameter = paste("\\%", lpjml_soiltypes)
)
for (resol in names(soildata)) {
  share_grid <- share_land <- share_n <- numeric(length(lpjml_soiltypes))
  for (tex in seq_along(lpjml_soiltypes)) {
    tex_cells <- which(soildata[[resol]] == tex)
    share_grid[tex] <- sum(gridarea[[resol]][tex_cells]) /
      sum(gridarea[[resol]]) * 100
    share_land[tex] <- sum(
      (gridarea[[resol]] * landfrac_data[[resol]])[tex_cells]
    ) / sum(gridarea[[resol]] * landfrac_data[[resol]]) * 100
    share_n[tex] <- length(tex_cells) / length(soildata[[resol]]) * 100
  }
  soil_texture_share <- cbind(
    soil_texture_share,
    share_n,
    share_grid,
    share_land
  )
  col <- seq(to = ncol(soil_texture_share), length.out = 3)
  colnames(soil_texture_share)[col] <- paste0(
    c("cell_share", "gridarea_share", "landarea_share"), "_", resol
  )
}

# Use grid area share for paper and round to 1 digit
col <- grep("gridarea_share", colnames(soil_texture_share))
paper_table <- data.frame(
  parameter = soil_texture_share$parameter,
  round(soil_texture_share[, col], 1)
)
## Number of cells gap-filled in each dataset                                 ##
gapfilled_table <- data.frame(parameter = "Number of cells gap-filled")
for (resol in names(soildata)) {
  if (!is.null(soil_aggregation_data[[resol]])) {
    gapfilled_table <- cbind(
      gapfilled_table,
      as.character(
        nrow(soil_aggregation_data[[resol]][["tmp_lpjml_soiltexture_gapfilling"]])
      )
    )
  } else {
    gapfilled_table <- cbind(gapfilled_table, "---")
  }
  colnames(gapfilled_table)[ncol(gapfilled_table)] <- grep(
    resol,
    colnames(paper_table),
    value = TRUE
  )
}
# Order soil types by decreasing area share
if (exists("order_table"))
  rm(order_table)
for (ds in grep("parameter", colnames(paper_table), value = TRUE, invert = TRUE)) {
  if (exists("order_table")) {
    order_table <- cbind(
      order_table,
      paper_table$parameter[order(paper_table[, ds], decreasing = TRUE)]
    )
  } else {
    order_table <- data.frame(
      paper_table$parameter[order(paper_table[, ds], decreasing = TRUE)]
    )
  }
  colnames(order_table)[ncol(order_table)] <- ds
}
# Save paper_table to file in LaTeX format
options(width = 160)
cat("Table written to paper_soil_table.txt:\n")
write.table(
  format(
    rbind(format(paper_table, nsmall = 1), gapfilled_table),
    justify = "right"
  ),
  file = "paper_soil_table.txt",
  sep = " & ", quote = FALSE, eol = " \\\\\n",
  col.names = FALSE, row.names = FALSE
)
print(
  format(
    rbind(format(paper_table, nsmall = 1), gapfilled_table),
    justify = "right"
  )
)
cat("Soil texture class ordered by decreasing area share:\n")
print(order_table)
################################################################################
