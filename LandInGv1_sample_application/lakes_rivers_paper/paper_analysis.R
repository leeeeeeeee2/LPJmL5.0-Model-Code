################################################################################
## Copyright (C) 2022 Potsdam Institute for Climate Impact Research (PIK),    ##
## see COPYRIGHT file.                                                        ##
##                                                                            ##
## This file is part of LandInG and licensed under GNU AGPL Version 3 or      ##
## later. See LICENSE file or go to http://www.gnu.org/licenses/              ##
## Contact: https://github.com/PIK-LPJmL/LandInG/                             ##
################################################################################

################################################################################
## Script to do analysis of generated lake and river inputs for documentation ##
## paper                                                                      ##
################################################################################

rm(list = ls())
library(udunits2)

################################################################################
## Set up file names                                                          ##
cft_gridname <- c(
  "5arcmin" = "../gadm_paper/grid_gadm_5arcmin.bin",
  "30arcmin" = "../gadm_paper/grid_gadm_30arcmin_predefined_grid.bin",
  "old_version" = file.path(
    "/p", "projects", "lpjml", "input", "historical", "input_VERSION2",
    "grid.bin"
  )
)
lake_mask_options <- c(
  default = "", polygon = "polygon-based",
  # test versions run on a different version of sf package:
  nos2 = "polygon-based_test_nos2", s2 = "polygon-based_test_s2"
)
lake_mask_basename <- "glwd_lakes_and_rivers_"
## Lake and river input from previous version
lake_mask_old <- file.path(
  "/p", "projects", "lpjml", "input", "historical", "input_VERSION2",
  "glwd_lakes_and_rivers.bin"
)
lake_mask_old_has_header <- FALSE
lake_mask_old_grid <- "old_version"
table_units <- "1e6 km2"
################################################################################

## Load helper functions to work with LPJmL format                            ##
if (file.exists(file.path("..", "lpjml_format_helper_functions.R"))) {
  source(file.path("..", "lpjml_format_helper_functions.R"))
} else {
  stop("Please update path to script with LPJmL input format helper function")
}

################################################################################
## Load data                                                                  ##
cft_griddata <- cft_gridheader <- lakefrac_data <- list()
for (resol in names(cft_gridname)) {
  tmpheader <- read_header(cft_gridname[resol])
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
  cft_gridheader[[resol]] <- tmpheader
  lakefrac_data[[resol]] <- list()
  for (opt in names(lake_mask_options)) {
    filename <- paste0(
      lake_mask_basename,
      ifelse(
        nchar(lake_mask_options[opt]) > 0,
        paste0(lake_mask_options[opt], "_"),
        ""
      ),
      resol,
      ".bin"
    )
    if (!file.exists(filename)) {
      warning("File ", filename, " does not exist.",
              call. = FALSE, immediate. = TRUE)
      next
    }
    tmpheader <- read_header(filename)
    zz <- file(filename, "rb")
    seek(zz, get_headersize(tmpheader))
    lakefrac_data[[resol]][[opt]] <- matrix(
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
}

if (lake_mask_old_has_header) {
  tmpheader <- read_header(lake_mask_old)
  zz <- file(lake_mask_old, "rb")
  seek(zz, get_headersize(tmpheader))
} else {
  tmpheader <- cft_gridheader[[lake_mask_old_grid]]
  tmpheader$header["nbands"] <- 1
  tmpheader$header["datatype"] <- 0
  tmpheader$header["scalar"] <- 1/100
  zz <- file(lake_mask_old, "rb")
}
if (typeof(get_datatype(tmpheader)$type) == "raw") {
  lakefrac_data_old <- matrix(
    as.integer(
      readBin(
        zz,
        what = get_datatype(tmpheader)$type,
        size = get_datatype(tmpheader)$size,
        n = prod(tmpheader$header[c("ncell", "nbands")])
      )
    ) * tmpheader$header["scalar"],
    ncol = tmpheader$header["nbands"],
    byrow = tmpheader$header["order"] != 4
  )
} else {
  lakefrac_data_old <- matrix(
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

global_lake_area <- data.frame(parameter = "Global lake and river area")
for (resol in names(lakefrac_data)) {
  for (opt in names(lakefrac_data[[resol]])) {
    global_lake_area <- cbind(
      global_lake_area,
      sum(lakefrac_data[[resol]][[opt]] * cellarea(
        cft_griddata[[resol]][, 2],
        cft_gridheader[[resol]]$header["cellsize_lon"],
        cft_gridheader[[resol]]$header["cellsize_lat"]
      )) * ud.convert(1, "m2", table_units)
    )
    colnames(global_lake_area)[ncol(global_lake_area)] <- paste0(opt, "_", resol)
  }
}
global_lake_area <- cbind(
  global_lake_area,
  sum(lakefrac_data_old * cellarea(
    cft_griddata[[resol]][, 2],
    cft_gridheader[[lake_mask_old_grid]]$header["cellsize_lon"],
    cft_gridheader[[lake_mask_old_grid]]$header["cellsize_lat"]
  )) * ud.convert(1, "m2", table_units)
)
colnames(global_lake_area)[ncol(global_lake_area)] <- ifelse(
  grepl("old", lake_mask_old_grid),
  lake_mask_old_grid,
  paste0("old_", lake_mask_old_grid)
)
latex_units <- ""
latex_units_mathmode <- FALSE
if (grepl("\\d+e\\d+", table_units)) {
  latex_units_mathmode <- TRUE
  latex_units <- " \\cdot "
  scalar <- regmatches(
    table_units,
    regexec("(\\d+)e(\\d+)", table_units)
  )[[1]][-1]
  if (scalar[1] != "1") {
    latex_units <- paste0(latex_units, scalar[1], " \\cdot ")
  }
  latex_units <- paste0(latex_units, "10^{", scalar[2], "}$~")
}
if (grepl("([[:alpha:]]+)(\\d*)$", table_units)) {
  unit_string <- regmatches(
    table_units,
    regexec("([[:alpha:]]+)(\\d*)$", table_units)
  )[[1]][-1]
#  if (length(unit_string) > 1)
#    latex_units_mathmode <- TRUE
  if (nchar(latex_units) > 0) {
    latex_units <- paste0(
      latex_units,
      "\\unit{", unit_string[1],
#      ifelse(
#        latex_units_mathmode,
#        paste0("\\mathrm{", unit_string[1], "}"),
#        unit_string[1]
#      ),
      ifelse(length(unit_string) > 1, paste0("^{", unit_string[2], "}"), ""),
      "}"
    )
  } else {
    latex_units <- paste0(
      "~",
      "\\unit{", unit_string[1],
#      ifelse(
#        latex_units_mathmode,
#        paste0("\\mathrm{", unit_string[1], "}"),
#        unit_string[1]
#      ),
      ifelse(length(unit_string) > 1, paste0("^{", unit_string[2], "}"), ""),
      "}"
    )
  }
}

paper_table <- data.frame(
  parameter = global_lake_area$parameter,
  matrix(
    paste0(
      ifelse(latex_units_mathmode, "$", ""),
      round(
        global_lake_area[, grep("default|old_", colnames(global_lake_area))],
        3
      ),
      latex_units#,
      #ifelse(latex_units_mathmode, "$", "")
    ),
    nrow = 1,
    dimnames = list(
      NULL, grep("default|old_", colnames(global_lake_area), value = TRUE)
    )
  )
)

write.table(
  paper_table,
  file = "paper_lakes_table.txt",
  sep = " & ", eol = "\\\\\n", quote = FALSE,
  row.names = FALSE, col.names = FALSE
)
