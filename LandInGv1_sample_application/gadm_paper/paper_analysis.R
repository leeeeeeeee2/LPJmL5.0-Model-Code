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
# old_version refers to Schaphoff18 dataset
gridname <- c(
  "5arcmin" = "grid_gadm_5arcmin.bin",
  "30arcmin" = "grid_gadm_30arcmin_predefined_grid.bin",
  old_version = file.path(
    "/p", "projects", "lpjml", "input", "historical", "input_VERSION2",
    "grid.bin"
  )
)
celllist_name <- c(
  "5arcmin" = "celllist_5arcmin.RData",
  "30arcmin" = "celllist_30arcmin_predefined_grid.RData",
  old_version = NULL
)
landfrac_name <- c(
  "5arcmin" = "landfrac_gadm_5arcmin.bin",
  "30arcmin" = "landfrac_gadm_30arcmin_predefined_grid.bin",
  # Not available for Schaphoff18 dataset, use 30arcmin
  old_version = "landfrac_gadm_30arcmin_predefined_grid.bin"
)
cowname <- c(
  "5arcmin" = "cow_gadm_5arcmin.bin",
  "30arcmin" = "cow_gadm_30arcmin_predefined_grid.bin",
  old_version = file.path(
    "/p", "projects", "lpjml", "input", "historical", "input_VERSION2",
    "cow_full_2018.bin"
  )
)
# Meta data not available for Schaphoff18 dataset
cowmetaname <- c(
  "5arcmin" = "cow_gadm_5arcmin_countries.csv",
  "30arcmin" = "cow_gadm_30arcmin_predefined_grid_countries.csv",
  old_version = NULL
)
regmetaname <- c(
  "5arcmin" = "cow_gadm_5arcmin_regions.csv",
  "30arcmin" = "cow_gadm_30arcmin_predefined_grid_regions.csv",
  old_version = NULL
)
## Is the grid pre-defined?
predefined_grid <- c("5arcmin" = FALSE, "30arcmin" = TRUE, old_version = TRUE)
## The country code input for LPJmL includes regions for some large           ##
## countries. Set up which countries (3-letter ISO code):                     ##
include_regions <- c("AUS", "BRA", "CAN", "CHN", "IND", "RUS", "USA")
## Country code to be used for cells without GADM coverage (if                ##
## force_grid == TRUE). Set variable value to country name to be used, and    ##
## name attribute to ISO code. The ISO code should consist of three capital   ##
## letters and start with an "X" to denote it is not an official ISO code.    ##
gadm_no_land <- c(XNL = "No land")
## Code in Schaphoff18 dataset corresponding to "No land" entry. Taken from   ##
## LPJmL source code.                                                         ##
old_version_no_land <- 146
## The default for global LPJmL setups is to exclude Antarctica (ATA).        ##
## You may add more countries to exclude by adding their ISO codes.           ##
## Set to NULL if you do not want to skip any countries.                      ##
skip_countries <- c("ATA")
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
griddata <- gridarea <- celllist <- landfrac_data <- cowdata <- cowmetadata <-
  regmetadata <- list()

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
  if (resol %in% names(cowname)) {
    tmpheader <- read_header(cowname[resol])
    zz <- file(cowname[resol], "rb")
    seek(zz, get_headersize(tmpheader))
    cowdata[[resol]] <- matrix(
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
  if (resol %in% names(cowmetaname[resol])) {
    cowmetadata[[resol]] <- read.csv(cowmetaname[resol])
  }
  if (resol %in% names(regmetaname[resol])) {
    regmetadata[[resol]] <- read.csv(regmetaname[resol])
  }
  if (resol %in% names(celllist_name)) {
    tmpenv <- new.env()
    load(celllist_name[resol], envir = tmpenv)
    if (predefined_grid[resol]) {
      tmpenv$predefined_ngrid <- nrow(griddata[[resol]])
    } else {
      tmpenv$predefined_ngrid <- Inf
    }
    celllist[[resol]] <- tmpenv
  }
}
################################################################################


################################################################################
## Comparison table                                                           ##
comparison_table <- data.frame(
  parameter = c(
    "Number of cells in grid",
    "Total grid area [$10^6$~\\unit{km^2}]",
    "Total land area [$10^6$~\\unit{km^2}]",
    "Number of countries assigned",
    "Number of cells assigned as `No land'",
    "Number of regions$^\\ddagger$ assigned"
  ),
  stringsAsFactors = FALSE
)
tot_countries <- tot_regions <- 0
for (resol in names(griddata)) {
  col_vals <- numeric(nrow(comparison_table))
  names(col_vals) <- c("ngrid", "gridarea", "landarea", "ncountries",
                       "n_noland", "nregions")
  col_vals["ngrid"] <- nrow(griddata[[resol]])
  col_vals["gridarea"] <- round(sum(gridarea[[resol]]) * 1e-12, 1)
  if (resol %in% names(landfrac_data)) {
    col_vals["landarea"] <- round(
      sum(gridarea[[resol]] * landfrac_data[[resol]]) * 1e-12,
      1
    )
  } else {
    col_vals["landarea"] <- NA
  }
  if (resol %in% names(cowdata)) {
    col_vals["ncountries"] <- length(unique(cowdata[[resol]][, 1]))
    if (resol %in% names(cowmetadata)) {
      tot_countries <- max(
        tot_countries,
        nrow(cowmetadata[[resol]]) + length(
          which(!skip_countries %in% cowmetadata[[resol]]$ISO)
        ) + length(
          which(!names(gadm_no_land) %in% cowmetadata[[resol]]$ISO)
        )
      )
      # Codes of countries without regions
      index <- which(!cowmetadata[[resol]]$ISO %in% include_regions)
      no_region_codes <- cowmetadata[[resol]]$ID[index]
      region_codes <- cowmetadata[[resol]]$ID[-index]
      # No land
      index <- which(cowmetadata[[resol]]$ISO %in% names(gadm_no_land))
      no_land_id <- cowmetadata[[resol]]$ID[index]
    } else {
      # If now meta data is available assume that countries with regions are
      # last in ID list. Assume that data uses same number of countries as
      # include_regions
      region_codes <- seq(
        to = max(cowdata[[resol]][, 1]),
        length.out = length(include_regions)
      )
      no_region_codes <- seq(0, min(region_codes) - 1)
      no_land_id <- old_version_no_land
    }
    col_vals["n_noland"] <- length(which(cowdata[[resol]][, 1] %in% no_land_id))
    col_vals["nregions"] <- length(
      unique(setdiff(cowdata[[resol]][, 2], no_region_codes))
    )
    if (resol %in% names(regmetadata)) {
      tot_regions <- max(
        tot_regions,
        length(unique(setdiff(regmetadata[[resol]]$ID, no_region_codes)))
      )
    }
  }
  if (col_vals["n_noland"] > 0) {
    col_vals["landarea"] <- paste0(col_vals["landarea"], "$^*$")
  }
  if (!resol %in% names(cowmetadata)) {
    col_vals[c("ncountries", "n_noland")] <- paste0(
      col_vals[c("ncountries", "n_noland")], "$^\\dagger$"
    )
  }
  if (!resol %in% names(regmetadata)) {
    col_vals["nregions"] <- paste0(col_vals["nregions"], "$^\\dagger$")
  }
  if (is.numeric(col_vals)) {
    col_vals[c("ngrid", "ncountries", "n_noland", "nregions")] <- paste0(
      col_vals[c("ngrid", "ncountries", "n_noland", "nregions")]
    )
    col_vals[c("gridarea", "landarea")] <- format(
      as.numeric(col_vals[c("gridarea", "landarea")]),
      nsmall = 1
    )
  }
  comparison_table <- cbind(comparison_table, unname(col_vals))
  colnames(comparison_table)[ncol(comparison_table)] <- resol
}
r <- grep("Number of countries", comparison_table$parameter)
comparison_table$parameter[r] <- paste0(
  comparison_table$parameter[r], " (total: ", tot_countries, ")"
)
r <- grep("Number of regions", comparison_table$parameter)
comparison_table$parameter[r] <- paste0(
  comparison_table$parameter[r], " (total: ", tot_regions, ")"
)
# Write data to file in LaTeX syntax
cat("Comparison table saved to paper_analysis_grid_admin.txt:\n")
write.table(
  format(comparison_table, justify = "right"),
  file = "paper_analysis_grid_admin.txt",
  sep = " & ", quote = FALSE, eol = " \\\\\n",
  row.names = FALSE, col.names = FALSE
)
print(comparison_table)
################################################################################


################################################################################
## Additional information about grid intersection                             ##
# Number of cells with land according to grid intersection
cat("Number of cells with land according to grid intersection:\n")
print(
  sapply(
    celllist,
    function(listdata) length(which(!sapply(listdata$cell_list, is.null)))
  )
)

landarea_celllist <- function(celldata) {
  if (is.null(celldata))
    return(NA)
  return(celldata["Landarea"])
}

gridid_celllist <- function(celldata) {
  if (is.null(celldata))
    return(NA)
  return(celldata["GridID"])
}

cat("Number of cells with < 1000 m2 land according to grid intersection:\n")
print(
  sapply(
    celllist,
    function(listdata) length(
      which(sapply(listdata$cell_list, landarea_celllist) < 1000)
    )
  )
)
cat("Total land area in cells with < 1000 m2 land area:\n")
print(
  sapply(
    celllist,
    function(listdata) sum(
      ifelse(
        sapply(listdata$cell_list, landarea_celllist) < 1000,
        sapply(listdata$cell_list, landarea_celllist),
        0
      ),
      na.rm = TRUE
    )
  )
)


cat("Number of cells with land according to grid intersection, but not",
    "included due to pre-defined grid:\n")
print(
  sapply(
    celllist,
    function(listdata) length(
      which(
        sapply(listdata$cell_list, gridid_celllist) > listdata$predefined_ngrid
      )
    )
  )
)
cat("Countries with regional distinction:\n")
cat(
  toString(
    unique(
      c(
        unlist(
          sapply(
            cowmetadata,
            function(cowmeta, include)
              cowmeta$country[which(cowmeta$ISO %in% include)],
            include = include_regions
          )
        )
      )
    )
  ),
  "\n"
)
################################################################################
