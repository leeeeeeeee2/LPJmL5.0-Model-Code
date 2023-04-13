################################################################################
## Copyright (C) 2022 Potsdam Institute for Climate Impact Research (PIK),    ##
## see COPYRIGHT file.                                                        ##
##                                                                            ##
## This file is part of LandInG and licensed under GNU AGPL Version 3 or      ##
## later. See LICENSE file or go to http://www.gnu.org/licenses/              ##
## Contact: https://github.com/PIK-LPJmL/LandInG/                             ##
################################################################################

################################################################################
## Script to do analysis of generated land use inputs for documentation paper ##
################################################################################

rm(list = ls())

library(geosphere)

################################################################################
## Datasets to compare                                                        ##
datasets <- list(
  "5arcmin" = c("hydrosheds_v1.1", "merit"),
  "30arcmin" = c("ddm30", "stn")
)
# River routing input filenames
drainname <- list()
for (resol in names(datasets)) {
  drainname[[resol]] <- paste0("drainage_", datasets[[resol]], "_", resol, ".bin")
  names(drainname[[resol]]) <- datasets[[resol]]
}
# Cell lists and upstream areas
celllists_RData <- upstream_RData <- list()
for (resol in names(datasets)) {
  celllists_RData[[resol]] <- paste0(
    "drainage_celllists_", datasets[[resol]], "_", resol, ".RData"
  )
  names(celllists_RData[[resol]]) <- datasets[[resol]]
  upstream_RData[[resol]] <- paste0(
    "drainage_upstreamarea_", datasets[[resol]], "_", resol, ".RData"
  )
  names(upstream_RData[[resol]]) <- datasets[[resol]]
}
# Options for irrigation neighbour setting. Include "self" for upstream areas
# of cells themselves.
neighbour_upstream_options <- c(
  "self",
  "adjacent",
  "adjacent_exclude_downstream_exclude_upstream",
  "75000m_radius_exclude_downstream_exclude_upstream",
  "75000m_radius_exclude_downstream_exclude_upstream_idw"
)
neighbour_upstream_option_names <- list(
  "Median upstream area of cell",
  "Adjacent",
  c("Adjacent", "(no upstream/downstream cell)"),
  c("75 km region", "(no upstream/downstream cell)"),
  c("75 km region", "(no u./d. cell, distance-weighted)")
)
names(neighbour_upstream_option_names) <- neighbour_upstream_options
## River basin outlets (rough coordinates, will be slightly different for     ##
## each dataset, used to assign outlet cells in river routing to basins)
river_outlet <- array(
  dim = c(5, 2),
  dimnames = list(
    c("Amazon Basin", "Congo Basin", "Mississippi Basin", "Nile Basin",
      "Rio de la Plata Basin"),
    c("lon", "lat")
  )
)
river_outlet["Amazon Basin", ] <- c(-50.09, 0.72)
river_outlet["Congo Basin", ] <- c(12.36, -6.05)
river_outlet["Mississippi Basin", ] <- c(-89.25, 29.15)
river_outlet["Nile Basin", ] <- c(31.15, 30.17)
river_outlet["Rio de la Plata Basin", ] <- c(-55.78, -35.67)
################################################################################

## Load helper functions to work with LPJmL format                            ##
if (file.exists(file.path("..", "lpjml_format_helper_functions.R"))) {
  source(file.path("..", "lpjml_format_helper_functions.R"))
} else {
  stop("Please update path to script with LPJmL input format helper function")
}

################################################################################
## Load river routing data                                                    ##
draindata <- list()
for (resol in names(datasets)) {
  draindata[[resol]] <- list()
  for (ds in datasets[[resol]]) {
    tmpheader <- read_header(drainname[[resol]][ds])
    zz <- file(drainname[[resol]][ds], "rb")
    seek(zz, get_headersize(tmpheader))
    draindata[[resol]][[ds]] <- matrix(
      readBin(
        zz,
        what = get_datatype(tmpheader)$type,
        size = get_datatype(tmpheader)$size,
        n = prod(tmpheader$header[c("ncell", "nbands")]),
        endian = tmpheader$endian
      ) * tmpheader$header["scalar"],
      ncol = tmpheader$header["nbands"],
      byrow = tmpheader$header["nbands"] != 4
    )
    close(zz)
    tmpenv <- new.env()
    load(celllists_RData[[resol]][ds], envir = tmpenv)
    # Reduce memory requirement by removing unneeded data
    tmpenv$drainage_dsclist <- NULL
    load(upstream_RData[[resol]][ds], envir = tmpenv)
    envname <- paste0("env_", ds, "_", resol)
    assign(envname, tmpenv)
  }
}
## % of grid area covered by DDM, median basin size, number of basins with    ##
## > 1 cells                                                                  ##
basic_attr_table <- data.frame(
  parameter = c("\\% grid area covered by DDM", "Median river basin size ",
                 "River basins with at least 2 cells")
)
for (resol in names(datasets)) {
  for (ds in datasets[[resol]]) {
    column <- character(3)
    non_missing <- which(draindata[[resol]][[ds]][, 1] != -9)
    tmparea <- get(paste0("env_", ds, "_", resol))$drainage_gridarea
    column[1] <- round(sum(tmparea[non_missing]) / sum(tmparea) * 100, 1)
    outlet <- which(draindata[[resol]][[ds]][, 1] < 0)
    column[2] <- paste0(
      round(
        median(
          get(paste0("env_", ds, "_", resol))$drainage_upstreamarea[outlet]
        ) * 1e-6,
      ),
      "~\\unit{km^2}"
    )
    column[3] <- length(
      which(
        sapply(get(paste0("env_", ds, "_", resol))$drainage_usclist[outlet],
               length) > 0
      )
    )
    basic_attr_table <- cbind(basic_attr_table, column)
    colnames(basic_attr_table)[ncol(basic_attr_table)] <- paste0(ds, "_", resol)
  }
}

## size and outlet of 10 largest river basins                                 ##
outlet_table <- array(
  dim = c(10, 3, length(unlist(datasets))),
  dimnames = list(
    NULL,
    c("lon", "lat", "upstream"),
    paste(
      unlist(datasets),
      rep(names(datasets), times = sapply(datasets, length)),
      sep = "_"
    )
  )
)
for (resol in names(datasets)) {
  for (ds in datasets[[resol]]) {
    tmpenv <- get(paste0("env_", ds, "_", resol))
    sequence <- intersect(
      order(tmpenv$drainage_upstreamarea, decreasing = TRUE),
      which(draindata[[resol]][[ds]][, 1] < 0)
    )[seq_len(10)]
    outlet_table[, , paste(ds, resol, sep = "_")] <- cbind(
      tmpenv$drainage_griddata,
      tmpenv$drainage_upstreamarea
    )[sequence, ]
  }
}
outlet_table_paper <- data.frame(
  parameter = c(rbind(rownames(river_outlet), ""))
)
for (ds in dimnames(outlet_table)[[3]]) {
  column <- character(nrow(outlet_table_paper))
  for (river in rownames(river_outlet)) {
    # Find closest match
    r <- which.min(rowSums(abs(outlet_table[, c("lon", "lat"), ds] - rep(
      river_outlet[river, ], each = 10))))
    column[match(river, outlet_table_paper[, "parameter"])] <- paste0(
      "$",
      round(outlet_table[r, "upstream", ds] * 1e-12, 2),
      " \\cdot 10^6$~\\unit{km^2}"
    )
    column[match(river, outlet_table_paper[, "parameter"]) + 1] <- paste0(
      "($", round(abs(outlet_table[r, "lon", ds]), 2), "^\\circ$",
      ifelse(outlet_table[r, "lon", ds] < 0, "W", "E"),
      ", $", round(abs(outlet_table[r, "lat", ds]), 2), "^\\circ$",
      ifelse(outlet_table[r, "lat", ds] < 0, "S", "N"),
      ")"
    )
  }
  outlet_table_paper <- cbind(
    outlet_table_paper,
    column
  )
  colnames(outlet_table_paper)[ncol(outlet_table_paper)] <- ds
}
################################################################################


################################################################################
## Load neighbour irrigation data                                             ##
neighbour_upstream_table <- array(
  dim = c(length(neighbour_upstream_options), length(unlist(datasets))),
  dimnames = list(
    neighbour_upstream_options,
    paste(unlist(datasets), rep(names(datasets), sapply(datasets, length)),
          sep = "_")
  )
)
neighbour_distance_table <- neighbour_bigger_table <- neighbour_upstream_table
for (opt in setdiff(dimnames(neighbour_upstream_table)[[1]], "self")) {
  for (ds in dimnames(neighbour_upstream_table)[[2]]) {
    filename <- paste0(
      "neighbour_irrig_", ds, "_", opt, ".bin"
    )
    if (file.exists(filename)) {
      tmpheader <- read_header(filename)
      zz <- file(filename, "rb")
      seek(zz, get_headersize(tmpheader))
      filedata <- matrix(
        readBin(
          zz,
          what = get_datatype(tmpheader)$type,
          size = get_datatype(tmpheader)$size,
          n = tmpheader$header["ncell"] * tmpheader$header["nbands"]
        ) * tmpheader$header["scalar"],
        ncol = tmpheader$header["nbands"],
        byrow = tmpheader$header["order"] == 1
      )
      assign(
        paste0("neighbourdata_", ds, "_", opt),
        filedata
      )
      close(zz)
      neighbour_upstream_table["self", ds] <- median(
        get(paste0("env_", ds))[["drainage_upstreamarea"]]
      ) * 1e-6
      neighbour_upstream_table[opt, ds] <- median(
        get(paste0("env_", ds))[["drainage_upstreamarea"]][filedata + 1]
      ) * 1e-6
      neighbour_distance_table[opt, ds] <- median(
        distHaversine(
          get(paste0("env_", ds))[["drainage_griddata"]],
          get(paste0("env_", ds))[["drainage_griddata"]][filedata + 1, ],
          r = earthradius
        )
      ) * 1e-3
      neighbour_bigger_table[opt, ds] <- length(
        which(
          get(paste0("env_", ds))[["drainage_upstreamarea"]][filedata + 1] >
            get(paste0("env_", ds))[["drainage_upstreamarea"]]
        )
      ) / length(filedata)
      rm(filedata)
    }
  }
}
neighbour_table_paper <- matrix(
  paste0(
    round(neighbour_upstream_table["self", ]),
    "~\\unit{km^2}"
  ),
  nrow = 1,
  dimnames = list(
    neighbour_upstream_option_names[["self"]][1],
    colnames(neighbour_upstream_table)
  )
)
if (length(neighbour_upstream_option_names[["self"]]) > 1) {
  for (r in seq(2, length(neighbour_upstream_option_names[["self"]]))) {
    neighbour_table_paper <- rbind(
      neighbour_table_paper,
      matrix("", nrow = 1, ncol = ncol(neighbour_table_paper))
    )
    rownames(neighbour_table_paper)[nrow(neighbour_table_paper)] <-
      neighbour_upstream_option_names[["self"]][r]
  }
}
for (opt in setdiff(dimnames(neighbour_upstream_table)[[1]], "self")) {
  neighbour_table_paper <- rbind(
    neighbour_table_paper,
    matrix(
      paste0(
        round(neighbour_upstream_table[opt, ]),
        "~\\unit{km^2} / ",
        round(neighbour_distance_table[opt, ]),
        "~\\unit{km} / ",
        round(neighbour_bigger_table[opt, ] * 100),
        "~\\%"
      ),
      nrow = 1,
      dimnames = list(
        neighbour_upstream_option_names[[opt]][1],
        colnames(neighbour_upstream_table)
      )
    )
  )
  if (length(neighbour_upstream_option_names[[opt]]) > 1) {
    for (r in seq(2, length(neighbour_upstream_option_names[[opt]]))) {
      neighbour_table_paper <- rbind(
        neighbour_table_paper,
        matrix("", nrow = 1, ncol = ncol(neighbour_table_paper))
      )
      rownames(neighbour_table_paper)[nrow(neighbour_table_paper)] <-
        neighbour_upstream_option_names[[opt]][r]
    }
  }
}
################################################################################


################################################################################
## Write tables to text file in LaTeX syntax                                  ##
write.table(
  format(
    rbind(
      basic_attr_table,
      outlet_table_paper,
      data.frame(
        parameter = rownames(neighbour_table_paper),
        neighbour_table_paper
      ),
      make.row.names = FALSE
    ),
    justify = "right"
  ),
  file = "paper_river_routing_table.txt",
  sep = " & ", eol = "\\\\\n", quote = FALSE,
  row.names = FALSE, col.names = FALSE
)
