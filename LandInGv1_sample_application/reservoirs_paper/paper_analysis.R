################################################################################
## Script to do analysis of generated reservoir inputs for documentation      ##
## paper                                                                      ##
################################################################################

rm(list = ls())

library(RColorBrewer)

################################################################################
## Basic set up of generated reservoir inputs                                 ##
# Reservoir inputs based on different drainage direction maps
src_datasets <- list(
  "5arcmin" = c("hydrosheds_v1.1", "merit"),
  "30arcmin" = c("ddm30", "stn")
)
# Number of parameter combinations for deviation penalty, distance penalty and
# sign penalty
num_combi <- 4
# Use only combinations 1-3, because reservoir inputs for combination 4 are
# identical to combination 1
show_combi <- seq_len(3)
src_resol_names <- c("5arcmin" = "5'", "30arcmin" = "30'")
src_ds_names <- list(
  "5arcmin" = c("hydrosheds_v1.1" = "HydroSHEDS", "merit" = "MERIT"),
  "30arcmin" = c("ddm30" = "DDM30", "stn" = "STN-30")
)
grand_area_unit <- "km2"
################################################################################

## Load helper functions to work with LPJmL format                            ##
if (file.exists(file.path("..", "lpjml_format_helper_functions.R"))) {
  source(file.path("..", "lpjml_format_helper_functions.R"))
} else {
  stop("Please update path to script with LPJmL input format helper function")
}

################################################################################
## Load reservoir diagnostics and reservoir input files created               ##
reservoir_data <- list()
for (resol in names(src_datasets)) {
  reservoir_data[[resol]] <- list()
  for (ds in src_datasets[[resol]]) {
    tmp_env <- new.env()
    diagname <- paste0(
      "reservoir_diagnostics_",
      ifelse(nchar(ds) > 0, paste0(ds, "_"), ""),
      resol, "_", num_combi, "_settings.RData"
    )
    if (!file.exists(diagname)) {
      warning(diagname, " does not exist. Skipping dataset ", sQuote(ds),
              " at resolution ", sQuote(resol),
              call. = FALSE, immediate. = TRUE)
      next
    }
    cat("Load", sQuote(diagname), "\n")
    load(diagname, envir = tmp_env)
    lpjml_filename <- paste0(
      "reservoir_",
      ifelse(nchar(ds) > 0, paste0(ds, "_"), ""),
      resol,
      "_par", seq_len(num_combi),
      ".bin"
    )
    present <- file.exists(lpjml_filename)
    if (any(present)) {
      tmpheader <- read_header(lpjml_filename[min(which(present))])
      reservoir_data[[resol]][[ds]] <- array(
        dim = c(tmpheader$header[c("ncell", "nbands")], num_combi),
        dimnames = list(
          NULL,
          c(
            "year",
            "capacity",
            "area",
            "inst_cap",
            "height",
            paste0("purpose", seq_len(5))
          ),
          NULL
        )
      )
      for (i in seq_len(num_combi)) {
        if (file.exists(lpjml_filename[i])) {
          tmpheader <- read_header(lpjml_filename[i])
          if (tmpheader$header["nbands"] != 10) {
            stop("Unexpected number of bands ", tmpheader$header["nbands"],
                 " in file ", sQuote(lpjml_filename[min(which(present))]))
          }
          if (tmpheader$header["ncell"] != nrow(tmp_env$griddata)) {
            stop("Mismatch between number of cells in file ",
                 sQuote(lpjml_filename[i]), " and reservoir diagnostics data",
                 " from ", sQuote(diagname))
          }
          cat("Load", sQuote(lpjml_filename[i]), "\n")
          zz <- file(lpjml_filename[i], "rb")
          seek(zz, get_headersize(tmpheader))
          for (cell in seq_len(tmpheader$header["ncell"])) {
            # Reservoir data hads mixed data type. Confirm with
            # process_reservoir.R
            reservoir_data[[resol]][[ds]][cell, 1, i] <- readBin(
              zz, what = integer(), size = 4, n = 1, endian = tmpheader$endian
            )
            reservoir_data[[resol]][[ds]][cell, c(2, 3), i] <- readBin(
              zz, what = numeric(), size = 4, n = 2, endian = tmpheader$endian
            )
            reservoir_data[[resol]][[ds]][cell, seq(4, 10), i] <- readBin(
              zz, what = integer(), size = 4, n = 7, endian = tmpheader$endian
            )
          }
          close(zz)
        } else {
          warning("Reservoir input file ", sQuote(lpjml_filename[i]),
                  " does not exist.",
                  call. = FALSE, immediate. = TRUE)
        }
      }
      env_name <- paste0("res_diag_", ds, "_", resol)
      assign(env_name, tmp_env)
    } else {
      warning("No reservoir input corresponding to dataset ", sQuote(ds),
              " at resolution ", sQuote(resol), " found.",
              call. = FALSE, immediate. = TRUE)
    }
  }
}
################################################################################


################################################################################
## Calculate measure for paper                                                ##
diag_table <- data.frame(
  parameter = c(
    "Number of cells with reservoirs", rep("", length(show_combi)),
    "Median distance [m]", rep("", length(show_combi)),
    "\\% of dams not grid-based", rep("", length(show_combi)),
    "\\% of dams strong underestimation", rep("", length(show_combi)),
    "\\% of dams medium underestimation", rep("", length(show_combi)),
    "\\% of dams good match", rep("", length(show_combi)),
    "\\% of dams medium overestimation", rep("", length(show_combi)),
    "\\% of dams strong overestimation", rep("", length(show_combi))
  )
)
for (resol in names(reservoir_data)) {
  for (ds in names(reservoir_data[[resol]])) {
    env_name <- paste0("res_diag_", ds, "_", resol)
    tmp_env <- get(env_name)
    if (!"assignment_strategy" %in% colnames(diag_table)) {
      diag_table <- cbind(
        diag_table,
        assignment_strategy = c(
          "grid-based",
          paste0(
            "$p_{\\mathrm{dis}}=", tmp_env$parameter_setting[show_combi, "distance_penalty"],
            "$, $p_{\\mathrm{dev}}=",
            tmp_env$parameter_setting[show_combi, "deviation_penalty"],
            "$"
          )
        )
      )
    }
    column <- rep(NA, nrow(diag_table))
    # Number of cells with reservoirs
    id_cols <- c("CELL_ID_GRID", paste0("CELL_ID_CATCH", show_combi))
    res_cells <- apply(tmp_env$granddata[, id_cols], 2,
                       function(indata) length(unique(indata)))
    ind <- c(1, seq_along(show_combi) + 1)
    column[ind] <- as.character(res_cells)
    # Median distance
    ind <- ind + length(show_combi) + 1
    dis <- apply(
      tmp_env$diagnostics_table[, c("dist_grid", paste0("dist_best", show_combi))],
      2,
      median
    )
    column[ind] <- as.character(round(dis))
    # % of cells not located in grid-based cell
    ind <- ind + length(show_combi) + 1
    rel <- apply(
      tmp_env$granddata[, id_cols],
      2,
      function(cell, cell2) length(which(cell != cell2)),
      cell2 = tmp_env$granddata[, "CELL_ID_GRID"]
    ) / nrow(tmp_env$granddata) * 100
    column[ind] <- as.character(round(rel, 1))
    # strong underestimation
    dev_cols <- c("deviation_grid", paste0("deviation_best", show_combi))
    dev_cols2 <- c(
      "deviation_grid_cellsize",
      paste0("deviation_best", show_combi, "_cellsize")
    )
    ind <- ind + length(show_combi) + 1
    est <- apply(
      tmp_env$diagnostics_table[, dev_cols],
      2,
      function(indata) length(which(indata < -0.5))
    ) / nrow(tmp_env$diagnostics_table) * 100
    column[ind] <- as.character(round(est, 1))
    # medium underestimation
    ind <- ind + length(show_combi) + 1
    est <- apply(
      tmp_env$diagnostics_table[, dev_cols],
      2,
      function(indata) length(which(indata < -0.1 & indata >= -0.5))
    ) / nrow(tmp_env$diagnostics_table) * 100
    column[ind] <- as.character(round(est, 1))
    # good match
    ind <- ind + length(show_combi) + 1
    est <- apply(
      tmp_env$diagnostics_table[, dev_cols],
      2,
      function(indata) length(which(abs(indata) <= 0.1))
    ) / nrow(tmp_env$diagnostics_table) * 100
    est2 <- apply(
      tmp_env$diagnostics_table[, dev_cols2],
      2,
      function(indata) length(which(abs(indata) <= 0.1))
    ) / nrow(tmp_env$diagnostics_table) * 100
    column[ind] <- paste0(round(est, 1), " (", round(est2, 1), ")")
    # Medium overestimation
    ind <- ind + length(show_combi) + 1
    est <- apply(
      tmp_env$diagnostics_table[, dev_cols],
      2,
      function(indata) length(which(indata > 0.1 & indata <= 1))
    ) / nrow(tmp_env$diagnostics_table) * 100
    est2 <- apply(
      tmp_env$diagnostics_table[, dev_cols2],
      2,
      function(indata) length(which(indata > 0.1 & indata <= 1))
    ) / nrow(tmp_env$diagnostics_table) * 100
    column[ind] <- paste0(round(est, 1), " (", round(est2, 1), ")")
    # strong overestimation
    ind <- ind + length(show_combi) + 1
    est <- apply(
      tmp_env$diagnostics_table[, dev_cols],
      2,
      function(indata) length(which(indata > 1))
    ) / nrow(tmp_env$diagnostics_table) * 100
    est2 <- apply(
      tmp_env$diagnostics_table[, dev_cols2],
      2,
      function(indata) length(which(indata > 1))
    ) / nrow(tmp_env$diagnostics_table) * 100
    column[ind] <-paste0(round(est, 1), " (", round(est2, 1), ")")
    diag_table <- cbind(diag_table, column)
    colnames(diag_table)[ncol(diag_table)] <- paste0(ds, "_", resol)
  }
}
## Write measures to file in LaTeX format                                     ##
write.table(
  format(diag_table, justify = "right"),
  file = "paper_reservoir_diag.txt",
  sep = " & ", eol = "\\\\\n", quote = FALSE,
  col.names = FALSE, row.names = FALSE
)

width <- 6.5
height <- width / max(sapply(src_datasets, length)) * length(src_datasets)
p <- 1
pdf(
  "paper_catchment_areas.pdf",
  width = width, height = height, paper = "special",
  pointsize = 8
)
par(mfrow = c(length(src_datasets), max(sapply(src_datasets, length))),
    mar = c(3, 3, 1, 0) + 0.3)
for (resol in names(src_datasets)) {
  for (ds in src_datasets[[resol]]) {
    env_name <- paste0("res_diag_", ds, "_", resol)
    tmp_env <- get(env_name)
    cols <- grep("CATCH_", colnames(tmp_env$granddata), value = TRUE)
    min_val <- min(
      as.double(
        tmp_env$granddata[, "CATCH_SKM"][which(tmp_env$granddata[, "CATCH_SKM"] > 0)]
      ),
      as.matrix(tmp_env$granddata[, setdiff(cols, "CATCH_SKM")])
    )
    plot(x = pmax(min_val, tmp_env$granddata[, "CATCH_SKM"]),
         y = pmax(min_val, tmp_env$granddata[, "CATCH_AREA_LPJ_GRID"]),
         xlab = "", ylab = "", asp = 1, type = "p", pch = 0, mgp = c(2, 0.75, 0),
         xlim = c(min_val, max(tmp_env$granddata[, cols], na.rm = TRUE)),
         ylim = c(min_val, max(tmp_env$granddata[, cols], na.rm = TRUE)),
         cex = 0.8, lwd = 0.25, log = "xy"
    )
    for (i in seq_along(show_combi)) {
      points(x = pmax(min_val, tmp_env$granddata[, "CATCH_SKM"]),
             y = pmax(min_val, tmp_env$granddata[, paste0("CATCH_AREA_LPJ_BEST", i)]),
             pch = 1 + i, col = brewer.pal(num_combi, "Set1")[i],
             cex = 0.8 - i * 0.1, lwd = 0.25
      )
    }
    title(
      main = paste0(
        letters[p], ") ",
        src_resol_names[resol], " (", src_ds_names[[resol]][ds], ")"
      ),
      line = 0.2, adj = 0
    )
    title(
      xlab = switch(
        grand_area_unit,
        "ha" = "Catchment area (GRanD) [ha]",
        "km2" = expression(paste("Catchment area (GRanD) [", km^2, "]")),
        "m2" = expression(paste("Catchment area (GRanD) [", m^2, "]")),
        substitute(
          paste("Catchment area (GRanD) [", un, "]",
          list(un = grand_area_unit))
        )
      ),
      ylab = switch(
        grand_area_unit,
        "ha" = "DDM upstream area [ha]",
        "km2" = expression(paste("DDM upstream area [", km^2, "]")),
        "m2" = expression(paste("DDM upstream area [", m^2, "]")),
        substitute(
          paste("DDM upstream area [", un, "]",
          list(un = grand_area_unit))
        )
      ),
      mgp = c(2, 1, 0)
    )
    abline(0, 1, lwd = 0.25, lty = 2)
    p <- p + 1
  }
  for (i in seq_len(
    length(src_datasets[[resol]]) - max(sapply(src_datasets, length))
  )) {
    plot.new()
  }
}
legend(
  "bottomright", 
  legend = c(
    "grid-based",
    expression(paste(p[dis], " = ", 1, ", ", p[dev], " = ", 1)),
    expression(paste(p[dis], " = ", 1.5, ", ", p[dev], " = ", 1)),
    expression(paste(p[dis], " = ", 1, ", ", p[dev], " = ", 1.5))
  ),
  pch = c(0, 1 + show_combi),
  pt.cex = c(0.8, 0.8 - 0.1 * show_combi), pt.lwd = 0.25,
  col = c("black", brewer.pal(num_combi, "Set1")[show_combi]),
)
dev.off()
