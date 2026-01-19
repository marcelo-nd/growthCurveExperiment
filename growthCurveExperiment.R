# Install and load packages
if (!require("readxl", quietly = TRUE))
  install.packages("readxl")

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")

if (!require("growthcurver", quietly = TRUE))
  install.packages("growthcurver")

library(readxl)
library(ggplot2)
library(dplyr)
library(growthcurver)

# Class "GrowthCurve"
GrowthCurve <- setRefClass("GrowthCurve",
                           fields = list(name = "character", data = "data.frame",
                                         means = "data.frame", sds = "data.frame"),
                           methods = list(
                             initialize = function(name = character(), data = data.frame(),
                                                   strain_plate_cols = list(), #contains the plate columns used (generally species)
                                                   strain_plate_rows = list()) # contains the plate rows used (generally replicates)
                               { 
                               .self$name <- name # contains the name of the sample
                               .self$data <- as.data.frame(data) # contains a DF with the growth curve data (rows: time, columns: wells)
                               
                               # Function to calculate stats (means and sd for each strain (columns))
                               # This data is stored in means and sds DF attributes in each GrowthCurve object
                               .self$calculate_stats()
                               },
                             calculate_stats = function(){
                               rowSDs <- function(x, na.rm=F) {
                                 # Vectorised version of variance filter
                                 return(apply(x, MARGIN = 1, FUN = sd))
                               }
                               .self$means <- dplyr::select(.self$data, "Time") # gets the "Time" column for the means DF
                               .self$sds <- dplyr::select(.self$data, "Time") # gets the "Time" column for the sds DF
                               .self$means <- cbind(.self$means, data.frame(rowMeans(.self$data[2:ncol(.self$data)]))) # calculates means
                               .self$sds <- cbind(.self$sds, data.frame(rowSDs(.self$data[2:ncol(.self$data)]))) # calculates sds
                               colnames(.self$means) <- c("Time", "mean") # renames columns in means DF
                               colnames(.self$sds) <- c("Time", "sd") # renames columns in sds DF
                             }
                           )
)

# Class "GrowthCurveExperiment"
GrowthCurveExperiment <- setRefClass("GrowthCurveExperiment",
                            fields = list(name = "character", strains_names = "vector",
                                          growthCurveObjects = "list",
                                          od_means = "data.frame",
                                          od_sds = "data.frame",
                                          gc_d2f = "data.frame"),
                            methods = list(
                              # Create GrowthCurve object
                              initialize = function(name = character())
                                {
                                .self$name <- name # the only required value to init. a GreowthCurveExperiment object is the NAME
                              },
                              # The most common way to populate a GrowthCurveExperiment object is from a table.
                              # In this kind of experiment rows are replicates and columns are species or samples.
                              # If blank is used it should be in one column as a sample
                              create_gc_objects_from_table = function(gc_df_path,# path to excel file
                                                                      gc_range = character(),     # range inside the excel file
                                                                      strains_names = c(),        # sample/strain names (one per plate column group)
                                                                      strains_plate_cols = list(),# plate columns used (e.g. "1".."12")
                                                                      strains_plate_rows = list(),# plate rows used (e.g. "A".."H")
                                                                      blank = logical(),          # subtract blank well mean? TRUE/FALSE
                                                                      blank_col = NULL,           # which plate column is the blank (e.g. 12)
                                                                      pr_correction = TRUE,       # keep for later (reader can use it)
                                                                      psheet = character(),       # excel sheet name
                                                                      time_format = c("excel_day_fraction")  # for now: only this format
                                                                      ) {
                                
                                # --- Store names ---
                                .self$strains_names <- strains_names
                                
                                # --- Defaults ---
                                if (length(psheet) == 0) psheet <- 1  # allow "1" or first sheet
                                time_format <- match.arg(time_format)
                                
                                # --- Read table using the new general reader (we'll build it next) ---
                                gc_df <- .self$read.gc.file(
                                  gc_path       = gc_df_path,
                                  gc_range      = gc_range,
                                  psheet2       = psheet,
                                  n_plate_cols  = length(strains_plate_cols),
                                  n_plate_rows  = length(strains_plate_rows),
                                  time_format   = time_format,
                                  pr_correction = pr_correction
                                )
                                
                                print(head(gc_df))
                                print(tail(gc_df))
                                
                                # --- Optional: subtract blank measurements (your existing method) ---
                                if (isTRUE(blank) && !is.null(blank_col)) {
                                  gc_df_blnk <- .self$substract_blank(
                                    gc_df,
                                    blank_col = blank_col,
                                    str_plate_rows = strains_plate_rows
                                  )
                                } else {
                                  gc_df_blnk <- gc_df
                                }

                                # --- Build GrowthCurve objects ---
                                .self$growthCurveObjects <- list()  # reset list each time

                                bacteria_count <- 1
                                for (colstr in strains_plate_cols) {
                                  strain_wells <- c()
                                  for (rowstr in strains_plate_rows) {
                                    strain_wells <- c(strain_wells, paste(rowstr, colstr, sep = ""))
                                  }

                                  .self$growthCurveObjects <- append(
                                    .self$growthCurveObjects,
                                    GrowthCurve$new(
                                      name = .self$strains_names[bacteria_count],
                                      data = dplyr::select(gc_df_blnk, all_of(c("Time", strain_wells)))
                                    )
                                  )

                                  bacteria_count <- bacteria_count + 1
                                }

                                # --- Calculate stats DF for each strain/sample ---
                                .self$get_strains_stats_df()
                              },
                              read.gc.file = function(gc_path, gc_range, psheet2,
                                                      n_plate_cols = NULL,
                                                      n_plate_rows = NULL,
                                                      time_format = c("excel_day_fraction"),
                                                      pr_correction = TRUE,
                                                      correction_fun = NULL,
                                                      correction_a = 3.17,
                                                      correction_b = -0.26,
                                                      ignore_cols = NULL,
                                                      verbose = FALSE
                                                      ) {
                                time_format <- match.arg(time_format)
                                
                                # ---- 1) Read excel range ----
                                rd_raw <- readxl::read_excel(
                                  path = gc_path,
                                  sheet = psheet2,
                                  range = gc_range,
                                  col_names = TRUE
                                )
                                
                                if (!is.null(ignore_cols)) {
                                  rd_raw <- dplyr::select(rd_raw, -dplyr::any_of(ignore_cols))
                                }
                                
                                # ---- 2) Identify Time column ----
                                nms <- names(rd_raw)
                                time_idx <- which(tolower(nms) == "time")
                                if (length(time_idx) == 0) {
                                  stop("Could not find a 'Time' column (case-insensitive). Check the header row in your selected range.")
                                }
                                if (length(time_idx) > 1) {
                                  warning("More than one 'Time' column found; using the first one.")
                                  time_idx <- time_idx[1]
                                }
                                time_col_name <- nms[time_idx]
                                
                                # ---- 3) Keep only well columns A1..H12 ----
                                well_idx <- which(grepl("^[A-H](?:[1-9]|1[0-2])$", nms, perl = TRUE))
                                if (length(well_idx) == 0) {
                                  stop("No well columns detected. Expected names like A1..H12. Check whether your header row contains those well IDs.")
                                }
                                
                                rd_gc_df <- dplyr::select(rd_raw, dplyr::all_of(c(time_col_name, nms[well_idx])))
                                names(rd_gc_df)[1] <- "Time"
                                
                                # ---- 4) Optional sanity check ----
                                if (!is.null(n_plate_cols) && !is.null(n_plate_rows)) {
                                  expected_wells <- n_plate_cols * n_plate_rows
                                  got_wells <- ncol(rd_gc_df) - 1
                                  if (got_wells != expected_wells) {
                                    warning(sprintf(
                                      "Expected %d well columns (%d cols x %d rows) but detected %d. Continuing anyway.",
                                      expected_wells, n_plate_cols, n_plate_rows, got_wells
                                    ))
                                  }
                                }
                                
                                # ---- 5) Coerce well columns to numeric ----
                                for (j in 2:ncol(rd_gc_df)) {
                                  if (!is.numeric(rd_gc_df[[j]])) {
                                    rd_gc_df[[j]] <- suppressWarnings(as.numeric(rd_gc_df[[j]]))
                                  }
                                }
                                
                                # ---- 6) Parse Time + UNWRAP across midnight (FIX) ----
                                if (time_format == "excel_day_fraction") {
                                  t <- rd_gc_df$Time
                                  
                                  # Helper: unwrap a vector of hours where time resets at midnight
                                  unwrap_hours <- function(hours_vec) {
                                    hours_vec <- as.numeric(hours_vec)
                                    if (length(hours_vec) <= 1 || all(is.na(hours_vec))) return(hours_vec)
                                    wrap <- c(FALSE, diff(hours_vec) < 0)   # TRUE when time goes backwards
                                    hours_vec + cumsum(wrap) * 24
                                  }
                                  
                                  # A) Character input (often "HH:MM:SS"; sometimes Excel shows "1:0:09:10" but readxl drops the day)
                                  if (is.character(t)) {
                                    parse_hms_like_to_hours <- function(x) {
                                      parts <- strsplit(x, ":", fixed = TRUE)[[1]]
                                      nums <- suppressWarnings(as.numeric(parts))
                                      if (any(is.na(nums))) return(NA_real_)
                                      
                                      # Most robust approach here: interpret the LAST 3 fields as HH:MM:SS
                                      # This handles both "23:59:10" and "1:0:09:10" (we treat it as 00:09:10,
                                      # then the unwrapping step adds +24 based on the previous timepoint).
                                      if (length(nums) >= 3) {
                                        tail3 <- tail(nums, 3)
                                        h <- tail3[1]; m <- tail3[2]; s <- tail3[3]
                                        return(h + (m/60) + (s/3600))
                                      } else {
                                        return(NA_real_)
                                      }
                                    }
                                    
                                    hours <- as.numeric(vapply(t, parse_hms_like_to_hours, numeric(1)))
                                    rd_gc_df$Time <- unwrap_hours(hours)
                                    
                                  } else if (inherits(t, c("POSIXct", "POSIXt"))) {
                                    lt <- as.POSIXlt(t)
                                    hours <- lt$hour + (lt$min / 60) + (lt$sec / 3600)
                                    rd_gc_df$Time <- unwrap_hours(hours)
                                    
                                  } else if (inherits(t, "Date")) {
                                    rd_gc_df$Time <- as.numeric(t - t[1]) * 24
                                    
                                  } else {
                                    # B) Numeric input: could be Excel fraction-of-day (0..1) OR already hours
                                    x <- suppressWarnings(as.numeric(t))
                                    
                                    if (!all(is.na(x))) {
                                      if (max(x, na.rm = TRUE) <= 1.5) {
                                        # day fractions -> hours within day, then unwrap across midnight
                                        hours <- x * 24
                                        rd_gc_df$Time <- unwrap_hours(hours)
                                      } else {
                                        # already hours; still unwrap in case it resets for some reason
                                        rd_gc_df$Time <- unwrap_hours(x)
                                      }
                                    } else {
                                      rd_gc_df$Time <- x
                                    }
                                  }
                                }
                                
                                # ---- 7) Apply correction (optional) ----
                                if (!isFALSE(pr_correction)) {
                                  fun <- NULL
                                  if (is.function(pr_correction)) {
                                    fun <- pr_correction
                                  } else if (is.function(correction_fun)) {
                                    fun <- correction_fun
                                  } else if (isTRUE(pr_correction)) {
                                    fun <- function(x) x * correction_a + correction_b
                                  }
                                  
                                  if (is.function(fun)) {
                                    rd_gc_df[, 2:ncol(rd_gc_df)] <- as.data.frame(lapply(rd_gc_df[, 2:ncol(rd_gc_df)], fun))
                                  }
                                }
                                
                                if (verbose) {
                                  print(utils::head(rd_gc_df))
                                  print(summary(rd_gc_df$Time))
                                }
                                
                                return(as.data.frame(rd_gc_df))
                              },
                              parse_time_in_hours_biotek = function(time_string){
                                # Parse Time from Biotek-style results (strings in the form: "00:09:10")
                                # Import as numeric reults in values from 0 to 1.0 where 1.0 is 24 hours
                                return(time_string * 24)
                              },
                              parse_time_in_hours_infinite = function(time_string){
                                # Parse Time from Infinite-style results (in seconds)
                                return((time_string/60)/60)
                              },
                              substract_blank = function(nob_gc_df, blank_col = 1, str_plate_rows){
                                blank_wells <- c()
                                # Iterate over rows provided to generate the wells vector
                                for (rowstr in str_plate_rows) {
                                  blank_wells <- c(blank_wells, paste(rowstr, blank_col, sep = ""))
                                }
                                # Select blank wells and calculate blank means
                                blank_means <- rowMeans(dplyr::select(nob_gc_df, all_of(c(blank_wells))))
                                
                                # Copy non blank DF to avoid modifying the original one
                                new_df <- nob_gc_df[,2:ncol(nob_gc_df)]
                                # Now iterate over dataframe to subtract blank value in each row
                                if (length(blank_means == nrow(nob_gc_df))) {
                                  for (gc_row in 1:nrow(new_df)) {
                                    new_df[gc_row, ] <- new_df[gc_row, ] - blank_means[gc_row]
                                  }
                                }else{
                                    print("Gc and blank means length are not equal")
                                }
                                # Add the Time column to the blanked DF and make sure names are correct
                                new_df <- cbind(nob_gc_df[, "Time"], new_df)
                                colnames(new_df) <- c("Time", colnames(new_df[, 2:ncol(new_df)]))
                                return(new_df)
                              },
                              read.gc.file.biotek = function(gc_path, gc_range, n_plate_cols, n_plate_rows, psheet2, btk_correction=TRUE, blanked = TRUE){
                                print("Started Biotek reading file")
                                # Calculate the amount of samples (number of columns, each well is a sample belonging to a strain and a replicate)
                                n_samples <- n_plate_cols * n_plate_rows
                                print(n_samples)
                                # Read excel file
                                if (blanked) {
                                  rd_gc_df <- readxl::read_excel(path = gc_path, sheet = psheet2, range = gc_range, col_names = TRUE, col_types = c("numeric", rep("numeric", n_samples)))
                                }else{
                                  rd_gc_df <- readxl::read_excel(path = gc_path, sheet = psheet2, range = gc_range, col_names = TRUE, col_types = c("numeric", "numeric", rep("numeric", n_samples)))
                                  #rd_gc_df <- readxl::read_excel(path = gc_path, sheet = psheet2, range = gc_range, col_names = TRUE, col_types = c("numeric", rep("numeric", n_samples)))
                                  
                                  rd_gc_df <- rd_gc_df[, -2]
                                }
                                
                                print(head(rd_gc_df))
                                # Parse time 
                                rd_gc_df$Time <- lapply(rd_gc_df$Time, .self$parse_time_in_hours_biotek)
                                rd_gc_df$Time <- as.numeric(rd_gc_df$Time)
                                # Adjust with Plate Readers conversion factor (od = preader.od * 4.4131 - 0.3588)
                                if (btk_correction == TRUE) {
                                  rd_gc_df <- cbind(rd_gc_df$Time, rd_gc_df[, 2:ncol(rd_gc_df)]*3.17 - 0.26)
                                }
                                colnames(rd_gc_df) <- c("Time", colnames(rd_gc_df)[2:length(rd_gc_df)])
                                print("Finished reading Biotek file")
                                return(rd_gc_df)
                              },
                              read.gc.file.spark = function(gc_path, gc_range, n_plate_cols, n_plate_rows, spk_correction=TRUE, psheet2){
                                # Calculate the amount of samples 
                                n_samples <- n_plate_cols * n_plate_rows
                                # Read excel file
                                .self$gc_df <- readxl::read_excel(path = gc_path, sheet = psheet2, range = gc_range, col_names = TRUE, col_types = c("numeric", rep("numeric", n_samples)))
                                # Parse time 
                                .self$gc_df$Time <- lapply(gc_df$Time, .self$parse_time_in_hours_biotek)
                                .self$gc_df$Time <- as.numeric(.self$gc_df$Time)
                                # Adjust with Plate Readers conversion factor (od = preader.od * 4.4131 - 0.3588)
                                if (spk_correction) {
                                  
                                }
                                .self$gc_df <- cbind(.self$gc_df$Time, .self$gc_df[, 2:ncol(.self$gc_df)]*3.18 - 0.28)
                                colnames(.self$gc_df) <- c("Time", colnames(.self$gc_df)[2:length(.self$gc_df)])
                                return(gc_df)
                              },
                              read.gc.file.infinite = function(gc_path, gc_range, n_plate_cols, n_plate_rows, inf_correction=TRUE, psheet2){
                                # Calculate the amount of samples 
                                n_samples <- n_plate_cols * n_plate_rows
                                # Read excel file
                                rd_gc_df <- readxl::read_excel(path = gc_path, sheet = psheet2, range = gc_range, col_names = TRUE, col_types = c("numeric", rep("numeric", n_samples+1)))
                                # Get rid of column with temperature values in Infinite PR result files
                                rd_gc_df <- subset(rd_gc_df, select = -c(`Temp. [°C]`))
                                colnames(rd_gc_df) <- c("Time", colnames(rd_gc_df)[2:length(rd_gc_df)])
                                # Parse time 
                                rd_gc_df$Time <- lapply(rd_gc_df$Time, .self$parse_time_in_hours_infinite)
                                rd_gc_df$Time <- as.numeric(rd_gc_df$Time)
                                # Adjust with Plate Readers conversion factor (od = preader.od * 4.4131 - 0.3588)
                                if (inf_correction) {
                                  rd_gc_df <- cbind(rd_gc_df$Time, rd_gc_df[, 2:ncol(rd_gc_df)]*4.4131 - 0.3588)
                                  colnames(rd_gc_df) <- c("Time", colnames(rd_gc_df)[2:length(rd_gc_df)])
                                }
                                return(rd_gc_df)
                              },
                              calculate_growth_curve_models = function(growth.curve.points){
                                # Based on https://rpubs.com/angelov/growthcurver.
                                summG <- function(x) {growthcurver::SummarizeGrowth(growth.curve.points$Time,x)}
                                # Get models
                                models.all <- lapply(growth.curve.points[2:ncol(growth.curve.points)], summG)
                                # Create predicted values DF
                                df.predicted <- data.frame(time = growth.curve.points$Time)
                                # Iterate over all samples/species
                                for (i in names(growth.curve.points[2:ncol(growth.curve.points)])) 
                                {
                                  df.predicted_val <- NA
                                  # Use tryCatch to avoid halting scripts execution when models cannot be calculated
                                  tryCatch(
                                    {
                                      # Just to highlight: if you want to use more than one
                                      # R expression in the "try" part then you'll have to
                                      # use curly brackets.
                                      # 'tryCatch()' will return the last evaluated expression
                                      # in case the "try" part was completed successfully
                                      #message("This is the 'try' part")
                                      
                                      #suppressWarnings(df.predicted[[i]] <- predict(models.all[[i]]$model))
                                      df.predicted_val <- predict(models.all[[i]]$model)
                                      # The return value of `readLines()` is the actual value
                                      # that will be returned in case there is no condition
                                      # (e.g. warning or error).
                                    },
                                    error = function(cond) {
                                      #message(paste("URL does not seem to exist:", url))
                                      message("Here's the original error message:")
                                      message(conditionMessage(cond))
                                      # Choose a return value in case of error
                                      df.predicted_val <- NA
                                    },
                                    finally = {
                                      # NOTE:
                                      # Here goes everything that should be executed at the end,
                                      # regardless of success or error.
                                      # If you want more than one expression to be executed, then you
                                      # need to wrap them in curly brackets ({...}); otherwise you could
                                      # just have written 'finally = <expression>' 
                                      #message("Some other message at the end")
                                      df.predicted[[i]] <- df.predicted_val
                                    }
                                  )
                                }
                                # Reassign names
                                colnames(df.predicted) <- c("Time", colnames(df.predicted)[2:ncol(df.predicted)])
                                # Return de DF with predicted values for all species/strains
                                return(df.predicted)
                              },
                              get_strains_stats_df = function(){
                                print("Getting stats for all GCO")
                                # Empty dfs for means and sds, first column is "Time" variable.
                                .self$od_means <- dplyr::select(.self$growthCurveObjects[[1]]$data, "Time")
                                .self$od_sds <- dplyr::select(.self$growthCurveObjects[[1]]$data, "Time")
                                for (gco in growthCurveObjects) {
                                  # Calculate the stats for each strain in GrowthCurve objects list
                                  gco$calculate_stats()
                                  # and append to each stats df
                                  gco_means <- dplyr::select(gco$means, "mean")
                                  colnames(gco_means) <- gco$name
                                  od_means <<- cbind(od_means, gco_means)
                                  
                                  gco_sds <- dplyr::select(gco$sds, "sd")
                                  colnames(gco_sds) <- gco$name
                                  od_sds <<- cbind(od_sds, gco_sds)
                                }
                              },
                              add_gco = function(gco_to_add){
                                # Add strain/sample name to names list
                                .self$strains_names <- c(.self$strains_names, gco_to_add[[1]]$name)
                                # Add GrowthCurve object to list
                                .self$growthCurveObjects <- c(.self$growthCurveObjects, gco_to_add)
                                # Recalculate strains/samples stats
                                get_strains_stats_df()
                              },
                              remove_gco = function(gco_to_remove_name){
                                # Iterate over all GrowthCurves object in list 
                                for (gco in 1:length(.self$growthCurveObjects)) {
                                  # If the name in gco from list corresponds to the one to be removed
                                  if (.self$growthCurveObjects[[gco]]$name == gco_to_remove_name) {
                                    # When the Growth curve object in list with the name to remove is found, delete it
                                    .self$growthCurveObjects <- .self$growthCurveObjects[.self$strains_names != gco_to_remove_name]
                                    break # exit loop
                                  }
                                }
                                # Also remove strain/species name from list
                                .self$strains_names <- .self$strains_names[.self$strains_names != gco_to_remove_name]
                                # Recalculate strains/samples stats
                                get_strains_stats_df()
                              },
                              merge_experiments = function(gcexperiment){},
                              plot_curves = function(yScalemin = NULL, yScalemax = NULL, calculate_model = TRUE){
                                # 
                                od_g <- tidyr::gather(.self$od_means, key = "Species", value = "OD", 2:(length(.self$strains_names) + 1))
                                # 
                                od_sds_g <- tidyr::gather(.self$od_sds, key = "Species", value = "SD", 2:(length(.self$strains_names) + 1))
                                # 
                                od_means_sds_preds <- cbind(od_g, od_sds_g$SD)
                                # 
                                colnames(od_means_sds_preds) <- c("Time", "Species", "od", "sd")
                                
                                if (calculate_model) {
                                  #
                                  df.predicted.plate <- calculate_growth_curve_models(.self$od_means)
                                  
                                  # 
                                  
                                  pred_g <- tidyr::gather(df.predicted.plate, key = "Species", value = "pred.od", 2:(length(.self$strains_names) + 1))
                                  
                                  od_means_sds_preds <- cbind(od_g, od_sds_g$SD, pred_g$pred.od)
                                  
                                  colnames(od_means_sds_preds) <- c("Time", "Species", "od", "sd", "pred.od")
                                }
                                
                                color_scale = c("#AD7BE9", "#3E54AC", "#658864", "#6D9886", "#E96479", "#E69F00", "#FC7300",  "#183A1D", "#635985", "#F99417", "#FEA1BF", "#5BC0F8", "#DC0000", "#495579")
                                
                                if (!is.null(yScalemin) && !is.null(yScalemax)) {
                                  curves_plot <- od_means_sds_preds %>%
                                    ggplot(aes(x = Time, y = od, group = Species, color = Species)) +
                                    geom_errorbar(aes(ymin = od - sd, ymax = od + sd), width= 0.1) +
                                    geom_point(alpha=0.7) +
                                    scale_color_manual(values=color_scale) +
                                    ylim(yScalemin, yScalemax) +
                                    labs(y= "OD (600nm)", x = "Time (h)") 
                                }else{
                                  curves_plot <- od_means_sds_preds %>%
                                    ggplot(aes(x = Time, y = od, group = Species, color = Species)) +
                                    geom_errorbar(aes(ymin = od - sd, ymax = od + sd), width= 0.1) +
                                    geom_point(alpha=0.7) +
                                    scale_color_manual(values=color_scale) +
                                    labs(y= "OD (600nm)", x = "Time (h)")
                                }
                                
                                if (calculate_model) {
                                  # Add model lines
                                  curves_plot <- curves_plot +
                                    geom_line(aes(x = Time, y=pred.od), linewidth = 0.5)
                                }
                                return(curves_plot)
                              },
                              remove_replicate_in_gco = function(gco_name = NULL, reps_2_remove = NULL){
                                print(reps_2_remove)
                                if (is.null(gco_name) | is.null(reps_2_remove)) {
                                  print("Missing arguments")
                                  break
                                }else{
                                  new_gco_data <- .self$growthCurveObjects[.self$strains_names == gco_name][[1]]$data
                                  #new_gco_data <- new_gco_data[, -(rep_2_remove + 1)]
                                  new_gco_data <- new_gco_data[, !(names(new_gco_data) %in% reps_2_remove)]
                                  .self$growthCurveObjects[.self$strains_names == gco_name][[1]]$data <- new_gco_data
                                  
                                  #.self$growthCurveObjects[.self$strains_names == gco_name][[1]]$data
                                  
                                  get_strains_stats_df()
                                }
                              }
                            )
)

# To do: merge experiments, packages calls, multipanel graphs.
