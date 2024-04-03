# Install and load packages
if (!require("readxl", quietly = TRUE))
  install.packages("readxl")

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")

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
                              create_gc_objects_from_table = function(gc_df_path, # path to table
                                                                      gc_range = character(), # The range inside the excel file
                                                                      plate_reader_type = character(), # One of the 3 types of plate readers
                                                                      strains_names = c(), # The strains or samples names 
                                                                      strains_plate_cols = list(), # the columns measured in the plate
                                                                      strains_plate_rows = list(), # The rows measured in the plate
                                                                      blank = logical(), # if blank has to be substracted, TRUE
                                                                      blank_col = NULL, # What is the column of the blank
                                                                      pr_correction = TRUE, # should plate reader correction be applied use TRUE
                                                                      psheet = "Plate 1 - Sheet1"
                              ){
                                .self$strains_names = strains_names
                                # Decide which plate reader type generated the table file and use respective function
                                if(plate_reader_type == "Biotek") {
                                  gc_df <- .self$read.gc.file.biotek(gc_path = gc_df_path, gc_range = gc_range,
                                                                           n_plate_cols = length(strains_plate_cols),
                                                                           n_plate_rows = length(strains_plate_rows),
                                                                           btk_correction = pr_correction,
                                                                           psheet = psheet)
                                }else if (plate_reader_type == "Spark"){
                                  gc_df <- .self$read.gc.file.spark(gc_path = gc_df_path, gc_range = gc_range,
                                                                          n_plate_cols = length(strains_plate_cols),
                                                                          n_plate_rows = length(strains_plate_rows),
                                                                          spk_correction = pr_correction)
                                }else if (plate_reader_type == "Infinite"){
                                  gc_df <- .self$read.gc.file.infinite(gc_path = gc_df_path, gc_range = gc_range,
                                                                             n_plate_cols = length(strains_plate_cols),
                                                                             n_plate_rows = length(strains_plate_rows),
                                                                             inf_correction = pr_correction)
                                }else{
                                  print("Plate reader type does not exist")
                                }
                                
                                # Subtract blank measurements from blank column
                                if (blank == TRUE && !is.null(blank_col)) {
                                  # call the substract_blank function contained in the GrowthCUrveExperiment Class
                                  gc_df_blnk <- substract_blank(gc_df, blank_col = blank_col, str_plate_rows = strains_plate_rows)
                                }else{
                                  # If no blank is going to be substracted just change the name pf the DF
                                  gc_df_blnk <- gc_df
                                }
                                
                                # Now we create the GrowthCurve Objects from the data in the table file
                                # The GrowthCurve objects are stored in a list in the GrowthCurveExperiment object
                                # First we generate the wells list for each strain/sample
                                bacteria_count = 1 # counter to keep track of the strain/sample number
                                # Each column is a strain.
                                for (colstr in strains_plate_cols) {
                                  strain_wells <- c()
                                  for (rowstr in strains_plate_rows) {
                                   # Each strain has nrow number of replicates
                                    strain_wells <- c(strain_wells, paste(rowstr, colstr, sep = ""))
                                  }
                                  # Generate and add the growth curve object to list
                                  # The generated well names are used to select columns from GC DF for each strain (column)
                                  # The name of the GrowthCurve object is taken from the list of strain provided so order of this list is important
                                  # The selected data is used to generate a GrowthCurve object which is appended to the Growthcurve objects list
                                  growthCurveObjects <<- append(growthCurveObjects, GrowthCurve$new(name=.self$strains_names[bacteria_count], data  = dplyr::select(gc_df_blnk, all_of(c("Time", strain_wells)))))
                                  bacteria_count = bacteria_count + 1 # counter to keep track of the strain/sample number
                                }
                                # Calculate stats DF for each strain/sample. These DFs are stored in the GrowthCurve objects
                                get_strains_stats_df()
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
                              read.gc.file.biotek = function(gc_path, gc_range, n_plate_cols, n_plate_rows, btk_correction=TRUE, psheet = psheet){
                                # Calculate the amount of samples (number of columns, each well is a sample belonging to a strain and a replicate)
                                n_samples <- n_plate_cols * n_plate_rows
                                # Read excel file
                                rd_gc_df <- readxl::read_excel(path = gc_path, sheet = psheet, range = gc_range, col_names = TRUE, col_types = c("numeric", rep("numeric", n_samples)))
                                # Parse time 
                                rd_gc_df$Time <- lapply(rd_gc_df$Time, .self$parse_time_in_hours_biotek)
                                rd_gc_df$Time <- as.numeric(rd_gc_df$Time)
                                # Adjust with Plate Readers conversion factor (od = preader.od * 4.4131 - 0.3588)
                                if (btk_correction == TRUE) {
                                  rd_gc_df <- cbind(rd_gc_df$Time, rd_gc_df[, 2:ncol(rd_gc_df)]*3.17 - 0.26)
                                }
                                colnames(rd_gc_df) <- c("Time", colnames(rd_gc_df)[2:length(rd_gc_df)])
                                return(rd_gc_df)
                              },
                              read.gc.file.spark = function(gc_path, gc_range, n_plate_cols, n_plate_rows, spk_correction=TRUE){
                                # Calculate the amount of samples 
                                n_samples <- n_plate_cols * n_plate_rows
                                # Read excel file
                                .self$gc_df <- readxl::read_excel(path = gc_path, sheet = "Plate 1 - Sheet1", range = gc_range, col_names = TRUE, col_types = c("numeric", rep("numeric", n_samples)))
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
                              read.gc.file.infinite = function(gc_path, gc_range, n_plate_cols, n_plate_rows, inf_correction=TRUE, p_sheet = "Sheet2"){
                                # Calculate the amount of samples 
                                n_samples <- n_plate_cols * n_plate_rows
                                # Read excel file
                                rd_gc_df <- readxl::read_excel(path = gc_path, sheet = p_sheet, range = gc_range, col_names = TRUE, col_types = c("numeric", rep("numeric", n_samples+1)))
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
                                # Empty dfs for means and sds, first column is "Time" variable.
                                .self$od_means <- dplyr::select(.self$growthCurveObjects[[1]]$data, "Time")
                                .self$od_sds <- dplyr::select(.self$growthCurveObjects[[1]]$data, "Time")
                                for (gco in growthCurveObjects) {
                                  # Calculate the stats for each strain in GrowthCurve objects list
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
                              }
                            )
)
# To do: merge experiments, packages calls, multipanel graphs.
