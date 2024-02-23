library(readxl)

# Class "GrowthCurve"
GrowthCurve <- setRefClass("GrowthCurve",
                           fields = list(name = "character", data = "data.frame",
                                         data_stats = "list", means = "data.frame", sds = "data.frame"),
                           methods = list(
                             initialize = function(name = character(), data = data.frame(),
                                                   strain_plate_cols = list(),
                                                   strain_plate_rows = list()){
                               .self$name <- name
                               .self$data <- as.data.frame(data)
                               #.self$data$Time <- lapply(.self$data$Time, .self$parse_time_in_hours)
                               
                               # Calculate stats
                               #.self$data_stats <- .self$calculate_stats()
                               #.self$means <- .self$data_stats[[1]]
                               #.self$sds <- .self$data_stats[[2]]
                               .self$calculate_stats()
                               },
                             calculate_stats = function(){
                               rowSDs <- function(x, na.rm=F) {
                                 # Vectorised version of variance filter
                                 return(apply(x, MARGIN = 1, FUN = sd))
                               }
                               .self$means <- dplyr::select(.self$data, "Time")
                               .self$sds <- dplyr::select(.self$data, "Time")
                               .self$means <- cbind(.self$means, data.frame(rowMeans(.self$data[2:ncol(.self$data)])))
                               .self$sds <- cbind(.self$sds, data.frame(rowSDs(.self$data[2:ncol(.self$data)])))
                               colnames(.self$means) <- c("Time", "mean")
                               colnames(.self$sds) <- c("Time", "sd")
                               #return(list(.self$means, .self$sds))
                             }
                           )
)

# Class "GrowthCurveExperiment"
GrowthCurveExperiment <- setRefClass("GrowthCurveExperiment",
                            fields = list(name = "character", strains_names = "vector",
                                          growthCurveObjects = "list",
                                          od_means = "data.frame",
                                          od_sds = "data.frame",
                                          gc_df = "data.frame"),
                            methods = list(
                              # Create GrowthCurve object
                              initialize = function(name = character())
                                {
                                .self$name <- name
                              },
                              parse_time_in_hours_biotek = function(time_string){
                                # Parse Time from Biotek-style results (strings in the form: "00:09:10")
                                return(time_string * 24)
                              },
                              read.gc.file.biotek = function(gc_path, gc_range, n_plate_cols, n_plate_rows)
                                {
                                # Calculate the amount of samples 
                                n_samples <- n_plate_cols * n_plate_rows
                                # Read excel file
                                .self$gc_df <- readxl::read_excel(path = gc_path, sheet = "Plate 1 - Sheet1", range = gc_range, col_names = TRUE, col_types = c("numeric", rep("numeric", n_samples)))
                                # Parse time 
                                .self$gc_df$Time <- lapply(gc_df$Time, .self$parse_time_in_hours_biotek)
                                .self$gc_df$Time <- as.numeric(.self$gc_df$Time)
                                return(gc_df)
                              },
                              create_gc_objects_from_table1 = function(gc_df_path,
                                                                       gc_range = character(),
                                                                       plate_reader_type = character(),
                                                                       strains_names = c(),
                                                                       strains_plate_cols = list(),
                                                                       strains_plate_rows = list(),
                                                                       blank = logical()
                                                                       ){
                                .self$strains_names = strains_names
                                
                                if(plate_reader_type == "Biotek") {
                                  .self$gc_df <- .self$read.gc.file.biotek(gc_path = gc_df_path, gc_range = gc_range,
                                                                     n_plate_cols = length(strains_plate_cols),
                                                                     n_plate_rows = length(strains_plate_rows))
                                }
                                
                                # Create list of GrowthCurve Objects.
                                # Get wells list for each strain
                                bacteria_count = 1
                                for (colstr in strains_plate_cols) {
                                  #strain_wells <- list()
                                  strain_wells <- c()
                                  for (rowstr in strains_plate_rows) {
                                    #print(paste(rowstr, colstr, sep = ""))
                                    #strain_wells <- append(strain_wells, paste(rowstr, colstr, sep = ""))
                                    strain_wells <- c(strain_wells, paste(rowstr, colstr, sep = ""))
                                  }
                                  # Generate and add the growth curve object to list
                                  #print(strains_names[[bacteria_count]])
                                  #print(strain_wells)
                                  #print(c("Time", strain_wells))
                                  #print(head(.self$data))
                                  growthCurveObjects <<- append(growthCurveObjects, GrowthCurve$new(name=.self$strains_names[bacteria_count], data  = dplyr::select(.self$gc_df, all_of(c("Time", strain_wells)))))
                                  bacteria_count = bacteria_count + 1
                                }
                                get_strains_stats_df()
                              },
                              get_strains_stats_df = function(){
                                # Empty dfs for means and sds, first column is "Time" variable.
                                .self$od_means <- dplyr::select(.self$gc_df, "Time")
                                .self$od_sds <- dplyr::select(.self$gc_df, "Time")
                                for (gco in growthCurveObjects) {
                                  #print(gco$name)
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
                                .self$strains_names <- c(.self$strains_names, gco_to_add[[1]]$name)
                                .self$growthCurveObjects[[1]] <- c(.self$growthCurveObjects[[1]], gco_to_add)
                                get_strains_stats_df()
                              },
                              remove_gco = function(gco_to_remove_name){
                                for (gco in 1:length(.self$growthCurveObjects)) {
                                  #print(gco)
                                  if (.self$growthCurveObjects[[gco]]$name == gco_to_remove_name) {
                                    #print(.self$growthCurveObjects[[gco]]$name != gco_to_remove_name)
                                    .self$growthCurveObjects <- .self$growthCurveObjects[.self$strains_names != gco_to_remove_name]
                                    break
                                  }
                                }
                                .self$strains_names <- .self$strains_names[.self$strains_names != gco_to_remove_name]
                                
                                get_strains_stats_df()
                              },
                              merge_experiments = function(gcexperiment){},
                              plot_curves = function(yScalemin = NULL, yScalemax = NULL, calculate_model = TRUE){
                                #print(length(.self$strains_names[[1]]))
                                od_g <- tidyr::gather(.self$od_means, key = "Species", value = "OD", 2:(length(.self$strains_names) + 1))
                                #print(head(od_g))
                                od_sds_g <- tidyr::gather(.self$od_sds, key = "Species", value = "SD", 2:(length(.self$strains_names) + 1))
                                #print(head(od_sds_g))
                                od_means_sds_preds <- cbind(od_g, od_sds_g$SD)
                                #print(head(od_means_sds_preds))
                                #print(summary(od_means_sds_preds))
                                #print(sapply(od_means_sds_preds, mode))
                                colnames(od_means_sds_preds) <- c("Time", "Species", "od", "sd")
                                
                                if (calculate_model) {
                                  #
                                  df.predicted.plate <- calculate_growth_curve_models(.self$od_means)
                                  
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
                                
                                if (calculate_model == TRUE) {
                                  # Add model lines
                                  curves_plot +
                                    geom_line(aes(x = Time, y=pred.od), linewidth = 0.5)
                                }
                              }
                            )
)
