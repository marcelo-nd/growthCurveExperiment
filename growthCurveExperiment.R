# Class "GrowthCurve"
GrowthCurve <- setRefClass("GrowthCurve",
                           fields = list(name = "character", data = "data.frame",
                                         data_stats = "list", means = "data.frame", sds = "data.frame"),
                           
                           methods = list(
                             initialize = function(name = character(), data = data.frame(),
                                                   strain_plate_cols = list(),strain_plate_rows = list()){
                               .self$name <- name
                               .self$data <- as.data.frame(data)
                               #.self$data$Time <- lapply(.self$data$Time, .self$parse_time_in_hours)
                               
                               # Calculate stats
                               .self$data_stats <- .self$calculate_stats()
                               .self$means <- .self$data_stats[[1]]
                               .self$sds <- .self$data_stats[[2]]
                               },
                             parse_time_in_hours = function(time_string){
                               time_split <- strsplit(time_string, split = ":")
                               hours <- as.numeric(time_split[[1]][1])
                               mins <- as.numeric(time_split[[1]][2]) + (as.numeric(time_split[[1]][3])/60)
                               return(hours + (mins/60))
                               },
                             parse_time_column = function(){
                               .self$data$Time <- lapply(.self$data$Time, .self$parse_time_in_hours)
                               },
                             calculate_stats = function(){
                               rowSDs <- function(x, na.rm=F) {
                                 # Vectorised version of variance filter
                                 return(apply(x, MARGIN = 1, FUN = sd))
                               }
                               
                               od_means <- dplyr::select(.self$data, "Time")
                               od_sds <- dplyr::select(.self$data, "Time")
                               od_means <- cbind(od_means, data.frame(rowMeans(.self$data[2:ncol(.self$data)])))
                               od_sds <- cbind(od_sds, data.frame(rowSDs(.self$data[2:ncol(.self$data)])))
                               colnames(od_means) <- c("Time", "mean")
                               colnames(od_sds) <- c("Time", "sd")
                               return(list(od_means, od_sds))
                             }
                           )
)

# Class "GrowthCurveExperiment"
GrowthCurveExperiment <- setRefClass("GrowthCurveExperiment",
                            fields = list(name = "character", data = "data.frame",
                                          strains_names = "list",
                                          strains_plate_cols = "list",
                                          strains_plate_rows = "list", 
                                          replicates = "numeric",
                                          blank = "logical",
                                          growthCurveObjects = "list",
                                          od_means = "data.frame",
                                          od_sds = "data.frame",
                                          parseTime = "logical"),
                            methods = list(
                              initialize = function(name = character()
                              ){
                                .self$name <- name
                              },
                              create_gc_objects_from_table1 = function(data = data.frame(),
                                                                       strains_names = list(),
                                                                       strains_plate_cols = list(),
                                                                       strains_plate_rows = list(),
                                                                       blank = logical(),
                                                                       parseTime = logical()){
                                # Set all fields with the values passed in the constructor
                                .self$strains_plate_rows <- strains_plate_rows
                                .self$data <- as.data.frame(data)
                                .self$od_means <- data.frame()
                                .self$od_sds <- data.frame()
                                .self$strains_names <- as.list(strains_names)
                                .self$parseTime <- parseTime
                                
                                # Get number of replicates
                                .self$replicates = length(.self$strains_plate_rows)
                                
                                # Parse time column in the data DF.
                                if (.self$parseTime) {
                                  .self$parse_time_column()
                                }
                                #print(head(.self$data))
                                
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
                                  # generate and add the growth curve object to list
                                  #print(strains_names[[bacteria_count]])
                                  growthCurveObjects <<- append(growthCurveObjects, GrowthCurve$new(name=strains_names[[1]][bacteria_count], data  = select(.self$data, all_of(c("Time", strain_wells)))))
                                  bacteria_count = bacteria_count + 1
                                }
                              },
                              parse_time_in_hours = function(time_string){
                                time_split <- strsplit(time_string, split = ":")
                                hours <- as.numeric(time_split[[1]][1])
                                mins <- as.numeric(time_split[[1]][2]) + (as.numeric(time_split[[1]][3])/60)
                                return(hours + (mins/60))
                              },
                              parse_time_column = function(){
                                .self$data$Time <- as.double(lapply(.self$data$Time, .self$parse_time_in_hours))
                              },
                              get_strains_stats_df = function(){
                                # Empty dfs for means and sds, first column is "Time" variable.
                                .self$od_means <- dplyr::select(.self$data, "Time")
                                .self$od_sds <- dplyr::select(.self$data, "Time")
                                for (gco in growthCurveObjects) {
                                  print(gco$name)
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
                                strains_names <<- append(strains_names, gco_to_add$name)
                                
                                gco_means <- dplyr::select(gco_to_add$means, "mean")
                                colnames(gco_means) <- gco_to_add$name
                                od_means <<- cbind(od_means, gco_means)
                                
                                gco_sds <- dplyr::select(gco_to_add$sds, "sd")
                                colnames(gco_sds) <- gco_to_add$name
                                od_sds <<- cbind(od_sds, gco_sds)
                              },
                              remove_gco = function(gco_to_remove){
                                od_means <<- dplyr::select(.self$od_means, -any_of(c(gco_to_remove)))
                                od_sds <<- dplyr::select(.self$od_sds, -any_of(c(gco_to_remove)))
                                .self$strains_names[[1]] <- .self$strains_names[[1]][.self$strains_names[[1]] != gco_to_remove]
                              },
                              merge_experiments = function(gco_to_remove){
                                #toDo
                              },
                              plot_curves = function(){
                                #print(length(.self$strains_names[[1]]))
                                od_g <- tidyr::gather(.self$od_means, key = "Species", value = "OD", 2:(length(.self$strains_names[[1]]) + 1))
                                print(od_g, 20)
                                od_sds_g <- tidyr::gather(.self$od_sds, key = "Species", value = "SD", 2:(length(.self$strains_names[[1]]) + 1))
                                
                                od_means_sds_preds <- cbind(od_g, od_sds_g$SD)
                                
                                colnames(od_means_sds_preds) <- c("Time", "Species", "od", "sd")
                                
                                color_scale = c("#AD7BE9", "#3E54AC", "#658864", "#6D9886", "#E96479", "#E69F00", "#FC7300",  "#183A1D", "#635985", "#F99417", "#FEA1BF", "#5BC0F8", "#DC0000", "#495579")
                                
                                curves_plot <- od_means_sds_preds %>%
                                  ggplot(aes(x = Time, y = od, group = Species, color = Species)) +
                                  geom_errorbar(aes(ymin = od - sd, ymax = od + sd), width= 0.1) +
                                  geom_point(alpha=0.7) +
                                  scale_color_manual(values=color_scale)
                                curves_plot
                              }
                            )
)
