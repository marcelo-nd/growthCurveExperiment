GrowthCurve <- setRefClass("GrowthCurve",
                           fields = list(name = "character", data = "data.frame"),
                           
                           methods = list(
                             initialize = function(name = character(), data = data.frame(),
                                                   strain_plate_cols = list(),strain_plate_rows = list()){
                               .self$name <- name
                               .self$data <- as.data.frame(data)
                               #.self$data$Time <- lapply(.self$data$Time, .self$parse_time_in_hours)
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


Growthcurve.group3 <- read.csv2("C:/Users/Marcelo/OneDrive - UT Cloud/Postdoc Tü/Sci/NoseSynComProject/6_Growth curves/Growthcurve group3 repeats.csv")

gc2 <- GrowthCurve$new(name="test_name", data  = Growthcurve.group3)

gc2$parse_time_column()

gc2$data

gc2$calculate_stats()

# create a class "GrowthCurveExperiment"
GrowthCurveExperiment <- setRefClass("GrowthCurveExperiment",
                            fields = list(name = "character", data = "data.frame",
                                          strains_names = "list",
                                          strains_plate_cols = "list",
                                          strains_plate_row_start = "list", 
                                          replicates = "numeric",
                                          blank = "logical",
                                          growthCurveObjects = "list"),
                            methods = list(
                              
                            )
                            )


# 
experiment1 <- GrowthCurveExperiment(name = "exp1", data = Growthcurve.group3,
                                     strains_names = list("R. mucilaginosa", "C. tuberculosteariscum",
                                                          "C. kroppenstedtii", "C. propinquum",
                                                          "C. accolens"),
                                     strains_plate_cols = list(4,5,6,7,8),
                                     strains_row_cols = list("A:H"),
                                     blank = FALSE
                                     )

# experiment1
experiment1



#strain_plate_cols = "list", strain_plate_rows = "list"
#.self$strain_plate_cols  <- as.list(strain_plate_cols)
#.self$strain_plate_rows  <- as.list(strain_plate_rows)
