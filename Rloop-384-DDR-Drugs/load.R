library(plyr)

#Clear out existing environment
rm(list = ls())

#Set the working directory to the local input folder
setwd('/Users/webbert/Documents/Work/Experiments/LateSubDriv Large files:folders (unsynced)/Gene workup/DDR/DDR IF/DDR IF/20181128 DDR 5EU 384 IR/R analysis')


#Get output file list
setwd('./input/Plates')
file.list = list.files(getwd())
View(file.list)
##Make a list of folder paths
directories <- sapply(file.list, normalizePath)
file.order <- c(1,2)
ordered.directories <- directories[file.order]
head(directories)
head(ordered.directories)
##Make a list of the header files (plate.csv) in these folders
files <- lapply(ordered.directories, list.files, pattern="Cell.csv", full.names = TRUE)

##empty list 
all.data.list <- list()
##Create list of csvs for data files and add a column for the plate ID (folder name)
all.data.list <- lapply(files, read.csv)

for(i in 1:length(all.data.list)){
  if(all.data.list[[i]][1,1] == 1){
    all.data.list[[i]][1] <- i
  } else{ 
    print("script already run")}
}
all.data <- rbindlist(all.data.list)

# setwd('./Input/') 
# df <- read.csv("Well.csv")
# df.cells <- read.csv("Cell.csv")
# head(df)

