library(data.table)

setwd('/Users/webbert/Documents/Work/Experiments/LateSubDriv Large files:folders (unsynced)/Gene workup/DDR/DDR IF/DDR IF/20181128 DDR 5EU 384 IR/R analysis/input')

#Add names of well treatments, stains and plate treatments
plate.key <- read.csv("PlateKey.csv")

df <- all.data
setorder(df, Row, Column)
df$WellId <- as.factor(gsub("\\s", "", df$WellId))
df <- merge(df, plate.key, by =)
# df$fixed.treatment <- as.factor(gsub("\\s", "", df$fixed.treatment))

## Cells
