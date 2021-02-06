library(data.table)

setwd('./R analysis/input')

#Add names of well treatments, stains and plate treatments
plate.key <- read.csv("PlateKey.csv")

df <- all.data
setorder(df, Row, Column)
df$WellId <- as.factor(gsub("\\s", "", df$WellId))
df <- merge(df, plate.key, by =)
# df$fixed.treatment <- as.factor(gsub("\\s", "", df$fixed.treatment))

## Cells
