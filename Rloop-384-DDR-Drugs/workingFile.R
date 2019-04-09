library(data.table)
library(plyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(gtools)

#Set the working directory to the local input folder
setwd('/Users/webbert/Documents/Work/Experiments/LateSubDriv Large files:folders (unsynced)/Gene workup/DDR/DDR IF/DDR IF/20181128 DDR 5EU 384 IR/R analysis')


## Check loaded files
remove.list <- c("CV_", "SD_", "X\\.", "SE_", "Ratio")
df.short <- df[, grep(paste(remove.list, collapse = "|"), colnames(df)):=NULL]


##Clean up dataframe to remove columns with all-zeroes:

print(colnames(df.short))
print(colnames(df))

#colnames to plot
graph.features = data.frame( channel = c("Ch1", "Ch1", "Ch1", 
                                         "Ch2", "Ch2", 
                                         "Ch3", "Ch3"), 
                                       measurement = c("SelectedObjectCount", "MEAN_ObjectAreaCh1", "MEAN_ObjectAvgIntenCh1",
                                                       "MEAN_RingAvgIntenCh2", "MEAN_CircAvgIntenCh2",
                                                       "MEAN_RingAvgIntenCh3", "MEAN_CircAvgIntenCh3",
                                                       "MEAN_CircAvgIntenChNorm"),
                                       title =       c("Nuclei count", "Nuclear area", "Avg. DAPI intensity",
                                                       "Mean Total Cyto. S9.6 intensity", "Mean Nuc. S9.6 intensity", 
                                                       "Mean Total Cyto. Nucleolin intensity", "Mean Nuc. Nucleolin intensity", "Mean Nuc. Specific S9.6 signal"))

graph.features.cells = data.frame(
  channel = c("Ch1", "Ch1", "Ch1", 
              "Ch2", 
              "Ch3", 
              "Ch4",
              "Ch5"), 
  channel.label = c("DAPI", "DAPI", "DAPI",
                    "GFP",
                    "RAD51",
                    "Cyclin A",
                    "5-EU"),
  measurement = c("CellNumber", "ObjectAreaCh1", "ObjectTotalIntenCh1",
                  "CircTotalIntenCh2",
                  "CircTotalIntenCh3",
                  "CircTotalIntenCh4", 
                  "CircTotalIntenCh5"),
  title =       c("Nuclei count", "Nuclear area", "Avg. DAPI intensity",
                  "Nuc. GFP intensity", 
                  "Nuc. RAD51", 
                  "Nuc. Cyclin A",
                  "Nuc. 5-EU intensity"),
  ylimits =     c(NA, NA, NA,
                  NA,
                  NA,
                  NA,
                  0.6e7))


drug.order <- c("DMSO",
                "DRB",
                "Caffeine",
                "VE-821 (ATRi)",
                "KU55933 (ATMi)",
                "Cisplatin",
                "MMC",
                "Doxorubicin")
siRNA.order <- c("siNT", "siGENE", "siprotein/2A", "siGENE + siprotein/2A")
#x's to plot
reorder(wellTreatment, Row)
reorder(siRNA, Column)

df.EU.hi <- setDT(df.treatments)[CircAvgIntenCh5 > 500]
df.EU.lo <- setDT(df.treatments)[CircAvgIntenCh5 <= 500]

#Iterate through this list to create labelled bar plots (and tiffs) for each stain
setwd('./Output/')  

save.image = FALSE

plot.title <- "(-/+ 6 Gy IR)"

Make.violin = function(save.image = FALSE, df, plot.title, graph.features, feature){
  # Use measurement from dataframe of features
  y_Axis = as.character(graph.features$title[graph.features$measurement == as.character(feature)])
  # Choice of label for plot title: 
  channel.label = as.character(graph.features$channel.title[graph.features$measurement == as.character(feature)])
  ylim = graph.features$ylimits[graph.features$measurement == as.character(feature)]
  # Graph title 
  graphTitle = paste(y_Axis, " (", channel.label, ") ", plot.title, sep = "")
  
  # # Colour of graph based on IR or UV
  # if(plot.title=="(IR treated)"){
  #   col.pal = brewer.pal(n=9, "Oranges")[4:8]
  # } else if(plot.title=="(UV treated)"){
  #   col.pal = brewer.pal(n=9, "Purples")[4:8]
  # }
  
  # Plot function, takes dataframe, plots siRNA (plate order) on x axis by treatment time, against measurement on y axis
  WorkingPlot = ggplot(df,
                       aes(x = factor(siRNA, levels = siRNA.order), fill = factor(well.treatment, levels = drug.order), y = df[[feature]])) + 
    # labels
    labs(x = "siRNA", y = y_Axis, fill = "Plate Treatment", title = graphTitle) + #S9.6 #pRPA2 S33    
    # graph type, stat+ position plot treatments alongside one another
    geom_violin(stat = "ydensity", position = "dodge", alpha = 0.3) + 
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "crossbar", 
                 width = 0.25,
                 position = position_dodge(width = .25)) +
    ## colour of bars
    scale_fill_brewer(palette = "Set1") +
    facet_wrap(~factor(Ch2, levels = c("M27_protein", "protein_D210N", "EGFP"))
               +reorder(plate.treatment, PlateNumber)
               +factor(well.treatment, levels = drug.order)
               , nrow = 3) +
    # text orientation
    theme_bw() +
    theme(panel.border = element_blank(), 
          axis.line = element_line(colour = "grey"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    #scale_y_continuous(limits = c(0, ylim))
  
  
  print(WorkingPlot)
  
  ## Saving in function
  if(save.image){
    ggsave(
      paste(gsub("/", "x", graphTitle), ".pdf", sep=""),
      WorkingPlot,
      width = 40,
      height = 20,
      scale = 1,
      dpi = 180
    )
    dev.off()
    print("saved")
    print(WorkingPlot)
  }

}

Make.bar = function(save.image = FALSE, df, plot.title, graph.features, feature){
  # Use measurement from dataframe of features
  y_Axis = as.character(graph.features$title[graph.features$measurement == as.character(feature)])
  # Choice of label for plot title: 
  if( as.character(graph.features$channel[graph.features$measurement == as.character(feature)]) == "Ch2"){
    channel.label <- df[[1,"Ch2"]]
  } else if (as.character(graph.features$channel[graph.features$measurement == as.character(feature)]) == "Ch3"){
    channel.label <- df[[1,"Ch3"]]
  } else if (as.character(graph.features$channel[graph.features$measurement == as.character(feature)]) == "Ch4"){
    channel.label <- "GFP"
  } else {
    channel.label <- "DAPI"
  }
  # Graph title 
  graphTitle = paste(y_Axis, " (", channel.label, ") ", plot.title, sep = "")
  
  # # Colour of graph based on IR or UV
  # if(plot.title=="(IR treated)"){
  #   col.pal = brewer.pal(n=9, "Oranges")[4:8]
  # } else if(plot.title=="(UV treated)"){
  #   col.pal = brewer.pal(n=9, "Purples")[4:8]
  # }
  
  # Plot function, takes dataframe, plots siRNA (plate order) on x axis by treatment time, against measurement on y axis
  WorkingPlot = ggplot(df,
                       aes(x = factor(siRNA, levels = siRNA.order), fill = factor(well.treatment, levels = drug.order), y = df[[feature]])) + 
    # labels
    labs(x = "siRNA", y = y_Axis, fill = "Plate Treatment", title = graphTitle) + #S9.6 #pRPA2 S33    
    # graph type, stat+ position plot treatments alongside one another
    geom_bar(stat = "identity", position = "dodge") + 
    ## colour of bars
    scale_fill_brewer(palette = "Set1") +
    facet_wrap(~reorder(Ch2, Row)+reorder(plate.treatment, PlateNumber)+reorder(well.treatment, Row), nrow = 4) +
    # text orientation
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    scale_y_continuous(limits = c(0, 2000))
  
  print(WorkingPlot)
  
  ###Saving in function
  #if(makeTiff){
  #  ggsave(
  #    paste(gsub("/", "x", graphTitle), ".png", sep=""),
  #    WorkingPlot,
  #    width = 6,
  #    height = 6,
  #    scale = 1,
  #    dpi = 180
  #  )
  #    dev.off()
  #    print("saved")
  #    print(WorkingPlot)
  
  #file_name = paste(gsub("/", "x", graphTitle), ".png", sep="")
  #png(file_name, width = 960, height = 960, bg = "transparent", pointsize = 12)
  #print(WorkingPlot)
  
  
}

Make.bar(df = df, plot.title = plot.title, graph.features = graph.features.cells, feature = "CircAvgIntenCh3")

Make.box = function(save.image = FALSE, df, plot.title, graph.features, feature){
  # Use measurement from dataframe of features
  y_Axis = as.character(graph.features$title[graph.features$measurement == as.character(feature)])
  # Choice of label for plot title: 
  channel.label = as.character(graph.features$channel.label[graph.features$measurement == as.character(feature)])
  ylim = graph.features$ylimits[graph.features$measurement == as.character(feature)]
  print(ylim)
  # Graph title 
  graphTitle = paste(y_Axis, " (", channel.label, ") ", plot.title, sep = "")
  
  # # Colour of graph based on IR or UV
  # if(plot.title=="(IR treated)"){
  #   col.pal = brewer.pal(n=9, "Oranges")[4:8]
  # } else if(plot.title=="(UV treated)"){
  #   col.pal = brewer.pal(n=9, "Purples")[4:8]
  # }
  
  # Plot function, takes dataframe, plots siRNA (plate order) on x axis by treatment time, against measurement on y axis
  WorkingPlot = ggplot(df,
                       aes(x = factor(siRNA, levels = siRNA.order), fill = factor(well.treatment, levels = drug.order), y = df[[feature]])) + 
    # labels
    labs(x = "siRNA", y = y_Axis, fill = "Plate Treatment", title = graphTitle) + #S9.6 #pRPA2 S33    
    # graph type, stat+ position plot treatments alongside one another
    geom_boxplot(position = "dodge", outlier.size = 0.03, outlier.alpha = 0.1, alpha = 0.3) + 
    ## colour of bars
    scale_fill_brewer(palette = "Set1") +
    facet_wrap(~factor(Ch2, levels = c("M27_protein", "protein_D210N", "EGFP"))
               +reorder(plate.treatment, PlateNumber)
               +factor(well.treatment, levels = drug.order)
               , nrow = 3) +
    # text orientation
    theme_bw() +
    theme(panel.border = element_blank(), 
          axis.line = element_line(colour = "grey"),
             axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    scale_y_continuous(limits = c(0, ylim))
    
  
  print(WorkingPlot)
  
  ## Saving in function
  if(save.image){
   ggsave(
     paste(gsub("/", "x", graphTitle), ".pdf", sep=""),
     WorkingPlot,
     width = 6,
     height = 9,
     scale = 1,
     dpi = 180
   )
     dev.off()
     print("saved")
     print(WorkingPlot)
  }
  

}

Make.box(save.image = FALSE, 
         df = 
           # df,
           subset(df, well.treatment %in% c("DMSO")
           & siRNA %in% c("siNT", "siGENE")
           ),
         plot.title = paste(plot.title),
                            # ,"subset Untreated DMSO", 
                            # sep = " "), 
         graph.features = graph.features.cells, 
         feature = "CircTotalIntenCh5")

for(i in 2:8){
  Make.box(save.image = TRUE, 
           df = subset(df, well.treatment %in% c(drug.order[1], drug.order[i])
                       & siRNA %in% c("siNT", "siGENE")),
           plot.title = paste(plot.title,
                              "subset", drug.order[i], "ARHG",
                              sep = " "),
           graph.features = graph.features.cells, 
           feature = "CircTotalIntenCh4")
}
