library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)
library(ggsci)

##### present BGC total number 3183
##### generate a map of 3183x262 abundance score with metagroup info and other mapping info ###

rm(list = ls())
setwd("~/Documents/Qiqi-2022-work/2022-04/#5 Mapping/ALL")
dt <- read.csv("BGC_mapping_result_all.csv", row.names = 1, header = T)

## ratio for heatmap
num = 3195
ratio = as.integer(num/262)

## mapping result for different isolation site
mx <- dt[,c(1:262)]

## discribing metrices
avg_ab <- dt[,263]
avg_corecov <- data.frame(dt[,264], 1-dt[,264])
avg_cov <- data.frame(dt[,265], 1-dt[,265])
high_ab <- dt[,266]
presence <- data.frame(dt[,267], 100-dt[,267])
mut_dt <- mx
mut_dt[mut_dt > 0] <- 1
mut_dt <- colSums(mut_dt,na.rm = TRUE)

## play with color
meta_group_col = c("#731D1C", "#A3593E", "#CBAC7C", "#AB7B4A", "#975B35", "#CA793F", "#FFEEC2",
                   "#E1B913", "#F7CB15", "#E9B419", "#DA9C1D", 
                   "#a8c256")
g = c("Corn_Rhizosphere","Switchgrass_Rhizosphere",
      "Sorghum_Rhizosphere","Miscanthus_Rhizosphere",
      "Maize_Rhizosphere", "Populus_Rhizospehre",
      "Arabidopsis_Rhizosphere",
      "Arabidopsis_Root", "Populus_Root", "M_polymorpha_Root", "Legume-Nodule",
      "Phyllosphere_VariousHosts")

## generate legend for meta groups
lgd = Legend(labels = g, title = "Meta sample groups", 
              legend_gp = gpar(fill = meta_group_col),
              labels_gp = gpar(fontsize = 18),
              title_gp = gpar(fontsize = 24, fontface = "bold"),
              grid_height = unit(0.5, "cm"), grid_width = unit(0.5, "cm"))
mapcol = colorRamp2(c(10,9,8,7,6,5,4,3,2,1,0), c("#10451D",
                                                 "#155D27",
                                                 "#1A7431",
                                                 "#208B3A",
                                                 "#25A244",
                                                 "#2DC653",
                                                 "#4AD66D",
                                                 "#6EDE8A",
                                                 "#92E6A7",
                                                 "#B7EFC5",
                                                 "white"))

## generate legend for heatmap
lgd2 = Legend(col_fun = mapcol, title = "Abundance Score", at = c(0, 2.5, 5, 7.5, 10),
             direction = "vertical", title_position = "topleft",
             title_gp = gpar(fontsize = 20, fontface = "bold"),
             labels_gp = gpar(fontsize = 14),
             legend_height = unit(5, "cm"), grid_width = unit(0.8, "cm"))

## add average coverage as annotation
cov_anno = rowAnnotation(
  Cov = anno_barplot(avg_cov,
                     gp = gpar(fill = c("#25A244", "white"), alpha = 1, col = c("#25A244", "white")), border = TRUE,
                     bar_width = 1, width = unit(ratio*2*0.1,"cm"),
                     axis_param = list(gp = gpar(fontsize = 18))), 
  CoreCov = anno_barplot(avg_corecov,
                         gp = gpar(fill = c("#10451D", "white"), alpha = 1, col = c("#10451D", "white")), border = TRUE,
                         bar_width = 1, width = unit(ratio*2*0.1,"cm"),axis_param = list(gp = gpar(fontsize = 18))),
  gap = unit(8, "mm"),
  annotation_name_side = "top",
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 18, fontface = "bold"))

space = rowAnnotation(space = anno_empty(border = FALSE), annotation_name_gp = gpar(fontsize = 0), width = unit(0.5,"cm"))
space_col = HeatmapAnnotation(space = anno_empty(border = FALSE), annotation_name_gp = gpar(fontsize = 0), height = unit(0.25,"cm"))


## add mark annotation to heatmap
pos_ctl = c("Pyrrolnitrin",
            "Bacterial gibberelin",
            "Nitrogen fixation (B. japonicum)",
            "DAPG",
            "Nodulation (S. meliloti)",
            "Nodulation (Rhizobium sp.)",
            "Nitrogen fixation (K. pneumoniae)",
            "Thanamycin",
            "Zwittermicin A",
            "Nodulation (Sinorhizobium sp.)",
            "Brabantamide A",
            "Bacillomycin D")
mark = rowAnnotation(foo = anno_mark(at = c(2,4,5,6,30,32,37,57,86,181,187,261), labels = pos_ctl,
                                     labels_gp = gpar(fontsize = 18)))

## add presence/distribution percentage annotation
pre_anno <- rowAnnotation(Presence = anno_barplot(presence, gp = gpar(fill = c("#008751", "white"), alpha = 1, 
                                                                      col = c("#008751", "white")), 
                                                  border = TRUE, bar_width = 1, width = unit(ratio*4*0.1,"cm"),
                                                  axis_param = list(gp = gpar(fontsize = 18))),
                          annotation_name_side = "top",
                          annotation_name_rot = 0,
                          annotation_name_gp = gpar(fontsize = 18, fontface = "bold"))

## add metagroup annotation
meta_groups <- NULL
fill_col <- NULL
for(i in (1:length(col_slice))){
  group_n = col_slice[i]
  sample_n = s[i]
  c = gcol[i]
  m = c(rep(group_n, sample_n))
  fc = c(rep(c, sample_n))
  meta_groups = c(meta_groups, m)
  fill_col = c(fill_col, fc)}

group_anno = HeatmapAnnotation(height = unit(0.75, "cm"),
                               foo = anno_block(gp = gpar(fill = gcol, col = gcol, show_annotation_name = FALSE, lwd = 1.5,
                                                          labels = g,labels_gp = gpar(col = "white", fontsize = 0))))

## add total BGC/metasample annotation(use mut data)
b <- HeatmapAnnotation('Detected BGC' = anno_barplot((mut_dt/10), bar_width = 1, baseline = 0,
                                           gp = gpar(fill = fill_col, col = fill_col),border = FALSE,
                                           axis_param = list(at = c(0, 20, 40, 60),
                                                             labels = c("0", "200", "400", "600"),
                                                             gp = gpar(fontsize = 18)),
                                           height = unit(4, "cm")),
                       annotation_name_side = "left",
                       annotation_name_rot = 0,
                       annotation_name_gp = gpar(fontsize = 18, fontface = "bold")
                       )

## draw heatmap (do not divide by isolation site)
h = 0.012*num
w = ratio*1*0.01*262
s = c(32,27,28,28,11,36,14,32,9,5,7,33)
col_slice = LETTERS[1:length(g)] ## slice heatmap by metagroup
gcol = meta_group_col[1:length(g)] ## use first 7 color
heatmap <- Heatmap(mx, name = "heatmap", col = mapcol, cluster_rows = F, 
                   column_names_gp = gpar(fontsize = 0),border_gp = gpar(col = "white", lwd = 1.5),
                   row_names_gp = gpar(fontsize = 0), 
                   cluster_columns = F,rect_gp = gpar(col = 0), 
                   width = unit(w, "cm"), height = unit(h, "cm"),
                   na_col = "white", top_annotation = c(b,space_col,group_anno),
                   column_split = meta_groups, column_title = NULL, left_annotation = c(cov_anno, space),
                   right_annotation = c(space, pre_anno, mark),show_heatmap_legend = FALSE)


svg(file = "heatmap.svg", width = (w + 4), height = (h + 4))
draw(heatmap)
draw(lgd2, x = unit(66, "cm"), y = unit(58, "cm"))
draw(lgd, x = unit(65, "cm"), y = unit(50, "cm"))
dev.off()

