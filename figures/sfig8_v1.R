library(dplyr)
library(stringr)
library(tidyr)
library(forcats)
library(ggplot2)
library(rstatix) #dplyr-compatible stats package
library(ComplexHeatmap) #heatmap generation
library(circlize) #part of ComplexHeatmap requirements
library(cluster) #required for heatmap generation


### Create dataframe with disease staging information}
dfAAb <- read.csv("./metadata/Donor_AAb.csv", header = TRUE, check.names = FALSE)

# bin donors by pseudo-staging, with 4 groups depending upon AAb number
# and/or years since T1D diagnosis
dfLineageAAbBinning <- dfLineageFilter %>%
  left_join(dfAAb, by = c("HPAP Donor", "Disease Status")) %>%
  mutate(onsetBinStaging = case_when(`Disease Status` == "ND" ~ "1",
                                     `Disease Status` == "AAb+" & Aabs == 1 ~ "2",
                                     `Disease Status` == "AAb+" & Aabs >= 2 ~ "3",
                                     `Disease Status` == "T1D" & `LOD (y)` <= 3 ~ "3",
                                     `Disease Status` == "T1D" & `LOD (y)` > 3 ~ "4")) 


### Significantly changing populations with "onset" binning

dfLineageOnsetStats <- dfLineageAAbBinning %>%
  dplyr::filter(!str_detect(metric, "CD. Tn CD127.CD27.")) %>% #filter out metric with no variance, will throw error
  dplyr::filter(!str_detect(metric,  "CD. Tn CD27.")) %>% #filter out metric with no variance, will throw error
  dplyr::filter(!str_detect(metric,  "CD. Tn CD127.")) %>% #filter out metric with no variance, will throw error
  dplyr::filter(!str_detect(metric, "CD. T.* CD127.CD27.")) %>% #redundant with single positive markers, hard to interpret
  dplyr::filter(!str_detect(metric, "^CD. Mem CD.*$")) %>% #difficult to interpret with memory subsets
  dplyr::filter(!str_detect(metric, "^CD. Mem .*\\+$")) %>% #difficult to interpret with memory subsets
  dplyr::filter(LN_type == "pLN") %>% #just looking at pLN
  group_by(metric) %>%
  dunn_test(value ~ onsetBinStaging, p.adjust.method = "holm")

# adding FDR p-value adjustment due to number of comparisons
dfLineageOnsetStats <- mutate(dfLineageOnsetStats, 
                              p.adj.fdr = p.adjust(dfLineageOnsetStats$p.adj, method = "BH"))

# filter out populations affected by cold ischemia, since group 3 has a higher cold ischemia time than all the others
coldIsPops3 <- (str_replace_all(coldIsPops2, "_", " "))

dfLineageOnsetStats <- dfLineageOnsetStats %>%
  dplyr::filter(!(metric %in% coldIsPops3))

onsetBinPopsSig <- dfLineageOnsetStats %>%
  dplyr::filter(p.adj.fdr < 0.1) %>% # filter for populations below 10% FDR
  dplyr::select(metric) %>%
  unique() # significant populations changing across stages


### Heatmap of significantly changing populations across pseudo-stages

# make df of z-scored tissues from significant populations in onset binning
dfNormTissuePopOnsetBinning <- dfLineageAAbBinning %>%
  dplyr::filter(metric %in% onsetBinPopsSig$metric) %>%
  dplyr::filter(LN_type == "pLN") %>%
  group_by(metric) %>%
  mutate(z_score = scale(value)) %>%
  mutate(onsetBinStaging = case_when(onsetBinStaging == "1" ~ "ND",
                                     onsetBinStaging == "2" ~ "Single AAb+",
                                     onsetBinStaging == "3" ~ "Onset",
                                     onsetBinStaging == "4" ~ "Late T1D"))

# whittle down to onset binning, make other edits for heatmap compatibility
dfsigPopsOnsetHeatmap <- dfNormTissuePopOnsetBinning %>%
  dplyr::select(`HPAP Donor`, Tissue, metric, value, onsetBinStaging, z_score) %>%
  group_by(metric, onsetBinStaging) %>%
  summarise(zMean = mean(z_score)) %>%
  pivot_wider(names_from = onsetBinStaging, values_from = zMean) %>% 
  rename("immune_pop" = "metric") %>%
  ungroup() %>%
  mutate(Sample_ID = paste("S", as.character(row_number()), sep = "")) %>%
  relocate(Sample_ID, .before = immune_pop) %>%
  as.data.frame()

rownames(dfsigPopsOnsetHeatmap) <- dfsigPopsOnsetHeatmap$Sample_ID
dfsigPopsOnsetHeatmap <- dplyr::select(dfsigPopsOnsetHeatmap, immune_pop:ncol(dfsigPopsOnsetHeatmap))
immune_pops_names <- structure(dfsigPopsOnsetHeatmap$immune_pop, names = rownames(dfsigPopsOnsetHeatmap))

dfsigPopsOnsetHeatmapPlot <- dfsigPopsOnsetHeatmap %>%
  dplyr::select(`Late T1D`:`Single AAb+`) %>%
  as.matrix()

col_fun_DiseaseHm <- colorRamp2(c(-0.75, 0, 0.75), c("blue", "white", "red")) 
set.seed(42)

OnsetHeatmap <- Heatmap(dfsigPopsOnsetHeatmapPlot,
                        name = "Mean Z-score",
                        col = col_fun_DiseaseHm,
                        column_order = c("ND", "Single AAb+", "Onset", "Late T1D"),
                        show_column_names = TRUE,
                        column_names_side = "top",
                        column_names_rot = 45,
                        column_names_centered = FALSE,
                        column_names_gp = gpar(fontsize = 8, fontface = "bold"),
                        row_km = 3,
                        row_gap = unit(0.1, "in"),
                        row_title_gp = gpar(fontsize = 10),
                        row_title_rot = 0,
                        row_labels = immune_pops_names,
                        row_names_gp = gpar(fontsize = 8),
                        row_dend_width = unit(0.3, "in"),
                        heatmap_legend_param = list(
                          legend_direction = "horizontal",
                          legend_width = unit(0.6, "in"),
                          title_position = "topcenter",
                          grid_height = unit(0.1, "in"),
                          labels_gp = gpar(fontsize = 8),
                          title_gp = gpar(fontsize = 10)
                        ))
draw(OnsetHeatmap, heatmap_legend_side = "bottom")
saveRDS(OnsetHeatmap, file = "./outs/rds/OnsetHeatmap.RDS")


### CD38 in populations of interest

# CD38 in CD4s and CD8s
diseaseColors4 <- c("ND" = "#949494", "Single AAb+" = "#7bccc4", "Onset" = "#43a2ca", "Late T1D" = "#0868ac")

dfDiseaseScalesOnset <- dfDiseaseScales %>%
  dplyr::select(`HPAP Donor`, LN_type, Tissue, metric, value, cd, tpop) %>%
  left_join(dfLineageAAbBinning, by = c("HPAP Donor", "LN_type", "Tissue", "metric", "value")) %>%
  mutate(onsetBinStaging = case_when(onsetBinStaging == "1" ~ "ND",
                                     onsetBinStaging == "2" ~ "Single AAb+",
                                     onsetBinStaging == "3" ~ "Onset",
                                     onsetBinStaging == "4" ~ "Late T1D")) %>%
  mutate(onsetBinStaging = factor(onsetBinStaging, levels = c("ND", "Single AAb+", "Onset", "Late T1D")))

onset_CD38_boxplot_list <- list("CD38+_CD4_CD4+ T cells_pLN", "CD38+_CD8_CD8+ T cells_pLN") #string will be split at _ for graph creation

onset_CD38_facetLabels <- c("Tn", "Tnl", "Tcm", "Tem", "Temra")

lapply(onset_CD38_boxplot_list, function(c) {
  strings <- unlist(strsplit(c, "_"))
  message(paste("Generating a graph for", strings[4], strings[1], strings[3]))
  print(strings[1:4])
  graph <- dfDiseaseScalesOnset %>%
    dplyr::filter(grepl(strings[1], metric)) %>%
    dplyr::filter(LN_type == strings[4]) %>%
    dplyr::filter(cd == strings[2]) %>%
    dplyr::filter(!str_detect(metric, "HLA")) %>% #HLA-DR+ CD38+ need to be filtered out
    ggplot(aes(onsetBinStaging, value, color = onsetBinStaging)) +
    ylab(paste("% of", strings[3])) +
    geom_boxplot(width = 0.8,
                 outlier.shape = NA,
                 show.legend = FALSE) +
    geom_point(
      size = 1,
      stroke = 0.2,
      alpha = 0.4,
      show.legend = TRUE,
      position = position_jitterdodge(jitter.width = 1, dodge.width = 0.8)) +
    scale_color_manual(values = diseaseColors4) +
    guides(color = guide_legend(title = "Pseudo-stage", override.aes = list(size = 4))) +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 8)) +
    facet_grid(cols = vars(tpop), switch = "y") +
    theme(strip.background = element_blank())
  saveRDS(graph, file = paste0("./outs/rds/onset_",strings[1], " ", strings[3], " boxplots.RDS"))

  
### Plotting specific markers of interest
  
# B cells CD27, CD25 CD4 Tem, CD4 Tn, NK CD56briCD16+

onset_BoxPlots_pops <- list("pLN_B cells CD27+_B cells", 
                            "pLN_CD4 Tem CD25+_CD4+ T cells", 
                            "pLN_NK CD56dimCD16+_NK cells", 
                            "pLN_CD4 Tn_CD4+ T cells") #underscore denotes where string split will occur
})

onset_BoxPlots_graphs <- lapply(onset_BoxPlots_pops, function(d) {
  message("Generating a graph for ", d)
  strings <- unlist(strsplit(d, "_"))
  graph <- dfLineageFilter %>%
    dplyr::select(`HPAP Donor`, LN_type, Tissue, metric, value) %>%
    left_join(dfLineageAAbBinning, by = c("HPAP Donor", "LN_type", "Tissue", "metric", "value")) %>%
    mutate(onsetBinStaging = case_when(onsetBinStaging == "1" ~ "ND",
                                       onsetBinStaging == "2" ~ "Single AAb+",
                                       onsetBinStaging == "3" ~ "Onset",
                                       onsetBinStaging == "4" ~ "Late T1D")) %>%
    mutate(onsetBinStaging = factor(onsetBinStaging, levels = c("ND", "Single AAb+", "Onset", "Late T1D"))) %>%
    dplyr::filter(LN_type == strings[1]) %>%
    dplyr::filter(metric == strings[2]) %>%
    ggplot(aes(onsetBinStaging, value, color = onsetBinStaging)) +
    geom_boxplot(width = 0.8,
                 outlier.shape = NA,
                 show.legend = FALSE) +
    geom_point(
      size = 1,
      stroke = 0.2,
      alpha = 0.4,
      show.legend = FALSE,
      position = position_jitterdodge(jitter.width = 1, 
                                      dodge.width = 0.8)) +
    ggtitle(strings[2]) +
    ylab(paste("% of", strings[3])) +
    scale_color_manual(values = diseaseColors4) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_text(size = 8)
    )
  saveRDS(graph, file = paste0("./outs/rds/onset_",strings[2], " ", strings[1], " boxplots.RDS"))
})
