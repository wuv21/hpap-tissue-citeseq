# note that this code is written to be run from the project base directory
# renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

source("figures/genericFigureSettings.R")
source("scripts/dimPlots.R")
source("scripts/deg.R")
library(cowplot)
library(ppcor)

set.seed(42)

################################################################################
# figure 7 derived from @Greg's code
################################################################################
#Correlation analysis for HLA genetic risk score and immune population frequency
grsHLAcorr <- function(Allpops, NApops, dfAllpops, dfNApops, disease_status_filter_out = "") {
  population <- c()
  statistic <- c()
  estimate <- c()
  p.value <- c()
  method <- c()
  statistic_corrected <- c()
  estimate_corrected <- c()
  p.value_corrected <- c()
  method_corrected <- c()
  
  dataFrame <- c() 
  
  dfAllpops <- dfAllpops %>%
    filter(!is.na(cold_ischemia))
  dfNApops <- dfNApops %>%
    filter(!is.na(cold_ischemia))
  
  for (q in Allpops) {
    if (q %in% NApops) {
      message(paste("Running correlation analyses on", q, "accounting for NA values"))
      dataFrame <- dfNApops %>%
        filter(metric == q) %>%
        filter(Disease_Status != disease_status_filter_out)
    } else {
      message(paste("Running correlation analyses on", q))
      dataFrame <- dfAllpops %>%
        filter(metric == q) %>%
        filter(Disease_Status != disease_status_filter_out)
    }
    
    message("Accounting for cold ischemia effect.")
    corrtestA <- cor.test(dataFrame$value, dataFrame$HLA_score, method = "kendall")
    corrtestB <- ppcor::pcor.test(dataFrame$value, dataFrame$HLA_score, dataFrame[,c("dummyDisease", "cold_ischemia")], method = "kendall")
    
    population <- append(population, q)
    statistic <- append(statistic, corrtestA$statistic)
    estimate <- append(estimate, corrtestA$estimate)
    p.value <- append(p.value, corrtestA$p.value)
    method <- append(method, "kendall")
    statistic_corrected <- append(statistic_corrected, corrtestB$statistic)
    estimate_corrected <- append(estimate_corrected, corrtestB$estimate)
    p.value_corrected <- append(p.value_corrected, corrtestB$p.value)
    method_corrected <- append(method_corrected, "kendall")
  }  
  dfCorrStats <- data.frame(population = population,
    statistic = statistic,
    estimate = estimate,
    p.value = p.value,
    method = method,
    statistic_corrected = statistic_corrected,
    estimate_corrected = estimate_corrected,
    p.value_corrected = p.value_corrected,
    method_corrected = method_corrected)
  
  dfCorrStats <- mutate(dfCorrStats, p.corrected.adjusted = p.adjust(dfCorrStats$p.value_corrected, method = "BH"))
  return(dfCorrStats)
}



################################################################################
# A - hla grs score
################################################################################
parentDir <- "figures/greg_flow_data"

dfLineageFilter <- readRDS("figures/greg_flow_data/rds/dfLineageFilter.rds")
dfGRShla <- read.csv(file = "metadata/HLAgrs.csv", header = TRUE, check.names = FALSE, row.names = 1)

dfGRShla <- dplyr::select(dfGRShla, c(donor, HLA_score)) %>%
  rename(`HPAP Donor` = donor)

dfLineageHLA <- dfLineageFilter %>%
  left_join(dfGRShla, by = c("HPAP Donor"))

dfAllPopsFreqStats <- dfLineageFilter %>%
  pivot_wider(names_from = metric, values_from = value)
colnames(dfAllPopsFreqStats) <- gsub("\\+ \\(\\%T ?cells\\)", "", colnames(dfAllPopsFreqStats))
colnames(dfAllPopsFreqStats) <- gsub(" ", "_", colnames(dfAllPopsFreqStats))

coldIsFreqStats <- filter(dfAllPopsFreqStats, !is.na(cold_ischemia)) %>%
  filter(Disease_Status == "ND") %>%
  filter(Tissue != "Spleen")
coldIsPops <- colnames(coldIsFreqStats[15:ncol(coldIsFreqStats)])


allpops <- colnames(dfAllPopsFreqStats[13:ncol(coldIsFreqStats)])
dfAllPopsFreqStatsGRS <- dfLineageHLA %>%
  group_by(metric, LN_type) %>%
  mutate(z_score = scale(value)) %>%
  ungroup() %>%
  relocate(HLA_score, .before = metric) %>%
  pivot_wider(id_cols = `HPAP Donor`:HLA_score, names_from = metric, values_from = z_score)
colnames(dfAllPopsFreqStatsGRS) <- gsub("\\+ \\(\\%T ?cells\\)", "", colnames(dfAllPopsFreqStatsGRS))
colnames(dfAllPopsFreqStatsGRS) <- gsub(" ", "_", colnames(dfAllPopsFreqStatsGRS))

dfAllPopsFreqStatsGRS <- dfAllPopsFreqStatsGRS %>%
  mutate(Disease_Status = factor(Disease_Status, levels = c("ND", "AAb+", "T1D")))

dfAllPopsFreqStatsGRSpln <- filter(dfAllPopsFreqStatsGRS, LN_type == "pLN")

HLAscoreDisease <- dfAllPopsFreqStatsGRSpln %>%
  dplyr::select(HPAP_Donor, Disease_Status, HLA_score) %>%
  unique() %>%
  ggplot(aes(Disease_Status, HLA_score, color = Disease_Status)) +
  geom_boxplot(width = 0.8,
    outlier.shape = NA,
    show.legend = FALSE) +
  geom_point(
    size = 1.5,
    stroke = 0.2,
    alpha = 0.4,
    show.legend = FALSE,
    position = position_jitterdodge(jitter.width = 1, 
      dodge.width = 0.8)) +
  scale_color_manual(values = COLORS$disease) +
  labs(y = "HLA-GRS score") +
  theme_classic() +
  subplotTheme +
  theme(
    legend.position = "blank",
    plot.title = element_text(size = BASEPTFONTSIZE, color = "#000000", hjust = 0.5, margin = margin(b = BASEPTFONTSIZE)),
    plot.title.position = "panel",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = BASEPTFONTSIZE),
    axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000"))

dfGRSstatsPLN <- dfAllPopsFreqStatsGRSpln %>%
  dplyr::select(HPAP_Donor, Disease_Status, HLA_score) %>%
  unique() %>%
  rstatix::dunn_test(HLA_score ~ Disease_Status, p.adjust.method = "holm")

################################################################################
# B - correlation 
################################################################################
dfAllPopsFreqStatsGRSplndummy <- dfAllPopsFreqStatsGRSpln %>%
  mutate(dummyDisease = case_when(Disease_Status == "ND" ~ 1,
    Disease_Status == "AAb+" ~ 2,
    Disease_Status == "T1D" ~ 3)) %>%
  relocate(dummyDisease, .after = Disease_Status)

# create dataframe without samples containing NAs
dfAllPopsFreqStatsGRSplndummyNoNAs <- dfAllPopsFreqStatsGRSplndummy %>%
  filter(!is.na(`CD4_Mem_PD1+`))

# populations with 0 variation and will fail in in some statistical tests
popsTotal <- dfAllPopsFreqStatsGRSpln %>%
  dplyr::select(`CD45+`:ncol(dfAllPopsFreqStatsGRSpln)) %>%
  colnames() %>%
  str_remove_all("CD._Tn_CD127.CD27.") %>%
  str_remove_all("CD._Tn_CD27.") %>%
  str_remove_all("CD._Tn_CD127.")

# all populations to be tested and some populations with NA values (will cause some correlation tests to fail)
popsTotal <- popsTotal[nzchar(popsTotal)]
popsNAs <- c(
  "CD4_Mem_PD1+",
  "CD4_Tcm_PD1+",
  "CD4_Tem_PD1+",
  "CD4_Temra_PD1+",
  "CD4_Tn_PD1+",
  "CD4_Tnl_PD1+", 
  "CD4_Tem_HLA-DR+",
  "CD4_Temra_HLA-DR+",
  "CD4_Tcm_HLA-DR+",
  "CD8_Temra_HLA-DR+",
  "CD4_Tem_HLA-DR+_CD38+",
  "CD4_Temra_HLA-DR+_CD38+",
  "CD4_Tcm_HLA-DR+_CD38+",
  "CD8_Temra_HLA-DR+_CD38+",
  "CD4_Tcm_CD25+",
  "CD4_Tnl_CD25+")

# df's for correlation analysis
dfAllPopsFreqStatsGRSplndummy <- dfAllPopsFreqStatsGRSplndummy %>%
  pivot_longer(cols = c(`CD45+`:ncol(.)), names_to = "metric", values_to = "value")

dfAllPopsFreqStatsGRSplndummyNoNAs <- dfAllPopsFreqStatsGRSplndummyNoNAs %>%
  pivot_longer(cols = c(`CD45+`:ncol(.)), names_to = "metric", values_to = "value")

###USE PPCOR package function to run a partial correlation controlling for the effect of disease state and use Kendall correlation (less sensitive to outliers)

#both disease state and cold ischemia (since it covaries with HLA score in T1D) are controlled for 
dfGRShlaCorrStats <- grsHLAcorr(popsTotal, popsNAs, dfAllPopsFreqStatsGRSplndummy, dfAllPopsFreqStatsGRSplndummyNoNAs)

#the above correlation has improper controlling for disease state, as immune populations that only change in AAb+ donors will not be controlled for. 
#create separate dfs for ND -> AAb correlation and ND -> T1D correlation, re-run analyses
dfGRShlaCorrStatsNDAAb <- grsHLAcorr(popsTotal, popsNAs, dfAllPopsFreqStatsGRSplndummy, dfAllPopsFreqStatsGRSplndummyNoNAs, "T1D")
dfGRShlaCorrStatsNDT1D <- grsHLAcorr(popsTotal, popsNAs, dfAllPopsFreqStatsGRSplndummy, dfAllPopsFreqStatsGRSplndummyNoNAs, "AAb+")
dfGRShlaCorrStatsAAbT1D <- grsHLAcorr(popsTotal, popsNAs, dfAllPopsFreqStatsGRSplndummy, dfAllPopsFreqStatsGRSplndummyNoNAs, "ND")


hlaGRS_NDAAb <- dfGRShlaCorrStatsNDAAb %>%
  mutate(population = str_replace_all(population, "_", " ")) %>%
  mutate(graphingP = -log10(p.corrected.adjusted)) %>%
  arrange(desc(graphingP)) %>%
  mutate(population = factor(population, levels = rev(population))) %>% 
  filter(graphingP >= 1) %>% 
  mutate(corrPosNeg = case_when(estimate_corrected >= 0 ~ "Positive Correlation",
    estimate_corrected < 0 ~ "Negative Correlation")) %>%
  mutate(corrPosNeg = factor(corrPosNeg, levels = c("Positive Correlation", "Negative Correlation"))) %>%
  {ggplot(., aes(x = graphingP, y = population, fill = estimate_corrected)) +
      geom_vline(xintercept = -log10(0.05), color = "#d3d4d3") +
      geom_segment(aes(x = 0, xend = graphingP, y = population, yend = population), 
        linetype = "dotted", color = "#555555") +
      geom_point(size = 1.3, shape = 21, color = "#000000") +
      labs(x = "-log10(p adj)",
        fill = "Corr.Estimate (τ)") + 
      ggtitle("ND and AAb+") +
      guides(fill = guide_colorbar(
        title.position = "top",
        title.hjust = 1,
        title.vjust = 1,
        barwidth = 0.5,
        barheight = 4)) +
      scale_fill_gradient(low = "#ffebea", high = "red") +
      scale_x_continuous(expand = c(0, 0), limits = c(0, 4.1)) +
      theme_classic() +
      subplotTheme +
      theme(
        axis.text = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.margin = margin(l = -10),
        legend.title = element_text(size = BASEPTFONTSIZE - 4),
        legend.text = element_text(size = BASEPTFONTSIZE - 4),
        plot.title = element_text(size = BASEPTFONTSIZE, hjust = 0.5, margin = margin(b = BASEPTFONTSIZE))) +
      facet_grid(corrPosNeg ~ ., scales = "free", space = "free") +
      theme(strip.background = element_blank(),
        strip.text = element_blank())
  }

################################################################################
# Figure CD/FG
################################################################################
GRSpops <- c("ND_AAb CD8_Tem_HLA-DR+_CD38+ CD8+_Tem_cells:_HLA-DR+_CD38+", "ND_AAb CD8_Tcm_HLA-DR+_CD38+ CD8+_Tcm_cells:_HLA-DR+_CD38+", "AAb_T1D B_cells_CD69+ B_cells:_CD69+", "AAb_T1D CD4_Temra_CD127-CD27+ CD4+_Temra:_CD127-_CD27-") #split by the sapce

HLAgrsGraphs <- lapply(GRSpops, function(x) {
  message("Generating a graph for ", x)
  strings <- unlist(strsplit(x, " "))
  
  if (strings[1] == "ND_AAb") {
    diseaseFilter <- "T1D"
  } else if (strings[1] == "AAb_T1D") {
    diseaseFilter <- "ND"
  } else if (strings[1] == "ND_T1D") {
    diseaseFilter <- "AAb+"
  } else {
    stop("Error: No disease state comparisons found.")
  }
  
  GRSpopGraph <- dfAllPopsFreqStatsGRSplndummy %>%
    filter(Disease_Status != diseaseFilter) %>%
    filter(metric == strings[2]) %>%
    mutate(metric = str_replace_all(metric, "_", " ")) %>%
    ggplot(aes(HLA_score, value)) +
    geom_point(size = 0.5) +
    geom_smooth(method = "lm", size = 0.75) +
    ggtitle(str_wrap(str_replace_all(strings[3], "_", " "), width = 15)) +
    xlab("HLA-GRS score") +
    ylab("z score") +
    theme_classic() +
    subplotTheme +
    theme(
      axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
      axis.title = element_text(size = BASEPTFONTSIZE, color = "#000000"),
      plot.title = element_text(size = BASEPTFONTSIZE, hjust = 0.5, margin = margin(b = BASEPTFONTSIZE)),
      plot.title.position = "panel")
})

################################################################################
# Figure E
################################################################################
hlaGRS_AAbT1D <- dfGRShlaCorrStatsAAbT1D %>%
  mutate(population = str_replace_all(population, "_", " ")) %>%
  mutate(graphingP = -log10(p.corrected.adjusted)) %>%
  arrange(desc(graphingP)) %>%
  mutate(population = factor(population, levels = (population))) %>% 
  filter(graphingP >= 1) %>% 
  mutate(corrPosNeg = case_when(estimate_corrected >= 0 ~ "Positive Correlation",
    estimate_corrected < 0 ~ "Negative Correlation")) %>%
  mutate(corrPosNeg = factor(corrPosNeg, levels = c("Positive Correlation", "Negative Correlation"))) %>%
  {ggplot(., aes(x = population, y = graphingP, fill = estimate_corrected)) +
      geom_hline(yintercept = -log10(0.05), color = "#d3d4d3") +
      geom_segment(aes(y = 0, yend = graphingP, x = population, xend = population), linetype = "dotted", color = "#555555") +
      geom_point(size = 1.3, shape = 21, color = "#000000") +
      labs(y = "-log10(p adj)",
        fill = "Corr.Estimate (τ)") + 
      guides(fill = guide_colorbar(
        title.position = "top",
        title.hjust = 1,
        title.vjust = 1,
        barwidth = 0.5,
        barheight = 4)) +
      scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red") +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 4.1)) +
      ggtitle("AAb+ and T1D") +
      facet_wrap(corrPosNeg ~ ., scales = "free") +
      theme_classic() +
      subplotTheme +
      theme(
        axis.text.x.bottom = element_text(angle = 45, size = BASEPTFONTSIZE - 4, hjust = 1, vjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 1),
        legend.position = "right",
        plot.margin = unit(c(0, 0, 0, 0.1), units = "in"),
        axis.text = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
        legend.margin = margin(l = -5),
        legend.title = element_text(size = BASEPTFONTSIZE - 4),
        legend.text = element_text(size = BASEPTFONTSIZE - 4),
        plot.title = element_text(size = BASEPTFONTSIZE, hjust = 0.5, margin = margin(b = BASEPTFONTSIZE)),
        strip.background = element_blank())
  }


################################################################################
# Final layout and plot all
################################################################################
layout <- c(
  patchwork::area(1, 1, 2, 2), # a
  patchwork::area(1, 3, 2, 6), # a
  patchwork::area(1, 7, 2, 8), # a
  patchwork::area(1, 9, 2, 10), # a
  patchwork::area(3, 1, 4, 6), # a
  patchwork::area(3, 7, 4, 8), # a
  patchwork::area(3, 9, 4, 10) # a
)

p <- wrap_elements(plot = HLAscoreDisease) +
  wrap_elements(plot = hlaGRS_NDAAb) +
  wrap_elements(plot = HLAgrsGraphs[[1]]) +
  wrap_elements(plot = HLAgrsGraphs[[2]]) +
  wrap_elements(plot = hlaGRS_AAbT1D) +
  wrap_elements(plot = HLAgrsGraphs[[3]]) +
  wrap_elements(plot = HLAgrsGraphs[[4]]) +
  plot_annotation(tag_levels = list(LETTERS[1:length(layout)])) +
  plot_layout(design = layout) &
  plotTagTheme


saveFinalFigure(
  plot = p,
  prefixDir = "figures/outs",
  fn = "fig7_v3_final",
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 8,
  gheight = 4)
