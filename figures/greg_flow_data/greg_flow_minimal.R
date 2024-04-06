# note this code is from @Greg but has been refactored by @Vincent
# generates all flow data and graphs for figures

library(dplyr)
library(stringr)
library(tidyr)
library(forcats)
library(ggplot2)
library(rstatix)
library(ggpubr)

parentDir <- "figures/greg_flow_data"

################################################################################
# load in data
################################################################################
df <- read.csv(
  file = paste0(parentDir, "/Lineage_MasterV3.csv"), 
  check.names = FALSE,
  stringsAsFactors = FALSE)

df <- dplyr::select(df, 1:`CD4 Tnaive1- CD127+`)

df <- df %>%
  mutate(across(.cols = "CD45+":ncol(df), .fns = ~ str_replace_all(., "%$", ""))) %>%
  mutate(across(.cols = "CD45+":ncol(df), .fns = ~ as.numeric(.))) %>%
  mutate(across(.cols = "CD45+":ncol(df), .fns = ~ case_when(. < 0.01 ~ 0.005, TRUE ~ as.numeric(.)))) %>%
  mutate("Age" = as.numeric(df$Age)) %>%
  mutate("Aabs" = as.numeric(df$Aabs)) %>%
  mutate(LN_type = case_when(str_detect(Tissue, "_pLN$") ~ "pLN",
    str_detect(Tissue, "_mLN$") ~ "mLN",
    Tissue == "Spleen" ~ "Spleen",
    Tissue == "PBMC" ~ "PBMC",
    Tissue == "Whole_Blood" ~ "Whole_Blood",
    Tissue == "Tonsil" ~ "Tonsil",
    Tissue == "Thymus" ~ "Thymus",
    Tissue == "Islet_Sup" ~ "Islet")) %>%
  relocate(LN_type, .after = Tissue)

################################################################################
# load in metadata
################################################################################
dfMeta <- read.csv(
  file = paste0(parentDir, "/Donor_Metadata.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE)

dfMeta <- dfMeta %>%
  filter(str_detect(donor, "HPAP"))

dfLong <- df %>%
  pivot_longer(cols = c(`CD45+`:ncol(.)), names_to = "metric", values_to = "value")

dfFresh <- filter(dfLong, `Fresh/Thawed` == "fresh")

################################################################################
# notes from @Greg
# a more specific naive T cell gating strategy was needed, so the Tnaive 
# population (CD45RA+CCR7+) was separated into two populations based upon 
# CD27 and CD127 expression: Tnaive1 (CD27+CD127+ naive T cells, 
# bona fide naive) and Tnaive 1- (all other CD45RA+CCR7+ T cells). 
# 
# The Tnaive population will be filtered out, since it encompasses all
# CD45RA+CCR7+ T cells and is redundant/less specific.
# Tnaive1 will be re-named "Tn" and Tnaive1- will be renamed "Tnl"
# 
# Thymus and mixed pLN are not standard HPAP tissues, 
# and therefore will be filtered out
################################################################################
dfFresh <- dfFresh %>%
  filter(
    !str_detect(metric, "Tnaive ") & 
    !str_detect(metric, "Tnaive$") &
    Tissue != "Thymus" &
    Tissue != "pLN_mix_pLN") %>%
  mutate(metric = str_replace_all(.$metric, "Tnaive1-", "Tnl")) %>%
  mutate(metric = str_replace_all(.$metric, "Tnaive1", "Tn"))

# Add additional patient metadata:
# cold ischemia to be explored as variable that impacts immune pop frequency
dfMetaSubset <- dfMeta %>%
  dplyr::select(c("donor", "cold_ischemia")) %>% 
  rename("HPAP Donor" = "donor")

dfFresh <- dfFresh %>%
  left_join(dfMetaSubset, by = "HPAP Donor") %>%
  relocate(cold_ischemia, .after = `Fresh/Thawed`)

saveRDS(dfFresh, file = paste0(parentDir, "/rds/dfFresh.rds"))

################################################################################
# notes from @Greg
# as HPAP is a study occurring over several years, different researchers have 
# led the project therefore, the effect of researcher on immune population 
# frequency was checked
################################################################################

# these researchers (D and A) did the vast majority of the data acquisition, 
# and are the only ones considered for this comparison
# 
# additional filters include for metrics where all filtered values are 0 or 100,
# which means that some stats test would fail
dfScientistStat <- dfFresh %>%
  filter(Staining_Flow_Researcher == "D" | Staining_Flow_Researcher == "A") %>% 
  filter(
    LN_type == "pLN" &
    `Disease Status` == "ND" &
    !grepl("Tn CD127.CD27.", metric) &
    !grepl("Tn CD127.", metric) &
    !grepl("Tn CD27.", metric)
  ) %>%
  group_by(metric) %>%
  wilcox_test(value ~ Staining_Flow_Researcher) %>%
  adjust_pvalue() %>%
  filter(p.adj < 0.05) %>%
  mutate(pvalSymbol = case_when(
    p.adj < 0.001 ~ "***",
    p.adj >= 0.001 & p.adj < 0.01 ~ "**",
    p.adj >= 0.01 & p.adj < 0.05 ~ "*")) %>%
  arrange(p.adj)

ScientistComp <- dfFresh %>%
  filter(Staining_Flow_Researcher != "C") %>%
  filter(metric %in% dfScientistStat$metric) %>%
  filter(`Disease Status` == "ND") %>%
  mutate(metric = fct_relevel(`metric`, dfScientistStat$metric)) %>%
  ggplot(aes(Staining_Flow_Researcher, value, color = Staining_Flow_Researcher)) +
  geom_boxplot(outlier.size = 0.3) +
  stat_pvalue_manual(data = dfScientistStat, label = "pvalSymbol", 
    y.position = 105, xmin = "group1", xmax = "group2",
    tip.length = 0, label.size = 2) +
  ylab("% of Parent") +
  guides(color = guide_legend(title = "Researcher")) +
  theme_classic() +
  scale_color_brewer(palette = "Dark2") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.length = unit(0, "in"),
  ) +
  facet_wrap(~factor(metric), nrow = 1, strip.position = "bottom") +
  theme(strip.text.x.bottom = element_text(angle = 90, hjust = 1, size = 6),
    strip.background = element_blank(),
    strip.placement = "outside")

saveRDS(ScientistComp, file = paste0(parentDir, "/rds/ScientistComparison.RDS"))

################################################################################
# notes from @Greg
# Results show that there are several parameters that differ between
# researchers, mainly in PD1, HLA-DR, and CD25 on certain CD4+ T cells. This 
# could be due to differential flow cytometry parameter settings on the 
# cytometer (antibody titrations or conjugations did not change
# between researchers). 
# 
# Populations that were significantly different between researchers were
# filtered (see below), keeping the researcher's results (Researcher D) who
# analyzed the data and compiled it for the manuscript.
################################################################################
ScientistFilter <- c("CD4 Mem PD1+", "CD4 Tcm PD1+", "CD4 Tem PD1+",
  "CD4 Temra PD1+", "CD4 Tn PD1+", "CD4 Tnl PD1+", "CD4 Tem HLA-DR+",
  "CD4 Temra HLA-DR+", "CD4 Tcm HLA-DR+", "CD8 Temra HLA-DR+",
  "CD4 Tem HLA-DR+ CD38+", "CD4 Temra HLA-DR+ CD38+", "CD4 Tcm HLA-DR+ CD38+",
  "CD8 Temra HLA-DR+ CD38+", "CD4 Tcm CD25+", "CD4 Tnl CD25+")

dfScientistFilter <- dfFresh %>%
  filter(Staining_Flow_Researcher == "A") %>%
  filter(metric %in% ScientistFilter)

#generate df filtering out non-standard HPAP tissues
dfLineageFilter <- dplyr::anti_join(dfFresh, dfScientistFilter) %>%
  filter(Tissue != "PBMC") %>% #removed since PBMCs have artifacts in the neutrophil counts due to long term blood storage
  filter(Tissue != "Islet_Sup") %>% #removed since it is not a standard HPAP tissue
  filter(Tissue != "Whole_Blood") %>% #removed since it is not a standard HPAP tissue
  mutate(metric = str_replace(metric, " \\(%CD45\\+\\)", "")) %>%
  mutate(metric = str_replace(metric, " \\(%CD45\\)", "")) %>%
  mutate(metric = str_replace(metric, " \\(%T cells\\)", " T cells"))

saveRDS(dfLineageFilter, file = paste0(parentDir, "/rds/dfLineageFilter.rds"))
