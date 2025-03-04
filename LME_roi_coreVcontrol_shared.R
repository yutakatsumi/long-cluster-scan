library(nlme)
library(lmerTest)
library(ggeffects)
library(emmeans)
library(tidyr)
library(dplyr)
library(tibble)
library(stringr)
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(freesurferformats)
library(fsbrain)
library(flextable)
library(htmltools)
library(lmeInfo)
library(effectsize)

# flextable options
set_flextable_defaults(
  font.family = "Arial", font.size = 11, 
  border.color = "gray")

subjects_dir = '/xxxx/xxx'
df <- read.csv('data_cluster_allsubs.csv', sep=',')
# Drop ROI-wise data
df <- df[c("subj_id", "subj_id_buckner", "sess_id", "tp", "day", "session")]

# motion data
dfm <- read.csv('meanrms.csv', sep=',')
dfm$sess_id <- paste(dfm$session,"_",sprintf("%03d", dfm$run), sep="")
dfm <- dfm[,c("sess_id", "mean.rms")]

# Merge two DFs based on the session IDs
dfx <- merge(df,dfm,by="sess_id")
# Drop any scans with suprathreshold motion
motion_thresh = 10
dfxx <- dfx[dfx$mean.rms < motion_thresh,]
dfxx %>% 
  group_by(subj_id) %>%
  summarise(n = n())

datadir = '/xxxx/xxx'
roidir = '/xxxx/xxx'

core_files <- list.files(roidir, pattern="*roidefined.label", full.names=TRUE)
control_files <- list.files(roidir, pattern="*control.label", full.names=TRUE)

# Specify subject ID
subj = 'xxxx'

df_subj = dfxx[dfxx$subj_id==subj,]

alldata_lh <- array(numeric(),c(nrow(df_subj),163842))
alldata_rh <- array(numeric(),c(nrow(df_subj),163842))

# Load core and control ROI masks
core_lh_file <- grep(paste('lh.*',subj, sep=''), core_files, value=TRUE)
core_rh_file <- grep(paste('rh.*',subj, sep=''), core_files, value=TRUE)
control_lh_file <- grep(paste('lh.*',subj, sep=''), control_files, value=TRUE)
control_rh_file <- grep(paste('rh.*',subj, sep=''), control_files, value=TRUE)
if (length(core_lh_file) > 0){
  core_lh_idx <- read.fs.label(core_lh_file)
  core_lh <- array(0,c(1,163842))
  core_lh[,core_lh_idx] = 1
}
if (length(core_rh_file) > 0){
  core_rh_idx <- read.fs.label(core_rh_file)
  core_rh <- array(0,c(1,163842))
  core_rh[,core_rh_idx] = 1
}
if (length(control_lh_file) > 0){
  control_lh_idx <- read.fs.label(control_lh_file)
  control_lh <- array(0,c(1,163842))
  control_lh[,control_lh_idx] = 1
}
if (length(control_rh_file) > 0){
  control_rh_idx <- read.fs.label(control_rh_file)
  control_rh <- array(0,c(1,163842))
  control_rh[,control_rh_idx] = 1
}

# Concatenate ROI masks
if (exists('core_lh') & exists('core_rh')){
  core_bh <- c(core_lh, core_rh)
} else if (exists('core_lh')) {
  core_bh <- c(core_lh, array(0,c(1,163842)))
} else if (exists('core_rh')) {
  core_bh <- c(array(0,c(1,163842)), core_rh)
} else {
  stop("Don't know what to do...")
}

if (exists('control_lh') & exists('control_rh')){
  control_bh <- c(control_lh, control_rh)
} else if (exists('control_lh')) {
  control_bh <- c(control_lh, array(0,c(1,163842)))
} else if (exists('control_rh')) {
  control_bh <- c(array(0,c(1,163842)), control_rh)
} else {
  stop("Don't know what to do...")
}

# Load vtx-wise data
for (i in 1:nrow(df_subj)){
  alldata_lh[i,] <- read.fs.morph(file.path(datadir, paste("lh.", subj, "_", df_subj$sess_id[i], ".thickness.fsaverage.sm15.mgz", sep="")))
  alldata_rh[i,] <- read.fs.morph(file.path(datadir, paste("rh.", subj, "_", df_subj$sess_id[i], ".thickness.fsaverage.sm15.mgz", sep="")))
}

# Convert vtx-wise data to ROI-wise data
alldata_bh <- cbind(alldata_lh,alldata_rh)
rm(alldata_lh)
rm(alldata_rh)
alldata_bh_core <- rowMeans(alldata_bh[,core_bh==1])
alldata_bh_control <- rowMeans(alldata_bh[,control_bh==1])
  
# Do modeling
df_subj$core <- alldata_bh_core
df_subj$control <- alldata_bh_control
df_subj_long <- gather(df_subj, key = "roi", value = "thickness", core, control)

# Reorder levels of ROI so that core is shown first
df_subj_long$roi <- factor(df_subj_long$roi, levels = c("core", "control"))

model <- lme(thickness ~ factor(tp)*factor(roi) + mean.rms, random = ~ 1 | day/session , data=df_subj_long, weights = varIdent(form = ~1|roi), method="REML")
modelSummary <- anova(model)
rownames(modelSummary) <- c("Intecept", "TP", "ROI", "Motion", "TP\u00D7ROI")

model.emm <- emmeans(model, pairwise ~ tp | roi, lmer.df = "satterthwaite", adjust = "none")
emm_df <- summary(model.emm)$emmeans  # Convert EMMs to a data frame
names(emm_df)[1] <- "TP"  # Rename the factor column appropriately
write.csv(emm_df, paste0('nlme_results_', subj, '.csv'), row.names=FALSE)
write.csv(df_subj_long, paste0('df_', subj, '.csv'), row.names=FALSE)

##
modelSummary
# Drop intercept
modelSummary <- modelSummary[-1, ]

ft <- modelSummary %>%
  rownames_to_column() %>%
  flextable() %>%
  set_header_labels(rowname = "Variable",
                    NumDF = "df1",
                    DenDF = "df2",
                    `F-value` = "F",
                    `p-value` = "p") %>%
  set_caption(caption = "ANOVA Table") %>%
  padding(padding.top = 0.5, padding.bottom = 0.5)
ft <- set_formatter(ft, `F-value` = function(x) format(round(x, 3), nsmall = 3),
                    `p-value` = function(x) format(round(x, 3), nsmall = 3))

interaction_contrast_df <- summary(contrast(emmeans(model, ~ tp * roi), interaction = TRUE, method = "pairwise"))
# Rename contrast names
interaction_contrast_df <- interaction_contrast_df %>%
  mutate(tp_pairwise = case_when(
    tp_pairwise == "1 - 2" ~ "BL - 3mo",
    tp_pairwise == "1 - 3" ~ "BL - 6mo",
    tp_pairwise == "2 - 3" ~ "3mo - 6mo",
    TRUE ~ tp_pairwise
  ))
ft_interaction_contrasts <- interaction_contrast_df %>%
  flextable() %>%
  set_header_labels(tp_pairwise = "TP contrast",
                    roi_pairwise = "ROI contrast",
                    estimate = "Estimate",
                    `t.ratio` = "T",
                    `p.value` = "p") %>%
  set_caption(caption = "Interaction contrasts") %>%
  padding(padding.top = 0.5, padding.bottom = 0.5)
ft_interaction_contrasts <- set_formatter(ft_interaction_contrasts,
                                          estimate = function(x) format(round(x, 3), nsmall = 3),
                                          SE = function(x) format(round(x, 3), nsmall = 3),
                                          `t.ratio` = function(x) format(round(x, 3), nsmall = 3),
                                          `p.value` = function(x) format(round(x, 3), nsmall = 3))
ft_interaction_contrasts <- width(ft_interaction_contrasts, j = "tp_pairwise", width = 1)
ft_interaction_contrasts <- width(ft_interaction_contrasts, j = "roi_pairwise", width = 1.2)

contrast_summary <- summary(model.emm)
ROI <- c("Core", "Core", "Core", "Control", "Control", "Control")
`TP contrast` <- c("BL - 3mo", "BL - 6mo", "3mo - 6mo", "BL - 3mo", "BL - 6mo", "3mo - 6mo")
Estimate <- contrast_summary$contrasts$estimate
SE <- contrast_summary$contrasts$SE
df <- contrast_summary$contrasts$df
`T` <- contrast_summary$contrasts$t.ratio
`p` <- contrast_summary$contrasts$p.value
contrast_df <- data.frame(ROI, `TP contrast`, Estimate, SE, df, `T`, `p`)
ft_contrasts <- contrast_df %>%
  flextable() %>%
  set_header_labels(`TP.contrast` = "TP contrast") %>%
  set_caption(caption = "Simple contrasts") %>%
  padding(padding.top = 0.5, padding.bottom = 0.5)
ft_contrasts <- set_formatter(ft_contrasts,
                              Estimate = function(x) format(round(x, 3), nsmall = 3),
                              SE = function(x) format(round(x, 3), nsmall = 3),
                              `T` = function(x) format(round(x, 3), nsmall = 3),
                              `p` = function(x) format(round(x, 3), nsmall = 3))
ft_contrasts <- width(ft_contrasts, j = "TP.contrast", width = 1)

# Combine and display
ft1_file <- tempfile(fileext = ".html")
ft2_file <- tempfile(fileext = ".html")
ft3_file <- tempfile(fileext = ".html")
save_as_html(ft, path = ft1_file)
save_as_html(ft_interaction_contrasts, path = ft2_file)
save_as_html(ft_contrasts, path = ft3_file)

# Read and combine the HTML content
combined_display <- browsable(
  tagList(
    includeHTML(ft1_file),
    includeHTML(ft2_file),
    includeHTML(ft3_file)
  )
)
combined_display