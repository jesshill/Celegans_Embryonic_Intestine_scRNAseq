#### Script for graphing RNAi quantification data ####

getwd()

library(readr)
library(dplyr)
library(ggplot2)

#### Read in data ####

emb_data <- read_csv("./RNAi_Data/RNAi_Data_for_analysis_emb.csv")
View(emb_data)

L4_data <- read_csv("./RNAi_Data/RNAi_Data_for_analysis_L4.csv")
View(L4_data)


#### Clean data ####
clean_data_emb <- emb_data %>% filter(!grepl("background", image_number, ignore.case = TRUE), !is.na(noramlized_intestine_flourescence), noramlized_intestine_flourescence >= 0)
View(clean_data_emb)


clean_data_L4 <- L4_data %>% filter(!grepl("background", image_number, ignore.case = TRUE), !is.na(noramlized_intestine_flourescence), noramlized_intestine_flourescence >= 0) %>% filter(!(RNAi_strain == "c14c6.5" & image_number == "10")) %>% filter(!(RNAi_strain == "pop-1"))
View(clean_data_L4)


#### Calculate L4440 control means (within each replicate, when applicable) ####
control_means_emb <- clean_data_emb %>% filter(RNAi_strain == "L4440") %>% group_by(worm_strain, rep) %>% summarise(control_mean = mean(noramlized_intestine_flourescence, na.rm = TRUE), n_control = n(), .groups = "drop")
View(control_means_emb)


control_means_L4 <- clean_data_L4 %>% filter(RNAi_strain == "L4440") %>% group_by(worm_strain) %>% summarise(control_mean = mean(noramlized_intestine_flourescence, na.rm = TRUE), n_control = n(), .groups = "drop")
View(control_means_L4)

#### Normalize individual measurements ####
clean_data_normalized_emb <- clean_data_emb %>% left_join(control_means_emb, by = c("worm_strain", "rep")) %>% filter(!is.na(control_mean), is.finite(control_mean), control_mean > 0) %>% mutate(normalized_to_L4440 = noramlized_intestine_flourescence / control_mean)
View(clean_data_normalized_emb)

clean_data_normalized_L4 <- clean_data_L4 %>% left_join(control_means_L4, by = c("worm_strain")) %>% filter(!is.na(control_mean), is.finite(control_mean), control_mean > 0) %>% mutate(normalized_to_L4440 = noramlized_intestine_flourescence / control_mean)
View(clean_data_normalized_L4)

#### Combine replicates and calculate final summary statistics ####
normalized_group_means_emb <- clean_data_normalized_emb %>% group_by(worm_strain, RNAi_strain) %>% summarise(mean_normalized = mean(normalized_to_L4440, na.rm = TRUE), sd_normalized = sd(normalized_to_L4440, na.rm = TRUE), n = sum(!is.na(normalized_to_L4440)), sem_normalized = sd_normalized / sqrt(n), .groups = "drop")
View(normalized_group_means_emb)

normalized_group_means_L4 <- clean_data_normalized_L4 %>% group_by(worm_strain, RNAi_strain) %>% summarise(mean_normalized = mean(normalized_to_L4440, na.rm = TRUE), sd_normalized = sd(normalized_to_L4440, na.rm = TRUE), n = sum(!is.na(normalized_to_L4440)), sem_normalized = sd_normalized / sqrt(n), .groups = "drop")
View(normalized_group_means_L4)

#### Set RNAi strain order ####
strain_order <- c("L4440", "ugt-14", "cpr-1", "y32f6a.5", "c14c6.5", "endu-2", "clec-56", "pbo-4")

clean_data_normalized_emb <- clean_data_normalized_emb %>% mutate(RNAi_strain = factor(RNAi_strain, levels = strain_order))
View(clean_data_normalized_emb)

clean_data_normalized_L4 <- clean_data_normalized_L4 %>% mutate(RNAi_strain = factor(RNAi_strain, levels = strain_order))
View(clean_data_normalized_L4)

normalized_group_means_emb <- normalized_group_means_emb %>% mutate(RNAi_strain = factor(RNAi_strain, levels = strain_order))
View(normalized_group_means_emb)

normalized_group_means_L4 <- normalized_group_means_L4 %>% mutate(RNAi_strain = factor(RNAi_strain, levels = strain_order))
View(normalized_group_means_L4)

#### Visualize data ####

# pdf("./RNAi_data/JM149_RNAi_emb_data.pdf", height = 6, width = 9)
# ggplot(normalized_group_means_emb %>% filter(worm_strain == "JM149"), aes(x = RNAi_strain, y = mean_normalized)) + geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
#   geom_errorbar(aes(ymin = mean_normalized - sem_normalized, ymax = mean_normalized + sem_normalized), width = 0.2, linewidth = 0.7) +
#   geom_jitter(data = clean_data_normalized_emb %>% filter(worm_strain == "JM149"),
#               aes(
#                 x = RNAi_strain,
#                 y = normalized_to_L4440,
#                 color = factor(rep)
#               ),
#               width = 0.15,
#               size = 2,
#               alpha = 0.8,
#               inherit.aes = FALSE
#   ) +
#   scale_color_manual(
#     values = c("1" = "#9932CC", "2" = "#FFB90F"),
#     labels = c("1" = "Rep 1", "2" = "Rep 2"),
#     name = "Replicate"
#   ) +
#   geom_hline(
#     yintercept = 1,
#     linetype = "dashed",
#     color = "red",
#     linewidth = 1) +
#   theme_minimal() +
#   labs(
#     x = "RNAi strain",
#     y = "Normalized fluorescence (to L4440)",
#     title = "Normalized Intestinal Fluorescence (JM149)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# dev.off()


# pdf("./RNAi_data/JM149_RNAi_emb_boxplot_rep.pdf", height = 6, width = 9)
# ggplot(
#   clean_data_normalized_emb %>% filter(worm_strain == "JM149"),
#   aes(x = RNAi_strain, y = normalized_to_L4440)
# ) +
#   geom_boxplot(
#     aes(fill = RNAi_strain),
#     alpha = 0.5,
#     outlier.shape = NA
#   ) +
#   geom_jitter(
#     aes(color = factor(rep)),
#     width = 0.15,
#     size = 2,
#     alpha = 0.8
#   ) +
#   scale_color_manual(
#     values = c("1" = "#9932CC", "2" = "#FFB90F"),
#     labels = c("1" = "Rep 1", "2" = "Rep 2"),
#     name = "Replicate"
#   ) +
#   scale_fill_brewer(palette = "Set2", guide = "none") +
#   geom_hline(
#     yintercept = 1,
#     linetype = "dashed",
#     color = "red",
#     linewidth = 1
#   ) +
#   theme_minimal() +
#   labs(
#     x = "RNAi strain",
#     y = "Normalized fluorescence (to L4440)",
#     title = "Normalized Intestinal Fluorescence (JM149)"
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# dev.off()


pdf("./RNAi_data/JM149_RNAi_emb_boxplot.pdf", height = 6, width = 9)
ggplot(
  clean_data_normalized_emb %>% filter(worm_strain == "JM149"),
  aes(x = RNAi_strain, y = normalized_to_L4440)
) +
  geom_boxplot(
    aes(fill = RNAi_strain),
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter() +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  theme_minimal() +
  labs(
    x = "RNAi strain",
    y = "Normalized fluorescence (to L4440)",
    title = "Normalized Intestinal Fluorescence (JM149)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# pdf("./RNAi_data/MR142_RNAi_emb_data.pdf", height = 6, width = 9)
# ggplot(normalized_group_means_emb %>% filter(worm_strain == "MR142"), aes(x = RNAi_strain, y = mean_normalized)) + geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
#   geom_errorbar(aes(ymin = mean_normalized - sem_normalized, ymax = mean_normalized + sem_normalized), width = 0.2, linewidth = 0.7) + 
#   geom_jitter(data = clean_data_normalized_emb %>% filter(worm_strain == "MR142"),
#               aes(
#                 x = RNAi_strain,
#                 y = normalized_to_L4440,
#                 color = factor(rep)
#               ),
#               width = 0.15,
#               size = 2,
#               alpha = 0.8,
#               inherit.aes = FALSE
#   ) +
#   scale_color_manual(
#     values = c("1" = "#9932CC", "2" = "#FFB90F"),
#     labels = c("1" = "Rep 1", "2" = "Rep 2"),
#     name = "Replicate"
#   ) +
#   geom_hline(
#     yintercept = 1,
#     linetype = "dashed",
#     color = "red",
#     linewidth = 1) +
#   theme_minimal() +
#   labs(
#     x = "RNAi strain",
#     y = "Normalized fluorescence (to L4440)",
#     title = "Normalized Intestinal Fluorescence (MR142)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# dev.off()


# pdf("./RNAi_data/MR142_RNAi_emb_boxplot_rep.pdf", height = 6, width = 9)
# ggplot(
#   clean_data_normalized_emb %>% filter(worm_strain == "MR142"),
#   aes(x = RNAi_strain, y = normalized_to_L4440)
# ) +
#   geom_boxplot(
#     aes(fill = RNAi_strain),
#     alpha = 0.5,
#     outlier.shape = NA
#   ) +
#   geom_jitter(
#     aes(color = factor(rep)),
#     width = 0.15,
#     size = 2,
#     alpha = 0.8
#   ) +
#   scale_color_manual(
#     values = c("1" = "#9932CC", "2" = "#FFB90F"),
#     labels = c("1" = "Rep 1", "2" = "Rep 2"),
#     name = "Replicate"
#   ) +
#   scale_fill_brewer(palette = "Set2", guide = "none") +
#   geom_hline(
#     yintercept = 1,
#     linetype = "dashed",
#     color = "red",
#     linewidth = 1
#   ) +
#   theme_minimal() +
#   labs(
#     x = "RNAi strain",
#     y = "Normalized fluorescence (to L4440)",
#     title = "Normalized Intestinal Fluorescence (MR142)"
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# dev.off()


pdf("./RNAi_data/MR142_RNAi_emb_boxplot.pdf", height = 6, width = 9)
ggplot(
  clean_data_normalized_emb %>% filter(worm_strain == "MR142"),
  aes(x = RNAi_strain, y = normalized_to_L4440)
) +
  geom_boxplot(
    aes(fill = RNAi_strain),
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter() +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  theme_minimal() +
  labs(
    x = "RNAi strain",
    y = "Normalized fluorescence (to L4440)",
    title = "Normalized Intestinal Fluorescence (MR142)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()





pdf("./RNAi_data/JM149_RNAi_L4_boxplot.pdf", height = 6, width = 9)
ggplot(
  clean_data_normalized_L4 %>% filter(worm_strain == "JM149"),
  aes(x = RNAi_strain, y = normalized_to_L4440)
) +
  geom_boxplot(
    aes(fill = RNAi_strain),
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter() + 
  scale_fill_brewer(palette = "Greys", guide = "none") +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  theme_minimal() +
  labs(
    x = "RNAi strain",
    y = "Normalized fluorescence (to L4440)",
    title = "Normalized Intestinal Fluorescence (JM149)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


pdf("./RNAi_data/MR142_RNAi_L4_boxplot.pdf", height = 6, width = 9)
ggplot(
  clean_data_normalized_L4 %>% filter(worm_strain == "MR142"),
  aes(x = RNAi_strain, y = normalized_to_L4440)
) +
  geom_boxplot(
    aes(fill = RNAi_strain),
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter() + 
  scale_fill_brewer(palette = "Greys", guide = "none") +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  theme_minimal() +
  labs(
    x = "RNAi strain",
    y = "Normalized fluorescence (to L4440)",
    title = "Normalized Intestinal Fluorescence (MR142)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


pdf("./RNAi_data/ERT60_RNAi_L4_boxplot.pdf", height = 6, width = 9)
ggplot(
  clean_data_normalized_L4 %>% filter(worm_strain == "ERT60"),
  aes(x = RNAi_strain, y = normalized_to_L4440)
) +
  geom_boxplot(
    aes(fill = RNAi_strain),
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter() + 
  scale_fill_brewer(palette = "Greys", guide = "none") +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  theme_minimal() +
  labs(
    x = "RNAi strain",
    y = "Normalized fluorescence (to L4440)",
    title = "Normalized Intestinal Fluorescence (ERT60)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



#### STATISTICS - Welch's t-test ####

#### extract data points ####
data_for_Ttest_emb <- clean_data_normalized_emb[, c("worm_strain", "RNAi_strain", "rep", "noramlized_intestine_flourescence")]
View(data_for_Ttest_emb)

data_for_Ttest_L4 <- clean_data_normalized_L4[, c("worm_strain", "RNAi_strain", "noramlized_intestine_flourescence")]
View(data_for_Ttest_L4)


#### check normality ####
install.packages("rstatix")
library(rstatix)

# Shapiro test - Null Hypothesis (H₀): The data is normally distributed vs. Alternative Hypothesis (H₁): The data is not normally distributed
normality_L4 <- data_for_Ttest_L4 %>% group_by(worm_strain, RNAi_strain) %>% shapiro_test(noramlized_intestine_flourescence)
View(normality_L4)


normality_emb <- data_for_Ttest_emb %>% group_by(worm_strain, RNAi_strain) %>% shapiro_test(noramlized_intestine_flourescence)
View(normality_emb)


# Levene Test - Null Hypothesis (H₀): The variances across the groups are equal (homogeneity of variance holds) vs. Alternative Hypothesis (H₁): At least one group has a significantly different variance
variance_L4 <- data_for_Ttest_L4 %>% group_by(worm_strain) %>% levene_test(noramlized_intestine_flourescence ~ RNAi_strain)
View(variance_L4)

variance_emb <- data_for_Ttest_emb %>% group_by(worm_strain) %>% levene_test(noramlized_intestine_flourescence ~ RNAi_strain)
View(variance_emb)


# Q-Q plot
qqnorm(data_for_Ttest_emb$noramlized_intestine_flourescence, main = "Normal Q-Q Plot", col = "blue", pch = 19)
qqline(data_for_Ttest_emb$noramlized_intestine_flourescence, col = "red", lwd = 2)

qqnorm(data_for_Ttest_L4$noramlized_intestine_flourescence, main = "Normal Q-Q Plot", col = "blue", pch = 19)
qqline(data_for_Ttest_L4$noramlized_intestine_flourescence, col = "red", lwd = 2)


# compare two groups
qqplot(L4440_JM149_L4$noramlized_intestine_flourescence, y32f6a.5_JM149_L4$noramlized_intestine_flourescence, xlab = "L4440", ylab = "y32f6a.5") 



# Q-Q plots for each group - display multiple
JM149_emb <- data_for_Ttest_emb[data_for_Ttest_emb$worm_strain == "JM149", ]
MR142_emb <- data_for_Ttest_emb[data_for_Ttest_emb$worm_strain == "MR142", ]

JM149_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$worm_strain == "JM149", ]
MR142_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$worm_strain == "MR142", ]
ERT60_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$worm_strain == "ERT60", ]



RNAi_list_JM149 <- unique(JM149_emb$RNAi_strain)
n <- length(RNAi_list_JM149)
# Choose 3 columns of plots
ncol_plot <- 3
nrow_plot <- ceiling(n / ncol_plot)

pdf("./RNAi_data/QQplots_emb_JM149.pdf", height = 10, width = 12)
# Arrange all plots on one page
par(mfrow = c(nrow_plot, ncol_plot),
  mar = c(3, 3, 3, 1),
  oma = c(0, 0, 2, 0))

for (RNAi in RNAi_list_JM149) {
  group_data <- JM149_emb %>% filter(RNAi_strain == RNAi)
  
  x <- group_data$noramlized_intestine_flourescence
  x <- x[is.finite(x)]
  
  qqnorm(
    x,
    main = as.character(RNAi),
    col = "blue",
    pch = 19,
    cex = 0.6
  )
  
  qqline(x, col = "red", lwd = 2)
}
mtext("JM149 Embryos: Normal Q-Q Plots",
      outer = TRUE, cex = 1.2)
dev.off()



RNAi_list_MR142 <- unique(MR142_emb$RNAi_strain)
n <- length(RNAi_list_MR142)
# Choose 3 columns of plots
ncol_plot <- 3
nrow_plot <- ceiling(n / ncol_plot)

pdf("./RNAi_data/QQplots_emb_MR142.pdf", height = 10, width = 12)
# Arrange all plots on one page
par(mfrow = c(nrow_plot, ncol_plot),
    mar = c(3, 3, 3, 1),
    oma = c(0, 0, 2, 0))

for (RNAi in RNAi_list_MR142) {
  group_data <- MR142_emb %>% filter(RNAi_strain == RNAi)
  
  x <- group_data$noramlized_intestine_flourescence
  x <- x[is.finite(x)]
  
  qqnorm(
    x,
    main = as.character(RNAi),
    col = "blue",
    pch = 19,
    cex = 0.6
  )
  
  qqline(x, col = "red", lwd = 2)
}
mtext("MR142 Embryos: Normal Q-Q Plots",
      outer = TRUE, cex = 1.2)
dev.off()




RNAi_list_JM149_L4 <- unique(JM149_L4$RNAi_strain)
n <- length(RNAi_list_JM149_L4)
# Choose 3 columns of plots
ncol_plot <- 3
nrow_plot <- ceiling(n / ncol_plot)

pdf("./RNAi_data/QQplots_L4_JM149.pdf", height = 10, width = 12)
# Arrange all plots on one page
par(mfrow = c(nrow_plot, ncol_plot),
    mar = c(3, 3, 3, 1),
    oma = c(0, 0, 2, 0))

for (RNAi in RNAi_list_JM149_L4) {
  group_data <- JM149_L4 %>% filter(RNAi_strain == RNAi)
  
  x <- group_data$noramlized_intestine_flourescence
  x <- x[is.finite(x)]
  
  qqnorm(
    x,
    main = as.character(RNAi),
    col = "blue",
    pch = 19,
    cex = 0.6
  )
  
  qqline(x, col = "red", lwd = 2)
}
mtext("JM149 L4s: Normal Q-Q Plots",
      outer = TRUE, cex = 1.2)
dev.off()



RNAi_list_MR142_L4 <- unique(MR142_L4$RNAi_strain)
n <- length(RNAi_list_MR142_L4)
# Choose 3 columns of plots
ncol_plot <- 3
nrow_plot <- ceiling(n / ncol_plot)

pdf("./RNAi_data/QQplots_L4_MR142.pdf", height = 10, width = 12)
# Arrange all plots on one page
par(mfrow = c(nrow_plot, ncol_plot),
    mar = c(3, 3, 3, 1),
    oma = c(0, 0, 2, 0))

for (RNAi in RNAi_list_MR142_L4) {
  group_data <- MR142_L4 %>% filter(RNAi_strain == RNAi)
  
  x <- group_data$noramlized_intestine_flourescence
  x <- x[is.finite(x)]
  
  qqnorm(
    x,
    main = as.character(RNAi),
    col = "blue",
    pch = 19,
    cex = 0.6
  )
  
  qqline(x, col = "red", lwd = 2)
}
mtext("MR142 L4s: Normal Q-Q Plots",
      outer = TRUE, cex = 1.2)
dev.off()



RNAi_list_ERT60_L4 <- unique(ERT60_L4$RNAi_strain)
n <- length(RNAi_list_ERT60_L4)
# Choose 3 columns of plots
ncol_plot <- 3
nrow_plot <- ceiling(n / ncol_plot)

pdf("./RNAi_data/QQplots_L4_ERT60.pdf", height = 10, width = 12)
# Arrange all plots on one page
par(mfrow = c(nrow_plot, ncol_plot),
    mar = c(3, 3, 3, 1),
    oma = c(0, 0, 2, 0))

for (RNAi in RNAi_list_ERT60_L4) {
  group_data <- ERT60_L4 %>% filter(RNAi_strain == RNAi)
  
  x <- group_data$noramlized_intestine_flourescence
  x <- x[is.finite(x)]
  
  qqnorm(
    x,
    main = as.character(RNAi),
    col = "blue",
    pch = 19,
    cex = 0.6
  )
  
  qqline(x, col = "red", lwd = 2)
}
mtext("ERT60 L4s: Normal Q-Q Plots",
      outer = TRUE, cex = 1.2)
dev.off()


#### Separate data by RNAi strain  ####
L4440_JM149_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "L4440" & data_for_Ttest_L4$worm_strain == "JM149", ]
L4440_ERT60_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "L4440" & data_for_Ttest_L4$worm_strain == "ERT60", ]
L4440_MR142_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "L4440" & data_for_Ttest_L4$worm_strain == "MR142", ]

y32f6a.5_JM149_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "y32f6a.5" & data_for_Ttest_L4$worm_strain == "JM149", ]
y32f6a.5_ERT60_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "y32f6a.5" & data_for_Ttest_L4$worm_strain == "ERT60", ]
y32f6a.5_MR142_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "y32f6a.5" & data_for_Ttest_L4$worm_strain == "MR142", ]

cpr.1_JM149_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "cpr-1" & data_for_Ttest_L4$worm_strain == "JM149", ]
cpr.1_ERT60_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "cpr-1" & data_for_Ttest_L4$worm_strain == "ERT60", ]
cpr.1_MR142_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "cpr-1" & data_for_Ttest_L4$worm_strain == "MR142", ]

ugt.14_JM149_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "ugt-14" & data_for_Ttest_L4$worm_strain == "JM149", ]
ugt.14_ERT60_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "ugt-14" & data_for_Ttest_L4$worm_strain == "ERT60", ]
ugt.14_MR142_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "ugt-14" & data_for_Ttest_L4$worm_strain == "MR142", ]

endu.2_JM149_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "endu-2" & data_for_Ttest_L4$worm_strain == "JM149", ]
endu.2_ERT60_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "endu-2" & data_for_Ttest_L4$worm_strain == "ERT60", ]
endu.2_MR142_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "endu-2" & data_for_Ttest_L4$worm_strain == "MR142", ]

clec.56_JM149_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "clec-56" & data_for_Ttest_L4$worm_strain == "JM149", ]
clec.56_ERT60_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "clec-56" & data_for_Ttest_L4$worm_strain == "ERT60", ]
clec.56_MR142_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "clec-56" & data_for_Ttest_L4$worm_strain == "MR142", ]

c14c6.5_JM149_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "c14c6.5" & data_for_Ttest_L4$worm_strain == "JM149", ]
c14c6.5_ERT60_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "c14c6.5" & data_for_Ttest_L4$worm_strain == "ERT60", ]
c14c6.5_MR142_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "c14c6.5" & data_for_Ttest_L4$worm_strain == "MR142", ]

pbo.4_JM149_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "pbo-4" & data_for_Ttest_L4$worm_strain == "JM149", ]
pbo.4_ERT60_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "pbo-4" & data_for_Ttest_L4$worm_strain == "ERT60", ]
pbo.4_MR142_L4 <- data_for_Ttest_L4[data_for_Ttest_L4$RNAi_strain == "pbo-4" & data_for_Ttest_L4$worm_strain == "MR142", ]




L4440_JM149_emb <- data_for_Ttest_emb[data_for_Ttest_emb$RNAi_strain == "L4440" & data_for_Ttest_emb$worm_strain == "JM149", ]
L4440_MR142_emb <- data_for_Ttest_emb[data_for_Ttest_emb$RNAi_strain == "L4440" & data_for_Ttest_emb$worm_strain == "MR142", ]

y32f6a.5_JM149_emb <- data_for_Ttest_emb[data_for_Ttest_emb$RNAi_strain == "y32f6a.5" & data_for_Ttest_emb$worm_strain == "JM149", ]
y32f6a.5_MR142_emb <- data_for_Ttest_emb[data_for_Ttest_emb$RNAi_strain == "y32f6a.5" & data_for_Ttest_emb$worm_strain == "MR142", ]

cpr.1_JM149_emb <- data_for_Ttest_emb[data_for_Ttest_emb$RNAi_strain == "cpr-1" & data_for_Ttest_emb$worm_strain == "JM149", ]
cpr.1_MR142_emb <- data_for_Ttest_emb[data_for_Ttest_emb$RNAi_strain == "cpr-1" & data_for_Ttest_emb$worm_strain == "MR142", ]

ugt.14_JM149_emb <- data_for_Ttest_emb[data_for_Ttest_emb$RNAi_strain == "ugt-14" & data_for_Ttest_emb$worm_strain == "JM149", ]
ugt.14_MR142_emb <- data_for_Ttest_emb[data_for_Ttest_emb$RNAi_strain == "ugt-14" & data_for_Ttest_emb$worm_strain == "MR142", ]

endu.2_JM149_emb <- data_for_Ttest_emb[data_for_Ttest_emb$RNAi_strain == "endu-2" & data_for_Ttest_emb$worm_strain == "JM149", ]
endu.2_MR142_emb <- data_for_Ttest_emb[data_for_Ttest_emb$RNAi_strain == "endu-2" & data_for_Ttest_emb$worm_strain == "MR142", ]

clec.56_JM149_emb <- data_for_Ttest_emb[data_for_Ttest_emb$RNAi_strain == "clec-56" & data_for_Ttest_emb$worm_strain == "JM149", ]
clec.56_MR142_emb <- data_for_Ttest_emb[data_for_Ttest_emb$RNAi_strain == "clec-56" & data_for_Ttest_emb$worm_strain == "MR142", ]

c14c6.5_JM149_emb <- data_for_Ttest_emb[data_for_Ttest_emb$RNAi_strain == "c14c6.5" & data_for_Ttest_emb$worm_strain == "JM149", ]
c14c6.5_MR142_emb <- data_for_Ttest_emb[data_for_Ttest_emb$RNAi_strain == "c14c6.5" & data_for_Ttest_emb$worm_strain == "MR142", ]

pbo.4_JM149_emb <- data_for_Ttest_emb[data_for_Ttest_emb$RNAi_strain == "pbo-4" & data_for_Ttest_emb$worm_strain == "JM149", ]
pbo.4_MR142_emb <- data_for_Ttest_emb[data_for_Ttest_emb$RNAi_strain == "pbo-4" & data_for_Ttest_emb$worm_strain == "MR142", ]



#### Set up Welch's t-test  ####

set.seed(42)

## for y32f6a-5 ##
t.test(L4440_JM149_L4$noramlized_intestine_flourescence, y32f6a.5_JM149_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
t.test(L4440_ERT60_L4$noramlized_intestine_flourescence, y32f6a.5_ERT60_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95) 
t.test(L4440_MR142_L4$noramlized_intestine_flourescence, y32f6a.5_MR142_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95) ##significant


t.test(L4440_JM149_emb$noramlized_intestine_flourescence, y32f6a.5_JM149_emb$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
t.test(L4440_MR142_emb$noramlized_intestine_flourescence, y32f6a.5_MR142_emb$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)


## for cpr-1 ##
t.test(L4440_JM149_L4$noramlized_intestine_flourescence, cpr.1_JM149_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
t.test(L4440_ERT60_L4$noramlized_intestine_flourescence, cpr.1_ERT60_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95) ##significant
t.test(L4440_MR142_L4$noramlized_intestine_flourescence, cpr.1_MR142_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95) 


t.test(L4440_JM149_emb$noramlized_intestine_flourescence, cpr.1_JM149_emb$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95) ##significant
t.test(L4440_MR142_emb$noramlized_intestine_flourescence, cpr.1_MR142_emb$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95) ##significant


## for ugt-14 ##
t.test(L4440_JM149_L4$noramlized_intestine_flourescence, ugt.14_JM149_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
t.test(L4440_ERT60_L4$noramlized_intestine_flourescence, ugt.14_ERT60_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
t.test(L4440_MR142_L4$noramlized_intestine_flourescence, ugt.14_MR142_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)


t.test(L4440_JM149_emb$noramlized_intestine_flourescence, ugt.14_JM149_emb$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
t.test(L4440_MR142_emb$noramlized_intestine_flourescence, ugt.14_MR142_emb$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)


## for endu-2 ##
t.test(L4440_JM149_L4$noramlized_intestine_flourescence, endu.2_JM149_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
t.test(L4440_ERT60_L4$noramlized_intestine_flourescence, endu.2_ERT60_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
t.test(L4440_MR142_L4$noramlized_intestine_flourescence, endu.2_MR142_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)


t.test(L4440_JM149_emb$noramlized_intestine_flourescence, endu.2_JM149_emb$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95) ## significant
t.test(L4440_MR142_emb$noramlized_intestine_flourescence, endu.2_MR142_emb$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95) ## significant


## for clec-56 ##
t.test(L4440_JM149_L4$noramlized_intestine_flourescence, clec.56_JM149_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
t.test(L4440_ERT60_L4$noramlized_intestine_flourescence, clec.56_ERT60_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
t.test(L4440_MR142_L4$noramlized_intestine_flourescence, clec.56_MR142_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)


t.test(L4440_JM149_emb$noramlized_intestine_flourescence, clec.56_JM149_emb$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
t.test(L4440_MR142_emb$noramlized_intestine_flourescence, clec.56_MR142_emb$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95) ## significant


## for c14c6.5 ##
t.test(L4440_JM149_L4$noramlized_intestine_flourescence, c14c6.5_JM149_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
t.test(L4440_ERT60_L4$noramlized_intestine_flourescence, c14c6.5_ERT60_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
t.test(L4440_MR142_L4$noramlized_intestine_flourescence, c14c6.5_MR142_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)


t.test(L4440_JM149_emb$noramlized_intestine_flourescence, c14c6.5_JM149_emb$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
t.test(L4440_MR142_emb$noramlized_intestine_flourescence, c14c6.5_MR142_emb$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95) ## significant


## for pbo-4 ##
t.test(L4440_JM149_L4$noramlized_intestine_flourescence, pbo.4_JM149_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
t.test(L4440_ERT60_L4$noramlized_intestine_flourescence, pbo.4_ERT60_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
t.test(L4440_MR142_L4$noramlized_intestine_flourescence, pbo.4_MR142_L4$noramlized_intestine_flourescence, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)


t.test(L4440_JM149_emb$noramlized_intestine_flourescence, pbo.4_JM149_emb$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(L4440_MR142_emb$noramlized_intestine_flourescence, pbo.4_MR142_emb$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95) ## significant
