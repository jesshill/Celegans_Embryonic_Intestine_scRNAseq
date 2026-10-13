### Script for graphing RNAi quantification data ###

getwd()

library(readr)
library(dplyr)
library(ggplot2)


### read in data ###
#data <- read_csv("./RNAi_Data/RNAi_Data_for_analysis_L4.csv")
data <- read_csv("./RNAi_Data/RNAi_Data_for_analysis_emb.csv")
View(data)


### clean data by removing background image readings ###
clean_data <- data %>% filter(!grepl("background", image_number)) 
View(clean_data)

  
### remove any outliers ###

#L4 data only
clean_data_no_outlier <- clean_data %>% filter(!(RNAi_strain == "c14c6.5" & image_number == "10"))
clean_data_no_outlier <- clean_data_no_outlier %>% filter(!(RNAi_strain == "pop-1"))
View(clean_data_no_outlier)

#emb data only
clean_data_no_outlier <- clean_data %>% filter(noramlized_intestine_flourescence >= 0)
View(clean_data_no_outlier)


### calc group means ###

# L4 data only
group_means <- clean_data_no_outlier %>% group_by(worm_strain, RNAi_strain) %>% summarise(
    mean_fluor = mean(noramlized_intestine_flourescence, na.rm = TRUE),
    sd_fluor = sd(noramlized_intestine_flourescence, na.rm = TRUE),
    n = n(),
    sem_fluor = sd_fluor / sqrt(n),
    .groups = "drop")
View(group_means)

# emb data only
group_means <- clean_data_no_outlier %>% group_by(worm_strain, RNAi_strain, rep) %>% summarise(
  mean_fluor = mean(noramlized_intestine_flourescence, na.rm = TRUE),
  sd_fluor = sd(noramlized_intestine_flourescence, na.rm = TRUE),
  n = n(),
  sem_fluor = sd_fluor / sqrt(n),
  .groups = "drop")
View(group_means)


### calc negative control (L4440) means ###

# L4 data only
control_means <- group_means %>% filter(RNAi_strain == "L4440") %>% dplyr::select(worm_strain, control_mean = mean_fluor)
View(control_means)

# emb data only
control_means_rep1 <- group_means %>% filter(rep == "1" & RNAi_strain == "L4440") %>% dplyr::select(worm_strain, control_mean = mean_fluor)
View(control_means_rep1)

control_means_rep2 <- group_means %>% filter(rep == "2" & RNAi_strain == "L4440") %>% dplyr::select(worm_strain, control_mean = mean_fluor)
View(control_means_rep2)


### normalize individual observations ###

# L4 data only
clean_data_normalized <- clean_data_no_outlier %>% left_join(control_means, by = "worm_strain") %>% mutate(normalized_to_L4440 = noramlized_intestine_flourescence / control_mean)
View(clean_data_normalized)

# emb data only
clean_data_normalized_rep1 <- clean_data_no_outlier %>% filter(rep == "1") %>% left_join(control_means_rep1, by = "worm_strain") %>% mutate(normalized_to_L4440 = noramlized_intestine_flourescence / control_mean)
View(clean_data_normalized_rep1)

clean_data_normalized_rep2 <- clean_data_no_outlier %>% filter(rep == "2") %>% left_join(control_means_rep2, by = "worm_strain") %>% mutate(normalized_to_L4440 = noramlized_intestine_flourescence / control_mean)
View(clean_data_normalized_rep2)


### summarize normalized replicates ###

# L4 data only
normalized_group_means <- clean_data_normalized %>% group_by(worm_strain, RNAi_strain) %>%
  summarise(mean_normalized = mean(normalized_to_L4440, na.rm = TRUE), sd_normalized = sd(normalized_to_L4440, na.rm = TRUE), n = n(), sem_normalized = sd_normalized / sqrt(n), .groups = "drop")
View(normalized_group_means)


# emb data only
normalized_group_means_rep1 <- clean_data_normalized_rep1 %>% group_by(worm_strain, RNAi_strain) %>%
  summarise(mean_normalized = mean(normalized_to_L4440, na.rm = TRUE), sd_normalized = sd(normalized_to_L4440, na.rm = TRUE), n = n(), sem_normalized = sd_normalized / sqrt(n), .groups = "drop")
View(normalized_group_means_rep1)

normalized_group_means_rep2 <- clean_data_normalized_rep2 %>% group_by(worm_strain, RNAi_strain) %>%
  summarise(mean_normalized = mean(normalized_to_L4440, na.rm = TRUE), sd_normalized = sd(normalized_to_L4440, na.rm = TRUE), n = n(), sem_normalized = sd_normalized / sqrt(n), .groups = "drop")
View(normalized_group_means_rep2)


### set factor order ###

strain_order <- c("L4440", "ugt-14", "cpr-1", "y32f6a.5", "c14c6.5", "endu-2", "clec-56", "pbo-4")

# L4 data only
normalized_group_means <- normalized_group_means %>% mutate(RNAi_strain = factor(RNAi_strain, levels = strain_order))
clean_data_normalized <- clean_data_normalized %>% mutate(RNAi_strain = factor(RNAi_strain, levels = strain_order))

# emb data only
normalized_group_means_rep1 <- normalized_group_means_rep1 %>% mutate(RNAi_strain = factor(RNAi_strain, levels = strain_order))
clean_data_normalized_rep1 <- clean_data_normalized_rep1 %>% mutate(RNAi_strain = factor(RNAi_strain, levels = strain_order))

normalized_group_means_rep2 <- normalized_group_means_rep2 %>% mutate(RNAi_strain = factor(RNAi_strain, levels = strain_order))
clean_data_normalized_rep2 <- clean_data_normalized_rep2 %>% mutate(RNAi_strain = factor(RNAi_strain, levels = strain_order))


### graph the data - for embryos ###

pdf("./RNAi_data/JM149_RNAi_emb_data_rep1.pdf", height = 6, width = 9)
ggplot(normalized_group_means_rep1 %>% filter(worm_strain == "JM149"), aes(x = RNAi_strain, y = mean_normalized)) + geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_normalized - sem_normalized, ymax = mean_normalized + sem_normalized), width = 0.2, linewidth = 0.7) + geom_jitter(
    data = clean_data_normalized_rep1 %>% filter(worm_strain == "JM149"),
    aes(x = RNAi_strain, y = normalized_to_L4440),
    width = 0.15,
    size = 2,
    alpha = 0.8,
    color = "black",
    inherit.aes = FALSE) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red",
    linewidth = 1) +
  theme_minimal() +
  labs(
    x = "RNAi strain",
    y = "Normalized fluorescence (to L4440)",
    title = "Normalized Intestinal Fluorescence (JM149)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("./RNAi_data/JM149_RNAi_emb_data_rep2.pdf", height = 6, width = 9)
ggplot(normalized_group_means_rep2 %>% filter(worm_strain == "JM149"), aes(x = RNAi_strain, y = mean_normalized)) + geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_normalized - sem_normalized, ymax = mean_normalized + sem_normalized), width = 0.2, linewidth = 0.7) + geom_jitter(
    data = clean_data_normalized_rep2 %>% filter(worm_strain == "JM149"),
    aes(x = RNAi_strain, y = normalized_to_L4440),
    width = 0.15,
    size = 2,
    alpha = 0.8,
    color = "black",
    inherit.aes = FALSE) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red",
    linewidth = 1) +
  theme_minimal() +
  labs(
    x = "RNAi strain",
    y = "Normalized fluorescence (to L4440)",
    title = "Normalized Intestinal Fluorescence (JM149)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("./RNAi_data/MR142_RNAi_emb_data_rep1.pdf", height = 6, width = 9)
ggplot(normalized_group_means_rep1 %>% filter(worm_strain == "MR142"), aes(x = RNAi_strain, y = mean_normalized)) + geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_normalized - sem_normalized, ymax = mean_normalized + sem_normalized), width = 0.2, linewidth = 0.7) + geom_jitter(
    data = clean_data_normalized_rep1 %>% filter(worm_strain == "MR142"),
    aes(x = RNAi_strain, y = normalized_to_L4440),
    width = 0.15,
    size = 2,
    alpha = 0.8,
    color = "black",
    inherit.aes = FALSE) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red",
    linewidth = 1) +
  theme_minimal() +
  labs(
    x = "RNAi strain",
    y = "Normalized fluorescence (to L4440)",
    title = "Normalized Intestinal Fluorescence (MR142)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("./RNAi_data/MR142_RNAi_emb_data_rep2.pdf", height = 6, width = 9)
ggplot(normalized_group_means_rep2 %>% filter(worm_strain == "MR142"), aes(x = RNAi_strain, y = mean_normalized)) + geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_normalized - sem_normalized, ymax = mean_normalized + sem_normalized), width = 0.2, linewidth = 0.7) + geom_jitter(
    data = clean_data_normalized_rep2 %>% filter(worm_strain == "MR142"),
    aes(x = RNAi_strain, y = normalized_to_L4440),
    width = 0.15,
    size = 2,
    alpha = 0.8,
    color = "black",
    inherit.aes = FALSE) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red",
    linewidth = 1) +
  theme_minimal() +
  labs(
    x = "RNAi strain",
    y = "Normalized fluorescence (to L4440)",
    title = "Normalized Intestinal Fluorescence (MR142)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


### graph the data - for L4s ###

pdf("./RNAi_data/JM149_RNAi_L4_data.pdf", height = 6, width = 9)
ggplot(normalized_group_means %>% filter(worm_strain == "JM149"), aes(x = RNAi_strain, y = mean_normalized)) + geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_normalized - sem_normalized, ymax = mean_normalized + sem_normalized), width = 0.2, linewidth = 0.7) + geom_jitter(
    data = clean_data_normalized %>% filter(worm_strain == "JM149"),
    aes(x = RNAi_strain, y = normalized_to_L4440),
    width = 0.15,
    size = 2,
    alpha = 0.8,
    color = "black",
    inherit.aes = FALSE) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red",
    linewidth = 1) +
  theme_minimal() +
  labs(
    x = "RNAi strain",
    y = "Normalized fluorescence (to L4440)",
    title = "Normalized Intestinal Fluorescence (JM149)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


pdf("./RNAi_data/ERT60_RNAi_L4_data.pdf", height = 6, width = 9)
ggplot(normalized_group_means %>% filter(worm_strain == "ERT60"), aes(x = RNAi_strain, y = mean_normalized)) + geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_normalized - sem_normalized, ymax = mean_normalized + sem_normalized), width = 0.2, linewidth = 0.7) + geom_jitter(
    data = clean_data_normalized %>% filter(worm_strain == "ERT60"),
    aes(x = RNAi_strain, y = normalized_to_L4440),
    width = 0.15,
    size = 2,
    alpha = 0.8,
    color = "black",
    inherit.aes = FALSE) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red",
    linewidth = 1) +
  theme_minimal() +
  labs(
    x = "RNAi strain",
    y = "Normalized fluorescence (to L4440)",
    title = "Normalized Intestinal Fluorescence (ERT60)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


pdf("./RNAi_data/MR142_RNAi_L4_data.pdf", height = 6, width = 9)
ggplot(normalized_group_means %>% filter(worm_strain == "MR142"), aes(x = RNAi_strain, y = mean_normalized)) + geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_normalized - sem_normalized, ymax = mean_normalized + sem_normalized), width = 0.2, linewidth = 0.7) + geom_jitter(
    data = clean_data_normalized %>% filter(worm_strain == "MR142"),
    aes(x = RNAi_strain, y = normalized_to_L4440),
    width = 0.15,
    size = 2,
    alpha = 0.8,
    color = "black",
    inherit.aes = FALSE) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red",
    linewidth = 1) +
  theme_minimal() +
  labs(
    x = "RNAi strain",
    y = "Normalized fluorescence (to L4440)",
    title = "Normalized Intestinal Fluorescence (MR142)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


### students t test ###
data_for_Ttest <- clean_data_normalized[, c("worm_strain", "RNAi_strain", "noramlized_intestine_flourescence")]
view(data_for_Ttest)


pop.1_JM149 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "pop-1" &
    data_for_Ttest$worm_strain == "JM149", ]
pop.1_ERT60 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "pop-1" &
    data_for_Ttest$worm_strain == "ERT60", ]
pop.1_MR142 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "pop-1" &
    data_for_Ttest$worm_strain == "MR142", ]

L4440_JM149 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "L4440" &
    data_for_Ttest$worm_strain == "JM149", ]
L4440_ERT60 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "L4440" &
    data_for_Ttest$worm_strain == "ERT60", ]
L4440_MR142 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "L4440" &
    data_for_Ttest$worm_strain == "MR142", ]

y32f6a.5_JM149 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "y32f6a.5" &
    data_for_Ttest$worm_strain == "JM149", ]
y32f6a.5_ERT60 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "y32f6a.5" &
    data_for_Ttest$worm_strain == "ERT60", ]
y32f6a.5_MR142 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "y32f6a.5" &
    data_for_Ttest$worm_strain == "MR142", ]

cpr.1_JM149 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "cpr-1" &
    data_for_Ttest$worm_strain == "JM149", ]
cpr.1_ERT60 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "cpr-1" &
    data_for_Ttest$worm_strain == "ERT60", ]
cpr.1_MR142 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "cpr-1" &
    data_for_Ttest$worm_strain == "MR142", ]

ugt.14_JM149 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "ugt-14" &
    data_for_Ttest$worm_strain == "JM149", ]
ugt.14_ERT60 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "ugt-14" &
    data_for_Ttest$worm_strain == "ERT60", ]
ugt.14_MR142 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "ugt-14" &
    data_for_Ttest$worm_strain == "MR142", ]

endu.2_JM149 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "endu-2" &
    data_for_Ttest$worm_strain == "JM149", ]
endu.2_ERT60 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "endu-2" &
    data_for_Ttest$worm_strain == "ERT60", ]
endu.2_MR142 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "endu-2" &
    data_for_Ttest$worm_strain == "MR142", ]

clec.56_JM149 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "clec-56" &
    data_for_Ttest$worm_strain == "JM149", ]
clec.56_ERT60 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "clec-56" &
    data_for_Ttest$worm_strain == "ERT60", ]
clec.56_MR142 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "clec-56" &
    data_for_Ttest$worm_strain == "MR142", ]

c14c6.5_JM149 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "c14c6.5" &
    data_for_Ttest$worm_strain == "JM149", ]
c14c6.5_ERT60 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "c14c6.5" &
    data_for_Ttest$worm_strain == "ERT60", ]
c14c6.5_MR142 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "c14c6.5" &
    data_for_Ttest$worm_strain == "MR142", ]

pbo.4_JM149 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "pbo-4" &
    data_for_Ttest$worm_strain == "JM149", ]
pbo.4_ERT60 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "pbo-4" &
    data_for_Ttest$worm_strain == "ERT60", ]
pbo.4_MR142 <- data_for_Ttest[
  data_for_Ttest$RNAi_strain == "pbo-4" &
    data_for_Ttest$worm_strain == "MR142", ]

# for pop-1
t.test(L4440_JM149$noramlized_intestine_flourescence, pop.1_JM149$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(L4440_ERT60$noramlized_intestine_flourescence, pop.1_ERT60$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(L4440_MR142$noramlized_intestine_flourescence, pop.1_MR142$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

wilcox.test(L4440_JM149$noramlized_intestine_flourescence, pop.1_JM149$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)
wilcox.test(L4440_ERT60$noramlized_intestine_flourescence, pop.1_ERT60$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)
wilcox.test(L4440_MR142$noramlized_intestine_flourescence, pop.1_MR142$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)


# for y32f6a-5
t.test(L4440_JM149$noramlized_intestine_flourescence, y32f6a.5_JM149$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(L4440_ERT60$noramlized_intestine_flourescence, y32f6a.5_ERT60$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(L4440_MR142$noramlized_intestine_flourescence, y32f6a.5_MR142$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95) ##significant

wilcox.test(L4440_JM149$noramlized_intestine_flourescence, y32f6a.5_JM149$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)
wilcox.test(L4440_ERT60$noramlized_intestine_flourescence, y32f6a.5_ERT60$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)
wilcox.test(L4440_MR142$noramlized_intestine_flourescence, y32f6a.5_MR142$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)  ##significant


# for cpr-1
t.test(L4440_JM149$noramlized_intestine_flourescence, cpr.1_JM149$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(L4440_ERT60$noramlized_intestine_flourescence, cpr.1_ERT60$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95) ##significant
t.test(L4440_MR142$noramlized_intestine_flourescence, cpr.1_MR142$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95) 

wilcox.test(L4440_JM149$noramlized_intestine_flourescence, cpr.1_JM149$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)
wilcox.test(L4440_ERT60$noramlized_intestine_flourescence, cpr.1_ERT60$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95) ##significant
wilcox.test(L4440_MR142$noramlized_intestine_flourescence, cpr.1_MR142$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95) 


# for ugt-14
t.test(L4440_JM149$noramlized_intestine_flourescence, ugt.14_JM149$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(L4440_ERT60$noramlized_intestine_flourescence, ugt.14_ERT60$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(L4440_MR142$noramlized_intestine_flourescence, ugt.14_MR142$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

wilcox.test(L4440_JM149$noramlized_intestine_flourescence, ugt.14_JM149$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)
wilcox.test(L4440_ERT60$noramlized_intestine_flourescence, ugt.14_ERT60$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)
wilcox.test(L4440_MR142$noramlized_intestine_flourescence, ugt.14_MR142$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)


# for endu-2
t.test(L4440_JM149$noramlized_intestine_flourescence, endu.2_JM149$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(L4440_ERT60$noramlized_intestine_flourescence, endu.2_ERT60$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(L4440_MR142$noramlized_intestine_flourescence, endu.2_MR142$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

wilcox.test(L4440_JM149$noramlized_intestine_flourescence, endu.2_JM149$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)
wilcox.test(L4440_ERT60$noramlized_intestine_flourescence, endu.2_ERT60$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)
wilcox.test(L4440_MR142$noramlized_intestine_flourescence, endu.2_MR142$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)


# for clec-56
t.test(L4440_JM149$noramlized_intestine_flourescence, clec.56_JM149$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(L4440_ERT60$noramlized_intestine_flourescence, clec.56_ERT60$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(L4440_MR142$noramlized_intestine_flourescence, clec.56_MR142$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

wilcox.test(L4440_JM149$noramlized_intestine_flourescence, clec.56_JM149$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)
wilcox.test(L4440_ERT60$noramlized_intestine_flourescence, clec.56_ERT60$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)
wilcox.test(L4440_MR142$noramlized_intestine_flourescence, clec.56_MR142$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)


# for c14c6.5
t.test(L4440_JM149$noramlized_intestine_flourescence, c14c6.5_JM149$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(L4440_ERT60$noramlized_intestine_flourescence, c14c6.5_ERT60$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(L4440_MR142$noramlized_intestine_flourescence, c14c6.5_MR142$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

wilcox.test(L4440_JM149$noramlized_intestine_flourescence, c14c6.5_JM149$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)
wilcox.test(L4440_ERT60$noramlized_intestine_flourescence, c14c6.5_ERT60$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)
wilcox.test(L4440_MR142$noramlized_intestine_flourescence, c14c6.5_MR142$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)


# for pbo-4
t.test(L4440_JM149$noramlized_intestine_flourescence, pbo.4_JM149$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(L4440_ERT60$noramlized_intestine_flourescence, pbo.4_ERT60$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(L4440_MR142$noramlized_intestine_flourescence, pbo.4_MR142$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

wilcox.test(L4440_JM149$noramlized_intestine_flourescence, pbo.4_JM149$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)
wilcox.test(L4440_ERT60$noramlized_intestine_flourescence, pbo.4_ERT60$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)
wilcox.test(L4440_MR142$noramlized_intestine_flourescence, pbo.4_MR142$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)
