# Script for graphing RNAi quantification data

getwd()

library(readr)
library(dplyr)
library(ggplot2)

## read in data
data <- read_csv("./RNAi_data/RNAi_Data_for_analysis.csv")
View(data)

## clean data by removing background image readings 
clean_data <- data %>% filter(!grepl("background", image_number))
view(clean_data)
dim(clean_data)
  
## remove an outlier in MR142 group
clean_data_no_outlier <- clean_data %>% filter(!(worm_strain == "MR142" & image_number == "10"))
view(clean_data_no_outlier)
dim(clean_data_no_outlier)


## calc group means
group_means <- clean_data_no_outlier %>% group_by(worm_strain, RNAi_strain) %>% summarise(
    mean_fluor = mean(noramlized_intestine_flourescence, na.rm = TRUE),
    sd_fluor = sd(noramlized_intestine_flourescence, na.rm = TRUE),
    n = n(),
    sem_fluor = sd_fluor / sqrt(n),
    .groups = "drop")
view(group_means)

## calc negative control (L4440) means
control_means <- group_means %>% filter(RNAi_strain == "L4440") %>% dplyr::select(worm_strain, control_mean = mean_fluor)
view(control_means)

## normalize indivudal observations
clean_data_normalized <- clean_data_no_outlier %>% left_join(control_means, by = "worm_strain") %>% mutate(normalized_to_L4440 = noramlized_intestine_flourescence / control_mean)
view(clean_data_normalized)

## summarize normalized replicates
normalized_group_means <- clean_data_normalized %>% group_by(worm_strain, RNAi_strain) %>%
  summarise(mean_normalized = mean(normalized_to_L4440, na.rm = TRUE), sd_normalized = sd(normalized_to_L4440, na.rm = TRUE), n = n(), sem_normalized = sd_normalized / sqrt(n), .groups = "drop")
view(normalized_group_means)


## set factor order
strain_order <- c("L4440", "pop-1", "ugt-14", "cpr-1", "y32f6a.5", "c14c6.5", "endu-2", "clec-56", "pbo-4")

normalized_group_means <- normalized_group_means %>% mutate(RNAi_strain = factor(RNAi_strain, levels = strain_order))
clean_data_normalized <- clean_data_normalized %>% mutate(RNAi_strain = factor(RNAi_strain, levels = strain_order))


## graph the data 

pdf("./RNAi_data/JM149_RNAi_data.pdf", height = 6, width = 9)
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


pdf("./RNAi_data/ERT60_RNAi_data.pdf", height = 6, width = 9)
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


pdf("./RNAi_data/MR142_RNAi_data.pdf", height = 6, width = 9)
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


### students t test
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
