# Script for graphing RNAi quantification data


## read in and clean the data
getwd()

library(readr)

data <- read_csv("./RNAi_data/RNAi_Data_for_analysis.csv")
View(data)

library(dplyr)
clean_data <- data %>% filter(!grepl("background", image_number))
view(clean_data)
  

group_means <- clean_data %>% group_by(worm_strain, RNAi_strain) %>% summarise(mean_fluor = mean(noramlized_intestine_flourescence, na.rm = TRUE), .groups = "drop")
view(group_means)

control_means <- group_means %>% filter(RNAi_strain == "L4440") %>% dplyr::select(worm_strain, control_mean = mean_fluor)
view(control_means)

normalized_group_means <- group_means %>% left_join(control_means, by = "worm_strain") %>% mutate(normalized_to_L4440 = mean_fluor / control_mean)
view(normalized_group_means)


## graph the data 

library(ggplot2)

### re-factor the RNAistrains for the graph

normalized_group_means_ordered <- normalized_group_means %>% mutate(RNAi_strain = factor(RNAi_strain, levels = c("L4440", "pop-1", "y32f6a.5", "cpr-1", "ugt-14", "endu-2", "clec-56", "c14c6.5", "pbo-4")))
view(normalized_group_means_ordered)


pdf("./RNAi_data/JM149_RNAi_data.pdf", height = 6, width = 9)
ggplot(normalized_group_means_ordered %>% filter(worm_strain == "JM149"), aes(x = RNAi_strain, y = normalized_to_L4440)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  theme_minimal() +
  labs(
    x = "RNAi strain",
    y = "Normalized fluorescence (to L4440)",
    title = paste("Normalized Intestinal Fluorescence (JM149)")
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


pdf("./RNAi_data/ERT60_RNAi_data.pdf", height = 6, width = 9)
ggplot(normalized_group_means_ordered %>% filter(worm_strain == "ERT60"), aes(x = RNAi_strain, y = normalized_to_L4440)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  theme_minimal() +
  labs(
    x = "RNAi strain",
    y = "Normalized fluorescence (to L4440)",
    title = paste("Normalized Intestinal Fluorescence (ERT60)")
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



pdf("./RNAi_data/MR142_RNAi_data.pdf", height = 6, width = 9)
ggplot(normalized_group_means_ordered %>% filter(worm_strain == "MR142"), aes(x = RNAi_strain, y = normalized_to_L4440)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  theme_minimal() +
  labs(
    x = "RNAi strain",
    y = "Normalized fluorescence (to L4440)",
    title = paste("Normalized Intestinal Fluorescence (MR142)")
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



pdf("./RNAi_data/all_worm_strains_RNAi_data.pdf", height = 8, width = 12)
ggplot(
  normalized_group_means_ordered,
  aes(x = RNAi_strain, y = normalized_to_L4440)
) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  facet_wrap(~ worm_strain, ncol = 3) +
  theme_minimal() +
  labs(
    x = "RNAi strain",
    y = "Normalized fluorescence (to L4440)",
    title = "Normalized Intestinal Fluorescence"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off()


### students t test

data_for_Ttest <- clean_data[, c("RNAi_strain", "noramlized_intestine_flourescence")]
view(data_for_Ttest)

pop.1 <- data_for_Ttest[data_for_Ttest$RNAi_strain == "pop-1", ]
L4440 <- data_for_Ttest[data_for_Ttest$RNAi_strain == "L4440", ]
y32f6a.5 <- data_for_Ttest[data_for_Ttest$RNAi_strain == "y32f6a.5", ]
cpr.1 <- data_for_Ttest[data_for_Ttest$RNAi_strain == "cpr-1", ]
ugt.14 <- data_for_Ttest[data_for_Ttest$RNAi_strain == "ugt-14", ]
endu.2 <- data_for_Ttest[data_for_Ttest$RNAi_strain == "endu-2", ]
clec.56 <- data_for_Ttest[data_for_Ttest$RNAi_strain == "clec-56", ]
c14c6.5 <- data_for_Ttest[data_for_Ttest$RNAi_strain == "c14c6.5", ]
pbo.4 <- data_for_Ttest[data_for_Ttest$RNAi_strain == "pbo-4", ]


# for pop-1
t.test(L4440$noramlized_intestine_flourescence, pop.1$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

wilcox.test(L4440$noramlized_intestine_flourescence, pop.1$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)


# for y32f6a-5
t.test(L4440$noramlized_intestine_flourescence, y32f6a.5$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

wilcox.test(L4440$noramlized_intestine_flourescence, y32f6a.5$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)


# for cpr-1
t.test(L4440$noramlized_intestine_flourescence, cpr.1$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

wilcox.test(L4440$noramlized_intestine_flourescence, cpr.1$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)


# for ugt-14
t.test(L4440$noramlized_intestine_flourescence, ugt.14$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

wilcox.test(L4440$noramlized_intestine_flourescence, ugt.14$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)


# for endu-2
t.test(L4440$noramlized_intestine_flourescence, endu.2$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

wilcox.test(L4440$noramlized_intestine_flourescence, endu.2$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)


# for clec-56
t.test(L4440$noramlized_intestine_flourescence, clec.56$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

wilcox.test(L4440$noramlized_intestine_flourescence, clec.56$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)


# for c14c6.5
t.test(L4440$noramlized_intestine_flourescence, c14c6.5$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

wilcox.test(L4440$noramlized_intestine_flourescence, c14c6.5$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)


# for pbo-4
t.test(L4440$noramlized_intestine_flourescence, pbo.4$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

wilcox.test(L4440$noramlized_intestine_flourescence, pbo.4$noramlized_intestine_flourescence, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, conf.level = 0.95)

