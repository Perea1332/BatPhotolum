# Glowing Green: A Quantitative Analysis of Photoluminescence in Six North American Bat Species
# Roberson, B., Perea, S., Derose-Broeckert, D., & Castleberry, S. (2025). Glowing green: A quantitative analysis of photoluminescence in 
# six North American bats. https://doi.org/10.1002/ece3.71885

# Packages 
library(dplyr)
library(data.table)
library(ggplot2)
library(vegan)
library(glmmTMB)

# importing and merging individual csv scans
setwd('C:/Users/Bri/OneDrive - University of Georgia/Desktop/Senior Thesis/data analysis/CSVs')
files <- list.files(pattern = ".csv")
temp <- lapply(files, fread, sep=",")
joined_data <- rbindlist(temp)


# subtracting out dark and reference scans to calculate photoluminescence
joined_data$Photoluminescence <- joined_data$Value - (joined_data$Dark + joined_data$Reference)
str(joined_data)

# check that all specimens are present
table(joined_data$Species)


# calculate peak wavelength & irradiance (subsetting out noise & outside of range of emmission)
max_wavelength <- joined_data %>%
  filter(Wavelength >= 450 & Wavelength <= 600) %>%
  group_by(Species, ID) %>%
  filter(Photoluminescence == max(Photoluminescence)) %>%
  mutate(Photoluminescence = if_else(Photoluminescence < 0, 0, Photoluminescence)) %>%  
  select(Species, ID, Wavelength, Photoluminescence, Sex) %>%
  ungroup()

max_wavelength <- max_wavelength %>%
  rename(MaxWavelength = Wavelength, MaxPhotoluminescence = Photoluminescence)


range(max_wavelength$MaxWavelength)


write.csv(max_wavelength, "max_wavelength.csv", row.names = FALSE)


mean_by_species <- max_wavelength %>%
  group_by(Species) %>%
  summarise(mean_MaxWavelength = mean(MaxWavelength, na.rm = TRUE))

# Plotting spectral characteristics --------------------------------------------
joined_data$Species <- factor(joined_data$Species, 
                              levels = c("EPFU", "TABR", "MYAU", "MYGR", "LABO", "LASE"))

ggplot(joined_data %>% filter(Wavelength != 620), aes(Wavelength, Photoluminescence, color = Species, fill = Species, group = Species)) +
  geom_smooth(method = "gam",
              alpha = 0.2,  # Transparency for CI
              se = TRUE,
              span = 1.5,
              level = 0.95,
              na.rm = TRUE) +
  scale_x_continuous(limits = c(480, 700)) +
  scale_y_continuous(limits = c(0, 0.0075)) +
  scale_color_manual(values = c("LASE" = "#481567FF", 
                                "LABO" = "#33638DFF", 
                                "MYGR" = "#1F968BFF", 
                                "MYAU" = "#3CBB75FF", 
                                "TABR" = "#95D840FF", 
                                "EPFU" = "#FDE725FF"), 
                     labels = c("LASE" = expression(italic("Lasiurus seminolus")),
                                "LABO" = expression(italic("Lasiurus borealis")),
                                "MYGR" = expression(italic("Myotis grisescens")),
                                "MYAU" = expression(italic("Myotis austroriparius")),
                                "TABR" = expression(italic("Tadarida brasiliensis")),
                                "EPFU" = expression(italic("Eptesicus fuscus")))) +
  scale_fill_manual(values = c("LASE" = "#481567FF", 
                               "LABO" = "#33638DFF", 
                               "MYGR" = "#1F968BFF", 
                               "MYAU" = "#3CBB75FF", 
                               "TABR" = "#95D840FF", 
                               "EPFU" = "#FDE725FF"), 
                    labels = c("LASE" = expression(italic("Lasiurus seminolus")),
                               "LABO" = expression(italic("Lasiurus borealis")),
                               "MYGR" = expression(italic("Myotis grisescens")),
                               "MYAU" = expression(italic("Myotis austroriparius")),
                               "TABR" = expression(italic("Tadarida brasiliensis")),
                               "EPFU" = expression(italic("Eptesicus fuscus")))) +
  xlab("Wavelength (nm)") +
  ylab(expression("Irradiance (μW/cm"^2*")")) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size=10, color = "#101010", hjust = 0.2),
        axis.text.x = element_text(size=10,margin = margin(t = 0, r = 0, b = 0, l = 0), vjust = -0.5, color = "#101010"),
        axis.title.y = element_text(size=12, color = "#101010", hjust = 0.5, l = 0, vjust = 5 , margin = margin(l = 20)),
        axis.title.x = element_text(size=12, vjust = -5, color = "#101010" , margin = margin(b = 20)),
        legend.text = element_text(size=10, color = "#101010"),
        legend.title = element_text(size=12, color = "#101010"), 
        plot.margin = margin(t = 20))



# Statistics -------------------------------------------------------------------
# Differences in peak wavelength 

# Variable = MaxWavelength, which is the wavelength at the max irradiance value 
# Test for normality: Shapiro-Wilkes test
class(max_wavelength)
max_wavelength_numeric <- as.numeric(max_wavelength$MaxWavelength)
peak_normality_wavelength <- shapiro.test(max_wavelength_numeric)
print(peak_normality_wavelength)
#W = 0.74795, p-value = 8.412e-09 -> Data is not normally distributed 


# Test for beta-diversity: PERMDISP , is there a significant difference in dispersion of wavelength between species?
distance_matrix_wavelength <- dist(max_wavelength$MaxWavelength)
print(distance_matrix_wavelength)
permdisp_result_wavelength <- betadisper(distance_matrix_wavelength, group = max_wavelength$Species)
summary(permdisp_result_wavelength)
str(permdisp_result_wavelength$vectors)

permdisp_result_wavelength
permutest_wavelength <- permutest(permdisp_result_wavelength, permutations = 999)
print(permutest_wavelength)
# F=1.8633, p=0.105 -> no significant differences in variation of wavelength
# no permdisp in final analysis 


# Test for significant interspecific differences in peak wavelength
# Kruskal-Wallis test
kruskal.test(MaxWavelength ~ Species, data = max_wavelength)
#chi-squared = 10.595, df = 5, p-value = 0.06002 -> there is not a significant difference in peak wavelength between species 

# Differences in sex
# Because there is not a significant difference of peak wavelength between species, we can group sexes of different species together 

# Kruskal-Wallis test
kruskal.test(MaxPhotoluminescence ~ Sex, data = max_wavelength)
# chi-squared = 0.80365, df = 1, p-value = 0.37 -> there is no significant difference in peak wavelength between sexes 

# Influence of specimen age 
specimen_data <- read.csv("C:/Users/Bri/OneDrive - University of Georgia/Desktop/Senior Thesis/data analysis/specimen_information.csv")
max_wavelength <- merge(max_wavelength, specimen_data[, c("ID", "Age")], by = "ID", all.x = TRUE)

# Test for normality: Shapiro-Wilkes test
Age_numeric <- as.numeric(max_wavelength$Age)
normality_age <- shapiro.test(Age_numeric)
print(normality_age)
#W = 0.92564, p-value = 0.001607, Age is not normally distributed 

# Remove specimens with "ND" (no date) listed 
max_wavelength_clean <- max_wavelength %>%
  filter(!is.na(Age), !is.na(MaxPhotoluminescence)) %>%
  mutate(MaxPhotoluminescence = if_else(MaxPhotoluminescence <= 0, 0.000001, MaxPhotoluminescence))
# check that all specimens are present
table(max_wavelength_clean$Species)

#Scatterplot of age vs irradiance - for all species 
max_wavelength_clean$Species <- factor(max_wavelength_clean$Species, 
                                       levels = c("EPFU", "TABR", "MYAU", "MYGR", "LABO", "LASE"))


age_scatterplot <- ggplot(max_wavelength_clean, aes(x = Age, y = MaxWavelength, color = Species)) +
  geom_point(size = 3) +
  geom_smooth(method = "glm", se = FALSE, color = "black") +
  facet_wrap(~ Species, labeller = label_parsed) +
  scale_color_manual(values = c("LASE" = "#481567FF", 
                                "LABO" = "#33638DFF", 
                                "MYGR" = "#1F968BFF", 
                                "MYAU" = "#3CBB75FF", 
                                "TABR" = "#95D840FF", 
                                "EPFU" = "#FDE725FF"), 
                     labels = c("LASE" = expression(italic("Lasiurus seminolus")),
                                "LABO" = expression(italic("Lasiurus borealis")),
                                "MYGR" = expression(italic("Myotis grisescens")),
                                "MYAU" = expression(italic("Myotis austroriparius")),
                                "TABR" = expression(italic("Tadarida brasiliensis")),
                                "EPFU" = expression(italic("Eptesicus fuscus")))) +
  scale_fill_manual(values = c("LASE" = "#481567FF", 
                               "LABO" = "#33638DFF", 
                               "MYGR" = "#1F968BFF", 
                               "MYAU" = "#3CBB75FF", 
                               "TABR" = "#95D840FF", 
                               "EPFU" = "#FDE725FF"), 
                    labels = c("LASE" = expression(italic("Lasiurus seminolus")),
                               "LABO" = expression(italic("Lasiurus borealis")),
                               "MYGR" = expression(italic("Myotis grisescens")),
                               "MYAU" = expression(italic("Myotis austroriparius")),
                               "TABR" = expression(italic("Tadarida brasiliensis")),
                               "EPFU" = expression(italic("Eptesicus fuscus")))) +
  xlab("Specimen Age") +
  ylab(expression("Wavelength")) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size=10, color = "#101010", hjust = 0.2),
        axis.text.x = element_text(size=10,margin = margin(t = 0, r = 0, b = 0, l = 0), vjust = -0.5, color = "#101010"),
        axis.title.y = element_text(size=12, face = "bold" , color = "#101010", hjust = 0.5, l = 0, vjust = 5 , margin = margin(l = 20)),
        axis.title.x = element_text(size=12, face = "bold" , vjust = -5, color = "#101010" , margin = margin(b = 20)),
        legend.text = element_text(size=10, color = "#101010"),
        legend.title = element_text(size=10, color = "#101010"), 
        plot.margin = margin(t = 20)) 

print(age_scatterplot)


# GLM null model analysis - using Gamma family due to continuous, positive, non normal data distribution

# null model
null_model_age_glm <- glm(MaxWavelength ~ 1, 
                                  data = max_wavelength_clean, 
                                  family = Gamma(link = "log"))
# age model 
predict_model_age_glm <- glm(MaxWavelength ~ Age, 
                                     data = max_wavelength_clean, 
                                     family = Gamma(link = "log"))
# comparison
AIC(null_model_age_glm, predict_model_age_glm)
#difference is 0.596, so there is not a significant difference in model fit (<2)

# Mean wavelength values 

# plotting mean peak wavelength values 
# summarize 
summary_max_wavelength <- max_wavelength %>%
  group_by(Species) %>%
  summarise(
    MeanWavelength = mean(MaxWavelength, na.rm = TRUE),
    SD = sd(MaxWavelength, na.rm = TRUE),
    SE = SD / sqrt(n()) 
  )

# plot


ggplot(summary_max_wavelength, aes(x = Species, y = MeanWavelength, color = Species)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = MeanWavelength - SE, ymax = MeanWavelength + SE), 
                width = 0.2) +
  coord_cartesian(ylim = c(500, 560)) +
  labs(x = "Species", y = "Peak Wavelength (nm)") +
  scale_color_manual(values = c("LASE" = "#481567FF", 
                                "LABO" = "#33638DFF", 
                                "MYGR" = "#1F968BFF", 
                                "MYAU" = "#3CBB75FF", 
                                "TABR" = "#95D840FF", 
                                "EPFU" = "#FDE725FF"),
                     labels = c("LASE" = expression(italic("Lasiurus seminolus")),
                                "LABO" = expression(italic("Lasiurus borealis")),
                                "MYGR" = expression(italic("Myotis grisescens")),
                                "MYAU" = expression(italic("Myotis austroriparius")),
                                "TABR" = expression(italic("Tadarida brasiliensis")),
                                "EPFU" = expression(italic("Eptesicus fuscus")))) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), 
    axis.title.y = element_text(size=12, color = "#101010", hjust = 0.5, l = 0, vjust = 5 , margin = margin(l = 20)),
    axis.title.x = element_text(size=12, vjust = -5, color = "#101010" , margin = margin(b = 20)), 
    axis.text.y = element_text(size=10, color = "#101010", hjust = 0.2),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
    panel.grid = element_blank(),
    legend.text = element_text(size=10, color = "#101010"),
    legend.title = element_text(size=12, color = "#101010"), 
    plot.margin = margin(t = 20)) 
    
