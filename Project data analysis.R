# 0 Set up ---------------------------------------------------------------------
library(ggplot2) # for creating plots
library(car) # for Levene's test
library(AICcmodavg) # for Akaike information criterion (AIC) test for models
library(MHTdiscrete) # for Sidak post-hoc p-value adjustment
library(gridExtra) # allows multiple graphs to be arranged nicely in one image

#setwd("C:/Users/El Richardson/OneDrive - Lancaster University/Biomedicine/Year 
  #3/387 Project/Lab work")


# 1 Data processing ------------------------------------------------------------
data <- read.csv("microglia_data.csv")
View(data)

# Adding in extra columns of information to the data frame-----
data$split <- data$file_name |>
  strsplit(split = " ")

# Brain region and animal number
data$animal <- NA
data$region <- NA

for(i in 1:nrow(data)){
  data$animal[[i]] <- data$split[[i]][1]
  data$region[[i]] <- data$split[[i]][2]
}

# sex and diet (unbound)
data$sex <- NA
data$diet <- NA 

for(i in 1:nrow(data)){
  if(data$animal[[i]] == "PC44" | data$animal[[i]] == "PC45" | 
     data$animal[[i]] == "PC71"){
    data$sex[[i]] <- "male"
    data$diet[[i]] <- "C/C"
  }
  else if(data$animal[[i]] == "PC1" | data$animal[[i]] == "PC33" | 
          data$animal[[i]] == "PC40"){
    data$sex[[i]] <- "male"
    data$diet[[i]] <- "HF/C"
  }
  else if(data$animal[[i]] == "PC20" | data$animal[[i]] == "PC49" | 
          data$animal[[i]] == "PC70"){
    data$sex[[i]] <- "female"
    data$diet[[i]] <- "C/C"
  }
  else if(data$animal[[i]] == "PC6" | data$animal[[i]] == "PC32" | 
          data$animal[[i]] == "PC57"){
    data$sex[[i]] <- "female"
    data$diet[[i]] <- "HF/C"
  }
}

# Quantifying variability in light between the two microscopy sessions
# This is done through comparison of data obtained from capturing an image of 
# the same section of the same sample in both sessions
day_1_control <- data[data$file_name == 
                        "PC44 ctx x10 a - day 2.jpg",]["percent_area"] |> 
  as.numeric()
day_2_control <- data[data$file_name == "PC44 ctx x10 a.jpg",]["percent_area"]|> 
  as.numeric()

session_variability <- day_2_control/day_1_control

# Adjusting % coverage from day 1 to be comparable to day 2 
data$percent_area_adjusted <- NA

for(i in 1:nrow(data)){
  
  if(data$animal[[i]] == "PC40" | data$animal[[i]] == "PC44" | 
     data$animal[[i]] == "PC71"){
    
    unadjusted <- as.numeric(data$percent_area[[i]])
    data$percent_area_adjusted[[i]] <- unadjusted*session_variability
  }
  else{
    data$percent_area_adjusted[[i]] <- data$percent_area[[i]]
  }
}

data$percent_area_adjusted <- as.numeric(data$percent_area_adjusted)

# Creating a new data frame to export and include in paper
export_df <- data[c("animal", "diet", "sex", "region", "count", "total_area", 
                    "average_size", "percent_area", "percent_area_adjusted")]
#write.csv(export_df, "processed data for paper.csv", row.names = FALSE)


# 3 Plot -----------------------------------------------------------------------
# Creating subsets of the data that I want to look at
subset_list <- list(
  "control" = data[data$diet == "C/C", "percent_area_adjusted"],
  "test" = data[data$diet == "HF/C", "percent_area_adjusted"],
  
  "control_male" = 
    data[data$diet == "C/C" & data$sex == "male", "percent_area_adjusted"],
  "contol_female" = 
    data[data$diet == "C/C" & data$sex == "female", "percent_area_adjusted"],
  "test_male" = 
    data[data$diet == "HF/C" & data$sex == "male", "percent_area_adjusted"],
  "test_female" = 
    data[data$diet == "HF/C" & data$sex == "female", "percent_area_adjusted"],
  
  "control_ctx" = 
    data[data$diet == "C/C" & data$region == "ctx", "percent_area_adjusted"],
  "control_cpu" = 
    data[data$diet == "C/C" & data$region == "cpu", "percent_area_adjusted"],
  "test_ctx" = 
    data[data$diet == "HF/C" & data$region == "ctx", "percent_area_adjusted"],
  "test_cpu" = 
    data[data$diet == "HF/C" & data$region == "cpu", "percent_area_adjusted"],
  
  "control_male_ctx" = data[data$diet == "C/C" & data$sex == "male" & 
                              data$region == "ctx", "percent_area_adjusted"],
  "control_female_ctx" = data[data$diet == "C/C" & data$sex == "female" & 
                                data$region == "ctx", "percent_area_adjusted"],
  "test_male_ctx" = data[data$diet == "HF/C" & data$sex == "male" & 
                           data$region == "ctx", "percent_area_adjusted"],
  "test_female_ctx" = data[data$diet == "HF/C" & data$sex == "female" &
                             data$region == "ctx", "percent_area_adjusted"],
  "control_male_cpu" = data[data$diet == "C/C" & data$sex == "male" & 
                              data$region == "cpu", "percent_area_adjusted"],
  "contol_female_cpu" = data[data$diet == "C/C" & data$sex == "female" & 
                               data$region == "cpu", "percent_area_adjusted"],
  "test_male_cpu" = data[data$diet == "HF/C" & data$sex == "male" & 
                           data$region == "cpu", "percent_area_adjusted"],
  "test_female_cpu" = data[data$diet == "HF/C" & data$sex == "female" & 
                             data$region == "cpu", "percent_area_adjusted"]
  )

# Creating a new data frame for statistical results
length(subset_list)
stats_df <- data.frame(matrix(ncol = 5, nrow = 18))
colnames(stats_df) <- c("test", "mean", "se", "shapiro_p", "h0_accept")

# Finding means, standard deviations, and normality of data subsets
se <- function(x) sd(x)/sqrt(length(x)) # creating a standard error function

for(i in 1:length(subset_list)){
  name <- subset_list[i]|> names()
  subset <- subset_list[[i]]
  subset_df <- data.frame(subset)
  
  stats_df[[1]][i] <- name
  stats_df[[2]][i] <- subset |> mean()
  stats_df[[3]][i] <- subset |> se()
  shapiro <- subset |> shapiro.test()
  stats_df[[4]][i] <- shapiro$p.value
  
  if(shapiro$p.value>0.05){
    stats_df[[5]][i] <- TRUE
  } else{
    stats_df[[5]][i] <- FALSE
  }
  
  colnames(subset_df) <- "percent" # (allows ggplot to plot varying columns)
  
  plot <- # Creating plots to display the distribution of the subset data
    ggplot(data = subset_df, mapping = aes(x = percent)) + 
    geom_density(colour = "cornflowerblue", fill = alpha("cornflowerblue", 0.3)) + 
    theme_minimal() + 
    ylab("Density") + xlab("Area Covered by Microglia %")
  #ggsave(plot = plot, filename = paste0("Graphs/", name, ".png"), width = 6.25, 
   #     height = 5)
}

# Creating and saving mean/standard error plots
oneway_diet_plot <- 
  stats_df[1:2,] |>
  ggplot(aes(x = test, y = mean)) +
  geom_bar(stat = "identity", width = 0.5, fill = "cornflowerblue") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_x_discrete(labels=c("C/C", "HF/C")) + 
  xlab("Paternal Diet") + ylab("Mean Microglia Coverage (%)") + 
  theme_bw()
#ggsave(plot = oneway_diet_plot, filename = "Graphs/comparison oneway diet.png", 
 #      width = 6.25, height = 5)

twoway_dietsex_plot <- 
  stats_df[3:6,] |>
  ggplot(aes(x = test, y = mean,)) +
  geom_bar(stat = "identity", width = 0.5, fill = "cornflowerblue") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_x_discrete(labels=c("C/C Female", "C/C Male", "HF/C Female", 
                            "HF/C Male")) + 
  xlab("Paternal Diet and Offspring sex") + ylab("Mean Microglia Coverage (%)")+
  theme_bw()
#ggsave(plot = twoway_dietsex_plot, filename = 
 #        "Graphs/comparison twoway diet sex.png", width = 6.25, height = 5)

twoway_dietregion_plot <- 
  stats_df[7:10,] |>
  ggplot(aes(x = test, y = mean)) +
  geom_bar(stat = "identity", width = 0.5, fill = "cornflowerblue") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_x_discrete(labels=c("C/C Striatum", "C/C Cortex", "HF/C Striatum", "HF/C
                            Cortex")) + 
  xlab("Paternal Diet and Offspring Brain Region") + ylab("Mean Microglia 
                                                          Coverage (%)") + 
  theme_bw()
#ggsave(plot = twoway_dietregion_plot, filename = "Graphs/comparison twoway diet
  #region.png", width = 6.25, height = 5)
 
threeway_dietsexregion_plot <-
  stats_df[11:18,] |>
  ggplot(aes(x = test, y = mean)) +
  geom_bar(stat = "identity", width = 0.5, fill = "cornflowerblue") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_x_discrete(labels=c("C/C Female Striatum", "C/C Female Cortex", 
                            "C/C Male Striatum", "C/C Male Cortex", 
                            "HF/C Female Striatum", "HF/C Female Cortex", 
                            "HF/C Male Striatum",  "HF/C Male Cortex")) + 
  xlab("Paternal Diet and Offspring Brain Region") + ylab("Mean Microglia 
                                                          Coverage (%)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
#ggsave(plot = threeway_dietsexregion_plot, filename = "Graphs/comparison 
 #      threeway diet sex region.png", width = 6.25, height = 5)


# 4 Significance Testing -------------------------------------------------------
# One way testing: diet
leveneTest(percent_area_adjusted ~ diet, data = data) # 0.272 equal variance

wilcox.test(subset_list[["control"]], subset_list[["test"]], 
            alternative = "two.sided", exact = FALSE)

# Two way testing (both additive and interactive models): diet and sex:
#leveneTest(percent_area_adjusted ~ diet*sex, data = data)
# Archived since apprently one/two way ANOVAs are robust to variance.
twoway_test_sex <- aov(percent_area_adjusted ~ diet+sex, data = data)
summary(twoway_test_sex)
twoway_test_sex_interaction <- aov(percent_area_adjusted ~ diet*sex, 
                                   data = data)
summary(twoway_test_sex)

# Two way testing (both additive and interactive models): diet and region
twoway_test_region <- aov(percent_area_adjusted ~ diet+region, data = data)
summary(twoway_test_region)
twoway_test_region_interaction <- aov(percent_area_adjusted ~ diet*region, 
                                      data = data)
summary(twoway_test_region)

# Three way testing (both additive and interactive models): diet, sex and region
threeway_test <- aov(percent_area_adjusted ~ diet+region+sex, data = data)
summary(threeway_test)
threeway_test_rsinteraction <- aov(percent_area_adjusted ~ diet+region*sex, 
                                   data = data)
summary(threeway_test)
threeway_test_drinteraction <- aov(percent_area_adjusted ~ diet*region+sex, 
                                   data = data)
summary(threeway_test)
threeway_test_interaction <- aov(percent_area_adjusted ~ diet*region*sex, 
                                 data = data)
summary(threeway_test)

# Comparing different models created by ANOVA to see which is the best fit
model_set <- list(twoway_test_sex, twoway_test_sex_interaction, 
                  twoway_test_region, twoway_test_region_interaction,
                  threeway_test, threeway_test_rsinteraction, 
                  threeway_test_drinteraction,threeway_test_interaction)

model_names <- c("twoway_test_sex", "twoway_test_sex_interaction", 
                 "twoway_test_region", "twoway_test_region_interaction",
                 "threeway_test", "threeway_test_rsinteraction", 
                 "threeway_test_drinteraction","threeway_test_interaction")

aictab(model_set, modnames = model_names) # "twoway_test_sex" is the best fit
# This means that sex and diet affect microglia coverage in additive manner


# 5 Post-hoc correction of significance values --=------------------------------
Sidak.p.adjust(p = c(0.000488, 0.000283), alpha = 0.05) # twoway_sex p-values

# 6 Re-plot for paper ----------------------------------------------------------
significant_box_plot <- 
  data |>
  ggplot(aes(x = diet, y = percent_area_adjusted, fill = sex)) +
  geom_boxplot() + geom_point(position=position_jitterdodge(0.05)) + 
  guides(fill=guide_legend(title="Offspring Sex")) +
  scale_fill_manual(values = c("#61B499", "#8E61B4")) + 
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))) + 
  xlab("Paternal Diet") + ylab("Microglia Coverage (%)")+
  theme_bw()
#ggsave(plot = significant_box_plot, filename = 
 #       "Graphs/significant box plot.png", width = 6.25, height = 5)

MinMeanSEMMax <- function(x) {
  v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + 
           sd(x)/sqrt(length(x)), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}

significant_box_plot_se <- 
  data |>
  ggplot(aes(x = diet, y = percent_area_adjusted, fill = sex)) +
  stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", 
               position=position_dodge(0.8)) +
  geom_point(position=position_jitterdodge(0.05)) + 
  guides(fill=guide_legend(title="Offspring Sex")) +
  scale_fill_manual(values = c("#61B499", "#8E61B4")) + 
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24), 
                     limits = c(0, 23), expand = expansion(mult = c(0, 0.05))) + 
  xlab("Paternal Diet") + ylab("Microglia Coverage (%)")+
  theme_bw()
#ggsave(plot = significant_box_plot_se, filename = 
 #      "Graphs/significant box plot se.png", width = 6.25, height = 5)

significant__plot_df <- stats_df[3:6,]
significant__plot_df$sex <- c("male", "female", "male", "female")
significant__plot_df$diet <- c("C/C", "C/C", "HF/C", "HF/C")

significant_bar_plot <- 
  significant__plot_df |>
  ggplot(aes(x = diet, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position=position_dodge(), width = 0.5) +
  guides(fill=guide_legend(title="Offspring Sex")) +
  scale_fill_manual(values = c("#61B499", "#8E61B4")) + 
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1), 
                position=position_dodge(.5)) + 
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14), 
                     expand = expansion(mult = c(0, 0.05))) + 
  xlab("Paternal Diet") + ylab("Mean Microglia Coverage (%)")+
  theme_bw()
#ggsave(plot = significant_bar_plot, filename = 
 #       "Graphs/significant bar plot.png", width = 6.25, height = 5)

# Arranging two plots in the same image
significant_box_plot <- significant_box_plot + theme(legend.position = "none")

combined <- grid.arrange(significant_box_plot, significant_bar_plot, 
          ncol = 2, nrow = 1)
#ggsave(plot = combined, filename = 
 #      "Graphs/significant combined plot.png", width = 10, height = 5)
