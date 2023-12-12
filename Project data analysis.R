# 0 Set up ---------------------------------------------------------------------
# Libraries for statistics
library(car) # for Levene's test
library(AICcmodavg) # for Akaike information criterion (AIC) test for models
library(MHTdiscrete) # for Sidak post-hoc p-value adjustment

# Libraries for plotting
library(ggplot2) # for creating plots
library(gridExtra) # allows multiple graphs to be arranged nicely in one image
library(ggsignif) # draws asterisks on ggplot plots to indicate significance

#setwd("C:/Users/El Richardson/OneDrive - Lancaster University/Biomedicine/Year 
  #3/387 Project/Lab work")


# 1 Data processing ------------------------------------------------------------
data <- read.csv("microglia_data.csv")
View(data)

# Adding in extra columns for brain region and animal number
data$region <- NA
data$animal <- NA

data$split <- data$file_name |>
  strsplit(split = " ")

for(i in 1:nrow(data)){
  data$animal[[i]] <- data$split[[i]][1]
  data$region[[i]] <- data$split[[i]][2]
}

# Adding in extra columns for sex and diet
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
day_1_control <- data[data$file_name == 
                        "PC44 ctx x10 a - day 2.jpg",]["percent_area"] |> 
  as.numeric()

day_2_control <- data[data$file_name == "PC44 ctx x10 a.jpg",]["percent_area"]|> 
  as.numeric()

session_variability <- day_2_control/day_1_control

# Removing the calibration duplicate data point
data <- data[-27,]
rownames(data) <- c(1:48)

# Adjusting % coverage from session 1 to be comparable to session 2
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
export_df <- data[c("animal", "diet", "sex", "region", 
                    "percent_area_adjusted")] #|> 
#  write.csv("processed data for paper.csv", row.names = FALSE)                  #


# 2 Creating subsets of the data that I want to look at-------------------------
subset_list <- list(
  "control" = data[data$diet == "C/C", "percent_area_adjusted"],
  "test" = data[data$diet == "HF/C", "percent_area_adjusted"],
  "male" = data[data$sex == "male", "percent_area_adjusted"],
  "female" = data[data$sex == "female", "percent_area_adjusted"],
  "ctx" = data[data$region == "ctx", "percent_area_adjusted"],
  "cpu" = data[data$region == "cpu", "percent_area_adjusted"],
  
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


# 3 Statistics -----------------------------------------------------------------
# Creating a new data frame for statistical results
length(subset_list) # seeing how long stats data frame needs to be 
stats_df <- data.frame(matrix(ncol = 5, nrow = 22))
colnames(stats_df) <- c("test", "mean", "se", "shapiro_p", "h0_accept")
View(stats_df)

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
    geom_density(colour = "cornflowerblue", 
                 fill = alpha("cornflowerblue", 0.3)) +
    theme_minimal() + 
    ylab("Density") + xlab("Area Covered by Microglia %")
  ggsave(plot = plot, filename = paste0("Graphs/", name, ".png"), width = 6.25, #
        height = 5)
}


# 4 Significance Testing -------------------------------------------------------
# One way testing: diet (non-parametric)
leveneTest(percent_area_adjusted ~ diet, data = data) # 0.348 equal variance

wilcox.test(subset_list[["control"]], subset_list[["test"]], 
            alternative = "two.sided", exact = FALSE) # p-value = 0.002698 **

# One way testing: sex (non-parametric)
leveneTest(percent_area_adjusted ~ sex, data = data) # 0.5587 equal variance

wilcox.test(subset_list[["female"]], subset_list[["male"]], 
            alternative = "two.sided", exact = FALSE) # p-value = 0.0008364 ***

# One way testing: region (non-parametric)
leveneTest(percent_area_adjusted ~ region, data = data) # 0.6869 equal variance

wilcox.test(subset_list[["ctx"]], subset_list[["cpu"]], 
            alternative = "two.sided", exact = FALSE) # p-value = 0.7649 NS

# Two way testing (both additive and interactive models): diet and sex:
twoway_test_sex <- aov(percent_area_adjusted ~ diet+sex, data = data)
summary(twoway_test_sex)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#diet         1  341.8   341.8   13.14 0.000732 ***
#  sex          1  390.8   390.8   15.03 0.000342 ***
#  Residuals   45 1170.4    26.0  

twoway_test_sex_interaction <- aov(percent_area_adjusted ~ diet*sex, 
                                   data = data)
summary(twoway_test_sex_interaction)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
#diet         1  341.8   341.8  13.386 0.000675 ***
#sex          1  390.8   390.8  15.303 0.000313 ***
#diet:sex     1   46.8    46.8   1.834 0.182537    
#Residuals   44 1123.6    25.5 

# Two way testing (both additive and interactive models): diet and region
twoway_test_region <- aov(percent_area_adjusted ~ diet+region, data = data)
summary(twoway_test_region)
#             Df Sum Sq Mean Sq F value  Pr(>F)   
#diet         1  341.8   341.8   9.904 0.00292 **
#region       1    8.1     8.1   0.235 0.63017   
#Residuals   45 1553.1    34.5  

twoway_test_region_interaction <- aov(percent_area_adjusted ~ diet*region, 
                                      data = data)
summary(twoway_test_region_interaction)
#             Df Sum Sq Mean Sq F value  Pr(>F)   
#diet         1  341.8   341.8   9.704 0.00323 **
#region       1    8.1     8.1   0.230 0.63370   
#diet:region  1    3.1     3.1   0.089 0.76710   
#Residuals   44 1550.0    35.2 

# Three way testing (both additive and interactive models): diet, sex and region
threeway_test <- aov(percent_area_adjusted ~ diet+region+sex, data = data)
summary(threeway_test)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
#diet         1  341.8   341.8  12.940 0.000810 ***
#region       1    8.1     8.1   0.307 0.582291    
#sex          1  390.8   390.8  14.793 0.000383 ***
#Residuals   44 1162.3    26.4  

threeway_test_rsinteraction <- aov(percent_area_adjusted ~ diet+region*sex, 
                                   data = data)
summary(threeway_test_rsinteraction)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
#diet         1  341.8   341.8  12.662 0.000925 ***
#region       1    8.1     8.1   0.300 0.586426    
#sex          1  390.8   390.8  14.475 0.000444 ***
#region:sex   1    1.5     1.5   0.054 0.817715    
#Residuals   43 1160.9    27.0 

threeway_test_drinteraction <- aov(percent_area_adjusted ~ diet*region+sex, 
                                   data = data)
summary(threeway_test_drinteraction)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
#diet         1  341.8   341.8  12.680 0.000918 ***
#region       1    8.1     8.1   0.301 0.586157    
#sex          1  390.8   390.8  14.496 0.000441 ***
#diet:region  1    3.1     3.1   0.116 0.735024    
#Residuals   43 1159.2    27.0  

threeway_test_interaction <- aov(percent_area_adjusted ~ diet*region*sex, 
                                 data = data)
summary(threeway_test_interaction)
#                 Df Sum Sq Mean Sq F value   Pr(>F)    
#diet             1  341.8   341.8  12.426 0.001077 ** 
#region           1    8.1     8.1   0.295 0.590133    
#sex              1  390.8   390.8  14.205 0.000529 ***
#diet:region      1    3.1     3.1   0.114 0.737715    
#diet:sex         1   46.8    46.8   1.703 0.199393    
#region:sex       1    1.5     1.5   0.053 0.819473    
#diet:region:sex  1   10.5    10.5   0.382 0.539913    
#Residuals       40 1100.4    27.5  

# Comparing different models created by ANOVA to see which is the best fit
model_set <- list(twoway_test_sex, twoway_test_sex_interaction, 
                  twoway_test_region, twoway_test_region_interaction,
                  threeway_test, threeway_test_rsinteraction, 
                  threeway_test_drinteraction,threeway_test_interaction)

model_names <- c("twoway_test_sex", "twoway_test_sex_interaction", 
                 "twoway_test_region", "twoway_test_region_interaction",
                 "threeway_test", "threeway_test_rsinteraction", 
                 "threeway_test_drinteraction","threeway_test_interaction")

aictab(model_set, modnames = model_names) 
#                               K   AICc Delta_AICc AICcWt Cum.Wt      LL
#twoway_test_sex                4 298.46       0.00   0.43   0.43 -144.76
#twoway_test_sex_interaction    5 298.99       0.54   0.33   0.77 -143.78
#threeway_test                  5 300.62       2.16   0.15   0.91 -144.60
#threeway_test_drinteraction    6 303.11       4.66   0.04   0.96 -144.53
#threeway_test_rsinteraction    6 303.18       4.72   0.04   1.00 -144.57
#threeway_test_interaction      9 309.30      10.84   0.00   1.00 -143.28
#twoway_test_region             4 312.03      13.58   0.00   1.00 -151.55
#twoway_test_region_interaction 5 314.44      15.98   0.00   1.00 -151.50
# twoway_test_sex has the highest model proability


# 5 Post-hoc test and correction of significance values-------------------------
data$sex_diet <- NA

for(i in 1:nrow(data)){ # Adding extra column to allow multiple comparisons 
  data$sex_diet[[i]] <- paste0(data$sex[[i]], data$diet[[i]]) |>
    as.character()
}

pairwise.t.test(data$percent_area_adjusted,  data$sex_diet, p.adj='bonferroni')
#             femaleC/C femaleHF/C maleC/C
#femaleHF/C   0.6622    -          -      
#  maleC/C    0.0033    1.8e-05    -      
#  maleHF/C   1.0000    0.4642     0.0057 


# 6 Plots ----------------------------------------------------------------------
# Creating and saving mean/standard error plots
oneway_diet_plot <- 
  stats_df[1:2,] |>
  ggplot(aes(x = test, y = mean)) +
  geom_bar(stat = "identity", width = 0.5, fill = "cornflowerblue") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1)) + 
  geom_signif(comparisons = list(c("control", "test")), annotation = "**") + 
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12),
                     expand = expansion(mult = c(0, 0.05))) + 
  scale_x_discrete(labels=c("C/C", "HF/C")) + 
  xlab("Paternal Diet") + ylab("Mean Microglia Coverage (%)") + 
  theme_bw()
#ggsave(plot = oneway_diet_plot, filename = "Graphs/comparison oneway diet.png", #
#      width = 6.25, height = 5)

oneway_sex_plot <- 
  stats_df[3:4,] |>
  ggplot(aes(x = test, y = mean)) +
  geom_bar(stat = "identity", width = 0.5, fill = "cornflowerblue") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1)) +
  geom_signif(comparisons = list(c("male", "female")), annotation = "***") + 
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12),
                     expand = expansion(mult = c(0, 0.05))) + 
  scale_x_discrete(labels=c("Female", "Male")) + 
  xlab("Offspring Sex") + ylab("Mean Microglia Coverage (%)") + 
  theme_bw()
#ggsave(plot = oneway_sex_plot, filename = "Graphs/comparison oneway sex.png",   #
#      width = 6.25, height = 5)

oneway_region_plot <- 
  stats_df[5:6,] |>
  ggplot(aes(x = test, y = mean)) +
  geom_bar(stat = "identity", width = 0.5, fill = "cornflowerblue") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1)) + 
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10),
                     expand = expansion(mult = c(0, 0.05))) + 
  scale_x_discrete(labels=c("Caudate Putamen", "Cortex")) + 
  xlab("Brain Region Imaged") + ylab("Mean Microglia Coverage (%)") + 
  theme_bw()
#ggsave(plot = oneway_region_plot, filename =                                    #
#         "Graphs/comparison oneway region.png", width = 6.25, height = 5)

twoway_dietsex_plot <- 
  stats_df[7:10,] |>
  ggplot(aes(x = test, y = mean,)) +
  geom_bar(stat = "identity", width = 0.5, fill = "cornflowerblue") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_x_discrete(labels=c("C/C Female", "C/C Male", "HF/C Female", 
                            "HF/C Male")) + 
  xlab("Paternal Diet and Offspring sex") + ylab("Mean Microglia Coverage (%)")+
  theme_bw()
#ggsave(plot = twoway_dietsex_plot, filename =                                   #
#        "Graphs/comparison twoway diet sex.png", width = 6.25, height = 5)

twoway_dietregion_plot <- 
  stats_df[11:14,] |>
  ggplot(aes(x = test, y = mean)) +
  geom_bar(stat = "identity", width = 0.5, fill = "cornflowerblue") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_x_discrete(labels=c("C/C Striatum", "C/C Cortex", "HF/C Striatum", "HF/C
                            Cortex")) + 
  xlab("Paternal Diet and Offspring Brain Region") + ylab("Mean Microglia 
                                                          Coverage (%)") + 
  theme_bw()
#ggsave(plot = twoway_dietregion_plot, filename =                                #
#"Graphs/comparison twoway diet region.png", width = 6.25, height = 5)

threeway_dietsexregion_plot <-
  stats_df[15:22,] |>
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
#ggsave(plot = threeway_dietsexregion_plot, filename =                           #
#"Graphs/comparison threeway diet sex region.png", width = 6.25, height = 5)


# 7 Plots for paractical report  -----------------------------------------------
# Box plot to show distribution
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
#ggsave(plot = significant_box_plot, filename =                                  #
#        "Graphs/significant box plot.png", width = 6.25, height = 5)

# Box plot with mean and standard error
MinMeanSEMMax <- function(x) { # creating a custom function to plot this
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
#ggsave(plot = significant_box_plot_se, filename =                               #
#       "Graphs/significant box plot se.png", width = 6.25, height = 5)

# Bar graph showing mean and standard error
significant__plot_df <- stats_df[3:6,]
significant__plot_df$sex <- c("male", "female", "male", "female")
significant__plot_df$diet <- c("C/C", "C/C", "HF/C", "HF/C")

significant_bar_plot_sex <- 
  significant__plot_df |>
  ggplot(aes(x = diet, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position=position_dodge(), width = 0.5) +
  geom_signif(stat = "identity",
              aes(x = 0.8, xend = 1.2, y = 12, yend = 12, annotation = "**")) + 
    guides(fill=guide_legend(title="Offspring Sex")) +
  scale_fill_manual(values = c("#61B499", "#8E61B4")) + 
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1), 
                position=position_dodge(.5)) + 
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14), 
                     expand = expansion(mult = c(0, 0.05))) + 
  xlab("Paternal Diet") + ylab("Mean Microglia Coverage (%)")+
  theme_bw()
#ggsave(plot = significant_bar_plot_sex, filename =                              #
#        "Graphs/significant bar plot sex.png", width = 6.25, height = 5)

significant_bar_plot_diet <- 
  significant__plot_df |>
  ggplot(aes(x = sex, y = mean, fill = diet)) +
  geom_bar(stat = "identity", position=position_dodge(), width = 0.5) +
  geom_signif(stat = "identity",
              aes(x = 1.8, xend = 2.2, y = 12, yend = 12, annotation = "**")) + 
  guides(fill=guide_legend(title="Paternal Diet")) +
  scale_fill_manual(values = c("#6185B4", "#B461AE")) + 
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1), 
                position=position_dodge(.5)) + 
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14), 
                     expand = expansion(mult = c(0, 0.05))) + 
  xlab("Offspring Sex") + ylab("Mean Microglia Coverage (%)")+
  theme_bw()
#ggsave(plot = significant_bar_plot_sex, filename =                              #
#         "Graphs/significant bar plot sex.png", width = 6.25, height = 5)

# Arranging two plots in the same image-----
# Sex and Region
oneway_sex_plot <- oneway_sex_plot + theme(axis.title.y = element_blank())
oneway_region_plot <- oneway_region_plot + theme(axis.title.y = element_blank())

combined_oneway <- grid.arrange(oneway_diet_plot, oneway_sex_plot, 
                                oneway_region_plot, ncol = 3, nrow = 1)
#ggsave(plot = combined_oneway, filename =                                       #
#      "Graphs/significant oneway plot.png", width = 10, height = 5)

# Sex and diet together
combined_twoway <- grid.arrange(
  significant_bar_plot_diet, significant_bar_plot_sex, ncol = 2, nrow = 1)
#ggsave(plot = combined_twoway, filename =                                       #
#       "Graphs/significant combined plot.png", width = 10, height = 5)
