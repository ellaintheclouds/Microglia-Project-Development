# 0 Set up ---------------------------------------------------------------------
# Libraries for statistics
library(car) # for Levene's test
library(AICcmodavg) # for Akaike information criterion (AIC) test for models

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
  write.csv("processed data for paper.csv", row.names = FALSE)                  #


# 2 Creating subsets of the data that I want to look at-------------------------
subset_list <- list(
  # Splitting into cortex and striatum
  "ctx" = data[data$region == "ctx", "percent_area_adjusted"],
  "cpu" = data[data$region == "cpu", "percent_area_adjusted"],
  
  # Splitting by diet (within cortex and striatum groups)
  "control_ctx" = 
    data[data$diet == "C/C" & data$region == "ctx", "percent_area_adjusted"],
  "test_ctx" = 
    data[data$diet == "HF/C" & data$region == "ctx", "percent_area_adjusted"],
  "control_cpu" = 
    data[data$diet == "C/C" & data$region == "cpu", "percent_area_adjusted"],
  "test_cpu" = 
    data[data$diet == "HF/C" & data$region == "cpu", "percent_area_adjusted"],
  
  # Splitting by sex (within cortex and striatum groups)
  "male_ctx" = 
    data[data$sex == "male" & data$region == "ctx", "percent_area_adjusted"],
  "female_ctx" = 
    data[data$sex == "female" & data$region == "ctx", "percent_area_adjusted"],
  "male_cpu" = 
    data[data$sex == "male" & data$region == "cpu", "percent_area_adjusted"],
  "female_cpu" = 
    data[data$sex == "female" & data$region == "cpu", "percent_area_adjusted"],
 
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
stats_df <- data.frame(matrix(ncol = 4, nrow = length(subset_list)))
colnames(stats_df) <- c("test", "mean", "se", "shapiro_result")
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
  
  if(shapiro$p.value>0.05){
    stats_df[[4]][i] <- "normal"
  } else{
    stats_df[[4]][i] <- "not normal"
  }
  
  colnames(subset_df) <- "percent" # (allows ggplot to plot varying columns)
  
  plot <- # Creating plots to display the distribution of the subset data
    ggplot(data = subset_df, mapping = aes(x = percent)) + 
    geom_density(colour = "cornflowerblue", 
                 fill = alpha("cornflowerblue", 0.3)) +
    theme_minimal() + 
    ylab("Density") + xlab("Area Covered by Microglia %")
  #ggsave(plot = plot, filename = paste0("Graphs/distribution density", 
  # name, ".png"), width = 6.25, height = 5)
}


# 4 One-way significance Testing ----------------------------------------------
# One way testing: diet, given region (non-parametric)
leveneTest(percent_area_adjusted ~ diet, 
           data = data[data$region == "ctx",]) # p = 0.2766 equal variance
wilcox.test(subset_list[["control_ctx"]], subset_list[["test"]], 
            alternative = "two.sided", exact = FALSE) # p-value = 0.002526 **

leveneTest(percent_area_adjusted ~ diet, 
           data = data[data$region == "cpu",]) # p = 0.7433 equal variance
wilcox.test(subset_list[["control_cpu"]], subset_list[["test"]], 
            alternative = "two.sided", exact = FALSE) # p-value = 0.002526 **

# One way testing: sex, given region (non-parametric)
leveneTest(percent_area_adjusted ~ sex, 
           data = data[data$region == "ctx",]) # p = 0.589 equal variance 
wilcox.test(subset_list[["female_ctx"]], subset_list[["male_ctx"]], 
            alternative = "two.sided", exact = FALSE) # p-value = 0.01529 *

leveneTest(percent_area_adjusted ~ sex, 
           data = data[data$region == "cpu",]) # p = 0.7737 equal variance 
wilcox.test(subset_list[["female_cpu"]], subset_list[["male_cpu"]], 
            alternative = "two.sided", exact = FALSE) # p-value = 0.02258 *

# One way testing: region (non-parametric)
leveneTest(percent_area_adjusted ~ region, 
           data = data) # p = 0.6869 equal variance
wilcox.test(subset_list[["ctx"]], subset_list[["cpu"]], 
            alternative = "two.sided", exact = FALSE) # p-value = 0.7649 NS


# 5 Two-way significance testing------------------------------------------------##########
# Two way testing: diet and sex, given region:
leveneTest(percent_area_adjusted ~ diet*sex, 
           data = data[data$region ==  "ctx",]) # p = 0.02996 * unequal variance
twoway_test_sex_ctx <- aov(percent_area_adjusted ~ diet+sex, 
                           data = data[data$region ==  "ctx",])
summary(twoway_test_sex_ctx)
        #             Df Sum Sq Mean Sq F value Pr(>F)   
        #diet         1  205.2  205.18   9.524 0.0056 **
        #sex          1  172.3  172.30   7.998 0.0101 * 
        #Residuals   21  452.4   21.54 
twoway_test_sex_ctx_interaction <- aov(percent_area_adjusted ~ diet*sex, 
                                       data = data[data$region ==  "ctx",])
summary(twoway_test_sex_ctx_interaction)
        #             Df Sum Sq Mean Sq F value  Pr(>F)   
        #diet         1  205.2  205.18   9.203 0.00656 **
        #sex          1  172.3  172.30   7.728 0.01156 * 
        #diet:sex     1    6.5    6.48   0.291 0.59563   
        #Residuals   20  445.9   22.30  

leveneTest(percent_area_adjusted ~ diet*sex, 
           data = data[data$region ==  "cpu",]) # p = 0.18 equal variance
twoway_test_sex_cpu <- aov(percent_area_adjusted ~ diet+sex, 
                           data = data[data$region ==  "cpu",])
summary(twoway_test_sex_cpu)
        #             Df Sum Sq Mean Sq F value Pr(>F)  
        #diet         1  139.8  139.78   4.162 0.0541 .
        #sex          1  219.9  219.94   6.548 0.0183 *
        #Residuals   21  705.3   33.59    
twoway_test_sex_cpu_interaction <- aov(percent_area_adjusted ~ diet*sex, 
                           data = data[data$region ==  "cpu",])
summary(twoway_test_sex_cpu_interaction)
        #             Df Sum Sq Mean Sq F value Pr(>F)  
        #diet         1  139.8  139.78   4.272 0.0519 .
        #sex          1  219.9  219.94   6.721 0.0174 *
        #diet:sex     1   50.9   50.87   1.555 0.2269  
        #Residuals   20  654.5   32.72    

# Two way testing: diet and region:
leveneTest(percent_area_adjusted ~ diet*region, 
           data = data) # p = 0.6656 equal variance
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

# Three way testing: diet, sex and region
leveneTest(percent_area_adjusted ~ diet*region*sex, 
           data = data) # p = 0.06312 equal variance
threeway_test <- aov(percent_area_adjusted ~ diet+region+sex, data = data)
summary(threeway_test)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
#diet         1  341.8   341.8  12.940 0.000810 ***
#region       1    8.1     8.1   0.307 0.582291    
#sex          1  390.8   390.8  14.793 0.000383 ***
#Residuals   44 1162.3    26.4  

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
model_set <- list(twoway_test_sex_ctx, twoway_test_sex_ctx_interaction,
                  twoway_test_sex_cpu, twoway_test_sex_cpu_interaction,
                  twoway_test_region, twoway_test_region_interaction, 
                  threeway_test, threeway_test_interaction)

model_names <- c("twoway_test_sex_ctx", "twoway_test_sex_ctx_interaction",
                 "twoway_test_sex_cpu", "twoway_test_sex_cpu_interaction",
                 "twoway_test_region", "twoway_test_region_interaction",
                 "threeway_test", "threeway_test_interaction")

aictab(model_set, modnames = model_names) 
#                               K   AICc Delta_AICc AICcWt Cum.Wt      LL
#twoway_test_sex_ctx             4 148.69       0.00   0.80   0.80  -69.29
#twoway_test_sex_ctx_interaction 5 151.57       2.88   0.19   0.99  -69.12
#twoway_test_sex_cpu             4 159.35      10.66   0.00   1.00  -74.62
#twoway_test_sex_cpu_interaction 5 160.78      12.09   0.00   1.00  -73.72
#threeway_test                   5 300.62     151.93   0.00   1.00 -144.60
#threeway_test_interaction       9 309.30     160.61   0.00   1.00 -143.28
#twoway_test_region              4 312.03     163.34   0.00   1.00 -151.55
#twoway_test_region_interaction  5 314.44     165.75   0.00   1.00 -151.50


# 6 Post-hoc test and correction of significance values-------------------------##########
data$sex_diet <- NA
data$region_diet <- NA
data$sex_region_diet <- NA

for(i in 1:nrow(data)){ # Adding extra column to allow multiple comparisons 
  data$sex_diet[[i]] <- paste0(data$sex[[i]], data$diet[[i]]) |>
    as.character()
  data$region_diet[[i]] <- paste0(data$region[[i]], data$diet[[i]]) |>
    as.character()
  data$sex_region_diet[[i]] <- paste0(data$sex[[i]], data$region[[i]], 
                                  data$diet[[i]]) |> 
    as.character()
}

# Two way testing: diet and sex, given region:
pairwise.wilcox.test(data[data$region == "ctx",]$percent_area_adjusted,  
                     data[data$region == "ctx",]$sex_diet, paired = FALSE,
                     p.adj = "bonferroni",
                     alternative = "two.sided", 
                     exact = FALSE)
#           femaleC/C femaleHF/C maleC/C
#femaleHF/C 1.00      -          -      
#maleC/C    0.03      0.03       -      
#maleHF/C   1.00      1.00       0.39

pairwise.wilcox.test(data[data$region == "cpu",]$percent_area_adjusted,  
                     data[data$region == "cpu",]$sex_diet, paired = FALSE,
                     p.adj = "bonferroni", alternative = "two.sided", 
                     exact = FALSE)
#           femaleC/C femaleHF/C maleC/C
#femaleHF/C 1.00      -          -      
#maleC/C    0.03      0.03       -      
#maleHF/C   1.00      1.00       0.39 

# Two way testing: diet and region:
pairwise.wilcox.test(data$percent_area_adjusted,  
                     data$region_diet, paired = FALSE,
                     p.adj = "bonferroni", alternative = "two.sided", 
                     exact = FALSE)
#       cpuC/C cpuHF/C ctxC/C
#cpuHF/C 0.414  -       -     
#ctxC/C  1.000  0.099   -     
#ctxHF/C 0.413  1.000   0.135 

# Three way testing: diet, sex and region
pairwise.wilcox.test(data$percent_area_adjusted,  
                     data$sex_region_diet, paired = FALSE,
                     p.adj = "bonferroni", alternative = "two.sided", 
                     exact = FALSE)
# Nothing was significant likely due to correction of many groups


# 6 Plots ----------------------------------------------------------------------
# One-way comparison -----
(# Diet (cortex)
  stats_df[3:4, ] |> 
  ggplot(aes(x = test, y = mean)) +
    geom_bar(stat = "identity", width = 0.5, fill = "cornflowerblue") +
    geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1)) +
    geom_signif(comparisons = list(c("control_ctx", "test_ctx")), 
                annotation = "**") + 
    scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12),
                       expand = expansion(mult = c(0, 0.05))) +
    ggtitle("Cortex") +
    xlab("Diet (Paternal/Offspring)") + ylab("Mean Microglia Coverage (%)") +
    scale_x_discrete(labels=c("C/C", "HF/C")) + 
    theme_bw()
  ) |>
  #ggsave(filename = "Graphs/oneway comparison/diet (cortex).png",              #
  #       width = 5, height = 5)

(# Diet (striatum)
  stats_df[5:6, ] |> 
    ggplot(aes(x = test, y = mean)) +
    geom_bar(stat = "identity", width = 0.5, fill = "cornflowerblue") +
    geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1)) +
    geom_signif(comparisons = list(c("control_cpu", "test_cpu")), 
                annotation = "**") + 
    scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12),
                       expand = expansion(mult = c(0, 0.05))) +
    ggtitle("Striatum") +
    xlab("Diet (Paternal/Offspring)") + ylab("Mean Microglia Coverage (%)") +
    scale_x_discrete(labels=c("C/C", "HF/C")) + 
    theme_bw()
) |>
  #ggsave(filename = "Graphs/oneway comparison/diet (striatum).png",              #
  #       width = 5, height = 5)
    
  (# Sex (cortex)
    stats_df[7:8, ] |> 
      ggplot(aes(x = test, y = mean)) +
      geom_bar(stat = "identity", width = 0.5, fill = "cornflowerblue") +
      geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1)) +
      geom_signif(comparisons = list(c("male_ctx", "female_ctx")), 
                  annotation = "*") + 
      scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14),
                         expand = expansion(mult = c(0, 0.05))) +
      ggtitle("Cortex") +
      xlab("Sex") + ylab("Mean Microglia Coverage (%)") +
      scale_x_discrete(labels=c("Female", "Male")) + 
      theme_bw()
  ) |>
  #ggsave(filename = "Graphs/oneway comparison/sex (cortex).png",              #
  #       width = 5, height = 5)
  
  (# Sex (striatum)
    stats_df[9:10, ] |> 
      ggplot(aes(x = test, y = mean)) +
      geom_bar(stat = "identity", width = 0.5, fill = "cornflowerblue") +
      geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1)) +
      geom_signif(comparisons = list(c("male_cpu", "female_cpu")), 
                  annotation = "*") + 
      scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14),
                         expand = expansion(mult = c(0, 0.05))) +
      ggtitle("Striatum") +
      xlab("Sex") + ylab("Mean Microglia Coverage (%)") +
      scale_x_discrete(labels=c("Female", "Male")) + 
      theme_bw()
  ) |>
  #ggsave(filename = "Graphs/oneway comparison/sex (striatum).png",              #
  #       width = 5, height = 5)
  
  (# Region
    stats_df[1:2, ] |> 
      ggplot(aes(x = test, y = mean)) +
      geom_bar(stat = "identity", width = 0.5, fill = "cornflowerblue") +
      geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1)) +
      scale_y_continuous(breaks = c(2, 4, 6, 8, 10),
                         expand = expansion(mult = c(0, 0.05))) +
      xlab("Region") + ylab("Mean Microglia Coverage (%)") +
      scale_x_discrete(labels=c("Striatium", "Cortex")) + 
      theme_bw()
  ) |>
  #ggsave(filename = "Graphs/oneway comparison/region.png",                     #
  #       width = 5, height = 5)  
  
# Two-way comparison -----
stats_df_twoway <- stats_df[11:18, ]
stats_df_twoway$sex <- c("male", "female", "male", "female", 
                         "male", "female", "male", "female")
stats_df_twoway$diet <- c("C/C", "C/C", "HF/C", "HF/C",
                          "C/C", "C/C", "HF/C", "HF/C")

(# Diet and sex (cortex) grouped by diet
  stats_df_twoway[1:4, ] |> 
    ggplot(aes(x = diet, y = mean, fill = sex)) +
    geom_bar(stat = "identity", position=position_dodge(), width = 0.5) +
    geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1), 
                  position=position_dodge(.5)) +
    guides(fill=guide_legend(title="Offspring Sex")) +
    scale_fill_manual(values = c("#61B499", "#8E61B4")) + 
    geom_signif(stat = "identity", aes(x = 0.8, xend = 1.2, y = 10, yend = 10, 
                annotation = "*")) + 
    scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16),
                       expand = expansion(mult = c(0, 0.05))) +
    ggtitle("Cortex") +
    xlab("Diet (Paternal/Offspring)") + ylab("Mean Microglia Coverage (%)") +
    scale_x_discrete(labels=c("C/C", "HF/C")) + 
    theme_bw()
) |>
  ggsave(filename = "Graphs/twoway comparison/diet and sex (cortex).png",              #
         width = 6.25, height = 5)
  
(# Diet and sex (cortex) grouped by sex
  stats_df_twoway[1:4, ] |> 
    ggplot(aes(x = sex, y = mean, fill = diet)) +
    geom_bar(stat = "identity", position=position_dodge(), width = 0.5) +
    geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1), 
                  position=position_dodge(.5)) +
    guides(fill=guide_legend(title="Diet (Paternal/Offspring)")) +
    scale_fill_manual(values = c("#6185B4", "#B461AE")) + 
    geom_signif(stat = "identity", aes(x = 1.8, xend = 2.2, y = 10, yend = 10, 
                                       annotation = "*")) + 
    scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16),
                       expand = expansion(mult = c(0, 0.05))) +
    ggtitle("Cortex") +
    xlab("Offspring Sex") + ylab("Mean Microglia Coverage (%)") +
    scale_x_discrete(labels=c("Female", "Male")) + 
    theme_bw()
) |>
  ggsave(filename = "Graphs/twoway comparison/diet and sex (cortex).png",              #
         width = 6.25, height = 5)

  
  
  


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
ggsave(plot = significant_box_plot, filename =                                  #
        "Graphs/significant box plot.png", width = 6.25, height = 5)

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
ggsave(plot = significant_box_plot_se, filename =                               #
       "Graphs/significant box plot se.png", width = 6.25, height = 5)

# Bar graph showing mean and standard error
significant__plot_df <- stats_df[7:10,]
significant__plot_df$sex <- c("male", "female", "male", "female")
significant__plot_df$diet <- c("C/C", "C/C", "HF/C", "HF/C")

significant_bar_plot_sex <- 
  significant__plot_df |>
  ggplot(aes(x = diet, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position=position_dodge(), width = 0.5) +
  geom_signif(stat = "identity",
              aes(x = 0.8, xend = 1.2, y = 12, yend = 12, annotation = "***")) + 
    #guides(fill=guide_legend(title="Offspring Sex")) +
  scale_fill_manual(values = c("#61B499", "#8E61B4")) + 
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1), 
                position=position_dodge(.5)) + 
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14), 
                     expand = expansion(mult = c(0, 0.05))) + 
  #xlab("Paternal Diet") + ylab("Mean Microglia Coverage (%)")+
  theme_bw()
ggsave(plot = significant_bar_plot_sex, filename =                              #
        "Graphs/significant bar plot sex.png", width = 6.25, height = 5)

significant_bar_plot_diet <- 
  significant__plot_df |>
  ggplot(aes(x = sex, y = mean, fill = diet)) +
  geom_bar(stat = "identity", position=position_dodge(), width = 0.5) +
  geom_signif(stat = "identity",
              aes(x = 1.8, xend = 2.2, y =11, yend = 11, annotation = "*")) + 
  guides(fill=guide_legend(title="Paternal Diet")) +
  scale_fill_manual(values = c("#6185B4", "#B461AE")) + 
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se, width = 0.1), 
                position=position_dodge(.5)) + 
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14), 
                     expand = expansion(mult = c(0, 0.05))) + 
  xlab("Offspring Sex") + ylab("Mean Microglia Coverage (%)")+
  theme_bw()
ggsave(plot = significant_bar_plot_sex, filename =                              #
         "Graphs/significant bar plot sex.png", width = 6.25, height = 5)


# Arranging two plots in the same image-----
# Sex and Region
oneway_sex_plot <- oneway_sex_plot + theme(axis.title.y = element_blank())
oneway_region_plot <- oneway_region_plot + theme(axis.title.y = element_blank())

combined_oneway <- grid.arrange(oneway_diet_plot, oneway_sex_plot, 
                                oneway_region_plot, ncol = 3, nrow = 1)
ggsave(plot = combined_oneway, filename =                                       #
      "Graphs/significant oneway plot.png", width = 10, height = 5)

# Sex and diet together
combined_twoway <- grid.arrange(
  significant_bar_plot_diet, significant_bar_plot_sex, ncol = 2, nrow = 1)
ggsave(plot = combined_twoway, filename =                                       #
       "Graphs/significant combined plot.png", width = 10, height = 5)
