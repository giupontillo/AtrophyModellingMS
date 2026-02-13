#===============================================================================
# demographics and statistical analyses
#=============================================================================== 
packages <- c(  "dplyr", "tidyr", "purrr", "stringr", "lubridate",  # data manipulation
                "stats", "rstatix", "coin",                         # statistics
                "ggplot2", "patchwork", "gghalves",                 # plotting
                "table1"                                            # summary tables
)

purrr::walk(packages, library, character.only = TRUE)

path="/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/data"
setwd(path)

#### create demographic, clinical and mri table 
my.render.cont <- function(x) {
  s <- summary(x, na.rm=TRUE)
  # compute median, Q1, Q3
  med <- median(x, na.rm=TRUE)
  q1  <- quantile(x, 0.25, na.rm=TRUE)
  q3  <- quantile(x, 0.75, na.rm=TRUE)
  
  # return string like "median (Q1 – Q3)"
  c("Median (IQR)" = sprintf("%.1f (%.1f – %.1f)", med, q1, q3))
}

df_cross <- readRDS("df_baseline.Rds")
df_cross$disease_duration[df_cross$disease_duration > 100] <- NA #remove above 100 y
df_cross <- df_cross %>%
  mutate(dmt_bin = ifelse(!is.na(dmt) & dmt != "" & dmt != "placebo", "yes", "no"))
df_cross$brain_seg_not_vent <- df_cross$brain_seg_not_vent / 10000000000 

#df_cross <- df_cross[grepl("amsterdam", df_cross$site, ignore.case = TRUE), ] #only amsterdam data
table1= table1(~ sex + age + disease_duration  + edss + clinical_phenotype + dmt_bin+ tlv + brain_seg_not_vent | group, data=df_cross, 
               render.continuous = my.render.cont
               )
table1

### test group difference demographic variables
perm_age <- oneway_test(age ~ group, data = df_cross, distribution = approximate(nresample = 10000))
perm_brain <- oneway_test(brain_seg_not_vent ~ group, data = df_cross, distribution = approximate(nresample = 10000))
sex_table <- table(df_cross$sex, df_cross$group)
chisq_result <- chisq.test(sex_table)

perm_age
perm_brain
chisq_result

## stages 
df_stages <- read.csv("stages.csv")
table2= table1(~ stage | stage, data=df_stages )
table2

#check number of session and FU time (mean + range)
df_long <- read.csv("data_atrophy_long.csv")
df_long_fil <- df_long[!(df_long$Session %in% c("ses-01", "ses-Y05") | df_long$Excluded == 1), ]
#df_long_graz <- subset(df_long, Site == "graz")
table2= table1(~ Sex + Age  + EDSS + FU_time | Group, data=df_long_fil )
table2

#check min and max years MRI 
df_year <- df_long %>% 
  mutate(year = parse_date_time(DateMRI, orders=c("dmy","my"), 
                                truncated = 1) %>%
           year())

df_year %>%
  summarise( 
    min_year <- min(as.numeric(year), na.rm = TRUE), 
    max_year <- max(as.numeric(year), na.rm = TRUE))

#####

##plot euler number distribution and cutt-off
df_euler <- df_long %>%
  select(mean_euler_bh) %>%
  filter(is.finite(mean_euler_bh))

binwidth <- 5
hist_counts <- hist(df_euler$mean_euler_bh, breaks = seq(floor(min(df_euler$mean_euler_bh)), ceiling(max(df_euler$mean_euler_bh)), by = binwidth), plot = FALSE)
ymax <- max(hist_counts$counts) * 1.05  # 5% higher than max

# Plot smooth density lines for both hemispheres
euler_dist <- ggplot(df_euler, aes(x = mean_euler_bh)) + 
  geom_histogram(aes(y = after_stat(count)), binwidth = binwidth, fill = "lightblue", color = "black") +  
  geom_density(aes(y = after_stat(count*(binwidth))),color  = "steelblue", linewidth = 1) + 
  geom_vline(xintercept = -120, color = "black", linetype = "dashed", linewidth = 0.5) +
  labs(
    title = "",
    x = "Euler Number",
    y = "Count"
  ) +
  theme_minimal(base_size = 30) +
  coord_cartesian(ylim=c(0, ymax), expand=FALSE) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.3 ),
    axis.line = element_line(colour = "black", linewidth =  0.3),
    legend.title = element_blank(),
    axis.text = element_text(size = 20),                                        
    axis.title = element_text(size = 20),  
    plot.title = element_text(hjust = 0.5), 
    panel.grid = element_blank()
  )

#Save image 
ggsave("../scripts/figures/euler_distribution.png", plot = euler_dist, width = 10, height = 5, dpi = 300)  


#####

#perform group comparisons
df_rvals <- read.csv("/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/output/june/cross_indv_results/2025.08.21_all_indv_rvals.csv",
                     header = TRUE, stringsAsFactors = FALSE) 

#check democlinical diff between assigment_2 (hydra group)
df_clean <- df_rvals[!duplicated(df_rvals$id), c("id", "assignment_2")]
df_merged <- merge(df_cross, df_clean, by ="id", all = TRUE)
df_merged <- df_merged[!is.na(df_merged$assignment_2),]
df_merged$assignment_2 <-factor(df_merged$assignment_2)
table3= table1(~ sex + age + disease_duration  + edss + clinical_phenotype + dmt_bin+ tlv + brain_seg_not_vent | assignment_2, data=df_merged, 
               #render.continuous = my.render.cont #uncomment for tlv
)
table3

### test group difference demographic variables
perm_age2 <- oneway_test(age ~ assignment_2, data = df_merged, distribution = approximate(nresample = 10000))
perm_brain2 <- oneway_test(brain_seg_not_vent ~ assignment_2, data = df_merged, distribution = approximate(nresample = 10000))
perm_sex2 <- independence_test(sex ~ assignment_2, data = df_merged, distribution = approximate(nresample = 10000))

# Print results
perm_age2
perm_brain2
perm_sex2

ggplot(df_merged, aes(x = brain_seg_not_vent, fill = factor(assignment_2))) +
  geom_histogram(binwidth = 50, alpha = 0.5, position = "identity", color = "black") +
  labs(title = "Age Distribution by Group", x = "Age", y = "Frequency", fill = "Group")

###


df_rvals <- df_rvals %>%
  mutate(category = recode(category,
                           "functional cortical hubs" = "functional cortico-cortical degree centrality",
                           "structural cortical hubs" = "structural cortico-cortical degree centrality",
                           "functional subcortical hubs" = "functional subcortico-cortical degree centrality",
                           "structural subcortical hubs" = "structural subcortico-cortical degree centrality",
                           "functional cortical neighbour atrophy" = "functional cortico-cortical neighbour atrophy",
                           "structural cortical neighbour atrophy" = "structural cortico-cortical neighbour atrophy",
                           "functional subcortical neighbour atrophy" = "functional subcortico-cortical neighbour atrophy",
                           "structural subcortical neighbour atrophy" = "structural subcortico-cortical neighbour atrophy",
                           "structural cortical disconnectome" = "structural cortico-cortical disconnection",
                           "structural subcortical disconnectome" = "structural subcortico-cortical disconnection", 
                           "cortical gene expression" = "cortico-cortical transcriptomic neighbour atrophy",
                           "subcortical gene expression" = "subcortico-cortical transcriptomic neighbour atrophy"))



categories <- c("functional cortico-cortical degree centrality",
                "structural cortico-cortical degree centrality",
                "functional subcortico-cortical degree centrality",
                "structural subcortico-cortical degree centrality",
                "functional cortico-cortical neighbour atrophy",
                "structural cortico-cortical neighbour atrophy",
                "functional subcortico-cortical neighbour atrophy",
                 "structural subcortico-cortical neighbour atrophy",
                "structural cortico-cortical disconnection",
                "structural subcortico-cortical disconnection", 
                "cortico-cortical transcriptomic neighbour atrophy",
                "subcortico-cortical transcriptomic neighbour atrophy")


# Filter dataset to relevant categories only
df_sub <- df_rvals %>%
  filter(category %in% categories)

# -------------------------------
# One-sample t-tests against 0
# -------------------------------
cohen_d_results <- df_sub %>%
  group_by(category) %>%
  cohens_d(z_value ~ 1, mu = 0, ci = TRUE)

one_sample_results <- df_sub %>%
  group_by(category) %>%
  t_test(z_value ~ 1, mu = 0, conf.level = 0.95) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  left_join(cohen_d_results, by="category")

# -------------------------------
# Two-sample Welch t-tests by assignment_2 (hydra)
# -------------------------------
two_sample_results <- df_sub %>%
  group_by(category) %>%
  t_test(z_value ~ assignment_2, var.equal = FALSE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()

# -------------------------------
# Combine results
# -------------------------------
df_stats <- list(
  OneSample = one_sample_results,
  TwoSample = two_sample_results
)

# Inspect
df_stats$OneSample
df_stats$TwoSample

final_table <- one_sample_results %>%
  rename(
    OneSample_t = statistic,
    OneSample_df = df,
    OneSample_p = p,
    OneSample_p.adj = p.adj,
    OneSample_sig = p.adj.signif, 
    OneSample_cohens_d = effsize
  ) %>%
  left_join(
    two_sample_results %>%
      rename(
        TwoSample_group1 = group1,
        TwoSample_group2 = group2,
        TwoSample_t = statistic,
        TwoSample_df = df,
        TwoSample_p = p,
        TwoSample_p.adj = p.adj,
        TwoSample_sig = p.adj.signif
      ),
    by = "category"
  )

# View results
View(final_table)
write.csv(final_table, "one_sample_two_sample_results_09092025.csv", row.names = FALSE)

# -------------------------------
# Violin plot for atrophy HYDRA (all subjects pooled)
# -------------------------------
df_all <- df_rvals %>%
  filter(!is.na(z_value), is.finite(z_value))  %>%
  mutate(category = factor(category, levels= categories))

wrapped_labels <- str_wrap(as.character(df_all$category), width = 20)
df_all$category_wrapped <- factor(wrapped_labels, levels = str_wrap(categories, width = 20))

effect_labels <- one_sample_results %>%
  as.data.frame() %>%
  select(category, magnitude) %>%
  distinct %>%
  mutate(category_wrapped = str_wrap(category, width = 20), 
         magnitude = as.character(magnitude),
         stars = case_when(
           magnitude == "small" ~ "*", 
           magnitude == "moderate" ~ "**", 
           magnitude == "large" ~ "***", 
           TRUE ~ ""
         ))

p <- ggplot(df_all, aes(x = category_wrapped, y = z_value, fill = category)) +
  geom_half_violin(side = "r", trim = FALSE, alpha = 0.6, color = "black") +
  geom_half_point(side = "l", width = 0.8, size = 0.7, alpha = 0.4, aes(color = factor(category))) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
  geom_text(data = effect_labels, aes(x=category_wrapped, y = 1.8, label = stars), 
            inherit.aes = FALSE, size =5, fontface ="bold") +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(title = "", 
       y = "Fisher's z", x = NULL) +
  theme_classic(base_size = 10) +
  theme(
    text= element_text(family= "Arial"),
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5, size=12, face = "bold"),
    axis.text.y = element_text(size=12, face = "bold"), 
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5)
  )

print(p)

ggsave(
  filename = "all_mechanisms.png", plot = p, width = 20, height = 4,  dpi = 300)


# -------------------------------
# plot violin for atrophy hydra 
# -------------------------------
plots <- lapply(categories, function(cat) {
  df_cat <- df_rvals %>% 
    filter(category == cat, !is.na(z_value), is.finite(z_value))
  
  # Wrap title
  cat_wrapped <- stringr::str_wrap(cat, width = 20)
  
  sig_label <- final_table %>%
    filter(category == cat) %>%
    pull(TwoSample_sig)
  
  df_cat <- df_cat %>%
    mutate(assignment_2_f = factor(assignment_2, levels = c(1, 2), labels = c("HYDRA1", "HYDRA2")))
  
  p <- ggplot(df_cat, aes(x = factor(assignment_2_f), y = z_value, fill = factor(assignment_2_f))) +
    geom_half_violin(side = "r", alpha = 0.6, color = "black") +
    geom_half_point(side = "l", width = 0.8, alpha = 0.7, aes(color = factor(assignment_2_f))) +
    stat_summary(fun.data = function(x) {
      r <- quantile(x, probs = c(0.25, 0.5, 0.75))
      data.frame(ymin = r[1], y = r[2], ymax = r[3])
    }, geom = "crossbar", width = 0.1, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
    coord_cartesian(ylim = c(-2, 2)) +
    labs(title = cat_wrapped, y = "Fisher's z", x = NULL) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      text= element_text(family= "Arial"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size=12, face = "bold"),
      axis.text.y = element_text(size=12, face = "bold")
    )
  
  if (!is.na(sig_label) && sig_label != "ns") {
    y_max <- max(df_cat$z_value, na.rm = TRUE)
    y_line <- y_max + 0.2
    y_text <- y_max + 0.25
    
    p <- p +
      geom_segment(aes(x = 1, xend = 2, y = y_line, yend = y_line), color = "black", size = 0.3) +
      geom_text(aes(x = 1.5, y = y_text, label = sig_label), size = 4.5)
  }
  
  return(p)
})

n_cols <- 6
plot_grid <- wrap_plots(plots, ncol = n_cols)

plot_grid
ggsave(
  filename = "hydra_plots.png", plot = plot_grid, width = 20, height = 6,  dpi = 300)


# -------------------------------
# ANOVA by stage
# -------------------------------
for(cat in categories) {
  cat("\n=== Category:", cat, "===\n")
  df_cat <- df_rvals[df_rvals$category == cat, ]
  z_clean <- df_cat$z_value[is.finite(df_cat$z_value)]
  
  # KS test (against normal distribution)
  if(length(z_clean) > 3) {
    ks_res <- ks.test(z_clean, "pnorm", mean = mean(z_clean), sd = sd(z_clean))
    print(ks_res)
  }
  
  # Shapiro-Wilk per stage within this category (only if group size <= 5000)
  sw_res <- by(df_cat$z_value, df_cat$stage, function(x) {
    x_clean <- x[is.finite(x)]
    if(length(x_clean) >= 3 & length(x_clean) <= 5000) {
      shapiro.test(x_clean)
    } else {
      NA
    }
  })
  print(sw_res)
  
  # Histogram per stage
  p_hist <- ggplot(df_cat, aes(x = z_value)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    facet_wrap(~stage) +
    ggtitle(paste("Histogram -", cat))
  print(p_hist)
  
  # Q-Q plot per stage
  p_qq <- ggplot(df_cat, aes(sample = z_value)) +
    stat_qq() +
    stat_qq_line() +
    facet_wrap(~stage) +
    ggtitle(paste("Q-Q plot -", cat))
  print(p_qq)
}

# Initialize empty data frames
anova_table <- data.frame()
tukey_table <- data.frame()

for(cat in categories) {
  
  df_cat <- df_rvals[df_rvals$category == cat, ]
  if(length(unique(df_cat$stage)) > 1) {
    
    anova_res <- aov(z_value ~ stage, data = df_cat)
    anova_summary <- summary(anova_res)[[1]]
    
    anova_df <- data.frame(
      Category = cat,
      Term = rownames(anova_summary),
      Df = anova_summary$Df,
      SumSq = anova_summary$`Sum Sq`,
      MeanSq = anova_summary$`Mean Sq`,
      F_value = anova_summary$`F value`,
      P_value = anova_summary$`Pr(>F)`
    )
    
    anova_table <- bind_rows(anova_table, anova_df)
    
    tukey_res <- TukeyHSD(anova_res)
    tukey_df <- as.data.frame(tukey_res$stage)
    tukey_df$Comparison <- rownames(tukey_res$stage)
    tukey_df$Category <- cat
    tukey_df <- tukey_df %>% select(Category, Comparison, everything())
    
    tukey_table <- bind_rows(tukey_table, tukey_df)
  }
}

p_to_sig <- function(p) {
  if(is.na(p)) return(NA)
  else if(p < 0.001) return("***")
  else if(p < 0.01) return("**")
  else if(p < 0.05) return("*")
  else return("ns")
}

anova_table$Significance <- sapply(anova_table$P_value, p_to_sig)
tukey_table$Significance <- sapply(tukey_table$`p adj`, p_to_sig)

# View the tables
anova_table
tukey_table


# -------------------------------
# Plot results 2
# -------------------------------
library(tibble)

stage_levels <- c("ER", "LR", "PP", "SP")
stage_colors <- RColorBrewer::brewer.pal(n = length(stage_levels), name = "Set2")

plots <- lapply(categories, function(cat) {
  df_cat <- df_rvals %>% 
    filter(category == cat, !is.na(z_value), is.finite(z_value))
  
  df_cat$stage <- factor(df_cat$stage, levels = stage_levels)
  cat_wrapped <- stringr::str_wrap(cat, width = 20)
  
  tukey_df <- tukey_table %>% filter(Category == cat)
  tukey_df <- tibble::rownames_to_column(tukey_df, var ="comparison")
  
  p <- ggplot(df_cat, aes(x = stage, y = z_value, fill = stage)) +
    geom_half_violin(side = "r", alpha = 0.6, color = "black") +
    geom_half_point(side = "l", size = 0.8, alpha = 0.7, aes(color = stage)) +
    stat_summary(fun.data = function(x) {
      r <- quantile(x, probs = c(0.25, 0.5, 0.75))
      data.frame(ymin = r[1], y = r[2], ymax = r[3])
    }, geom = "crossbar", width = 0.1, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
    coord_cartesian(ylim = c(-2, 2)) +
    labs(title = cat_wrapped, y = "Fisher's z", x = NULL) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      text= element_text(family= "Arial"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size=12, face = "bold"),
      axis.text.y = element_text(size=12, face = "bold")
    )
  
  if (!is.null(tukey_df)) {
    
    tukey_df <- tukey_df %>%
      mutate(
        x1 = sapply(strsplit(comparison, "-"), function(x) which(stage_levels == x[1])),
        x2 = sapply(strsplit(comparison, "-"), function(x) which(stage_levels == x[2])),
        span = abs(x2 - x1)
      ) %>%
      arrange(desc(span))   
    
    y_base <- max(df_cat$z_value, na.rm = TRUE) + 0.2  
    offset <- 0.2
    n_lines <- 0
    
    for (i in 1:nrow(tukey_df)) {
      label <- tukey_df$Significance[i]
      
      if (label != "ns") {
        # positions of the groups
        x1 <- tukey_df$x1[i]
        x2 <- tukey_df$x2[i]
        
        if (length(x1) == 1 && length(x2) == 1) {
          y_line <- y_base + n_lines * offset
          y_text <- y_line + 0.05
        
          p <- p +
            geom_segment(x = x1, xend = x2, y = y_line, yend = y_line, size = 0.3) +
            geom_text(x = (x1 + x2)/2, y = y_text, label = label, size = 4)
          
          n_lines <- n_lines + 1 
          }
      }
    }
  }
  
  return(p)
})

n_cols <- 4
plot_grid <- wrap_plots(plots, ncol = n_cols)

plot_grid
ggsave(
  filename = "clinical_phenotypes_plots.png", plot = plot_grid, width = 20, height = 12,  dpi = 300)

