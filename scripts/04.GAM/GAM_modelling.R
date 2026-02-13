## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling")
packages <- c("dplyr", "ggplot2", "here", "nlme", "mgcv", "MuMIn", "gratia", "rsample", "gamclass", "purrr", "broom", "mgcViz", "Metrics")
lapply(packages, library, character.only = TRUE)




## ----GAM modelling------------------------------------------------------------
## import and transform
df_GAM <- read.csv(here::here("data/df_GAM.csv"), header = TRUE, sep = ",")

df_GAM$node <- as.factor(df_GAM$node)
df_GAM$id <- as.factor(df_GAM$id)
df_GAM$subcortical <- as.factor(df_GAM$subcortical)
df_GAM$sex <- as.factor(df_GAM$sex)
df_GAM$site <- as.factor(df_GAM$site)

## model fitting and comparison
model_gam <- bam(atrophy_rate ~ s(atrophy_bl) + s(sc_centrality) + s(fc_centrality) + s(disc) + s(neighbourhood_atrophy_sc) + s(neighbourhood_atrophy_fc) + s(neighbourhood_atrophy_cge) + s(x,y,z) + s(age) + s(fu_time_duration) + sex + subcortical + site + s(id, bs = 're') + s(node, bs = 're') + s(atrophy_bl, node, bs = 're'), data = df_GAM, na.action = na.fail, method = "fREML", select = TRUE, discrete = TRUE,nthreads = 8)

model_gam_bl <- bam(atrophy_rate ~ s(atrophy_bl) + s(x,y,z) + s(age) + s(fu_time_duration) + sex + subcortical + site + s(id, bs = 're') + s(node, bs = 're') + s(atrophy_bl, node, bs = 're'), data = df_GAM, na.action = na.fail, method = "fREML", select = TRUE, discrete = TRUE, nthreads = 8)

anova_test <- anova(model_gam_bl, model_gam, test = "LRT")
AIC_test <- AIC(model_gam_bl, model_gam)

# save
save(model_gam, model_gam_bl, anova_test, AIC_test, file = here("data/GAM_modelling.RData"))

## LOOCV
custom_cvparts <- group_vfold_cv(df_GAM, group = id)

mod_form <- as.formula(atrophy_rate ~ s(atrophy_bl) + s(sc_centrality) + s(fc_centrality) + s(disc) + s(neighbourhood_atrophy_sc) + s(neighbourhood_atrophy_fc) + s(neighbourhood_atrophy_cge) + s(x,y,z) + s(age) + s(fu_time_duration) + sex + subcortical + site + s(id, bs = 're') + s(node, bs = 're') + s(atrophy_bl, node, bs = 're')) 

mod_form_bl <- as.formula(atrophy_rate ~ s(atrophy_bl) + s(x,y,z) + s(age) + s(fu_time_duration) + sex + subcortical + site + sex + subcortical + s(id, bs = 're') + s(node, bs = 're') + s(atrophy_bl, node, bs = 're')) 
  
holdout_results <- function(splits, ...){
  # fit the model
  mod <- bam(..., data = analysis(splits), na.action = na.fail, method = "fREML", select = TRUE, discrete = TRUE, nthreads = 8)
  # apply the model and assess performance
  holdout <- assessment(splits)
  res <- broom::augment(mod, newdata = holdout)
  performance <- data.frame(matrix(nrow = 1, ncol = 3))
colnames(performance) <- c("pearsonr", "RMSE", "MAE")
  performance$pearsonr <- cor(res$.fitted, holdout$atrophy_rate, method = "pearson")
  performance$RMSE <- rmse(holdout$atrophy_rate, res$.fitted)
  performance$MAE <- mae(holdout$atrophy_rate, res$.fitted)
  output <- list(res, performance)
  # return output
  return(output)
}

# example <- holdout_results(custom_cvparts$splits[[1]], mod_form)

LOOCV <- purrr::map(custom_cvparts$splits, holdout_results, mod_form)

LOOCV_bl <- purrr::map(custom_cvparts$splits, holdout_results, mod_form_bl)

# save
save(LOOCV, LOOCV_bl, file = here("data/GAM_LOOCV.RData"))

## LOOCV CVgam
# custom_cvparts <- rep(1:565, each = 114)
# 
# LOOCV <- CVgam(formula = atrophy_rate ~ s(atrophy_bl) + s(sc_centrality) + s(fc_centrality) + s(disc) + s(neighbourhood_atrophy_sc) + s(neighbourhood_atrophy_fc) + s(neighbourhood_atrophy_cge) + s(x,y,z) + s(age) + sex + subcortical + s(id, bs = 're') + s(node, bs = 're') + s(atrophy_bl, node, bs = 're'), data = df_GAM, printit = TRUE, method = "GCV.Cp", cvparts = custom_cvparts, seed = 42)

# LOOCV_bl <- CVgam(formula = atrophy_rate ~ s(atrophy_bl) + s(x,y,z) + s(age) + sex + subcortical + s(id, bs = 're') + s(node, bs = 're') + s(atrophy_bl, node, bs = 're'), data = df_GAM, printit = TRUE, method = "GCV.Cp", cvparts = custom_cvparts, seed = 42)
# 
# # save
# save(LOOCV_bl, file = here("data/LOOCV_bl.RData"))

## model selection with MuMIn
# model_selection <- dredge(model_gam, evaluate = TRUE, fixed = ~ s(atrophy_bl) + s(x,y,z) + s(age) + sex + subcortical + s(id, bs = 're') + s(node, bs = 're') + s(node, atrophy_bl, bs = 're'))

# best_model <- get.models(model_selection, 1)[[1]]


