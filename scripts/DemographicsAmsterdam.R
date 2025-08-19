#===============================================================================
# demographics and statistical analyses
#=============================================================================== 

library(table1)
path="/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/data"
setwd(path)

df_cross <- readRDS("df_baseline.Rds")
df_cross$disease_duration[df_cross$disease_duration > 100] <- NA #remove above 100 y
#df_cross$brain_seg_not_vent <- df_cross$brain_seg_not_vent / 10000000000 
table1= table1(~ sex + age + disease_duration  + edss + clinical_phenotype + tlv + brain_seg_not_vent | group, data=df_cross )
table1

## stages 
df_stages <- read.csv("stages.csv")
table2= table1(~ stage | stage, data=df_stages )
table2

df_long <- read.csv("data_atrophy_long.csv")
df_long_graz <- subset(df_long, Site == "graz")
table1= table1(~ Sex + Age  + EDSS | Group, data=df_long_graz )
table1

#####

  ## import atrophy and disc data 
atrophydir <- file.path(path, "atrophy_cross")
atrophy_cross_group_100 <- read.csv(file.path(atrophydir, "PrograMS_atrophy_baseline_SchaeferSubcortex114Parcels_group_cohend.csv"))
atrophy_cross_group_400 <- read.csv(file.path(atrophydir, "PrograMS_atrophy_baseline_SchaeferSubcortex414Parcels_group_cohend.csv"))

discdir <- file.path(path, "disc")

disc_group_100_ctx <- read.csv(file.path(discdir, "PrograMS_disc_SchaeferSubcortex114Parcels_group_ctx.csv"), header = FALSE)
disc_group_100_sctx <- read.csv(file.path(discdir, "PrograMS_disc_SchaeferSubcortex114Parcels_group_sctx.csv"), header = FALSE)

disc_group_400_ctx <- read.csv(file.path(discdir, "PrograMS_disc_SchaeferSubcortex414Parcels_group_ctx.csv"), header = FALSE)
disc_group_400_sctx <- read.csv(file.path(discdir, "PrograMS_disc_SchaeferSubcortex414Parcels_group_sctx.csv"), header = FALSE)

compute_dc <- function(fc_ctx, sc_ctx, fc_sctx, sc_sctx) {
  # Computes degree centrality for fMRI & dMRI
  fc_ctx_dc <- colSums(fc_ctx)  #  (axis = 0)
  sc_ctx_dc <- colSums(sc_ctx)  #  (axis = 0)
  fc_sctx_dc <- rowSums(fc_sctx)  # (axis = 1)
  sc_sctx_dc <- rowSums(sc_sctx)  #  axis = 1)
  
  return(list(fc_ctx_dc = fc_ctx_dc, sc_ctx_dc = sc_ctx_dc, 
              fc_sctx_dc = fc_sctx_dc, sc_sctx_dc = sc_sctx_dc))
}

dc_results <-compute_dc(disc_group_100_ctx, disc_group_400_ctx, disc_group_100_sctx, disc_group_400_sctx)

disc_ctx_dc_100 <- dc_results$fc_ctx_dc
disc_ctx_dc_400 <- dc_results$sc_ctx_dc  
disc_sctx_dc_100 <- dc_results$fc_sctx_dc
disc_sctx_dc_400 <- dc_results$sc_sctx_dc


