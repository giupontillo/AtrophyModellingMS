## ----setup, include=FALSE----------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling")
system("export LD_LIBRARY_PATH=/opt/aumc-apps-eb/software/Anaconda3/2024.02-1/lib")
packages <- c("dplyr", "ggplot2", "here", "corrplot", "nlme", "psych", "effsize", "janitor", "VertexWiseR", "R.matlab", "reticulate","rstudioapi")
lapply(packages, library, character.only = TRUE)
use_condaenv("atrophymodelling")
# openProject(path = "/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/AtrophyModelling.Rproj", newSession = FALSE)

## ----compute baseline neighbour atrophy weighted by SC/FC/CGE----------
### SchaeferSubcortex 114 parcels
## load normative functional/structural/CGE matrices
sc_ctx <- as.matrix(read.csv(here::here("data/hcp_connectivity/strucMatrix_ctx_schaefer_100.csv"), header = FALSE, sep = ","))
sc_sctx <- as.matrix(read.csv(here::here("data/hcp_connectivity/strucMatrix_sctx_schaefer_100.csv"), header = FALSE, sep = ","))
sc_sctx <- sc_sctx[,1:100] # the sctx matrix contains also subcortico-subcortical connections
fc_ctx <- as.matrix(read.csv(here::here("data/hcp_connectivity/funcMatrix_ctx_schaefer_100.csv"), header = FALSE, sep = ","))
fc_sctx <- as.matrix(read.csv(here::here("data/hcp_connectivity/funcMatrix_sctx_schaefer_100.csv"), header = FALSE, sep = ","))
cge_ctx <- as.matrix(read.csv(here::here("data/ge/CGE_SchaeferSubcortex114Parcels_ctx.csv"), header = FALSE, sep = ","))
cge_sctx <- as.matrix(read.csv(here::here("data/ge/CGE_SchaeferSubcortex114Parcels_sctx.csv"), header = FALSE, sep = ","))
## group-level
# load atrophy map
atrophy_mapping <- read.csv(here::here("data/atrophy_cross/atrophy_baseline_SchaeferSubcortex114Parcels_group_cohend.csv"), header = TRUE, sep = ",")
atrophy_ctx <- atrophy_mapping[,1:100]
atrophy_sctx <- atrophy_mapping[,101:114]
# compute neighbourhood atrophy sc ctx
neighbourhood_atrophy_sc_ctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_ctx)))
colnames(neighbourhood_atrophy_sc_ctx) <- paste0(colnames(atrophy_ctx), "_neighbours")
for (i in 1:length(atrophy_ctx)){
  neighbours <- which(sc_ctx[i,]!=0)
  weight_vec <- sc_ctx[i,c(neighbours)]
  neighbourhood_atrophy_sc_ctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy sc sctx
neighbourhood_atrophy_sc_sctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_sctx)))
colnames(neighbourhood_atrophy_sc_sctx) <- paste0(colnames(atrophy_sctx), "_neighbours")
for (i in 1:length(atrophy_sctx)){
  neighbours <- which(sc_sctx[i,]!=0)
  weight_vec <- sc_sctx[i,c(neighbours)]
  neighbourhood_atrophy_sc_sctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy fc ctx
neighbourhood_atrophy_fc_ctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_ctx)))
colnames(neighbourhood_atrophy_fc_ctx) <- paste0(colnames(atrophy_ctx), "_neighbours")
for (i in 1:length(atrophy_ctx)){
  neighbours <- which(fc_ctx[i,]!=0)
  weight_vec <- fc_ctx[i,c(neighbours)]
  neighbourhood_atrophy_fc_ctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy fc sctx
neighbourhood_atrophy_fc_sctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_sctx)))
colnames(neighbourhood_atrophy_fc_sctx) <- paste0(colnames(atrophy_sctx), "_neighbours")
for (i in 1:length(atrophy_sctx)){
  neighbours <- which(fc_sctx[i,]!=0)
  weight_vec <- fc_sctx[i,c(neighbours)]
  neighbourhood_atrophy_fc_sctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy cge ctx
neighbourhood_atrophy_cge_ctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_ctx)))
colnames(neighbourhood_atrophy_cge_ctx) <- paste0(colnames(atrophy_ctx), "_neighbours")
for (i in 1:length(atrophy_ctx)){
  neighbours <- which(cge_ctx[i,]!=0 & !is.na(cge_ctx[i,]))
  weight_vec <- cge_ctx[i,c(neighbours)]
  neighbourhood_atrophy_cge_ctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy cge sctx
neighbourhood_atrophy_cge_sctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_sctx)))
colnames(neighbourhood_atrophy_cge_sctx) <- paste0(colnames(atrophy_sctx), "_neighbours")
for (i in 1:length(atrophy_sctx)){
  neighbours <- which(cge_sctx[i,]!=0 & !is.na(cge_sctx[i,]))
  weight_vec <- cge_sctx[i,c(neighbours)]
  neighbourhood_atrophy_cge_sctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# write out
write.csv(neighbourhood_atrophy_sc_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_sc_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_sc_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_sc_sctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_fc_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_fc_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_fc_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_fc_sctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_cge_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_cge_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_cge_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_cge_sctx.csv"), row.names = FALSE)
## group-level long
# load atrophy map
atrophy_mapping <- read.csv(here::here("data/atrophy_cross/atrophy_baseline_SchaeferSubcortex114Parcels_group_long_cohend.csv"), header = TRUE, sep = ",")
atrophy_ctx <- atrophy_mapping[,1:100]
atrophy_sctx <- atrophy_mapping[,101:114]
# compute neighbourhood atrophy sc ctx
neighbourhood_atrophy_sc_ctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_ctx)))
colnames(neighbourhood_atrophy_sc_ctx) <- paste0(colnames(atrophy_ctx), "_neighbours")
for (i in 1:length(atrophy_ctx)){
  neighbours <- which(sc_ctx[i,]!=0)
  weight_vec <- sc_ctx[i,c(neighbours)]
  neighbourhood_atrophy_sc_ctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy sc sctx
neighbourhood_atrophy_sc_sctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_sctx)))
colnames(neighbourhood_atrophy_sc_sctx) <- paste0(colnames(atrophy_sctx), "_neighbours")
for (i in 1:length(atrophy_sctx)){
  neighbours <- which(sc_sctx[i,]!=0)
  weight_vec <- sc_sctx[i,c(neighbours)]
  neighbourhood_atrophy_sc_sctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy fc ctx
neighbourhood_atrophy_fc_ctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_ctx)))
colnames(neighbourhood_atrophy_fc_ctx) <- paste0(colnames(atrophy_ctx), "_neighbours")
for (i in 1:length(atrophy_ctx)){
  neighbours <- which(fc_ctx[i,]!=0)
  weight_vec <- fc_ctx[i,c(neighbours)]
  neighbourhood_atrophy_fc_ctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy fc sctx
neighbourhood_atrophy_fc_sctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_sctx)))
colnames(neighbourhood_atrophy_fc_sctx) <- paste0(colnames(atrophy_sctx), "_neighbours")
for (i in 1:length(atrophy_sctx)){
  neighbours <- which(fc_sctx[i,]!=0)
  weight_vec <- fc_sctx[i,c(neighbours)]
  neighbourhood_atrophy_fc_sctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy cge ctx
neighbourhood_atrophy_cge_ctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_ctx)))
colnames(neighbourhood_atrophy_cge_ctx) <- paste0(colnames(atrophy_ctx), "_neighbours")
for (i in 1:length(atrophy_ctx)){
  neighbours <- which(cge_ctx[i,]!=0 & !is.na(cge_ctx[i,]))
  weight_vec <- cge_ctx[i,c(neighbours)]
  neighbourhood_atrophy_cge_ctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy cge sctx
neighbourhood_atrophy_cge_sctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_sctx)))
colnames(neighbourhood_atrophy_cge_sctx) <- paste0(colnames(atrophy_sctx), "_neighbours")
for (i in 1:length(atrophy_sctx)){
  neighbours <- which(cge_sctx[i,]!=0 & !is.na(cge_sctx[i,]))
  weight_vec <- cge_sctx[i,c(neighbours)]
  neighbourhood_atrophy_cge_sctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# write out
write.csv(neighbourhood_atrophy_sc_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_long_sc_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_sc_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_long_sc_sctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_fc_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_long_fc_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_fc_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_long_fc_sctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_cge_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_long_cge_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_cge_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_long_cge_sctx.csv"), row.names = FALSE)
## individual-level
# load atrophy map
atrophy_mapping <- read.csv(here::here("data/atrophy_cross/atrophy_baseline_SchaeferSubcortex114Parcels_individuals_wscore.csv"), header = TRUE, sep = ",")
atrophy_ctx <- atrophy_mapping[,2:101]
atrophy_sctx <- atrophy_mapping[,102:115]
# compute neighbourhood atrophy sc ctx
neighbourhood_atrophy_sc_ctx <- data.frame(matrix(nrow = length(atrophy_mapping$id), ncol = (1+length(atrophy_ctx))))
colnames(neighbourhood_atrophy_sc_ctx) <- append("id", paste0(colnames(atrophy_ctx), "_neighbours"))
neighbourhood_atrophy_sc_ctx$id <- atrophy_mapping$id
for (l in atrophy_mapping$id){
  for (i in 1:length(atrophy_ctx)){
    neighbours <- which(sc_ctx[i,]!=0)
    weight_vec <- sc_ctx[i,c(neighbours)]
    neighbourhood_atrophy_sc_ctx[neighbourhood_atrophy_sc_ctx$id==l,i+1] <- sum(atrophy_ctx[atrophy_mapping$id==l,neighbours]*weight_vec)*(1/length(neighbours))
  }
  }
# compute neighbourhood atrophy sc sctx
neighbourhood_atrophy_sc_sctx <- data.frame(matrix(nrow = length(atrophy_mapping$id), ncol = (1+length(atrophy_sctx))))
colnames(neighbourhood_atrophy_sc_sctx) <- append("id", paste0(colnames(atrophy_sctx), "_neighbours"))
neighbourhood_atrophy_sc_sctx$id <- atrophy_mapping$id
for (l in atrophy_mapping$id){
  for (i in 1:length(atrophy_sctx)){
    neighbours <- which(sc_sctx[i,]!=0)
    weight_vec <- sc_sctx[i,c(neighbours)]
    neighbourhood_atrophy_sc_sctx[neighbourhood_atrophy_sc_sctx$id==l,i+1] <- sum(atrophy_ctx[atrophy_mapping$id==l,neighbours]*weight_vec)*(1/length(neighbours))
  }
  }
# compute neighbourhood atrophy fc ctx
neighbourhood_atrophy_fc_ctx <- data.frame(matrix(nrow = length(atrophy_mapping$id), ncol = (1+length(atrophy_ctx))))
colnames(neighbourhood_atrophy_fc_ctx) <- append("id", paste0(colnames(atrophy_ctx), "_neighbours"))
neighbourhood_atrophy_fc_ctx$id <- atrophy_mapping$id
for (l in atrophy_mapping$id){
  for (i in 1:length(atrophy_ctx)){
    neighbours <- which(fc_ctx[i,]!=0)
    weight_vec <- fc_ctx[i,c(neighbours)]
    neighbourhood_atrophy_fc_ctx[neighbourhood_atrophy_fc_ctx$id==l,i+1] <- sum(atrophy_ctx[atrophy_mapping$id==l,neighbours]*weight_vec)*(1/length(neighbours))
  }
  }
# compute neighbourhood atrophy fc sctx
neighbourhood_atrophy_fc_sctx <- data.frame(matrix(nrow = length(atrophy_mapping$id), ncol = (1+length(atrophy_sctx))))
colnames(neighbourhood_atrophy_fc_sctx) <- append("id", paste0(colnames(atrophy_sctx), "_neighbours"))
neighbourhood_atrophy_fc_sctx$id <- atrophy_mapping$id
for (l in atrophy_mapping$id){
  for (i in 1:length(atrophy_sctx)){
    neighbours <- which(fc_sctx[i,]!=0)
    weight_vec <- fc_sctx[i,c(neighbours)]
    neighbourhood_atrophy_fc_sctx[neighbourhood_atrophy_fc_sctx$id==l,i+1] <- sum(atrophy_ctx[atrophy_mapping$id==l,neighbours]*weight_vec)*(1/length(neighbours))
  }
  }
# compute neighbourhood atrophy cge ctx
neighbourhood_atrophy_cge_ctx <- data.frame(matrix(nrow = length(atrophy_mapping$id), ncol = (1+length(atrophy_ctx))))
colnames(neighbourhood_atrophy_cge_ctx) <- append("id", paste0(colnames(atrophy_ctx), "_neighbours"))
neighbourhood_atrophy_cge_ctx$id <- atrophy_mapping$id
for (l in atrophy_mapping$id){
  for (i in 1:length(atrophy_ctx)){
    neighbours <- which(cge_ctx[i,]!=0 & !is.na(cge_ctx[i,]))
    weight_vec <- cge_ctx[i,c(neighbours)]
    neighbourhood_atrophy_cge_ctx[neighbourhood_atrophy_cge_ctx$id==l,i+1] <- sum(atrophy_ctx[atrophy_mapping$id==l,neighbours]*weight_vec)*(1/length(neighbours))
  }
  }
# compute neighbourhood atrophy cge sctx
neighbourhood_atrophy_cge_sctx <- data.frame(matrix(nrow = length(atrophy_mapping$id), ncol = (1+length(atrophy_sctx))))
colnames(neighbourhood_atrophy_cge_sctx) <- append("id", paste0(colnames(atrophy_sctx), "_neighbours"))
neighbourhood_atrophy_cge_sctx$id <- atrophy_mapping$id
for (l in atrophy_mapping$id){
  for (i in 1:length(atrophy_sctx)){
    neighbours <- which(cge_sctx[i,]!=0 & !is.na(cge_sctx[i,]))
    weight_vec <- cge_sctx[i,c(neighbours)]
    neighbourhood_atrophy_cge_sctx[neighbourhood_atrophy_cge_sctx$id==l,i+1] <- sum(atrophy_ctx[atrophy_mapping$id==l,neighbours]*weight_vec)*(1/length(neighbours))
  }
  }
# write out
write.csv(neighbourhood_atrophy_sc_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_individuals_sc_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_sc_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_individuals_sc_sctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_fc_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_individuals_fc_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_fc_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_individuals_fc_sctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_cge_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_individuals_cge_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_cge_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_individuals_cge_sctx.csv"), row.names = FALSE)

### SchaeferSubcortex 414 parcels
## load normative functional/structural/CGE matrices
sc_ctx <- as.matrix(read.csv(here::here("data/hcp_connectivity/strucMatrix_ctx_schaefer_400.csv"), header = FALSE, sep = ","))
sc_sctx <- as.matrix(read.csv(here::here("data/hcp_connectivity/strucMatrix_sctx_schaefer_400.csv"), header = FALSE, sep = ","))
sc_sctx <- sc_sctx[,1:400] # the sctx matrix contains also subcortico-subcortical connections
fc_ctx <- as.matrix(read.csv(here::here("data/hcp_connectivity/funcMatrix_ctx_schaefer_400.csv"), header = FALSE, sep = ","))
fc_sctx <- as.matrix(read.csv(here::here("data/hcp_connectivity/funcMatrix_sctx_schaefer_400.csv"), header = FALSE, sep = ","))
cge_ctx <- as.matrix(read.csv(here::here("data/ge/CGE_SchaeferSubcortex414Parcels_ctx.csv"), header = FALSE, sep = ","))
cge_sctx <- as.matrix(read.csv(here::here("data/ge/CGE_SchaeferSubcortex414Parcels_sctx.csv"), header = FALSE, sep = ","))
## group-level
# load atrophy map
atrophy_mapping <- read.csv(here::here("data/atrophy_cross/atrophy_baseline_SchaeferSubcortex414Parcels_group_cohend.csv"), header = TRUE, sep = ",")
atrophy_ctx <- atrophy_mapping[,1:400]
atrophy_sctx <- atrophy_mapping[,401:414]
# compute neighbourhood atrophy sc ctx
neighbourhood_atrophy_sc_ctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_ctx)))
colnames(neighbourhood_atrophy_sc_ctx) <- paste0(colnames(atrophy_ctx), "_neighbours")
for (i in 1:length(atrophy_ctx)){
  neighbours <- which(sc_ctx[i,]!=0)
  weight_vec <- sc_ctx[i,c(neighbours)]
  neighbourhood_atrophy_sc_ctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy sc sctx
neighbourhood_atrophy_sc_sctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_sctx)))
colnames(neighbourhood_atrophy_sc_sctx) <- paste0(colnames(atrophy_sctx), "_neighbours")
for (i in 1:length(atrophy_sctx)){
  neighbours <- which(sc_sctx[i,]!=0)
  weight_vec <- sc_sctx[i,c(neighbours)]
  neighbourhood_atrophy_sc_sctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy fc ctx
neighbourhood_atrophy_fc_ctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_ctx)))
colnames(neighbourhood_atrophy_fc_ctx) <- paste0(colnames(atrophy_ctx), "_neighbours")
for (i in 1:length(atrophy_ctx)){
  neighbours <- which(fc_ctx[i,]!=0)
  weight_vec <- fc_ctx[i,c(neighbours)]
  neighbourhood_atrophy_fc_ctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy fc sctx
neighbourhood_atrophy_fc_sctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_sctx)))
colnames(neighbourhood_atrophy_fc_sctx) <- paste0(colnames(atrophy_sctx), "_neighbours")
for (i in 1:length(atrophy_sctx)){
  neighbours <- which(fc_sctx[i,]!=0)
  weight_vec <- fc_sctx[i,c(neighbours)]
  neighbourhood_atrophy_fc_sctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy cge ctx
neighbourhood_atrophy_cge_ctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_ctx)))
colnames(neighbourhood_atrophy_cge_ctx) <- paste0(colnames(atrophy_ctx), "_neighbours")
for (i in 1:length(atrophy_ctx)){
  neighbours <- which(cge_ctx[i,]!=0 & !is.na(cge_ctx[i,]))
  weight_vec <- cge_ctx[i,c(neighbours)]
  neighbourhood_atrophy_cge_ctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy cge sctx
neighbourhood_atrophy_cge_sctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_sctx)))
colnames(neighbourhood_atrophy_cge_sctx) <- paste0(colnames(atrophy_sctx), "_neighbours")
for (i in 1:length(atrophy_sctx)){
  neighbours <- which(cge_sctx[i,]!=0 & !is.na(cge_sctx[i,]))
  weight_vec <- cge_sctx[i,c(neighbours)]
  neighbourhood_atrophy_cge_sctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# write out
write.csv(neighbourhood_atrophy_sc_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_sc_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_sc_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_sc_sctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_fc_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_fc_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_fc_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_fc_sctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_cge_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_cge_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_cge_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_cge_sctx.csv"), row.names = FALSE)
## group-level long
# load atrophy map
atrophy_mapping <- read.csv(here::here("data/atrophy_cross/atrophy_baseline_SchaeferSubcortex414Parcels_group_long_cohend.csv"), header = TRUE, sep = ",")
atrophy_ctx <- atrophy_mapping[,1:400]
atrophy_sctx <- atrophy_mapping[,401:414]
# compute neighbourhood atrophy sc ctx
neighbourhood_atrophy_sc_ctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_ctx)))
colnames(neighbourhood_atrophy_sc_ctx) <- paste0(colnames(atrophy_ctx), "_neighbours")
for (i in 1:length(atrophy_ctx)){
  neighbours <- which(sc_ctx[i,]!=0)
  weight_vec <- sc_ctx[i,c(neighbours)]
  neighbourhood_atrophy_sc_ctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy sc sctx
neighbourhood_atrophy_sc_sctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_sctx)))
colnames(neighbourhood_atrophy_sc_sctx) <- paste0(colnames(atrophy_sctx), "_neighbours")
for (i in 1:length(atrophy_sctx)){
  neighbours <- which(sc_sctx[i,]!=0)
  weight_vec <- sc_sctx[i,c(neighbours)]
  neighbourhood_atrophy_sc_sctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy fc ctx
neighbourhood_atrophy_fc_ctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_ctx)))
colnames(neighbourhood_atrophy_fc_ctx) <- paste0(colnames(atrophy_ctx), "_neighbours")
for (i in 1:length(atrophy_ctx)){
  neighbours <- which(fc_ctx[i,]!=0)
  weight_vec <- fc_ctx[i,c(neighbours)]
  neighbourhood_atrophy_fc_ctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy fc sctx
neighbourhood_atrophy_fc_sctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_sctx)))
colnames(neighbourhood_atrophy_fc_sctx) <- paste0(colnames(atrophy_sctx), "_neighbours")
for (i in 1:length(atrophy_sctx)){
  neighbours <- which(fc_sctx[i,]!=0)
  weight_vec <- fc_sctx[i,c(neighbours)]
  neighbourhood_atrophy_fc_sctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy cge ctx
neighbourhood_atrophy_cge_ctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_ctx)))
colnames(neighbourhood_atrophy_cge_ctx) <- paste0(colnames(atrophy_ctx), "_neighbours")
for (i in 1:length(atrophy_ctx)){
  neighbours <- which(cge_ctx[i,]!=0 & !is.na(cge_ctx[i,]))
  weight_vec <- cge_ctx[i,c(neighbours)]
  neighbourhood_atrophy_cge_ctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# compute neighbourhood atrophy cge sctx
neighbourhood_atrophy_cge_sctx <- data.frame(matrix(nrow = 1, ncol = length(atrophy_sctx)))
colnames(neighbourhood_atrophy_cge_sctx) <- paste0(colnames(atrophy_sctx), "_neighbours")
for (i in 1:length(atrophy_sctx)){
  neighbours <- which(cge_sctx[i,]!=0 & !is.na(cge_sctx[i,]))
  weight_vec <- cge_sctx[i,c(neighbours)]
  neighbourhood_atrophy_cge_sctx[,i] <- sum(atrophy_ctx[neighbours]*weight_vec)*(1/length(neighbours))
  }
# write out
write.csv(neighbourhood_atrophy_sc_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_long_sc_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_sc_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_long_sc_sctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_fc_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_long_fc_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_fc_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_long_fc_sctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_cge_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_long_cge_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_cge_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_long_cge_sctx.csv"), row.names = FALSE)
## individual-level
# load atrophy map
atrophy_mapping <- read.csv(here::here("data/atrophy_cross/atrophy_baseline_SchaeferSubcortex414Parcels_individuals_wscore.csv"), header = TRUE, sep = ",")
atrophy_ctx <- atrophy_mapping[,2:401]
atrophy_sctx <- atrophy_mapping[,402:415]
# compute neighbourhood atrophy sc ctx
neighbourhood_atrophy_sc_ctx <- data.frame(matrix(nrow = length(atrophy_mapping$id), ncol = (1+length(atrophy_ctx))))
colnames(neighbourhood_atrophy_sc_ctx) <- append("id", paste0(colnames(atrophy_ctx), "_neighbours"))
neighbourhood_atrophy_sc_ctx$id <- atrophy_mapping$id
for (l in atrophy_mapping$id){
  for (i in 1:length(atrophy_ctx)){
    neighbours <- which(sc_ctx[i,]!=0)
    weight_vec <- sc_ctx[i,c(neighbours)]
    neighbourhood_atrophy_sc_ctx[neighbourhood_atrophy_sc_ctx$id==l,i+1] <- sum(atrophy_ctx[atrophy_mapping$id==l,neighbours]*weight_vec)*(1/length(neighbours))
  }
  }
# compute neighbourhood atrophy sc sctx
neighbourhood_atrophy_sc_sctx <- data.frame(matrix(nrow = length(atrophy_mapping$id), ncol = (1+length(atrophy_sctx))))
colnames(neighbourhood_atrophy_sc_sctx) <- append("id", paste0(colnames(atrophy_sctx), "_neighbours"))
neighbourhood_atrophy_sc_sctx$id <- atrophy_mapping$id
for (l in atrophy_mapping$id){
  for (i in 1:length(atrophy_sctx)){
    neighbours <- which(sc_sctx[i,]!=0)
    weight_vec <- sc_sctx[i,c(neighbours)]
    neighbourhood_atrophy_sc_sctx[neighbourhood_atrophy_sc_sctx$id==l,i+1] <- sum(atrophy_ctx[atrophy_mapping$id==l,neighbours]*weight_vec)*(1/length(neighbours))
  }
  }
# compute neighbourhood atrophy fc ctx
neighbourhood_atrophy_fc_ctx <- data.frame(matrix(nrow = length(atrophy_mapping$id), ncol = (1+length(atrophy_ctx))))
colnames(neighbourhood_atrophy_fc_ctx) <- append("id", paste0(colnames(atrophy_ctx), "_neighbours"))
neighbourhood_atrophy_fc_ctx$id <- atrophy_mapping$id
for (l in atrophy_mapping$id){
  for (i in 1:length(atrophy_ctx)){
    neighbours <- which(fc_ctx[i,]!=0)
    weight_vec <- fc_ctx[i,c(neighbours)]
    neighbourhood_atrophy_fc_ctx[neighbourhood_atrophy_fc_ctx$id==l,i+1] <- sum(atrophy_ctx[atrophy_mapping$id==l,neighbours]*weight_vec)*(1/length(neighbours))
  }
  }
# compute neighbourhood atrophy fc sctx
neighbourhood_atrophy_fc_sctx <- data.frame(matrix(nrow = length(atrophy_mapping$id), ncol = (1+length(atrophy_sctx))))
colnames(neighbourhood_atrophy_fc_sctx) <- append("id", paste0(colnames(atrophy_sctx), "_neighbours"))
neighbourhood_atrophy_fc_sctx$id <- atrophy_mapping$id
for (l in atrophy_mapping$id){
  for (i in 1:length(atrophy_sctx)){
    neighbours <- which(fc_sctx[i,]!=0)
    weight_vec <- fc_sctx[i,c(neighbours)]
    neighbourhood_atrophy_fc_sctx[neighbourhood_atrophy_fc_sctx$id==l,i+1] <- sum(atrophy_ctx[atrophy_mapping$id==l,neighbours]*weight_vec)*(1/length(neighbours))
  }
  }
# compute neighbourhood atrophy cge ctx
neighbourhood_atrophy_cge_ctx <- data.frame(matrix(nrow = length(atrophy_mapping$id), ncol = (1+length(atrophy_ctx))))
colnames(neighbourhood_atrophy_cge_ctx) <- append("id", paste0(colnames(atrophy_ctx), "_neighbours"))
neighbourhood_atrophy_cge_ctx$id <- atrophy_mapping$id
for (l in atrophy_mapping$id){
  for (i in 1:length(atrophy_ctx)){
    neighbours <- which(cge_ctx[i,]!=0 & !is.na(cge_ctx[i,]))
    weight_vec <- cge_ctx[i,c(neighbours)]
    neighbourhood_atrophy_cge_ctx[neighbourhood_atrophy_cge_ctx$id==l,i+1] <- sum(atrophy_ctx[atrophy_mapping$id==l,neighbours]*weight_vec)*(1/length(neighbours))
  }
  }
# compute neighbourhood atrophy cge sctx
neighbourhood_atrophy_cge_sctx <- data.frame(matrix(nrow = length(atrophy_mapping$id), ncol = (1+length(atrophy_sctx))))
colnames(neighbourhood_atrophy_cge_sctx) <- append("id", paste0(colnames(atrophy_sctx), "_neighbours"))
neighbourhood_atrophy_cge_sctx$id <- atrophy_mapping$id
for (l in atrophy_mapping$id){
  for (i in 1:length(atrophy_sctx)){
    neighbours <- which(cge_sctx[i,]!=0 & !is.na(cge_sctx[i,]))
    weight_vec <- cge_sctx[i,c(neighbours)]
    neighbourhood_atrophy_cge_sctx[neighbourhood_atrophy_cge_sctx$id==l,i+1] <- sum(atrophy_ctx[atrophy_mapping$id==l,neighbours]*weight_vec)*(1/length(neighbours))
  }
  }
# write out
write.csv(neighbourhood_atrophy_sc_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_individuals_sc_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_sc_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_individuals_sc_sctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_fc_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_individuals_fc_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_fc_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_individuals_fc_sctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_cge_ctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_individuals_cge_ctx.csv"), row.names = FALSE)
write.csv(neighbourhood_atrophy_cge_sctx, here::here("data/neighbourhood_atrophy/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_individuals_cge_sctx.csv"), row.names = FALSE)

