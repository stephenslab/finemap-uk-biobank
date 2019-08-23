library(susieR)
library(Matrix)

load('/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/height.chr3.141027000.141163045.RData')
info = readRDS('/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/height.chr3.141027000.141163045.info.rds')
maf = info$maf
rm(info)

time_attr = system.time(X.attr <- susieR:::set_X_attributes(X, center=T, scale=T))
cat(sprintf("Computing X attributes takes %d seconds.\n",round(time_attr["elapsed"])))
rm(X.attr)

time_mod = system.time(mod <- susie(X=X, Y=y, verbose=T))
cat(sprintf("Fitting SuSiE model takes %d seconds.\n",round(time_mod["elapsed"])))
rm(mod)

time_mod1 = system.time(mod1 <- susie(X=X, Y=y, verbose=T, L=1, max_iter=1, estimate_residual_variance=F, estimate_prior_variance=F))
cat(sprintf("Fitting SuSiE model (L=1, max_iter=1, estimate prior = F, estimate resid = F) takes %d seconds.\n",round(time_mod1["elapsed"])))
rm(mod1)

## Filter MAF
time_attr_maf = system.time(X.attr <- susieR:::set_X_attributes(X[, which(maf > 0.01)], center=T, scale=T))
cat(sprintf("MAF filter: Computing X attributes takes %d seconds.\n",round(time_attr_maf["elapsed"])))
rm(X.attr)

time_mod_maf = system.time(mod.maf <- susie(X=X[,which(maf > 0.01)], Y=y, verbose=T))
cat(sprintf("MAF filter: Fitting SuSiE model toobks %d seconds.\n",round(time_mod_maf["elapsed"])))
rm(mod.maf)

time_mod1_maf = system.time(mod1.maf <- susie(X=X[,which(maf>0.01)], Y=y, verbose=T, L=1, max_iter=1, estimate_residual_variance=F, estimate_prior_variance=F))
cat(sprintf("MAF filter: Fitting SuSiE model (L=1, max_iter=1, estimate prior = F, estimate resid = F) takes %d seconds.\n",round(time_mod1_maf["elapsed"])))
rm(mod1.maf)

cat('\n')
cat('Sparse:\n')
X.sparse = as(X, 'dgCMatrix')
time_attr = system.time(X.attr <- susieR:::set_X_attributes(X.sparse, center=T, scale=T))
cat(sprintf("Computing X attributes takes %d seconds.\n",round(time_attr["elapsed"])))
rm(X.attr)

time_mod = system.time(mod <- susie(X=X.sparse, Y=y, verbose=T))
cat(sprintf("Fitting SuSiE model takes %d seconds.\n",round(time_mod["elapsed"])))

time_mod1 = system.time(mod1 <- susie(X=X.sparse, Y=y, verbose=T, L=1, max_iter=1, estimate_residual_variance=F, estimate_prior_variance=F))
cat(sprintf("Fitting SuSiE model (L=1, max_iter=1, estimate prior = F, estimate resid = F) takes %d seconds.\n",round(time_mod1["elapsed"])))
rm(mod1)

## Filter MAF
time_attr_maf = system.time(X.attr <- susieR:::set_X_attributes(X.sparse[, which(maf > 0.01)], center=T, scale=T))
cat(sprintf("MAF filter: Computing X attributes takes %d seconds.\n",round(time_attr_maf["elapsed"])))
rm(X.attr)

time_mod_maf = system.time(mod.maf <- susie(X=X.sparse[,which(maf > 0.01)], Y=y, verbose=T))
cat(sprintf("MAF filter: Fitting SuSiE model toobks %d seconds.\n",round(time_mod_maf["elapsed"])))

time_mod1_maf = system.time(mod1.maf <- susie(X=X.sparse[,which(maf>0.01)], Y=y, verbose=T, L=1, max_iter=1, estimate_residual_variance=F, estimate_prior_variance=F))
cat(sprintf("MAF filter: Fitting SuSiE model (L=1, max_iter=1, estimate prior = F, estimate resid = F) takes %d seconds.\n",round(time_mod1_maf["elapsed"])))
rm(mod1.maf)

save(mod, mod.maf, file = "height.chr3.141027000.141163045.susie.model.RData")

