g6 <- strsplit(Gname,split = ";")
gdomain <- unlist(lapply(g6, function(x){x[1]}))
gkingdom <- unlist(lapply(g6, function(x){x[2]}))
gphylum <- unlist(lapply(g6, function(x){x[3]}))
gclass <- unlist(lapply(g6, function(x){x[4]}))
gorder <- unlist(lapply(g6, function(x){x[5]}))
gfamily <- unlist(lapply(g6, function(x){x[6]}))

datafamily <- read.csv("/Users/wangzhuofan/Project/coalescence/realdata/Family.csv",header = FALSE)
ds <- as.character(datafamily[1,-1])
rownames(datafamily) <- c("",datafamily$V1[-1])
colnames(datafamily) <- c("",ds)
datafamily <- datafamily[-1,-1]
metadata <- read.csv("/Users/wangzhuofan/Project/coalescence/realdata/metadata.csv")
ns <- metadata$Sample
ns_oral <- ns[metadata$Material=="Oral"]
ns_stool <- ns[metadata$Material=="Stool"]
meta_oral <- subset(metadata,Material=="Oral")
meta_stool <- subset(metadata,Material=="Stool")
pairsm <- data.frame(oral=ns_oral,stool=NA)
templ <- list()
for (i in 1:768) {
  pairsm$oral[i] <- ns_oral[i]
  temp <- which((meta_stool$Patient==metadata$Patient[i])&(meta_stool$Timepoint==metadata$Timepoint[i]))
  templ[[i]] <- temp
  if(length(temp)==1){
    pairsm$stool[i] <- ns_stool[temp]
  }
  if(length(temp)==2){
    pairsm$stool[i] <- ns_stool[min(temp)]
  }
}
pairsm <- na.omit(pairsm)
pairs_record <- pairsm
pairs_record$oral <- match(pairsm$oral,ds)
pairs_record$stool <- match(pairsm$stool,ds)
pairs_record <- na.omit(pairs_record)
#oral_index <- na.omit(match(pairsm$oral,ds))
#stool_index <- na.omit(match(pairsm$stool,ds))
family_oral <- datafamily[,pairs_record$oral]
family_stool <- datafamily[,pairs_record$stool]
ratio0 <- data.frame(oral=apply(family_oral,1,function(x){sum(x==0)/ncol(family_oral)}),stool=apply(family_stool,1,function(x){sum(x==0)/ncol(family_stool)}))
library(dplyr)
select_ratio_oral <- which(ratio0$oral<0.5)
select_ratio_stool <- which(ratio0$stool<0.5)
family_oral <- family_oral[select_ratio_oral,]
family_stool <- family_stool[select_ratio_stool,]
oral_mat <- matrix(as.numeric(as.matrix(family_oral)),nrow = length(select_ratio_oral))
stool_mat <- matrix(as.numeric(as.matrix(family_stool)),nrow = length(select_ratio_stool))
del <- which(colSums(stool_mat)==0)
oral_mat <- oral_mat[,-del]
stool_mat <- stool_mat[,-del]
library(robCompositions)
imp_oral <- as.matrix(impRZilr(t(oral_mat))$x)
imp_stool <- as.matrix(impRZilr(t(stool_mat))$x)
oral_scale <- t((imp_oral)/rowSums(imp_oral))
stool_scale <- t((imp_stool)/rowSums(imp_stool))
heatmap(oral_scale,Rowv = NA, Colv = NA)
heatmap(stool_scale,Rowv = NA, Colv = NA)
dd_data <- list(n = ncol(stool_scale), px = length(select_ratio_stool), py = length(select_ratio_oral), x = stool_scale, y = oral_scale)
dd_data <- list(n = 515, px = length(select_ratio), py = length(select_ratio), x = oral_scale, y = stool_scale)
library(rstan)
model_dd = stan_model('coalescencezhuofan.stan')
fit = sampling(model_dd,dd_data,iter=5000,chains=2,cores = parallel::detectCores())
 # heatmap(x,Rowv = NA, Colv = NA)
# h [eatmap(y,Rowv = NA, Colv = NA)
params = extract(fit)
sam_mest <- apply(params$mt,c(2,3),mean)
heatmap(t(sam_mest), Rowv = NA, Colv = NA,cexCol = 1,cexRow = 0.6)

oral_name = rownames(ratio0)[select_ratio_oral]
oraln = unlist(lapply(strsplit(oral_name,split = ";"), function(x){x[5]}))
stool_name = rownames(ratio0)[select_ratio_stool]
stooln = unlist(lapply(strsplit(stool_name,split = ";"), function(x){x[5]}))

rownames(sam_mest) <- oraln
colnames(sam_mest) <- stooln
