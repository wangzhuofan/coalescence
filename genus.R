g6 <- strsplit(Gname,split = ";")
gdomain <- unlist(lapply(g6, function(x){x[1]}))
gkingdom <- unlist(lapply(g6, function(x){x[2]}))
gphylum <- unlist(lapply(g6, function(x){x[3]}))
gclass <- unlist(lapply(g6, function(x){x[4]}))
gorder <- unlist(lapply(g6, function(x){x[5]}))
gfamily <- unlist(lapply(g6, function(x){x[6]}))

datagenus <- read.csv("/Users/wangzhuofan/Project/coalescence/realdata/Genus.csv",header = FALSE)
ds <- as.character(datagenus[1,-1])
rownames(datagenus) <- c("",datagenus$V1[-1])
colnames(datagenus) <- c("",ds)
datagenus <- datagenus[-1,-1]
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
genus_oral <- datagenus[,pairs_record$oral]
genus_stool <- datagenus[,pairs_record$stool]
ratio0 <- data.frame(oral=apply(genus_oral,1,function(x){sum(x==0)/ncol(genus_oral)}),stool=apply(genus_stool,1,function(x){sum(x==0)/ncol(genus_stool)}))
library(dplyr)
select_ratio <- which(ratio0$oral<0.9&ratio0$stool<0.9)
genus_oral <- genus_oral[select_ratio,]
genus_stool <- genus_stool[select_ratio,]
oral_mat <- matrix(as.numeric(as.matrix(genus_oral)),nrow = 15)
stool_mat <- matrix(as.numeric(as.matrix(genus_stool)),nrow = 15)
oral_mat <- oral_mat[,-del]
stool_mat <- stool_mat[,-del]
imp_oral <- as.matrix(impRZilr(t(oral_mat))$x)
imp_stool <- as.matrix(impRZilr(t(stool_mat))$x)
oral_scale <- t((imp_oral)/rowSums(imp_oral))
stool_scale <- t((imp_stool)/rowSums(imp_stool))
heatmap(t(oral_scale),Rowv = NA, Colv = NA)
heatmap(stool_scale,Rowv = NA, Colv = NA)
dd_data <- list(n = 514, px = 15, py = 15, x = stool_scale, y = oral_scale)
dd_data <- list(n = 514, px = 15, py = 15, x = oral_scale, y = stool_scale)
model_dd = stan_model('coalescencezhuofan.stan')
fit = sampling(model_dd,dd_data,iter=2000,chains=4,cores = parallel::detectCores())

# heatmap(x,Rowv = NA, Colv = NA)
# heatmap(y,Rowv = NA, Colv = NA)
params = extract(fit)
sam_mest <- apply(params$mt,c(2,3),mean)
heatmap(t(sam_mest), Rowv = NA, Colv = NA)

for (i in 1:15) {
  for (j in 1:15) {
    plot(x[i,],y[j,])
  }
}
testind <- c(1,3,5,6,7,8,10,11,12)
test_oral <- imp_oral[,testind]
test_stool <- imp_stool[,testind]
t_oral_scale <- t((test_oral)/rowSums(test_oral))
t_stool_scale <- t((test_stool)/rowSums(test_stool))
heatmap(t(t_oral_scale),Rowv = NA, Colv = NA)
heatmap(t(t_stool_scale),Rowv = NA, Colv = NA)
dd_data <- list(n = 514, px = 9, py = 9, x = t_stool_scale, y = t_oral_scale)

