g6 <- strsplit(Gname,split = ";")
gdomain <- unlist(lapply(g6, function(x){x[1]}))
gkingdom <- unlist(lapply(g6, function(x){x[2]}))
gphylum <- unlist(lapply(g6, function(x){x[3]}))
gclass <- unlist(lapply(g6, function(x){x[4]}))
gorder <- unlist(lapply(g6, function(x){x[5]}))
gfamily <- unlist(lapply(g6, function(x){x[6]}))

dataorder <- read.csv("/Users/wangzhuofan/Project/coalescence/realdata/Order.csv",header = FALSE)
ds <- as.character(dataorder[1,-1])
rownames(dataorder) <- c("",dataorder$V1[-1])
colnames(dataorder) <- c("",ds)
dataorder <- dataorder[-1,-1]
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
order_oral <- dataorder[,pairs_record$oral]
order_stool <- dataorder[,pairs_record$stool]
ratio0 <- data.frame(oral=apply(order_oral,1,function(x){sum(x==0)/ncol(order_oral)}),stool=apply(order_stool,1,function(x){sum(x==0)/ncol(order_stool)}))
library(dplyr)
select_ratio <- which(ratio0$oral<0.9&ratio0$stool<0.9)
order_oral <- order_oral[select_ratio,]
order_stool <- order_stool[select_ratio,]
oral_mat <- matrix(as.numeric(as.matrix(order_oral)),nrow = 14)
stool_mat <- matrix(as.numeric(as.matrix(order_stool)),nrow = 14)
imp_oral <- as.matrix(impRZilr(t(oral_mat))$x)
imp_stool <- as.matrix(impRZilr(t(stool_mat))$x)
oral_scale <- t((imp_oral)/rowSums(imp_oral))
stool_scale <- t((imp_stool)/rowSums(imp_stool))
heatmap((oral_scale),Rowv = NA, Colv = NA)
heatmap((stool_scale),Rowv = NA, Colv = NA)
dd_data <- list(n = ncol(oral_scale), px = nrow(oral_scale), py = nrow(stool_scale), x = stool_scale, y = oral_scale)
dd_data <- list(n = 515, px = 14, py = 14, x = oral_scale, y = stool_scale)
model_dd = stan_model('coalescencezhuofan.stan')
fit = sampling(model_dd,dd_data,iter=2000,chains=4,cores = parallel::detectCores())
# heatmap(x,Rowv = NA, Colv = NA)
# heatmap(y,Rowv = NA, Colv = NA)
params = extract(fit)
sam_mest <- apply(params$mt,c(2,3),mean)
heatmap(t(sam_mest), Rowv = NA, Colv = NA)
