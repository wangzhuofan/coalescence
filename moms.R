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