# Load data, setwd() in data dir with all idat to create csv file of pD define in Samples Selections
## Effectifs par technologies
ChAMP_csv_v2 <- split(ChAMP_csv_v2, ChAMP_csv_v2$Sample_Plate)
PlateEPICv1 <- table(ChAMP_csv_v2$EPICv1$Sample_Group)
PlateEPICV1
PlateEPICv2 <- table(ChAMP_csv_v2$EPICv2$Sample_Group)
#write.table(ChAMP_csv_v2,"pD_MD_CSV_V2.csv", row.names=F, quote=F, sep=",")

myLoads_list <- setNames(lapply(names(ChAMP_csv_v2), function(arrayType){
  write.table(ChAMP_csv_v2[[arrayType]], paste0(getwd(),"/pD_ChAMP.csv"), row.names=F, quote=F, sep=",")
  myLoad <- champ.load(getwd(),method="ChAMP",
                       methValue="B",
                       autoimpute=TRUE,
                       filterDetP=TRUE,
                       ProbeCutoff=0,
                       SampleCutoff=0.1,
                       detPcut=0.01,
                       filterBeads=TRUE,
                       beadCutoff=0.05,
                       filterNoCG=TRUE,
                       filterSNPs=TRUE,
                       population=NULL,
                       filterMultiHit=TRUE,
                       filterXY=TRUE,
                       force=FALSE,
                       arraytype=arrayType)
  return(myLoad)
}),names(ChAMP_csv_v2))
# Keep probes in common between V1 and V2
probes_in_com <- intersect(rownames(myLoads_list$EPICv1$beta),rownames(myLoads_list$EPICv2$beta))
probes_V2<-unlist(strsplit(rownames(myLoads_list$EPICv2$beta), "_.*"))
probes_V2
probes_in_com <- intersect(rownames(myLoads_list$EPICv1$beta),probes_V2)
length(rownames(myLoads_list[["EPICv1"]]$beta))
length(rownames(myLoads_list[["EPICv2"]]$beta))
length(probes_in_com)
#change rownames of EPICv2
myLoads_list2<-myLoads_list
rownames(myLoads_list2[["EPICv2"]]$beta)<-probes_V2
rownames(myLoads_list2[["EPICv2"]]$intensity)<-probes_V2
# Create a new myLoad object
myLoad <- NULL
myLoad$beta <- cbind(myLoads_list[["EPICv1"]]$beta[probes_in_com,],myLoads_list2[["EPICv2"]]$beta[probes_in_com,])
myLoad$intensity <- cbind(myLoads_list[["EPICv1"]]$intensity[probes_in_com,],myLoads_list2[["EPICv2"]]$intensity[probes_in_com,])
myLoad$pd <- rbind(myLoads_list[["EPICv1"]]$pd,myLoads_list2[["EPICv2"]]$pd)
dim(myLoad$beta)
myLoad$pd
