library(data.table)

####PART 0: definition of the list of aptamers per target#####

###load list of targets for assay version
candia2<-fread("/exchange/healthds/pQTL/Reference_datasets_for_QC_proteomics/Candia/Candia_supp1.csv",skip=2)
candia2<-candia2[,c(1,10,19)]
colnames(candia2)<-c("k7","k5","k1")
candia2<-lapply(candia2, function(x) paste("seq.",gsub("-",".",x),sep=""))
list_k7<-candia2$k7
list_k5<-candia2$k5
list_k1<-candia2$k1
k4<-read.csv2("/exchange/healthds/pQTL/Reference_datasets_for_QC_proteomics/Candia/202406_table_aptamers_sun.csv")
list_k4<-gsub("-",".",k4$V1)
list_k4<-paste("seq.",list_k4,sep="")

####PART 1: annotation of a dataset containing a "target" variable####
##the function defines in this part directly take as input a dataset containing
##a "target" column containing the seqid ("seq.XXXX.X"), and the 4 list of aptamers
##corresponding to each of the assay versions.
##It returns the input dataset with 2 additional variables:
## "new_target" which gives the information about the versions of the assay
##containing the aptamer
## "new" which gives the information about the aptamer being new in 7k or not

#create function to map assay versions with dataset input
assay_annotation_on_dataset<-function(lit, list_k7,list_k5,list_k4,list_k1){
lit$new_target<-rep("new",nrow(lit))
lit$new_target<-ifelse(lit$study_id%in%list_k5&lit$study_id%in%list_k1,"already_in_5k_and_1k",lit$new_target)
lit$new_target<-ifelse(lit$study_id%in%list_k5&!lit$study_id%in%list_k1,"already_in_5k",lit$new_target)
lit$new_target<-ifelse(lit$study_id%in%list_k1&!lit$study_id%in%list_k5,"already_in_1k",lit$new_target)
lit$new_target<-ifelse(lit$new_target=="already_in_5k"&lit$study_id%in%list_k4,"already_in_5k_and_4k",lit$new_target)
lit$new_target<-ifelse(lit$new_target=="already_in_1k"&lit$study_id%in%list_k4,"already_in_4k_and_1k",lit$new_target)
lit$new_target<-ifelse(lit$new_target=="already_in_5k_and_1k"&lit$study_id%in%list_k4,"already_in_5k_4k_and_1k",lit$new_target)
lit$new_target<-ifelse(lit$new_target=="new"&lit$study_id%in%list_k4,"already_in_4k",lit$new_target)
lit$new<-ifelse(lit$new_target=="new",TRUE,FALSE)
return(lit)
}


##load file to annotate
lit<-fread("")
##apply function
lit<-assay_annotation_on_dataset(lit,list_k7,list_k5,list_k4,list_k1)

####PART 2: annotation function with aptamer as an input
##the function defines in this part directly take as input an aptamer (format:
##"seq.XXXX.X" and the 4 list of aptamers corresponding to each of the assay versions.
##it returns a list of characters with each of the assay version containing
##this aptamer.

annotation_aptamer<-function(target, list_k7,list_k5,list_k4,list_k1){
  a<-list()
  if (target%in%list_k7){a<-c(a,list("in_7k_assay"))}
  if (target%in%list_k5){a<-c(a,list("in_5k_assay"))}
  if (target%in%list_k4){a<-c(a,list("in_4k_assay"))}
  if (target%in%list_k1){a<-c(a,list("in_1k_assay"))}
  return(a)
}
