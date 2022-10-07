#### load libraries
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(VariantAnnotation)
library(tidyr)
source("pks_context_selection_functions.R")

flank <- 0

##### list the vcfs you want to classify 
vcfs = list.files(path="path/to/vcfs", full.names = T)

names_vcfs = gsub(vcfs, ...) ### whatever pattern gets you the sample name 

#### function for loading vcfs
vcfs_list = function(vcf_file_list, names) {
  
  vcfs_grl = list()
  
  for (i in 1:length(vcf_file_list)) {
    vcf_name = vcf_file_list[i]
    print(vcf_name)
    name = names[i]
    print(name)
    vcf = readVcf(vcf_name)
    gr = granges(vcf)
    gr = gr[gr$FILTER == "PASS"]
    vcfs_grl[[name]] = gr
  }
  
  vcfs_grl = GRangesList(vcfs_grl)
  seqlevelsStyle(vcfs_grl) = "UCSC"
  
  return(vcfs_grl)
}

#### load vcfs in right format 
grl_samples = vcfs_list(vcfs, names_samples)


###### function to select snvs and 10 basepair context 
filter_pks = function(grl, output_folder, n_cores = 1){
  grl_snvs = mclapply(grl, FUN = select_snvs, mc.cores = n_cores)
  grls_pks_snv = mclapply(grl_snvs, FUN= select_context_snv, mc.cores = n_cores)


  for (name in names(grl)) {
    gr = grls_pks_snv[[name]]
    df <- data.frame(seqnames=gsub("chr","",seqnames(gr)),
                      starts=start(gr)-1-flank,
                      ends=end(gr)+flank,
                      names=gr$context,
                      scores=c(rep(".", length(gr))),
                      strands=strand(gr)
     )
     write.table(df, file= paste0(output_folder, "/snvs/", name, "_snvs.bed"), quote=F, sep="\t", row.names=F, col.names=F)
   }
}


#### provide grl and whatever output folder you want. Output folder should contain folder "snvs"
filter_pks(grl_samples, "samples")

