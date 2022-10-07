# ----- Functions to determine pks-patterns in whole-genome cancer driver data
library(GenomicRanges)
mutationDataToGranges = function(mutation_data) {
  # function to convert cBioPortal mutation tables to granges format
  mut_table = data.frame(chr = paste0("chr",mutation_data$chr), 
                         start = mutation_data$start_position, 
                         end = mutation_data$end_position,
                         REF = mutation_data$reference_allele, 
                         ALT = mutation_data$variant_allele)
  
  granges = makeGRangesFromDataFrame(mut_table, ignore.strand = T, keep.extra.columns = T)
  seqlevelsStyle(granges) = "UCSC"
  seqlevels(granges, pruning.mode ="coarse") = seqlevels(granges)[seqlevels(granges) != "chrNA"]
  
  return(granges)
}

select_snvs = function(mutation_gr) {
  # select only snvs in the indel
  mutation_gr = remove_mult_alts(mutation_gr)
  mutation_gr_snv = mutation_gr[nchar(mutation_gr$REF) == 1 &
                                  nchar(unlist(mutation_gr$ALT)) == 1 &
                                  !is.na(seqnames(mutation_gr)) & start(ranges(mutation_gr)) != -1 & seqnames(mutation_gr) != "chrNA" &
                                  mutation_gr$REF != "-" &
                                  unlist(mutation_gr$ALT) != "-", ]
  
  strand(mutation_gr_snv) = ifelse(mutation_gr_snv$REF == "C" | mutation_gr_snv$REF == "T", "+", "-")
  
  return(mutation_gr_snv)
}

#### find context around mutation
select_context_snv = function(gr) {
  gr = remove_mult_alts(gr)
  chr = seqnames(gr)
  start = start(gr) - 10 %>% as.integer()
  end = end(gr) + 10 %>% as.integer()
  strand = strand(gr) %>% as.character()
  contexts = getSeq(BSgenome.Hsapiens.UCSC.hg38, names = chr, start =start , end = end ,strand = strand, as.character = T )
  
  gr$context = contexts
  
  gr <- gr[gr$REF == 'T' | gr$REF == 'A']
  
  return(gr)
}

remove_mult_alts = function(gr) {
  
  mult_alts = elementNROWS(gr$ALT) > 1
  nr_mult_alts = sum(mult_alts)
  if (nr_mult_alts > 0){
    gr = gr[!mult_alts]
    warning(paste0("Removed ", nr_mult_alts, " variant(s), because they had multiple alternative alleles."), call. = F)
  }
  return(gr)
}
