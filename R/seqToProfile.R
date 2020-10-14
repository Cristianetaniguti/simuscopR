#' Create SimusCop profile
#' 
#' @param bam.file charater defining the BAM file
#' @param bed.file charater defining the BED file
#' @param vcf.file charater defining the VCF file
#' @param reference charater defining the reference genome FASTA file
#' @param out.profile charater defining the profile output name  
#' 
#' @export
seqToProfile <- function(bam.file = NULL, bed.file = NULL, vcf.file = NULL, 
                         reference = NULL, out.profile = NULL){
  system(
    paste0(system.file(package="simuscopR"),"/SimusCop/seqToProfile",
           " -b ", bam.file, " -t ", bed.file, " -v ", vcf.file, " -r ", reference, 
           " > ", out.profile)
  )
}