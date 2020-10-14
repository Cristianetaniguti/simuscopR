#' Converts variants in VCF to simuscop input file.
#'
#' Exclusively for biallelic SNVs, indels and insertions. Only one sample at the time and one chromosome.
#'
#' @param vcfR.object object of class vcfR
#' @param sample character defining sample ID to be evaluated
#' @param chrom character defining the chromosome to be evaluated
#'
#' @return list with data.frames for SNVs, indels and insertions
#'
#' @import dplyr
#'
#' @export
vcf2variants <- function(vcfR.object=NULL,
                         sample=NULL,
                         chrom=NULL){

  if(!is(vcfR.object, "vcfR")) stop("Input object is not of class vcfR")

  ALT <- vcfR.object@fix[,5]

  # Remove markers with more than one alternative
  rm_mks <- grep(ALT, pattern = ",")
  GT <- vcfR.object@gt

  # Genotypes
  GT <- GT %>% as.data.frame %>% select(all_of(sample))  %>%
    apply(.,2,as.character) %>% strsplit(., ":") %>% sapply(., "[[",1) %>%
    matrix(., ncol=length(sample)) %>%
    gsub(pattern = "[|]", replacement = "/")

  # remove missing data
  rm_mks <- unique(c(rm_mks, grep(pattern = "[.]", GT)))
  # remove homozygous for reference allele
  rm_mks <- unique(c(rm_mks, which(GT=="0/0")))

  ALT <- ALT[-rm_mks]
  REF <- vcfR.object@fix[,4][-rm_mks]
  CHROM <- vcfR.object@fix[,1][-rm_mks]
  POS <- vcfR.object@fix[,2][-rm_mks]
  GT <- vcfR.object@gt[-rm_mks,]

  if(!is.null(chrom)){
    chr <- which(CHROM %in% chrom)
    ALT <- ALT[chr]
    REF <- REF[chr]
    CHROM <- CHROM[chr]
    POS <- POS[chr]
    GT <- GT[chr,]
  }

  ALT_len <- str_length(ALT)
  REF_len <- str_length(REF)

  # Deletions
  del <- which(REF_len > ALT_len)
  del_sp <- REF_len[del] - ALT_len[del]

  # Insertions
  ins <- which(REF_len < ALT_len)
  ins_sp <- ALT_len[ins] - REF_len[ins]
  ins_sp <- str_sub(ALT[ins],-ins_sp,-1)

  GT <- GT %>% as.data.frame %>% select(all_of(sample))  %>%
    apply(.,2,as.character) %>% strsplit(., ":") %>% sapply(., "[[",1) %>%
    matrix(., ncol=length(sample)) %>%
    gsub(pattern = "[|]", replacement = "/")

  GT_simu <- sapply(strsplit(GT, "/"), function(x) if(x[1] == x[2]){
    "homo"
  } else {
    "het"
  } )

  GT_simu_matrix <- matrix(GT_simu, ncol=length(sample))

  # SNVs
  GT_SNVs <- as.vector(GT_simu_matrix[-c(ins, del),])

  variants_SNVs <- tibble(
    V1 = "s",
    V2 = rep(sample, each= length(REF[-c(ins, del)])),
    V3 = rep(CHROM[-c(ins, del)], length(sample)),
    V4 = rep(POS[-c(ins, del)], length(sample)),
    V5 = rep(REF[-c(ins, del)], length(sample)),
    V6 = rep(ALT[-c(ins, del)], length(sample)),
    V7 = GT_SNVs
  )

  # insertions
  GT_ins <- as.vector(GT_simu_matrix[ins,])

  variants_ins <- tibble(
    V1 = "i",
    V2 = rep(sample, each= length(ins)),
    V3 = rep(CHROM[ins], length(sample)),
    V4 = rep(POS[ins], length(sample)),
    V5 = rep(ins_sp, length(sample)),
    V6 = GT_ins
  )

  # deletions
  GT_del <- as.vector(GT_simu_matrix[del,])

  variants_del <- tibble(
    V1 = "d",
    V2 = rep(sample, each= length(del)),
    V3 = rep(CHROM[del], length(sample)),
    V4 = rep(POS[del], length(sample)),
    V5 = rep(del_sp, length(sample)),
    V6 = GT_del
  )

  variants <- list(SNVs = variants_SNVs, indels = variants_del,
                   insertions = variants_ins)

  return(variants)
}
