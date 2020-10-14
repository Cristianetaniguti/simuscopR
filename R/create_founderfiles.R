#' Create PedigreeSim founder file
#'
#' @param vcfR.object object of class vcfR
#' @param parents_id vector of character defining parents ID
#' @param chrom character defining the chromosome to be evaluated
#'
#' @import dplyr
#' @export
create_founderfiles <- function(vcfR.object, parents_id, chrom){
  if(!is(vcfR.object, "vcfR")) stop("Input object is not of class vcfR")

  ALT <- vcfR.object@fix[,5]

  # Remove markers with more than one alternative
  rm_mks <- grep(ALT, pattern = ",")
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

  # Genotypes
  GT <- GT %>% as.data.frame %>% select(all_of(parents_id))  %>%
    apply(.,2,as.character) %>% strsplit(., ":") %>% sapply(., "[[",1) %>%
    matrix(., ncol=length(samples)) %>%
    gsub(pattern = "[|]", replacement = "/")

  colnames(GT) <- parents_id

  P1 <- GT %>% as.data.frame %>% select(.,parents_id[1])
  P1 <- strsplit(as.character(P1[,1]), "/")
  P1_1 <- sapply(P1, "[[", 1)
  P1_2 <- sapply(P1, "[[", 2)
  P1_1[P1_1 == "0"] <- REF[P1_1 == "0"]
  P1_1[P1_1 == "1"] <- ALT[P1_1 == "1"]
  P1_2[P1_2 == "0"] <- REF[P1_2 == "0"]
  P1_2[P1_2 == "1"] <- ALT[P1_2 == "1"]

  P2 <- GT %>% as.data.frame %>% select(.,parents_id[2])
  P2 <- strsplit(as.character(P2[,1]), "/")
  P2_1 <- sapply(P2, "[[", 1)
  P2_2 <- sapply(P2, "[[", 2)
  P2_1[P2_1 == "0"] <- REF[P2_1 == "0"]
  P2_1[P2_1 == "1"] <- ALT[P2_1 == "1"]
  P2_2[P2_2 == "0"] <- REF[P2_2 == "0"]
  P2_2[P2_2 == "1"] <- ALT[P2_2 == "1"]

  # Phases in parents will be the same as the VCF
  founderfile <- tibble(
    marker = str_c(CHROM, POS, sep = "_"),
    P1_1,
    P1_2,
    P2_1,
    P2_2
  )

  return(founderfile)
}
