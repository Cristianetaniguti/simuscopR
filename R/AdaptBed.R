#' Adapt BED positions to RADseq simulation
#'
#' @param bed.file character defining BED file
#' @param out.file character defining output BED file name
#' @param int.size interger defining threshold of interval size to be kept in dataset
#'
#' @import intervals
#' @export
AdaptBed <- function(bed.file, out.file, int.size){
  df <- read.table(bed.file)
  int <- Intervals(df[,2:3])
  df_edit <- interval_intersection(int)

  diff <- df_edit@.Data[,2] - df_edit@.Data[,1]

  start <- df_edit@.Data[,1] - mean(diff) # Realocate by the mean of differences between start and end of the interval
  end <- df_edit@.Data[,2] - mean(diff)

  idx <- which(diff < int.size) # Keep only intervals smaller than threshold

  new <- data.frame(unique(df$V1), round(start[idx],0), round(end[idx],0))

  write.table(new, out.file, quote = F, col.names = F, row.names = F, sep = "\t")
}
