map <- function(seqId = character(), chr, start, end, mapping_file) {
  map_temp <- mapping_file[mapping_file$target == seqId, ]
  if (nrow(map_temp) == 0) {
    print("seqId not in mapping file")
    return("trans")
  } else{
    map_temp <- map_temp[map_temp$chromosome == chr, ]
    if (nrow(map_temp) == 0) {
      return("trans")
    } else{
      ir1 <- IRanges(start = map_temp$cis_start, end = map_temp$cis_end)
      ir2 <- IRanges(start = start, end = end)
      ov <- countOverlaps(ir1, ir2)
      if (TRUE %in% ov >= 1) {
        return("cis")
      } else{
        return("trans")
      }
    }
  }
}
