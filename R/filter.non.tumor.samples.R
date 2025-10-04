filter.non.tumor.samples <- function(raw.datum, only.primary=TRUE) {
  # 01 is primary, 06 is metastatic, 03 is blood derived cancer
  if (!only.primary)
    return(raw.datum[,substring(colnames(raw.datum), 14, 15) %in% c('01', '03', '06')])
  else
    return(raw.datum[,substring(colnames(raw.datum), 14, 15) %in% c('01')])
}
normalize.matrix <- function(data.matrix) {
  temp = data.matrix - rowMeans(data.matrix)
  should.keep = (apply(temp, 1, sd) != 0)
  return ((temp / apply(temp, 1, sd))[should.keep, ])
}
