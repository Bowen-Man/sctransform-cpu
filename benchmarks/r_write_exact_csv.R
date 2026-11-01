# CSV is convenient to inspect, but write.csv's numeric formatting need not
# round-trip all binary64 values even with options(digits=17). Format each
# double explicitly before writing; R is used only for oracle generation.
write_exact_csv <- function(data, path, row.names=FALSE) {
  frame <- as.data.frame(data, check.names=FALSE)
  for (i in seq_along(frame)) {
    column <- frame[[i]]
    if (is.double(column)) {
      formatted <- sprintf('%.17g', column)
      formatted[is.na(column) & !is.nan(column)] <- NA_character_
      frame[[i]] <- formatted
    }
  }
  write.csv(frame, path, row.names=row.names, na='NA')
}
