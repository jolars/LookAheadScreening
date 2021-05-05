readLIBSVM <- function(file) {
  dat <- data.table::fread(file, sep = "\n", header = FALSE)[[1]]

  l <- stringi::stri_split_regex(dat, "[ ]+")

  y <- as.numeric(vapply(l, "[", 1, FUN.VALUE = character(1)))

  vals <- do.call(rbind, lapply(l, function(x) {
    do.call(rbind, stringi::stri_split_fixed(x[-1], ":"))
  }))

  row_ind <- rep(seq_len(length(l)), times = lengths(l) - 1)
  col_ind <- as.integer(vals[, 1])

  X <- Matrix::sparseMatrix(row_ind, col_ind, x = as.numeric(vals[, 2]))

  density <- Matrix::nnzero(X) / length(X)

  if (density > 0.5) {
    X <- as.matrix(X)
  }

  if (length(unique(y)) == 2) {
    y <- as.numeric(as.factor(y)) - 1
  }

  list(X = X, y = y)
}

getLibsvmData <- function(url) {
  is_bz2 <- grepl(".bz2$", url)

  tmp_file <- tempfile(fileext = if (is_bz2) ".bz2" else "")

  download.file(url, tmp_file, quiet = TRUE)

  readLIBSVM(tmp_file)
}

colon_cancer <- getLibsvmData(
  "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/colon-cancer.bz2"
)
cpusmall <- getLibsvmData(
  "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/regression/cpusmall"
)
heart <- getLibsvmData(
  "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/heart"
)
leukemia <- getLibsvmData(
  "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/leu.t.bz2"
)

usethis::use_data(colon_cancer, cpusmall, heart, leukemia, overwrite = TRUE)
