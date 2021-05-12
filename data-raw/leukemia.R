readLIBSVM <- function(file) {
  dat <- data.table::fread(file,
    sep = "\n",
    header = FALSE,
    showProgress = FALSE
  )[[1]]

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

datafiles <- c(
  "leukemia-train" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/leu.bz2"
)

datanames <- names(datafiles)

for (i in seq_along(datanames)) {
  cat(i, "/", length(datanames), ": ", datanames[i], "..", sep = "")

  file <- file.path("data", paste0(datanames[i], ".rds"))

  if (file.exists(file)) {
    cat(" already downloaded, skipping!\n")
  } else {
    cat(" downloading and parsing .. ")
    d <- readLIBSVM(datafiles[i])
    cat(" done!\n")

    density <- Matrix::nnzero(d$X) / length(d$X)

    if (density > 0.25) {
      d$X <- as.matrix(d$X)
    }

    saveRDS(d, file)
  }
}
