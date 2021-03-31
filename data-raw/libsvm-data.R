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

# getLibsvmData <- function(url) {
#   is_bz2 <- grepl(".bz2$", url)

#   tmp_file <- tempfile(fileext = if (is_bz2) ".bz2" else "")

#   download.file(url, tmp_file, quiet = TRUE)

#   readLIBSVM(tmp_file)
# }

datafiles <- c(
  "abalone" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/regression/abalone",
  "cadata" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/regression/cadata",
  "colon-cancer" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/colon-cancer.bz2",
  "covtype" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/covtype.libsvm.binary.bz2",
  "diabetes" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/diabetes",
  "duke-breast-cancer" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/duke.bz2",
  "e2006-log1p-test" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/regression/log1p.E2006.test.bz2",
  "e2006-log1p-train" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/regression/log1p.E2006.train.bz2",
  "e2006-tfidf-test" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/regression/E2006.test.bz2",
  "e2006-tfidf-train" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/regression/E2006.train.bz2",
  # "epsilon-train" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/epsilon_normalized.bz2",
  # "epsilon-test" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/epsilon_normalized.t.bz2",
  "german-numer" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/german.numer",
  "gisette-test" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/gisette_scale.t.bz2",
  "gisette-train" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/gisette_scale.bz2",
  "heart" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/heart",
  "ijcnn1-train" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/ijcnn1.tr.bz2",
  "ijcnn1-test" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/ijcnn1.t.bz2",
  "ionosphere" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/ionosphere_scale",
  "leukemia-test" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/leu.t.bz2",
  "leukemia-train" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/leu.bz2",
  "madelon-test" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/madelon.t",
  "madelon-train" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/madelon",
  "mushrooms" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/mushrooms",
  "news20" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/news20.binary.bz2",
  "phishing" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/phishing",
  "pyrim-scaled-expanded5" = "https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/QJEUKR/GVA3LP",
  "rcv1-train" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/rcv1_train.binary.bz2",
  "real-sim" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/real-sim.bz2",
  "sonar" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/sonar_scale",
  "splice-train" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/splice",
  "splice-test" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/splice.t",
  "YearPredictionMSD-test" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/regression/YearPredictionMSD.t.bz2",
  "YearPredictionMSD-train" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/regression/YearPredictionMSD.bz2"
  # "epsilon-test" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/epsilon_normalized.t.bz2",
  # "epsilon-train" = "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/epsilon_normalized.bz2",
  # NOTE(jolars): epsilon datasets are very large
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
