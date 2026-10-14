# Robust helper functions for GSE43292 / GSE163154 validation
# FIXED: column-name normalization now lowercases BEFORE removing non-alphanumeric
# characters, so "Gene" -> "gene" and "Module" -> "module".

clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "NULL", "-", "---")] <- NA_character_
  x
}

read_module_table_auto <- function(path) {

  attempts <- list(
    function() read.csv(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    function() read.delim(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    function() read.table(
      path,
      header = TRUE,
      sep = ";",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  )

  out <- NULL

  for (fun in attempts) {
    tmp <- tryCatch(fun(), error = function(e) NULL)

    if (!is.null(tmp) && ncol(tmp) >= 2L) {
      out <- tmp
      break
    }
  }

  if (is.null(out)) {
    stop(
      "Could not parse module assignment file as comma-, tab-, or semicolon-delimited text."
    )
  }

  out
}

read_brown_blue <- function(path) {

  x <- read_module_table_auto(path)

  cat("\nModule-assignment file detected:\n")
  cat("Path:", path, "\n")
  cat("Rows:", nrow(x), " Columns:", ncol(x), "\n")
  cat("Column names:\n")
  print(names(x))

  # IMPORTANT FIX:
  # tolower() must happen BEFORE gsub(), otherwise uppercase G/M in
  # "Gene"/"Module" are removed by the [a-z0-9] pattern.
  clean_names <- gsub(
    "[^a-z0-9]",
    "",
    tolower(names(x))
  )

  gene_idx <- which(
    clean_names %in%
      c(
        "gene",
        "genes",
        "symbol",
        "genesymbol",
        "geneid",
        "genename",
        "feature",
        "featureid"
      )
  )

  module_idx <- which(
    clean_names %in%
      c(
        "module",
        "modulecolor",
        "modulecolour",
        "color",
        "colour",
        "moduleassignment",
        "modulelabel"
      )
  )

  # Fallback: find the column that actually contains "brown" and "blue".
  if (length(module_idx) == 0L) {

    module_score <- vapply(
      x,
      function(v) {

        vals <- tolower(trimws(as.character(v)))

        sum(
          vals %in% c("brown", "blue"),
          na.rm = TRUE
        )
      },
      numeric(1)
    )

    if (max(module_score) > 0L) {
      module_idx <- which.max(module_score)
    }
  }

  if (length(module_idx) == 0L) {
    stop(
      "Could not identify the module column. ",
      "No column contains module labels such as 'brown' and 'blue'."
    )
  }

  module_idx <- module_idx[1]

  # Fallback for unusual gene headers.
  if (length(gene_idx) == 0L) {

    candidates <- setdiff(
      seq_len(ncol(x)),
      module_idx
    )

    gene_scores <- vapply(
      candidates,
      function(j) {

        vals <- trimws(as.character(x[[j]]))

        mean(
          grepl(
            "^[A-Za-z][A-Za-z0-9._-]{1,30}$",
            vals
          ),
          na.rm = TRUE
        )
      },
      numeric(1)
    )

    gene_idx <- candidates[
      which.max(gene_scores)
    ]
  }

  gene_idx <- gene_idx[1]

  cat(
    "Using gene column:",
    names(x)[gene_idx],
    "\n"
  )

  cat(
    "Using module column:",
    names(x)[module_idx],
    "\n"
  )

  gene <- clean_gene(
    x[[gene_idx]]
  )

  module <- tolower(
    trimws(
      as.character(
        x[[module_idx]]
      )
    )
  )

  brown <- unique(
    gene[
      module == "brown" &
        !is.na(gene)
    ]
  )

  blue <- unique(
    gene[
      module == "blue" &
        !is.na(gene)
    ]
  )

  cat(
    "Brown genes detected:",
    length(brown),
    "\n"
  )

  cat(
    "Blue genes detected:",
    length(blue),
    "\n\n"
  )

  if (length(brown) != 4262L) {
    warning(
      "Brown module size is ",
      length(brown),
      "; expected 4262 from the frozen Idea-2 WGCNA."
    )
  }

  if (length(blue) != 5197L) {
    warning(
      "Blue module size is ",
      length(blue),
      "; expected 5197 from the frozen Idea-2 WGCNA."
    )
  }

  list(
    brown = brown,
    blue = blue
  )
}

collapse_to_gene <- function(expr_probe, symbols) {

  symbols <- clean_gene(symbols)

  keep <- !is.na(symbols)

  limma::avereps(
    expr_probe[
      keep,
      ,
      drop = FALSE
    ],
    ID = symbols[
      keep
    ]
  )
}

row_zscore <- function(mat) {

  mat <- as.matrix(mat)

  z <- t(
    apply(
      mat,
      1,
      function(v) {

        s <- sd(
          v,
          na.rm = TRUE
        )

        if (is.na(s) || s == 0) {
          rep(
            0,
            length(v)
          )
        } else {
          (
            v -
              mean(
                v,
                na.rm = TRUE
              )
          ) /
            s
        }
      }
    )
  )

  rownames(z) <- rownames(mat)
  colnames(z) <- colnames(mat)

  z
}

score_gene_set <- function(expr_gene, genes) {

  genes <- intersect(
    unique(
      clean_gene(
        genes
      )
    ),
    rownames(
      expr_gene
    )
  )

  if (length(genes) < 2L) {

    return(
      list(
        score = rep(
          NA_real_,
          ncol(expr_gene)
        ),
        present = genes
      )
    )
  }

  sub <- expr_gene[
    genes,
    ,
    drop = FALSE
  ]

  vv <- apply(
    sub,
    1,
    var,
    na.rm = TRUE
  )

  use <- names(vv)[
    is.finite(vv) &
      vv > 0
  ]

  if (length(use) < 2L) {

    return(
      list(
        score = rep(
          NA_real_,
          ncol(expr_gene)
        ),
        present = genes
      )
    )
  }

  z <- row_zscore(
    sub[
      use,
      ,
      drop = FALSE
    ]
  )

  list(
    score = colMeans(
      z,
      na.rm = TRUE
    ),
    present = genes
  )
}

find_module_file <- function(preferred) {

  if (file.exists(preferred)) {
    return(preferred)
  }

  hits <- list.files(
    ".",
    pattern = "(module.*assign|wgcna.*module).*\\.csv$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )

  if (length(hits) == 1L) {
    return(hits[1])
  }

  NA_character_
}

map_probes_with_annotation_package <- function(
    probe_ids,
    annotation_package
) {

  if (!requireNamespace(annotation_package, quietly = TRUE)) {
    return(NULL)
  }

  db <- get(
    annotation_package,
    envir = asNamespace(annotation_package)
  )

  if (!"PROBEID" %in%
      AnnotationDbi::keytypes(db)) {
    return(NULL)
  }

  unname(
    AnnotationDbi::mapIds(
      db,
      keys = as.character(probe_ids),
      column = "SYMBOL",
      keytype = "PROBEID",
      multiVals = "first"
    )
  )
}

map_probes_with_gpl <- function(
    probe_ids,
    gpl_id
) {

  message(
    "Falling back to GEO platform annotation: ",
    gpl_id
  )

  tab <- GEOquery::Table(
    GEOquery::getGEO(
      gpl_id
    )
  )

  id_candidates <- c(
    "ID",
    "ID_REF",
    "Probe_Id",
    "PROBE_ID"
  )

  id_col <- id_candidates[
    id_candidates %in%
      names(tab)
  ][1]

  if (
    length(id_col) == 0L ||
    is.na(id_col)
  ) {
    stop(
      "Could not identify probe-ID column in ",
      gpl_id
    )
  }

  symbol_candidates <- grep(
    "gene.?symbol|symbol",
    names(tab),
    ignore.case = TRUE,
    value = TRUE
  )

  if (length(symbol_candidates) > 0L) {

    symbols <- as.character(
      tab[[symbol_candidates[1]]]
    )

  } else {

    assign_candidates <- grep(
      "gene.*assignment",
      names(tab),
      ignore.case = TRUE,
      value = TRUE
    )

    if (length(assign_candidates) == 0L) {
      stop(
        "Could not locate gene symbol information in ",
        gpl_id
      )
    }

    raw <- as.character(
      tab[[assign_candidates[1]]]
    )

    symbols <- vapply(
      raw,
      function(z) {

        if (
          is.na(z) ||
          !nzchar(z)
        ) {
          return(
            NA_character_
          )
        }

        bits <- trimws(
          unlist(
            strsplit(
              z,
              " /// | // ",
              perl = TRUE
            )
          )
        )

        hit <- bits[
          grepl(
            "^[A-Za-z0-9][A-Za-z0-9._-]{1,20}$",
            bits
          )
        ]

        if (length(hit) == 0L) {
          NA_character_
        } else {
          hit[1]
        }
      },
      character(1)
    )
  }

  symbols[
    match(
      as.character(probe_ids),
      as.character(
        tab[[id_col]]
      )
    )
  ]
}

inflammasome_65 <- c(
  "ABHD17A","ABHD8","ADCY7","ATAT1","BRCC3","BTK","CARD8","CASR",
  "CCDC92","CD36","CPTP","CSNK1A1","DDX3X","DHX33","EIF2AK2","ELP6",
  "FBXL2","GBP2","GBP5","GKN2","HSPA8","IFI16","IRAK1","IRAK2",
  "IRAK4","IRGM","KCNJ8","KCNK13","KCNK6","LAMP2","LATS1","LATS2",
  "MAPK8","MARCHF5","MARK4","MAVS","MEFV","MYD88","NEK7","NLRC3",
  "NLRP1","NLRP3","NLRP6","OGT","P2RX7","PLCG2","PPP2CA","PRKD1",
  "PTPN22","PYCARD","PYDC1","PYDC2","SIRT2","STMP1","TLR4","TLR6",
  "TLR9","TREM2","TRIM11","TRIM31","TRIM65","USP50","ZDHHC1",
  "ZDHHC12","ZDHHC5"
)
