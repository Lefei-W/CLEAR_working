# Download credible sets and locus-to-gene predictions for GCST studies.

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(data.table)
  library(GenomicRanges)
  library(S4Vectors)
  library(IRanges)
})

PROJECT_ROOT <- "/working/lab_jonathb/lefeiW/projects/ATAC_BCAC"
STUDIES_PATH <- file.path(
  PROJECT_ROOT,
  "data/opentarget_trait/_all_studies/studies.tsv"
)

OUT_ROOT <- "/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/API_opentargets/data/GCST_download"
STUDY_DIR <- file.path(OUT_ROOT, "per_study")
LOG_DIR <- file.path(OUT_ROOT, "logs")

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(STUDY_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

STATUS_FILE <- file.path(OUT_ROOT, "GCST_processing_status.tsv")

BASE_URL <- "https://api.platform.opentargets.org/api/v4/graphql"
PAGE_IDS <- 500L
PAGE_LOCUS <- 500L
SLEEP_SEC <- 0.12

get_proxy_credentials <- function() {
  proxy <- c(Sys.getenv("https_proxy"), Sys.getenv("http_proxy"))
  proxy <- proxy[nzchar(proxy)][1]
  x <- sub("^https?://", "", proxy)
  
  cred <- sub("@.*$", "", x)
  host_port <- sub("^.*@", "", x)
  host_port <- sub("/+$", "", host_port)
  
  user <- sub(":.*$", "", cred)
  pass <- sub("^[^:]*:", "", cred)
  
  host <- sub(":.*$", "", host_port)
  port <- as.integer(sub("^.*:", "", host_port))
  
  list(user = user, pass = pass, host = host, port = port)
}

PX <- get_proxy_credentials()

append_status <- function(study_id, status, n_cs = NA_integer_, n_variants = NA_integer_,
                          n_l2g = NA_integer_, message_text = NA_character_) {
  dt <- data.table(
    timestamp  = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    studyId    = study_id,
    status     = status,
    n_cs       = n_cs,
    n_variants = n_variants,
    n_l2g      = n_l2g,
    message    = gsub("[\r\n\t]+", " ", as.character(message_text))
  )
  
  fwrite(
    dt,
    STATUS_FILE,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    append = file.exists(STATUS_FILE),
    col.names = !file.exists(STATUS_FILE)
  )
}

collapse_rsids <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}

safe_p <- function(mantissa, exponent) {
  suppressWarnings(as.numeric(mantissa) * 10^as.numeric(exponent))
}

safe_minus_log10_p <- function(mantissa, exponent) {
  mantissa <- suppressWarnings(as.numeric(mantissa))
  exponent <- suppressWarnings(as.numeric(exponent))
  out <- -log10(mantissa) - exponent
  out[is.na(mantissa) | is.na(exponent) | mantissa <= 0] <- NA_real_
  out
}

empty_variants_table <- function() {
  data.table(
    studyId = character(),
    studyLocusId = character(),
    credibleSetIndex = integer(),
    finemappingMethod = character(),
    confidence = character(),
    lead_id = character(),
    lead_rs = character(),
    var_id = character(),
    rsids = character(),
    chr = character(),
    pos = integer(),
    A1 = character(),
    A2 = character(),
    PIP = numeric(),
    beta = numeric(),
    p_mantissa = numeric(),
    p_exponent = numeric(),
    p = numeric(),
    minus_log10_p = numeric(),
    logBF = numeric()
  )
}

empty_l2g_table <- function() {
  data.table(
    studyId = character(),
    studyLocusId = character(),
    gene_id = character(),
    gene_symbol = character(),
    l2g_score = numeric()
  )
}

gql_post <- function(query, variables, max_tries = 5L, sleep_base = 0.7) {
  body <- jsonlite::toJSON(
    list(query = query, variables = variables),
    auto_unbox = TRUE,
    null = "null",
    digits = 22
  )
  
  ua <- httr::user_agent("CLEAR-GCST-batch/1.0")
  attempt <- 1L
  
  while (attempt <= max_tries) {
    res <- httr::POST(
      BASE_URL,
      body = body,
      encode = "raw",
      ua,
      httr::add_headers(`Content-Type` = "application/json"),
      httr::use_proxy(
        url      = PX$host,
        port     = PX$port,
        username = PX$user,
        password = PX$pass
      ),
      httr::timeout(60)
    )
    
    txt <- tryCatch(
      httr::content(res, as = "text", encoding = "UTF-8"),
      error = function(e) ""
    )
    
    code <- httr::status_code(res)
    ok <- !httr::http_error(res)
    
    obj <- if (nzchar(txt)) {
      tryCatch(jsonlite::fromJSON(txt, flatten = TRUE), error = function(e) list(raw = txt))
    } else {
      list()
    }
    
    if (!ok && code %in% c(429, 500, 502, 503, 504) && attempt < max_tries) {
      wait <- sleep_base * attempt
      message("HTTP retry ", attempt, "/", max_tries, " after HTTP ", code)
      Sys.sleep(wait)
      attempt <- attempt + 1L
      next
    }
    
    return(obj)
  }
}

Q_IDS <- '
query ($studyIds:[String!], $index:Int!, $size:Int!) {
  credibleSets(studyIds:$studyIds, page:{index:$index, size:$size}) {
    count
    rows {
      studyLocusId
    }
  }
}'

Q_CS_SINGLE <- '
query ($id:String!, $index:Int!, $size:Int!) {
  credibleSet(studyLocusId:$id) {
    studyId
    studyLocusId
    credibleSetIndex
    finemappingMethod
    confidence

    variant {
      id
      rsIds
    }

    l2GPredictions(page:{index:0, size:50}) {
      rows {
        score
        target {
          id
          approvedSymbol
        }
      }
    }

    locus(page:{index:$index, size:$size}) {
      count
      rows {
        posteriorProbability
        beta
        pValueMantissa
        pValueExponent
        logBF

        variant {
          id
          rsIds
          chromosome
          position
          referenceAllele
          alternateAllele
        }
      }
    }
  }
}'

studies <- fread(STUDIES_PATH)
gcst_studies <- studies[grepl("^GCST", studyId)]

message("Total studies in studies.tsv: ", nrow(studies))
message("GCST studies retained: ", nrow(gcst_studies))

f_gcst_studies <- file.path(OUT_ROOT, "GCST_studies_full.tsv.gz")
fwrite(gcst_studies, f_gcst_studies, sep = "\t", quote = FALSE, na = "NA")
message("Wrote full GCST study table: ", f_gcst_studies)

get_cs_ids <- function(study_id, page_size = PAGE_IDS) {
  out <- character()
  i <- 0L
  total <- Inf
  
  while (i * page_size < total) {
    d <- gql_post(
      Q_IDS,
      list(
        studyIds = list(study_id),
        index = i,
        size = page_size
      )
    )
    
    cs <- d$data$credibleSets
    ids <- as.character(cs$rows$studyLocusId)
    ids <- ids[!is.na(ids) & nzchar(ids)]
    out <- c(out, ids)
    
    total <- suppressWarnings(as.integer(cs$count))
    i <- i + 1L
    Sys.sleep(SLEEP_SEC)
  }
  
  unique(out)
}

pull_one_cs <- function(cs_id, page_size = PAGE_LOCUS) {
  j <- 0L
  total <- Inf
  chunks <- list()
  meta <- NULL
  l2g_tbl <- NULL
  
  while (j * page_size < total) {
    d <- gql_post(
      Q_CS_SINGLE,
      list(
        id = cs_id,
        index = j,
        size = page_size
      )
    )
    
    cs <- d$data$credibleSet

    if (is.null(meta)) {
      lead_rs <- collapse_rsids(cs$variant$rsIds)
      
      meta <- list(
        studyId = cs$studyId,
        studyLocusId = cs$studyLocusId,
        credibleSetIndex = cs$credibleSetIndex,
        finemappingMethod = cs$finemappingMethod,
        confidence = cs$confidence,
        lead_id = cs$variant$id,
        lead_rs = lead_rs
      )
      
      if (!is.null(cs$l2GPredictions$rows) && NROW(cs$l2GPredictions$rows) > 0L) {
        L2 <- as.data.frame(cs$l2GPredictions$rows)
        
        l2g_tbl <- data.table(
          studyId = cs$studyId,
          studyLocusId = cs$studyLocusId,
          gene_id = L2$target.id,
          gene_symbol = L2$target.approvedSymbol,
          l2g_score = as.numeric(L2$score)
        )
      }
    }
    
    loc <- cs$locus
    loc_rows <- loc$rows
    total <- suppressWarnings(as.integer(loc$count))
    
    if (is.data.frame(loc_rows) && nrow(loc_rows) > 0L) {
      chunks[[length(chunks) + 1L]] <- as.data.table(loc_rows)
    }

    j <- j + 1L
    Sys.sleep(SLEEP_SEC)
  }
  
  v <- rbindlist(chunks, fill = TRUE)
  
  rsids_vec <- vapply(v$variant.rsIds, collapse_rsids, character(1))
  
  vars <- data.table(
    studyId = meta$studyId,
    studyLocusId = meta$studyLocusId,
    credibleSetIndex = meta$credibleSetIndex,
    finemappingMethod = meta$finemappingMethod,
    confidence = meta$confidence,
    
    lead_id = meta$lead_id,
    lead_rs = meta$lead_rs,
    
    var_id = v$variant.id,
    rsids = rsids_vec,
    chr = as.character(v$variant.chromosome),
    pos = as.integer(v$variant.position),
    A1 = as.character(v$variant.alternateAllele),
    A2 = as.character(v$variant.referenceAllele),
    
    PIP = as.numeric(v$posteriorProbability),
    beta = as.numeric(v$beta),
    
    p_mantissa = as.numeric(v$pValueMantissa),
    p_exponent = as.numeric(v$pValueExponent),
    p = safe_p(v$pValueMantissa, v$pValueExponent),
    minus_log10_p = safe_minus_log10_p(v$pValueMantissa, v$pValueExponent),
    
    logBF = as.numeric(v$logBF)
  )
  
  list(
    variants = vars,
    l2g = l2g_tbl
  )
}

to_clear_gr <- function(df) {
  df <- as.data.table(df)
  
  if (nrow(df) == 0L) {
    return(GRanges())
  }
  
  rs_first <- ifelse(
    is.na(df$rsids) | !nzchar(df$rsids),
    "",
    sub(";.*$", "", df$rsids)
  )
  
  fallback <- paste0(df$chr, ":", df$pos, "_", df$A2, "_", df$A1)
  nm <- ifelse(nzchar(rs_first), rs_first, fallback)
  
  chr <- ifelse(grepl("^chr", df$chr), df$chr, paste0("chr", df$chr))
  pos <- as.integer(df$pos)
  
  keep <- !is.na(chr) & !is.na(pos) & pos > 0
  
  df <- df[keep]
  chr <- chr[keep]
  pos <- pos[keep]
  nm <- nm[keep]
  
  if (nrow(df) == 0L) {
    return(GRanges())
  }
  
  gr <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = pos, end = pos)
  )
  
  mcols(gr)$names <- nm
  mcols(gr)$signal <- df$studyLocusId
  
  mcols(gr)$studyId <- df$studyId
  mcols(gr)$studyLocusId <- df$studyLocusId
  mcols(gr)$credibleSetIndex <- df$credibleSetIndex
  mcols(gr)$finemappingMethod <- df$finemappingMethod
  mcols(gr)$confidence <- df$confidence
  
  mcols(gr)$lead_id <- df$lead_id
  mcols(gr)$lead_rs <- df$lead_rs
  mcols(gr)$var_id <- df$var_id
  mcols(gr)$rsids <- df$rsids
  
  mcols(gr)$A1 <- df$A1
  mcols(gr)$A2 <- df$A2
  
  mcols(gr)$PIP <- df$PIP
  mcols(gr)$beta <- df$beta
  mcols(gr)$p_mantissa <- df$p_mantissa
  mcols(gr)$p_exponent <- df$p_exponent
  mcols(gr)$p <- df$p
  mcols(gr)$minus_log10_p <- df$minus_log10_p
  mcols(gr)$logBF <- df$logBF
  
  if ("top_gene_id" %in% names(df)) {
    mcols(gr)$top_gene_id <- df$top_gene_id
    mcols(gr)$top_gene_symbol <- df$top_gene_symbol
    mcols(gr)$top_l2g_score <- df$top_l2g_score
  }
  
  gr
}

gr_to_tsv <- function(gr) {
  if (length(gr) == 0L) {
    return(data.table(
      seqnames = character(),
      start = integer(),
      end = integer(),
      names = character(),
      signal = character()
    ))
  }
  
  dt <- data.table(
    seqnames = as.character(seqnames(gr)),
    start = start(gr),
    end = end(gr)
  )
  
  mc <- as.data.table(as.data.frame(mcols(gr)))
  cbind(dt, mc)
}

process_study <- function(study_id, study_row = NULL, skip_if_done = TRUE) {
  sid_dir <- file.path(STUDY_DIR, study_id)
  dir.create(sid_dir, recursive = TRUE, showWarnings = FALSE)
  
  f_cs_ids <- file.path(sid_dir, paste0(study_id, "_studyLocusIds.txt"))
  f_variants <- file.path(sid_dir, paste0(study_id, "_credible_set_variants_full.tsv.gz"))
  f_l2g <- file.path(sid_dir, paste0(study_id, "_credible_set_l2g.tsv.gz"))
  f_vtop <- file.path(sid_dir, paste0(study_id, "_credible_set_variants_with_top_l2g.tsv.gz"))
  f_gr_rds <- file.path(sid_dir, paste0(study_id, "_CLEAR_credible_set_variants.gr.rds"))
  f_gr_tsv <- file.path(sid_dir, paste0(study_id, "_CLEAR_credible_set_variants.tsv.gz"))
  f_study_meta <- file.path(sid_dir, paste0(study_id, "_study_metadata.tsv"))
  f_failed_cs <- file.path(sid_dir, paste0(study_id, "_failed_studyLocusIds.txt"))
  
  if (skip_if_done && all(file.exists(f_variants, f_l2g, f_vtop, f_gr_rds, f_gr_tsv))) {
    message("  already done: ", study_id)
    append_status(study_id, "skipped_already_done")
    return(invisible(TRUE))
  }
  
  message("Processing: ", study_id)
  
  if (!is.null(study_row)) {
    fwrite(as.data.table(study_row), f_study_meta, sep = "\t", quote = FALSE, na = "NA")
  }
  
  cs_ids <- get_cs_ids(study_id)
  
  writeLines(cs_ids, f_cs_ids)
  message("  credible sets: ", length(cs_ids))
  
  if (!length(cs_ids)) {
    fwrite(empty_variants_table(), f_variants, sep = "\t", quote = FALSE, na = "NA")
    fwrite(empty_l2g_table(), f_l2g, sep = "\t", quote = FALSE, na = "NA")
    fwrite(empty_variants_table(), f_vtop, sep = "\t", quote = FALSE, na = "NA")
    saveRDS(GRanges(), f_gr_rds)
    fwrite(gr_to_tsv(GRanges()), f_gr_tsv, sep = "\t", quote = FALSE, na = "NA")
    
    append_status(study_id, "no_credible_sets", n_cs = 0L, n_variants = 0L, n_l2g = 0L)
    return(invisible(TRUE))
  }
  
  parts <- vector("list", length(cs_ids))
  failed_cs <- character()
  
  for (k in seq_along(cs_ids)) {
    csid <- cs_ids[k]
    
    parts[[k]] <- tryCatch(
      pull_one_cs(csid),
      error = function(e) {
        msg <- conditionMessage(e)
        message("    failed CS: ", csid, " | ", msg)
        failed_cs <<- c(failed_cs, paste(csid, msg, sep = "\t"))
        NULL
      }
    )
    
    Sys.sleep(SLEEP_SEC)
  }
  
  parts <- Filter(Negate(is.null), parts)
  
  if (length(failed_cs) > 0L) {
    writeLines(failed_cs, f_failed_cs)
  } else if (file.exists(f_failed_cs)) {
    file.remove(f_failed_cs)
  }
  
  variants <- rbindlist(lapply(parts, `[[`, "variants"), fill = TRUE)
  
  l2g_list <- Filter(Negate(is.null), lapply(parts, `[[`, "l2g"))
  l2g <- if (length(l2g_list)) {
    rbindlist(l2g_list, fill = TRUE)
  } else {
    empty_l2g_table()
  }
  
  fwrite(variants, f_variants, sep = "\t", quote = FALSE, na = "NA")
  fwrite(l2g, f_l2g, sep = "\t", quote = FALSE, na = "NA")
  
  if (nrow(l2g) > 0L) {
    top_l2g <- l2g[
      order(studyLocusId, -l2g_score),
      .SD[1],
      by = studyLocusId
    ][
      ,
      .(
        studyLocusId,
        top_gene_id = gene_id,
        top_gene_symbol = gene_symbol,
        top_l2g_score = l2g_score
      )
    ]
    
    variants_top <- merge(
      variants,
      top_l2g,
      by = "studyLocusId",
      all.x = TRUE
    )
  } else {
    variants_top <- copy(variants)
    variants_top[, `:=`(
      top_gene_id = NA_character_,
      top_gene_symbol = NA_character_,
      top_l2g_score = NA_real_
    )]
  }
  
  fwrite(variants_top, f_vtop, sep = "\t", quote = FALSE, na = "NA")
  
  gr <- to_clear_gr(variants_top)
  
  saveRDS(gr, f_gr_rds)
  
  gr_tsv <- gr_to_tsv(gr)
  fwrite(gr_tsv, f_gr_tsv, sep = "\t", quote = FALSE, na = "NA")
  
  status <- if (length(failed_cs) > 0L) "done_partial_failed_cs" else "done"
  
  append_status(
    study_id,
    status,
    n_cs = length(cs_ids),
    n_variants = nrow(variants),
    n_l2g = nrow(l2g),
    message_text = if (length(failed_cs) > 0L) paste0(length(failed_cs), " credible sets failed") else NA_character_
  )
  
  message("  wrote:")
  message("    ", f_variants)
  message("    ", f_l2g)
  message("    ", f_vtop)
  message("    ", f_gr_rds)
  message("    ", f_gr_tsv)
  
  invisible(TRUE)
}

message("Starting GCST processing...")

for (i in seq_len(nrow(gcst_studies))) {
  sid <- gcst_studies$studyId[i]

  process_study(
    study_id = sid,
    study_row = gcst_studies[i],
    skip_if_done = TRUE
  )
  
  Sys.sleep(SLEEP_SEC)
}

message("Done.")
