
#### Version 2 ####
library(httr)
library(jsonlite)
library(data.table)
library(dplyr)
library(GenomicRanges)
library(S4Vectors)

setwd('/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/')

proxy_url <- "http://lefei.wang:Mrhuai990116&7@webproxy.adqimr.ad.lan:8080/"

Sys.setenv(
  http_proxy  = proxy_url,
  https_proxy = proxy_url,
  HTTP_PROXY  = proxy_url,
  HTTPS_PROXY = proxy_url
)

r <- GET(
  "https://api.platform.opentargets.org/api/v4/graphql",
  use_proxy(
    url = "webproxy.adqimr.ad.lan",
    port = 8080,
    username = "lefei.wang",
    password = "Mrhuai990116&7"
  ),
  timeout(30)
)

status_code(r)
content(r, "text", encoding = "UTF-8")

#  config 
BASE_URL    <- "https://api.platform.opentargets.org/api/v4/graphql"
OUT_ROOT    <- "opentarget_trait"
PAGE_IDS    <- 200      # page size for set IDs
PAGE_LOCUS  <- 200      # page size for variants
SLEEP_SEC   <- 0.15

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

study_map <- c(
  "GCST90002344" = "monocyte_count_2020_Chen_cs374",
  "GCST004988"   = "breast_carcinoma_2017_Michailidou_cs210",
  "GCST010099"   = "breast_carcinoma_2020_Zhang_cs10",
  "GCST010100"   = "tnbc_brca1_2020_Zhang_cs18",
  "GCST010797"   = "bc_ovarian_prostate_pleiotropy_2016_Kar_cs18",
  "GCST90011732" = "mammographic_density_nonDenseArea_2020_Sieh_cs11",
  "GCST90011731" = "mammographic_density_DenseArea_2020_Sieh_cs10",
  "GCST90027158" = "alzheimer_2022_bellenguez_cs83",
  "GCST90479813" = "brain_and_nervous_cancer_2024_verma_cs52",
  "GCST90479812" = "brain_cancer_2024_verma_cs84",
  "GCST90435268" = "colorectal_cancer_2024_tian_2024_cs70",
  "GCST90455657" = "heart_failure_2025_lee_cs165",
  "GCST90320054" = "kidney_cancer_2024_purdue_2024_cs89",
  "GCST90002357" = "platelet_count_2020_chen_cs1514",
  "GCST90492734" = "t2d_2024_suzuki_cs687",
  "FINNGEN_R12_HEIGHT_IRN" = "height_finngen_cs1903",
  "GCST90310178" = "sitting_height_2025_Hu_cs686",
  "GCST90016665" = "ovarian_cancer_2022_Dareng_cs27",
  "FINNGEN_R12_CD2_INSITU_BREAST_INTRADUCTAL_EXALLC" = "dcis_2023_Kurki_cs8",
  "GCST90292538" = "IBD_2024_Liu_cs293",
  "GCST90255621" = "BMI_2022_Huang_cs1071",
  "GCST90310294" = "systolic_blood_pressure_2024_Keaton_cs1162",
  "GCST90310295" = "diastolic_blood_pressure_2024_Keaton_cs1108",
  "GCST90132314" = "coronary_artery_disease_2022_Aragam_cs250",
  "GCST90105038" = "educational_attainment_2022_Okbay_cs1599",
  "GCST90503210" = "schizophrenia_2025_Dang_cs247"
  
)

slug <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", tolower(x))
  sub("^_+|_+$", "", x)
}

#  graphql 
gql_post <- function(query, variables) {
  body <- toJSON(
    list(query = query, variables = variables),
    auto_unbox = TRUE, null = "null", digits = 22
  )
  
  res <- POST(
    BASE_URL,
    body = body,
    encode = "raw",
    add_headers(`Content-Type` = "application/json"),
    use_proxy(
      url = "webproxy.adqimr.ad.lan",
      port = 8080,
      username = "lefei.wang",
      password = "Mrhuai990116&7"
    ),
    timeout(60)
  )
  
  txt <- content(res, "text", encoding = "UTF-8")
  obj <- tryCatch(fromJSON(txt, flatten = TRUE), error = function(e) list(raw = txt))
  
  if (http_error(res)) stop(sprintf("HTTP %s\n%s", status_code(res), txt))
  if (!is.null(obj$errors)) {
    msgs <- vapply(obj$errors, function(e) e$message, character(1))
    stop(paste(msgs, collapse = "\n"))
  }
  
  obj
}

Q_IDS <- '
query ($studyIds:[String!], $index:Int!, $size:Int!) {
  credibleSets(studyIds:$studyIds, page:{index:$index, size:$size}) {
    count
    rows { studyLocusId }
  }
}'

Q_VARIANTS <- '
query ($ids:[String!]!, $index:Int!, $size:Int!) {
  credibleSets(studyLocusIds:$ids) {
    rows {
      studyId
      studyLocusId
      credibleSetIndex
      finemappingMethod
      confidence
      variant { id rsIds }

      l2GPredictions(page:{index:0, size:50}) {
        rows { score target { id approvedSymbol } }
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
            id rsIds chromosome position referenceAllele alternateAllele
          }
        }
      }
    }
  }
}'

#  helpers 
get_set_ids <- function(study_id, page_size = PAGE_IDS) {
  out <- character()
  i <- 0L
  repeat {
    d <- gql_post(Q_IDS, list(studyIds = list(study_id), index = i, size = page_size))
    cs <- d$data$credibleSets
    out <- c(out, unlist(cs$rows$studyLocusId))
    if ((i + 1L) * page_size >= cs$count) break
    i <- i + 1L
    Sys.sleep(SLEEP_SEC)
  }
  unique(out)
}

pull_set <- function(study_locus_id, page_size = PAGE_LOCUS) {
  j <- 0L
  chunks <- list()
  meta <- NULL
  l2g <- NULL
  
  repeat {
    d <- gql_post(Q_VARIANTS, list(ids = list(study_locus_id), index = j, size = page_size))
    rows <- d$data$credibleSets$rows
    if (is.null(rows) || NROW(rows) == 0L) return(NULL)
    r <- rows[1, , drop = FALSE]
    
    # once per set
    if (is.null(meta)) {
      lead_rs <- if (!is.null(r$variant.rsIds[[1]]) && length(r$variant.rsIds[[1]]) > 0)
        paste(r$variant.rsIds[[1]], collapse = ";") else NA_character_
      
      meta <- list(
        studyId           = r$studyId,
        studyLocusId      = r$studyLocusId,
        credibleSetIndex  = r$credibleSetIndex,
        finemappingMethod = r$finemappingMethod,
        confidence        = r$confidence,
        lead_id           = r$variant.id,
        lead_rs           = lead_rs
      )
      
      if (!is.null(r$l2GPredictions.rows[[1]]) && NROW(r$l2GPredictions.rows[[1]]) > 0L) {
        L2 <- r$l2GPredictions.rows[[1]]
        l2g <- tibble::tibble(
          studyId      = r$studyId,
          studyLocusId = r$studyLocusId,
          gene_id      = L2$target.id,
          gene_symbol  = L2$target.approvedSymbol,
          l2g_score    = L2$score
        )
      }
    }
    
    loc_rows <- r$locus.rows[[1]]
    if (is.data.frame(loc_rows) && nrow(loc_rows) > 0L) {
      chunks[[length(chunks) + 1L]] <- loc_rows
    }
    
    total <- suppressWarnings(as.integer(r$locus.count))
    enough <- (!is.na(total) && ((j + 1L) * page_size >= total)) ||
      (is.data.frame(loc_rows) && nrow(loc_rows) < page_size) ||
      is.null(loc_rows)
    if (enough) break
    
    j <- j + 1L
    Sys.sleep(SLEEP_SEC)
  }
  
  if (!length(chunks)) return(NULL)
  v <- dplyr::bind_rows(chunks)
  
  vars <- tibble::tibble(
    studyId          = meta$studyId,
    studyLocusId     = meta$studyLocusId,
    credibleSetIndex = meta$credibleSetIndex,
    finemappingMethod= meta$finemappingMethod,
    confidence       = meta$confidence,
    lead_id          = meta$lead_id,
    lead_rs          = meta$lead_rs,
    var_id           = v$variant.id,
    rsids            = vapply(v$variant.rsIds, function(x) if (length(x)) paste(x, collapse=";") else NA_character_, ""),
    chr              = v$variant.chromosome,
    pos              = v$variant.position,
    A1               = v$variant.alternateAllele,
    A2               = v$variant.referenceAllele,
    PIP              = v$posteriorProbability,
    beta             = v$beta,
    p                = v$pValueMantissa * 10^v$pValueExponent,
    logBF            = v$logBF
  )
  
  list(variants = vars, l2g = l2g)
}

to_clear_gr <- function(df) {
  rs_first <- ifelse(is.na(df$rsids), "", sub(";.*$", "", df$rsids))
  fallback <- paste0(df$chr, ":", df$pos, "_", df$A2, "_", df$A1)
  nm <- ifelse(nzchar(rs_first), rs_first, fallback)
  
  chr <- ifelse(grepl("^chr", df$chr), df$chr, paste0("chr", df$chr))
  pos <- as.integer(df$pos)
  
  # control for NA positions: replace with 1
  pos[is.na(pos)] <- 1L
  
  gr <- GRanges(seqnames = chr,
                ranges = IRanges(pmax(pos - 1L, 1L), pos))
  mcols(gr)$names  <- nm
  mcols(gr)$signal <- df$lead_id
  gr
}


is_study_done <- function(study_id, name_readable) {
  tag     <- slug(name_readable)
  out_dir <- file.path(OUT_ROOT, paste0(tag, "__", study_id))
  f_vars_top <- file.path(out_dir, sprintf("%s__%s_credible_set_variants_with_top_l2g.tsv", tag, study_id))
  f_gr_rds   <- file.path(out_dir, sprintf("%s__%s_CLEAR_credible_set_variants.gr.rds", tag, study_id))
  # consider it "done" if the key downstream files exist
  dir.exists(out_dir) && file.exists(f_vars_top) && file.exists(f_gr_rds)
}
#  one study 
process_study <- function(study_id, name_readable, skip_if_done = TRUE) {
  tag <- slug(name_readable)
  out_dir <- file.path(OUT_ROOT, paste0(tag, "__", study_id))
  
  if (skip_if_done && is_study_done(study_id, name_readable)) {
    message("== ", study_id, " | ", name_readable, " already extracted. Skipping.")
    return(invisible(list(skipped = TRUE, out_dir = out_dir)))
  }
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("== ", study_id, " | ", name_readable)
  
  ids <- get_set_ids(study_id)
  message("  credible sets: ", length(ids))
  
  parts <- lapply(ids, function(x) {
    tryCatch(pull_set(x), error = function(e) { message("  failed set: ", x); NULL })
  })
  
  variants <- dplyr::bind_rows(lapply(parts, `[[`, "variants"))
  l2g      <- dplyr::bind_rows(lapply(parts, `[[`, "l2g"))
  
  f_vars <- file.path(out_dir, sprintf("%s__%s_credible_set_variants.tsv", tag, study_id))
  fwrite(variants, f_vars, sep = "\t", quote = FALSE, na = "NA")
  
  f_l2g  <- file.path(out_dir, sprintf("%s__%s_credible_set_l2g.tsv", tag, study_id))
  if (!is.null(l2g) && nrow(l2g) > 0) {
    fwrite(l2g, f_l2g, sep = "\t", quote = FALSE, na = "NA")
  } else {
    fwrite(data.frame(studyId=character(), studyLocusId=character(),
                      gene_id=character(), gene_symbol=character(), l2g_score=numeric()),
           f_l2g, sep = "\t", quote = FALSE, na = "NA")
  }
  
  if (!is.null(l2g) && nrow(l2g) > 0) {
    top <- l2g %>%
      dplyr::group_by(studyLocusId) %>%
      dplyr::slice_max(l2g_score, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(studyLocusId,
                    top_gene_id = gene_id,
                    top_gene_symbol = gene_symbol,
                    top_l2g_score = l2g_score)
    variants_top <- dplyr::left_join(variants, top, by = "studyLocusId")
  } else {
    variants_top <- dplyr::mutate(variants,
                                  top_gene_id = NA_character_,
                                  top_gene_symbol = NA_character_,
                                  top_l2g_score = NA_real_)
  }
  
  f_vars_top <- file.path(out_dir, sprintf("%s__%s_credible_set_variants_with_top_l2g.tsv", tag, study_id))
  fwrite(variants_top, f_vars_top, sep = "\t", quote = FALSE, na = "NA")
  
  writeLines(ids, file.path(out_dir, sprintf("%s__%s_studyLocusIds.txt", tag, study_id)))
  
  gr <- to_clear_gr(variants)
  f_gr_rds <- file.path(out_dir, sprintf("%s__%s_CLEAR_credible_set_variants.gr.rds", tag, study_id))
  saveRDS(gr, f_gr_rds)
  
  f_gr_tsv <- file.path(out_dir, sprintf("%s__%s_CLEAR_credible_set_variants.tsv", tag, study_id))
  write.table(
    data.frame(
      seqnames = as.character(GenomeInfoDb::seqnames(gr)),
      start    = IRanges::start(gr),
      end      = IRanges::end(gr),
      names    = S4Vectors::mcols(gr)$names,
      signal   = S4Vectors::mcols(gr)$signal,
      stringsAsFactors = FALSE
    ),
    file = f_gr_tsv, sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  message("  wrote:\n    ", f_vars,
          "\n    ", f_l2g,
          "\n    ", f_vars_top,
          "\n    ", f_gr_rds,
          "\n    ", f_gr_tsv)
  
  invisible(list(variants = variants, l2g = l2g, gr = gr, ids = ids, out_dir = out_dir))
}

# run all
for (sid in names(study_map)) {
  process_study(sid, study_map[[sid]])
}







