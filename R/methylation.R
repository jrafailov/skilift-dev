#' lift_methylation
#'
#' Create methylation tracks for each sample in a cohort.
#' Expects a cohort input column (default: "methylation") containing either
#' a GRanges object or a path to an RDS file that loads to a GRanges.
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param methylation_col Cohort input column containing methylation GRanges or RDS path
#' @param beta_field Field to plot for methylation beta values
#' # @param intensity_field Optional field for methylation intensity track
#' @param cores Number of cores for parallel processing
#' @return Modified cohort object (invisibly)
#' @export
lift_methylation <- function(
  cohort,
  output_data_dir,
  methylation_col = "methylation",
  beta_field = "p_methyl",
  # intensity_field = NULL,
  cores = 1
) {
  if (!inherits(cohort, "Cohort")) {
    stop("Input must be a Cohort object")
  }
  if (!dir.exists(output_data_dir)) {
    dir.create(output_data_dir, recursive = TRUE)
  }
  if (!methylation_col %in% names(cohort$inputs)) {
    warning("Missing methylation column in cohort inputs: ", methylation_col)
    return(invisible(cohort))
  }
  
  # Get reference from cohort
  reference_name <- cohort$reference_name
  if (is.null(reference_name)) {
    stop("Reference name not found in cohort object")
  }

  # choose_intensity_field <- function(gr, beta_field) {
  #   candidate_fields <- c("ratio", "intensity", "signal", "coverage", beta_field)
  #   candidate_fields <- candidate_fields[candidate_fields %in% names(S4Vectors::mcols(gr))]
  #   if (length(candidate_fields) == 0) return(NA_character_)
  #   candidate_fields[[1]]
  # }

  iterate_fun <- function(i) {
    main <- function() {
      row <- cohort$inputs[i, ]
      pair_dir <- file.path(output_data_dir, row$pair)
      if (!dir.exists(pair_dir)) {
        dir.create(pair_dir, recursive = TRUE)
      }

      methyl_source <- row[[methylation_col]]
      if (is.null(methyl_source) || (length(methyl_source) == 1 && is.na(methyl_source))) {
        return(NULL)
      }

      gr <- methyl_source
      if (is.character(methyl_source) && length(methyl_source) == 1) {
        if (!file.exists(methyl_source)) {
          warning("Methylation RDS path does not exist for ", row$pair, ": ", methyl_source)
          return(NULL)
        }
        gr <- readRDS(methyl_source)
      }

      if (!inherits(gr, "GRanges")) {
        stop("Methylation input must be GRanges for ", row$pair)
      }

      if (!beta_field %in% names(S4Vectors::mcols(gr))) {
        stop("Missing beta field '", beta_field, "' for ", row$pair)
      }

      # # Get metadata for purity/ploidy (only needed for intensity track)
      # metadata_path <- file.path(pair_dir, "metadata.json")
      # js <- NULL
      # purity <- NULL
      # ploidy <- NULL
      # if (file.exists(metadata_path)) {
      #   js <- jsonlite::fromJSON(metadata_path, simplifyVector = FALSE)
      #   purity <- js[[1]]$purity
      #   ploidy <- js[[1]]$ploidy
      # }

      # # Create intensity track if field exists
      # intensity_field_use <- intensity_field
      # if (is.null(intensity_field_use) || is.na(intensity_field_use)) {
      #   intensity_field_use <- choose_intensity_field(gr, beta_field)
      # }
      # 
      # if (!is.na(intensity_field_use) && intensity_field_use %in% names(S4Vectors::mcols(gr))) {
      #   # Compute rel2abs params if purity/ploidy available
      #   params <- NULL
      #   if (!is.null(purity) && !is.null(ploidy) && is.finite(purity) && is.finite(ploidy)) {
      #     tryCatch({
      #       params <- skitools::rel2abs(
      #         gr,
      #         field = intensity_field_use,
      #         purity = purity,
      #         ploidy = ploidy,
      #         return.params = TRUE
      #       )
      #     }, error = function(e) {
      #       warning("Failed to compute rel2abs params for ", row$pair, ": ", e$message)
      #     })
      #   }
      #   
      #   methylation_intensity <- Skilift:::granges_to_arrow_scatterplot(
      #     gr,
      #     field = intensity_field_use,
      #     ref = reference_name,
      #     bin.width = NA,
      #     mask = FALSE,
      #     arrow_x_type = arrow::float64
      #   )
      #   
      #   arrow::write_feather(
      #     methylation_intensity,
      #     file.path(pair_dir, "methylation_intensity.arrow"),
      #     compression = "uncompressed"
      #   )
      #   
      #   # Update metadata with slope/intercept
      #   if (!is.null(js) && !is.null(params)) {
      #     js[[1]]$methylation_intensity_cov_slope <- params["slope"]
      #     js[[1]]$methylation_intensity_cov_intercept <- params["intercept"]
      #   }
      # }

      # Compute beta colors

      #browser()

      beta_vals <- S4Vectors::mcols(gr)[[beta_field]]
      ramp <- grDevices::colorRamp(c("blue", "red"))
      normvals <- pmin(pmax(beta_vals, 0), 1)
      colmat <- ramp(normvals)
      cols <- grDevices::rgb(colmat, maxColorValue = 255)
      S4Vectors::mcols(gr)$cols <- cols
      GenomicRanges::ranges(gr) <- IRanges::IRanges(start = GenomicRanges::start(gr), width = 1)
      gr <- gr %Q% (seqnames %in% c(1:22, "X", "Y"))

      if (reference_name == "hg19") {
        cpg.sites = system.file("extdata/data", "cpg_positions_hg19_sorted.bed", package = "Skilift")
        cpg.sites = data.table::fread(cpg.sites, header = FALSE)
        setnames(cpg.sites, c("chrom", "start", "end", "name", "hell", "strand", "cpg_code"))
        cpg.sites = cpg.sites %>% dt2gr()

        gr = gr %&% cpg.sites

      }


      # Create beta track
      methylation_beta <- Skilift:::granges_to_arrow_scatterplot(
        gr,
        field = beta_field,
        ref = reference_name,
        bin.width = NA,
        cov.color.field = "cols",
        mask = FALSE,
        arrow_x_type = arrow::float64
      )

      arrow::write_feather(
        methylation_beta,
        file.path(pair_dir, "methylation_beta.arrow"),
        compression = "uncompressed"
      )

      # # Write updated metadata (only needed for intensity track)
      # if (!is.null(js)) {
      #   jsonlite::write_json(
      #     js,
      #     metadata_path,
      #     auto_unbox = TRUE,
      #     pretty = TRUE,
      #     null = "null"
      #   )
      # }

      return(NULL)
    }

    futile.logger::flog.threshold("ERROR")
    tryCatchLog::tryCatchLog(
      main(),
      error = function(e) {
        pair_name <- cohort$inputs$pair[i]
        print(sprintf("Error processing %s: %s", pair_name, e$message))
        NULL
      }
    )
  }

  parallel::mclapply(seq_len(nrow(cohort$inputs)), iterate_fun, mc.cores = cores, mc.preschedule = TRUE)
  invisible(cohort)
}
