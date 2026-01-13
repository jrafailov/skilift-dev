#' Converts a color from hexadecimal to numeric format
#'
#' @param color A character vector of colors in hexadecimal format (e.g., "#FF0000").
#' @return A numeric vector representing the colors.
color2numeric <- function(color) {
    # Remove the '#' character if present
    color <- gsub("#", "", color)
    # Convert hexadecimal to numeric
    return(as.numeric(strtoi(color, base = 16L)))
}


#'
#' @param ref_seqinfo_json Path to the JSON file containing metadata.
#' @param ref Optional reference name. If not provided, the default reference is used.
#' @return A data.table containing the sequence information for the specified reference.
get_ref_metadata <- function(ref_seqinfo_json, ref = NULL) {
    meta <- jsonlite::read_json(ref_seqinfo_json)

    if (!("coordinates" %in% names(meta))) {
        stop("Input meta file is not a proper settings.json format.")
    }

    coord <- meta$coordinates

    if (is.null(ref)) {
        if (!("default" %in% names(coord))) {
            stop(ref_seqinfo_json, " is missing a default coordinates attribute")
        }
        ref <- coord$default
        message('No reference name provided so using default: "', ref, '".')
    }

    if (!("sets" %in% names(coord))) {
        stop("Could not find reference sets in ", ref_seqinfo_json, "in attribute 'sets'.")
    }

    sets <- coord$sets

    if (!(ref %in% names(sets))) {
        stop(ref, " does not appear in set of references in ", ref_seqinfo_json)
    }

    seq.info <- data.table::rbindlist(sets[[ref]])
    return(seq.info)
}

#' @name granges_to_arrow_scatterplot
#' @description
#'
#' Converts GRanges into an arrow file
#'
#' @param gr_path input GRanges
#' @param field name of GRanges column to use as arrow Y axis
#' @param ref the name of the reference to use.
#' @param cov.color.field name of GRanges column containing color of each point
#' @param overwrite (logical) by default, if the output path already exists, it will not be overwritten.
#' @param ref_seqinfo_json path to JSON file with ref seqlengths
#' @param bin.width (integer) bin width for rebinning the coverage (default: 1e4)
#' @param mask (logical) whether to mask the coverage data (default: TRUE)
#' @param mask.path path to mask file (default: system.file("extdata", "data", "maskA_re.rds", package = "Skilift"))
#' @param ... additional arguments to pass to cov2cov.js
#' @param arrow_x_type arrow type function for x values (default: arrow::float32)
#' @param arrow_y_type arrow type function for y values (default: arrow::float32)
#' @return arrow table
#' @author Alon Shaiber, Shihab Dider
granges_to_arrow_scatterplot <- function(
    gr_path,
    field = "foreground",
    ref = "hg19",
    cov.color.field = NULL,
    ref_seqinfo_json = Skilift:::default_settings_path,
    bin.width = 1e4,
    mask = TRUE,
    mask.path = system.file("extdata", "data", "maskA_re.rds", package = "Skilift"),
    arrow_x_type = arrow::float32,
    arrow_y_type = arrow::float32,

    ...) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
        stop('You must have the package "arrow" installed in order for this function to work. Please install it.')
    }

    #' gr_path must be a path or an actual gRanges object    
    if (is.character(gr_path)) {
        if (!file.exists(gr_path)) {
            stop("Please provide a valid path to a GRanges object.")
        }
    } else if (!inherits(gr_path, "GRanges")) {
        stop("Please provide a valid GRanges object.")
    }

    gr <- gr_path
    if (is.character(gr_path)){
        gr <- readRDS(gr_path)
    } 

    if (mask) {
        mask <- readRDS(mask.path)
        gr <- gr.val(gr, mask, "mask")
        gr <- gr %Q% (is.na(mask))
        gr$mask <- NULL
    }

    message("Preparing GRanges for conversion to arrow format")
    dat <- gGnome::cov2cov.js(
        gr,
        meta.js = ref_seqinfo_json,
        js.type = "PGV",
        field = field,
        ref = ref,
        cov.color.field = cov.color.field,
        bin.width = bin.width,
        ...
    )

    if (!is.null(cov.color.field)) {
        dat[, color := color2numeric(get(cov.color.field))]
    } else {
        if (!is.null(ref_seqinfo_json)) {
            ref_meta <- get_ref_metadata(ref_seqinfo_json, ref)
            data.table::setkey(ref_meta, "chromosome")
            dat$color <- color2numeric(ref_meta[dat$seqnames]$color)
        } else {
            dat$color <- 0
        }
    }

    outdt <- dat[, .(x = new.start, y = get(field), color)]

    # if there are any NAs for colors then set those to black
    outdt[is.na(color), color := 0]

    # remove NAs
    outdt <- outdt[!is.na(y)]

    # sort according to x values (that is what front-end expects)
    outdt <- outdt[order(x)]

    ## outdt$x = bit64::as.integer64(outdt$x)
    arrow_table <- arrow::Table$create(
        outdt,
        schema = arrow::schema(
            # x = arrow::float32(),
            # y = arrow::float32(),
            x = arrow_x_type(),
            y = arrow_y_type(),
            color = arrow::float32()
        )
    )
    return(arrow_table)
}

#' @name make_allelic_hets
#' @title make_allelic_hets
#'
#' @description
#'
#' returns melted hetsnps gRanges with allelic (major/minor) counts
#'
#' @param het_pileups (character) path to sites.txt with hetsnps
#' @param min_normal_freq (numeric) in [0, 1] min frequency in normal to count as het site
#' @param max_normal_freq (numeric) in [0, 1] max frequency in normal to count as het site
#'
#' @return GRanges with major and minor counts
#' @author Jonathan Rafailov, Shihab Dider
make_allelic_hetsnps <- function(
    het_pileups = NULL,
    min_normal_freq = 0.2,
    max_normal_freq = 0.8) {
    if (is.null(het_pileups) || !file.exists(het_pileups)) {
        stop("Please provide a valid path to a hetsnps file.")
    }

    ## prepare and filter
    hetsnps_dt <- fread(het_pileups)

    if (all(c("alt.frac.n", "alt.frac.t") %in% colnames(hetsnps_dt))) {
        hetsnps_dt <- hetsnps_dt[
            alt.frac.n > min_normal_freq & alt.frac.n < max_normal_freq,
        ]
    }
    # hetsnps_dt <- fread(hetsnps_dt)[alt.frac.n > min_normal_freq & alt.frac.n < max_normal_freq, ]
    ## add major and minor
    hetsnps_dt[, which.major := ifelse(alt.count.t > ref.count.t, "alt", "ref")]
    hetsnps_dt[, major.count := ifelse(which.major == "alt", alt.count.t, ref.count.t)]
    hetsnps_dt[, minor.count := ifelse(which.major == "alt", ref.count.t, alt.count.t)]

    ## melt the data frame
    hetsnps_melted <- rbind(
        hetsnps_dt[, .(seqnames, start, end, count = major.count, allele = "major")],
        hetsnps_dt[, .(seqnames, start, end, count = minor.count, allele = "minor")]
    )

    ## make GRanges
    allelic_hetsnps <- dt2gr(hetsnps_melted[, .(seqnames, start, end, count, allele)])

    return(allelic_hetsnps)
}


#' @name subsample_hetsnps
#' @title subsample_hetsnps
#' @description subsets the hetsnps to masked unique sites and colors by major/minor allele
#'
#' @param het_pileups sites.txt from het pileups
#' @param mask rds file with masked regions
#' @param sample_size number of snps to randomly sample
#' @export
#' @author Jonathan Rafailov, Shihab Dider
subsample_hetsnps <- function(
    het_pileups,
    mask = NULL,
    sample_size = 100000,
    min_normal_freq = 0.2,
    max_normal_freq = 0.8) {
    if (is.null(het_pileups)) {
        stop("Please provide a valid path to a hetsnps file.")
    }

    if (is.null(mask)) {
        warning("No mask provided, using default mask.")
        maska_path <- system.file("extdata", "data", "maskA_re.rds", package = "Skilift")
        maska <- readRDS(maska_path)
    } else {
        maska <- readRDS(mask)
    }

    # apply mask
    allelic_hetsnps <- make_allelic_hetsnps(
        het_pileups,
        min_normal_freq = 0.2,
        max_normal_freq = 0.8
    )
    allelic_hetsnps <- gr.val(allelic_hetsnps, maska, "mask")
    allelic_hetsnps <- allelic_hetsnps %Q% (is.na(mask))
    allelic_hetsnps$mask <- NULL

    # color by allele
    allelic_hetsnps$col <- ifelse(allelic_hetsnps$allele == "major", "red", "blue")

    # subset to unique sites
    unique_snps <- unique(gr2dt(allelic_hetsnps)[, .(seqnames, start, end)])
    n_snps <- nrow(unique_snps)
    message(paste(n_snps, "snps found"))

    # subsample to reduce number of points drawn on front-end scatterplot
    if (!is.na(sample_size) && n_snps > sample_size) {
        message(paste("subsampling", sample_size, "points..."))
        set.seed(42) ## should set seed for reproducibility
        snps_to_include <- unique_snps[sample(n_snps, sample_size)] %>% dt2gr()
        subsampled_allelic_hetsnps <- allelic_hetsnps %&% snps_to_include
    } else {
        snps_to_include <- unique_snps %>% dt2gr()
        subsampled_allelic_hetsnps <- allelic_hetsnps %&% snps_to_include
    }

    return(subsampled_allelic_hetsnps)
}

#' @name lift_denoised_coverage
#' @title lift_denoised_coverage
#' @description
#' Create denoised coverage arrow files for all samples in a cohort
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @param mask Logical indicating whether to mask the coverage data (default: TRUE)
#' @return None
#' @export
lift_denoised_coverage <- function(
    cohort,
    output_data_dir,
    cores = 1
) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }

    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }

    # Validate required columns exist
    required_cols <- c("pair", "tumor_coverage")
    missing_cols <- required_cols[!required_cols %in% names(cohort$inputs)]
    if (length(missing_cols) > 0) {
        stop("Missing required columns in cohort: ", paste(missing_cols, collapse = ", "))
    }

    # Process each sample in parallel
    mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        row <- cohort$inputs[i, ]
        pair_dir <- file.path(output_data_dir, row$pair)

        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }

        out_file <- file.path(pair_dir, "coverage.arrow")

        futile.logger::flog.threshold("ERROR")
        tryCatchLog(
            {
                if (!is.null(row$tumor_coverage) && file.exists(row$tumor_coverage)) {

                    cov <- row$tumor_coverage %>% readRDS()
					coverage_values = base::get(row$denoised_coverage_field, as.environment(as.list(mcols(cov))))
					lst_cov_bool = Skilift::test_coverage_normalized(coverage_values)
					is_cov_likely_normalized = lst_cov_bool$is_cov_likely_normalized
					is_cov_near_one = lst_cov_bool$is_cov_near_one
					if (!is_cov_likely_normalized && !is_cov_near_one) {
						message("Assuming coverage is in read coverage per bin and paired end, 151 bp reads, rescaling")
						mcols(cov)[[row$denoised_coverage_field]] = coverage_values * PAIRED_READS_FACTOR * READ_LENGTH / width(cov)
					} else {
						message("Assuming coverage is mean-normalized, ignoring rescaling to base coverage")
					}
                    # mcols(cov)[[row$denoised_coverage_field]] <- mcols(cov)[, row$denoised_coverage_field] * 2 * 151 / width(cov)
                                        
                    # Create arrow table
                    arrow_table <- granges_to_arrow_scatterplot(
                        gr_path = cov,
                        field = row$denoised_coverage_field,
                        ref = cohort$reference_name,
                        cov.color.field = row$denoised_coverage_color_field,
                        bin.width = row$denoised_coverage_bin_width,
                        mask = row$denoised_coverage_apply_mask
                    )

                    # Write arrow table
                    message(sprintf("Writing coverage arrow file for %s", row$pair))
                    arrow::write_feather(arrow_table, out_file, compression = "uncompressed")
                } else {
                    warning(sprintf("Tumor coverage file missing for %s", row$pair))
                }
            },
            error = function(e) {
                print(sprintf("Error processing %s: %s", row$pair, e$message))
                NULL
            }
        )
    }, mc.cores = cores, mc.preschedule = TRUE)

    invisible(NULL)
}

#' @name lift_hetsnps
#' @title lift_hetsnps
#' @description
#' Create hetsnps arrow files for all samples in a cohort
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @return None
#' @export
lift_hetsnps <- function(cohort, output_data_dir, cores = 1) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }

    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }

    # Validate required columns exist
    required_cols <- c("pair", "het_pileups")
    missing_cols <- required_cols[!required_cols %in% names(cohort$inputs)]
    if (length(missing_cols) > 0) {
        stop("Missing required columns in cohort: ", paste(missing_cols, collapse = ", "))
    }

    # Process each sample in parallel
    mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        row <- cohort$inputs[i, ]
        pair_dir <- file.path(output_data_dir, row$pair)

        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }

        out_file <- file.path(pair_dir, "hetsnps.arrow")

        futile.logger::flog.threshold("ERROR")
        tryCatchLog(
            {
                if (!is.null(row$het_pileups) && file.exists(row$het_pileups)) {
                    # Subsample hetsnps
                    hetsnps_gr <- subsample_hetsnps(
                        het_pileups = row$het_pileups,
                        mask = row$hetsnps_mask,
                        sample_size = row$hetsnps_subsample_size,
                        min_normal_freq = row$hetsnps_min_normal_freq,
                        max_normal_freq = row$hetsnps_max_normal_freq
                    )

                    # col2numeric() in granges_to_arrow_scatterplot
                    # expects a hexadecimal color value, not a name
                    hetsnps_gr$col <- grDevices::rgb(
                        t(
                            grDevices::col2rgb(hetsnps_gr$col)
                        ),
                        maxColorValue = 255
                    )

                    # Create arrow table
                    arrow_table <- granges_to_arrow_scatterplot(
                        gr_path = hetsnps_gr,
                        field = row$hetsnps_field,
                        ref = cohort$reference_name,
                        cov.color.field = row$hetsnps_color_field,
                        bin.width = row$hetsnps_bin_width
                    )

                    # Write arrow table
                    message(sprintf("Writing hetsnps arrow file for %s", row$pair))
                    arrow::write_feather(arrow_table, out_file, compression = "uncompressed")
                } else {
                    warning(sprintf("Het pileups file missing for %s", row$pair))
                }
            },
            error = function(e) {
                print(sprintf("Error processing %s: %s", row$pair, e$message))
                NULL
            }
        )
    }, mc.cores = cores, mc.preschedule = TRUE)

    invisible(NULL)
}
