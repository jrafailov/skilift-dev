#' @name initialize_metadata_columns
#' @title Initialize Metadata Columns
#' @description Creates a data.table with all expected columns initialized to NA
#' @param pair The sample pair identifier
#' @return A data.table with all possible columns initialized to NA
initialize_metadata_columns <- function(pair) {
    if (is.null(pair)) {
        stop("pair must be a character string")
    }
    if (!is.character(pair) && is.integer(pair)) {
        warning("pair is an integer, converting to character")
        pair <- as.character(pair)
    }
    if (!is.character(pair)) {
        stop("pair must be a character string")
    }
    if (nchar(pair) == 0) {
        stop("pair cannot be empty")
    }
    
    dt <- data.table(
        # Basic metadata
        pair = pair,
        tumor_type = NA_character_,
        disease = NA_character_,
        primary_site = NA_character_,
        inferred_sex = NA_character_,
        
        # Coverage QC
        coverage_qc = list(NULL),
        
        # Variant counts
        snv_count = NA_integer_,
        snv_count_normal_vaf_greater0 = NA_integer_,
        
        # Coverage parameters
        cov_slope = NA_real_,
        cov_intercept = NA_real_,
        
        # Signature lists
        deconstructsigs_sbs_fraction = list(NULL),
        sigprofiler_indel_fraction = list(NULL),
        sigprofiler_sbs_count = list(NULL),
        sigprofiler_sbs_fraction = list(NULL),
        
        # HRD scores
        hrd_score = NA_real_,
        hrd = list(NULL),
        b1_2 = NA_real_,
        b1 = NA_real_,
        b2 = NA_real_,
        tmb = NA_real_  # Add this line
    )
    return(dt)
}

#' @name add_basic_metadata
#' @title Add Basic Metadata
#' @description
#' Adds basic metadata information such as tumor type, disease, and primary site.
#'
#' @param metadata A data.table containing metadata.
#' @param tumor_type The type of tumor.
#' @param disease The disease associated with the sample.
#' @param primary_site The primary site of the tumor.
#' @return Updated metadata with basic information added.
add_basic_metadata <- function(
    metadata,
    input_tumor_type,
    input_disease,
    input_primary_site
) {
    # Input validation
    if (is.null(metadata)) {
        stop("metadata cannot be NULL")
    }
    if (!is.data.table(metadata)) {
        stop("metadata must be a data.table")
    }
    
    # Validate tumor_type if provided
    if (!is.null(input_tumor_type)) {
        if (!is.character(input_tumor_type)) {
            stop("tumor_type must be NULL or a character string")
        }
        metadata[, tumor_type := input_tumor_type]
    }
    
    # Validate disease if provided
    if (!is.null(input_disease)) {
        if (!is.character(input_disease)) {
            stop("disease must be NULL or a character string")
        }
        metadata[, disease := input_disease]
    }
    
    # Validate primary_site if provided
    if (!is.null(input_primary_site)) {
        if (!is.character(input_primary_site)) {
            stop("primary_site must be NULL or a character string")
        }
        metadata[, primary_site := input_primary_site]
    }
    
    return(metadata)
}

#' @name add_sex_information
#' @title Add Sex Information
#' @description
#' Adds inferred sex information to the metadata based on various inputs.
#'
#' @param metadata A data.table containing metadata.
#' @param input_inferred_sex The inferred sex of the sample.
#' @param jabba_gg A file path to the jabba_gg data.
#' @param tumor_coverage Coverage data for the tumor.
#' @return Updated metadata with inferred sex information added.
add_sex_information <- function(
    metadata,
    input_inferred_sex = NULL,
    jabba_gg = NULL,
    tumor_coverage = NULL 
) {
    # Direct sex specification
    if (!is.null(input_inferred_sex)) {
        if (!input_inferred_sex %in% c("Male", "Female")) {
            warning("inferred sex must be one of `Male` or `Female`")
        }
        metadata[, inferred_sex := input_inferred_sex]
        return(metadata)
    }
    
    # Infer from jabba_gg if available
    if (!is.null(jabba_gg)) {
        gg <- process_jabba(jabba_gg)
        ncn.x <- gg$nodes$dt[
            (seqnames == "X" | seqnames == "chrX"),
            weighted.mean(cn, w = end - start + 1, na.rm = TRUE)
        ]
        metadata[, inferred_sex := ifelse(ncn.x < 1.4, "Male", "Female")]
        return(metadata)
    }
    
    # Infer from tumor coverage if available
    if (!is.null(tumor_coverage)) {
        t_cov <- gGnome:::readCov(tumor_coverage)
        ncn.x <- unique(gr2dt(t_cov)[, foreground.chr := mean(foreground), by = seqnames][, .(seqnames, foreground.chr)])[seqnames %in% c("X", "chrX")]$foreground.chr
        if (length(ncn.x) == 0) {
            warning("Could not extract X chromosome coverage from tumor coverage data")
            return(metadata)
        }
        metadata[, inferred_sex := ifelse(ncn.x < 0.7, "Male", "Female")]
        return(metadata)
    }
    
    return(metadata)
}

#' @name dlrs
#' @title derivative log ratio spread
#' @description
#'
#' function to get the dlrs for coverage data
#' used in meta_data_json
#'
#' @param x foreground from dryclean
#' @return dlrs
#' @export
#' @author Joel Rosiene
dlrs <- function(x) {
    nx <- length(x)
    if (nx < 3) {
        stop("Vector length>2 needed for computation")
    }
    tmp <- embed(x, 2)
    diffs <- tmp[, 2] - tmp[, 1]
    dlrs <- IQR(diffs, na.rm = TRUE) / (sqrt(2) * 1.34)
    return(dlrs)
}

#' @name extract_metrics
#' @title Extract Metrics from QC File
#' @description Helper function to extract specified metrics from a QC file
#' @param qc_data Data.table containing QC metrics
#' @param metrics Named vector of column names to extract
#' @param pair Sample pair identifier
#' @return Data.table with extracted metrics
extract_metrics <- function(qc_data, metrics, pair) {
    # Verify all metrics exist in data
    missing_cols <- setdiff(metrics, names(qc_data))
    if (length(missing_cols) > 0) {
        stop(sprintf("Missing columns in QC data: %s", paste(missing_cols, collapse = ", ")))
    }
    
    # Select and rename columns
    result <- qc_data[, metrics, with = FALSE]
    setnames(result, old = metrics, new = names(metrics))
    
    # Add pair identifier
    result[, pair := pair]
    
    return(result)
}

#' @name process_qc_metrics
#' @title Process Quality Control Metrics
#' @description Processes quality control metrics from multiple QC files and combines them
#' @param estimate_library_complexity Path to the estimate_library_complexity_metrics file
#' @param alignment_summary_metrics Path to the alignment_summary_metrics file
#' @param insert_size_metrics Path to the insert_size_metrics file
#' @param wgs_metrics Path to the wgs_metrics file
#' @param pair Sample pair identifier
#' @return A list containing processed QC metrics
process_qc_metrics <- function(
    estimate_library_complexity,
    alignment_summary_metrics,
    insert_size_metrics,
    tumor_wgs_metrics,
    normal_wgs_metrics,
    pair
) {


    # Define metric mappings for each file type
    complexity_metrics_cols <- c(
        read_pairs_examined = "READ_PAIRS_EXAMINED",
        read_pair_duplicates = "READ_PAIR_DUPLICATES",
        read_pair_optical_duplicates = "READ_PAIR_OPTICAL_DUPLICATES",
        percent_duplication = "PERCENT_DUPLICATION"
    )
    
    alignment_metrics_cols <- c(
        total_reads = "TOTAL_READS",
        pf_reads_aligned = "PF_READS_ALIGNED",
        pf_aligned_bases = "PF_ALIGNED_BASES",
        mean_read_length = "MEAN_READ_LENGTH"
    )
    
    insert_metrics_cols <- c(
        # median_insert_size = "MEDIAN_INSERT_SIZE"
        insert_size = "MEDIAN_INSERT_SIZE"
    )
    
    tumor_wgs_metrics_cols <- c(
        # median_coverage = "MEAN_COVERAGE"
        tumor_median_coverage = "MEDIAN_COVERAGE",
        # pct_30x = "PCT_30X",
        greater_than_or_equal_to_30x = "PCT_30X",
        # pct_50x = "PCT_50X"
        greater_than_or_equal_to_50x = "PCT_50X",
        fraction_excluded = "PCT_EXC_TOTAL"
    )

    normal_wgs_metrics_cols <- c(
        normal_median_coverage = "MEDIAN_COVERAGE",
        fraction_excluded = "PCT_EXC_TOTAL"
    )
    
    test_file_is_present = function(x) {
        (
            !is.null(x)
            && is.character(x)
            && NROW(x) == 1
            && file.exists(x)
        )
    }

    # Read and extract metrics from each file
    complexity_data = data.table(pair = character(0))
    if (test_file_is_present(estimate_library_complexity)) {
        complexity_data <- extract_metrics(
            fread(estimate_library_complexity),
            complexity_metrics_cols,
            pair
        )
    }
    
    alignment_data = data.table(pair = character(0))
    if (test_file_is_present(alignment_summary_metrics)) {
        alignment_data <- extract_metrics(
            fread(alignment_summary_metrics)[CATEGORY=="PAIR"],
            alignment_metrics_cols,
            pair
        )
    }
    
    insert_data = data.table(pair = character(0))
    if (test_file_is_present(insert_size_metrics)) {
        insert_data <- extract_metrics(
            fread(insert_size_metrics)[PAIR_ORIENTATION == "FR"],
            insert_metrics_cols,
            pair
        )
    }

    tumor_wgs_data = data.table(pair = character(0))
    if (test_file_is_present(tumor_wgs_metrics)) {
        tumor_wgs_data <- extract_metrics(
            fread(tumor_wgs_metrics),
            tumor_wgs_metrics_cols,
            pair
        )
        tumor_wgs_data$tumor_median_coverage = round(
            tumor_wgs_data$tumor_median_coverage / (1 - tumor_wgs_data$fraction_excluded)
        )
        tumor_wgs_data$fraction_excluded = NULL
    }

    normal_wgs_data = data.table(pair = character(0))
    if (test_file_is_present(normal_wgs_metrics)) {
        normal_wgs_data <- extract_metrics(
            fread(normal_wgs_metrics),
            normal_wgs_metrics_cols,
            pair
        )
        normal_wgs_data$normal_median_coverage = round(
            normal_wgs_data$normal_median_coverage / (1 - normal_wgs_data$fraction_excluded)
        )
        normal_wgs_data$fraction_excluded = NULL
    }
    
    # c("pair", "total_reads", "pf_reads_aligned", "pf_aligned_bases", 
    # "mean_read_length", "median_insert_size", "median_coverage", 
    # "pct_30x", "pct_50x", "percent_optical_duplication", 
    # "percent_aligned", "percent_optical_dups_of_dups")
    # Merge all metrics on pair
    lst_to_merge = list(
        complexity_data, 
        alignment_data, 
        insert_data, 
        tumor_wgs_data,
        normal_wgs_data
    )
    qc_metrics <- Reduce(function(x, y) {
        data.table::merge.data.table(
            x, y, by = "pair", 
            all.x = TRUE, 
            all.y = TRUE,
            suffixes = c("_x", "_y")
        )
    }   , lst_to_merge)
    
    # Calculate derivative metrics
    # Separating out and writing long way for robustness to
    # missing data.
    qc_metrics$m_reads = qc_metrics$total_reads / 1e6
    qc_metrics$m_reads_mapped = qc_metrics$pf_reads_aligned / 1e6
    qc_metrics$percent_optical_duplication = (
        qc_metrics$read_pair_optical_duplicates /
        qc_metrics$read_pairs_examined
    )
    qc_metrics$percent_aligned = (
        qc_metrics$pf_aligned_bases / 
        ( qc_metrics$total_reads * qc_metrics$mean_read_length )
    )
    qc_metrics$percent_optical_dups_of_dups = (
        qc_metrics$read_pair_optical_duplicates / 
        qc_metrics$read_pair_duplicates
    )
    # # FIXME, need to account for tumor and normal
    # qc_metrics$normal_median_coverage = NA_integer_

    # Remove any metrics that are NA
    # This can happen if any derivative metrics
    # are calculated from qc inputs not provided
    # in above lines.
    for (colnm in names(qc_metrics)) {
        is_all_na = all(is.na(qc_metrics[[colnm]]))
        if (is_all_na) qc_metrics[[colnm]] = NULL
    }

    # qc_metrics[, `:=`(
    #     m_reads = total_reads / 1e6,
    #     m_reads_mapped = pf_reads_aligned / 1e6,
    #     percent_optical_duplication = read_pair_optical_duplicates / read_pairs_examined,
    #     percent_aligned = pf_aligned_bases / (total_reads * mean_read_length),
    #     percent_optical_dups_of_dups = read_pair_optical_duplicates / read_pair_duplicates
    # )]
    
    return(as.list(qc_metrics))
}

#' @name add_coverage_metrics
#' @title Add Coverage Metrics
#' @description
#' Adds coverage metrics to the metadata based on tumor coverage and QC metrics.
#'
#' @param metadata A data.table containing metadata.
#' @param tumor_coverage Coverage data for the tumor.
#' @param foreground_col_name Name of the column in the coverage data to use for foreground coverage
#' @param estimate_library_complexity Path to library complexity metrics file
#' @param alignment_summary_metrics Path to alignment summary metrics file
#' @param insert_size_metrics Path to insert size metrics file
#' @param wgs_metrics Path to WGS metrics file
#' @return Updated metadata with coverage metrics added.
add_coverage_metrics <- function(
    metadata,
    tumor_coverage = NULL,
    foreground_col_name = "foreground",
    estimate_library_complexity = NULL,
    alignment_summary_metrics = NULL,
    insert_size_metrics = NULL,
    tumor_wgs_metrics = NULL,
    normal_wgs_metrics = NULL
) {
    coverage_variance <- NULL
    if (!is.null(tumor_coverage)) {
        foreground = as.data.table(readRDS(tumor_coverage))[[foreground_col_name]]
        cov_var <- dlrs(foreground) /100
        coverage_variance <- list(coverage_variance = cov_var)
        metadata$coverage_qc <- list(list(coverage_variance))
    }
    
    processed_metrics <- NULL
    if (
        !is.null(estimate_library_complexity) ||
        !is.null(alignment_summary_metrics) ||
        !is.null(insert_size_metrics) ||
        !is.null(tumor_wgs_metrics) ||
        !is.null(normal_wgs_metrics)
        ) {
        
        processed_metrics <- process_qc_metrics(
            estimate_library_complexity,
            alignment_summary_metrics,
            insert_size_metrics,
            tumor_wgs_metrics,
            normal_wgs_metrics,
            metadata$pair
        )
        
        metadata$coverage_qc <- list(as.list(c(processed_metrics, coverage_variance)))
        ## FIXME: hardcoding this for now.
        ## median and normal coverages are at top level..
        metadata$tumor_median_coverage = processed_metrics$tumor_median_coverage
        metadata$normal_median_coverage = processed_metrics$normal_median_coverage
    }

    return(metadata)
}

#' @name vcf_count
#' @title vcf_count
#' @description
#' takes in a vcf and returns a data.table with total count and count of variants with normal vaf greater than 0
#'
#' @param vcf_path Path to the VCF file
#' @param genome Reference genome name (e.g., "hg19", "hg38")
#' @return data.table
#' @export
#' @author Shihab Dider, Tanubrata Dey
vcf_count <- function(
    vcf_path,
    genome) {
    if (!file.exists(vcf_path)) {
        stop("VCF file does not exist.")
    }
    is_len_one = NROW(vcf_path) == 1
    is_character = is.character(vcf_path)
    is_na = is_len_one && is_character && (is.na(vcf_path) || vcf_path %in% c("NA", "na"))
    is_vcf = is_len_one && is_character && !is_na && grepl("vcf(.gz|.bgz)?$", vcf_path)
    is_bcf = is_len_one && is_character && !is_na && grepl("bcf(.gz|.bgz)?$", vcf_path)
    is_valid_path = is_vcf || is_bcf
    if (!is_valid_path) stop("either bcf or vcf path must be provided to vcf_count()")
    if (is_bcf) {
        temp_vcf = tempfile(fileext = ".vcf")
        on.exit(unlink(temp_vcf, force = TRUE))
        system2("bcftools", c("view", vcf_path, "-Ov", "-o", temp_vcf))
        vcf_path = temp_vcf
    }
    vcf <- readVcf(vcf_path, genome)

    # Filter for PASS variants
    pass_variants <- rowRanges(vcf)$FILTER == "PASS"
    vcf <- vcf[pass_variants, ]

    snv_count <- length(vcf)

    return(
        data.table(
            category = c("snv_count", "snv_count_normal_vaf_greater0"),
            counts = c(snv_count, NA)
        )
    )
}
    
#' @name add_variant_counts
#' @title Add Variant Counts
#' @description
#' Adds variant counts to the metadata based on somatic and germline SNVs.
#'
#' @param metadata A data.table containing metadata.
#' @param somatic_snvs Path to somatic SNV VCF file
#' @param genome The genome reference used.
#' @return Updated metadata with variant counts added.
add_variant_counts <- function(
    metadata,
    somatic_snvs = NULL,
    genome = "hg19"
) {
    if (!is.null(somatic_snvs)) {
        is_path_character = is.character(somatic_snvs)
        is_length_one = NROW(somatic_snvs) == 1
        is_snvs_exists = is_path_character && is_length_one && file.exists(somatic_snvs)
        is_rds = is_snvs_exists && grepl("rds$", somatic_snvs)
        is_vcf = is_snvs_exists && grepl("(v|b)cf(.gz|.bgz)?$", somatic_snvs)
        is_txt = is_snvs_exists && grepl("(txt|maf|(t|c)sv)(.gz|.bgz)?$", somatic_snvs)
        is_other = is_txt || is_rds
        snv_counts_dt = data.table()
        if (is_vcf) {
            snv_counts_dt <- vcf_count(somatic_snvs, genome = genome)
        }
        if (is_rds) snvs = readRDS(somatic_snvs)
        if (is_txt) snvs = fread(somatic_snvs)
        if (is_other) {
            snv_counts_dt = data.table(
                category = c("snv_count", "snv_count_normal_vaf_greater0"),
                counts = c(NROW(snvs), NA_integer_)
            )
        }
        
        metadata$snv_count <- snv_counts_dt[category == "snv_count", ]$counts
        ## count up filtered variants in VCF 

        metadata$snv_count_normal_vaf_greater0 <- snv_counts_dt[category == "snv_count_normal_vaf_greater0", ]$counts
    }
    return(metadata)
}

#' @name add_sv_counts
#' @title Add Structural Variant Counts
#' @description Adds junction, loose end, and total SV counts from JaBbA graph
#' @param meta_dt A data.table containing metadata
#' @param jabba_gg Path to JaBbA graph RDS file
#' @return Updated metadata with SV counts
add_sv_counts <- function(metadata, jabba_gg = NULL) {
    if (is.null(jabba_gg)) {
        return(metadata)
    }
    
    gg <- process_jabba(jabba_gg)
    
    # Count junctions
    if (nrow(gg$junctions$dt) > 0) {
        metadata$junction_count <- nrow(gg$junctions$dt[type != "REF", ])
    } else {
        warning("No junctions found in JaBbA graph, skipping SV counts...")
        metadata$junction_count <- 0
    }
    
    # Count loose ends
    if (length(gg$loose) > 0) {
        metadata$loose_count <- nrow(as.data.table(gg$loose)[terminal == FALSE, ])
    } else {
        warning("No loose ends found in JaBbA graph, skipping loose counts...")
        metadata$loose_count <- 0
    }
    
    # Calculate total SV count
    metadata[, sv_count := (junction_count)]
    
    return(metadata)
}

#' @name add_purity_ploidy
#' @title Add Purity and Ploidy Metrics
#' @description Adds purity and ploidy information from JaBbA graph
#' @param metadata A data.table containing metadata
#' @param jabba_gg Path to JaBbA graph RDS file
#' @return Updated metadata with purity and ploidy information
add_purity_ploidy <- function(metadata, purple_pp_bestFit = NULL, jabba_gg = NULL, tumor_coverage = NULL) {
    is_null_or_na_jabba_gg = is.null(jabba_gg) || any(is.na(jabba_gg))
    is_null_or_na_purple = is.null(purple_pp_bestFit) || any(is.na(purple_pp_bestFit))
    is_jabba_gg_and_purple_absent = is_null_or_na_jabba_gg && is_null_or_na_purple

    if (is_jabba_gg_and_purple_absent) {
        return(metadata)
    }
    
    if (!is_null_or_na_jabba_gg) {
        gg <- process_jabba(jabba_gg)

        purity = base::get("purity", gg$meta) # Errors out if not found
        ploidy = base::get("ploidy", gg$meta) # Errors out if not found
    } else if (!is_null_or_na_purple) {
        purple_best_fit = fread(purple_pp_bestFit)
        purity = purple_best_fit$purity
        ploidy = purple_best_fit$ploidy
    }

    metadata$purity <- purity
    metadata$ploidy <- ploidy
    metadata$beta = purity / (purity * ploidy + 2*(1 - purity)) # from Multiplicity
    metadata$gamma = 2*(1 - purity) / (purity * ploidy + 2*(1 - purity)) # from Multiplicity
    
    return(metadata)
}

#' @name add_loh
#' @title Add Loss of Heterozygosity Metrics
#' @description Adds LOH metrics from JaBbA graph
#' @param metadata A data.table containing metadata
#' @param jabba_gg Path to JaBbA graph RDS file
#' @param seqnames_loh Vector of sequence names to consider for LOH
#' @return Updated metadata with LOH metrics
add_loh <- function(
    metadata,
    jabba_gg = NULL,
    seqnames_loh = c(1:22)
) {
    if (is.null(jabba_gg) || is.null(seqnames_loh)) {
        return(metadata)
    }
    
    gg <- process_jabba(jabba_gg)
    nodes.dt <- gg$nodes$dt
    nodes.dt[, seqnames := gsub("chr", "", seqnames)]
    nodes.dt <- nodes.dt[seqnames %in% seqnames_loh]
    totalseglen <- nodes.dt$width %>% sum()
    
    if ("cn.low" %in% names(nodes.dt)) {
        lohsegs <- nodes.dt[cn.low == 0, ] %>%
            .[cn.high > 0] %>%
            .$width %>%
            sum()
        loh_frc <- lohsegs / totalseglen
        metadata[, `:=`(
            loh_fraction = loh_frc,
            loh_seglen = lohsegs,
            loh_total_genome = totalseglen
        )]
    } else {
        metadata$loh_fraction <- "Not Allelic Jabba"
    }
    
    return(metadata)
}

#' @name add_fga
#' @title Add Fraction Genome Altered (FGA)
#' @param metadata A data.table containing metadata
#' @param jabba_gg Path to JaBbA graph RDS file
#' @param seqnames Vector of sequence names to consider
#' @return Updated metadata with LOH metrics
add_fga <- function(
    metadata,
    jabba_gg = NULL,
    seqnames = c(1:22)
) {
    if (is.null(jabba_gg) || is.null(seqnames)) {
        return(metadata)
    }
    
    gg <- process_jabba(jabba_gg)
    nodes.dt <- gg$nodes$dt
    nodes.dt[, seqnames := gsub("chr", "", seqnames)]
    nodes.dt <- base::subset(nodes.dt, subset = nodes.dt$seqnames %in% seqnames)
    totalseglen <- nodes.dt$width %>% sum()

    
    fga_width = sum(nodes.dt[cn != 2]$width)
    fga <- fga_width / totalseglen
    metadata$loh_fraction = fga ## FIXME: This is FGA, this is to maintain compatibility with gOS frontend
    
    return(metadata)
}

#' @name add_genome_length
#' @title Add Total Genome Length
#' @description Adds total genome length from JaBbA graph
#' @param metadata A data.table containing metadata
#' @param jabba_gg Path to JaBbA graph RDS file
#' @param seqnames_genome_width_or_genome_length Vector of sequence names to include in genome length or genome length as a numeric
#' @return Updated metadata with genome length information
add_genome_length <- function(
    metadata,
    jabba_gg = NULL,
    seqnames_genome_width_or_genome_length = c(1:22, "X", "Y")
) {

    is_null_jabba = is.null(jabba_gg)

    is_genome_length_already_provided = is.numeric(seqnames_genome_width_or_genome_length) && NROW(seqnames_genome_width_or_genome_length) == 1

    # handle targeted panels/whole exome
    if (is_genome_length_already_provided) {
        metadata$total_genome_length <- seqnames_genome_width_or_genome_length
        return(metadata)
    }

    # otherwise, handle genome-wide
    is_no_jabba_or_genome_length_provided = is_null_jabba && is.null(seqnames_genome_width_or_genome_length)
    if (is_no_jabba_or_genome_length_provided) {
        return(metadata)
    }

    if (is_null_jabba) {
         ## FIXME: hardcoding to hg19!!
        genome_length_hg19 = 3095677412
        message("NO JABBA OR GENOME LENGTH PROVIDED, hardcoding to hg19: ", genome_length_hg19)
        metadata$total_genome_length = genome_length_hg19
        return(metadata)
    }

    gg <- process_jabba(jabba_gg)
    nodes.gr <- gg$nodes$gr
    seqlengths.dt <- suppressWarnings(
        as.data.table(
            seqinfo(nodes.gr),
            keep.rownames = "seqnames"
        )
    )
    seqlengths.dt[, seqnames := gsub("chr", "", seqnames)]
    if (!is.null(seqnames_genome_width_or_genome_length)) {
        seqlengths.dt <- seqlengths.dt[seqnames %in% seqnames_genome_width_or_genome_length, ]
    }
    
    metadata$total_genome_length <- sum(seqlengths.dt$seqlengths)
    
    return(metadata)
}

#' @name add_sv_types
#' @title Add Structural Variant Types
#' @description Adds counts of different SV types from JaBbA graph and complex events
#' @param metadata A data.table containing metadata
#' @param jabba_gg Path to JaBbA graph RDS file
#' @param complex Path to complex events RDS file
#' @return Updated metadata with SV type counts
add_sv_types <- function(metadata, jabba_gg = NULL, complex = NULL) {
    sv_like_types_count <- NULL
    if (!is.null(jabba_gg)) {
        gg <- process_jabba(jabba_gg)
        if (!is.null(gg$edges$dt) && nrow(gg$edges$dt) > 0) {
            sv_like_types_count <- table(gg$edges$dt[type == "ALT" & class != "REF"]$class)
            metadata$sv_types_count <- list(list(list(sv_like_types_count)))
        } else {
            warning("No edges found in JaBbA graph, skipping SV-like types count...")
        }
    }
    
    qrp_counts <- data.table(qrpmin = 0, qrpmix = 0, qrppos = 0)
    
    sv_types_count <- NULL
    if (!is.null(complex)) {
        complex_data <- process_jabba(complex)
        if (!is.null(complex_data$meta$events$type)) {
            sv_types_count <- table(complex_data$meta$events$type)
            
            # Update QRP counts
            qrp_types <- c("qrpmin", "qrpmix", "qrppos")
            for (qrp_type in qrp_types) {
                if (qrp_type %in% names(sv_types_count)) {
                    qrp_counts[[qrp_type]] <- sv_types_count[[qrp_type]]
                }
            }
            
        }
    } else {
        warning("Complex events and JaBbA graph not found as inputs, skipping SV types count...")
    }
    
    metadata$sv_types_count <- list(as.list(c(sv_types_count, sv_like_types_count)))
    metadata$hrd <- list(as.list(qrp_counts))
    return(metadata)
}

#' @name add_coverage_parameters
#' @title Add Coverage Parameters
#' @description Adds coverage parameters to the metadata based on tumor coverage
#' @param metadata A data.table containing metadata
#' @param tumor_coverage Coverage data for the tumor
#' @param field The field to use for coverage calculation
#' @return Updated metadata with coverage parameters added
add_coverage_parameters <- function(metadata, tumor_coverage, field = "foreground") {
    if (!is.null(tumor_coverage)) {
        if (is.null(metadata$purity) || is.null(metadata$ploidy)) {
            warning("Purity and ploidy not found in metadata, cov_slope and cov_intercept will not be calculated")
        }
        cov <- tumor_coverage %>% readRDS()
		coverage_values = base::get0(field, as.environment(as.list(mcols(cov))), ifnotfound = NULL)
		lst_cov_bool = Skilift::test_coverage_normalized(coverage_values)
		is_cov_likely_normalized = lst_cov_bool$is_cov_likely_normalized
		is_cov_near_one = lst_cov_bool$is_cov_near_one
		if (!is_cov_likely_normalized && !is_cov_near_one) {
			message("Assuming coverage is in read coverage per bin and paired end, 151 bp reads, rescaling")
			mcols(cov)[[field]] = coverage_values * PAIRED_READS_FACTOR * READ_LENGTH / width(cov)
		} else {
			message("Assuming coverage is mean-normalized, ignoring rescaling to base coverage")
		}
        # mcols(cov)[[field]] <- mcols(cov)[, field] * 2 * 151 / width(cov)
        rel2abs.cov <- skitools::rel2abs(cov,
            field = field,
            purity = metadata$purity,
            ploidy = metadata$ploidy,
            return.params = TRUE
        )
        metadata$cov_slope <- rel2abs.cov[1] %>% unname()
        metadata$cov_intercept <- rel2abs.cov[2] %>% unname()
    }
    return(metadata)
}

#' @name add_het_pileups_parameters
#' @title Add Heterozygous Pileups Parameters
#' @description Adds slope and intercept parameters calculated from heterozygous pileups
#' @param metadata A data.table containing metadata
#' @param het_pileups_wgs Path to heterozygous pileups WGS data
#' @return Updated metadata with heterozygous pileups parameters
add_het_pileups_parameters <- function(metadata, het_pileups) {
    if (!is.null(het_pileups)) {
        hets.read <- grab.hets(het_pileups) %>% gr2dt()
        hets.read <- dt2gr(hets.read[, .(count = sum(count)), by = c("seqnames", "start", "end")])
        if (is.null(metadata$purity) || is.null(metadata$ploidy)) {
            warning("Purity and ploidy not found in metadata, het_slope and het_intercept will not be calculated")
        }
        rel2abs.hets <- skitools::rel2abs(hets.read,
            field = "count",
            purity = metadata$purity,
            ploidy = metadata$ploidy,
            return.params = TRUE
        )
        metadata$hets_slope <- rel2abs.hets[1] %>% unname()
        metadata$hets_intercept <- rel2abs.hets[2] %>% unname()
    }
    return(metadata)
}

#' @name add_tmb
#' @title Add Tumor Mutation Burden
#' @description
#' Calculates and adds tumor mutation burden (TMB) to the metadata based on SNV count and genome length
#'
#' @param metadata A data.table containing metadata
#' @param somatic_snvs Path to somatic SNV VCF file
#' @param jabba_gg Path to JaBbA graph RDS file
#' @param genome The genome reference used
#' @param seqnames_genome_width_or_genome_length Sequence names for genome width calculation
#' @return Updated metadata with TMB value added
add_tmb <- function(
    metadata,
    somatic_snvs = NULL,
    jabba_gg = NULL,
    genome = "hg19",
    seqnames_genome_width_or_genome_length = c(1:22, "X", "Y")
) {
    meta_dt <- copy(metadata)
    # First ensure we have snv_count
    if (is.null(meta_dt$snv_count) || is.na(meta_dt$snv_count)) {
        meta_dt <- add_variant_counts(meta_dt, somatic_snvs, genome)
    }
    
    # Then ensure we have total_genome_length
    if (is.null(meta_dt$total_genome_length) || is.na(meta_dt$total_genome_length)) {
        meta_dt <- add_genome_length(meta_dt, jabba_gg, seqnames_genome_width_or_genome_length)
    }
    
    
    # Calculate TMB if we have both required values
    is_tmb_computable <- !is.na(meta_dt$snv_count) && !is.na(meta_dt$total_genome_length) && !is.null(meta_dt$snv_count) && !is.null(meta_dt$total_genome_length)
    if (is_tmb_computable) {
        tmb_vals = meta_dt$snv_count / (as.numeric(meta_dt$total_genome_length) / 1e6)
        metadata[, tmb := tmb_vals]
        metadata[, tmb := round(tmb, digits = 3)]
    } else {
        warning("Cannot calculate TMB without both snv_count and total_genome_length")
    }
    
    return(metadata)
}

#' @name compute_signature_averages
#' @title Compute Signature Averages
#' @description
#' Computes signature averages from various signature files
#'
#' @param sig_file Path to signature file
#' @param sample Sample or pair name needed to subset to the correct sample
#' @param is_indel Boolean indicating if processing indel signatures
#' @param is_deconstruct_sigs Boolean indicating if using deconstructSigs format
#' @return List of signature averages
compute_signature_averages <- function(
    sig_file,
    is_indel = FALSE,
    is_deconstruct_sigs = FALSE
) {
    if (is_deconstruct_sigs) {
        sig_data <- readRDS(sig_file)
        weights <- as.data.table(sig_data$weights)
        weights <- t(weights)
        sigs.dt <- as.data.table(weights, keep.rownames = "Signature") %>% 
            setnames(., c("Signature", "weights"))
        sigs.dt[, Signature := gsub("Signature.", "SBS", Signature)]
        sigs.vect <- sigs.dt$weights
        names(sigs.vect) <- sigs.dt$Signature
        return(sigs.vect)
    } else {
        sig.dt <- fread(sig_file)
        if (nrow(sig.dt) == 0) {
            stop("No signatures found in file")
        }
        pair <- sig.dt$Samples[1]
        sig.dt[, pair := pair]
        sig.dt[, Samples := NULL]
        sig.dt_avg <- copy(sig.dt)
        sig.dt_avg[, pair := NULL]
        row_sum <- sig.dt_avg %>% rowSums()
        sigs.dt <- melt.data.table(sig.dt_avg, measure.vars = names(sig.dt_avg)) %>% 
            setnames(., c("signature", "value"))
        sigs.dt[, avg_value := value / row_sum]
        sigs.vect <- sigs.dt$avg_value
        names(sigs.vect) <- sigs.dt$signature
        sigs_counts.vect <- sigs.dt$value
        names(sigs_counts.vect) <- sigs.dt$signature
        
        if (is_indel) {
            return(list(
                indel_fraction = sigs.vect,
                indel_count = sigs_counts.vect
            ))
        } else {
            return(list(
                sbs_fraction = sigs.vect,
                sbs_count = sigs_counts.vect
            ))
        }
    }
}

#' @name compute_cosine_similarity
#' @title Compute Cosine Similarities   
#' @description compute consine similarities between reference and attributed signatures
#'
#' @param metadata A data.table containing metadata (must have sigprofiler_sbs_count or sigprofiler_indel_count)
#' @param probabilities Path to probabilities file
#' @param matrix_file Path to matrix file
#' @param is_indel Boolean indicating if processing indel signatures
#' @param is_deconstruct_sigs Boolean indicating if using deconstructSigs format
#' @param reference Reference genome name (e.g., "hg19", "hg38")
#' @return List of signature averages
compute_cosine_similarity <- function(
    metadata,
    probabilities,
    matrix_file,
    is_indel = FALSE,
    is_deconstruct_sigs = FALSE,
    reference = "hg19"
) {

    if (reference == "hg19") {
        if (is_indel) {
            sig_file <- system.file("extdata", "/data/COSMIC_v3.4_ID_GRCh37.txt", package = "Skilift")
        } else {
            sig_file <- system.file("extdata", "/data/COSMIC_v3.4_SBS_GRCh37.txt", package = "Skilift")
        }
    } else {
        if (is_indel) {
            sig_file <- system.file("extdata", "/data/COSMIC_v3.4_ID_GRCh38.txt", package = "Skilift")
        } else {
            sig_file <- system.file("extdata", "/data/COSMIC_v3.4_SBS_GRCh38.txt", package = "Skilift")
        }
    }
    
    weight_reads <- fread(sig_file, fill = TRUE)
    expected_reads <- copy(weight_reads)

    row_names <- expected_reads[[1]]
    expected_reads <- expected_reads[, 2:ncol(expected_reads)]
    rownames(expected_reads) <- row_names
    
    ### EXPECTED SIGNATURES ###
    activities <- if(is_indel) {
        metadata$sigprofiler_indel_count %>% unlist()
    } else {
        metadata$sigprofiler_sbs_count %>% unlist()
    }
    matching_signatures <- intersect(names(activities), colnames(expected_reads))
    
    if(length(matching_signatures) == 0) {
        warning("No matching signature names found between expected_reads and activities")
        return(NULL)
    }
    
    weighted_expected <- as.matrix(expected_reads)
    
    weighted_expected[, matching_signatures] <- sweep(weighted_expected[, matching_signatures], 2, activities[matching_signatures], "*")

    ### ATTRIBUTED SIGNATURES ###     
    probabilities = fread(probabilities, fill = TRUE)
    prob.matrix <- as.matrix(probabilities[, 3:ncol(probabilities)])

    matrix_file <- fread(matrix_file, fill = TRUE)
    if (ncol(matrix_file) < 2) {
        stop("matrix_file must have at least two columns")
    }
    multiplier <- as.numeric(matrix_file[[2]])
    if (length(multiplier) != nrow(prob.matrix)) {
        stop("Length of multiplier vector must match the number of rows in prob.matrix")
    }
    
    ### COMPUTING COSINE SIMILARITY ###
    weighted_prob_matrix <- sweep(prob.matrix, 1, multiplier, "*")

    cosine_similarities <- numeric(ncol(weighted_prob_matrix))
    names(cosine_similarities) <- colnames(weighted_prob_matrix)
    
    # Ensure matching columns only
    common_sigs <- intersect(colnames(weighted_prob_matrix), colnames(weighted_expected))
    
    if (length(common_sigs) == 0) {
        warning("No common signatures between matrices")
        return(NULL)
    }
    
    # Subset to common columns
    prob_subset <- weighted_prob_matrix[, common_sigs, drop = FALSE]
    exp_subset <- weighted_expected[, common_sigs, drop = FALSE]
    
    # Vectorized cosine similarity calculation
    dot_products <- colSums(prob_subset * exp_subset)
    norms_prob <- sqrt(colSums(prob_subset^2))
    norms_exp <- sqrt(colSums(exp_subset^2))
    
    # Handle zero norm cases
    zero_norms <- (norms_prob == 0) | (norms_exp == 0)
    cosine_similarities <- ifelse(zero_norms, 0, dot_products / (norms_prob * norms_exp))
    names(cosine_similarities) <- common_sigs

    # Replace NaN values with 0
    cosine_similarities[is.nan(cosine_similarities)] <- 0

    return(list(cosine_similarities = cosine_similarities, 
                weighted_expected = weighted_expected, 
                weighted_prob_matrix = weighted_prob_matrix))

}

#' @name add_signatures
#' @title Add Signatures
#' @description
#' Adds mutational signatures to the metadata based on various signature matrices.
#'
#' @param metadata A data.table containing metadata.
#' @param signatures_pair_name The name of the signature pair.
#' @param activities_sbs_signatures Activities of SBS signatures.
#' @param decomposed_sbs_signatures Probabilities of SBS signatures.
#' @param matrix_sbs_signatures Matrix of SBS signatures.
#' @param activities_indel_signatures Activities of indel signatures.
#' @param decomposed_indel_signatures Matrix of indel signatures.
#' @param matrix_indel_signatures Matrix of indel signatures.
#' @param deconstructsigs_sbs_signatures DeconstructSigs SBS signatures.
#' @return Updated metadata with mutational signatures added.
add_signatures <- function(
    metadata,
    activities_sbs_signatures,
    decomposed_sbs_signatures,
    matrix_sbs_signatures,
    activities_indel_signatures,
    decomposed_indel_signatures,
    matrix_indel_signatures,
    deconstructsigs_sbs_signatures
) {

    # Initialize list columns with proper length
    n <- nrow(metadata)
    metadata[, `:=`(
        deconstructsigs_sbs_fraction = list(list()),
        deletionInsertion = list(list()),
        sigprofiler_indel_fraction = list(list()),
        sigprofiler_indel_count = list(list()),
        sigprofiler_sbs_fraction = list(list()),
        sigprofiler_sbs_count = list(list()),
        signatures = list(list())
    )]

    if (!is.null(activities_sbs_signatures)) {
        
        signatures <- compute_signature_averages(
            sig_file = activities_sbs_signatures,
            is_indel = FALSE,
            is_deconstruct_sigs = FALSE
        )

        metadata$sigprofiler_sbs_fraction <- list(as.list(signatures[["sbs_fraction"]]))
        metadata$sigprofiler_sbs_count <- list(as.list(signatures[["sbs_count"]]))
        metadata$signatures <- list(as.list(signatures[["sbs_fraction"]]))
        
        if (!is.null(decomposed_sbs_signatures) && !is.null(matrix_sbs_signatures)) {
            tryCatch({

                cosine_similarities <- compute_cosine_similarity(
                metadata = metadata,
                probabilities = decomposed_sbs_signatures,
                matrix_file = matrix_sbs_signatures,
                is_indel = FALSE,
                is_deconstruct_sigs = FALSE #' TODO ADD FUNCTIONALITY FOR different references
                )

                metadata$sigprofiler_sbs_cosine_similarity <- list(as.list(cosine_similarities[["cosine_similarities"]]))
            }, error = function(e) {
                warning(sprintf("Error computing SBS cosine similarities: %s", e$message))
            })
        }
    }

    if (!is.null(activities_indel_signatures)) {
        deletionInsertion <- compute_signature_averages(
            sig_file = activities_indel_signatures,
            is_indel = TRUE,
            is_deconstruct_sigs = FALSE
        )
        metadata$deletionInsertion <- list(as.list(deletionInsertion[["indel_fraction"]]))
        metadata$sigprofiler_indel_fraction <- list(as.list(deletionInsertion[["indel_fraction"]]))
        metadata$sigprofiler_indel_count <- list(as.list(deletionInsertion[["indel_count"]]))

        if(!is.null(decomposed_indel_signatures) && !is.null(matrix_indel_signatures)) {
            

            tryCatch({
                cosine_similarities <- compute_cosine_similarity(
                    metadata = metadata,
                    probabilities = decomposed_indel_signatures,
                    matrix_file = matrix_indel_signatures,
                    is_indel = TRUE,
                    is_deconstruct_sigs = FALSE #' TODO ADD FUNCTIONALITY FOR different references
                )

                metadata$sigprofiler_indel_cosine_similarity <- list(as.list(cosine_similarities[["cosine_similarities"]]))
            }, error = function(e) {
                warning(sprintf("Error computing indel cosine similarities: %s", e$message))
            })
        }
    }

    if (!is.null(deconstructsigs_sbs_signatures)) {
        signatures <- compute_signature_averages(
            sig_file = deconstructsigs_sbs_signatures,
            is_indel = FALSE,
            is_deconstruct_sigs = TRUE
        )
        metadata$deconstructsigs_sbs_fraction <- list(as.list(signatures))
    }

    return(metadata)
}

#' @name add_hrd_scores
#' @title Add HRD Scores
#' @description
#' Adds homologous recombination deficiency (HRD) scores to the metadata.
#'
#' @param metadata A data.table containing metadata.
#' @param hrdetect HRDetect scores.
#' @param onenesstwoness Oneness and twoness scores.
#' @return Updated metadata with HRD scores added.
add_hrd_scores <- function(metadata, hrdetect, onenesstwoness) {
    if (!is.null(onenesstwoness)) {
        onetwo <- readRDS(onenesstwoness)
        hrd_values <- data.table(
            b1_2_score = onetwo$ot_scores[, "BRCA1"] + onetwo$ot_scores[, "BRCA2"],
            b1_score = onetwo$ot_scores[, "BRCA1"],
            b2_score = onetwo$ot_scores[, "BRCA2"],
            wt_score = onetwo$ot_scores[, "OTHER"],
            DUP_1kb_100kb = onetwo$expl_variables[, "DUP_1kb_100kb"],
            SNV3 = onetwo$expl_variables[, "SNV3"],
            SNV8 = onetwo$expl_variables[, "SNV8"],
            RS3 = onetwo$expl_variables[, "RS3"],
            RS5 = onetwo$expl_variables[, "RS5"],
            del_mh_prop = onetwo$expl_variables[, "del.mh.prop"],
            ihdel = onetwo$expl_variables[, "ihdel"],
            loh_score = onetwo$expl_variables[, "hrd"],
            qrppos = onetwo$expl_variables[, "qrppos"],
            qrpmin = onetwo$expl_variables[, "qrpmin"],
            qrpmix = onetwo$expl_variables[, "qrpmix"],
            qrdup = onetwo$expl_variables[, "qrdup"],
            qrdel = onetwo$expl_variables[, "qrdel"],
            tib = onetwo$expl_variables[, "tib"]
        )

        if (!is.null(hrdetect)) {
            hrd <- readRDS(hrdetect)
            hrd_values <- hrd_values[, hrd_score := hrd$hrdetect_output[1, "Probability"]]
        } else {
            warning("HRDetect scores not found, skipping HRD scores...")
        }

        metadata$hrd <- list(as.list(hrd_values))

    } else {
        warning("Oneness and twoness scores not found, skipping...")
    }
    
    return(metadata)
}

#' @name add_msisensor_score
#' @title Add MSIsensor Score
#' @description
#' Adds MSIsensor score to the metadata.    
#' 
#' @param metadata A data.table containing metadata.
#' @param msisensor_pro Path to MSIsensor profile file.
#' @return Updated metadata with MSIsensor score added.
add_msisensor_score <- function(metadata, msisensorpro) {
    tryCatch({
        if (!is.null(msisensorpro)) {
            msisensor_data <- fread(msisensorpro)
            if (nrow(msisensor_data) != 0) {
                score <- msisensor_data[[3]][1]
                label.msi <- ifelse(score < 10, "MSS",
                                    ifelse(score < 20, "MSI-Low", "MSI-High"))
                dt <- data.table(
                    score = score / 100,
                    n_unstable = msisensor_data[[2]][1],
                    n_evaluated = msisensor_data[[1]][1],
                    label = label.msi
                )
                metadata$msisensor <- list(as.list(dt))
            } else {
                warning("MSIsensor profile is empty, skipping MSIsensor score...")
            }
        } else {
            warning("MSIsensor profile not provided, skipping MSIsensor score...")
        }
    }, error = function(e) {
        warning(sprintf("Error processing MSIsensor score: %s", e$message))
    })
    return(metadata)
}

#' @title Create Metadata for a Sample
#' @description
#' Creates a comprehensive metadata object for a single sample pair by aggregating various data inputs.
#'
#' @param pair The sample pair identifier.
#' @param tumor_type The type of tumor.
#' @param disease The disease associated with the sample.
#' @param primary_site The primary site of the tumor.
#' @param inferred_sex The inferred sex of the sample.
#' @param jabba_gg A file path to the jabba_gg data.
#' @param events Structural variant events.
#' @param somatic_snvs Somatic single nucleotide variants.
#' @param germline_snvs Germline single nucleotide variants.
#' @param tumor_coverage Coverage data for the tumor.
#' @param estimate_library_complexity Path to the estimate_library_complexity_metrics file.
#' @param alignment_summary_metrics Path to the alignment_summary_metrics file.
#' @param insert_size_metrics Path to the insert_size_metrics file.
#' @param wgs_metrics Path to the wgs_metrics file.
#' @param het_pileups Heterozygous pileups data.
#' @param signatures_pair_name The name of the signature pair.
#' @param matrix_indel_signatures Matrix of indel signatures.
#' @param matrix_sbs_signatures Matrix of SBS signatures.
#' @param decomposed_sbs_signatures Decomposed SBS signatures.
#' @param decomposed_indel_signatures Decomposed indel signatures.
#' @param deconstructsigs_sbs_signatures DeconstructSigs SBS signatures.
#' @param activities_indel_signatures Activities of indel signatures.
#' @param activities_sbs_signatures Activities of SBS signatures.
#' @param hrdetect HRDetect scores.
#' @param onenesstwoness Oneness and twoness scores.
#' @param msisensorpro MSIsensor profile file.
#' @param genome The genome reference used.
#' @param seqnames_loh Sequence names for loss of heterozygosity.
#' @param seqnames_genome_width_or_genome_length Sequence names and genome width in list or genome length as a numeric
#' @return A data.table containing the metadata for a single sample pair.
#' @export
create_metadata <- function(
    pair,
    tumor_type = NULL,
    disease = NULL, 
    primary_site = NULL,
    inferred_sex = NULL,
    purple_pp_bestFit = NULL,
    jabba_gg = NULL,
    events = NULL,
    somatic_snvs = NULL,
    germline_snvs = NULL,
    tumor_coverage = NULL,
    foreground_col_name = "foreground",
    estimate_library_complexity = NULL,
    alignment_summary_metrics = NULL,
    insert_size_metrics = NULL,
    tumor_wgs_metrics = NULL,
    normal_wgs_metrics = NULL,
    het_pileups = NULL,
    activities_indel_signatures = NULL,
    deconstructsigs_sbs_signatures = NULL,
    activities_sbs_signatures = NULL,
    decomposed_sbs_signatures = NULL,
    decomposed_indel_signatures = NULL,
    matrix_indel_signatures = NULL,
    matrix_sbs_signatures = NULL,
    hrdetect = NULL,
    onenesstwoness = NULL,
    msisensorpro = NULL,
    genome = "hg19",
    seqnames_loh = c(1:22),
    seqnames_autosomes = c(1:22),
    seqnames_genome_width_or_genome_length = c(1:22, "X", "Y"),
    denoised_coverage_field = "foreground",
    is_visible = TRUE,
	summary = NULL,
    cohort_type = NULL
) {

    # Initialize metadata with all possible columns
    metadata <- initialize_metadata_columns(pair)
    # change NA to NULL
    fix_entries = c("tumor_type", "disease", "primary_site", "inferred_sex", "jabba_gg", "events", "somatic_snvs", "germline_snvs", "tumor_coverage", "estimate_library_complexity", "alignment_summary_metrics", "insert_size_metrics", "wgs_metrics", "het_pileups", "activities_indel_signatures", "decomposed_indel_signatures", "matrix_indel_signatures", "deconstructsigs_sbs_signatures", "activities_sbs_signatures", "decomposed_sbs_signatures","matrix_sbs_signatures","hrdetect", "onenesstwoness", "msisensorpro", "denoised_coverage_field", "summary")
    for (x in fix_entries) {
        if (!exists(x) || is.null(get(x)) || is.na(get(x))) {
            assign(x, NULL)
        }
    }    

    # Add each component sequentially
    metadata <- add_basic_metadata(metadata, tumor_type, disease, primary_site)
    metadata <- add_sex_information(metadata, inferred_sex, jabba_gg, tumor_coverage)
    # Add coverage metrics
    metadata <- add_coverage_metrics(
        metadata = metadata,
        tumor_coverage = tumor_coverage,
        foreground_col_name = foreground_col_name,
        estimate_library_complexity = estimate_library_complexity,
        alignment_summary_metrics = alignment_summary_metrics,
        insert_size_metrics = insert_size_metrics,
        tumor_wgs_metrics = tumor_wgs_metrics,
        normal_wgs_metrics = normal_wgs_metrics
    )
    metadata <- add_variant_counts(metadata, somatic_snvs, genome)
    
    # New SV-related function calls
    metadata <- add_sv_counts(metadata, jabba_gg)
    metadata <- add_purity_ploidy(metadata, purple_pp_bestFit = purple_pp_bestFit, jabba_gg = jabba_gg, tumor_coverage = tumor_coverage)
    # metadata <- add_loh(metadata, jabba_gg, seqnames_loh)
    metadata <- add_fga(metadata, jabba_gg, seqnames_autosomes)
    metadata <- add_genome_length(metadata, jabba_gg, seqnames_genome_width_or_genome_length)
    metadata <- add_sv_types(metadata, jabba_gg, events)
    metadata <- add_coverage_parameters(metadata, tumor_coverage, denoised_coverage_field)
    metadata <- add_het_pileups_parameters(metadata, het_pileups)
    
    # Add TMB calculation
    if (!cohort_type == "heme")
        metadata <- add_tmb(metadata, somatic_snvs, jabba_gg, genome, seqnames_genome_width_or_genome_length)
    
    metadata <- add_signatures(
        metadata,
        activities_sbs_signatures,
        decomposed_sbs_signatures,
        matrix_sbs_signatures,
        activities_indel_signatures,
        decomposed_indel_signatures,
        matrix_indel_signatures,
        deconstructsigs_sbs_signatures
    )
    
    # Add HRD scores
    metadata <- add_hrd_scores(metadata, hrdetect, onenesstwoness)

    # Add MSIsensor score
    metadata <- add_msisensor_score(metadata, msisensorpro)    

    if (!as.logical(is_visible)) {
        metadata$visible <- FALSE
    }

	metadata$summary = summary
    
    return(metadata)
}

#' @name lift_metadata
#' @title lift_metadata
#' @description
#' Create metadata files for all samples in a cohort
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @param genome_length Genome length for the samples (for targeted panels or WES data)
#' @return None
#' @export
lift_metadata <- function(cohort, output_data_dir, cores = 1, genome_length = c(1:22, "X", "Y"), do_lift_datafiles_json = TRUE) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }

    lift_inputs = pairify_cohort_inputs(cohort, tumor_normal_columns = c("wgs_metrics"), keep_remaining = TRUE, sep_cast = "_")

    jabba_column = Skilift::DEFAULT_JABBA(object = cohort)
    # Define all possible columns
    all_cols <- c(
        "pair", "tumor_type", "disease", "primary_site", "inferred_sex",
        # "jabba_gg", 
		jabba_column,
		"events", "oncokb_snv", "somatic_snvs", "germline_snvs", "tumor_coverage",
        "estimate_library_complexity", "alignment_summary_metrics",
        "insert_size_metrics", "tumor_wgs_metrics", "normal_wgs_metrics",
        "het_pileups", "activities_sbs_signatures", "activities_indel_signatures",
        "hrdetect", "onenesstwoness", "msisensorpro", "string_summary"
    )
    
    # Check for required column
    if (!"pair" %in% names(lift_inputs)) {
        stop("Missing required column 'pair' in cohort")
    }
    
    # Warn about missing optional columns
    missing_cols <- all_cols[!all_cols %in% names(lift_inputs)]
    if (length(missing_cols) > 0) {
        warning("Missing optional columns in cohort: ", paste(missing_cols, collapse = ", "))
    }

    if(is.null(genome_length)) {
        warning("No genome length provided, assuming WGS data")
        genome_length <- c(1:22, "X", "Y")
    }
    
    cohort_type = cohort$type
    # Process each sample in parallel
    list_metadata = mclapply(seq_len(nrow(lift_inputs)), function(i) {
        row <- lift_inputs[i,]
        pair_dir <- file.path(output_data_dir, row$pair)
        
        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }
        
        out_file <- file.path(pair_dir, "metadata.json")

        # prefer oncokb_snv over somatic_snvs if available
        is_oncokb_present = !is.null(row$oncokb_snv) && !is.na(row$oncokb_snv)
        snvs_column = row$somatic_snvs
        if (is_oncokb_present) {
            snvs_column = row$oncokb_snv
        }
		
		inferred_sex_field = row$inferred_sex

		purple_qc_path_for_fread = row$purple_qc
		is_purple_qc_null = is.null(purple_qc_path_for_fread) 
		is_purple_pp_range_null = is.null(row$purple_pp_range)
		extracted_purple_qc_path = character(0)
		if (!is_purple_pp_range_null) {
			extracted_purple_qc_path = dir(dirname(as.character(row$purple_pp_range)), full.names = TRUE, pattern = ".qc$")
		}
		if (is_purple_qc_null && NROW(extracted_purple_qc_path) > 0) {
			purple_qc_path_for_fread = extracted_purple_qc_path[1]
		}
		
		is_purple_qc_path_valid = NROW(purple_qc_path_for_fread) == 1 && is.character(purple_qc_path_for_fread) && file.exists(purple_qc_path_for_fread)
		if (is_purple_qc_path_valid) {
			inferred_sex_field = fread(purple_qc_path_for_fread, header = FALSE)[V1 == "AmberGender"]$V2
			inferred_sex_field = tools::toTitleCase(tolower(inferred_sex_field))
		}
        
        futile.logger::flog.threshold("ERROR")
        tryCatchLog({

            # Create metadata object

            metadata <- create_metadata(
                pair = row$pair,
                tumor_type = row$tumor_type,
                disease = row$disease,
                primary_site = row$primary_site,
                inferred_sex = inferred_sex_field,
                purple_pp_bestFit = row$purple_pp_bestFit,
                jabba_gg = row[[jabba_column]],
                events = row$events,
                somatic_snvs = snvs_column,
                germline_snvs = row$germline_snvs,
                tumor_coverage = row$tumor_coverage,
                estimate_library_complexity = row$estimate_library_complexity,
                alignment_summary_metrics = row$alignment_summary_metrics,
                insert_size_metrics = row$insert_size_metrics,
                tumor_wgs_metrics = row$tumor_wgs_metrics,
                normal_wgs_metrics = row$normal_wgs_metrics,
                het_pileups = row$het_pileups,
                activities_sbs_signatures = row$activities_sbs_signatures,
                activities_indel_signatures = row$activities_indel_signatures,
                decomposed_sbs_signatures = row$decomposed_sbs_signatures,
                decomposed_indel_signatures = row$decomposed_indel_signatures,
                matrix_indel_signatures = row$matrix_indel_signatures,
                matrix_sbs_signatures = row$matrix_sbs_signatures,
                hrdetect = row$hrdetect,
                onenesstwoness = row$onenesstwoness,
                msisensorpro = row$msisensorpro,
                seqnames_genome_width_or_genome_length = genome_length,
                denoised_coverage_field = row$denoised_coverage_field,
                is_visible = row$metadata_is_visible,
				summary = row$string_summary,
                cohort_type = cohort_type
            )

            if (is.null(metadata)) {
                print(sprintf("No metadata generated for %s", row$pair))
                return()
            }
            
            # Write to JSON
            jsonlite::write_json(
                metadata,
                out_file,
                auto_unbox = TRUE,
                pretty = TRUE,
                null = "null"
            )

			return(metadata)
            
        }, error = function(e) {
            print(sprintf("Error processing %s: %s", row$pair, e$message))
            NULL
        })
    }, mc.cores = cores, mc.preschedule = TRUE)

	metadata_tbls = rbindlist(list_metadata, fill = TRUE)
	cohort$inputs = Skilift::merge.repl(cohort$inputs, metadata_tbls, by = "pair", prefer_x = TRUE, prefer_y = FALSE)

    # invisible(NULL)

    do_lift_datafiles_json = identical(do_lift_datafiles_json, TRUE)

    if (do_lift_datafiles_json && TRUE) ## excessive checking to make sure nothing dumb happens
        Skilift::lift_datafiles_json(output_data_dir = output_data_dir, cores = cores)
    
	return(cohort)
}

#' @name lift_datafiles_json
#' @title lift_datafiles_json
#' @description
#' Create a combined JSON file for all data files in a directory
#'
#' @param data_dir Directory containing data files
#' @return None
#' @export
lift_datafiles_json <- function(output_data_dir, cores = 1) {
  if (!dir.exists(output_data_dir)) {
    stop("Data directory does not exist.")
  }
  
  # Recursively look for all files named "metadata.json"
  metadata_files <- list.files(
    path = output_data_dir,
    pattern = "metadata\\.json$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(metadata_files) == 0) {
    stop("No metadata.json files found in the specified directory.")
  }
  
  # Read each JSON file and combine them into a list
combined_data <- mclapply(metadata_files, function(file) {
    
    if (!file.exists(file)) {
        warning(sprintf("File does not exist: %s", file))
        return(NULL)
    }
    
    if (!file.access(file, 4) == 0) {
        warning(sprintf("File is not readable: %s", file))
        return(NULL)
    }
    
    if (!grepl("\\.json$", file, ignore.case = TRUE)) {
        warning(sprintf("File is not a JSON file: %s", file))
        return(NULL)
    }
    
    file_info <- file.info(file)
    if (is.na(file_info$size) || file_info$size == 0) {
        warning(sprintf("JSON file is empty: %s", file))
        return(NULL)
    }
    json_content <- jsonlite::fromJSON(file)
    return(jsonlite::unbox(json_content))
    
}, mc.cores = cores)

# Remove NULL entries from failed reads
combined_data <- combined_data[!sapply(combined_data, is.null)]

  # Write the combined JSON list to "datafiles.json" in the data directory
  output_file <- file.path(output_data_dir, "datafiles.json")
  jsonlite::write_json(combined_data, output_file, auto_unbox = TRUE, pretty = TRUE, null = "null")
}
