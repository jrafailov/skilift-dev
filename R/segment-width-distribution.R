#' @title create_ppfit_json
#' @description
#'
#' function to create segmentation plots in case reports
#'
#' @param balanced_jabba_gg the balanced_gg ggraph or the path to it
#' @param tumor_coverage path to associated coverage, if null, will try to pull path from pairs table
#' @param coverage_field field in coverage file to use for segstats
#' @param settings path to settings.json (from PGV)
#' @param ref reference name
#' @param max_na max.na used with JaBbA
#' @param cores cores for JaBbA:::segstats
#' @return segstats table
#' @export
#' @author Stanley Clarke, Tanubrata Dey, Shihab Dider
get_segstats <- function(
    balanced_jabba_gg,
    tumor_coverage,
    coverage_field = "foreground",
    settings = Skilift:::default_settings_path,
    ref = "hg19",
    max_na = 0.9,
    cores = 1) {
    #' drcln_cov = readRDS(thisp$decomposed_cov[i])
    ## make this work with complex where the cov file was not an input and with jabba_gg
    ## x = path_obj %>% sniff %>% inputs %>% select(CovFile, maxna) #get coverage that was used for the jabba run
    if (is.null(balanced_jabba_gg)) {
        stop("Please provide a valid path to a non-integer balanced gGraph file.")
    }

  if (!is.null(tumor_coverage)) {
    if (is.character(tumor_coverage)) {
        cov <- readRDS(tumor_coverage)
    } else {
        cov <- tumor_coverage
    }
  } else {
      stop("Please provide a valid path to a coverage file.")
  }

    ## need to replace NaN with NA or JaBbA:::segstats breaks
    # mcols(cov)[[coverage_field]] <- ifelse(is.nan(mcols(cov)[[coverage_field]]),
    #     NA,
    #     mcols(cov)[[coverage_field]]
    # )
    # mcols(cov)[[coverage_field]] <- as.numeric(mcols(cov)[[coverage_field]])

    ## Modifying coverage field in place to coerce mean = 1
    signal = mcols(cov)[[coverage_field]]
    signal[is.na(signal)] = NA_real_
    signal[is.nan(signal)] = NA_real_
    signal = as.numeric(signal)
    mcols(cov)[[coverage_field]] = signal
    # cov_w = as.numeric(width(cov))
    # cov_w[is.na(signal)] = NA
    # cov_sw = sum(cov_w, na.rm = T)
    # signal_mutl = sum(signal * cov_w, na.rm = T)
    # mcols(cov)[[coverage_field]] = as.numeric(signal * (cov_sw / signal_mutl))

  if (is.character(balanced_jabba_gg)) {
    balanced_jabba_gg = process_jabba(balanced_jabba_gg)
  }
  balanced_gg_gr <- balanced_jabba_gg$nodes$gr
    segstats <- JaBbA:::segstats(
        balanced_gg_gr,
        cov,
        field = coverage_field,
        prior_weight = 1,
        max.chunk = 1e8,
        ## subsample = subsample,
        mc.cores = cores,
        verbose = FALSE,
        max.na = max_na,
        lp = FALSE
    )
    mu = segstats$mean
    mu[is.infinite(mu)] = NA
    w = as.numeric(width(segstats))
    w[is.na(mu)] = NA
    sw = sum(w, na.rm = T)
    mutl = sum(mu * w, na.rm = T)
    segstats$mean = mu * (sw / mutl)
    segstats_dt <- gr2dt(segstats)
    names(segstats_dt) <- gsub("\\.", "_", names(segstats_dt))
    return(segstats_dt)
}

#' @name lift_segment_width_distribution
#' @title lift_segment_width_distribution
#' @description
#' Create segment width distribution JSON files for all samples in a cohort
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param annotations (list) Optional annotations to include in JSON output (default: NULL)
#' @param cores Number of cores for parallel processing (default: 1)
#' @return None
#' @export
lift_segment_width_distribution <- function(
    cohort,
    output_data_dir,
    annotations = NULL,
    cores = 1) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }

    # Validate required columns exist
	jabba_column = Skilift::DEFAULT_JABBA(object = cohort)
    required_cols <- c("pair", jabba_column, "tumor_coverage")
    missing_cols <- required_cols[!required_cols %in% names(cohort$inputs)]
    if (length(missing_cols) > 0) {
        stop("Missing required columns in cohort: ", paste(missing_cols, collapse = ", "))
    }

    ## TODO: Deal with Hypersegmentated data

    # Process each sample in parallel
    mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        row <- cohort$inputs[i, ]
        pair_dir <- file.path(output_data_dir, row$pair)

        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }

        out_file <- file.path(pair_dir, "ppfit.json")

        futile.logger::flog.threshold("ERROR")
        tryCatchLog(
            {
                if (!file.exists(row[[jabba_column]])) {
                    print(sprintf("Balanced JaBbA file not found for %s: %s", row$pair, row[[jabba_column]]))
                    return(NULL)
                }

                # Read the gGraph
                ggraph <- process_jabba(row[[jabba_column]])
                if (!any(class(ggraph) == "gGraph")) {
                    print(sprintf("File is not a gGraph for %s: %s", row$pair, row[[jabba_column]]))
                    return(NULL)
                }

                # Get sequence lengths from the gGraph
                seq_lengths <- seqlengths(ggraph$nodes$gr)

                # Check for required fields
                colnames_check <- c(
                    "start_ix", "end_ix", "eslack_in", "eslack_out",
                    "edges_in", "edges_out", "tile_id", "snode_id",
                    "loose_left", "loose_right", "loose_cn_left",
                    "loose_cn_right", "node_id", "raw_mean", "raw_var",
                    "nbins", "nbins_tot", "nbins_nafrac", "wbins_nafrac",
                    "wbins_ok", "mean", "bad", "max_na", "loess_var",
                    "tau_sq_post", "post_var", "var", "sd"
                )

                if (all(colnames_check %in% names(ggraph$nodes$dt))) {
                    gg_w_segstats <- ggraph
                    fields.keep <- c(colnames_check, "seqnames", "start", "end", "strand", "width", "loose", "index")
                } else {
                    # Get segstats information
                    segstats.dt <- get_segstats(
                        balanced_jabba_gg = row[[jabba_column]],
                        tumor_coverage = row$tumor_coverage
                    )
                    segstats.gr <- GRanges(segstats.dt, seqlengths = seq_lengths) %>% trim()
                    gg_w_segstats <- gG(nodes = segstats.gr, edges = ggraph$edges$dt)
                    fields.keep <- names(segstats.dt) %>% grep("cn", ., invert = TRUE, value = TRUE)
                }

                gg_w_segstats$set(y.field = "cn")

                # Check for sequence name overlap
                ggraph.reduced <- gg_w_segstats[seqnames %in% names(seq_lengths)]
                if (length(ggraph.reduced) == 0) {
                    print(sprintf("No overlap in sequence names for %s", row$pair))
                    return(NULL)
                }

                annotations = row$segment_width_distribution_annotations
                annotations = if (!is.null(annotations)) unlist(annotations) else NULL

                # Create JSON
                gGnome::refresh(ggraph.reduced)$json(
                    filename = out_file,
                    verbose = TRUE,
                    annotations = annotations,
                    maxcn = 500,
                    nfields = fields.keep,
                    save = TRUE
                )
            },
            error = function(e) {
                print(sprintf("Error processing %s: %s", row$pair, e$message))
            }
        )
    }, mc.cores = cores, mc.preschedule = TRUE)

    invisible(NULL)
}

#' @name lift_allelic_pp_fit
#' @title lift_allelic_pp_fit
#' @description This function will create PNG files for allelic ploidy/purity fits
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @param file.name Name of the output file (default: "hetsnps_major_minor.png")
#' @export 
#' @author Johnathan Rafailov
lift_allelic_pp_fit <- function(cohort, 
                                output_data_dir,
                                cores = 1,
                                save_png = TRUE,
                                save_data = TRUE, 
                                binwidth = 1e4,
                                file.name = "hetsnps_major_minor.png") {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }

    # Validate required columns exist
	jabba_column = Skilift::DEFAULT_JABBA(object = cohort)
    required_cols <- c("pair", jabba_column, "het_pileups")
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

        #out_file <- file.path(pair_dir, "allelic_pp_fit.json")
        out_file_png <- file.path(pair_dir, file.name)

        futile.logger::flog.threshold("ERROR")
        tryCatchLog(
            {
                ppplot <- Skilift::pp_plot(jabba_rds = row[[jabba_column]],
                        hets.fname = row$het_pileups,
                        allele = TRUE,
                        scatter = TRUE,
                        binwidth = binwidth,
                        save = FALSE,
                        field = "count",
                        verbose = T) 

                if(save_png){
                    ggsave(file = out_file_png, plot = ppplot$plot, width = 6, height = 6, dpi = 300)
                }
        
                if(save_data){
                    write_json(ppplot$data[seqnames %in% c(1:22, "X", "Y"), .(major.cn, minor.cn, jabba_cn, color)], gsub(".png", ".json", out_file_png), pretty = TRUE)
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


#' @name lift_multiplicity_fits
#' @title lift_multiplicity_fits
#' @description This function will create JSON files for somatic, germline, and hetsnps multiplicity fits
#' -somatic: altered_copies and total_snv_copies
#' -germline: altered_copies and total_snv_copies
#' -hetsnps: major_snv_copies and minor_snv_copies
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @export
lift_multiplicity_fits <- function(cohort,
                                   output_data_dir,
                                   cores = 1) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }

    # Validate required columns exist
    required_cols <- c("pair", "multiplicity")
    optional_cols <-  c("germline_multiplicity", "hetsnps_multiplicity")
    required_cols <- required_cols[!required_cols %in% names(cohort$inputs)]
    iter_cols <- c("multiplicity", optional_cols[which(optional_cols %in% names(cohort$inputs))])

    if (length(required_cols) > 0) {
        stop("Missing required columns in cohort: ", paste(missing_cols, collapse = ", "))
    }

    ## TODO: update me once new plot is ready; right now using Zi plot from skitools
    ## also hacking and just getting rid of the purple sunrise plot to make it work for now
    ## eventually tell CX to update the png we are using on the frontend

    mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        row <- cohort$inputs[i, ]
        pair_dir <- file.path(output_data_dir, row$pair)

        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }

        for (col in iter_cols) {
			is_file_found = ! is.null(row[[col]]) && ! is.na(row[[col]]) && file.exists(row[[col]])
			if (!is_file_found) next
            message(sprintf("Processing %s for %s", col, row$pair))

            out_files <- switch(col,
                multiplicity = file.path(pair_dir, c("somatic_altered_hist.json", "somatic_total_hist.json")),
                germline_multiplicity = file.path(pair_dir, c("germline_altered_hist.json", "germline_total_hist.json")),
                hetsnps_multiplicity = file.path(pair_dir, c("hetsnps_major_hist.json", "hetsnps_minor_hist.json"))
            )

            field_to_use <- switch(col,
                multiplicity = c("altered_copies", "total_snv_copies"),
                germline_multiplicity = c("altered_copies", "total_snv_copies"),
                hetsnps_multiplicity = c("major_snv_copies", "minor_snv_copies")
            )

            futile.logger::flog.threshold("ERROR")
            tryCatchLog(
                {
                    # if (! is_file_found) {
                    #     stop(sprintf("Multiplicity file not found for %s: %s is %s", row$pair, col, row[[col]]))
                    # }

                    mapply(function(out_file, field_to_use) {
                        process_multiplicity_fit(row[[col]],
                            field = field_to_use,
                            out_file = out_file,
                            save_png = F,
                            save_html = F
                        )
                    }, out_files, field_to_use)

                    # if(col == "hetsnps_multiplicity"){
                    #     zi_plot(row[[col]], out_file = file.path(pair_dir, "hetsnps_major_minor.png"))
                    # }

                },
                error = function(e) {
                    print(sprintf("Error processing %s: %s", row$pair, e$message))
                }
            )
        }
    }, mc.cores = cores, mc.preschedule = TRUE)
}


#' @name process_multiplicity_fit
#' @title process_multiplicity_fit
#' @description a method to lift the mutations as a histogram to JSON for viewing in gOS
#' @param variants GRanges object containing the variants
#' @param field field to use from the input data; default is altered_copies. other often used paramaters are ["total_snv_copies","altered_copies", "major_snv_copies", "minor_snv_copies"]
#' @param mask logical value to mask the data or not; default is TRUE
#' @param mask_gr GRanges object containing the mask; default is maskA as provided in the package
#' @param bins  number of bins for histogram; should specify for lower limit to avoid performance issues; default = 100000
#' @param out_file output file path
process_multiplicity_fit <- function(
    variants,
    field = "altered_copies",
    mask = TRUE,
    mask_gr = system.file("extdata", "data", "maskA_re.rds", package = "Skilift"),
    save_data = TRUE,
    save_png = TRUE,
    save_html = FALSE,
    bins = 1e6,
    out_file
) {
    if (is.character(variants)) {
        variants <- readRDS(variants)
    }

    # Read the multiplicity
    if (!any(class(variants) == "GRanges")) {
        stop("input variants is not a GRanges object")
    }

    variants <- variants %Q% (!is.na(cn))
    variants$cn <- round(variants$cn, 0)
    if (mask) {
        mask_gr <- readRDS(mask_gr)
        variants$masked <- variants %^% mask_gr
        # variants <- variants %Q% (masked == TRUE) ## Should be vice versa..
        variants = variants[!variants$masked %in% TRUE]
    }

    # create histogram data
    hist_data <- gr2dt(variants)[, .(count = .N), by = .(mult_cn = get(field), jabba_cn = cn)]
    hist_data[, bin := cut(mult_cn,
        breaks = seq(min(mult_cn, na.rm = T), max(mult_cn, na.rm = T), length.out = bins), include.lowest = TRUE
    )]

    # create binned histogram data
    binned_hist_data <- hist_data[, .(mult_cn = mean(mult_cn, na.rm = T), count = .N), by = .(bin, jabba_cn)][order(mult_cn)][, bin := NULL]

    # write json to file
    if(save_data) {
        write_json(binned_hist_data, out_file, pretty = TRUE)
    }

    # create png
    if(save_png) {

        ggplot(binned_hist_data, aes(x = mult_cn, fill = factor(jabba_cn))) +
            geom_vline(xintercept = seq(0, 10, by = 1), color = "gray", linetype = "dashed") +
            geom_histogram(binwidth = 0.2, color = "black", linewidth = 0.01) +
            scale_fill_manual(values = colors.for.plot) +
            labs(
                x = "Multiplicity",
                y = "Count",
                fill = "JaBbA CN"
            ) +
            theme_bw() +
            xlim(0, 10)

        ggsave(file = gsub(".json", ".png", out_file), width = 6, height = 6, dpi = 300)
    }
    

    invisible(NULL)
}

#' @name lift_coverage_jabba_cn
#' @title lift_coverage_jabba_cn
#' @description This function will create a boxplot and calculate the linear correlation between the foreground from coverage GRanges and JaBbA graph CN.
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @param width tile width of the coverage GRanges  (default: 10kbp)
#' @param mask logical value to mask the data or not; default is TRUE
#' @param save_html Logical value to save HTML files (default: TRUE)
#' @param save_png Logical value to save PNG files (default: TRUE)
#' @param save_data Logical value to save data files (default: TRUE)
#' @return None
#' @author Johnathan Rafailov
#' @importFrom data.table fread setDT
#' @importFrom ggplot2 ggplot geom_boxplot geom_point aes theme_bw xlab ylab ggtitle ggsave
#' @importFrom plotly ggplotly subplot
#' @importFrom htmlwidgets saveWidget
#' @importFrom jsonlite write_json
#' @importFrom gUtils gr.tile si2gr
#' @importFrom GenomicRanges GRanges seqlengths
#' @importFrom grDevices png dev.off
#' @importFrom ggpubr stat_cor theme_bw
#' @export
lift_coverage_jabba_cn <- function(
    cohort,
    output_data_dir,
    cores = 1,
    width = 10000,
    mask_path = system.file("extdata", "data", "maskA_re.rds", package = "Skilift"),
    mask = TRUE,
    save_html = TRUE,
    save_png = TRUE,
    save_data = TRUE) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }

    # Validate required columns exist
	jabba_column = Skilift::DEFAULT_JABBA(object = cohort)
    required_cols <- c("pair", jabba_column, "tumor_coverage")
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

        out_file <- file.path(pair_dir, "coverage_cn_boxplot.json")
        out_file_denoised_png <- file.path(pair_dir, "coverage_cn_boxplot_denoised.png")
        out_file_original_png <- file.path(pair_dir, "coverage_cn_boxplot_original.png")
        out_file_html <- file.path(pair_dir, "coverage_cn_boxplot.html")
                
        futile.logger::flog.threshold("ERROR")
        tryCatchLog(
            {
                if (!file.exists(row[[jabba_column]])) {
                    print(sprintf("Balanced JaBbA file not found for %s: %s", row$pair, row[[jabba_column]]))
                    return(NULL)
                }

                if (!file.exists(row$tumor_coverage)) {
                    print(sprintf("Tumor coverage file not found for %s: %s", row$pair, row$tumor_coverage))
                    return(NULL)
                }

                # Read the gGraph and coverage
                ggraph <- process_jabba(row[[jabba_column]])
                cov <- readRDS(row$tumor_coverage)

                if (!any(class(ggraph) == "gGraph")) {
                    print(sprintf("File is not a gGraph for %s: %s", row$pair, row[[jabba_column]]))
                    return(NULL)
                }

                # Get sequence lengths from the gGraph
                seq_lengths <- seqlengths(ggraph$nodes$gr)

                # Prepare data for plotting
                tiles.gr <- gUtils::gr.tile(gUtils::si2gr(si = seq_lengths), width = width)
                tiles.gr <- tiles.gr %$% cov[, c("foreground", "input.read.counts")] %$% ggraph$nodes$gr[, c("cn")]

                if (mask) {
                    suppressWarnings({
                        mask_gr <- readRDS(mask_path)
                        tiles.gr$masked <- tiles.gr %^% mask_gr
                    })
                } else {
                    tiles.gr$masked <- FALSE
                }

                tiles.gr <- tiles.gr %Q% (!is.na(cn)) %Q% (!is.na(foreground))
                tiles.dt <- gr2dt(tiles.gr)
                tiles.dt$cn <- round(tiles.gr$cn)
                tiles.dt$foreground <- round(tiles.gr$foreground, 1)
                tiles.dt$input.read.counts <- round(tiles.gr$input.read.counts, 1)

                # Select only necessary columns for JSON
                tiles.dt <- tiles.dt[, .(cn, cov = foreground, og = input.read.counts, masked)]
                
                # Save data to JSON if save_data is TRUE
                if (save_data) {
                    write_json(tiles.dt, out_file, pretty = TRUE)
                }

                if (save_png) {
                    # Create PNG files
                    save_coverage_jabba_cn_png(tiles.dt, out_file_denoised_png, out_file_original_png)
                }

                if (save_html) {
                    # Create HTML file only if Cairo library is available
                    if (requireNamespace("Cairo", quietly = TRUE)) {
                        save_coverage_jabba_cn_html(tiles.dt, out_file_html)
                    }
                }   

            },
            error = function(e) {
                print(sprintf("Error processing %s: %s", row$pair, e$message))
            }
        )
    }, mc.cores = cores, mc.preschedule = TRUE)

    invisible(NULL)
}

save_coverage_jabba_cn_html <- function(tiles.dt, out_file_html) {
    # Create ggplot
    p <- ggplot(tiles.dt[masked == F]) +
        geom_boxplot(aes(x = cn, y = cov, fill = masked), alpha = 0.5) +
        geom_point(aes(x = cn, y = cov, color = masked), alpha = 0.5) +
        theme_bw() +
        xlab("JaBbA CN") +
        ylab("Coverage") +
        ggtitle("Coverage vs JaBbA CN")
    
    q <- ggplot(tiles.dt[masked == F]) +
        geom_boxplot(aes(x = cn, y = og, fill = masked), alpha = 0.5) +
        geom_point(aes(x = cn, y = og, color = masked), alpha = 0.5) +
        theme_bw() +
        xlab("JaBbA CN") +
        ylab("Original Coverage") +
        ggtitle("Original Coverage vs JaBbA CN")
    
    # Save as HTML
    p_plotly <- ggplotly(p, width = 800, height = 800)
    q_plotly <- ggplotly(q, width = 800, height = 800)
    subplot(p_plotly, q_plotly, nrows = 1, shareX = FALSE, shareY = FALSE) %>%
    plotly::layout(
        xaxis = list(title = "CN"),
        yaxis = list(title = "Coverage"),
        xaxis2 = list(title = "CN"),
        yaxis2 = list(title = "Original Coverage"),
        autosize = FALSE,
        margin = list(l = 50, r = 50, b = 50, t = 50, pad = 4),
        plot_bgcolor = 'rgba(0,0,0,0)',
        paper_bgcolor = 'rgba(0,0,0,0)'
    ) %>%
    saveWidget(
        file = out_file_html
        , 
        # selfcontained = TRUE # Is this necessary?? You have to have pandoc loaded.
        selfcontained = FALSE
    )
}
save_coverage_jabba_cn_png <- function(tiles.dt, out_file_denoised_png, out_file_original_png) {
    # Create ggplot
    p <- ggplot(tiles.dt[masked == F], aes(x = cn, y = cov)) +
        geom_jitter(size = 0.002, color = "black") +
        geom_boxplot(aes(x = cn, group = cn), color = "red", outlier.shape = NA, width = 0.65, notchwidth = 0.4, alpha = 0) +
        geom_abline(slope = coef(lm(cov ~ cn, data = tiles.dt[masked == F]))[2], intercept = coef(lm(cov ~ cn, data = tiles.dt[masked == F]))[1], color = "blue") +
        geom_abline(slope = 1, intercept = 0, color = "lightblue") +
        stat_cor(method = "pearson") +
        labs(
            title = "Denoised Coverage vs JaBbA CN",
            x = "JaBbA CN",
            y = "Coverage"
        ) +
        theme_bw() +
        theme(legend.position = "none")
    
    q <- ggplot(tiles.dt[masked == F], aes(x = cn, y = og)) +
        geom_jitter(size = 0.002, color = "black") +
        geom_boxplot(aes(x = cn, group = cn), color = "red", outlier.shape = NA, width = 0.65, notchwidth = 0.4, alpha = 0) +
        geom_abline(slope = coef(lm(og ~ cn, data = tiles.dt[masked == F]))[2], intercept = coef(lm(og ~ cn, data = tiles.dt[masked == F]))[1], color = "blue") +
        geom_abline(slope = 1, intercept = 0, color = "lightblue") +
        stat_cor(method = "pearson") +
        labs(
            title = "Original Coverage vs JaBbA CN",
            x = "JaBbA CN",
            y = "Original Coverage"
        ) +
        theme_bw() +
        theme(legend.position = "none")
    
    # Save as PNG
    ggsave(file = out_file_denoised_png, plot = p, width = 6, height = 6, dpi = 1000)
    ggsave(file = out_file_original_png, plot = q, width = 6, height = 6, dpi = 1000)
}

#' @name lift_purple_sunrise_plot
#' @title lift_purple_sunrise_plot
#' @description This function will create JSON files for purple purity/ploidy sunrise plots
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @param save_pngs Logical value to save PNG files (default: TRUE)
#' @param save_html Logical value to save HTML files (default: TRUE)
#' @param save_data Logical value to save data files (default: TRUE)
#' @importFrom data.table fread setDT
#' @importFrom ggplot2 ggplot geom_raster aes scale_fill_scico geom_point geom_segment theme_bw scale_y_continuous scale_x_continuous xlab ylab ggsave
#' @importFrom plotly ggplotly subplot
#' @importFrom htmlwidgets saveWidget
#' @importFrom jsonlite write_json
#' @importFrom scico scale_fill_scico scale_color_scico
#' @export
#' 
#' @references Code adapted from:
#'  - https://github.com/hartwigmedical/hmftools/tree/642436265858083a0bfc81b793a51ccde42edd02/purple/src/main/resources/r/copyNumberPlots.R
#' @author Johnathan Rafailov
lift_purple_sunrise_plot <- function(
    cohort,
    output_data_dir, cores = 1, 
    save_pngs = TRUE, 
    save_html = TRUE, 
    save_data = TRUE
) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }

    # Validate required columns exist
    required_cols <- c("pair", "purple_pp_range", "purple_pp_bestFit")
    if (!all(required_cols %in% names(cohort$inputs))) {
        stop("Missing required columns in cohort: ", paste(missing_cols, collapse = ", "))
    }

    # Process each sample in parallel
    mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        row <- cohort$inputs[i, ]
        pair_dir <- file.path(output_data_dir, row$pair)

        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }

        out_file <- file.path(pair_dir, "purple_sunrise.json")
        out_file_png <- file.path(pair_dir, "purple_sunrise_pp.png")
        out_file_beta_gamma_png <- file.path(pair_dir, "purple_sunrise_beta_gamma.png")
        out_file_html <- file.path(pair_dir, "combined_plot.html")

        futile.logger::flog.threshold("ERROR")
        tryCatchLog(
            {
                if (!file.exists(row$purple_pp_range)) {
                    print(sprintf("Purple purity range file not found for %s: %s", row$pair, row$purple_pp_range))
                    return(NULL)
                }

                if (!file.exists(row$purple_pp_bestFit)) {
                    print(sprintf("Purple best fit file not found for %s: %s", row$pair, row$purple_pp_bestFit))
                    return(NULL)
                }
                range <- fread(row$purple_pp_range)
                fit <- fread(row$purple_pp_bestFit)

                bestPurity <- fit[, purity]
                bestPloidy <- fit[, ploidy]
                bestScore <- fit[, score]

                setorder(range, purity, ploidy)
                range[, absScore := pmin(4, score)]
                range[, score := pmin(1, abs(score - bestScore) / score)]
                range[, leftPloidy := shift(ploidy, type = "lag"), by = purity]
                range[, rightPloidy := shift(ploidy, type = "lead"), by = purity]
                range[, xmin := ploidy - (ploidy - leftPloidy) / 2]
                range[, xmax := ploidy + (rightPloidy - ploidy) / 2]
                range[, ymin := purity - 0.005]
                range[, ymax := purity + 0.005]
                range[, xmin := ifelse(is.na(xmin), ploidy, xmin)]
                range[, xmax := ifelse(is.na(xmax), ploidy, xmax)]
                maxPloidy <- range[, .SD[.N], by = purity][, max(xmax)]
                minPloidy <- range[, .SD[1], by = purity][, min(xmin)]
                minPurity <- min(range$purity)
                maxPurity <- max(range$purity)
                range <- range[xmin <= maxPloidy & xmax >= minPloidy]
                range[, xmax := pmin(xmax, maxPloidy)]
                range[, xmin := pmax(xmin, minPloidy)]

                range[, beta := ploidy / (purity * ploidy + 2 * (1 - purity))]
                range[, gamma := 2 * (1 - purity) / (purity * ploidy + 2 * (1 - purity))]
                bestGamma <- 2 * (1 - bestPurity) / (bestPurity * bestPloidy + 2 * (1 - bestPurity))
                bestBeta <- bestPloidy / (bestPurity * bestPloidy + 2 * (1 - bestPurity))

                range[, color := scales::col_numeric(
                    palette = scico::scico(256, palette = "batlow"),
                    domain = range(1 - score, na.rm = TRUE)
                )(1 - score)]

                # Write the processed data to JSON if save_data is TRUE
                if (save_data) {
                    write_json(range[, .(purity, ploidy, score, xmin, xmax, ymin, ymax, color)], out_file, pretty = TRUE)
                }

                p <- create_purity_ploidy_plot(range, bestPloidy, bestPurity, minPurity, maxPurity, minPloidy, maxPloidy, use_geom_rect = TRUE)
                q <- create_beta_gamma_plot(range, bestBeta, bestGamma)

                if (save_html) {
                    if (requireNamespace("Cairo", quietly = TRUE)) {
                        p_html <- create_purity_ploidy_plot(range, bestPloidy, bestPurity, minPurity, maxPurity, minPloidy, maxPloidy, use_geom_rect = FALSE)
                        save_purple_sunrise_html(p_html, q, out_file_html)
                    }
                }

                if (save_pngs) {
                    save_purple_sunrise_pngs(p, q, out_file_png, out_file_beta_gamma_png)
                }

            },
            error = function(e) {
                print(sprintf("Error processing %s: %s", row$pair, e$message))
            }
        )
    }, mc.cores = cores, mc.preschedule = TRUE)
}

create_purity_ploidy_plot <- function(purple_purity_range, bestPloidy, bestPurity, minPurity, maxPurity, minPloidy, maxPloidy, use_geom_rect = TRUE) {
    if (use_geom_rect) {
        p <- ggplot(purple_purity_range) +
            geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = color)) +
            #scale_fill_scico(palette = "batlow", limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), direction = 1, name = "Relative\nScore", guide = "none") +
            geom_point(aes(x = bestPloidy, y = bestPurity), color = "red", size = 5, shape = 4, stroke = 0.5) +
            geom_segment(aes(x = bestPloidy, xend = bestPloidy, y = minPurity - 0.005, yend = maxPurity + 0.005), color = "red", linetype = "dotted") +
            geom_segment(aes(x = minPloidy, xend = maxPloidy, y = bestPurity, yend = bestPurity), color = "red", linetype = "dotted") +
            theme_bw() +
            theme(panel.background = element_rect(fill = "#f0f2f5", color = NA)) +
            scale_y_continuous(limits = c(minPurity - 0.005, maxPurity + 0.005), labels = c(paste0(minPurity * 100, "%"), "25%", "50%", "75%", "100%"), breaks = c(minPurity, 0.25, 0.5, 0.75, 1), expand = c(0, 0)) +
            scale_x_continuous(limits = c(minPloidy, maxPloidy), breaks = seq(minPloidy, maxPloidy, 0.5), labels = seq(minPloidy, maxPloidy, 0.5) %>% as.character(), expand = c(0, 0)) +
            xlab("Ploidy") +
            ylab("Purity")
        } else {
        p <- ggplot(purple_purity_range) +
            geom_raster(aes(x = ploidy, y = purity, fill = color)) +
            #scale_fill_scico(palette = "batlow", limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), direction = 1, name = "Relative\nScore", guide = "none") +
            geom_point(aes(x = bestPloidy, y = bestPurity), color = "red", size = 5, shape = 4, stroke = 0.5) +
            geom_segment(aes(x = bestPloidy, xend = bestPloidy, y = minPurity - 0.005, yend = maxPurity + 0.005), color = "red", linetype = "dotted") +
            geom_segment(aes(x = minPloidy, xend = maxPloidy, y = bestPurity, yend = bestPurity), color = "red", linetype = "dotted") +
            theme_bw() +
            theme(panel.background = element_rect(fill = "#f0f2f5", color = NA)) +
            scale_y_continuous(limits = c(minPurity - 0.005, maxPurity + 0.005), labels = c(paste0(minPurity * 100, "%"), "25%", "50%", "75%", "100%"), breaks = c(minPurity, 0.25, 0.5, 0.75, 1), expand = c(0, 0)) +
            scale_x_continuous(limits = c(minPloidy, maxPloidy), breaks = seq(minPloidy, maxPloidy, 0.5), labels = seq(minPloidy, maxPloidy, 0.5) %>% as.character(), expand = c(0, 0)) +
            xlab("Ploidy") +
            ylab("Purity")
        }
    return(p)
}

create_beta_gamma_plot <- function(purple_purity_range, bestBeta, bestGamma) {
    minBeta <- min(purple_purity_range$beta)
    maxBeta <- max(purple_purity_range$beta)
    minGamma <- min(purple_purity_range$gamma)
    maxGamma <- max(purple_purity_range$gamma)

    p <- ggplot(purple_purity_range) +
        geom_point(aes(x = beta, y = gamma, fill = score, color = 1 - score), size = 0.5) +
        scale_fill_scico(palette = "batlow", limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), direction = 1, name = "Relative\nScore") +
        scale_color_scico(palette = "batlow", limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), direction = 1, name = "Relative\nScore") +
        theme_bw() +
        theme(panel.background = element_rect(fill = "#f0f2f5", color = NA)) +
        scale_x_continuous(limits = c(minBeta, maxBeta), breaks = seq(minBeta, maxBeta, length.out = 5) %>% round(3)) +
        scale_y_continuous(limits = c(minGamma, maxGamma), breaks = seq(minGamma, maxGamma, length.out = 5) %>% round(3)) +
        geom_point(aes(x = bestBeta, y = bestGamma), color = "red", size = 4, shape = 4, stroke = 0.5) +
        xlab("Beta") +
        ylab("Gamma")
    return(p)
}

save_purple_sunrise_html <- function(p, q, out_file_html) {
    p_plotly <- ggplotly(p, width = 800, height = 800)
    q_plotly <- ggplotly(q, width = 800, height = 800)

    subplot(p_plotly, q_plotly, nrows = 1, shareX = FALSE, shareY = FALSE) %>%
    plotly::layout(
        xaxis = list(title = "Ploidy"),
        yaxis = list(title = "Purity"),
        xaxis2 = list(title = "Beta"),
        yaxis2 = list(title = "Gamma"),
        autosize = FALSE,
        margin = list(l = 50, r = 50, b = 50, t = 50, pad = 4),
        plot_bgcolor = 'rgba(0,0,0,0)',
        paper_bgcolor = 'rgba(0,0,0,0)'
    ) %>%
    saveWidget(
        file = out_file_html
        , 
        # selfcontained = TRUE # is this absolutely necessary?
        selfcontained = FALSE
    )
}

save_purple_sunrise_pngs <- function(p, q, out_file_png, out_file_beta_gamma_png) {
    ggsave(file = out_file_png, plot = p, width = 6, height = 6, dpi = 1000)
    ggsave(file = out_file_beta_gamma_png, plot = q, width = 6, height = 6, dpi = 1000)
}

#' Purity Ploidy Plot
#'
#' Lorem Ipsem
#'
#' @export
#' @author Aditya Deshpande
create_pp_plot <- function(jabba_gg = NA,
                           het_pileups = NULL,
                           seg.fname = NA,
                           field = "count",
                           purity = NA,
                           ploidy = NA,
                           plot.min = -2,
                           plot.max = 2,
                           bins = 500,
                           height = 800,
                           width = 800,
                           output.fname = "./",
                           is.wgs = FALSE) {
    suppressWarnings({
        if (is.na(purity) || is.wgs) {
            if (is.na(jabba_gg) || !file.exists(jabba_gg)) {
                stop("jabba_gg does not exist and no purity and ploidy provided")
            } else {
                jab <- process_jabba(jabba_gg)
                if (is.null(jab$meta$purity) & is.null(jab$purity)) {
                    stop("jabba_gg does not have purity and ploidy values")
                } else {
                    purity <- if (!is.null(jab$meta$purity)) {
                        jab$meta$purity
                    } else {
                        jab$purity
                    }
                    ploidy <- if (!is.null(jab$meta$ploidy)) {
                        jab$meta$ploidy
                    } else {
                        jab$ploidy
                    }
                }
            }
        }

        if (is.null(het_pileups) || !file.exists(het_pileups)) {
            stop("het_pileups not supplied")
        }
        hets <- Skilift:::grab.hets(het_pileups)
        field <- "count"
        if (!field %in% names(values(hets))) {
            stop("hets missing required field")
        }
        if (!field %in% names(mcols(hets))) {
            stop("hets missing required field")
        }
        hets$cn <- skitools::rel2abs(hets, field = field, purity = purity, ploidy = ploidy, allele = TRUE) # added allele == T since this are hets
        eqn <- skitools::rel2abs(hets, field = field, purity = purity, ploidy = ploidy, allele = TRUE, return.params = TRUE) # Trying to keep it as close to Zi's

        seqlevelsStyle(hets) <- "NCBI"

        # if (is.wgs) {
        #     segs = gr.stripstrand(jab$segstats %Q% (strand(jab$segstats)=="+"))[, c()]
        # } else {
        #     segs = readRDS(seg.fname)
        # }

        segs <- jab$nodes$gr[, c()]

        major.segs <- gr.val(segs, hets %Q% (allele == "major"), val = "cn", mean = TRUE, na.rm = TRUE)
        minor.segs <- gr.val(segs, hets %Q% (allele == "minor"), val = "cn", mean = TRUE, na.rm = TRUE)
        if (is.wgs) {
            tiles <- gr.tile(gr = segs, width = 1e4)
            major.tiles <- gr.val(tiles, major.segs, val = "cn", mean = TRUE, na.rm = TRUE)
            minor.tiles <- gr.val(tiles, minor.segs, val = "cn", mean = TRUE, na.rm = TRUE)
            major.tiles.subs <- as.data.table(major.tiles)[sample(.N, 500)]
            minor.tiles.subs <- as.data.table(minor.tiles)[sample(.N, 500)]
            dt <- cbind(
                major.tiles.subs[, .(seqnames, start, end, major.cn = cn)],
                minor.tiles.subs[, .(minor.cn = cn)]
            )
        } else {
            dt <- cbind(
                as.data.table(major.segs)[, .(seqnames, start, end, major.cn = cn)],
                as.data.table(minor.segs)[, .(minor.cn = cn)]
            )
        }

        maxval <- plot.max * ploidy # max dosage
        minval <- plot.min ## min dosage

        dt <- dt[major.cn < maxval & minor.cn < maxval &
            major.cn > minval & minor.cn > minval &
            grepl("[0-9]", seqnames) == TRUE, ]

        dt[, ratio := major.cn / minor.cn]
        return(
            list(
                pp_plot_data = dt,
                maxval = maxval,
                minval = minval,
                purity = purity,
                ploidy = ploidy,
                eqn = eqn
            )
        )
    })
}


#' Lift 2d purity ploidy plot
#' Lift 2d purity ploidy plot
#'
#'
#' @author Aditya Deshpande
#' @export
lift_pp_plot <- function(cohort, output_data_dir, cores = 1) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }

	jabba_column = Skilift::DEFAULT_JABBA(object = cohort)

    # Validate required columns exist
    required_cols <- c("pair", jabba_column, "het_pileups")

    missing_cols <- required_cols[!required_cols %in% names(cohort$inputs)]
    if (length(missing_cols) > 0) {
        stop("Missing required columns in cohort: ", paste(missing_cols, collapse = ", "))
    }

    # Process each sample in parallel
    iterate_function <- function(i) {
        row <- cohort$inputs[i, ]
        pair_dir <- file.path(output_data_dir, row$pair)

        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }

        # out_file <- file.path(pair_dir, "ppfit.json")
        png_path <- paste0(normalizePath(pair_dir), "/pp_plot.png")

        futile.logger::flog.threshold("ERROR")
        tryCatchLog(
            {
                pp_plot_list <- create_pp_plot(
                    jabba_gg = row[[jabba_column]],
                    het_pileups = row$het_pileups
                )
                pp_plot_data <- pp_plot_list$pp_plot_data
                maxval <- pp_plot_list$maxval
                minval <- pp_plot_list$minval
                purity <- pp_plot_list$purity
                ploidy <- pp_plot_list$ploidy
                eqn <- pp_plot_list$eqn

                pt <- ggplot(pp_plot_data, aes(x = major.cn, y = minor.cn)) +
                    scale_x_continuous(
                        breaks = 0:floor(maxval),
                        labels = 0:floor(maxval) %>% as.character(),
                        sec.axis = sec_axis(
                            trans = ~ (. - eqn["intercept"]) / eqn["slope"], # Trying to keep it as close to Zi's
                            name = "Major count"
                        )
                    ) +
                    scale_y_continuous(
                        breaks = 0:floor(maxval),
                        labels = 0:floor(maxval) %>% as.character(),
                        sec.axis = sec_axis(
                            trans = ~ (. - eqn["intercept"]) / eqn["slope"], # Trying to keep it as close to Zi's
                            name = "Minor count"
                        )
                    ) +
                    labs(x = "Major CN", y = "Minor CN") +
                    theme_bw() +
                    theme(
                        legend.position = "none",
                        legend.title = element_text(size = 10, family = "sans"),
                        legend.text = element_text(size = 10, family = "sans"),
                        axis.title = element_text(size = 10, family = "sans"),
                        axis.text.x = element_text(size = 10, family = "sans"),
                        axis.text.y = element_text(size = 10, family = "sans")
                    ) +
                    stat_density_2d(
                        geom = "polygon", contour = TRUE, aes(alpha = 0.5, fill = after_stat(level)),
                        bins = 10
                    ) +
                    scale_fill_distiller(palette = "Blues", direction = 1) +
                    geom_point(size = 2, shape = 4, alpha = 0.3) +
                    ggtitle(paste0("Purity: ", signif(purity, 2), " Ploidy: ", signif(ploidy, 2)))
                grDevices::png(png_path, units = "in", height = 5, width = 5, res = 600)
                print(pt)
                grDevices::dev.off()
            },
            error = function(e) {
                print(sprintf("Error processing %s: %s", row$pair, e$message))
                NULL
            }
        )
    }

    mclapply(
        seq_len(nrow(cohort$inputs)),
        iterate_function,
        mc.cores = cores,
        mc.preschedule = TRUE
    )

    invisible(NULL)
}


#' @name pp_plot
#' @title pp_plot
#'
#' @details
#'
#' create histogram of estimated copy number given purity and ploidy
#'
#' @param jabba_rds (character) JaBbA output
#' @param cov.fname (character) path to coverage GRanges (supply if allele = FALSE)
#' @param hets.fname (character) path to sites.txt (supply if allele = TRUE)
#' @param allele (logical) allelic CN? default FALSE
#' @param field (character) default ratio if allele is FALSE and count if allele is TRUE
#' @param plot.min (numeric) minimum CN default -2
#' @param plot.max (numeric) max CN (factor times ploidy) default 2
#' @param bins (numeric) number of histogram bins default 500
#' @param scatter (logical) default FALSE
#' @param height (numeric) plot height
#' @param width (numeric) plot width
#' @param output.fname (character) path of output directory
#' @param verbose (logical)
#'
#' @return output.fname
#' @author Zi-Ning Choo
#' @export
pp_plot = function(jabba_rds = NULL,
                   cov.fname = NULL,
                   hets.fname = NULL,
                   binwidth = 1e5,
                   allele = FALSE,
                   field = NULL,
                   plot.min = -2,
                   plot.max = 2,
                   scatter = FALSE,
                   bins = 500,
                   height = 800,
                   width = 800,
                   units = "px",
                   output.fname = "./plot.png",
                   save = TRUE,
                   verbose = FALSE) {
  
  if (is.null(jabba_rds) || !file.exists(jabba_rds)) {
    stop("jabba_rds does not exist")
  }

  grab.hets = function(agt.fname = NULL,
                       min.frac = 0.2,
                       max.frac = 0.8)
  {
    if (is.null(agt.fname) || !file.exists(agt.fname)) {
      stop("agt.fname does not exist")
    }

    ## prepare and filter
    agt.dt = fread(agt.fname)
    #only filter on normal fraction if the columns are present within the data table
    if ("alt.frac.n" %in% colnames(agt.dt)) {
      agt.dt = agt.dt[alt.frac.n > min.frac & alt.frac.n < max.frac,]
    } else {
      agt.dt = agt.dt[alt.frac.t > min.frac & alt.frac.t < max.frac,]
    }

    ## add major and minor
    agt.dt[, which.major := ifelse(alt.count.t > ref.count.t, "alt", "ref")]
    agt.dt[, major.count := ifelse(which.major == "alt", alt.count.t, ref.count.t)]
    agt.dt[, minor.count := ifelse(which.major == "alt", ref.count.t, alt.count.t)]

    ## melt the data frame
    agt.melted = rbind(agt.dt[, .(seqnames, start, end, count = major.count, allele = "major")],
                       agt.dt[, .(seqnames, start, end, count = minor.count, allele = "minor")]
                       )

    ## make GRanges
    agt.gr = dt2gr(agt.melted[, .(seqnames, start, end, count, allele)])

    return (agt.gr)
  }
  
  jab = process_jabba(jabba_rds)
  purity <- if (!is.null(jab$purity)) jab$purity else jab$meta$purity
  ploidy <- if (!is.null(jab$ploidy)) jab$ploidy else jab$meta$ploidy
  segstats <- if(!is.null(jab$segstats)) jab$segstats else jab$gr
  
  if (!allele) {
    if (is.null(cov.fname) || !file.exists(cov.fname)) {
      stop("cov.fname not supplied and allele = TRUE")
    }
    if (!grepl(pattern = "rds", x = cov.fname)) {
      stop("cov.fname must be .rds file containing GRanges object")
    }
    cov = readRDS(cov.fname)
    if (!inherits(cov, "GRanges")) {
      stop("cov is not GRanges")
    }
    if (length(cov) == 0) {
      stop("empty GRanges")
    }
    if (is.null(field)) {
      field = "ratio"
    }
    if (!(field %in% names(values(cov)))) {
      stop("cov missing required field")
    }
    if (verbose) {
      message("Grabbing coverage and converting rel2abs")
    }
    cov$cn = rel2abs(cov, field = field, purity = purity, ploidy = ploidy, allele = FALSE)
    cov$jabba_cn = gr.val(cov, segstats, val = "cn", mean = TRUE, na.rm = TRUE)$cn %>% round()
    ## get mean CN over JaBbA segments
    if (verbose) {
      message("computing mean over jabba segments")
    }
    segs = gr.stripstrand(segstats %Q% (strand(segstats)=="+"))[, c()]
    segs = gr.val(segs, cov, val = c("cn", "jabba_cn"), mean = TRUE, na.rm = TRUE)
    if (verbose) {
      message("tiling")
    }
    tiles = gr.tile(gr = segs, width = binwidth)
    tiles = gr.val(tiles, segs[, "cn"], val = c("cn", "jabba_cn"), mean = TRUE, na.rm = TRUE)
    tiles$jabba_cn <- round(tiles$jabba_cn)
    if (verbose) {
      message("Grabbing transformation slope and intercept")
    }
    eqn = rel2abs(cov, field = field, purity = purity, ploidy = ploidy, allele = FALSE, return.params = TRUE)
    dt = as.data.table(tiles)
  } else {
    if (is.null(hets.fname) || !file.exists(hets.fname)) {
      stop("hets.fname not supplied")
    }
    hets = grab.hets(hets.fname)
    if (is.null(field)) {
      field = "count"
    }
    if (!field %in% names(values(hets))) {
      stop("hets missing required field")
    }
    if (verbose) {
      message("Grabbing hets and converting rel2abs")
    }
    hets$cn = rel2abs(hets, field = field, purity = purity, ploidy = ploidy, allele = TRUE)
    hets$jabba_cn = gr.val(hets, segstats, val = "cn", mean = TRUE, na.rm = TRUE)$cn %>% round()
    eqn = rel2abs(hets, field = field, purity = purity, ploidy = ploidy, allele = TRUE, return.params = TRUE)
    if (verbose) {
      message("computing mean over jabba segments")
    }
    segs = gr.stripstrand(segstats %Q% (strand(segstats)=="+"))[, c()]
    major.segs = gr.val(segs, hets %Q% (allele == "major"), val = c("cn", "jabba_cn"), mean = TRUE, na.rm = TRUE)
    minor.segs = gr.val(segs, hets %Q% (allele == "minor"), val = c("cn", "jabba_cn"), mean = TRUE, na.rm = TRUE)
    major.segs$jabba_cn <- round(major.segs$jabba_cn); minor.segs$jabba_cn <- round(minor.segs$jabba_cn)
    if (verbose) {
      message("Tiling")
    }
    tiles = gr.tile(gr = segs, width = binwidth)
    major.tiles = gr.val(tiles, major.segs, val = c("cn", "jabba_cn"), mean = TRUE, na.rm = TRUE)
    minor.tiles = gr.val(tiles, minor.segs, val = c("cn", "jabba_cn"), mean = TRUE, na.rm = TRUE)
    dt = rbind(as.data.table(major.tiles)[, .(seqnames, start, end, allele = "major", cn, jabba_cn)],
               as.data.table(minor.tiles)[, .(seqnames, start, end, allele = "minor", cn, jabba_cn)])
  }

  maxval = plot.max * ploidy # max dosage
  minval = plot.min ## min dosage

  ## remove things with weird ploidy
  dt = dt[cn < maxval & cn > minval & grepl("[0-9]", seqnames)==TRUE]

  if (verbose) {
    message("Making plot for ", nrow(dt), " points")
  }
  
  if (!allele) {

    pt = ggplot(dt, aes(x = cn)) +
      geom_histogram(fill = "black", bins = bins, alpha = 0.8) +
      scale_x_continuous(breaks = 0:floor(maxval),
                         labels = 0:floor(maxval) %>% as.character,
                         sec.axis = sec_axis(trans = ~(. - eqn["intercept"])/eqn["slope"],
                                             name = field)) +
      geom_vline(xintercept = 0:floor(maxval), color = "red", linetype = "longdash") +
      labs(x = "Estimated CN", y = "count") +
      theme_bw() +
      theme(legend.position = "none",
            axis.title = element_text(size = 20, family = "sans"),
            axis.text.x = element_text(size = 20, family = "sans"),
            axis.text.y = element_text(size = 14, family = "sans"))

  } else {


    if (scatter) {
      dt = cbind(as.data.table(major.tiles)[, .(seqnames, start, end, major.cn = cn)],
                 as.data.table(minor.tiles)[, .(minor.cn = cn, jabba_cn = jabba_cn)])
      dt = dt[major.cn < maxval & minor.cn < maxval &
                major.cn > minval & minor.cn > minval &
                grepl("[0-9]", seqnames)==TRUE,]


      # Add a 'color' column with hex codes using Tableau10 color palette
      tableau10 <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
                   "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC")
      dt[, color := { u <- sort(unique(jabba_cn));
        setNames(tableau10[((seq_along(u) - 1) %% length(tableau10)) + 1], u)[as.character(jabba_cn)] }]
      dt <- dt[seqnames %in% c(1:22, "X", "Y")]
    
    pt = ggplot(dt, aes(x = major.cn, y = minor.cn, group = jabba_cn, fill = color, color = color)) +
        #stat_density_2d(geom = "polygon", alpha = 0.7, color = "white") +
        geom_point(size = 1, alpha = 0.2) +
        scale_fill_identity() +
        scale_x_continuous(breaks = 0:floor(maxval),
                           labels = as.character(0:floor(maxval)),
                           sec.axis = sec_axis(~ (. - eqn["intercept"]) / eqn["slope"],
                                               name = paste("Major", field))) +
        scale_y_continuous(breaks = 0:floor(maxval),
                           labels = as.character(0:floor(maxval)),
                           sec.axis = sec_axis(~ (. - eqn["intercept"]) / eqn["slope"],
                                               name = paste("Minor", field))) +
        labs(x = "Major CN", y = "Minor CN", fill = "Group Color", color = "Group Color") +
        theme_bw() +
        theme(legend.position = "none",
              axis.title = element_text(size = 20, family = "sans"),
              axis.text.x = element_text(size = 20, family = "sans"),
              axis.text.y = element_text(size = 14, family = "sans"))
              
    pt = ggExtra::ggMarginal(pt, type = "histogram",
                             xparams = list(bins = bins),
                             yparams = list(bins = bins))
      
    } else {

      pt = ggplot(dt, aes(x = cn)) +
        geom_histogram(fill = "gray", bins = bins, alpha = 0.8) +
        scale_x_continuous(breaks = 0:floor(maxval),
                           labels = 0:floor(maxval) %>% as.character,
                           sec.axis = sec_axis(trans = ~(. - eqn["intercept"])/eqn["slope"],
                                               name = field)) +
        geom_vline(xintercept = 0:floor(maxval), color = "red", linetype = "longdash") +
        labs(x = "Estimated CN", y = "count") +
        facet_grid(row = vars(allele)) +
        theme_bw() +
        theme(legend.position = "none",
              axis.title = element_text(size = 20, family = "sans"),
              axis.text.x = element_text(size = 20, family = "sans"),
              axis.text.y = element_text(size = 14, family = "sans"),
              strip.text.y = element_text(size = 20, family = "sans"))
    }

  }

  if(save) {
    if(verbose) {
      message("Saving results to: ", normalizePath(output.fname))
    }
    grDevices::png(
        filename = normalizePath(output.fname),
        height = height, 
        width = width,
        units = units
    )
    print(pt)
    grDevices::dev.off()
    # ppng(print(pt), 
    #      filename = normalizePath(output.fname), 
    #      height = height, 
    #      width = width,
    #      units = units)
  }

    return(list(
        plot = pt,
        data = dt,
        eqn = eqn
    ))
}
