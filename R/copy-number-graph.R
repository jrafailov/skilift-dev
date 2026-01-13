#' @name lift_copy_number_graph
#' @title lift_copy_number_graph
#' @description
#' Create copy number graph JSON files for all samples in a cohort
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param is_allelic Boolean indicating whether to process allelic (TRUE) or total (FALSE) copy number
#' @param cores Number of cores for parallel processing (default: 1)
#' @return None
#' @export
lift_copy_number_graph <- function(
    cohort,
    output_data_dir,
    settings = Skilift:::default_settings_path,
    is_allelic = FALSE,
    cores = 1
) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }
    
    # Determine which column to use based on is_allelic
    if (is_allelic) {
        message("Processing allelic copy number")
        cn_column <- "allelic_jabba_gg"
        out_filename <- "allelic.json"
    } else {
        message("Processing total copy number")
        cn_column <- Skilift::DEFAULT_JABBA(object = cohort)
        out_filename <- "complex.json"
    }
    
    # Validate required columns exist
    required_cols <- c("pair", cn_column)
    missing_cols <- required_cols[!required_cols %in% names(cohort$inputs)]
    if (length(missing_cols) > 0) {
        stop("Missing required columns in cohort: ", paste(missing_cols, collapse = ", "))
    }
    # Process each sample in parallel
    mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        row <- cohort$inputs[i,]
        pair_dir <- file.path(output_data_dir, row$pair)
        
        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }
        
        out_file <- file.path(pair_dir, out_filename)
        
        futile.logger::flog.threshold("ERROR")
        tryCatchLog({
            ggraph_path <- row[[cn_column]]
            
            if (!is.null(ggraph_path) && file.exists(ggraph_path)) {
                message(sprintf("Reading gGraph for %s", row$pair))
                ggraph <- process_jabba(ggraph_path)

				if (is_allelic) {
                    ggraph = Skilift::collapse_allelic_ggraph(ggraph)
					allele_var = base::get0(
						"allele",
						as.environment(as.list(ggraph$nodes$dt)),
						ifnotfound = rep_len(NA_character_, NROW(ggraph))
					)

					allele_col = ifelse(
						allele_var == "major",
						"#FF000080",
						ifelse(allele_var == "minor", "#0000FF80", NA_character_)
					)

					ggraph$nodes$mark(col = allele_col)
				}
                
                if (!any(class(ggraph) == "gGraph")) {
                    warning(sprintf("Input for %s is not a gGraph object", row$pair))
                    return(NULL)
                }
                
                # Check sequence names overlap with reference
                seq_lengths <- gGnome::parse.js.seqlengths(
                    settings,
                    js.type = "PGV",
                    ref = cohort$reference_name
                )
                
                # Reduce gGraph to only sequences that overlap with reference
                ggraph.reduced <- ggraph[seqnames %in% names(seq_lengths)]
                if (length(ggraph.reduced) == 0) {
                    warning(sprintf(
                        "No overlap between reference sequences and gGraph sequences for %s",
                        row$pair
                    ))
                    return(NULL)
                }
                
                # Set parameters for json export
                mc = mcols(ggraph.reduced$nodes$gr)
                nfields_arg = character(0)
                if ("col" %in% names(mc)) nfields_arg = c(nfields_arg, "col")
                nfields_present = names(mc)[names(mc) %in% c(unlist(row$copy_number_graph_annotations))]
                ggraph.reduced$nodes$mark(Events = NA_character_)
                if (NROW(nfields_present) > 0) {
                    events_string = gGnome:::.dtstring(
                        dt = mc[,nfields_present,drop=FALSE],
                        field_val_sep = ": ",
                        delimiter = "; "
                    )
                    ggraph.reduced$nodes$mark(Events = events_string)
                    ## nfields_arg = c("col", "Events")
                    nfields_arg = c(nfields_arg, "Events")
                }

                ## If still nothing in nfields_arg, set to NULL
                is_nfields_still_empty = NROW(nfields_arg) == 0
                if (is_nfields_still_empty) nfields_arg = NULL
                
                # if("col" %in% names(mcols(ggraph$nodes$gr))) "col" else NULL
                params <- list(
                    filename = out_file,
                    verbose = TRUE,
                    maxcn = row$copy_number_graph_max_cn,
                    # nfields = if("col" %in% names(mcols(ggraph$nodes$gr))) "col" else NULL,
                    nfields = nfields_arg,
                    annotations = unlist(row$copy_number_graph_annotations),
                    node_field_val_sep = ":",
                    node_delimiter = ";"
                )
                
                # Generate and write JSON
                message(sprintf("Writing copy number graph JSON for %s", row$pair))
                do.call(gGnome::refresh(ggraph.reduced)$json, params)
            } else {
                warning(sprintf("Copy number graph file missing for %s", row$pair))
            }
        }, error = function(e) {
            print(sprintf("Error processing %s: %s", row$pair, e$message))
        })
    }, mc.cores = cores, mc.preschedule = TRUE)
    
    invisible(NULL)
}


#' Collapse Allelic graph
#'
#' Reduce number of segments by ignoring phase switches.
#' And redundant segments -- adjacent segments with the same
#' major and minor allele copy numbers
#' 
#' @export
collapse_allelic_ggraph = function(gg) {
    gr = gg$nodes$gr
    gr_major = unname(gr[gr$allele == "major"])
    gr_minor = unname(gr[gr$allele == "minor"])

    is_valid = all(gr_major == gr_minor)
    if (!is_valid) stop("Allelic graph is not collapsible!")

    gr_major$cn_minor = gr_minor$cn
    gr_major$cn_major = gr_major$cn
    gr_minor$cn_major = gr_major$cn
    gr_minor$cn_minor = gr_minor$cn

    gr_by_cn = gUtils::gr_construct_by(c(gr_major, gr_minor), by = c("allele", "cn_minor", "cn_major"))

    gr_collapsed = (
        GenomicRanges::reduce(gr_by_cn, with.revmap = TRUE)
        %>% gUtils::gr_deconstruct_by(meta = TRUE, by = c("allele", "cn_minor", "cn_major"))
        %>% GenomeInfoDb::sortSeqlevels()
        %>% sort()
    )

    mc = mcols(gr_collapsed)
    for (col in c("cn_minor", "cn_major")) {
        mc[[col]] = as.numeric(mc[[col]])
    }
    mcols(gr_collapsed) = mc

    dt_collapsed = setDT(as.data.frame(gr_collapsed))

    nodeid_firstlast = (
        Skilift:::dunlist(dt_collapsed$revmap)
        [, .SD[order(V1)], by = listid]
        [, .(nodeid_left = V1[1], nodeid_right = V1[.N]), keyby = listid]
    )

    mc = mcols(gr_collapsed)
    for (col in names(nodeid_firstlast)) {
        mc[[col]] = nodeid_firstlast[[col]]
    }
    mcols(gr_collapsed) = mc

    gr_collapsed = gr_collapsed %Q% (order(allele, seqnames, start, end))

    gr_collapsed$nodeid_collapsed = seq_len(NROW(gr_collapsed))

    dt_collapsed = setDT(as.data.frame(gr_collapsed))

    dt_collapsed_melted = dt_collapsed %>% data.table::melt.data.table(measure.vars = c("nodeid_left", "nodeid_right"))

    dt_collapsed_melted[, variable := gsub("nodeid_", "", variable)][]

    dtedges = gg$edges$dt
    nr_edges = NROW(dtedges)
    is_any_edges = nr_edges > 0
    edges_collapsed = NULL
    if (is_any_edges) {
        alledges = (
            dtedges[
                (type == "ALT" & class != "REF")
                | (type == "REF")
            ]
        )
        edges_n1 = merge(
            alledges,
            dt_collapsed_melted[, .(allele, variable, value, n1_collapsed = nodeid_collapsed)],
            by.x = c("n1.allele", "n1.side", "n1"), by.y = c("allele", "variable", "value"), 
            all.x = TRUE
        )

        edges_n12 = merge(
            edges_n1,
            dt_collapsed_melted[, .(allele, variable, value, n2_collapsed = nodeid_collapsed)],
            by.x = c("n2.allele", "n2.side", "n2"), by.y = c("allele", "variable", "value"), all.x = TRUE
        )

        edges_collapsed = edges_n12[!is.na(n1_collapsed) & !is.na(n2_collapsed)]

        edges_collapsed$n1___parent = edges_collapsed$n1
        edges_collapsed$n2___parent = edges_collapsed$n2

        edges_collapsed$n1 = edges_collapsed$n1_collapsed
        edges_collapsed$n2 = edges_collapsed$n2_collapsed

    }

    gr_collapsed$col = ifelse(gr_collapsed$allele == "major", "#FF000080", "#0000FF80")
    gr_collapsed$ywid = 0.8

    gr_collapsed$cn = ifelse(gr_collapsed$allele == "major", gr_collapsed$cn_major, gr_collapsed$cn_minor)
    gg_collapsed = gG(nodes = gr_collapsed, edges = edges_collapsed)
    gg_collapsed$set(y.field = "cn")

    return(gg_collapsed)

}
