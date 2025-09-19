pick_first_snpeff = function(
    mutations.dt) {
  snpeff_annotations = mutations.dt$Consequence ## parsed by oncokb
  snpeff_annotations_delim = "," ## parsed by oncokb
  annotationsplit = strsplit(snpeff_annotations, snpeff_annotations_delim)
  annotationsplit = gGnome::dunlist(annotationsplit)
  annotationsplit[, ix := seq_len(.N), by = listid]
  annotationsplit[, num := .N, by = listid]
  first_annotation = annotationsplit[ix == 1]
  mutations.dt.ix = first_annotation$listid
  mutations.dt$snpeff_annotation = NA_character_
  mutations.dt[mutations.dt.ix]$snpeff_annotation = ifelse(
    is.na(first_annotation$V1),
    mutations.dt[mutations.dt.ix]$snpeff_annotation, ## use the original snpeff annotation, not the oncokb parsed snpeff annotation which could map to a different transcript/annotation combo
    first_annotation$V1
  )
  return(mutations.dt)
}

annotate_multihit = function(oncokb_mult) { ## oncokb multiplicity merged output from Skilift::merge_oncokb_multiplicity(oncokb, multiplicity)
  oncokb_mult = Skilift:::pick_first_snpeff(oncokb_mult)
	## Parse Multi-hit status
  is_protein_coding = (
	oncokb_mult$snpeff_annotation %in%
		Skilift:::snpeff_protein_coding_annotations
		# c("frameshift_variant", "stop_lost", "start_lost", "stop_gained", 
		# "chromosome", "exon_loss_variant", "feature_ablation", "duplication", 
		# "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", 
		# "inframe_insertion", "disruptive_inframe_insertion", "inframe_deletion", 
		# "disruptive_inframe_deletion", "coding_sequence_variant", "missense_variant", 
		# "protein_protein_contact", "structural_interaction_variant", 
		# "rare_amino_acid_variant")
	)

	oncokb_mult$is_protein_coding = is_protein_coding

	oncokb_pc_hits = (
		oncokb_mult[, .(
			num_protein_coding_hits_per_gene =
			sum(
			(
				(altered_copies > 0.25 & !is.na(altered_copies))
				| is.na(altered_copies) ## allowing a bit of lenience in case there's edge cases, e.g. PACT data sometimes provides no purity.
			)
			& is_protein_coding
			,
			na.rm = TRUE
			),
			.I
		),
		by = Hugo_Symbol
		]
	)

  setkey(oncokb_pc_hits, I)



	oncokb_mult$num_protein_coding_hits_per_gene = oncokb_pc_hits$num_protein_coding_hits



	## 
	is_vaf_gr65 = oncokb_mult[, VAF > 0.65 & !is.na(VAF)]
	is_loh = oncokb_mult[, segment_cn_low == 0 & !is.na(segment_cn_low)]
	is_hetdel = oncokb_mult[, segment_cn == 1 & !is.na(segment_cn)]

	oncokb_mult$is_vaf_gr65 = is_vaf_gr65
	oncokb_mult$is_loh = is_loh
	oncokb_mult$is_hetdel = is_hetdel


	## || does short-circuit so each element should only be length 1..
	## if this errors out you have a problem.
	oncokb_multi_hits = (
		oncokb_mult
		[, .(
			is_multi_hit_per_gene = (
			any(is_vaf_gr65 & is_protein_coding)
			|| num_protein_coding_hits_per_gene[1] > 1
			|| any(is_hetdel & is_protein_coding)
			|| any(
				(altered_copies[is_protein_coding] / segment_cn[is_protein_coding]) > 0.8
				& !is.na(altered_copies[is_protein_coding])
			)
			),
			.I ## expands query back out to the per variant level
		),
		by = Hugo_Symbol]
	)
  setkey(oncokb_multi_hits, I)


	oncokb_mult$is_multi_hit_per_gene = oncokb_multi_hits$is_multi_hit_per_gene

	return(oncokb_mult)
}

#' @description
#' Creates a multiplicity data table from SNV copy number data, handling both somatic and germline cases.
#'
#' @param snv_cn Path to SNV copy number RDS file or data.table object
#' @param is_germline Logical indicating if mutations are germline (default: FALSE)
#' @param field Column name to use for copy number values (default: "total_copies")
#'
#' @return data.table containing processed mutation data
#' @export
create_multiplicity <- function(snv_cn, oncokb_snv=NULL, is_germline = FALSE, field = "altered_copies", show_only_oncokb = FALSE) {
	is_character_snv = is.character(snv_cn)
	nrows_snv = NROW(snv_cn)
	is_empty_snv = nrows_snv == 0
	is_populated_snv = nrows_snv > 0
	is_lenone_snv = nrows_snv == 1
	is_possible_path_snv = is_character_snv && is_lenone_snv
    is_possible_path_non_na = is_possible_path_snv && !is_loosely_na(snv_cn, other_nas = base::nullfile())
	is_exists_snv = is_possible_path_non_na && file.exists(snv_cn)
	is_notexists_snv = is_possible_path_non_na && !file.exists(snv_cn)
	is_rds_snv = is_exists_snv && grepl("\\.rds$", snv_cn)
	is_other_snv = is_exists_snv && !grepl("\\.rds$", snv_cn)

	mutations.gr = snv_cn
	if (is_rds_snv) {
		mutations.gr = readRDS(snv_cn)
	} else if (is_other_snv) {
		mutations.gr = fread(snv_cn)
	} else if (is_notexists_snv) {
		stop("Please provide a valid path to snv_cn (multiplicity output)")
	}


	#   mutations.gr <- tryCatch({
	#     if (!grepl("\\.rds$", snv_cn)) {
	#       message("Expected .rds ending for mutations. Attempting to read anyway: ", snv_cn)
	#     }
	#     readRDS(snv_cn)
	#   }, error = function(e) {
	#     message(paste0("Input was not .rds; failed with '", e$message, "'\nAssuming input is .maf"))
	#     return(fread(snv_cn) %>% dt2gr())
	#   }, finally = {
	#     message("Finished attempting to load input.")
	#   })

	if (!inherits(mutations.gr, "GRanges")) {
		mutations.gr.tmp = tryCatch(
		{
			dt2gr(mutations.gr)
		}, error = function(e) tryCatch(
		{
			as(mutations.gr, "GRanges")
		}, error = function(e) NULL)
		)
		if (is.null(mutations.gr.tmp)) {
		stop("input snv_cn must be a maf/snv file coercible into GRanges")
		}
		mutations.gr = mutations.gr.tmp
		rm(mutations.gr.tmp)
	}

	if (is.null(mutations.gr)) {
		stop("Failed to assign a valid value to mutations.gr.")
	} else {
		message("Successfully loaded input snv_cn.")
	}


	mcols(mutations.gr)$snpeff_annotation = mutations.gr$annotation
	annotationsplit = strsplit(mcols(mutations.gr)$snpeff_annotation, "&")
	annotationsplit = gGnome::dunlist(annotationsplit)
	annotationsplit[, ix := seq_len(.N), by = listid]
	annotationsplit[, num := .N, by = listid]

	# Normalize everything to just the 1st variant
	# type that appears if we get
	# splice&intron_variant nonsense.
	mcols(mutations.gr)$snpeff_annotation = annotationsplit[ix == 1]$V1
	rm("annotationsplit")

	mut_ann <- ""

	# Annotation fields that should always be present
	# We're using 
	annotation_fields <- list(
		# Variant_Classification = "Type"
		# Let's just use the snpeff annotations here because
		# it should still be possible to populate with just SnpEff/Multiplicity.
		snpeff_annotation = "Type",
		Gene = "Gene",
		# HGVSc = "Variant",
		variant.c = "Variant",
		# HGVSp = "Protein_variant",
		variant.p = "Protein_variant",
		variant.g = "Genomic_variant",
		vaf = "VAF",
		alt = "Alt_count",
		ref = "Ref_count",
		normal.alt = "Normal_alt_count",
		normal.ref = "Normal_ref_count",
		FILTER = "Filter"
	)

	mutations.dt <- gr2dt(mutations.gr)
    
	is_null_oncokb_snv = is.null(oncokb_snv)
	is_path_character = is.character(oncokb_snv)
	is_length_one = NROW(oncokb_snv) == 1
	is_possible_path = is_path_character && is_length_one
	is_file_exists = is_possible_path && file.exists(oncokb_snv)
	is_rds = is_file_exists && grepl("rds$", oncokb_snv)
	is_txt = is_file_exists && grepl("maf$|(c|t)sv$|txt$", oncokb_snv)
	is_oncokb_snv_empty = TRUE

	if (!is_null_oncokb_snv) {
		message("oncokb_snv provided, processing input")
		if (is_rds) {
		oncokb_snv = readRDS(oncokb_snv)
		} else if (is_txt) {
		oncokb_snv = fread(oncokb_snv)
		} else if (is_file_exists) {
		stop("Provided oncokb_snv file extension not supported")
		}

		is_data_frame = inherits(oncokb_snv, "data.frame")
		if (is_data_frame) {
		oncokb_snv_tmp = tryCatch(
			{
			dt2gr(oncokb_snv)
			}, error = function(e) tryCatch(
			{
			as(oncokb_snv, "GRanges")
			}, error = function(e) NULL)
		)
		if (is.null(oncokb_snv_tmp)) {
			stop("oncokb_snv must be coercible to GRanges")
		}
		oncokb_snv = oncokb_snv_tmp
		}
		if (!inherits(oncokb_snv, "GRanges")) {
		stop("final oncokb_snv not a GRanges object")
		}

		nrows_oncokb_snv = NROW(oncokb_snv)
		is_oncokb_snv_empty = nrows_oncokb_snv == 0
	}

	is_show_only_oncokb_irrelevant = is_null_oncokb_snv && identical(show_only_oncokb, TRUE)

	if (is_show_only_oncokb_irrelevant) {
		message("WARNING: show_only_oncokb set to TRUE, but no oncokb provided to show")
	}
  


  

	## FIXME: merge_oncokb_multiplicity is a left join, so if
	## oncokb_snv is empty, mutations.gr.annotated will be empty
	## and here, we need multiplicity to be populated.
	## ideally should be an outer join
	is_oncokb_present_and_populated = !is_null_oncokb_snv && !is_oncokb_snv_empty
	if (is_oncokb_present_and_populated) {
		oncokb.mutations.gr.annotated = merge_oncokb_multiplicity(
		oncokb_snv,
		mutations.gr,
		overwrite = TRUE,
		other.cols.keep = c("snpeff_annotation")
		)

		
		oncokb.mutations.gr.annotated$gene = oncokb.mutations.gr.annotated$Hugo_Symbol

		## If show_only_oncokb is FALSE,
		## Use all multiplicity values, but match oncokb
        
		if (show_only_oncokb) {
			mutations.dt = gr2dt(oncokb.mutations.gr.annotated)
		} else {
			
			nms_oncokb = names(oncokb.mutations.gr.annotated)
			nms_mult = names(mutations.dt)
			cols_to_overwrite = c(
				nms_oncokb[!nms_oncokb %in% nms_mult]
			)
			multiplicity_id_match___ = oncokb.mutations.gr.annotated$multiplicity_id_match
			for (col in cols_to_overwrite) {
				col_class = class(oncokb.mutations.gr.annotated[[col]])
				mutations.dt[[col]] = rep_len(as(NA, col_class), NROW(mutations.dt))
				mutations.dt[multiplicity_id_match___,][[col]] = oncokb.mutations.gr.annotated[[col]]
			}
		} 

		


		## Overwrite with SnpEff annotations pulled out from OncoKB
		## Held in Consequence column.
		
		test_column_valid = function(dt, column) {
			colvals = dt
			is_column_missing = missing(column)
			is_column_null = !is_column_missing && (is.null(column) || NROW(column) == 0)
			is_column_present = !is_column_missing && !is_column_null
			is_dt_scalar = is.null(dim(dt))
			if (!is_dt_scalar && is_column_present) {
				colvals = get0(column, as.environment(as.list(dt)), ifnotfound = NULL)
			} else if (!is_dt_scalar && !is_column_present) {
				stop("A tabular object was provided to dt, please provide column")
			}
			is_null = is.null(colvals)
			is_empty = NROW(colvals) == 0
			is_all_na = !is_null && !is_empty && all(is.na(colvals))
			is_column_valid = !is_null && !is_empty && !is_all_na
			return(is_column_valid)
		}
		snpeff_annotations = mutations.dt$Consequence ## parsed by oncokb
		snpeff_annotations_delim = "," ## parsed by oncokb
		annotationsplit = strsplit(snpeff_annotations, snpeff_annotations_delim)
		annotationsplit = gGnome::dunlist(annotationsplit)
		annotationsplit[, ix := seq_len(.N), by = listid]
		annotationsplit[, num := .N, by = listid]
		first_annotation = annotationsplit[ix == 1]
		mutations.dt.ix = first_annotation$listid
		mutations.dt[mutations.dt.ix]$snpeff_annotation = ifelse(
            is.na(first_annotation$V1),
            mutations.dt[mutations.dt.ix]$snpeff_annotation, ## use the original snpeff annotation, not the oncokb parsed snpeff annotation which could map to a different transcript/annotation combo
            first_annotation$V1
        )
		rm("annotationsplit")

		## Process mutations
		setnames(mutations.dt, old = "VAF", new = "vaf", skip_absent = TRUE)
		mutations.dt <- mutations.dt[!is.na(get(field)), ]
		mutations.dt[start == end, end := end + 1]
		mutations.dt[, vaf := round(vaf, 3)] ## round for frontend legibility

		## FIXME: dealing with the left join issue from merge_oncokb_multiplicity.
		## ideally should be an outer join
		columns_from_oncokb = c("ONCOGENIC", "MUTATION_EFFECT", "HIGHEST_LEVEL")
		for (col in columns_from_oncokb) {
			is_column_valid = test_column_valid(mutations.dt, col)
			col_class = class(oncokb.mutations.gr.annotated[[col]])
			colassign = rep_len(as(NA, col_class), NROW(mutations.dt))
			if (!is_column_valid) mutations.dt[[col]] = colassign
		}
		
		mutations.dt[, ONCOGENIC := fcase(
		is.na(ONCOGENIC), "",
		grepl("Unknown", ONCOGENIC), "", ## necessitated by frontend implementation
		default = ONCOGENIC
		)]
		mutations.dt[, MUTATION_EFFECT := fcase(
		is.na(MUTATION_EFFECT), "",
		grepl("Unknown", MUTATION_EFFECT), "", ## extraneous string
		default = MUTATION_EFFECT
		)]
		mutations.dt[, HIGHEST_LEVEL := fcase(
		HIGHEST_LEVEL == "" | is.na(HIGHEST_LEVEL), "",
		default = gsub("LEVEL_", "", HIGHEST_LEVEL) ## extraneous string
		)]

		mutations.dt <- mutations.dt[FILTER == "PASS"] #### TEMPORARY before implementation of fast coverage

		if ("strand" %in% colnames(mutations.dt)) {
			mutations.dt[, strand := NULL]
		}

		## Add OncoKB specific annotation fields
		annotation_fields = c(
		annotation_fields,
		list(
			ONCOGENIC = "Oncogenicity",
			MUTATION_EFFECT = "Effect",
			HIGHEST_LEVEL = "Level"
		)
		)

		# Converting Gene annotation Hugo_Symbol
		# Meaning don't use ENSG* ids.
		names(annotation_fields)[
		names(annotation_fields) == "Gene"
		] = "Hugo_Symbol"

		# Converting to the OncoKB HGVSc and p variants.
		# for internal consistency
		names(annotation_fields)[
		names(annotation_fields) == "variant.c"
		] = "HGVSc"

		names(annotation_fields)[
		names(annotation_fields) == "variant.p"
		] = "HGVSp"
	} 

  if (!any(class(mutations.dt) == "data.table")) {
    stop("Input must be a data.table.")
  }

  ## parse mut_ann variable
  for (col in names(annotation_fields)) {
    if (col %in% colnames(mutations.dt)) {
      mut_ann <- paste0(mut_ann, annotation_fields[[col]], ": ", mutations.dt[[col]], "; ")
    }
  }

  mutations.dt[, annotation := mut_ann]
  return(mutations.dt[])
}

#' @title Convert Multiplicity Data to Intervals
#' @description
#' Converts multiplicity data into a list of intervals and settings for visualization
#'
#' @param multiplicity data.table containing mutation data with seqnames, start, and end columns
#' @param field Column name to use for y-axis values (default: "total_copies")
#' @param settings Path to settings JSON file
#' @param node_metadata Additional columns to include in node data (default: NULL)
#' @param reference_name Reference genome name (default: "hg19")
#'
#' @return List containing settings and intervals for visualization
#' @export
multiplicity_to_intervals <- function(
    multiplicity,
    field = "altered_copies",
    settings = Skilift:::default_settings_path,
    node_metadata = NULL,
    reference_name = "hg19",
    cohort_type,
	subsample_size = 1e3
) {
    # Load chromosome lengths from settings
    settings_data <- jsonlite::fromJSON(settings)
    chrom_lengths <- as.data.table(
        settings_data$coordinates$sets[[reference_name]])[
        , .(chromosome, startPoint, endPoint)
    ]
    setnames(chrom_lengths, c("seqnames", "start", "end"))

  # Ensure chromosome naming consistency
  if (nrow(chrom_lengths[grepl("chr", seqnames), ]) > 0) {
    chrom_lengths[
      !grepl("chr", seqnames),
      seqnames := paste0("chr", seqnames)
    ]
  }

  # Set y-values from specified field
  multiplicity[, y_value := get(field)]

    # Convert to GRanges
    gr = dt2gr(multiplicity[order(seqnames, start), ]) %>% sortSeqlevels()
    if (nrow(chrom_lengths[grepl("chr", seqnames), ]) > 0) {
        GenomeInfoDb::seqlevelsStyle(gr) = "UCSC"
    } else {
        GenomeInfoDb::seqlevelsStyle(gr) = "NCBI"
    }

    if (is.null(cohort_type)) stop("Cohort type is missing")

    # if (cohort_type == "heme") {
    #     hemedb = readRDS(Skilift:::HEMEDB())
    #     # FIXME: HARDCODED PATH!
    #     # gencode = "/gpfs/data/imielinskilab/DB/GENCODE/gencode.v19.annotation.gtf.nochr.rds"
    #     # gencode <- Skilift:::process_gencode(gencode)
    #     # genes = gencode[gencode$type == "gene"]
    #     # gr_heme_genes = genes[na.omit(match(hemedb$GENE, genes$gene_name))]
    #     message("Filtering multiplicity to heme relevant genes")
    #     # is_heme = (gr %^% gr_heme_genes) & (gr$gene %in% gr_heme_genes$gene_name)
    #     ## Removing dependency on gencode here. 
    #     is_heme = gr$gene %in% hemedb$GENE
    #     gr_in_query = gr[is_heme]
    #     gr_other = gr[!is_heme]
    #     remaining = 1e4 - NROW(gr_in_query)
	# 	is_other_more_than_remaining = NROW(gr_other) > remaining
	# 	is_other_less_eq_than_remaining = !is_other_more_than_remaining
	# 	is_subsampling_required_for_other = remaining > 0 && is_other_more_than_remaining
	# 	is_include_all_other = remaining > 0 && is_other_less_eq_than_remaining
    #     if (is_subsampling_required_for_other) {
    #         otherix = 1:NROW(gr_other)
    #         set.seed(42)
    #         sampled_otherix = sample(otherix, size = remaining, replace = FALSE)
    #         gr_subsampled_other = gr_other[sampled_otherix]
    #         gr = c(gr_in_query, gr_subsampled_other)
    #     } else if (is_include_all_other) {
	# 		gr = c(gr_in_query, gr_other)
    #     } else {
	# 		gr = gr_in_query
	# 	}   
    # } 

	do_subsample = (
		!is.null(subsample_size) 
		&& is.numeric(subsample_size) 
		&& !any(is.na(subsample_size))
		&& subsample_size > 0 
		&& subsample_size < NROW(gr)
	)

	is_cohort_heme = cohort_type == "heme"
	gene_query = character(0)
	
	if (do_subsample) {
		message("Filtering multiplicity to protein-coding genes")
		genc = Skilift::process_gencode()
		transcripts = genc[genc$transcript_type == "protein_coding" & genc$type == "transcript"]
		# transcripts_plus_prom = GenomicRanges::resize(transcripts, width = width(transcripts) + 1000, fix = "end")
		prom = GenomicRanges::flank(transcripts, width = 1000, start = TRUE)
		is_in_query = (gr %^% prom | gr$snpeff_annotation %in% snpeff_protein_coding_annotations)
		if (is_cohort_heme) {
			message("Filtering multiplicity to heme relevant genes")
			hemedb = readRDS(Skilift:::HEMEDB())
			gene_query = hemedb$GENE
			is_in_query = is_in_query & gr$gene %in% gene_query
		}

        gr_in_query = gr[is_in_query]
        gr_other = gr[!is_in_query]
		
        remaining = subsample_size - NROW(gr_in_query)
		is_other_more_than_remaining = NROW(gr_other) > remaining
		is_other_less_eq_than_remaining = !is_other_more_than_remaining
		is_subsampling_required_for_other = remaining > 0 && is_other_more_than_remaining
		is_include_all_other = remaining > 0 && is_other_less_eq_than_remaining
        if (is_subsampling_required_for_other) {
            otherix = 1:NROW(gr_other)
            set.seed(42)
            sampled_otherix = sample(otherix, size = remaining, replace = FALSE)
            gr_subsampled_other = gr_other[sampled_otherix]
            gr = c(gr_in_query, gr_subsampled_other)
        } else if (is_include_all_other) {
			gr = c(gr_in_query, gr_other)
        } else {
			gr = gr_in_query
		}  
	}

  # Validate ranges
  if (any(gr@seqinfo@seqlengths >
    chrom_lengths[seqnames %in% names(seqlengths(gr))]$end)) {
    stop(paste("Ranges exceed chromosome lengths in", reference_name))
  }

  # Create graph object and convert to data.table
  intervals = data.table()
  if (NROW(gr) > 0) {
	jab <- gG(nodes = gr)
	node_cols <- c("snode.id", "y_value", "annotation") ### node_metadata)
	node_dt <- gr2dt(jab$nodes$gr[, node_cols])

	# Create final node data
	intervals <- node_dt[, .(
		chromosome = seqnames,
		startPoint = start,
		endPoint = end,
		iid = snode.id,
		title = snode.id,
		type = "interval",
		y = y_value,
		annotation = node_dt$annotation
	)]
  }

  # Construct return list
  list(
    settings = list(
      y_axis = list(
        title = "copy number",
        visible = TRUE
      )
    ),
    intervals = intervals,
    connections = data.table()
  )
}


#' @name lift_multiplicity
#' @title lift_multiplicity
#' @description
#' Create multiplicity JSON files for all samples in a cohort
#'
#' @param cohort Cohort object containing sample information
#' @param is_germline Logical indicating if mutations are germline (default: FALSE)
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @return None
#' @export
lift_multiplicity <- function(
    cohort,
    output_data_dir,
    is_germline = FALSE,
    node_metadata = c("gene", "feature_type", "annotation", "REF", "ALT", "variant.c", "variant.p", "vaf", "transcript_type", "impact", "rank"),
    field = "altered_copies",
	show_only_oncokb = TRUE,
    cores = 1
) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }
    
    # Determine which column to use based on is_germline
    snv_cn_col <- if(is_germline) "germline_multiplicity" else "multiplicity"

    oncokb_snv_col <- if(is_germline) NULL else "oncokb_snv"
    
    # Validate required column exists
    if (!snv_cn_col %in% names(cohort$inputs)) {
        stop(sprintf("Missing required column in cohort: %s", snv_cn_col))
    }
    
    # Get reference name from cohort
    reference_name <- cohort$reference_name
    if (is.null(reference_name)) {
        stop("Reference name not found in cohort object")
    }
    
    # Process each sample in parallel
    mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        row <- cohort$inputs[i,]
        pair_dir <- file.path(output_data_dir, row$pair)
        
        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }
        
        # Determine output filename based on is_germline
        out_file <- file.path(
            pair_dir,
            if(is_germline) "germline_mutations.json" else "mutations.json"
        )
        
        snv_cn_path = row[[snv_cn_col]]
        oncokb_snv_path = NULL
        if (!is.null(oncokb_snv_col)) oncokb_snv_path = row[[oncokb_snv_col]]
        futile.logger::flog.threshold("ERROR")
        tryCatchLog({
            # Create multiplicity data.table
            mult_dt <- create_multiplicity(
                snv_cn = snv_cn_path,
                oncokb_snv = oncokb_snv_path,
                is_germline = is_germline,
                field = field,
				show_only_oncokb = show_only_oncokb
            )
            
            # Convert to intervals
            intervals_list <- multiplicity_to_intervals(
                multiplicity = mult_dt,
                reference_name = reference_name,
                node_metadata = node_metadata,
                field = field,
                cohort_type = cohort$type
            )
            
            # Write to JSON
            jsonlite::write_json(
                intervals_list,
                out_file,
                pretty = TRUE,
                auto_unbox = TRUE,
                digits = 4
            )
            
        }, error = function(e) {
            print(sprintf("Error processing %s: %s", row$pair, e$message))
            NULL
        })
    }, mc.cores = cores, mc.preschedule = TRUE)
    
    invisible(NULL)
}
