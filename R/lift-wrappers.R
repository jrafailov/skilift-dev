#' lift wrapper helper function
#'
#' Helper function to check if required columns exist
has_required_columns <- function(cohort, columns, any = FALSE, verbose = TRUE) {
#   is_non_allelic_jabba_column = columns %in% Skilift:::priority_columns_jabba
  is_non_allelic_jabba_column = columns %in% Skilift:::priority_columns_jabba_og
  ix_non_allelic_jabba_column = which(is_non_allelic_jabba_column)
  if (any(is_non_allelic_jabba_column)) {
	replace_column = Skilift::DEFAULT_JABBA(object = cohort)
	original_column = columns[ix_non_allelic_jabba_column[1]]
    is_not_same = !identical(replace_column, original_column)
	if (is_not_same) {
		if (verbose) {
            message(
                "Required jabba column specified: ", original_column, "\n",
                "Replacement jabba column: ", replace_column
            )
        }
        columns[ix_non_allelic_jabba_column[1]] = replace_column
	}
	
  }	
  if (any) {
    return(any(columns %in% names(cohort$inputs)))
  }
  all(columns %in% names(cohort$inputs))
}


  #' Define required columns for each lifter
  required_columns <- list(
    allelic_copy_number_graph = c("allelic_jabba_gg"),
    total_copy_number_graph = c("events"),
    denoised_coverage = c("tumor_coverage"),
    oncotable = c(
      "somatic_variant_annotations",
      "fusions",
      "jabba_gg",
      "karyograph",
      "events",
      "signature_counts",
      "oncokb_snv",
      "oncokb_cna",
      "oncokb_fusions"
    ),
    filtered_events = c("oncotable", "jabba_gg"),
    aggregated_events = "aggregated_events",
    highlighted_events = "highlighted_events",
    hetsnps = c("het_pileups"),
    germline_multiplicity = c("germline_multiplicity"),
    multiplicity = c("multiplicity"),
    segment_width = c("balanced_jabba_gg", "tumor_coverage"),
    allelic_pp_fit = c("jabba_gg", "het_pileups"),
    multiplicity_fits = c("multiplicity"),
    coverage_jabba_cn = c("jabba_gg", "tumor_coverage"),
    purple_sunrise_plot = c("purple_pp_range", "purple_pp_bestFit"),
    pp_plot = c("jabba_gg", "het_pileups"),
    signatures = c(
      "matrix_sbs_signatures",
      "decomposed_sbs_signatures",
      "matrix_indel_signatures",
      "decomposed_indel_signatures"
    ),
    variant_qc = c("somatic_snvs"),
    metadata = c(
      "tumor_type",
      "disease",
      "primary_site",
      "inferred_sex",
      "jabba_gg",
      "events",
      "somatic_snvs",
      "germline_snvs",
      "tumor_coverage",
      "estimate_library_complexity",
      "alignment_summary_metrics",
      "insert_size_metrics",
      "wgs_metrics",
      "het_pileups",
      "activities_sbs_signatures",
      "activities_indel_signatures",
      "hrdetect",
      "onenesstwoness"
    ),
    bam = c(
      "tumor_bam",
      "normal_bam",
      "bam"
    ),
    multiplicity_fits = c(
       "multiplicity", 
       "germline_multiplicity", 
       "hetsnps_multiplicity"
    ),

    allelic_pp_fit = c("jabba_gg", "het_pileups"),

    coverage_jabba_cn = c(
      "jabba_gg", "tumor_coverage"
    ),
    
    purple_sunrise_plot = c(
      "purple_pp_range",
      "purple_pp_bestFit"
    )

  )


#' Run all lift methods on a cohort
#'
#' This function runs all available lift methods on a cohort to generate the complete
#' set of output files needed for visualization.
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @param ... Additional arguments
#'
#' @return None
#' @export
lift_all <- function(
    cohort,
    output_data_dir,
    cores = 1,
    genome_length = NULL,
    width = 10000,
    gencode = Skilift::get_default_gencode(),
    ...) {
  
  list2env(list(...), envir = environment())

  # make sure pair is in cohort
  if (!"pair" %in% names(cohort$inputs)) {
    stop("Missing required column in cohort: pair")
  }

  # make sure output_data_dir is passed and exists, create if not
  if (missing(output_data_dir)) {
    stop("output_data_dir is required")
  } else if (!dir.exists(output_data_dir)) {
    dir.create(output_data_dir, recursive = TRUE)
  }

  if (!is.null(cohort$nextflow_results_path)) {
    oncotable_dir <- file.path(cohort$nextflow_results_path, "oncotable")
  } else {
    oncotable_dir <- file.path(output_data_dir, "oncotable")
    message("Oncotable outputs will be placed in: ", oncotable_dir)
  }


  message("Uploading in ", cohort$type, " mode")

  if (cohort$type == "paired") {
    cohort_mod = lift_paired(
      cohort = cohort,
      output_data_dir = output_data_dir,
      oncotable_dir = oncotable_dir,
      cores = cores,
      genome_length = genome_length,
      width = width,
      gencode = gencode,
      ... = ...
    )
  } else if (cohort$type == "heme") {
    cohort_mod = lift_heme(
      cohort = cohort,
      output_data_dir = output_data_dir,
      oncotable_dir = oncotable_dir,
      cores = cores,
      genome_length = genome_length,
      width = width,
      gencode = gencode,
      ... = ...
    )
  } else if (cohort$type == "tumor_only") {
    cohort_mod = lift_tumor_only(
      cohort = cohort,
      output_data_dir = output_data_dir,
      oncotable_dir = oncotable_dir,
      cores = cores,
      genome_length = genome_length,
      width = width,
      gencode = gencode,
      ... = ...
    )
  }

  datafiles_json_path <- file.path(dirname(output_data_dir), "datafiles.json")
  if (!file.exists(datafiles_json_path)) {
    warning("Creating datafiles.json directory")
  }

  return(cohort_mod)
}

#' Run MVP (Minimum Viable Product) lift methods
#'
#' Helper function containing the common lift methods used across different modes
#'
#' @inheritParams lift_all
#' @param oncotable_dir Directory for oncotable outputs
#'
#' @return Modified cohort object if oncotable is created, otherwise NULL
lift_mvp <- function(
  cohort,
  output_data_dir,
  oncotable_dir,
  cores,
  genome_length = c(1:22, "X", "Y"),
  ...
) {
  list2env(list(...), envir = environment())

  jabba_column = Skilift::DEFAULT_JABBA(object = cohort)

  if (has_required_columns(cohort, required_columns$total_copy_number_graph)) {
    Skilift::skimessage("Uploading total copy number segment gGraph via Skilift::lift_copy_number_graph()")
    Skilift::shutup({
        lift_copy_number_graph(
          cohort = cohort,
          output_data_dir = output_data_dir,
          cores = cores,
          is_allelic = FALSE
        )
    })
  }

  if (has_required_columns(cohort, required_columns$allelic_copy_number_graph)) {
    Skilift::skimessage("Uploading allelic copy number segment gGraph via Skilift::lift_copy_number_graph()")
    Skilift::shutup({
      lift_copy_number_graph(
        cohort = cohort,
        output_data_dir = output_data_dir,
        cores = cores,
        is_allelic = TRUE
      )
    })
  }

  if (has_required_columns(cohort, required_columns$denoised_coverage)) {
    Skilift::skimessage("Uploading binned coverage via Skilift::lift_denoised_coverage()")
    Skilift::shutup({
      lift_denoised_coverage(
        cohort = cohort,
        output_data_dir = output_data_dir,
        cores = cores
      )
    })
  }

  if (has_required_columns(cohort, required_columns$hetsnps)) {
    Skilift::skimessage("Uploading hetsnp pileup via Skilift::lift_hetsnps()")
    Skilift::shutup({
      lift_hetsnps(
        cohort = cohort,
        output_data_dir = output_data_dir,
        cores = cores
      )
    })
  }

  if (has_required_columns(cohort, required_columns$oncotable, any = TRUE) && jabba_column %in% names(cohort$inputs)) {
    Skilift::skimessage("Creating oncotable rds via Skilift::create_oncotable()")
    Skilift::shutup({
      cohort <- create_oncotable(
        cohort = cohort,
        outdir = oncotable_dir,
        cores = cores,
        gencode = gencode
      )
    })
  }

  if (has_required_columns(cohort, required_columns$filtered_events)) {
    Skilift::skimessage("Uploading filtered events table via Skilift::lift_filtered_events()")
    Skilift::shutup({
      cohort <- lift_filtered_events(
        cohort = cohort,
        output_data_dir = output_data_dir,
        cores = cores
      )
    })
  }

    
  # } else if (has_required_columns(cohort, required_columns$filtered_events)) {
  #   Skilift::skimessage("Uploading filtered events table via Skilift::lift_filtered_events()")
  #   Skilift::shutup({
  #     cohort <- lift_filtered_events(
  #       cohort = cohort,
  #       output_data_dir = output_data_dir,
  #       cores = cores
  #     )
  #   })
  # } 

  if (has_required_columns(cohort, required_columns$signatures)) {
    Skilift::skimessage("Uploading signatures via Skilift::lift_signatures()")
    Skilift::shutup({
      lift_signatures(
        cohort = cohort,
        output_data_dir = output_data_dir,
        cores = cores
      )
    })
  }

  if (has_required_columns(cohort, required_columns$multiplicity)) {
    Skilift::skimessage("Uploading somatic multiplicity via Skilift::lift_multiplicity()")
    Skilift::shutup({
      lift_multiplicity(
        cohort = cohort,
        output_data_dir = output_data_dir,
        cores = cores
      )
    })
  }

  if (has_required_columns(cohort, required_columns$segment_width)) {
    Skilift::skimessage("Uploading purity ploidy fit histogram via Skilift::lift_segment_width_distribution()")
    Skilift::shutup({
      lift_segment_width_distribution(
        cohort = cohort,
        output_data_dir = output_data_dir,
        cores = cores
      )
    })
  }

  if (has_required_columns(cohort, required_columns$variant_qc)) {
    Skilift::skimessage("Uploading Variant QC via Skilift::lift_variant_qc()")
    Skilift::shutup({
      lift_variant_qc(
        cohort = cohort,
        output_data_dir = output_data_dir,
        cores = cores
      )
    })
  }

  if (has_required_columns(cohort, Skilift:::required_columns$metadata, any = TRUE)) {
    Skilift::skimessage("Uploading metadata (tumor type, coverage, other QC, HRD score, MSI score, etc) for available fields via Skilift::lift_metadata()")
    Skilift::shutup({
      lift_metadata(
        cohort = cohort,
        output_data_dir = output_data_dir,
        cores = cores,
        genome_length = genome_length
      )
    })

    Skilift::skimessage("Generating datafiles.json via Skilift::lift_datafiles_json()")
    Skilift::shutup({
      lift_datafiles_json(
        output_data_dir = output_data_dir
      )
    })
  }

  # if (has_required_columns(cohort, required_columns$pp_plot)) {
  #   Skilift::skimessage("Uploading  via Skilift::lift_datafiles_json()")
  #   Skilift::shutup({
  #     lift_pp_plot(
  #       cohort = cohort,
  #       output_data_dir = output_data_dir
  #     )
  #   })
  # }

  if (has_required_columns(cohort, Skilift:::required_columns$bam, any = TRUE)) {
    
    Skilift::skimessage("Hard linking bam to data directory via Skilift::lift_bam()")
    Skilift::shutup({
      lift_bam(
        cohort = cohort,
        output_data_dir = output_data_dir,
        cores = cores
      )
    })
  }

  if (has_required_columns(cohort, Skilift:::required_columns$multiplicity_fits, any = TRUE)) {
    Skilift::skimessage("Uploading multiplicity fit histograms via Skilift::lift_multiplicity_fits()")
    Skilift::shutup({
      lift_multiplicity_fits(
        cohort,
        output_data_dir = output_data_dir,
        cores = cores
      )
    })
  }

  if (has_required_columns(cohort, Skilift:::required_columns$coverage_jabba_cn)) {
    Skilift::skimessage("Uploading multiplicity coverage vs copy number scatter+boxplot via Skilift::lift_coverage_jabba_cn()")
    Skilift::shutup({
      lift_coverage_jabba_cn(
        cohort,
        output_data_dir = output_data_dir,
        cores = cores,
        width = width
      )
    })
  }

  if (has_required_columns(cohort, Skilift:::required_columns$purple_sunrise_plot)) {
    Skilift::skimessage("Uploading sunrise plot via Skilift::lift_purple_sunrise_plot()")
    Skilift::shutup({
      lift_purple_sunrise_plot(
        cohort,
        output_data_dir = output_data_dir,
        cores = cores
      )
    })
  }

  if (has_required_columns(cohort, required_columns$allelic_pp_fit)) {
    Skilift::skimessage("Uploading Zi-Ning (Allelic purity ploidy fit) plot via Skilift::lift_allelic_pp_fit()")
    Skilift::shutup({
      lift_allelic_pp_fit(
        cohort = cohort,
        output_data_dir = output_data_dir,
        cores = cores
      )
    })
  }

  return(cohort)

  
}

lift_tumor_only <- function(cohort, output_data_dir, oncotable_dir, cores, ...) {
  cohort <- lift_mvp(cohort, output_data_dir, oncotable_dir, cores, ...)

  return(cohort)
}

lift_heme <- function(cohort, output_data_dir, oncotable_dir, cores, ...) {
  cohort <- lift_mvp(cohort, output_data_dir, oncotable_dir, cores, ...)

  return(cohort)
}

lift_paired <- function(cohort, output_data_dir, oncotable_dir, cores, ...) {
  cohort <- lift_mvp(cohort, output_data_dir, oncotable_dir, cores, ...)

  if (has_required_columns(cohort, required_columns$germline_multiplicity)) {
    Skilift::skimessage("Uploading germline multiplicity via Skilift::lift_multiplicity()")
    Skilift::shutup({
      lift_multiplicity(
        cohort = cohort,
        output_data_dir = output_data_dir,
        cores = cores,
        is_germline = TRUE
      )
    })
  }

  return(cohort)

}

