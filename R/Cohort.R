default_settings_path <- system.file("extdata", "test_data", "settings.json", package = "Skilift")


priority_columns_jabba_og = c(
	"events",
	"balanced_jabba_gg",
	"jabba_gg"
)

#' Jabba priority
#' 
#' Non-allelic jabbas
#' @export
priority_columns_jabba = priority_columns_jabba_og

#' Default Jabba Type 
#' 
#' Get correct jabba column to use globally across all methods
#' 
#' @export
DEFAULT_JABBA = function(object, priority_columns_jabba = getOption("skilift_jabba_columns"), get_object = "cohort", test_mode = all, verbose = TRUE) {

	is_priority_columns_jabba_valid = all(priority_columns_jabba %in% Skilift:::priority_columns_jabba_og)
	if (! is_priority_columns_jabba_valid) {
		stop(
			"Provided jabba column names: ", paste(priority_columns_jabba, collapse = ", "), "\n",
			"These must be one of: ", paste(Skilift:::priority_columns_jabba_og, collapse = ", ")
		)
	}
	
	if (!is.function(test_mode)) stop("test_mode must be a function, like all() or any()")
	is_object_not_provided = missing(object) || is.null(object)
	if (is_object_not_provided) {
		if (verbose) message("Pulling Cohort object from environment stack")
		cohort_object = base::dynGet(get_object, ifnotfound = NULL)
	} else {
		cohort_object = object
	}
	inp = NULL
	col = priority_columns_jabba[1]
	use_column = TRUE
	if (!is.null(cohort_object)) {
		inp = cohort_object$inputs
	} else {
		if (verbose) {
			message("Cohort object not found!")
			message("Returning: ", col)
		}
		return(col)
	}
	for (col in priority_columns_jabba) {
		Skilift::test_paths(inp[[col]], verbose = FALSE)
		use_column = test_mode(is_existent_path)
		if (use_column) break
	}
	return(col)
}


#' Default Jabba Type 
#' 
#' Get correct jabba column to use globally across all methods
#' 
#' @export
set_jabba_column = function(column) {
	is_provided_column = (
		!missing(column) 
		&& !is.null(column) 
		&& !(NROW(column) == 1 && Skilift::is_loosely_na(column, other_nas = base::nullfile()))
	)
	is_provided_column_valid = is_provided_column && all(column %in% Skilift:::priority_columns_jabba_og)
	if (is_provided_column && is_provided_column_valid) {
		message("Setting jabba column to: ", column)
	} else if (is_provided_column && !is_provided_column_valid) {
		stop(
			"Provided jabba column names: ", column, "\n",
			"These must be one of: ", paste(Skilift:::priority_columns_jabba_og, collapse = ", ")
		)
	} else if (!is_provided_column) {
		column = Skilift:::priority_columns_jabba_og
		message("Setting default jabba columns to: ", paste(column, collapse = ", "))
	}
	Skilift::assign_in_namespace("priority_columns_jabba", column, ns = asNamespace("Skilift"))
	options("skilift_jabba_columns" = column)
	return(column)
}

#' @export
Cohort <- R6Class("Cohort",
  public = list(
    #' @field inputs data.table containing cohort information
    inputs = NULL,

    #' @field reference_name character string specifying genome reference, defaults to hg19
    reference_name = "hg19",

    #' @field cohort_cols_to_x_cols mapping of cohort columns to possible input columns
    cohort_cols_to_x_cols = NULL,

    #' @field unified list of cohort column names and nf-casereport map
    unified_map = NULL,

    #' @field unified list of cohort column names and nf-casereport map
    path_patterns = NULL,

    #' @field enum("paired", "heme", "tumor_only")
    type = NULL,

    #' @field character
    nextflow_results_path = NULL,

    #' @description
    #' Initialize a new Cohort object
    #' @param x Either a data.table or path to pipeline output directory
    #' @param reference_name character string specifying genome reference
    #' @param col_mapping Optional list mapping cohort columns to possible input columns
    initialize = function(x,
                          reference_name = "hg19",
                          settings = Skilift:::default_settings_path,
                          col_mapping = NULL,
                          path_patterns = Skilift::nf_path_patterns,
                          cohort_type = "paired") {
      self$reference_name <- reference_name

      default_cohort_types <- c("paired", "heme", "tumor_only")
      if (!cohort_type %in% default_cohort_types) {
        stop("cohort_type must be one of: ", paste(default_cohort_types, collapse = ", "))
      }
      self$type <- cohort_type

      # Merge user-provided mapping with default mapping
      default_col_mapping <- Skilift::default_col_mapping
      if (!is.null(col_mapping)) {
        for (col_name in names(col_mapping)) {
          col_mapping
          if (col_name %in% names(default_col_mapping)) {
            # col_mapping goes first to override default_mapping
            default_col_mapping[[col_name]] <- unique(
              c(
                col_mapping[[col_name]],
                default_col_mapping[[col_name]]
              )
            )
          } else {
            default_col_mapping[[col_name]] <- col_mapping[[col_name]]
          }
        }
      }

      self$path_patterns <- path_patterns

      self$cohort_cols_to_x_cols <- default_col_mapping

      if (is.character(x) && length(x) == 1) {
        if (grepl("\\.csv$", x)) {
          self$inputs <- private$construct_from_datatable(data.table(read.csv(x, colClasses = c("patient_id" = "character", "pair" = "character", "patient" = "character"))))[]
        } else {
          self$inputs <- private$construct_from_path(x)[]
          self$nextflow_results_path <- x
          warning("Cohort initialized from path: ", x, "\n",
            "This is deprecated! You should use the gosh-cli to generate an outputs.csv file and then read it in using Cohort$new(path = 'outputs.csv')")
        }
      } else if (is.data.table(x)) {
        self$inputs <- private$construct_from_datatable(x)[]
      } else {
        stop("Input must be either a path (character) or data.table")
      }

      inp = self$inputs
      dm = dim(inp)
      has_dimensions = !is.null(dm)
      jx = integer()
      if (has_dimensions) jx = seq_len(NCOL(inp))
      for (j in jx) {
        v = inp[[j]]
        is_character = is.character(v)
        if (!is_character) next
        trimmedv = trimws(v)
        self$inputs[[j]] = ifelse(Skilift::is_loosely_na(trimmedv) | !nzchar(trimmedv), NA_character_, v)
      }
    },

    #' Print method for cohort object
    print = function(...) {
      cat("Cohort methods:\n")
      cat(format(self, ...), sep = "\n")
      cat("\n")
      base::print(self$inputs)
    },
    #' @description
    #' Validate inputs data.table for missing values and files
    #' @return NULL if all inputs are valid, or a data.table with details about missing data
    validate_inputs = function() {
      if (is.null(self$inputs) || nrow(self$inputs) == 0) {
        stop("No inputs data available to validate")
      }

      # Initialize results data.table
      missing_data <- data.table(
        pair = character(),
        column = character(),
        reason = character()
      )

      # Define metadata fields and configuration parameter fields that should only check for NA/NULL
      metadata_fields <- c("tumor_type", "disease", "primary_site", "inferred_sex", "metadata_is_visible")
      config_fields <- Skilift:::config_parameter_names # Use the existing config parameter names

      # Combine metadata and config fields
      na_check_only_fields <- c(metadata_fields, config_fields)

      # Check each column for missing values
      for (col in names(self$inputs)) {
        # Skip pair column as it's our identifier
        if (col == "pair") next

        # Check for NULL or NA values
        na_pairs <- self$inputs[is.na(get(col)) | is.null(get(col)), pair]
        if (length(na_pairs) > 0) {
          missing_data <- rbindlist(list(
            missing_data,
            data.table(
              pair = na_pairs,
              column = col,
              reason = "NULL or NA value"
            )
          ))
        }

        # For non-metadata/non-config columns that should contain file paths, check if files exist
        if (col %in% names(self$cohort_cols_to_x_cols) && !(col %in% na_check_only_fields)) {
          file_paths <- self$inputs[!is.na(get(col)), get(col)]
          for (path in file_paths) {
            if (!file.exists(path)) {
              missing_data <- rbindlist(list(
                missing_data,
                data.table(
                  pair = self$inputs[get(col) == path, pair],
                  column = col,
                  reason = "File does not exist"
                )
              ))
            }
          }
        }
      }

      # Return results
      if (nrow(missing_data) == 0) {
        message("All inputs are valid - no missing values or files found")
        return(NULL)
      } else {
        # Sort by pair and column for better readability
        setorder(missing_data, pair, column)
        return(missing_data)
      }
    }
  ),
  private = list(
    construct_from_path = function(pipeline_outdir) {
      if (!dir.exists(pipeline_outdir)) {
        stop("Pipeline directory does not exist: ", pipeline_outdir)
      }

      # Get sample metadata first - this contains our patient IDs
      sample_metadata <- private$get_pipeline_samples_metadata(pipeline_outdir)
      sample_metadata$pair_original = sample_metadata$pair

      is_castable = (
          all(!is.na(sample_metadata$sample_type))
      )
      is_paired = is_castable && all(c("normal", "tumor") %in% sample_metadata$sample_type)


      id_to_parse = "pair"
      if (is_castable) {
          sample_metadata = Skilift::dcastski(sample_metadata, id_columns = c("pair", "tumor_type", "disease", "primary_site", "inferred_sex", "pair_original"), type_columns = "sample_type", cast_columns = c("sample", "bam"), sep = "_")
      }

      if (is_paired) {
          sample_metadata$realpair = paste(sample_metadata$tumor_sample, "_vs_", sample_metadata$normal_sample, sep = "")
          id_to_parse = c("realpair", "tumor_sample", "normal_sample")
      }

      is_null_tumor_bam = is.null(sample_metadata$tumor_bam)
      is_all_na_tumor_bam = !is_null_tumor_bam && all(is.na(sample_metadata$tumor_bam))
      is_null_normal_bam = is.null(sample_metadata$normal_bam)
      is_all_na_normal_bam = !is_null_normal_bam && all(is.na(sample_metadata$normal_bam))
      if (is_all_na_tumor_bam) {
          sample_metadata$tumor_bam = NULL
      }
      if (is_all_na_normal_bam) {
          sample_metadata$normal_bam = NULL
      }
      
      

      if (is.null(sample_metadata)) {
        stop("Could not get sample metadata - this is required for patient IDs")
      }

      # Add default configuration parameters to sample_metadata
      for (cohort_col in names(self$cohort_cols_to_x_cols)) {
        possible_cols <- self$cohort_cols_to_x_cols[[cohort_col]]
        default_value <- attr(possible_cols, "default")

        if (!is.null(default_value)) {
          if (is.list(default_value) || length(default_value) > 1) {
            sample_metadata[, (cohort_col) := list(list(default_value))]
          } else {
            sample_metadata[, (cohort_col) := default_value]
          }
        }
      }

      # Get all file paths recursively
      pipeline_output_paths <- list.files(pipeline_outdir, recursive = TRUE, full.names = TRUE)

      outputs_lst = mclapply(id_to_parse, function(id_name) {
          sample_metadata$pair = sample_metadata[[id_name]]
          parse_pipeline_paths(
              pipeline_output_paths,
              initial_dt = sample_metadata,
              path_patterns = self$path_patterns,
              id_name = id_name
          )
      }, mc.cores = length(id_to_parse))
    #   }, mc.cores = 1)


      nm = names(sample_metadata)

      remove_cols = nm[!nm %in% c("pair_original")]

      for (i in seq_len(NROW(outputs_lst))) {
          nm_outputs = names(outputs_lst[[i]])
          outputs_lst[[i]] = base::subset(outputs_lst[[i]], select = !nm_outputs %in% remove_cols)
      }

      outputs = outputs_lst[[1]]
      if (length(outputs_lst) > 1) {
		outputs = Reduce(
			function(x,y) {
				out = merge(
					x, y, 
					by = "pair_original", 
					all = TRUE,
					suffixes = c("", "___THROWOUT")
				)
				out = base::subset(
					out,
					select = ! grepl("___THROWOUT$", names(out))
				)
				return(out)
			}, 
			outputs_lst
		)
	  }

      # outputs = Reduce(function(x,y) merge(x,y, by = "pair_original", all = TRUE), outputs_lst)
      outputs = merge(sample_metadata, outputs, by = "pair_original", all = TRUE)

      ## outputs <- parse_pipeline_paths(
      ##     pipeline_output_paths,
      ##     initial_dt = sample_metadata,
      ##     path_patterns = self$path_patterns
      ## )

      if (nrow(outputs) == 0) {
        warning("No data could be extracted from pipeline directory")
      }

      return(outputs)
    },
    get_pipeline_samples_metadata = function(pipeline_outdir) {
      report_path <- file.path(pipeline_outdir, "pipeline_info/pipeline_report.txt")
      launch_dir = character(0)
      samplesheet_path = character(0)
      if (!file.exists(report_path)) {
        message("Pipeline report not found: ", report_path)
        message("Attempting to read from nextflow logs")
        nflogs = list.files(pipeline_outdir, pattern = "\\.nextflow.*\\.log.*", all.files = TRUE, full.names = TRUE)
        if (NROW(nflogs) == 0) {
          pipeline_outdir_upone = dirname(normalizePath(pipeline_outdir))
          nflogs = list.files(
            pipeline_outdir_upone, 
            pattern = "\\.nextflow.*\\.log.*", 
            all.files = TRUE, 
            full.names = TRUE
          )
        }
        parsed_meta = lapply(nflogs, private$read_nf_log)
        parsed_meta = do.call(Map, c(f = c, parsed_meta))
        launch_dir = unique(parsed_meta$launchdir)
        launch_dir = launch_dir[nzchar(launch_dir)]
        samplesheet_path = unique(parsed_meta$samplesheet)
        samplesheet_path = samplesheet_path[nzchar(samplesheet_path)]
        if (NROW(launch_dir) > 1) stop("More than one launch_dir found - something's wrong")
        if (NROW(samplesheet_path) > 1) stop("More than one samplesheet_path found in pipeline directory")
      } else {
        # Read pipeline report
        report_lines <- readLines(report_path)

        # Extract launch directory and samplesheet path
        launch_dir <- grep("launchDir:", report_lines, value = TRUE)
        samplesheet_path <- grep("input:", report_lines, value = TRUE)

      }

      #### Edge case: multiple runs (with different samplesheets) point to the same results
      #### this needs to be handled a bit more carefully as the nflogs have a cap
      #### solution: coerce unique results directory per run (with common results directory)
      ####  merge the samplesheets with the same parent
      
      # nflogs = list.files(pipeline_outdir, pattern = "\\.nextflow.*\\.log.*", all.files = TRUE, full.names = TRUE)
      # parsed_meta = lapply(nflogs, private$read_nf_log)
      # ## transpose the list
      # parsed_meta = do.call(Map, c(f = c, parsed_meta))
      # samplesheets = unique(parsed_meta$samplesheet)
      # samplesheets = samplesheets[nzchar(samplesheets)]
      # launchdir = unique(parsed_meta$launchdir)
      # launchdir = launchdir[nzchar(launchdir)]
      #
      # # unique_patients = Reduce(union, lapply(samplesheets, function(x) fread(x)$patient))
      #
      # samplesheets = samplesheets[file.exists(samplesheets)]
      # samplesheet = data.table::rbindlist(lapply(samplesheets, fread), fill = TRUE)
      #
      # metadata <- data.table(
      #   pair = samplesheet$patient,
      #   sample = samplesheet$sample,
      #   status = samplesheet$status,
      #   sex = samplesheet$sex,
      #   tumor_type = samplesheet$tumor_type,
      #   disease = samplesheet$disease,
      #   primary_site = samplesheet$primary_site
      # )
      #
      # # remove duplicates
      # metadata <- unique(metadata)
      #
      # return(metadata)

      

      if (length(launch_dir) == 0 || length(samplesheet_path) == 0) {
        warning("Could not find launch directory or samplesheet path in pipeline report")
        return(NULL)
      }


      # Clean up paths
      launch_dir <- gsub(".*launchDir: ", "", launch_dir)
      samplesheet_filename <- gsub(".*input: ", "", samplesheet_path)
      samplesheet_filename <- gsub("^\\./", "", samplesheet_filename)
      if (grepl("^/", samplesheet_filename)) { # TRUE if the samplesheet is already a full path
        samplesheet_path <- samplesheet_filename
      } else {
        samplesheet_path <- file.path(launch_dir, gsub("^\\./", "", samplesheet_filename))
      }

      if (!file.exists(samplesheet_path)) {
        warning("Samplesheet not found: ", samplesheet_path)
        return(NULL)
      }

      

      # Read samplesheet and extract metadata
      samplesheet <- fread(samplesheet_path, colClasses = c("patient" = "character"))
      metacols = c("patient", "sample", "tumor_type", "status", "disease", "primary_site", "sex", "bam")
      metavars = base::mget(
        metacols, 
        as.environment(as.list(samplesheet)),
        ifnotfound = rep_len(
          list(rep_len(NA_character_, NROW(samplesheet))),
          NROW(metacols)
        )
      )
      metadata <- data.table(
        pair = metavars$patient,
        sample = metavars$sample,
        tumor_type = metavars$tumor_type,
        status = metavars$status,
        disease = metavars$disease,
        primary_site = metavars$primary_site,
        inferred_sex = metavars$sex,
        bam = metavars$bam        
      )
      metadata$sample_type = ifelse(metadata$status == 0, "normal", ifelse(metadata$status == 1, "tumor", NA_character_))

      # remove duplicates
     
      metadata <- unique(metadata)

      return(metadata)
    },

    # Determine column based on file path
    construct_from_datatable = function(dt) {
      if (!is.data.table(dt)) {
        stop("Input must be a data.table")
      }

      result_dt <- data.table()
      remaining_dt_cols <- names(dt) # Track available columns

      # First pass: exact matches only
      for (cohort_col in names(self$cohort_cols_to_x_cols)) {
        possible_cols <- self$cohort_cols_to_x_cols[[cohort_col]]

        # Get default value if it exists
        default_value <- attr(possible_cols, "default")

        # Try exact matches first
        found_col <- intersect(possible_cols, remaining_dt_cols)[1]

        if (!is.null(found_col)) {
          result_dt[, (cohort_col) := dt[[found_col]]]
          # Remove matched column from candidates
          remaining_dt_cols <- setdiff(remaining_dt_cols, found_col)
        }
      }

      # Second pass: prefix matches for columns that weren't exactly matched
      for (cohort_col in names(self$cohort_cols_to_x_cols)) {
        if (!is.null(result_dt[[cohort_col]]) && !all(is.na(result_dt[[cohort_col]]))) next # Skip if already matched

        possible_cols <- self$cohort_cols_to_x_cols[[cohort_col]]
        default_value <- attr(possible_cols, "default")

        found_col <- NULL
        # Try prefix matches
        for (col in possible_cols) {
          for (name in remaining_dt_cols) {
            if (startsWith(name, col)) {
              found_col <- name
              break
            }
          }
          if (!is.null(found_col)) break
        }

        if (!is.null(found_col)) {
          result_dt[, (cohort_col) := dt[[found_col]]]
          # Remove matched column from candidates
          remaining_dt_cols <- setdiff(remaining_dt_cols, found_col)
        } else if (!is.null(default_value) && nrow(dt) > 0) {
          # Only add default values if the input data.table is not empty
          if (is.list(default_value) || length(default_value) > 1) {
            result_dt[, (cohort_col) := list(replicate(nrow(dt), list(default_value), simplify = FALSE))]
          } else {
            result_dt[, (cohort_col) := default_value]
          }
        } else if (nrow(dt) > 0 && !cohort_col %in% Skilift:::config_parameter_names) {
          warning(sprintf(
            "No matching column found for '%s'. Expected one of: %s",
            cohort_col, paste(possible_cols, collapse = ", ")
          ))
        }
      }

      if (nrow(result_dt) == 0) {
        warning("No data could be extracted from input data.table")
      }

      if ("pair" %in% names(result_dt) && inherits(result_dt, "data.table")) {
        data.table::setkey(result_dt, pair)
      }

      return(result_dt)
    },
    filter_inputs = function() {
      stop("Method not implemented yet")
    }
  ),
  active = list()
)

Cohort$private_methods[["length"]] = function() {
  return(base::NROW(self$inputs))
}

Cohort$private_methods[["key"]] = function() {
  return(data.table::key(self$inputs))
}

Cohort$active[["length"]] = function() {
  return(private$length())
}

Cohort$active[["getkeys"]] = function() {
  input_key = private$key()
  return(self$inputs[[input_key]])
}

Cohort$public_methods[["setkey"]] = function(..., verbose = getOption("datatable.verbose"), physical = TRUE) {
  Skilift::setkey.Cohort(x = self, ... = ..., verbose = verbose, physical = physical)
  return(self)
}

Cohort$active[["key"]] = function() {
  input_key = private$key()
  return(input_key)
}

#' Cohort length
#' 
#' Get Cohort length (use inputs)
#' @export 
#' @export length.Cohort
length.Cohort = function(x) {
  return(x$length)
}

#' @export 
setkey <- function(x, ...) {
  UseMethod("setkey")
}

#' @export 
setkey.default = function(x, ..., verbose = getOption("datatable.verbose"), physical = TRUE) {
  return(
    data.table::setkey(x, ..., verbose = getOption("datatable.verbose"), physical = TRUE)
  )
}

#' Cohort setkey
#' 
#' setkey on Cohort inputs
#' @export 
#' @export setkey.Cohort
setkey.Cohort = function(x, ..., verbose = getOption("datatable.verbose"), physical = TRUE) {
  return(
    data.table::setkey(x$inputs, ..., verbose = verbose, physical = physical)
  )
}

#' @export 
key <- function(x, ...) {
  UseMethod("key")
}

#' @export 
key.default = function(x) {
  return(data.table::key(x))
}

#' Cohort key
#' 
#' Get key from Cohort inputs
#' @export 
key.Cohort = function(x) {
  return(
    data.table::key(x$inputs)
  )
}

# Cohort$private_fields = list()
# Cohort$private_fields[[".inputs"]] = data.table::data.table()


#' Expected nf-casereport patterns
#'
#' Encoding the expected path variables
#'
#' This variable contains the paths to be expected from
#' the nf-casereport directory structure.
#' You can output the default template to modify it using
#' base::dput(Skilift::nf_path_patterns)
#' @export
nf_path_patterns <- list(
  balanced_jabba_gg = "non_integer_balance/.*/non_integer.balanced.gg.rds$",
  tumor_coverage = "dryclean_tumor/.*/drycleaned.cov.rds$",
  het_pileups = "(hetpileups|amber)/.*/sites.txt$",
  purple_pp_range = "purple/.*purple.purity.range.tsv",
  purple_pp_bestFit = "purple/.*purple.purity.tsv",
  jabba_gg = "jabba/.*/jabba.simple.gg.rds$",
  events = "events/.*/complex.rds$",
  fusions = "fusions/.*/fusions.rds$",
  fusions_junctions = "fusions/.*altedge.annotations.tsv$",
  fragcounter_normal = "fragcounter_normal/.*/.*cov.rds$",
  fragcounter_tumor = "fragcounter_tumor/.*/.*cov.rds$",
  segments_cbs = "cbs/.*/.*(?<!n)seg.rds$",
  structural_variants = c("gridss.*/.*/.*high_confidence_somatic.vcf.bgz$", "tumor_only_junction_filter/.*/.*somatic.filtered.sv.rds$"),
  structural_variants_unfiltered = "gridss.*/.*.gridss.filtered.vcf.gz$",
  karyograph = "jabba/.*/karyograph.rds$",
  allelic_jabba_gg = "lp_phased_balance/.*/lp_phased.balanced.gg.rds$",
  somatic_snvs = c("sage/somatic/tumor_only_filter/.*/.*.sage.pass_filtered.tumoronly.vcf.gz$", "sage/somatic/.*/.*.sage.pass_filtered.vcf.gz$"),
  somatic_snvs_unfiltered = c("sage/somatic/.*/.*sage.somatic.vcf.gz$"),
  somatic_variant_annotations = "snpeff/somatic/.*/.*ann.bcf$",
  multiplicity = c("snv_multiplicity/.*/.*est_snv_cn_somatic.rds", "snv_multiplicity3/.*/.*est_snv_cn_somatic.rds"),
  germline_multiplicity = c("snv_multiplicity/.*/.*est_snv_cn_germline.rds", "snv_multiplicity3/.*/.*est_snv_cn_germline.rds"), ### TO DO FIX ME
  hetsnps_multiplicity = c("snv_multiplicity/.*/.*est_snv_cn_hets.rds", "snv_multiplicity3/.*/.*est_snv_cn_hetsnps.rds"), ### TO DO FIX ME
  ## discard for heme/tumor-only
  activities_sbs_signatures = "signatures/sigprofilerassignment/somatic/.*/sbs_results/Assignment_Solution/Activities/sbs_Assignment_Solution_Activities.txt",
  matrix_sbs_signatures = "signatures/sigprofilerassignment/somatic/.*/SBS/sigmat_results.SBS96.all",
  decomposed_sbs_signatures = "signatures/sigprofilerassignment/somatic/.*/sbs_results/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities*.txt",
  activities_indel_signatures = "signatures/sigprofilerassignment/somatic/.*/indel_results/Assignment_Solution/Activities/indel_Assignment_Solution_Activities.txt",
  matrix_indel_signatures = "signatures/sigprofilerassignment/somatic/.*/ID/sigmat_results.ID83.all",
  decomposed_indel_signatures = "signatures/sigprofilerassignment/somatic/.*/indel_results/.*/Decomposed_MutationType_Probabilities.txt",
  hrdetect = "hrdetect/.*/hrdetect_results.rds",
  ## ^^ discard for heme/tumor-only
  estimate_library_complexity = "qc_reports/gatk/.*/.*metrics",
  alignment_summary_metrics = "qc_reports/picard/.*/.*alignment_summary_metrics",
  insert_size_metrics = "qc_reports/picard/.*/.*insert_size_metrics",
  wgs_metrics = "qc_reports/picard/.*/.*coverage_metrics",
  msisensorpro = c("msisensorpro/(?!.*(_dis|_germline|_somatic|_scan.list))", "msisensorpro/.*/.*msisensor_pro_results.tsv"), ## TODO FILL ME
  oncokb_snv = "oncokb/.*/merged_oncokb.maf",
  oncokb_cna = "oncokb/.*/merged_oncokb_cna.tsv",
  oncokb_fusions = "oncokb/.*/merged_oncokb_fusions.tsv"
)

#' Default column mappings
#'
#' Defines naming convention cohort inputs table
#'
#' This variable contains the possible regex patterns
#' from an input data.table with columns corresponding to processed data.
#' You can output the default template to modify it using
#' base::dput(Skilift::default_col_mapping)
#' @export
default_col_mapping <- list(
  pair = c("patient_id", "pair", "pair_id", "sample"),
  tumor_type = c("tumor_type", "status"),
  tumor_sample = c("tumor_sample"),
  normal_sample = c("normal_sample"),
  disease = c("disease"),
  primary_site = c("primary_site"),
  inferred_sex = c("inferred_sex", "sex"),
  tumor_bam = c("tumor_bam", "bam_tumor"),
  normal_bam = c("normal_bam", "bam_normal"),
  structural_variants = c("structural_variants", "gridss_somatic", "gridss_sv", "svaba_sv", "sv", "svs", "vcf"),
  structural_variants_retiered = c("structural_variants_retiered"),
  structural_variants_unfiltered = "structural_variants_unfiltered",
  structural_variants_raw = "structural_variants_raw",
  tumor_coverage = c("tumor_coverage", "coverage_tumor", "dryclean_tumor", "tumor_dryclean_cov", "dryclean_cov"),
  somatic_snvs = c("somatic_snvs", "snvs_somatic", "sage_somatic_vcf", "strelka_somatic_vcf", "strelka2_somatic_vcf", "somatic_snv", "snv_vcf", "somatic_snv_vcf", "snv_somatic_vcf"),
  somatic_snvs_unfiltered = c("somatic_snvs_unfiltered", "snvs_somatic_unfiltered"),
  germline_snvs = c("germline_snvs", "snvs_germline", "sage_germline_vcf", "germline_snv", "germline_snv_vcf"),
  fragcounter_normal = c("fragcounter_normal", "frag_cov_normal"),
  fragcounter_tumor = c("fragcounter_tumor", "frag_cov_tumor"),
  segments_cbs = c("cbs_seg_rds", "seg_rds", "cbs_seg", "seg"),
  het_pileups = c("het_pileups", "hets", "sites_txt", "hets_sites"),
  multiplicity = c("multiplicity", "somatic_snv_cn", "snv_multiplicity"),
  germline_multiplicity = c("germline_multiplicity", "multiplicity_germline", "germline_snv_cn"),
  hetsnps_multiplicity = c("hetsnps_multiplicity", "multiplicity_hetsnps", "hets_snv_cn"),
  somatic_variant_annotations = c("somatic_variant_annotations", "variant_annotations_somatic", "annotated_bcf", "variant_somatic_ann"),
  germline_variant_annotations = c("germline_variant_annotations",  "variant_annotations_germline","annotated_vcf_germline"),
  oncokb_snv = c("oncokb_snv", "oncokb_maf", "maf"),
  oncokb_cna = c("oncokb_cna", "cna"),
  oncokb_fusions = c("oncokb_fusions", "oncokb_fusion", "fusion_maf"),
  itdseek_vcf = c("itdseek_vcf"),
  itdseek_rds = c("itdseek_rds"),
  jabba_gg = c("jabba_gg", "jabba_simple", "jabba_rds", "jabba_simple_gg"),
  karyograph = c("karyograph"),
  balanced_jabba_gg = c("balanced_jabba_gg", "jabba_gg_balanced", "non_integer_balance", "balanced_gg", "ni_balanced_gg"),
  events = c("events", "complex"),
  fusions = c("fusions"),
  fusions_junctions = c("fusions_junctions"),
  allelic_jabba_gg = c("allelic_jabba_gg", "jabba_gg_allelic", "lp_phased_balance", "allelic_gg", "lp_balanced_gg"),
  activities_sbs_signatures = c("activities_sbs_signatures", "signatures_activities_sbs", "sbs_activities"),
  matrix_sbs_signatures = c("matrix_sbs_signatures", "signatures_matrix_sbs", "sbs_matrix"),
  decomposed_sbs_signatures = c("decomposed_sbs_signatures", "signatures_decomposed_sbs", "sbs_decomposed"),
  sbs_post_prob_signatures = c("signatures_post_prob_sbs"),
  activities_indel_signatures = c("activities_indel_signatures", "signatures_activities_indel", "indel_activities"),
  matrix_indel_signatures = c("matrix_indel_signatures", "signatures_matrix_indel", "indel_matrix"),
  decomposed_indel_signatures = c("decomposed_indel_signatures", "signatures_decomposed_indel", "indel_decomposed"),
  indel_post_prob_signatures = c("signatures_post_prob_indel"),
  ffpe_impact_vcf = "ffpe_impact_vcf",
  ffpe_impact_vcf_filtered = "ffpe_impact_vcf_filtered",
  hrdetect = c("hrdetect", "hrd"),
  onenesstwoness = c("onenesstwoness","oneness_twoness"),
  oncotable = c("oncotable"),
  estimate_library_complexity = c("estimate_library_complexity", "qc_dup_rate", "library_complexity_metrics", "est_lib_complex", "qc_dup_rate_tumor"),
  estimate_library_complexity_tumor = c("estimate_library_complexity", "estimate_library_complexity_tumor", "qc_dup_rate", "qc_dup_rate_tumor", "library_complexity_metrics_tumor", "est_lib_complex_tumor"),
  estimate_library_complexity_normal = c("estimate_library_complexity_normal", "qc_dup_rate_normal", "library_complexity_metrics_normal", "est_lib_complex_normal"),
  alignment_summary_metrics = c("alignment_summary_metrics", "alignment_summary_metrics_tumor",  "qc_alignment_summary", "alignment_metrics", "qc_alignment_summary_tumor"),
  alignment_summary_metrics_tumor = c("alignment_summary_metrics", "alignment_summary_metrics_tumor", "qc_alignment_summary_tumor", "alignment_metrics_tumor"),
  alignment_summary_metrics_normal = c("alignment_summary_metrics_normal", "qc_alignment_summary_normal", "alignment_metrics_normal"),
  insert_size_metrics = c("insert_size_metrics", "qc_insert_size", "qc_insert_size_tumor", "insert_metrics", "insert_size_metrics_tumor"),
  insert_size_metrics_tumor = c("insert_size_metrics_tumor", "qc_insert_size_tumor", "insert_metrics_tumor", "insert_metrics"),
  insert_size_metrics_normal = c("insert_size_metrics_normal", "qc_insert_size_normal", "insert_metrics_normal"),
  tumor_wgs_metrics = c("tumor_wgs_metrics", "qc_coverage_metrics", "wgs_metrics", "qc_coverage_metrics_tumor", "tumor_wgs_stats"),
  wgs_metrics = c("wgs_metrics", "qc_coverage_metrics", "tumor_wgs_metrics", "qc_coverage_metrics_tumor","wgs_stats"),
  normal_wgs_metrics = c("normal_wgs_metrics", "qc_coverage_metrics_normal",  "normal_wgs_stats"),
  purple_pp_range = c("purple_pp_range", "purple_range"),
  purple_qc = c("purple_qc"),
  purple_pp_bestFit = c("purple_pp_bestFit", "purple_pp_best_fit", "purple_bestFit", "purple_solution"),
  msisensorpro = c("msisensorpro", "msisensor_pro", "msisensor_pro_results", "msisensor_results"),
  conpair_concordance = c("conpair_concordance"),
  conpair_contamination = c("conpair_contamination"),
  # Configuration parameters with default values
  metadata_is_visible = structure(c("metadata_is_visible"), default = TRUE),
  copy_number_graph_max_cn = structure(c("copy_number_graph_max_cn"), default = 100),
  copy_number_graph_annotations = structure(c("copy_number_graph_annotations"), default = list(c("bfb", "chromoplexy", "chromothripsis", "del", "dm", "cpxdm", "dup", "pyrgo", "rigma", "simple", "tic", "tyfonas"))),
  multiplicity_node_metadata = structure(c("multiplicity_node_metadata"), default = c("gene", "feature_type", "annotation", "REF", "ALT", "variant.c", "variant.p", "vaf", "transcript_type", "impact", "rank")),
  multiplicity_field = structure(c("multiplicity_field"), default = "altered_copies"),
  denoised_coverage_apply_mask = structure(c("denoised_coverage_apply_mask"), default = TRUE),
  denoised_coverage_field = structure(c("denoised_coverage_field"), default = "foreground"),
  denoised_coverage_color_field = structure(c("denoised_coverage_color_field"), default = NULL),
  denoised_coverage_bin_width = structure(c("denoised_coverage_bin_width"), default = 1e3L), ## More headaches with testing setting to NA than worth changing
  hetsnps_field = structure(c("hetsnps_field"), default = "count"),
  hetsnps_color_field = structure(c("hetsnps_color_field"), default = "col"),
  hetsnps_bin_width = structure(c("hetsnps_bin_width"), default = NA),
  hetsnps_mask = structure(c("hetsnps_mask"), default = system.file("extdata", "data", "maskA_re.rds", package = "Skilift")),
  hetsnps_subsample_size = structure(c("hetsnps_subsample_size"), default = 100000),
  hetsnps_min_normal_freq = structure(c("hetsnps_min_normal_freq"), default = 0.2),
  hetsnps_max_normal_freq = structure(c("hetsnps_max_normal_freq"), default = 0.8),
  segment_width_distribution_annotations = structure(c("segment_width_distribution_annotations"), default = NULL)
)

#' Configuration parameter names
#'
#' List of configuration parameter column names
#'
#' This variable contains the names of columns that are used
#' for configuration parameters in the Cohort inputs data.table
#' @export
config_parameter_names <- c(
  "copy_number_graph_max_cn",
  "copy_number_graph_annotations",
  "multiplicity_node_metadata",
  "multiplicity_field",
  "denoised_coverage_apply_mask",
  "denoised_coverage_field",
  "denoised_coverage_color_field",
  "denoised_coverage_bin_width",
  "hetsnps_field",
  "hetsnps_color_field",
  "hetsnps_bin_width",
  "hetsnps_mask",
  "hetsnps_subsample_size",
  "hetsnps_min_normal_freq",
  "hetsnps_max_normal_freq",
  "segment_width_distribution_annotations"
)

#' Create unified column and nf casereports map
#'
#' Parse column mapping and nf-casereport path mapping
#'
#' Takes in a list of named column regex patterns and
#' a list of named nf-casereport path regex patterns and
#' coerces them into a unified object for later use.
#' To be added as an active binding.
#' Not exported.
unify_colmap_nfmap <- function(col_mapping, path_patterns) {
  if (!inherits(col_mapping, "list") || !inherits(path_patterns, "list")) {
    stop("col_mapping and path_patterns objects must be lists")
  }
  expected_map_fields <- c(
    "column_name_regex",
    "nf_path_regex"
  )

  unified_map <- list()

  all_col_names <- union(
    names(col_mapping),
    names(path_patterns)
  )

  for (col in all_col_names) {
    unified_map[[col]] <- list(
      column_name_regex = default_col_mapping[[col]],
      nf_path_regex = nf_path_patterns[[col]]
    )
  }
  return(unified_map)
}

Cohort$private_methods[["unify_colmap_nfmap"]] <- unify_colmap_nfmap


#' nf-casereport Parser
#'
#' Parses outputs of the nf-casereport directory structure
#'
#' This is a function that extracts paths from a
#' pipeline. The default structure is that the
#' path patterns expected follow that of nf-casereport run.
#' The columns it extracts are defined by
#' variable Skilift::nf_column_map.
#' You can output the default template to modify it using
#' base::dput(Skilift::nf_path_patterns)
#' @export
parse_pipeline_paths <- function(
    paths,
    initial_dt = data.table(pair = character(0)),
    path_patterns = Skilift::nf_path_patterns,
    id_name = "pair") {
  if (length(paths) == 1 && dir.exists(paths)) {
    paths <- list.files(
      paths,
      recursive = TRUE,
      full.names = TRUE
    )
  }

  # Get the list of pairs we're looking for
  pairs_to_match <- unique(initial_dt[[id_name]])
  # Map of regex patterns to column names
  for (col_name in names(path_patterns)) {
    
    patterns = path_patterns[[col_name]]
    for (pattern in patterns) {
      present_paths <- grep(pattern, paths, value = TRUE, perl = TRUE)
      if (!length(present_paths) > 0) next

      # Create mapping between paths and pairs
      dt <- data.table()
      for (pair in pairs_to_match) {
        # Find paths that contain this pair name
        # FIXME: hack for snpeff and other paths that have "snpeff/somatic/PAIR-lane_X/..." directory structure
        # Seems not to be localized to just snpeff
        # balanced_jabba_gg, cbs, jabba_gg, somatic_variant_annotations, events, allelic_jabba_gg
          pair_paths <- grep(paste0("/", pair, "(-lane_.*)?", "/?"), present_paths, value = TRUE, perl = TRUE)
          pair_paths = pair_paths[which.max(file.mtime(pair_paths))]
        # if (col_name %in% c("somatic_variant_annotations", "jabba_gg")) {
        #   pair_paths <- grep(paste0("/", pair, "(-lane.*)?", "/"), present_paths, value = TRUE, perl = TRUE)
        # } else {
        #   pair_paths <- grep(paste0("/", pair, "/"), present_paths, value = TRUE, perl = TRUE)
        # }
        
        if (length(pair_paths) > 0) {
          dt <- rbindlist(list(dt, data.table(
            pair = pair,
            path = pair_paths
          )))
        }
      }

      if (nrow(dt) > 0) {
        setnames(dt, "path", col_name)
        initial_dt <- merge.data.table(
          initial_dt,
          dt,
          by = "pair",
          all.x = TRUE,
          all.y = FALSE, # Only keep pairs that were in initial_dt
          allow.cartesian = FALSE,
		  suffixes = c("", "___THROWOUT")
        )
		# initial_dt = base::subset(
		# 	initial_dt,
		# 	select = !grepl("___THROWOUT$", names(initial_dt))
		# )
        break
      }
    }
  }

  return(initial_dt)
}


#' deparse1 copy
#'
#' Robustly deparse object
#'
deparse1 <- function(expr, collapse = " ", width.cutoff = 500L, ...) {
  paste(deparse(expr, width.cutoff, ...), collapse = collapse)
}

cohort_attributes <- c(
  "inputs",
  "reference_name",
  "type",
  "nextflow_results_path"
)


#' Subset Cohort object
#'
#' Overloads subset operator for Cohort
#'
#' @export
"[.Cohort" <- function(obj, i = NULL, j = NULL, with = TRUE, ...) {
  expri <- deparse1(substitute(i))
  is_i_given <- any(expri != "NULL")
  vector_of_column_names <- deparse1(substitute(j))
  is_j_given <- any(vector_of_column_names != "NULL")
  tbl <- data.table::copy(obj$inputs)
  tblj <- tbl
  if (is_j_given) {
    if (any(grepl("\"|\'", vector_of_column_names))) vector_of_column_names <- j
    selectj <- unique(c("pair", vector_of_column_names))
    tblj <- base::subset(tblj, select = selectj)
  }
  colmap <- as.list(names(tblj))
  names(colmap) <- names(tblj)
  tbli <- tblj[]
  if (is_i_given) {
      expr = substitute(i)
      if (is.call(expr)) {
          tbli <- tblj[eval(expr), , with = with][]
      } else {
          tbli <- tblj[expr, , with = with][]
      }
      
  }
  # obj_out$inputs = tbli
  suppressWarnings({
    obj_out <- Skilift::Cohort$new(
      x = data.table()
    )
  })
  attributes_to_copy <- Skilift:::cohort_attributes
  attributes_to_copy <- attributes_to_copy[!attributes_to_copy %in% "inputs"]
  for (attribute in attributes_to_copy) {
    obj_out[[attribute]] <- obj[[attribute]]
  }
  obj_out$inputs <- tbli[]
  invisible(obj_out$inputs[])
  return(obj_out)
}


#' Parse .nextflow.log
#'
#' Parse nextflow log files for launchdir and input options
#' to get samplesheet
#'
#' @author Kevin Hadi
read_nf_log <- function(nflog_path, max_lines = 1000000L) {
  f <- file(nflog_path, open = "r")
  on.exit({
    close(f)
  })
  line <- readLines(f, n = 1L)
  is_line_at_options <- grepl("Input/output options", line)
  is_line_at_launchdir <- grepl("launchDir", line)
  are_all_lines_parsed <- is_line_at_options && is_line_at_launchdir
  counter <- 1
  ansi_pattern <- "\\\033\\[((?:\\d|;)*)([a-zA-Z])"
  input_samplesheet <- ""
  launchdir <- ""
  while (!are_all_lines_parsed && counter <= max_lines) {
    line <- readLines(f, n = 1L)
    if (!identical(is_line_at_options, TRUE)) {
      is_line_at_options <- length(line) && grepl("Input/output options", line)
      if (identical(is_line_at_options, TRUE)) {
        line <- readLines(f, n = 1L)
        input_samplesheet <- gsub(ansi_pattern, "", line)
        input_samplesheet <- gsub("input[[:space:]]+:[[:space:]]+", "", trimws(input_samplesheet))
      }
    }
    if (!identical(is_line_at_launchdir, TRUE)) {
      is_line_at_launchdir <- length(line) && grepl("launchDir", line)
      if (identical(is_line_at_launchdir, TRUE)) {
        launchdir_line <- gsub(ansi_pattern, "", line)
        launchdir <- gsub("launchDir[[:space:]]+:[[:space:]]+", "", trimws(launchdir_line))
      }
    }
    are_all_lines_parsed <- is_line_at_options && is_line_at_launchdir
    counter <- counter + 1
  }
  if (counter == max_lines && !are_all_lines_parsed) {
    stop("Could not find Input/output options!")
  }
  return(
    list(
      samplesheet = input_samplesheet,
      launchdir = launchdir
    )
  )
}


## Assign new methods
Cohort$private_methods[["read_nf_log"]] <- read_nf_log


## Already - we need to deal with backwards compatibility
## It makes sense to refer to Cohort$type and in the arguments to functions, specify "cohort_type"
## But there are cached cohorts where we've already instantiated the cohort_type attribute.
## Skilift::refresh() can restore this, but we would need to make sure that the right attributes are
## propagated forward.
refresh_attributes = list(
  list("type", "cohort_type") ## attributes that could be present, but need to be updated
)

#' Refresh Cohort object
#'
#' Reinstantiate Cohort object
#'
#' @export
refresh_cohort <- function(cohort) {
  former_inputs = data.table::copy(cohort$inputs) # Create a deep copy of the inputs
  former_inputs_colnames = names(former_inputs)
  obj_out <- Skilift::Cohort$new(
    x = former_inputs 
  )
  for (attribute in setdiff(Skilift:::cohort_attributes, "inputs")) {
    obj_out[[attribute]] <- cohort[[attribute]]
  }
  for (attribute_lst in Skilift:::refresh_attributes) {
    if ( !is.null(cohort[[ attribute_lst[[2]] ]]) ) {
      obj_out[[ attribute_lst[[1]] ]] <- cohort[[ attribute_lst[[2]] ]]
    }
  }
  for (remaining_col in former_inputs_colnames[! former_inputs_colnames %in% names(obj_out$inputs)]) {
    obj_out$inputs[[remaining_col]] = former_inputs[[remaining_col]]
  }
  return(obj_out)
}


#' pairify_cohort_inputs
#' 
#' Create paired cohort inputs if tumor and normal provided
#' for paired outputs
#' 
#' @return paired table from cohort inputs
#' @export
pairify_cohort_inputs = function(cohort, tumor_status = 1L, normal_status = 0L, tumor_normal_columns, sep_cast = "___", keep_remaining = FALSE, ...) {
  nextflow_results_path = NULL
  if (R6::is.R6(cohort) && !is.null(cohort$inputs)) {
    return_inputs = inputs = data.table::copy(cohort$inputs)
    nextflow_results_path = cohort$nextflow_results_path
  } else if (inherits(cohort, "data.table")) {
    return_inputs = inputs = data.table::copy(cohort)
  } else {
    stop("Cohort object or inputs data.table must be provided")
  }

  default_tumor_type_cols = Skilift::default_col_mapping[["tumor_type"]]

  tumor_type_vec = character(NROW(inputs))
  is_tumor_type_found = FALSE
  for (i in seq_along(default_tumor_type_cols)) {
    tumor_type_col = default_tumor_type_cols[i]
    tumor_type_vals = inputs[[tumor_type_col]]
    is_character = is.character(tumor_type_vals)
    is_numeric = is.integer(tumor_type_vals) || is.numeric(tumor_type_vals)
    is_0or1 = is_numeric && all(tumor_type_vals %in% c(0, 1))  
    if (is_character) {
      tumor_type_vec[
        grepl("tumor", tumor_type_vals, ignore.case = TRUE)
      ] = "tumor"
      tumor_type_vec[
        grepl("normal", tumor_type_vals, ignore.case = TRUE)
      ] = "normal"
    } else if (is_0or1) {
      tumor_type_vec[
        tumor_type_vals == tumor_status
      ] = "tumor"
      tumor_type_vec[
        tumor_type_vals == normal_status
      ] = "normal"
    }
    if (all(tumor_type_vec %in% c("tumor", "normal"))) is_tumor_type_found = TRUE
    if (is_tumor_type_found) {
      break
    }
  }

  # if (!is_tumor_type_found) stop("Tumor type was not found!")
  inputs$tumor_type = tumor_type_vec

	# is_status_in_cohort = !is.null(inputs$status)

	is_nextflow_results_path_present = !is.null(nextflow_results_path)
	# is_tumor_bam_in_cohort = !is.null(inputs$tumor_bam)
	# is_normal_bam_in_cohort = !is.null(inputs$normal_bam)

	# is_unpaired = is_status_in_cohort && is_nextflow_results_path_present
	# is_paired = is_tumor_bam_in_cohort && is_normal_bam_in_cohort

  is_unpaired = is_tumor_type_found ## TODO: test if this is sufficient.

  tumor_normal_columns = tumor_normal_columns[tumor_normal_columns %in% names(inputs)]

	if (is_unpaired && NROW(tumor_normal_columns) > 0) {
		# inputs$tumor_type = ifelse(inputs$status == tumor_status, "tumor", "normal")
		return_inputs = Skilift::dcastski(
			inputs,
			id_columns = "pair",
			type_columns = "tumor_type",
			cast_columns = tumor_normal_columns,
			sep_cast = sep_cast,
			keep_remaining = keep_remaining,
      ...
		)
	}
	return(return_inputs)
}
