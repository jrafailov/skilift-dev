
#' Print useful messages to console
#' 
#' For skilift lift-wrappers.R.
#' 
#' @export
skimessage = function(..., pre = 'Skilift'){
    message(pre, ' ', paste0(as.character(Sys.time()), ': '), ...)
}

#' Silence outputs.
#' 
#' Capture error messages
#' but silence superfluous intermediate messages
#' printed from lift_all().
#' 
#' @export
shutup = function(top_level_expr, capture_output_type = "message") {
  nr = NROW(capture_output_type)
  do_capture_output = identical("output", capture_output_type)
  do_capture_message =identical("message", capture_output_type)
  do_capture_both = all(c("message", "output") %in% capture_output_type)
  tryCatch(
    expr = {
      if (do_capture_output)
        capture.output(... = top_level_expr, type = "output")
      
      if (do_capture_message)
        capture.output(... = top_level_expr, type = "message")
      
      if (do_capture_both) {
        capture.output(
          ... = capture.output(... = top_level_expr, type = "message"),
          type = "output"
        )
      }
    },
    error = function(e) {
      stop(e)
    }
  )
}


#' Cast table
#' 
#' Cast using base R.
#' 
#' Using base R for robustness and flexibility.
#' 
#' @export
dcastski = function(
	tbl, 
	id_columns, 
	type_columns, 
	cast_columns, 
	keep_remaining = FALSE, 
	sep_cast = "___", 
	prefix_type = TRUE,
	drop = FALSE,
	use_regex = FALSE
) {
  if (identical(use_regex, TRUE)) {
	.NotYetImplemented()
  }
  if (base::anyDuplicated(names(tbl)) > 0) {
    stop("Duplicated names present in table! Dedup first: names(tbl) = base::make.unique(names(tbl))")
  }
  remaining_cols = character(0)
  columns_to_process = unique(c(id_columns, type_columns))
  if (identical(keep_remaining, TRUE)) {
    remaining_cols = names(tbl)[!names(tbl) %in% c(columns_to_process, cast_columns)]
    columns_to_process = c(columns_to_process, remaining_cols)
  }
  columns_to_process = c(columns_to_process, cast_columns)
  ## A very unfortunate circumstance.
  ## Have to use unique.data.frame due to lists present in cohort inputs
  ## The data.frame base R method is slow, but works on list elements in columns.
  ## If tbl is a data.table, skeleton will remain a data.table
  skeleton = unique.data.frame(
    base::subset(tbl, select = names(tbl) %in% c(id_columns, remaining_cols))
  )
  reduced_tbl = base::subset(tbl, select = names(tbl) %in% columns_to_process)
  reduced_tbl$types___reduced___l9vuK2F6lw = reduced_tbl[[type_columns[1]]]
  types = unique(reduced_tbl$types___reduced___l9vuK2F6lw)
  if (is.factor(reduced_tbl$types___reduced___l9vuK2F6lw) && !identical(drop, TRUE)) {
    types = base::levels(reduced_tbl$types___reduced___l9vuK2F6lw)
  }
  if (length(type_columns) > 1) {
    .NotYetImplemented()
    lst = as.list(base::subset(reduced_tbl, select = names(reduced_tbl) %in% type_columns))
    types = do.call(interaction, lst)
    skeleton$types___reduced___l9vuK2F6lw = types
  }
  for (type in types) {
    type___l9vuK2F6lw = type
    merge_type = reduced_tbl[reduced_tbl$types___reduced___l9vuK2F6lw == type___l9vuK2F6lw,]
    merge_type = base::subset(merge_type, select = names(merge_type) %in% c(id_columns, cast_columns))
    colnms = names(merge_type)
    names_to_change = colnms[colnms %in% cast_columns]
	if (prefix_type) {
		names_to_change = paste(type, names_to_change, sep = sep_cast)
	} else {
		names_to_change = paste(names_to_change, type, sep = sep_cast)
	}
    
    colnms[colnms %in% cast_columns] = names_to_change
    names(merge_type) = colnms
    skeleton = merge(skeleton, merge_type, by = id_columns, all.x = TRUE)
  }
  return(skeleton)
}


#' Melt table
#' 
#' Melt using base R.
#' 
#' Using base R for robustness and flexibility, not speed necessarily.
#' Since this will be used primarily for sample/pairs table manipulation.
#' 
#' @export
meltski = function(
    tbl,
    id.vars,
    measure.vars,
    is_measure_regex = FALSE,
    variable.name = "variable",
    value.name = "value",
    group_regex = "",
    group_select = "\\1",
    replacement_regex = "",
    replacement = "",
    drop = FALSE,
    keep_first_variable_col = FALSE,
    keep_remaining_cols = FALSE,
    keep_original_cols = FALSE,
    return_as_data_table = FALSE
) {
  keep_first_variable_col = identical(keep_first_variable_col, TRUE)
  keep_remaining_cols = identical(keep_remaining_cols, TRUE)
  keep_original_cols = identical(keep_original_cols, TRUE)
  
  is_table = is.table(tbl)
  
  dimensions = dim(tbl)
  is_dim_null = is.null(dimensions)
  dim_NR = NROW(dimensions)
  is_table_matrix_like = is_table && !is_dim_null && dim_NR == 2
  is_table_1d = is_table && !is_dim_null && dim_NR == 1
  is_table_array = is_table && !is_dim_null && dim_NR > 2
  if (is_table_array) stop("melting more than 2d not supported")
  if (is_table_1d) {
    tbl = as.data.frame.matrix(t(as.matrix(tbl)))
  }
  if (is_table_matrix_like) {
    tbl = as.data.frame.matrix(tbl)
  }

  if (is_table) measure.vars = colnames(tbl)

  is_data_frame = inherits(tbl, "data.frame")

  if (!is_data_frame) stop("Input tbl must be a data.table/data.frame or coercible to a data.frame-like object")

  all_names = names(tbl)

  is_measure_vars_list = is.list(measure.vars)
  all_measure_vars = measure.vars
  if (is_measure_vars_list) all_measure_vars = unlist(measure.vars)
    
  if (missing(id.vars)) {
    id.vars = all_names
    id.vars = id.vars[!id.vars %in% all_measure_vars]
  }
  skel = NULL
  collected_var_names = character(0)
  if (is_measure_regex) {
    message("!!WARNING!!: Regex using measure.vars does not guarantee ordered outputs.")
    message("!!WARNING!!: You must ensure that input columns are ordered with respect to the regex that you provide")
    message("!!WARNING!!: Otherwise, you should check that your output value columns are not mixed.")
    message("!!WARNING!!: Attempting to help by sorting column names first.")
    all_names = sort(all_names)
  }
  for (var_i in seq_along(measure.vars)) {
    var = measure.vars[var_i]
    if (is_measure_regex) var = grep(var, all_names, value = TRUE)
    if (is.list(var)) var = unlist(var)
    tbl_to_rbind = base::subset(tbl, select = names(tbl) %in% c(id.vars, var))
    if (drop) {
      is_none_na = base::complete.cases(
        base::subset(tbl_to_rbind, select = names(tbl_to_rbind) %in% var)
      )
      tbl_to_rbind = base::subset(tbl_to_rbind, is_none_na)
    }
    nm = names(tbl_to_rbind)
    for (i in seq_along(var)) {
      len_value_name = NROW(value.name)
      if (!(len_value_name > 1 && len_value_name == length(var))) {
        value_col_name = paste(value.name, "_", i, sep = "")
      } else {
        value_col_name = value.name[i]
      }
      nm[nm %in% var[i]] = value_col_name
    }
    names(tbl_to_rbind) = nm
    for (i in seq_along(var)) {
      var.col = paste(variable.name, "_", i, sep = "")
      tbl_to_rbind[[var.col]] = var[i]
      collected_var_names = c(collected_var_names, var.col)
    }    
    skel = data.table:::rbind.data.table(skel, tbl_to_rbind, fill = TRUE)
  }
  collected_var_names = unique(collected_var_names)
  if (nzchar(group_regex)) {
    for (variable in collected_var_names) {
      skel[[variable]] = gsub(group_regex, group_select, skel[[variable]])
    }
  }
  if (nzchar(replacement_regex)) {
    for (variable in collected_var_names) {
      skel[[variable]] = gsub(replacement_regex, replacement, skel[[variable]])
    }    
  }
  remaining_var_names = collected_var_names[-1]
  remaining_names = names(skel)[!names(skel) %in% remaining_var_names]
  if (keep_first_variable_col && NROW(remaining_names) > 0) {
    skel = base::subset(skel, select = names(skel) %in% remaining_names)
    names(skel)[names(skel) == collected_var_names[1]] = gsub("_[0-9]+$", "", collected_var_names[1])
  }
  is_keep_remaining_flag_on = keep_remaining_cols || keep_original_cols
  
  if (is_keep_remaining_flag_on) {
    is_id_var = all_names %in% id.vars
    is_selection = rep_len(TRUE, length(all_names))
    if (!keep_original_cols) {
      is_selection = all_names %in% id.vars | (!all_names %in% c(all_measure_vars, value.name))
    }      
    remaining_tbl = base::subset(tbl, select = is_selection)
    skel = merge(skel, remaining_tbl, by = id.vars, all.x = TRUE)
  }
  return(skel)
}


#' Merge Cohort objects
#'
#' Combines two or more Cohort objects into a single Cohort
#'
#' @param x Cohort object to merge
#' @param y Cohort object to merge
#' @param ... Arbitrary remaining Cohort objects to merge
#' @param warn_duplicates Logical indicating whether to warn about duplicate pairs (default TRUE)
#' @param rename_duplicates Logical indicating whether to rename duplicate pairs instead of overwriting (default FALSE)
#' @return A new Cohort object containing all data from input Cohorts
#' 
#' @export merge.Cohort 
#' @export
merge.Cohort <- function(
	x, 
	y, 
	..., 
	prefer_x = TRUE, 
	prefer_y = FALSE, 
	prefer_y_na = FALSE, 
	warn_duplicates = TRUE, 
	rename_duplicates = FALSE
) {
  cohorts_dots <- list(...)
  is_x_missing = missing(x)
  is_y_missing = missing(y)
  is_x_or_y_missing = is_x_missing || is_y_missing
  number_cohorts_dots = NROW(cohorts_dots)
  length_cohorts = sum(as.integer(!is_x_missing), as.integer(!is_y_missing), number_cohorts_dots)
  # Validate inputs are all Cohorts
  if (length_cohorts < 2) {
    stop("At least two Cohort objects must be provided")
  }
  if (is_x_or_y_missing) {
    if (is_x_missing) print("x argument not provided!")
    if (is_y_missing) print("y argument not provided!")
    stop("Both x and y must be provided as Cohort objects")
  }
  cohorts = c(list(x), list(y), cohorts_dots)
  for (cohort in cohorts) {
    if (!inherits(cohort, "Cohort")) {
      stop("All arguments must be Cohort objects")
    }
  }

  cohorts_type = unlist(lapply(cohorts, function(cohort) cohort$type))
  is_cohort_type_same_in_all = all(cohorts_type == cohorts_type[1])
  if (!is_cohort_type_same_in_all) {
	stop(
		"Cohorts to be merged have different types!\n", 
		"Cohort type is defined at cohort object level, not at a sample level, currently.\n"
	)
  }

  # Get first cohort to initialize merged result
  merged_dt <- data.table::copy(cohorts[[1]]$inputs)

  # Merge remaining cohorts
  for (i in 2:length(cohorts)) {
    current_dt <- data.table::copy(cohorts[[i]]$inputs)

    # Check for duplicate pairs
    duplicate_pairs <- intersect(merged_dt$pair, current_dt$pair)
    if (length(duplicate_pairs) > 0) {
      if (warn_duplicates) {
        message(sprintf(
          "Found %d duplicate pair(s): %s",
          length(duplicate_pairs),
          paste(duplicate_pairs, collapse = ", ")
        ))
      }

      if (rename_duplicates) {
        # Add incrementing suffix to duplicate pairs
        for (dup_pair in duplicate_pairs) {
          suffix <- 1
          while (paste0(dup_pair, "_", suffix) %in% merged_dt$pair) {
            suffix <- suffix + 1
          }
          current_dt[pair == dup_pair, pair := paste0(pair, "_", suffix)]
        }
	  }
    }

    # Merge inputs data.tables
    merged_dt <- Skilift::merge.repl(
      merged_dt, 
      current_dt,
      by="pair", 
      all = TRUE, ## Need these two on to do outer join, which should be expected behavior, you don't want to lose anything,
	  prefer_x = prefer_x,
	  prefer_y = prefer_y,
	  prefer_y_na = prefer_y_na
    )
  }

  # Create new Cohort with merged data
  result = Skilift::copy(cohorts[[1]])
  result$inputs = merged_dt

  return(result)
}

#' Test path
#' 
#' Test if path exists robustly
#' 
#' @export
test_path = function(
  object, 
  rds_regex = ".rds$",
  gff_regex = ".gtf(.gz){0,}$|.gff([0-9]){0,}(.gz){0,}$",
  bcf_regex = ".bcf(.bgz|.gz){0,}$",
  vcf_regex = ".vcf(.bgz|.gz){0,}$",
  verbose = TRUE
) {
  is_null = is.null(object)
  is_character = is.character(object)
  is_len_one = NROW(object) == 1
  is_not_valid = is_character && ! NROW(object) == 1
  is_na = is_len_one && all(Skilift::is_loosely_na(object, other_nas = base::nullfile()))
  is_possible_path = is_character && is_len_one && !is_na
  is_existent_path = is_possible_path && file.exists(object)
  is_rds = is_possible_path && grepl(rds_regex, object)
  is_vcf = is_possible_path && grepl(vcf_regex, object)
  is_bcf = is_possible_path && grepl(bcf_regex, object)
  is_gff = is_possible_path && grepl(gff_regex, object)
  logicals = as.list(data.frame(
      is_null,
      is_character,
      is_len_one,
      is_not_valid,
      is_na,
      is_possible_path,
      is_existent_path,
      is_rds,
      is_vcf,
      is_bcf,
      is_gff
  ))
  if (verbose) {
    nms = names(logicals)
    message("Assigning:")
    for (nm in nms) {
      message(nm)
    }
  }
  list2env(logicals, envir = parent.frame())
  return(logicals)
}

#' Test path vector
#' 
#' Test if paths exist robustly
#' 
#' @export
test_paths = function(
  objects, 
  rds_regex = ".rds$",
  gff_regex = ".gtf(.gz){0,}$|.gff([0-9]){0,}(.gz){0,}$",
  bcf_regex = ".bcf(.bgz|.gz){0,}$",
  vcf_regex = ".vcf(.bgz|.gz){0,}$",
  verbose = TRUE
) {
  is_null = is.null(objects)
  is_character = is.character(objects)
  is_na = Skilift::is_loosely_na(objects, other_nas = base::nullfile())
  is_possible_path = is_character & !is_na
  if (!is_character) {
	objects = as.character(objects)
  }
  is_robust_file_exists = tryCatch(file.exists(objects), error = function(e) rep_len(FALSE, NROW(objects)))
  is_existent_path = is_possible_path & is_robust_file_exists
  is_rds = is_possible_path & grepl(rds_regex, objects)
  is_vcf = is_possible_path & grepl(vcf_regex, objects)
  is_bcf = is_possible_path & grepl(bcf_regex, objects)
  is_gff = is_possible_path & grepl(gff_regex, objects)
  logicals = list(
      is_null = is_null,
      is_character = is_character,
      is_na = is_na,
      is_possible_path = is_possible_path,
      is_existent_path = is_existent_path,
      is_rds = is_rds,
      is_vcf = is_vcf,
      is_bcf = is_bcf,
      is_gff = is_gff
  )
  if (verbose) {
    nms = names(logicals)
    message("Assigning:")
    for (nm in nms) {
      message(nm)
    }
  }
  list2env(logicals, envir = parent.frame())
  return(logicals)
}

#' Read jabba
#' 
#' Function to flexibly read in jabba files whether providing raw jabba list or ggraph
#' 
#' @export
process_jabba = function(jabba, allow_return_null = FALSE) {
  logicals = test_path(jabba, verbose = FALSE)
  do_stop = FALSE
  if (is_existent_path && is_rds) {
    jabba <- readRDS(jabba)
  } else if (is_existent_path) {
    msg = sprintf("Path exists but is not rds: %s", jabba)
    do_stop = TRUE
  } else if (is_character && (is_len_one || is_na)) {
    msg = sprintf("Path provided does not exist %s", jabba)
    do_stop = TRUE
  } else if (is_not_valid) {
    msg = sprintf("Path provided must be a length one string %s", jabba)
    do_stop = TRUE
  } else if (is_null) {
    msg = sprintf("Path provided is NULL")
    do_stop = TRUE
  }

  do_return_null = do_stop && allow_return_null
  do_error = do_stop && !allow_return_null
  if (do_return_null) return(NULL)
  if (do_error) stop(msg)

  is_jabba_list = is.list(jabba) && all(c("segstats", "purity", "ploidy", "junctions") %in% names(jabba))
  
  gg = jabba
  if (is_jabba_list) gg = gG(jabba = jabba)
  is_jabba_gg = inherits(gg, "gGraph")
  if (!is_jabba_gg) {
    stop("jabba must be a gGraph object or jabba like object")
  }
  return(gg)
}


#' robust name()
#'
#' gives back character vector same length of input regardless whether named or not
#'
#' @param str a path string
#' @return a string with multiple parentheses replaced with a single parenthesis
#' @export
names2 = function(x) {
    nm = names(x)
    if (is.null(nm))
        return(rep("", length.out = length(x)))
    else
        return(nm)
}

#' convert columns with NA to false
#'
#' coerce NA in columns of class "logical" to FALSE
#'
#' @param dt data.table
#' @param these_cols NULL by default, will select columns of class logical, otherwise will be specified
#' @return A data.table
#' @export
dt_na2false = function(dt, these_cols = NULL) {
    na2false = function(v)
    {
        ## v = ifelse(is.na(v), v, FALSE)
        v[is.na(v)] = FALSE
        as.logical(v)
    }
    if (is.null(these_cols)) {
        these_cols = which(sapply(dt, class) == "logical")
    }
    for (this_col in these_cols) {
        ## this_val = as.data.frame(dt[, this_col, with = FALSE])[,1]
        this_val = dt[[this_col]]
        data.table::set(dt, j = this_col, value = na2false(this_val))
    }
    return(dt)
}

#' Kind of NA
#' 
#' Test for NAs of different types
#' 
#' NA may not be the only type that we want to test for
#' 
#' If it's a character, we may want to loosely test for other "NAs"
#' @export
is_loosely_na = function(values, character_nas = c("NA", "na", "N/A", "n/a", "NULL", "Null", "null"), other_nas = NULL) {
	is_proper_na = is.na(values)
	is_values_character = is.character(values)
	is_empty = is_other_na_type = rep_len(FALSE, NROW(values))
	if (is.null(other_nas) || NROW(other_nas) == 0 || all(is.na(other_nas))) other_nas = character(0)
	if (!is.character(other_nas)) other_nas = as.character(other_nas)
	if (is_values_character) {
		is_other_na_type = tolower(values) %in% c(character_nas, other_nas)
		is_empty = nchar(values) == 0
	}
	is_na = is_proper_na | is_other_na_type | is_empty
	return(is_na)
}


#' merging data tables with collapsing columns with the same name
#'
#' Merge two data tables with various replacing strategies
#' for columns common between x and y that are not used to merge
#' (i.e. not specified in the "by" argument)
#'
#' @param replace_NA logical, only use values in dt.y, any dt.x not in dt.y is clobbered (NA)
#' @param delete_xiny logical, delete columns in x that intersect with y columns (remove previous columns)
#' @param force_y logical, should x and y common columns be merged?
#' @param prefer_y logical, for x and y entries sharing common by keys, prefer y (logic to control whether to ignore na or not encoded by prefer_y_na)
#' @param prefer_y_na logical, for x and y entries sharing common by keys, prefer y and ignore NAs
#' @param overwrite_x logical, if force_y = TRUE, should NA values in y replace x?
#' @return A data.table
#' @export merge.repl
merge.repl = function(dt.x,
                      dt.y,
                      sep = "_",
                      replace_NA = TRUE,
					  delete_xiny,
                      force_y = TRUE,
					  prefer_y,
					  prefer_x,
                      overwrite_x = FALSE,
					  prefer_y_na,
                      keep_order = FALSE,
                      keep_colorder = TRUE,
                      keep_factor = TRUE,
					  suffixes = c(".x", ".y"),
                      ...) {
    
    arg_lst = as.list(match.call())
    by.y = eval(arg_lst[['by.y']], parent.frame())
    by.x = eval(arg_lst[['by.x']], parent.frame())
    by = eval(arg_lst[['by']], parent.frame())
    all.x = eval(arg_lst[['all.x']], parent.frame())
    all.y = eval(arg_lst[['all.y']], parent.frame())
    all = eval(arg_lst[['all']], parent.frame())
    allow.cartesian = eval(arg_lst[['allow.cartesian']])
	
    key_x = key(dt.x)
    is_all.x_not_provided = is.null(all.x)
    is_all.y_not_provided = is.null(all.y)
    is_all.y_provided = !is.null(all.y)
    is_all_provided = !is.null(all)
    if (is_all.x_not_provided && is_all.y_not_provided) {
        all.x = TRUE
        all.y = FALSE
    }
    if (is_all.x_not_provided && is_all.y_provided) {
        all.x = FALSE
    }
    if (is_all_provided && identical(all, TRUE)) {
        all.y = TRUE
        all.x = TRUE
    }
    if (is.null(allow.cartesian)) {
        allow.cartesian = FALSE
    }
    if (!inherits(dt.x, "data.table")) {
        dt.x = as.data.table(dt.x)
    }
    if (!inherits(dt.y, "data.table")) {
        dt.y = as.data.table(dt.y)
    }
    if (keep_order == TRUE) {
        dt.x[['tmp.2345098712340987']] = seq_len(nrow(dt.x))
    }

	call_replace_NA = eval(arg_lst[['replace_NA']], parent.frame())
	was_replace_NA_provided = !is.null(call_replace_NA)


	call_delete_xiny = eval(arg_lst[['delete_xiny']], parent.frame())
	was_delete_xiny_provided = !is.null(call_delete_xiny)

	if (was_replace_NA_provided || !was_delete_xiny_provided) {
		message("replace_NA is deprecated, please provide delete_xiny instead next time..")
		message("using !replace_NA value provided: ", !replace_NA)
		delete_xiny = !replace_NA
	}


	call_force_y = eval(arg_lst[['force_y']], parent.frame())
	# call_replace_na_x = eval(arg_lst[['replace_na_x']], parent.frame())
	was_force_y_provided = !is.null(call_force_y)

    call_prefer_y = eval(arg_lst[['prefer_y']], parent.frame())
	was_prefer_y_provided = !is.null(call_prefer_y)

    
	if (was_force_y_provided || !was_prefer_y_provided) {
		message("force_y is deprecated, please provide prefer_y instead next time..")
		message("using force_y value provided: ", force_y)
		prefer_y = force_y
	}


	call_overwrite_x = eval(arg_lst[['overwrite_x']], parent.frame())
	was_overwrite_x_provided = !is.null(call_overwrite_x)

    call_prefer_y_na = eval(arg_lst[['prefer_y_na']], parent.frame())
	was_prefer_y_na_provided = !is.null(call_prefer_y_na)

    
	if (was_overwrite_x_provided || !was_prefer_y_na_provided) {
		message("overwrite_x is deprecated, please provide prefer_y_na instead next time..")
		message("using prefer_y_na value provided: ", overwrite_x)
		prefer_y_na = overwrite_x
	}

	call_prefer_x = eval(arg_lst[['prefer_x']], parent.frame())
	was_prefer_x_provided = !is.null(call_prefer_x)


    if (was_prefer_x_provided && identical(prefer_x, prefer_y)) {
        prefer_y = !prefer_x
        prefer_y_na = !prefer_x
    }

  # was_prefer_y_and_not_prefer_x_provided = was_prefer_y_provided && !was_prefer_x_provided
  # if (was_prefer_y_provided && !was_prefer_x_provided) {
  #   prefer_x = !prefer_y
  # }

  if (!was_prefer_x_provided) {
    prefer_x = !prefer_y
  }

  if (was_prefer_y_na_provided && !was_prefer_x_provided) {
    prefer_x = !prefer_y_na
  }

	is_preference_invalid = identical(prefer_y, prefer_x)
	if (is_preference_invalid) {
		stop("prefer_y and prefer_x are mutually exlusive")
	}


	is_only_delete_provided_and_true = was_delete_xiny_provided && !was_prefer_y_provided && identical(delete_xiny, TRUE)

	is_delete_and_prefer_choice_invalid = (
		identical(delete_xiny, TRUE) 
		&& (
			identical(delete_xiny, prefer_y) 
			|| identical(delete_xiny, prefer_y_na)
		)
	)
	if ( ! is_only_delete_provided_and_true && is_delete_and_prefer_choice_invalid) {
		message("delete_xiny: ", delete_xiny)
		message("prefer_y: ", prefer_y)
		message("prefer_y_na: ", prefer_y_na)
		stop("delete_x must be mutually exclusive with both prefer_y and prefer_y_na")
	}

	is_prefer_y_na_invalid = identical(prefer_y_na, TRUE) && !identical(prefer_y_na, prefer_y)
	if (is_prefer_y_na_invalid) {
		message("prefer_y_na and prefer_y args are incompatible")
		message("If prefer_y_na is set to TRUE, specify explicitly that you would prefer y values for joined columns by setting prefer_y to TRUE also")
		stop("prefer_y_na set to FALSE, but prefer_y set to TRUE")
	}
	

	is_suffixes_unique = identical(anyDuplicated(suffixes), 0L)
	is_suffixes_all_nonempty_string = is.character(suffixes) && all(nzchar(suffixes))
	is_suffixes_len_2 = NROW(suffixes) == 2
	if (!is_suffixes_len_2 || !is_suffixes_unique || !is_suffixes_all_nonempty_string) {
		stop("suffixes argument must be provided as two non empty, unique strings, e.g. suffixes = c('.x', '.y')")
	}

    dt.x[['in.x.2345098712340987']] = rep(TRUE, length.out = nrow(dt.x))
    dt.y[['in.y.2345098712340987']] = rep(TRUE, length.out = nrow(dt.y))

    new_ddd_args = list(
      by = by, by.x = by.x, by.y = by.y,
      all.x = all.x, all.y = all.y,
      allow.cartesian = allow.cartesian,
	  suffixes = suffixes
    )

    if (is.null(by.y) && is.null(by.x) && is.null(by)) {

        if (
			length(attributes(dt.x)[['sorted']]) > 0 
			&& length(attributes(dt.y)[['sorted']]) > 0
		) {
            k.x = data.table::key(dt.x)
            k.y = data.table::key(dt.y)
        } else {
            k.y = k.x = intersect(names2(dt.x), names2(dt.y))
            if (length(k.x) == 0)
                stop("no common columns to merge by!")
            message("intersecting by: ", paste(k.x, collapse = ", "))
            new_ddd_args[['by']] = k.x
        }
        if ((!is.null(k.x) && !is.null(k.y)) && !identical(k.x, k.y)) {
          stop(
            "neither by.x/by.y nor by are supplied, ",
            "keys of dt.x and dt.y ",
            "must be identical and non NULL"
          )
        }
        x.cols = setdiff(names(dt.x), k.x)
        y.cols = setdiff(names(dt.y), k.y)

    } else if (!is.null(by.x) && !is.null(by.y)) {

        x.cols = setdiff(names(dt.x), by.x)
        y.cols = setdiff(names(dt.y), by.y)
        new_ddd_args = new_ddd_args[setdiff(names(new_ddd_args), c("by"))]

    } else if (!is.null(by)) {

      x.cols = setdiff(names(dt.x), by)
      y.cols = setdiff(names(dt.y), by)
      if (! all(by %in% colnames(dt.x)) || ! all(by %in% colnames(dt.y))) {
        stop(
          "column ",
          by,
          " does not exist in one of the tables supplied",
          "\nCheck the column names"
        )
      }
      new_ddd_args = new_ddd_args[setdiff(names(new_ddd_args), c("by.y", "by.x"))]

    }
    intersecting_colnames = intersect(x.cols, y.cols)
    ## if (replace_in_x) {
    # if (!replace_NA) {
	dt.x.tomerge = dt.x
	if (delete_xiny) {
		dt.x.tomerge = data.table::copy(dt.x.tomerge)
		for (this_col in intersecting_colnames) {
			data.table::set(dt.x.tomerge, i = NULL, j = this_col, value = NULL)
		}
	}
	dt.repl = suppressWarnings(
		do.call(
			"merge",
			args = c(list(x = dt.x.tomerge, y = dt.y), new_ddd_args)
		)
	)
	dt_na2false(dt.repl, c("in.x.2345098712340987", "in.y.2345098712340987"))
	in.x = which(dt.repl[["in.x.2345098712340987"]])
	in.y = which(dt.repl[["in.y.2345098712340987"]])
	this_env = environment()
	
	for (this_col in intersecting_colnames) {
		x_cname = paste0(this_col, suffixes[1])
		y_cname = paste0(this_col, suffixes[2])
		x_col = dt.repl[[x_cname]]
		y_col = dt.repl[[y_cname]]
		is_x_factor = inherits(x_col, "factor")
		is_y_factor = inherits(y_col, "factor")
		is_either_xy_factor = is_x_factor || is_y_factor
		if ( (is_either_xy_factor) && keep_factor) {
			if (!is_x_factor) { x_col = factor(x_col); is_x_factor = TRUE }
			if (!is_y_factor) { y_col = factor(y_col); is_y_factor = TRUE }
		}
		if (is_x_factor && !keep_factor) { x_col = as.character(x_col); is_x_factor = FALSE } 
		if (is_y_factor && !keep_factor) { y_col = as.character(y_col); is_y_factor = FALSE }
		is_either_xy_factor = is_x_factor || is_y_factor
		
		if (prefer_y) {
			## overwrite_x means that y NAs in keys that overlap between the
			## two x and y tables will be preferred 
			if (!prefer_y_na) {
				if (is_either_xy_factor) {
					new_col = factor(y_col, forcats::lvls_union(list(y_col, x_col)))
					new_col[is.na(new_col)] = x_col[is.na(new_col)]
				} else {
					new_col = ifelse(!Skilift::is_loosely_na(y_col), y_col, x_col)
				}
			} else {
				if (is_either_xy_factor) {
					new_col = factor(x_col, forcats::lvls_union(list(y_col, x_col)))
				} else {
					new_col = x_col
				}
				new_col[in.y] = y_col[in.y]
			}
		} else if (prefer_x) {
			## Only take X column if 
			if (is_either_xy_factor) {
				new_col = factor(x_col, forcats::lvls_union(list(x_col, y_col)))
				new_col[Skilift::is_loosely_na(new_col) & !Skilift::is_loosely_na(y_col)] = y_col[is.na(new_col) & !is.na(y_col)]
			} else {
				new_col = ifelse(Skilift::is_loosely_na(x_col) & !Skilift::is_loosely_na(y_col), y_col, x_col)
			}
		} else {
			stop("How did this happen?")
		}
		data.table::set(
			dt.repl,
			j = c(x_cname, y_cname, this_col),
			value = list(NULL, NULL, this_env[["new_col"]])
		)
	}
    if (identical(keep_order, TRUE)) {
        data.table::setorderv(dt.repl, "tmp.2345098712340987")
        dt.repl[['tmp.2345098712340987']] = NULL
    }
    data.table::set(
      dt.repl,
      j = c("in.y.2345098712340987", "in.x.2345098712340987"),
      value = list(NULL, NULL)
    )
    if (keep_colorder) {
        x_cols = colnames(dt.x)
        ## get the order of columns in dt.repl in order of X with
        ## additional columns tacked on end
        data.table::setcolorder(
          dt.repl,
          intersect(
            union(
              colnames(dt.x),
              colnames(dt.repl)),
            colnames(dt.repl)
          )
        )
    }
    return(dt.repl)
}

#' make deep copy, recursively
#'
#' useful for dev
#' makes deep copy of R6 object, S4 object, or anything else really
#'
#' @name copy
#' @export copy
copy = function (x, recurse_list = TRUE) {
    if (inherits(x, "R6")) {
        x2 = rlang::duplicate(x$clone(deep = T))
        for (name in intersect(names(x2$.__enclos_env__), c("private", 
                                                            "public", "self"))) {
            for (nname in names(x2$.__enclos_env__[[name]])) {
                tryCatch({
                    x2$.__enclos_env__[[name]][[nname]] = Skilift::copy(x2$.__enclos_env__[[name]][[nname]])
                }, error = function(e) NULL)
            }
        }
        return(x2)
    } else if (base::isS4(x)) {
        x2 = rlang::duplicate(x)
        slns = slotNames(x2)
        for (sln in slns) {
            tryCatch({
                slot(x2, sln) = Skilift::copy(slot(x2, sln))
            }, error = function(e) NULL)
        }
        return(x2)
    } else if (inherits(x, c("list"))) {
        x2 = rlang::duplicate(x)
        x2 = rapply(x2, Skilift::copy, how = "replace")
        return(x2)
    } else {
        x2 = rlang::duplicate(x)
        return(x2)
    }
}


#' Version of utils::assignInNamespace
#'
#' can be used to reassign function into a namespace
#' USE WITH CAUTION
#'
#' @export
assign_in_namespace = function (x, value, ns, pos = -1, envir = as.environment(pos)) {
    nf <- sys.nframe()
    if (missing(ns)) {
        nm <- attr(envir, "name", exact = TRUE)
        if (is.null(nm) || substr(nm, 1L, 8L) != "package:")
            stop("environment specified is not a package")
        ns <- asNamespace(substring(nm, 9L))
    }
    else ns <- asNamespace(ns)
    ns_name <- getNamespaceName(ns)
    ## if (nf > 1L) {
    ##     if (ns_name %in% tools:::.get_standard_package_names()$base)
    ##         stop("locked binding of ", sQuote(x), " cannot be changed",
    ##             domain = NA)
    ## }
    if (bindingIsLocked(x, ns)) {
        in_load <- Sys.getenv("_R_NS_LOAD_")
        if (nzchar(in_load)) {
            if (in_load != ns_name) {
                msg <- gettextf("changing locked binding for %s in %s whilst loading %s",
                  sQuote(x), sQuote(ns_name), sQuote(in_load))
                if (!in_load %in% c("Matrix", "SparseM"))
                  warning(msg, call. = FALSE, domain = NA, immediate. = TRUE)
            }
        }
        else if (nzchar(Sys.getenv("_R_WARN_ON_LOCKED_BINDINGS_"))) {
            warning(gettextf("changing locked binding for %s in %s",
                sQuote(x), sQuote(ns_name)), call. = FALSE, domain = NA,
                immediate. = TRUE)
        }
        unlockBinding(x, ns)
        assign(x, value, envir = ns, inherits = FALSE)
        w <- options("warn")
        on.exit(options(w))
        options(warn = -1)
        lockBinding(x, ns)
    }
    else {
        assign(x, value, envir = ns, inherits = FALSE)
    }
    if (!isBaseNamespace(ns)) {
        S3 <- .getNamespaceInfo(ns, "S3methods")
        if (!length(S3))
            return(invisible(NULL))
        S3names <- S3[, 3L]
        if (x %in% S3names) {
            i <- match(x, S3names)
            genfun <- get(S3[i, 1L], mode = "function", envir = parent.frame())
            if (.isMethodsDispatchOn() && methods::is(genfun,
                "genericFunction"))
                genfun <- methods::slot(genfun, "default")@methods$ANY
            defenv <- if (typeof(genfun) == "closure")
                environment(genfun)
            else .BaseNamespaceEnv
            S3Table <- get(".__S3MethodsTable__.", envir = defenv)
            remappedName <- paste(S3[i, 1L], S3[i, 2L], sep = ".")
            if (exists(remappedName, envir = S3Table, inherits = FALSE))
                assign(remappedName, value, S3Table)
        }
    }
    invisible(NULL)
}

#' Test coverage values
#' 
#' Test for mean normalized coverage values
#' 
#' Dryclean, fragcounter, cobalt, or some other
#' binned coverage counter may emit coverage values in  
#' units of reads or mean normalized, gc- and/or mappability-bias
#' corrected coverage units.
#'
#' @export
test_coverage_normalized = function(coverage_values, tolerance = 0.1, fraction_near_1 = 0.9) {
	is_mean_all_equal = identical(all.equal(
		target = 1, 
		current = mean(coverage_values, na.rm = TRUE), 
		tolerance = tolerance
	), TRUE)
	is_median_all_equal = identical(all.equal(
		target = 1, 
		current = median(coverage_values, na.rm = TRUE), 
		tolerance = tolerance
	), TRUE)
	is_cov_na = is.na(coverage_values)
	ix_cov_na = which(is_cov_na)
	if (NROW(ix_cov_na)) {
		coverage_values[ix_cov_na] = -3e9
	}
	histdt = Skilift::test_hist(coverage_values, integer_breaks = -1:5)
	fraction_of_values_inside_integer_breaks = histdt[histdt$is_in_integer_breaks == TRUE, max(cfrac, na.rm = TRUE)]
	is_cov_near_one = fraction_of_values_inside_integer_breaks > fraction_near_1
	is_cov_likely_normalized = is_median_all_equal || is_mean_all_equal
	return(
		list(
			is_cov_likely_normalized = is_cov_likely_normalized,
			is_cov_near_one = is_cov_near_one
		)
	)
}


#' Test histogram
#' 
#' Get histogram of values + custom integer breaks
#' 
#' @export
test_hist = function(values, integer_breaks = -1:5, tolerance = 1e-12) {
	is_integer_breaks_empty = is.null(integer_breaks) || NROW(integer_breaks) == 0
	is_integer_breaks_na = any(is.na(integer_breaks)) 
	# (NROW(integer_breaks == 1) && is.na(integer_breaks)) || any(is.na(integer_breaks))
	if (is_integer_breaks_empty || is_integer_breaks_na) {
		integer_breaks = integer(0)
	}
	hist_breaks = c(
		integer_breaks, 
		min(values, na.rm = TRUE),
		max(values, na.rm = TRUE)
	)
	hist_breaks = unique(hist_breaks)
	histobj = graphics::hist(
		values, 
		breaks = hist_breaks,
    plot = FALSE
	)
	histdt = data.table(from = histobj$breaks[-NROW(histobj$breaks)], to = histobj$breaks[-1])
	histdt$counts = histobj$counts
	is_min_integer_break = abs(integer_breaks - min(integer_breaks, na.rm = TRUE)) <= tolerance
	integer_breaks_no_min = integer_breaks[!is_min_integer_break]
	is_break_min_int = abs(histdt$from - min(integer_breaks_no_min, na.rm = TRUE)) <= tolerance
	is_break_max_int = abs(histdt$to - max(integer_breaks, na.rm = TRUE)) <= tolerance
	ix_is_break_min_int = which(is_break_min_int)
	ix_is_break_max_int = which(is_break_max_int)
	ix_histdt_select = seq(from = ix_is_break_min_int, to = ix_is_break_max_int, by = 1)
	histdt[, csum := cumsum(counts)]
	histdt[, total := sum(counts)]
	histdt[, cfrac := csum / total]
	histdt[, is_in_integer_breaks := FALSE]
	histdt[ix_histdt_select, is_in_integer_breaks := TRUE]
	return(histdt)
}
