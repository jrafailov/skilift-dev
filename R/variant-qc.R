#' @name create_variant_qc
#' @title create_variant_qc
#' @description
#'
#' Function to create quality control metrics for somatic variants
#'
#' @param oncokb_path Path to oncokb maf file to be used 
#' @param multiplicity_path Path to multiplitcity output
#' @param sbs_path Path to posterior probs for sbs signatures
#' @param indel_path Path to posterior probs for ID signatures
#' @return data.table containing variant QC metrics along with highest signature posterior
#' @export
#' @author Kevin Hadi, Aditya Deshpande
create_variant_qc = function(oncokb_path, multiplicity_path, sbs_path, indel_path, cohort_type){
  oncokb = fread(oncokb_path)
  oncokb = Skilift:::parse_oncokb_tier(oncokb)
  oncokb = Skilift::merge_oncokb_multiplicity(oncokb, multiplicity_path)
  message("read in the oncokb maf")
  ## Columns to keep
  columns_to_keep = c("Hugo_Symbol",
                      "Chromosome",
                      "Start_Position",
                      "End_Position",
                      "Strand",
                      "Variant_Classification",
                      "Variant_Type",
                      "Reference_Allele",
                      "Tumor_Seq_Allele1",
                      "Tumor_Seq_Allele2",
                      "HGVSc",
                      "HGVSp_Short",
                      "Consequence",
                      "t_depth",
                      "t_ref_count",
                      "t_alt_count",
                      "FILTER",
                      "ONCOGENIC",
                      "tier",
                      "tier_factor",
                      "Role",
                      "cn",
                      "altered_copies")

  tumor_normal_columns = c("n_depth",
                           "n_ref_count",
                           "n_alt_count")

  if (cohort_type == "paired"){
    columns_to_keep = c(columns_to_keep, tumor_normal_columns)
  }
  
  oncokb = oncokb[, ..columns_to_keep]
  oncokb[,  key_col := paste0(Chromosome, "_", as.character(Start_Position))]       

  ## Get signature posteriors for sbs
  sbs = fread(sbs_path)
  num_cols = names(sbs)[sapply(sbs, is.numeric)]
  keep_cols = num_cols[colSums(sbs[, ..num_cols]) != 0]
  sbs = sbs[, c(keep_cols, setdiff(names(sbs), num_cols)), with = FALSE]
  sbs_m = as.data.table(melt(sbs, id.vars = c("Sample Names", "Chr", "Pos", "MutationType")))
  sbs_m[, key_col := paste0(Chr, "_", as.character(Pos))]
  setkey(sbs_m, "key_col")
  sbs_m = sbs_m[, .SD[which.max(value)], by = key_col]
  sbs_m[, sbs_signature := paste0(variable, ": ", signif(value, 0.2))]
  message("read in sbs signatures")

  ## Get signature posteriors for indels
  indels = fread(indel_path)
  num_cols = names(indels)[sapply(indels, is.numeric)]
  keep_cols = num_cols[colSums(indels[, ..num_cols]) != 0]
  indels = indels[, c(keep_cols, setdiff(names(indels), num_cols)), with = FALSE]
  indels_m = as.data.table(melt(indels, id.vars = c("Sample Names", "Chr", "Pos", "MutationType")))
  indels_m[, Pos := ifelse(grepl("Del", MutationType), Pos+1, Pos)]
  indels_m[, key_col := paste0(Chr, "_", as.character(Pos))]
  setkey(indels_m, "key_col")
  indels_m = indels_m[, .SD[which.max(value)], by = key_col]
  indels_m[, indel_signature := paste0(variable, ": ", signif(value, 0.2))]
  message("read in indel signatures")
  all_sigs = rbind(sbs_m, indels_m, fill =  T)
  all_sigs[, signature := ifelse(!is.na(sbs_signature), sbs_signature, indel_signature)]
  all_sigs[, signature_max := variable]

  ## Remove SBS sigs from DNP and TNP
  oncokb_merged = as.data.table(merge(oncokb, all_sigs[, .(key_col, signature, signature_max)], by = "key_col", all.x = T))
  oncokb_merged[, signature := ifelse(Variant_Type %in% c("DNP", "TNP"), NA, signature)]
  message("calculating stats")

  ## Calculate and aggregate 
  oncokb_merged[, tumor_vaf := t_alt_count/t_depth]
  if (cohort_type == "paired"){
    oncokb_merged[, normal_vaf := n_alt_count/(n_ref_count + n_alt_count)]
  }
  is_oncokb_varp_populated = nzchar(oncokb_merged$HGVSp_Short) & !is.na(oncokb_merged$HGVSp_Short)
  is_oncokb_varc_populated = nzchar(oncokb_merged$HGVSc) & !is.na(oncokb_merged$HGVSc)
  use_oncokb_varp = is_oncokb_varp_populated & !is_oncokb_varc_populated
  use_oncokb_varc = !is_oncokb_varp_populated & is_oncokb_varc_populated
  use_oncokb = use_oncokb_varp | use_oncokb_varc
  vardisplay = data.table::fcase(
    use_oncokb_varp, oncokb_merged$HGVSp_Short,
    use_oncokb_varc, oncokb_merged$HGVSc
    )
  oncokb_merged$variant = vardisplay
  
  lst_strings_t = with(as.list((oncokb_merged)), {
    list(
      tumor_string = as.character(stringr::str_c(t_alt_count, " / ", t_depth, " / ", signif(tumor_vaf, 2)))
    )
  })
  oncokb_merged$tumor_alt_dp_vaf = lst_strings_t$tumor_string

  if (cohort_type == "paired"){
    lst_strings_n = with(as.list((oncokb_merged)), {
      list(
        normal_string = as.character(stringr::str_c(n_alt_count, " / ", n_depth, " / ", signif(normal_vaf, 2)))
      )
    })
    oncokb_merged$normal_alt_dp_vaf = lst_strings_n$normal_string
  }

  cols_to_pass = c(
    "Chromosome" = "Chromosome",
    "Start_Position" = "position",
    "Reference_Allele" = "reference",
    "Tumor_Seq_Allele2" = "alternate",
    "Hugo_Symbol" = "gene",
    "variant" = "variant",
    "ONCOGENIC" = "annotation",
    "tumor_alt_dp_vaf" = "tumor_alt_dp_vaf",
    "t_depth" = "tumor_depth",
    "t_alt_count" = "tumor_alt_counts",
    "tumor_vaf" = "tumor_vaf",
    "altered_copies" = "altered_copies",
    "cn" = "total_cn",
    "tier" = "oncokb_tier"
  )

  if (cohort_type == "paired"){
    cols_to_pass = c(cols_to_pass,
                     c("n_depth" = "normal_depth",
                       "n_alt_count" = "normal_alt_counts",
                       "normal_vaf" = "normal_vaf",
                       "normal_alt_dp_vaf" = "normal_alt_dp_vaf"
                       )
                     )
  }
  
  meta_gr_specials = c("seqnames", "start")
  dt_cols_to_pass = names(cols_to_pass)[!names(cols_to_pass) %in% meta_gr_specials]
  sq = Skilift:::change_names(oncokb_merged, cols_to_pass)

  consider_numeric <- c(
    "tumor_depth",
    "tumor_alt_counts",
    "tumor_vaf",
    "altered_copies",
    "total_cn",
    "oncokb_tier"
  )

  if (cohort_type == "paired"){
    consider_numeric = c(consider_numeric,
                     c("normal_depth",
                       "normal_alt_counts",
                       "normal_vaf"
                       )
                     )
  }
  
  for (col in consider_numeric) {
    data.table::set(sq, j = col, value = as.numeric(sq[[col]]))
  }
  sq[, id := .I]

  ## Convert all logical fields to character
  jx = seq_len(NCOL(sq))
  for (j in jx) {
    val = sq[[j]]
    is_null_or_other = is.null(val) || !is.logical(val)
    if (is_null_or_other) next
    data.table::set(sq, j = j, value = as.character(sq[[j]]))
  }

  return(sq)
}





#' @name lift_variant_qc
#' @title lift_variant_qc
#' @description
#' Creates variant QC metrics for all samples in a cohort and writes them to JSON files
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @return None
#' @export
lift_variant_qc <- function(cohort, output_data_dir, cores = 1) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }
    
    # Validate required columns exist
    if (any(!c("oncokb_snv", "multiplicity", "indel_post_prob_signatures", "sbs_post_prob_signatures") %in% names(cohort$inputs))) {
        stop("Missing required columns  in cohort")
    }
    
    # Process each sample in parallel
    mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        row <- cohort$inputs[i,]
        pair_dir <- file.path(output_data_dir, row$pair)
        
        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }
        
        out_file <- file.path(pair_dir, "sage.qc.json")
        
        futile.logger::flog.threshold("ERROR")
        tryCatchLog({
            # Generate QC metrics
          qc_data <- create_variant_qc(
            oncokb_path = row$oncokb_snv,
            multiplicity_path = row$multiplicity,
            sbs_path = row$sbs_post_prob_signatures,
            indel_path = row$indel_post_prob_signatures,
            cohort_type = cohort$type
            )
            # Write to JSON
            jsonlite::write_json(qc_data, out_file, pretty = TRUE)
        }, error = function(e) {
            print(sprintf("Error processing %s: %s", row$pair, e$message))
        })
    }, mc.cores = cores, mc.preschedule = TRUE)
    
    invisible(NULL)
}
