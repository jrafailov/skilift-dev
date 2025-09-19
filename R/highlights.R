#' create_heme_highlights
#'
#' Create highlighted heme events for a single sample. 
#' This code sits downstream of create_filtered_events
#' and is used within the lift_filtered_events method.
#' 
#' @export
create_heme_highlights = function(
  events_tbl, ## filtered events R output
  jabba_gg,
  out_file,
  hemedb_path = Skilift:::HEMEDB(),
  duncavage_path = Skilift:::DUNCAVAGEDB(),
  tumor_type = NULL,
  cohorttuple
) {
  ## {
  ##   "karotype": <string> | null,
  ##   "gene_mutations": [
  ##     {
  ##       "gene_name": <string>,
  ##       "variant_p": <string> | null,
  ##       "vaf": <float> | null,
  ##       "altered_copies": <float> | null,
  ##       "total_copies": <float> | null,
  ##       "alteration_type": <enum("small", "cna", "rearrangement")> | null,
  ##       "aggregate_label": <string> | null
  ##     }
  ##     ] | null
  ##   }
  
  
  
  ## TODO: remove everything below except for karyotype_string
  karyotype_string = ""
  if (NROW(jabba_gg) == 1 && is.character(jabba_gg) && file.exists(jabba_gg)) {
    lst = annotate_karyotype(Skilift::process_jabba(jabba_gg))
    karyotype_string = lst$karyotype_string
    karyotype_list = lst$karyotype_list 
  }
    

  emptyDfForJson = structure(list(gene_name = character(0), variant = character(0), 
    vaf = numeric(0), tier = integer(0), altered_copies = numeric(0), 
    total_copies = numeric(0), alteration_type = character(0), 
    aggregate_label = character(0), indication = list()), class = c("data.table", 
  "data.frame"), row.names = c(NA, 0L), sorted = "gene_name")

  hemedb = readRDS(hemedb_path)
  hemedb_guideline = hemedb[, .(GENE, GUIDELINE, DISEASE)][GUIDELINE==TRUE] %>% unique()
  hemedb_guideline = hemedb_guideline[, .(GUIDELINE = GUIDELINE[1], DISEASE = list(DISEASE)), by = GENE]

  ## events_tbl = events_tbl[[1]]

  small_muts = events_tbl[events_tbl$vartype == "SNV",]
  smallForJson = data.table::copy(emptyDfForJson)
  criterias = list(
    is_small_in_guidelines = small_muts$gene %in% hemedb_guideline$GENE, ## accounts for merge step
    is_frequent = small_muts$gene %in% hemedb[hemedb$FREQ >= 5]$GENE,
    is_tier2_or_better = small_muts$Tier <= 2
  )
  is_small_mutation_heme_relevant = base::Reduce("&", criterias)
  if (
    NROW(small_muts) > 0
    && any(
      is_small_mutation_heme_relevant
    )
  ) {
    small_guidelines = small_muts[is_small_mutation_heme_relevant]
    ## Inner join
    small_guidelines = data.table::merge.data.table(small_guidelines, hemedb_guideline, all.x = TRUE, by.x = "gene", by.y = "GENE")

    small_guidelines$variant = ifelse(is.na(small_guidelines$Variant), small_guidelines$type %>% paste(., "Variant"), small_guidelines$Variant)  %>% tools::toTitleCase()
    small_guidelines$alteration_type = "Small"

    small_guidelines$aggregate_label = NA_character_

    small_guidelines$aggregate_label = glue::glue(
      '{small_guidelines[, signif(estimated_altered_copies, 2)]} out of {small_guidelines$segment_cn} copies mutated (VAF: {signif(small_guidelines$VAF, 2)})' 
      ) %>% as.character()

    ## Special cases to highlight biallelic/multi-hit events

    hits_multi = (
      small_guidelines[
        is_multi_hit_per_gene == TRUE,
        .(
          variant = NA_character_,
          VAF = NA_real_,
          Tier = min(Tier, na.rm = TRUE),
          estimated_altered_copies = NA_integer_,
          segment_cn = NA_real_,
          alteration_type = "Multi-Hit",
          aggregate_label = NA_character_,
          DISEASE = DISEASE[1]
        ),
        by = gene
      ]
    )


    changemap = c(
      "gene" = "gene_name",
      "variant" = "variant",
      "VAF" = "vaf",
      "Tier" = "tier",
      "estimated_altered_copies" = "altered_copies",
      "segment_cn" = "total_copies",
      "alteration_type" = "alteration_type",
      "aggregate_label" = "aggregate_label",
      "DISEASE" = "indication"
    )


    smallForJson = base::subset(
      change_names(
        data.table:::rbind.data.table(
          small_guidelines,
          hits_multi,
          fill = TRUE
        ),
        changemap
      ),
      select = changemap
      )
  }

  ## FLT3-ITD special annotation
  is_flt3itd = identical(tolower(readRDS(cohorttuple$itdseek_rds)), "positive")

  flt3ForJson = data.table::copy(emptyDfForJson)
  if (is_flt3itd) {
    dt_flt3 = data.table(
        gene = "FLT3",
        variant = NA_character_,
        VAF = NA_real_,
        Tier = NA_real_,
        estimated_altered_copies = NA_integer_,
        segment_cn = NA_real_,
        alteration_type = "ITD",
        aggregate_label = NA_character_,
        DISEASE = NA
    )


    changemap = c(
      "gene" = "gene_name",
      "variant" = "variant",
      "VAF" = "vaf",
      "Tier" = "tier",
      "estimated_altered_copies" = "altered_copies",
      "segment_cn" = "total_copies",
      "alteration_type" = "alteration_type",
      "aggregate_label" = "aggregate_label",
      "DISEASE" = "indication"
    )


    flt3ForJson = base::subset(
      Skilift:::change_names(
          dt_flt3,
          changemap
      ),
      select = changemap
    )
  }





  dunc = readRDS(duncavage_path)
  ## base::dput(dunc[Gene == "BCR::ABL1"][1])
  dunc = rbind(
    dunc,
    (
      structure(list(Indication = "MPN", Type = structure(2L, .Label = c("SNV", 
        "Fusions/SV", "SV"), class = "factor"), Gene = "BCR::ABL1", DISEASE = "MPN"), row.names = c(NA, 
          -1L), class = c("data.table", "data.frame"))
      [, DISEASE := "ALL"]
      [, Indication := "ALL"]
    ),
    fill = TRUE
  )
  dunc[, DISEASE := gsub("^B-|^T-", "", Indication)]

  hemedb_fusions = dunc[Type == "Fusions/SV", .(Type = Type[1], DISEASE = list(DISEASE)),  by = .(Gene)]

  hemedb_fusions_rev = data.table::copy(hemedb_fusions)
  fg_split = data.table::tstrsplit(hemedb_fusions_rev$Gene, "::")
  hemedb_fusions_rev$Gene = paste(fg_split[[2]], fg_split[[1]], sep = "::")

  hemedb_fusions_fr = rbind(hemedb_fusions, hemedb_fusions_rev)

  svs = events_tbl[grepl("fusion", events_tbl$type, ignore.case = TRUE)]
  svsForJson = data.table::copy(emptyDfForJson)
  fg_exon = svs$fusion_genes ## this still works if svs is empty
  fg_exon = gsub("@.*$", "", fg_exon)
  fg = gsub("\\([0-9]+\\)", "", fg_exon)
  any_svs_in_guidelines = length(intersect(fg, hemedb_fusions_fr$Gene)) > 0 ## accounts for merge step later
  
  if (NROW(svs) > 0 && any_svs_in_guidelines) {
    lst_exons = lapply(strsplit(fg_exon, "::"), function(x) {
      exons = gsub(".*(\\([0-9]+\\)).*", "\\1", x)
      exons = gsub("\\(|\\)", "", exons)
      exons
    })

    svs$lst_exons = lst_exons
    svs$fg = fg


    svs_guidelines = data.table::merge.data.table(svs, hemedb_fusions_fr, by.x = "fg", by.y = "Gene")
    
    exons = do.call(Map, c(f = c, svs_guidelines$lst_exons))

    element_type = dplyr::case_when(
      grepl("UTR", exons[[1]]) & grepl("UTR", exons[[2]]) ~ "UTRs",
      !grepl("UTR", exons[[1]]) & !grepl("UTR", exons[[2]]) ~ "exons",
      TRUE ~ "UTR/exon"
    )

    svs_guidelines$variant = glue::glue('Fusion involving {element_type} {exons[[1]]} <> {exons[[2]]}') %>% as.character()
    svs_guidelines$vaf = NA_real_
    svs_guidelines$alteration_type = "Rearrangement"
    svs_guidelines$aggregate_label = glue::glue('{svs_guidelines$estimated_altered_copies} fusion cop{ifelse(svs_guidelines$estimated_altered_copies == 1, "y", "ies")}') %>% as.character()


    changemap = c(
      "fg" = "gene_name",
      "variant" = "variant",
      "VAF" = "vaf",
      "Tier" = "tier",
      "estimated_altered_copies" = "altered_copies",
      "segment_cn" = "total_copies",
      "alteration_type" = "alteration_type",
      "aggregate_label" = "aggregate_label",
      "DISEASE" = "indication"
    )

    svsForJson = base::subset(
      change_names(svs_guidelines, changemap),
      select = changemap
    )
  }
  

  cna = events_tbl[grepl("SCNA", events_tbl$type, ignore.case = TRUE),]
  cnaForJson = data.table::copy(emptyDfForJson)
  # any_cna_in_guidelines = length(intersect(cna$gene, hemedb_guideline$GENE)) > 0 ## accounts for merge stelater
  criterias = list(
    is_cna_in_guidelines = cna$gene %in% hemedb_guideline$GENE,
    is_tier2_or_better = cna$Tier <= 2
  )
  is_cna_heme_relevant = Reduce("|", criterias)
  if (NROW(cna) > 0 && any(is_cna_heme_relevant)) {
    is_amp = grepl("amp", cna$vartype, ignore.case = TRUE)
    is_homdel = grepl("del", cna$vartype, ignore.case = TRUE)
    is_one_cn_change = cna$estimated_altered_copies == 1
    suffix = ifelse(is_one_cn_change, "y", "ies")
    aggregate_label = data.table::fcase(
      is_amp,
      glue::glue("{cna$estimated_altered_copies} cop{suffix} gained"),
      is_homdel,
      glue::glue("{cna$estimated_altered_copies} cop{suffix} lost")
    ) %>% as.character()
    # cna$aggregate_label = glue::glue('{cna$estimated_altered_copies} out of {cna$estimated_altered_copies} copies altered') %>% as.character()
    cna$aggregate_label = aggregate_label
    cna$alteration_type = "CN"
    cna_guidelines = cna[is_cna_heme_relevant,]
    cna_guidelines = data.table::merge.data.table(cna_guidelines, hemedb_guideline, by.x = "gene", by.y = "GENE", all.x = TRUE)
    cna_guidelines$vartype = tools::toTitleCase(tolower(cna_guidelines$vartype))
    cna_guidelines$DISEASE = lapply(
        cna_guidelines$DISEASE,
        function(x) {
            if (is.null(x)) {
                "MULTIPLE"
            } else {
                x
            }
        }
    )

    changemap = c(
      "gene" = "gene_name",
      "vartype" = "variant",
      "VAF" = "vaf",
      "Tier" = "tier",
      "estimated_altered_copies" = "altered_copies",
      "segment_cn" = "total_copies",
      "alteration_type" = "alteration_type",
      "aggregate_label" = "aggregate_label",
      "DISEASE" = "indication"
    )

    cnaForJson = base::subset(
      change_names(cna_guidelines, changemap),
      select = changemap
    )
  }
  risk_level_eln = NULL

  is_aml = (
    !is.null(tumor_type)
    && (grepl("AML", tumor_type, ignore.case = TRUE)
    || grepl("acute myeloid", tumor_type, ignore.case = TRUE)) 
    && TRUE ## this ensures that this is length one boolean, otherwise errors
  )

  if (is_aml) {
    risk_level_eln = create_heme_highlights_eln_risk_score(
      small_muts = small_muts,
      svs = svs,
      cna = cna,
      karyotype_list = karyotype_list,
      jabba_gg = jabba_gg
    )
  }

  allOutputsForJson = rbind(
    smallForJson,
    cnaForJson,
    svsForJson,
    flt3ForJson
  )

  risk_json_array = list()

  is_eln_risk_populated = !is.null(risk_level_eln) && !any(is.na(risk_level_eln)) && NROW(risk_level_eln) > 0
  if (is_eln_risk_populated) {
    risk_json_array = c(
      risk_json_array,
      ## below syntax is to ensure that the key:value pair gets added as a json record
      list(list("ELN" = jsonlite::unbox(risk_level_eln)))
    )
  }
  

  highlights_output = list(
    karyotype = jsonlite::unbox(karyotype_string),
    ## risk_level_eln = jsonlite::unbox(risk_level_eln),
    risk_score = risk_json_array,
    gene_mutations = allOutputsForJson
  )

  jsonlite::write_json(
    highlights_output, 
    path = out_file,
    na = "null", 
    null = "list",
    pretty = TRUE
  )
  # jsonlite::toJSON(list(gene_mutations = allOutputsForJson), na = "null", null = "list")s

}


#' ELN Risk Score
#' 
#' Logic adapted from Fig 2 of Lachowiez et al. Blood Adv 2022
#' 
#' Heuristic for risk score based on small mutations, 
#' fusion genes (not taken from karyotypes, but the gene level annotations), 
#' and arm level CNA (taken from karyotypes).
create_heme_highlights_eln_risk_score = function(small_muts, svs, cna, karyotype_list, jabba_gg) {
  gg = Skilift:::process_jabba(jabba_gg)

  levels_eln = c("Favorable" = "Favorable", "Intermediate" = "Intermediate", "Adverse" = "Adverse")
  sv_mat = stringr::str_split_fixed(svs$gene, "::", 2)
  sv_5p = sv_mat[,1]
  sv_3p = sv_mat[,2]

  is_eln_risk_favorable = any(svs$gene %in% c("CBFB::MYH11", "RUNX1::RUNX1T1"))
  if (is_eln_risk_favorable) return(levels_eln["Favorable"])

  is_eln_risk_adverse_fusions = any(svs$gene %in% c("DEK::NUP214", "BCR::ABL1", "KAT6A::CREBBP", "GATA2::MECOM"))
  is_eln_risk_adverse_rearranged_5p = any(sv_5p %in% c("KMT2A"))
  is_eln_risk_adverse_rearranged_3p = any(sv_3p %in% c("MECOM"))
  is_eln_risk_adverse = (
    is_eln_risk_adverse_fusions 
    || is_eln_risk_adverse_rearranged_5p
    || is_eln_risk_adverse_rearranged_3p
  )
  if (is_eln_risk_adverse) return(levels_eln["Adverse"])

  is_eln_risk_adverse = (
    any(grepl("-5", karyotype_list$annotated_chr_cna))
    || any(grepl("-7", gsub("\\(|\\)", "", karyotype_list$annotated_chr_cna)))
    || any(grepl("-17", gsub("\\(|\\)", "", karyotype_list$annotated_chr_cna)))
    || any(grepl("del5p", gsub("\\(|\\)", "", karyotype_list$annotated_arm_cna)))
    || any(grepl("del17p", gsub("\\(|\\)", "", karyotype_list$annotated_arm_cna)))
  )
  if (is_eln_risk_adverse) return(levels_eln["Adverse"])

  ploidy_genome = gg$meta$ploidy
  gr_nodes = gg$nodes$gr
  GenomeInfoDb::seqlevelsStyle(gr_nodes) = "NCBI"
  gr_nodes = GenomeInfoDb::sortSeqlevels(gr_nodes)
  gr_nodes = sort(gr_nodes, ignore.strand = TRUE)
  ## FIXME: Need to account for Y once finished.
  gr_nodes = (
      gr_nodes[
          as.character(seqnames(gr_nodes))
          %in% c(
            1:22,
            "X"
           ## ,
            ## "Y"
          )
      ]
  )
  dt_nodes = gr2dt(gr_nodes)
  ploidy_chrom_per = dt_nodes[!is.na(cn), .(ploidy_chrom = sum(cn * width) / sum(width)), by = seqnames]

  ploidy_chrom = sum(round(ploidy_chrom_per$ploidy_chrom))

  ## FIXME: Y
  cytomap = Skilift:::create_cytomap()
  MIN_DIST=5e6
  PAD=round(MIN_DIST / 2)
  gr_aberrant_seg = GenomicRanges::reduce(gr_nodes[
    gr_nodes$cn != 2 & !as.logical(seqnames(gr_nodes) == "Y")
  ] + PAD) - PAD
  seqdf = as.data.frame(GenomeInfoDb::seqlengths(gr_nodes))
  gr_aberrant_seg$chromwidth = seqdf[as.character(seqnames(gr_aberrant_seg)),]
  is_chrom = as.integer(width(gr_aberrant_seg)) / gr_aberrant_seg$chromwidth > 0.9
  is_large = width(gr_aberrant_seg) > MIN_DIST
  gr_aberrant_seg = gr_aberrant_seg[!is_chrom & is_large]
  dtov = gr2dt(gr_aberrant_seg %*% dt2gr(cytomap))
  dtov = dtov[, .(fracov = sum(width) / category_width[1], band = list(band)), by = .(seqnames, query.id, category)]
  dtov = dtov[!category == "chromosome" & fracov > 0.99]
  dtov = dtov[!category == "arm" & fracov > 0.99]
  dtovbigbands = dtov

  chromosomal_abnormalities = c(
    karyotype_list$annotated_chr_cna,
    karyotype_list$annotated_arm_cna
  )

  
  

  is_hyperdiploid = sum(ploidy_chrom) > 46 && ploidy_genome > 2.7
  is_hypodiploid = sum(ploidy_chrom) <= 24 && ploidy_genome < 2
  is_complex = (
      ! is_hyperdiploid
      & NROW(chromosomal_abnormalities) + NROW(dtovbigbands) >= 3
  )
  is_eln_risk_adverse = is_complex || is_hypodiploid

  if (is_eln_risk_adverse) return(levels_eln["Adverse"])

  is_eln_risk_intermediate = any(svs$gene %in% c("MLLT3::KMT2A", "KMT2A::MLLT3"))
  if (is_eln_risk_intermediate) return(levels_eln["Intermediate"])
  
  is_gene_tp53 = small_muts$gene == "TP53"
  is_gene_tp53_vaf_gr10 = is_gene_tp53 & (small_muts$VAF >= 0.1 & !is.na(small_muts$VAF))
  is_eln_risk_adverse = any(is_gene_tp53_vaf_gr10)

  if (is_eln_risk_adverse) return(levels_eln["Adverse"])

  is_flt3 = svs$gene == "FLT3"
  is_itd = grepl("ITD", svs$Variant)
  is_flt3_itd = is_flt3 & is_itd
  is_any_flt3_itd = any(is_flt3_itd)

  epigene_list = c("ASXL1", "BCOR", "EZH2", "RUNX1", "SF3B1", "STAG2", "U2AF1", "ZRSR2")
  is_gene_epigene_small_mut = small_muts$gene %in% epigene_list
  is_gene_epigene_cna = cna[vartype == "HOMDEL"]$gene %in% epigene_list
  is_any_gene_epigene_altered = any(is_gene_epigene_small_mut) || any(is_gene_epigene_cna)
  
  
  is_eln_risk_adverse = is_any_flt3_itd && is_any_gene_epigene_altered
  
  if (is_eln_risk_adverse) return(levels_eln["Adverse"])

  is_eln_risk_intermediate = is_any_flt3_itd && !is_any_gene_epigene_altered

  if (is_eln_risk_intermediate) return(levels_eln["Intermediate"])  

  is_any_npm1_mutated = any(small_muts$gene %in% "NPM1")
  is_eln_risk_favorable = is_any_flt3_itd &&  is_any_npm1_mutated
  
  if (is_eln_risk_favorable) return(levels_eln["Favorable"])
  
  is_any_cebpa_multihit = any(small_muts$gene == "CEBPA" & small_muts$is_multi_hit_per_gene)

  variantp = regmatches(small_muts$Variant, gregexpr("(p\\.[a-zA-Z0-9_]+)", small_muts$Variant))
  variantp[base::lengths(variantp) == 0] = ""
  variantp = dunlist(variantp)

  variantp_coord = regmatches(variantp$V1, gregexpr("[0-9]+", variantp$V1))
  variantp_coord[base::lengths(variantp_coord) == 0] = ""
  variantp$coord = variantp_coord

  dt_is_in_bzip = variantp[, any(as.integer(coord[[1]]) %in% c(278:345)), keyby = listid]

  is_sanity_check_invalid = !(NROW(dt_is_in_bzip) == NROW(small_muts))
  if (is_sanity_check_invalid) {
      stop('small muts parsing is incorrect for eln risk scoring')
  }
  is_in_bzip = dt_is_in_bzip$V1
  is_any_cebpa_bzip = any(small_muts$gene == "CEBPA" & is_in_bzip)

  is_eln_risk_favorable = is_any_cebpa_multihit && is_any_cebpa_bzip
  

  if (is_eln_risk_favorable) return(levels_eln["Favorable"])

  is_eln_risk_adverse = is_any_gene_epigene_altered
  
  if (is_eln_risk_adverse) return(levels_eln["Adverse"])
  
  return(levels_eln["Intermediate"]) 

}



#' Change names of data.frame/matrix
#' 
#' Mean to be an easy wrapper that uses either two
#' equal length character vectors, or a single named
#' character vector to rename the columns of a data frame
#' or matrix.
#' 
#' 
change_names = function(obj, old, new) {
  onames = names(obj)
  if (missing(new) && length(names(old)) > 0 && length(old) > 0) {
    new = old
    old = names(old)
  }
  newnames = onames
  if (length(old) != length(new)) stop("lengths of old and new vectors must match")
  for (i in 1:NROW(old)) {
    if (!nzchar(old[i])) next
    if (old[i] == new[i]) next
    newnames[newnames %in% old[i]] = new[i]
  }
  names(obj) = newnames
  return(obj)
}




#' create_summary
#'
#' Create high level summary of
#' mutations that will get passed
#' to metadata.json.
#' 
#' @export
create_summary = function(
  events_tbl, ## filtered events R output
  altered_copies_threshold = 0.9,
  cohort_type
) {

  small_muts = events_tbl[events_tbl$vartype == "SNV",]
  hemedb = readRDS(Skilift:::HEMEDB())
  hemedb_guideline = hemedb[, .(GENE, GUIDELINE, DISEASE)][GUIDELINE==TRUE] %>% unique()
  hemedb_guideline = hemedb_guideline[, .(GUIDELINE = GUIDELINE[1], DISEASE = list(DISEASE)), by = GENE]

  criterias = list(
    is_tier_or_better = small_muts$Tier <= 1
   ,
    is_clonal = is.na(small_muts$estimated_altered_copies) | small_muts$estimated_altered_copies >= altered_copies_threshold
	, # allow NA's through
	is_small_in_guidelines = small_muts$gene %in% hemedb_guideline$GENE, ## heme relevant only
	is_frequent = small_muts$gene %in% hemedb[hemedb$FREQ >= 5]$GENE ## heme relevant only

  )
  if (!cohort_type == "heme") {
	criterias = criterias[1:2]
  }
  is_small_mutation_relevant = base::Reduce("&", criterias)
  small_muts_parsed = ""
  if (
    NROW(small_muts) > 0
    && any(
      is_small_mutation_relevant
    )
  ) {
    small_muts_out = small_muts[is_small_mutation_relevant]
    hits_multi = (
      small_muts_out[
        is_multi_hit_per_gene == TRUE,
        .(
          type = "Multi-Hit"
        ),
        by = gene
      ]
    )
    small_muts_out = data.table:::rbind.data.table(
      small_muts_out,
      hits_multi,
      fill = TRUE
    )

    small_muts_tally = (
      base::subset(small_muts_out)[, .(gene, type)]
      [, table(type, gene)]
      %>% reshape2::melt()
      %>% as.data.table()
    )
    small_muts_tally = small_muts_tally[value > 0]
    # small_muts_parsed = paste(
    #   small_muts_tally[, paste(type, ": ", paste(gene, collapse = ","), sep = ""), by = type]$V1,
    #   collapse = "; "
    # )

	small_muts_parsed = small_muts_tally[, paste(type, ": ", gene, sep = ""), by = type]$V1
    
  }
  
  dunc = readRDS(Skilift:::DUNCAVAGEDB())
  dunc = rbind(
    dunc,
    (
      structure(list(Indication = "MPN", Type = structure(2L, .Label = c("SNV", 
        "Fusions/SV", "SV"), class = "factor"), Gene = "BCR::ABL1", DISEASE = "MPN"), row.names = c(NA, 
          -1L), class = c("data.table", "data.frame"))
      [, DISEASE := "ALL"]
      [, Indication := "ALL"]
    ),
    fill = TRUE
  )
  dunc[, DISEASE := gsub("^B-|^T-", "", Indication)]
  hemedb_fusions = dunc[Type == "Fusions/SV", .(Type = Type[1], DISEASE = list(DISEASE)),  by = .(Gene)]
  hemedb_fusions_rev = data.table::copy(hemedb_fusions)
  fg_split = data.table::tstrsplit(hemedb_fusions_rev$Gene, "::")
  hemedb_fusions_rev$Gene = paste(fg_split[[2]], fg_split[[1]], sep = "::")

  hemedb_fusions_fr = rbind(hemedb_fusions, hemedb_fusions_rev)

  svs = events_tbl[grepl("fusion", events_tbl$type, ignore.case = TRUE)]
  #   svsForJson = data.table::copy(emptyDfForJson)
  fg_exon = svs$fusion_genes ## this still works if svs is empty
  fg_exon = gsub("@.*$", "", fg_exon)
  fg = gsub("\\([0-9]+\\)", "", fg_exon)

  sv_recurrent_map = hemedb_fusions_fr$Gene
  cosmic = fread(Skilift:::COSMIC_FUSIONS())
  sv_recurrent_map = c(
    sv_recurrent_map, 
    unique(cosmic[, paste(FIVE_PRIME_GENE_SYMBOL, THREE_PRIME_GENE_SYMBOL, sep = "::")])
  )

  any_svs_in_guidelines = length(intersect(fg, sv_recurrent_map)) > 0 ## accounts for merge step later
  criterias = list(
    is_tier2_or_better = svs$Tier <= 2,
	  is_svs_in_guidelines = fg %in% sv_recurrent_map
  )
  if (!cohort_type == "heme") {
	criterias = criterias[1]
  }
  is_svs_relevant = Reduce("&", criterias)

  svs_parsed = ""
  if (
    NROW(svs) > 0
    && any(
      is_svs_relevant
    )
  ) {
    # svs_parsed = (
    #   base::subset(svs, is_svs_relevant)[, .(gene, type)]
    #     [, paste(type, ": ", paste(gene, collapse = ","), sep = ""), by = type]
    #     $V1
    # )

	svs_parsed = (
      base::subset(svs, is_svs_relevant)[, .(gene, type)]
        [, paste(type, ": ", gene, sep = ""), by = type]
        $V1
    )

  }
  

  cna = events_tbl[grepl("SCNA", events_tbl$type, ignore.case = TRUE),]
  cna_parsed = ""
  criterias = list(
    is_tier_or_better = cna$Tier <= 2,
	is_cna_in_guidelines = cna$gene %in% hemedb_guideline$GENE
  )
  if (!cohort_type == "heme") {
	criterias = criterias[1]
  }
  is_cna_relevant = Reduce("&", criterias)
  if (NROW(cna) > 0 && any(is_cna_relevant)) {
    cna_out = cna[is_cna_relevant,]
    cna_out$vartype = tools::toTitleCase(tolower(cna_out$vartype))


    cna_tally = (
      base::subset(cna_out)[, .(gene, vartype)]
      [, table(vartype, gene)]
      %>% reshape2::melt()
      %>% as.data.table()
    )
    cna_tally = cna_tally[value > 0]
    # cna_parsed = paste(
    #   cna_tally[, paste(vartype, ": ", paste(gene, collapse = ","), sep = ""), by = vartype]$V1,
    #   collapse = "; "
    # )
	cna_parsed = (
      cna_tally[, paste(vartype, ": ", gene, sep = ""), by = vartype]$V1
    )
  }

  summary_string = paste(
    paste(small_muts_parsed, collapse = "\n"),
    paste(cna_parsed, collapse = "\n"),
    paste(svs_parsed, collapse = "\n"),
    sep = "\n"
  )
  ## summary_string = paste(small_muts_parsed, cna_parsed, svs_parsed, sep = "\n")
  summary_string = trimws(summary_string)
  summary_string = gsub("\n{2,}", "\n", summary_string)

  return(summary_string)
  

}
