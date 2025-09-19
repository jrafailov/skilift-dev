#' SV in GRangesList to bedpe
#' 
#' Convert GRangesList representation of SVs to bedpe-like format
#' 
#' @export
grl2bedpe = function (
	grl, 
	add_breakend_mcol = FALSE, 
	flip = FALSE, 
	as.data.table = TRUE, 
    zerobased = TRUE
) {
    grpiv = gUtils::grl.pivot(grl)
    if (zerobased) {
      GenomicRanges::start(grpiv[[1]]) = GenomicRanges::start(grpiv[[1]]) - 1
      GenomicRanges::start(grpiv[[2]]) = GenomicRanges::start(grpiv[[2]]) - 1
    }
    mcgrl = as.data.frame(mcols(grl))
    df1 = as.data.frame(grpiv[[1]])[, c(1:3, 5), drop = F]
    df2 = as.data.frame(grpiv[[2]])[, c(1:3, 5), drop = F]
    colnames(df1) = c("chrom1", "start1", "end1", "strand1")
    colnames(df2) = c("chrom2", "start2", "end2", "strand2")
    mc1 = data.frame()[seq_len(NROW(df1)), , drop = F]
    mc2 = data.frame()[seq_len(NROW(df2)), , drop = F]
    if (identical(add_breakend_mcol, TRUE)) {
        mc1 = as.data.frame(grpiv[[1]])[, -c(1:5), drop = F]
        mc2 = as.data.frame(grpiv[[2]])[, -c(1:5), drop = F]
        colnames(mc1) = paste0("first.", colnames(mc1))
        colnames(mc2) = paste0("second.", colnames(mc2))
    }
    out = cbind(df1, df2, name = as.character(seq_len(NROW(df1))), 
        score = rep_len(0, NROW(df1)), mcgrl, mc1, mc2)
    canon_col = c(1, 2, 3, 5, 6, 7, 9, 10, 4, 8)
    nix = seq_len(ncol(out))
    out = out[, c(canon_col, nix[!nix %in% canon_col]), drop = F]
    if (flip) {
        out$strand1 = c(`+` = "-", `-` = "+")[out$strand1]
        out$strand2 = c(`+` = "-", `-` = "+")[out$strand2]
    }
    if (as.data.table) 
        return(data.table::as.data.table(out))
    else return(out)
}

#' Annotate Karyotype
#' 
#' Convert gGraph to karyotype
#' 
#' @export

annotate_karyotype = function(
	gg, 
	segment_size_threshold = 5e6, 
	chr_fraction_threshold = 0.9, 
	arm_fraction_threshold = 0.9, 
	band_segment_ratio_threshold = 0.1,
    event_count_thresh = 5
) {
  cyto = Skilift:::process_cytoband()
  cyto$arm = gsub("^(p|q).*", "\\1", cyto$band)
  cyto$band_width = as.integer(width(cyto))

  cyto_by_arm = gUtils::gr_construct_by(cyto, "arm")
  rcyto_by_arm = GenomicRanges::reduce(cyto_by_arm)
  rcyto_arm = gUtils::gr_deconstruct_by(rcyto_by_arm, by = "arm", meta = TRUE)
  rcyto_arm_dt = gUtils::gr2dt(rcyto_arm)

  annotated_tra_inv_raw = annotated_tra_inv = ""

  ge = gg$edges
  geAlt = ge[ge$dt$type == "ALT"]
  grl_alt = GRangesList()
  is_any_alt_edge_present = NROW(geAlt) > 0
  # Maybe unnecessary but still..
  if (is_any_alt_edge_present) {
    grl_alt = geAlt$grl
    ## Pre-filtering grl_alt by cytobands
    ## to get rid of odd junctions
    ## e.g. Autosomal :: MT junctions
    ## e.g. Junctions mapping to decoy chromosomes
    is_in_cyto = gUtils::grl.in(
      grl = grl_alt,
      windows = cyto,
      only = TRUE ## means returns TRUE for each grl element if ALL windows match granges in grl.
    )
    grl_alt = grl_alt[is_in_cyto]
  }
  if (length(grl_alt) > 0) {
    sortedGrl = BiocGenerics::sort(GenomeInfoDb::sortSeqlevels(grl_alt), ignore.strand = TRUE)

    widths = width(range(sortedGrl, ignore.strand = TRUE))
    widths[base::lengths(widths) == 2] = list(NA_integer_)
    mcols(sortedGrl)$widths = unlist(widths)
    breakpoints = gUtils::grl.unlist(sortedGrl)
    breakpoints = gUtils::gr.end(breakpoints, ignore.strand = FALSE)
    breakpoints = breakpoints %*% cyto

    grl = S4Vectors::split(breakpoints, breakpoints$grl.ix)
  	bedpe = grl2bedpe(grl, add_breakend_mcol = TRUE, zerobased = FALSE)
    
    isTranslocation = !(bedpe$first.chrom_name == bedpe$second.chrom_name)
    k_type = dplyr::case_when(
      bedpe$first.class == "TRA-like" ~ "t",
      bedpe$first.class == "INV-like" ~ "inv",
      TRUE ~ NA_character_
    )
    k_chrom = ifelse(
      bedpe$first.chrom_name == bedpe$second.chrom_name,
      bedpe$first.chrom_name,
      paste(bedpe$first.chrom_name, bedpe$second.chrom_name, sep = ";")
    )
    isBig = (
		(
			isTranslocation | 
			(bedpe$first.widths > segment_size_threshold)
		)
		& !is.na(k_type)
	)
    k_bands = paste(bedpe$first.band, bedpe$second.band, sep = ";")

    if (any(isBig))
      annotated_tra_inv = glue::glue('{k_type}({k_chrom})({k_bands})')[isBig]

    ## If the length is too big, don't show. Set to 5 for now
    annotated_tra_inv_raw = annotated_tra_inv
    if(length(annotated_tra_inv) > event_count_thresh){
      annotated_tra_inv = "multiple tra/inv events"  
    }

  }

  segments = gg$nodes$gr
  segments = BiocGenerics::sort(
	GenomeInfoDb::sortSeqlevels(segments), 
	ignore.strand = TRUE
  )
  segments = GenomeInfoDb::keepStandardChromosomes(
	segments, 
	pruning.mode = "coarse"
  )
  genome_sl = GenomeInfoDb::seqlengths(segments)
  genome_sl = data.table::data.table(
    seqnames = names(genome_sl),
    chr_length = genome_sl
  )

  annotated_band_cna_raw = annotated_cna = ""

  events = segments[
    width(segments) > segment_size_threshold ## Removing small segs, will collapse eventually
    & seqnames(segments) != "Y"
  ]

  if (length(events) > 0) {
    ## Preserving cn for the chromosome level annotations
    events$label = ifelse(events$cn > 2, paste0("dup_", events$cn),
                          ifelse(events$cn == 2, paste0("dip_", events$cn),
                          paste0("del_", events$cn)))
    events_by_label = gUtils::gr_construct_by(events, "label")
    revents_by_label = GenomicRanges::reduce(events_by_label, ignore.strand = TRUE)
    revents = gUtils::gr_deconstruct_by(revents_by_label, by = "label", meta = TRUE)
    revents = revents %*% cyto
    revents_dt = gUtils::gr2dt(revents)
    revents_dt[, segment_arm_id := .GRP, by = .(query.id, arm)]
    revents_dt = data.table::merge.data.table(revents_dt, genome_sl, by = "seqnames")

    ## Some maneuvering to get the chromosomes where entire chr
    ## is gained and more events are overlayed
    revents_dt[, c("label", "cn") := tstrsplit(label, "_")]
    revents_dt[, cn := as.numeric(cn)]
    ## This is to assume a base level CN state for each CHR
    min_cn_dt = unique(revents_dt[, .(chr_lv_cn = if (all(cn > 2)) min(cn[cn > 2]) else (if (any(cn == 2)) 2 else max(cn))), by = seqnames])
    revents_dt = merge(revents_dt, min_cn_dt, by = "seqnames", all.x = T)
    revents_dt[, chr_level_group := paste0(label, "_", chr_lv_cn)]

    ## CHR LEVEL
    chr_stats = revents_dt[, sum(width) / chr_length[1],
                           by = .(seqnames, chr_level_group)]
    chr_events = data.table::copy(chr_stats[V1 > chr_fraction_threshold])
    chr_events = chr_events[!(grepl("dip", chr_level_group))]

    if (nrow(chr_events) > 0){
        chr_events[, c("label", "cn") := tstrsplit(chr_level_group, "_")]
        chr_events[, cn := as.numeric(cn)]

    ## This ensures that if there are 2 or more extra copies that gets
    ## propogated as per the convention
        chr_events[, suffix := ifelse(grepl("dup", label),
                                      ifelse(cn > 3, paste0("x", as.numeric(cn)-2), ""), "")]
        chr_events[, label := ifelse(grepl("dup", label), "+", "-")]
        annotated_chr_cna = with(chr_events, glue::glue('{label}{seqnames}{suffix}'))
    } else {
      annotated_chr_cna = ""
    }

    ## ARM LEVEL
    revents_arm_dt = data.table::merge.data.table(revents_dt, rcyto_arm_dt[, .(seqnames, arm, arm_width = width)], by = c("seqnames", "arm"))
    ## No need to exclude events that were in whole chr, since we can now remove segs
    ## that matched whole chr lvele cn events above
    if(nrow(chr_events) > 0){
        revents_arm_dt = revents_arm_dt[cn != chr_lv_cn]
    }
    revents_arm_dt = revents_arm_dt[!(grepl("dip", chr_level_group))] 
    arm_stats = revents_arm_dt[, sum(width) / arm_width[1],
                                                by = .(seqnames, label, arm)]
    arm_events = arm_stats[V1 > arm_fraction_threshold]
    annotated_arm_cna = with(arm_events, glue::glue('{label}({seqnames})({arm})'))

    ## BAND LEVEL
    ## Same for band level
    band_events_prefiltered = data.table::copy(revents_arm_dt)
    band_events_prefiltered = band_events_prefiltered[stain != "acen"]
    band_events_prefiltered[, band_segment_ratio := sum(width) / band_width[1], by = .(segment_arm_id)]
    band_events_prefiltered[, total_segment_width := sum(width), by = segment_arm_id]
    ## Original segment overlapping with band must be > segment size threshold, and also must cover > 10% of the coincident cytoband width.
    band_events = band_events_prefiltered[total_segment_width > segment_size_threshold & band_segment_ratio >= band_segment_ratio_threshold]

    band_events_collapsed = band_events[, .(
      seqnames = seqnames[1],
      min_band = band[1],
      max_band = band[.N],
      arm = arm[1],
      label = label[1]
    ), by = .(segment_arm_id)]

   ## Additional logic to collapse adjacent bands that have the same copy states 
   band_events_collapsed[, diff_next := abs(segment_arm_id - data.table::shift(segment_arm_id, type = "lag"))]
   band_events_collapsed[, diff_prev := abs(segment_arm_id - data.table::shift(segment_arm_id, type = "lead"))]
   band_events_collapsed[, diff := min(diff_next, diff_prev, na.rm = T), by = segment_arm_id]
   band_events_collapsed[, group_adj := paste0(seqnames, "_", arm, "_", label, "_", diff)]

   band_events_collapsed = band_events_collapsed[, .(
      seqnames = seqnames[1],
      min_band = min_band[1],
      max_band = max_band[.N],
      arm = arm[1],
      label = label[1]
    ), by = .(group_adj)]

    band_events_collapsed[, first_band := ifelse(arm == "p", max_band, min_band)]
    band_events_collapsed[, last_band := ifelse(arm == "p", min_band, max_band)]
    
    ## Filtering out arm level events
    if(nrow(arm_events) > 0){
      arm_events[, long_label := paste0(seqnames, "_", label, "_", arm)]
      band_events_collapsed[, long_label := paste0(seqnames, "_", label, "_", arm)]
      band_events_collapsed = band_events_collapsed[!long_label %in% arm_events$long_label]
    }

    annotated_band_cna = with(band_events_collapsed, {
      glue::glue('{label}({seqnames})({first_band},{last_band})')
    })

    ## If the length is too big, don't show. Set to 5 for now
    annotated_band_cna_raw = annotated_band_cna
    if(length(annotated_band_cna) > event_count_thresh){
      annotated_band_cna = "multiple band events"  
    }
  }

  karyotype_list = list(
    annotated_chr_cna = annotated_chr_cna,
    annotated_tra_inv = annotated_tra_inv_raw,
    annotated_arm_cna = annotated_arm_cna,
    annotated_band_cna = annotated_band_cna_raw
  )
  
  karyotype_string = paste(c(
    annotated_chr_cna,
    annotated_tra_inv,
    annotated_arm_cna,
    annotated_band_cna
  ), collapse = ",")

  karyotype_string = gsub(",{2,}", ",", karyotype_string)
  karyotype_string = gsub("^,|,$", "", karyotype_string)
  
  karyotype_out = list(karyotype_string = karyotype_string, karyotype_list = karyotype_list)
  return(karyotype_out)
}


