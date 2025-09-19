#' @title process_gencode
#' @description
#'
#' Helper script to process gencode parameter
#'
#' @param gencode path to gencode file. Gencode file must be either rds or some format accepted by rtracklayer::import (e.g. GTF) with just a single entry for each gene (so gencode entries for each gene are collapse to a single range). The input could be .gtf or .rds with GRanges object, or a GRanges object i.e. resulting from importing the (appropriate) GENCODE .gtf via rtracklayer, note: this input is only used in CNA to gene mapping.
#' @return gencode_gr GRanges
#' @author Marcin Imielinski
#' @export
process_gencode <- function(gencode = NULL, assembly = NULL, seqlevelsstyle = "NCBI") {
  is_null = is.null(gencode)
  
  if (is_null) {
    gencode = get_default_gencode(assembly = assembly)
    is_null = is.null(gencode)
  }
  
  if (is_null) stop("No default gencode found, file must be provided")

  is_character = is.character(gencode)
  is_len_one = NROW(gencode) == 1
  is_not_valid = is_character && ! NROW(gencode) == 1
  is_na = is_len_one && (is.na(gencode) || gencode %in% c("NA", base::nullfile()))
  is_possible_path = is_character && is_len_one && !is_na
  is_existent_path = is_possible_path && file.exists(gencode)
  is_rds = is_possible_path && grepl(".rds$", gencode)
  is_gff = is_possible_path && grepl(".gtf(.gz){0,}$|.gff([0-9]){0,}(.gz){0,}$", gencode)
  
  if (is_existent_path && is_rds) {
    gencode <- readRDS(gencode)
  } else if (is_existent_path && is_gff) {
    gencode <- rtracklayer::import(gencode)
    
  } else if (is_existent_path) {
    gencode = data.table::fread(gencode)
    gencode = tryCatch(
      gUtils::dt2gr(gencode),
      error = function(e) "ERR"
    )
    if (identical(gencode, "ERR")) {
      gencode = tryCatch(
        as(gencode, "GRanges"),
        error = function(e) NULL
      )
    }
  } else if (is_character && is_len_one) {
    stop("Path provided does not exist": gencode)
  } else if (is_not_valid) {
    stop("Path provided must be a length one string")
  }

  is_granges = inherits(gencode, "GRanges")
  if (!is_granges) stop("gencode must be read in as, provided as, or coercible to a GRanges")

  GenomeInfoDb::seqlevelsStyle(gencode) = seqlevelsstyle

  return(gencode)
}

#' Default gencode paths
#' 
#' Provided as an exported variable
#' 
#' Providing gencode defaults as a variable exposed to user
#' @export 
GENCODE_DEFAULTS = list(
  hg19 = "~/DB/GENCODE/gencode.v19.annotation.gtf.nochr.rds",
  hg38 = "~/DB/GENCODE/hg38/v29/gencode.v29.annotation.nochr.rds"
)

#' Get default gencode
#' 
#' Function to grab default gencode
#' 
#' Use environment variables
#' GENCODE_PATH
#' DEFAULT_ASSEMBLY
#' To get appropriate gencode and
#' assign it to a package level variable.
#' @export 
get_default_gencode = function(gencode_path_defaults = Skilift::GENCODE_DEFAULTS, assembly = NULL) {
  gencode_path_env = Sys.getenv("GENCODE_PATH")
  gencode_dir_env = Sys.getenv("GENCODE_DIR")
  
  is_assembly_arg_provided = !is.null(assembly)
  is_assembly_arg_character = is.character(assembly)
  is_assembly_arg_legit = (
    is_assembly_arg_provided
    && is_assembly_arg_character 
    && NROW(assembly) == 1 
  )
  if (is_assembly_arg_legit) assembly = tolower(assembly)
  is_assembly_arg_valid = is_assembly_arg_legit && all(assembly %in% c("hg19", "hg38", "grch37", "grch38"))
  default_assembly = assembly
  if (!is_assembly_arg_valid) {
    default_assembly = tolower(Sys.getenv("DEFAULT_ASSEMBLY"))
  }
  
  is_assembly_valid = all(default_assembly %in% c("", "hg19", "hg38", "grch37", "grch38"))
  is_assembly_provided = nzchar(default_assembly)
  
  if (!is_assembly_valid) {
    stop("Assembly must be one of: hg19, hg38, GRCh37, or GRCh38")
  }

  if (!is_assembly_provided) {
    message("No default assembly found, using GRCh37/hg19")
    default_assembly = "hg19"
  }

  remap_assembly = list(
    c("grch37", "hg19"),
    c("grch38", "hg38")
  )

  for (x in remap_assembly) {
    default_assembly[default_assembly == x[1]] = x[2]
  }
  
	is_gencode_path_provided = nzchar(gencode_path_env)
	is_gencode_path_existent = file.exists(gencode_path_env)

  is_gencode_dir_provided = nzchar(gencode_dir_env)
	is_gencode_dir_existent = dir.exists(gencode_dir_env)

  gencode_paths_to_parse = gencode_path_defaults[[default_assembly]]

	is_gencode_default_existent_locally = file.exists(gencode_paths_to_parse)
	gencode_paths = gencode_paths_to_parse[is_gencode_default_existent_locally]
	is_any_gencode_default_existent_locally = NROW(gencode_paths) > 0

	output_gencode_path = NULL

	if (is_gencode_path_provided && is_gencode_path_existent) {
		output_gencode_path = gencode_path_env
		message("GENCODE_PATH environment variable found: ", output_gencode_path)
	}

  is_output_null = is.null(output_gencode_path)
  
  if (is_output_null && is_gencode_dir_provided && is_gencode_dir_existent) {
    output_gencode_path = file.path(gencode_dir_env, basename(gencode_paths_to_parse))
    is_gencode_dir_path_existent = file.exists(output_gencode_path)
    if (!is_gencode_dir_path_existent) {
      message("GENCODE_DIR provided, but default gencode path could not be found: ", output_gencode_path)
      output_gencode_path = NULL
    } else {
      message("GENCODE_DIR provided, and default path found: ", output_gencode_path)
    }
  }

  is_output_null = is.null(output_gencode_path)

  if (is_output_null && is_any_gencode_default_existent_locally) {
		output_gencode_path = gencode_paths[1]
		message("Using gencode path as default: ", output_gencode_path)
	}

  return(output_gencode_path)
}

#' @title process_cytoband
#' @description
#'
#' Helper script to process cytoband
#'
#' @param cytoband path to cytoband file provided by UCSC
#' @return GenomicRanges::GRanges object
#' @author Kevin Hadi
process_cytoband <- function(cyto = NULL, coarse = FALSE) {
  if (is.null(cyto)) {
    cyto <- system.file("extdata", "data", "cytoband.rds", package = "Skilift")
    warning("Using default cytoband, hg19")
    cyto <- readRDS(cyto)
    return(cyto)
  } else if (is.character(cyto)) {
    if (grepl(".rds$", cyto)) {
      cyto <- readRDS(cyto)
      return(cyto)
    } else {
      cyto <- data.table::fread(cyto)
    }
  }

  names(cyto) <- c("seqnames", "start", "end", "band", "stain")
  isZeroStart <- any(cyto[, length(intersect(start, end)), by = seqnames]$V1 > 0) || any(cyto$start == 0)
  if (isZeroStart) cyto$start <- cyto$start + 1
  cyto <- gUtils::dt2gr(cyto)
  GenomeInfoDb::seqlevelsStyle(cyto) <- "NCBI"
  cyto$chrom_name <- as.character(seqnames(cyto))
  pasteband <- cyto$band
  if (coarse) {
    pasteband <- gsub("\\.[0-9]+", "", cyto$band)
  }
  # cyto$rough_band = gsub("\\.[0-9]+", "", cyto$band)
  cyto$chromband <- paste(
    cyto$chrom_name,
    pasteband,
    sep = ""
  )
  return(cyto)
}

eNROW = function (x) {
    return(vapply(x, NROW, integer(1)))
}

coalesce = function (..., bads = NA, bads.x = NULL, bads.y = NULL, r2l = FALSE, 
    fromLast = FALSE, opposite = TRUE, comparefun = NULL, remove.empty = TRUE) 
{
    lst = list(...)
    lens = eNROW(lst)
    maxlen = max(lens)
    if (length(unique(lens)) > 1) 
        lst = lapply(lst, function(x) rep(x, length.out = maxlen))
    if (remove.empty) 
        lst = lst[eNROW(lst) > 0]
    if (!length(bads) && !length(bads.x) && !length(bads.y)) 
        stop("You gotta set one of bads, bads.x, or bads.y")
    if ({
        length(bads.x) && length(bads.y)
    }) {
        message("bads.x and bads.y both set explicitly")
        message("setting opposite to FALSE")
        opposite = FALSE
    }
    anytrue = function(vec) rep(TRUE, length.out = length(vec))
    if (isTRUE(bads) || !length(bads)) {
        message("bads set to NULL or TRUE")
        message("setting opposite to FALSE")
        bads = anytrue
        opposite = FALSE
    }
    if (opposite) {
        yfun = get("!", mode = "function")
    }
    else {
        yfun = get("identity", mode = "function")
    }
    if (!length(bads.x)) 
        bads.x = bads
    if (!length(bads.y)) 
        bads.y = bads
    dofun = function(x, y) {
        if (is.function(bads.x)) 
            badsx = which(bads.x(x))
        else badsx = which(x %in% bads.x)
        if (is.function(bads.y)) 
            nbadsy = which(yfun(bads.y(y)))
        else nbadsy = which(yfun(y %in% bads.y))
        ix = intersect(badsx, nbadsy)
        return(replace(x, ix, rep(y[ix], length.out = length(ix))))
    }
    if (is.null(comparefun)) {
        if (!r2l) 
            return(Reduce(function(x, y) dofun(x, y), lst, right = fromLast))
        else return(Reduce(function(x, y) dofun(y, x), lst, right = fromLast))
    }
    else {
        yfun = get("identity", mode = "function")
        if (!r2l) {
            return(Reduce(function(x, y) {
                if (is.function(bads.x)) badsx = which(bads.x(x)) else badsx = which(x %in% 
                  bads.x)
                if (is.function(bads.y)) nbadsy = which(yfun(bads.y(y))) else nbadsy = which(yfun(y %in% 
                  bads.y))
                lg = which(comparefun(x, y))
                lg = setdiff(lg, nbadsy)
                out = x
                out[badsx] = y[badsx]
                out[lg] = y[lg]
                out
            }, lst, right = fromLast))
        }
        else {
            return(Reduce(function(x, y) {
                if (is.function(bads.x)) badsx = which(bads.x(x)) else badsx = which(x %in% 
                  bads.x)
                if (is.function(bads.y)) nbadsy = which(yfun(bads.y(y))) else nbadsy = which(yfun(y %in% 
                  bads.y))
                lg = which(comparefun(x, y))
                lg = setdiff(lg, nbadsy)
                out = y
                out[nbadsy] = x[nbadsy]
                out[lg] = x[lg]
                out
            }, lst, right = fromLast))
        }
    }
}

create_cytomap = function() {
  cytos = Skilift:::process_cytoband()

  digits2 = gsub("([pq][[:alnum:]]{2,}).[[:alnum:]]+", "\\1", cytos$band)

  digits3 = gsub("([pq][[:alnum:]]{2,}.[[:alnum:]]{1,1})[[:alnum:]]+", "\\1", cytos$band, perl = TRUE)
  digits3 = ifelse(nchar(digits3) != 5, NA_character_, digits3)

  digits4 = gsub("([pq][[:alnum:]]{2,}.[[:alnum:]]{2,2})[[:alnum:]]+", "\\1", cytos$band, perl = TRUE)
  digits4 = ifelse(nchar(digits4) != 6, NA_character_, digits4)
  digits5 = gsub("([pq][[:alnum:]]{2,}.[[:alnum:]]{3,3})[[:alnum:]]+", "\\1", cytos$band,  perl = TRUE)
  digits5 = ifelse(nchar(digits5) != 7, NA_character_, digits5)

  ## digits6 = gsub("([pq][[:alnum:]]{2,}.[[:alnum:]]{4,})[[:alnum:]]+", "\\1", cytos$band)
  ## digits6 = ifelse(nchar(digits6) != 8, NA_character_, digits6)


  cytos$digits2 = digits2
  cytos$digits3 = digits3
  cytos$digits4 = digits4
  cytos$digits5 = digits5

  minres_band23 = Skilift:::coalesce(cytos$digits3, cytos$digits2)
  cytos$terband = minres_band23

  cytos$arm = gsub("(p|q).*", "\\1", cytos$band)

  grcyto = (
    GenomicRanges::reduce(gUtils::gr_construct_by(cytos, c("chrom_name", "arm", "terband")))
    %>% gUtils::gr_deconstruct_by(c("chrom_name", "arm", "terband"), meta = TRUE)
  )
  GenomeInfoDb::seqlevelsStyle(grcyto) = "NCBI"
  grcyto = GenomeInfoDb::sortSeqlevels(grcyto)
  grcyto = GenomicRanges::sort(grcyto, ignore.strand = FALSE)

  dtcytoter = (
    grcyto
    %>% as.data.frame()
    %>% data.table::setDT()
  )


  dtcytoter = dtcytoter[, .SD[c(1, .N)], by = seqnames]


  gr_chr = (
    GenomicRanges::reduce(gUtils::gr_construct_by(cytos, c("chrom_name")))
    %>% gUtils::gr_deconstruct_by(c("chrom_name"), meta = TRUE)
  )
  GenomeInfoDb::seqlevelsStyle(gr_chr) = "NCBI"
  gr_chr = GenomeInfoDb::sortSeqlevels(gr_chr)
  gr_chr = GenomicRanges::sort(gr_chr, ignore.strand = FALSE)

  dtchr = (
    gr_chr
    %>% as.data.frame()
    %>% data.table::setDT()
  )


  ranges_2 = cytos[!is.na(cytos$digits2)]
  ranges_2$band = ranges_2$digits2

  ranges_3 = cytos[!is.na(cytos$digits3)]
  ranges_3$band = ranges_3$digits3

  ranges_4 = cytos[!is.na(cytos$digits4)]
  ranges_4$band = ranges_4$digits4

  ranges_5 = cytos[!is.na(cytos$digits5)]
  ranges_5$band = ranges_5$digits5

  ranges_arm = cytos[!is.na(cytos$arm)]
  ranges_arm$band = ranges_arm$arm

  ranges_ter = dtcytoter[, .(seqnames, start, end, band = ifelse(arm == "p", "pter", "qter"))] %>% as("GRanges")

  ranges_chr = dtchr[, .(seqnames, start, end, band = "chrom")] %>% as("GRanges")

  allcytolevels = c(
    ranges_2,
    ranges_3,
    ranges_4,
    ranges_5,
    ranges_arm,
    ranges_ter,
    ranges_chr
  )


  cytos_by_level = GenomicRanges::reduce(
    gUtils::gr_construct_by(allcytolevels, "band"),
    with.revmap = TRUE
  )

  allcytolevels = gUtils::gr_deconstruct_by(cytos_by_level, meta = TRUE, by = "band")

  cytomap = setDT(as.data.frame(allcytolevels))[]
  cytomap$chrom = as.character(cytomap$seqnames)
  setkeyv(cytomap, c("chrom", "band"))

  cytomap$category = data.table::fcase(
    cytomap$band == "chrom", "chromosome",
    cytomap$band %in% c("p", "q"), "arm",
    default = "band"
  )
  cytomap$category_width = as.numeric(cytomap$width)
  return(cytomap)
    

  ## saveRDS(cytomap, "~/DB/UCSC/hg19.cytoband.map.rds")
  ## normalizePath("~/DB/UCSC/hg19.cytoband.map.rds")
  ## "/gpfs/data/imielinskilab/DB/UCSC/hg19.cytoband.map.rds"

}



#' Global HemeDB path
HEMEDB = function() {
  system.file("extdata", "data", "hemedb.rds", package = "Skilift")
} 

#' Duncavage DB 
#' 
#' Duncavage et al. Blood 2022.
DUNCAVAGEDB = function() {
  system.file("extdata", "data", "Duncavage_Blood_22.rds", package = "Skilift")
}

#' COSMIC Fusions
#' 
#' https://cancer.sanger.ac.uk/cosmic/download/cosmic/v101/fusion
#' 
COSMIC_FUSIONS = function() {
  system.file("extdata", "data", "Cosmic_Fusion_v101_GRCh37.tsv.gz", package = "Skilift")
}

#' COSMIC fusions
get_curated_fusions <- function() {

  fusions_recurrent = character(0)
  fusions_recurrent = c(
    fusions_recurrent,
    readRDS(Skilift:::DUNCAVAGEDB())[grepl("::", Gene)]$Gene
  )

  fusions_cosmic = fread(Skilift:::COSMIC_FUSIONS())
  fusions_cosmic[, fusion_gene := paste(FIVE_PRIME_GENE_SYMBOL, "::", THREE_PRIME_GENE_SYMBOL, sep = "")]

  fusions_recurrent = c(
    fusions_recurrent,
    (
      fusions_cosmic[, .(
        num_samples = length(unique(COSMIC_SAMPLE_ID))
      ), by = fusion_gene]
      [num_samples > 1]
    )$fusion_gene
  )
  fusions_recurrent = unique(fusions_recurrent)
  
  ## Blacklist are fusions that absolutely should not
  ## appear. AFF1::KMT2A and RUNX1::ETV6 are both present
  ## as recurrent COSMIC fusions but are not themselves
  ## the drivers.
  blacklist = c(
    "AFF1::KMT2A",
    "RUNX1::ETV6"
  )

  return(
    list(
      fusions_recurrent = fusions_recurrent,
      blacklist = blacklist
    )
  )
}




#' @title collect_gene_fusions
#' @description
#' Collects gene fusion data from a specified file and processes it.
#'
#' @param fusions Path to the fusions.rds file.
#' @param pge GRanges object with gencode annotations.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed gene fusion information.
collect_gene_fusions <- function(fusions, pge, verbose = TRUE) {
  if (is.null(fusions) || !file.exists(fusions)) {
    if (verbose) message("Fusions file is missing or does not exist.")
    return(data.table(type = NA, source = "fusions"))
  }

  if (verbose) message("pulling fusions")
  fus <- readRDS(fusions)$meta

  if (nrow(fus) == 0) {
    if (verbose) message("No fusions found in the file.")
    return(data.table(type = NA, source = "fusions"))
  }

  non_silent_fusions <- fus[silent == FALSE, ]
  unique_fusions <- non_silent_fusions[!duplicated(genes), ]
  unique_fusions[, vartype := ifelse(in.frame == TRUE, "fusion", "outframe_fusion")]

  split_genes <- function(genes) unlist(strsplit(genes, ","))
  gene_lengths <- function(genes) sapply(strsplit(genes, ","), length)

  fus <- unique_fusions[, .(
    gene = split_genes(genes),
    vartype = rep(vartype, gene_lengths(genes)),
    fusion_genes = rep(genes, gene_lengths(genes))
  )][, `:=`(track = "variants", type = vartype, source = "fusions")]

  get_gene_coords <- function(genes) {
    coords <- lapply(genes, function(gene) {
      gene_ranges <- pge[mcols(pge)$gene_name == gene]
      paste0(seqnames(gene_ranges), ":", start(gene_ranges), "-", end(gene_ranges))
    })
    paste(unlist(coords), collapse = ",")
  }

  fus[, fusion_gene_coords := unlist(lapply(strsplit(fusion_genes, ","), get_gene_coords))]

  return(fus)
}
#' @description
#' Collects complex events from a specified file and processes them.
#'
#' @param complex Path to the complex.rds file.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed complex event information.
collect_complex_events <- function(complex, verbose = TRUE) {
  if (is.null(complex) || !file.exists(complex)) {
    if (verbose) message("Complex events file is missing or does not exist.")
    return(data.table(type = NA, source = "complex"))
  }

  if (verbose) message("pulling complex events")
  sv <- readRDS(complex)$meta$events

  if (nrow(sv) == 0) {
    if (verbose) message("No complex events found in the file.")
    return(data.table(type = NA, source = "complex"))
  }

  sv_summary <- sv[, .(value = .N), by = type]
  simple_sv_types <- c("del", "dup", "invdup", "tra", "inv")
  sv_summary[, track := ifelse(type %in% simple_sv_types, "simple sv", "complex sv")]
  sv_summary[, source := "complex"]

  return(sv_summary)
}

check_GRanges_compatibility = function (gr1, gr2, name1 = "first", name2 = "second") 
{
    non_overlapping_seqnames1 = setdiff(seqlevels(gr1), seqlevels(gr2))
    non_overlapping_seqnames2 = setdiff(seqlevels(gr2), seqlevels(gr1))
    overlap = intersect(seqlevels(gr1), seqlevels(gr2))
    message("The following seqnames are only in the ", name1, 
        " GRanges, but not in the ", name2, " GRanges: ", paste(non_overlapping_seqnames1, 
            collapse = ", "))
    message("The following seqnames are only in the ", name2, 
        " GRanges, but not in the ", name1, " GRanges: ", paste(non_overlapping_seqnames2, 
            collapse = ", "))
    message("The follosing seqnames are in both GRanges objects: ", 
        paste(overlap, collapse = ", "))
    if (length(non_overlapping_seqnames1) > 0 | length(non_overlapping_seqnames2) > 
        0) {
        return(FALSE)
    }
    return(TRUE)
}

get_gene_copy_numbers <- function(
  gg, 
  gene_ranges, 
  nseg = NULL, 
  gene_id_col = "gene_id", 
  simplify_seqnames = FALSE, 
  mfields = c("gene_name", "source", "gene_id", "gene_type", "level", "hgnc_id", "havana_gene"), 
  output_type = "data.table",
  min_width = 1e3,
  min_cn_quantile_threshold = 0.1,
  max_cn_quantile_threshold = 0.9
) {
    if (is.character(gg)) {
        gg = readRDS(gg)
    }
    if (!inherits(gene_ranges, "GRanges")) {
        gene_ranges = rtracklayer::import(gene_ranges)
    }
    if (!(output_type %in% c("GRanges", "data.table"))) {
        stop("Invalid output_type: ", output_type, ". outputtype must be either \"GRanges\" or \"data.table\".")
    }
    ngr = gg$nodes$gr
    if (simplify_seqnames) {
        ngr = gr.sub(ngr)
        gene_ranges = gr.sub(gene_ranges)
    }
    GRanges_are_compatible = check_GRanges_compatibility(ngr, 
        gene_ranges, "gGraph segments", "genes")
    if (!is.null(nseg)) {
        ngr = ngr %$% nseg[, c("ncn")]
    }
    else {
        message("No normal copy number segmentation was provided so assuming CN = 2 for all seqnames.")
        ngr$ncn = 2
    }
    ndt = gr2dt(ngr)
    seq_widths = as.numeric(width(ngr))
    normal_ploidy = round(sum(seq_widths * ngr$ncn, na.rm = T)/sum(seq_widths, 
        na.rm = T))
    ndt[, `:=`(normalized_cn, cn * normal_ploidy/(gg$meta$ploidy * 
                                                    ncn))]
    gene_ranges$gene_width = width(gene_ranges)
    gene_cn_segments = dt2gr(ndt, seqlengths = seqlengths(gg)) %*% 
      gene_ranges %>% gr2dt
    
    # gene_cn_stats = gene_cn_segments[, .(
    #   max_normalized_cn = max(normalized_cn, na.rm = TRUE),
    #   max_cn = max(cn, na.rm = TRUE),
    #   min_normalized_cn = min(normalized_cn, na.rm = TRUE),
    #   min_cn = min(cn, na.rm = TRUE),
    #   avg_normalized_cn = sum(normalized_cn * width, na.rm = TRUE) / sum(width),
    #   avg_cn = sum(cn * width, na.rm = TRUE) / sum(width),
    #   total_node_width = sum(width, na.rm = TRUE),
    #   number_of_cn_segments = .N,
    #   ncn = ncn[1],
    #   list_of_segs = list(.SD[, list(gene_name = gene_name[1], normalized_cn, cn, width)])
    #   ## gene_width = gene_width[1]
    # ), by = gene_name]

    gene_cn_segments[, ix := seq_len(.N)]

	# Order by copy number 
	reord_gene_cn_segments = data.table::copy(gene_cn_segments[order(cn)])
	reord_gene_cn_segments[
		,
		c("weight", "total_weight", "cweight", "from_cfrac", "to_cfrac", "is_at_min_quantile_threshold", "is_at_max_quantile_threshold") := {
			weight = width # Should the weight just be width? This would make it a true quantile
			# weight = width * (cn + 1e-9) # take care of 0's with a fudge factor
			total_weight = sum(weight)
			cweight = cumsum(weight)
			to_cfrac = cweight / total_weight
            ## setting lowest to something ridiculous to make sure argument of 0 works
			from_cfrac = c(-1e9, to_cfrac[-.N])
            ## setting highest to something ridiculous to make sure argument of 1 works
            to_cfrac[.N] = 1e9
			## interval is semi-closed - (from_cfrac, to_cfrac] (inclusive of to_cfrac, but not from_cfrac)
			## so any interval included where from_cfrac is greater than or equal to threshold should be excluded
			is_at_min_quantile_threshold = (
				data.table::between(min_cn_quantile_threshold, from_cfrac, to_cfrac) 
				& !from_cfrac >= min_cn_quantile_threshold
			)
			is_at_max_quantile_threshold = (
				data.table::between(max_cn_quantile_threshold, from_cfrac, to_cfrac) 
				& !from_cfrac >= max_cn_quantile_threshold
			)
			list(weight, total_weight, cweight, from_cfrac, to_cfrac, is_at_min_quantile_threshold, is_at_max_quantile_threshold)
		}
		,
		by = gene_name
	]
    
	gene_cn_segments = reord_gene_cn_segments[order(ix)]

	null_out_columns = c("ix")
	for (col in null_out_columns) {
		gene_cn_segments[[col]] = NULL
	}


    gene_cn_table = gene_cn_segments[, `:=`(
      max_normalized_cn = max(normalized_cn, na.rm = TRUE),
      max_cn = max(cn, na.rm = TRUE),
	  max_quantile_cn = cn[is_at_max_quantile_threshold],
	  max_quantile_normalized_cn = normalized_cn[is_at_max_quantile_threshold],
      min_normalized_cn = min(normalized_cn, na.rm = TRUE),
      min_cn = min(cn, na.rm = TRUE),
	  min_quantile_cn = cn[is_at_min_quantile_threshold],
	  min_quantile_normalized_cn = normalized_cn[is_at_min_quantile_threshold],
      avg_normalized_cn = sum(normalized_cn * width, na.rm = TRUE) / sum(width),
      avg_cn = sum(cn * width, na.rm = TRUE) / sum(width),
      # total_node_width = sum(width, na.rm = TRUE),
      number_of_cn_segments = .N,
      ncn = ncn[1],
      # list_of_segs = list(.SD[, list(gene_name = gene_name[1], normalized_cn, cn, width)])
      gene_width = gene_width[1]
    ), by = gene_name]
    
    if (output_type == "data.table") {
        return(gene_cn_table)
    }
    return(dt2gr(gene_cn_table, seqlengths = seqlengths(gene_ranges)))
}

#' Get gene amplifications and deletions from jabba object
#'
#' Copied from github::mskilab-org/skitools.
#' This function takes in a jabba object and returns
#' amplifications and deletions
#'
#' @param jab character path to jabba file or gGraph object
get_gene_ampdels_from_jabba <- function(jab, pge, amp.thresh = 4, del.thresh = 0.5, nseg = NULL, min_cn_quantile_threshold = 0.1, max_cn_quantile_threshold = 0.9) {
  gg <- jab
  if (!inherits(gg, "gGraph")) {
    gg <- gG(jabba = jab)
  }
  gene_CN <- Skilift:::get_gene_copy_numbers(
	gg, 
	gene_ranges = pge, 
	nseg = nseg, 
	min_cn_quantile_threshold = min_cn_quantile_threshold, 
	max_cn_quantile_threshold = max_cn_quantile_threshold
)
  gene_CN[, `:=`(type, NA_character_)]
  gencode = dynGet("gencode")
  exons_merged = GenomicRanges::reduce(gUtils::gr_construct_by(gencode[gencode$type == "exon"], by = "gene_name"))
  ## exons_merged = gUtils::gr_deconstruct_by(exons_merged, meta = TRUE, "gene_name")
  exons_merged$exon_width = width(exons_merged)
  gr_gene_cn = dt2gr(gene_CN)
  gr_gene_cn = gUtils::gr_construct_by(gr_gene_cn, by = "gene_name")
  exons_merged = dt2gr(gr2dt(exons_merged)[, total_exon_width := sum(exon_width), by = seqnames])
  ## gr_gene_cn = gUtils::gr.val(gr_gene_cn, exons_merged, val = c("exon_width"), mean = FALSE)
  exon_cn = gr2dt(gr_gene_cn %*% exons_merged)
  exon_cn_map = exon_cn[, .(overlapped_exon_width = sum(width), total_exon_width = total_exon_width[1]), keyby = query.id]
  exon_cn_map = exon_cn_map[list(seq_len(NROW(gene_CN)))]
  #   gene_CN[min_normalized_cn >= amp.thresh, `:=`(type, "amp")]

  gene_CN$overlapped_exon_width = exon_cn_map$overlapped_exon_width
  gene_CN$total_exon_width = exon_cn_map$total_exon_width

  
  gene_CN[min_quantile_normalized_cn >= amp.thresh & cn >= amp.thresh, `:=`(type, "amp")]
  gene_CN[min_cn > 1 & cn > 1 & min_normalized_cn < del.thresh, `:=`(
    type,
    "del"
    )]

  gene_CN[min_cn == 1 & cn == min_cn & min_cn < ncn, `:=`(type, "hetdel")]
  gene_CN[min_cn == 0 & cn == min_cn, `:=`(type, "homdel")]

  gene_CN[type == "amp", min_cn := min_quantile_cn]

  # scna_result = gene_CN[!is.na(type)]
  scna_result = (
      gene_CN[, 
      .(
        max_normalized_cn = max_normalized_cn[1],
        max_cn = max_cn[1],
		max_quantile_cn = max_quantile_cn[1],
		max_quantile_normalized_cn = max_quantile_normalized_cn[1],
        min_normalized_cn = min_normalized_cn[1],
        min_cn = min_cn[1],
		min_quantile_cn = min_quantile_cn[1],
		min_quantile_normalized_cn = min_quantile_normalized_cn[1],
        avg_normalized_cn = avg_normalized_cn[1],
        avg_cn = avg_cn[1],
        number_of_cn_segments = number_of_cn_segments[1],
        gene_width = gene_width[1],
        total_node_width = sum(width), # This must be calculated after nominating SCNA type, not before.,
        exon_frac = sum(overlapped_exon_width / total_exon_width, na.rm = TRUE)
      ), 
      by = .(gene_name, type)]
  )

  scna_result = base::subset(
    scna_result,
    !is.na(scna_result$type)
  )

  return(scna_result)
}

#' @title collect_copy_number_jabba
#' @description
#' Collects copy number and jabba data from a specified file and processes it.
#'
#' @param jabba_rds Path to the jabba.simple.rds file.
#' @param pge GRanges object with gencode annotations.
#' @param amp.thresh SCNA amplification threshold.
#' @param del.thresh SCNA deletion threshold.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @param karyograph Optional path to the karyograph.rds file.
#' @return A data.table containing processed copy number and jabba information.
collect_copy_number_jabba <- function(
    jabba_rds,
    pge,
    amp.thresh,
    del.thresh,
    verbose = TRUE,
    karyograph = NULL,
	min_cn_quantile_threshold = 0.1,
	max_cn_quantile_threshold = 0.9) {
  if (is.null(jabba_rds) || !file.exists(jabba_rds)) {
    if (verbose) message("Jabba RDS file is missing or does not exist.")
    return(data.table(type = NA, source = "jabba_rds"))
  }

  if (verbose) message("pulling jabba_rds to get SCNA and purity / ploidy")
  jab <- readRDS(jabba_rds)
  jabpurity = base::get0("purity", as.environment(jab$meta), ifnotfound = NULL)
  if (is.null(jabpurity)) jabpurity = base::get("purity", jab) ## will error out if not found
  jabploidy = base::get0("ploidy", as.environment(jab$meta), ifnotfound = NULL)
  if (is.null(jabploidy)) jabploidy = base::get("ploidy", jab) ## will error out if not found

  result <- data.table(
    value = c(jabpurity, jabploidy),
    type = c("purity", "ploidy"),
    track = "pp"
  )

  # get the ncn data from jabba
  nseg <- NULL

  # Don't fail out if karyograph isn't found
  # Just provide nseg as NULL, and get_gene_ampdels_from_jabba
  # assumes ncn = 2 (pretty safe assumption)
  if (!is.null(karyograph) && file.exists(karyograph)) {
    nseg <- readRDS(karyograph)$segstats[, c("ncn")]
  }

  scna <- get_gene_ampdels_from_jabba(
    jab,
    amp.thresh = amp.thresh,
    del.thresh = del.thresh,
    pge = pge,
    nseg = nseg,
	min_cn_quantile_threshold = min_cn_quantile_threshold,
	max_cn_quantile_threshold = max_cn_quantile_threshold
  )

  if (nrow(scna)) {
    scna[, track := "variants"][, source := "jabba_rds"][, vartype := "scna"]
    result <- rbind(
      result,
      scna[, .(value = min_cn, type, track, gene = gene_name)],
      fill = TRUE,
      use.names = TRUE
    )
  }

  return(result)
}

#' @description
#' Collects signature data from a specified file and processes it.
#'
#' @param signature_counts Path to the signature_counts.txt file.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed signature information.
collect_signatures <- function(signature_counts, verbose = TRUE) {
  # TODO: update this function to use sigprofiler
  if (!is.null(signature_counts) && file.exists(signature_counts)) {
    if (verbose) message("pulling signature_counts")
    sig <- fread(signature_counts)
    sig <- sig[, .(value = num_events, type = Signature, etiology = Etiology, frac = frac.events, track = "signature", source = "signature_counts")]
    return(sig)
  } else {
    return(data.table(type = NA, source = "signature_counts"))
  }
}

#' @title collect_gene_mutations
#' @description
#' Collects gene mutation data from a specified file and processes it.
#'
#' @param annotated_bcf Path to the annotated.bcf file.
#' @param jabba_rds Path to the jabba.simple.rds file.
#' @param filter Filter to apply to the variants.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed gene mutation information.
collect_gene_mutations <- function(
    annotated_bcf,
    jabba_rds,
    filter,
    verbose = TRUE) {
  if (is.null(annotated_bcf) || !file.exists(annotated_bcf)) {
    if (verbose) message("Annotated BCF file is missing or does not exist.")
    return(data.table(type = NA, source = "annotated_bcf"))
  }

  if (verbose) message("pulling annotated_bcf using FILTER=", filter)

  local_bcftools_path <- Sys.which("bcftools")
  if (local_bcftools_path == "") {
    stop("bcftools not found in the system PATH. Please install or module load bcftools.")
  }
  if (verbose) message("bcftools found at: ", local_bcftools_path)

  bcf <- skitools::grok_bcf(
    annotated_bcf,
    label = "sample",
    long = TRUE,
    filter = filter,
    bpath = local_bcftools_path
  )
  if (verbose) message(length(bcf), " variants pass filter")

  genome.size <- sum(GenomeInfoDb::seqlengths(bcf), na.rm = TRUE) / 1e6
  if (is.na(genome.size)) {
    genome.size <- sum(
      GenomeInfoDb::seqlengths(gGnome::gG(jabba = jabba_rds)),
      na.rm = TRUE
    ) / 1e6
  }

  nmut <- unique(data.table(
    as.character(seqnames(bcf)),
    start(bcf),
    end(bcf),
    bcf$REF,
    bcf$ALT
  ))[, .N]
  mut.density <- data.table(
    value = c(nmut, nmut / genome.size),
    type = c("count", "density"),
    track = "tmb",
    source = "annotated_bcf"
  )

  keepeff <- c("trunc", "cnadel", "cnadup", "complexsv", "splice", "inframe_indel", "fusion", "missense", "promoter", "regulatory", "mir")
  bcf <- bcf[bcf$short %in% keepeff]
  if (verbose) message(length(bcf), " variants pass keepeff")

  if (length(bcf) == 0) {
    return(mut.density)
  }

  bcf$variant.g <- paste0(seqnames(bcf), ":", start(bcf), "-", end(bcf), " ", bcf$REF, ">", bcf$ALT)
  vars <- gr2dt(bcf)[, .(gene, vartype, variant.g, variant.p, distance, annotation, type = short, track = "variants", source = "annotated_bcf")]
  setkey(vars, variant.g)
  vars <- vars[, .SD[1], by = variant.g]

  return(rbind(mut.density, vars, fill = TRUE, use.names = TRUE))
}

#' @title collect_oncokb_cna
#' @description
#' Collects OncoKB CNA data from a specified file and processes it.
#'
#' @param oncokb_cna Path to the oncokb CNA file.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed OncoKB CNA information.
collect_oncokb_cna <- function(oncokb_cna, jabba_gg, pge, amp.thresh, del.thresh, karyograph = NULL, verbose = TRUE, min_cn_quantile_threshold = 0.1, max_cn_quantile_threshold = 0.9) {
  if (is.null(oncokb_cna) || !file.exists(oncokb_cna)) {
    if (verbose) message("OncoKB CNA file is missing or does not exist.")
    return(data.table(type = NA, source = "oncokb_cna"))
  }

  oncokb_cna <- data.table::fread(oncokb_cna)

  nseg = NULL
  if (!is.null(karyograph) && file.exists(karyograph)) {
      nseg <- readRDS(karyograph)$segstats[, c("ncn")]
  }

  scna <- get_gene_ampdels_from_jabba(
      jabba_gg,
      amp.thresh = amp.thresh,
      del.thresh = del.thresh,
      pge = pge,
      nseg = nseg,
	  min_cn_quantile_threshold = 0.1, 
	  max_cn_quantile_threshold = 0.9
  ) 

  matches = list(
      c("Deletion", "homdel"),
      c("Amplification", "amp")
  )
  is_valid_oncokb_cna = rep(FALSE, NROW(oncokb_cna))
  for (match_lst in matches) {
      valid_gene = scna[
          ## scna$total_node_width > 1e3
          scna$exon_frac > 0.9
          & scna$type == match_lst[2],

      ]$gene_name
      is_valid_oncokb_cna = (
          is_valid_oncokb_cna
          | (oncokb_cna$HUGO_SYMBOL %in% valid_gene & oncokb_cna$ALTERATION %in% match_lst[1])
      )
  }
  oncokb_cna = oncokb_cna[is_valid_oncokb_cna,]

  if (NROW(oncokb_cna) > 0) {
    oncokb_cna <- parse_oncokb_tier(
      oncokb_cna,
      tx_cols = c("LEVEL_1", "LEVEL_2"),
      rx_cols = c("LEVEL_R1"),
      dx_cols = c("LEVEL_Dx1"),
      px_cols = c("LEVEL_Px1")
    )
    return(oncokb_cna[, .(
      gene = HUGO_SYMBOL,
      gene_summary = GENE_SUMMARY,
      role = Role,
      value = min_cn,
      type = ifelse(
        ALTERATION == "Amplification",
        "amp",
        ifelse(ALTERATION == "Deletion", "homdel", NA_character_)
      ),
      tier = tier,
      tier_description = tier_factor,
      variant_summary = VARIANT_SUMMARY,
      therapeutics = tx_string,
      resistances = rx_string,
      diagnoses = dx_string,
      prognoses = px_string,
      effect = MUTATION_EFFECT,
      effect_description = MUTATION_EFFECT_DESCRIPTION,
      track = "scna",
      source = "oncokb_cna"
    )])
  }
  return(data.table(type = NA, source = "oncokb_cna"))
}


#' @title collect_oncokb_fusions
#' @description
#' Collects OncoKB Fusion data from a specified file and processes it.
#'
#' @param oncokb_fusions Path to the oncokb fusions file.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed OncoKB Fusion information.
#' @author Kevin Hadi
collect_oncokb_fusions <- function(oncokb_fusions, pge, cytoband, verbose = TRUE) {
  out <- data.table(vartype = NA, source = "oncokb_fusions")
  if (is.null(oncokb_fusions) || !file.exists(oncokb_fusions)) {
    if (verbose) message("OncoKB Fusions file is missing or does not exist.")
    return(out)
  }

  oncokb_fusions <- data.table::fread(oncokb_fusions)

  if (NROW(oncokb_fusions) > 0) {
    oncokb_fusions <- parse_oncokb_tier(
      oncokb_fusions,
      tx_cols = c("LEVEL_1", "LEVEL_2"),
      rx_cols = c("LEVEL_R1"),
      dx_cols = c("LEVEL_Dx1"),
      px_cols = c("LEVEL_Px1")
    )

    non_silent_fusions <- oncokb_fusions[silent == FALSE, ] # already de-duped
    non_silent_fusions[, vartype := ifelse(in.frame == TRUE, "fusion", "outframe_fusion")]

    if (!NROW(non_silent_fusions) > 0) {
      return(out)
    }

    genes <- strsplit(non_silent_fusions$FUSION, "-")
    genes_matrix <- do.call(rbind, genes)
    ixA <- match(genes_matrix[, 1], pge$gene_name)
    ixB <- match(genes_matrix[, 2], pge$gene_name)
    # remove NA from ixA and ixB (needed for tests that use subset of gencode genes)
    is_na_A = is.na(ixA)
    is_na_B = is.na(ixB)
    is_either_na = is_na_A | is_na_B
    na.index <- which(is_either_na)
    # ixA <- ixA[!is.na(ixA)]
    # ixB <- ixB[!is.na(ixB)]
    grA <- pge[ixA[!is_na_A]]
    grB <- pge[ixB[!is_na_B]]
    coordB = coordA = character(NROW(non_silent_fusions))
    coordA[!is_na_A] <- gUtils::gr.string(grA)
    coordB[!is_na_B] <- gUtils::gr.string(grB)

    # Leaving cytoband query code in for now, because we will want to re-incorporate
    # at a later date
    grovA <- gUtils::gr.findoverlaps(grA, cytoband, scol = "chromband")
    grovB <- gUtils::gr.findoverlaps(grB, cytoband, scol = "chromband")

    
    if (length(grovA) > 0) {
      grovA <- GenomicRanges::sort(grovA) %Q% (order(query.id, ifelse(grepl("p", chromband), -1, 1) * start))
      non_silent_fusions$cytoA <- ifelse(!1:nrow(non_silent_fusions) %in% na.index,
        gUtils::gr2dt(grovA)[, paste(unique(chromband[chromband %in% c(chromband[1], tail(chromband, 1))]), collapse = "-"), by = query.id]$V1, NA
      )
    } else {
      non_silent_fusions$cytoA <- ""
    }

    if (length(grovB) > 0) {
      grovB <- GenomicRanges::sort(grovB) %Q% (order(query.id, ifelse(grepl("p", chromband), -1, 1) * start))
      non_silent_fusions$cytoB <- ifelse(!1:nrow(non_silent_fusions) %in% na.index,
        gUtils::gr2dt(grovB)[, paste(unique(chromband[chromband %in% c(chromband[1], tail(chromband, 1))]), collapse = "-"), by = query.id]$V1, NA
      )
    } else {
      non_silent_fusions$cytoB <- ""
    }

    # non_silent_fusions$fusion_genes <- paste0(
    #   genes_matrix[, 1], "(", non_silent_fusions$exonA, ")::", 
    #   genes_matrix[, 2], "(", non_silent_fusions$exonB, ")@", 
    #   non_silent_fusions$cytoA, "::", non_silent_fusions$cytoB
    # )

    non_silent_fusions$fusion_genes <- paste0(
      genes_matrix[, 1], "::", genes_matrix[, 2]
    )
    non_silent_fusions$fusion_gene_coords <- ifelse(!1:nrow(non_silent_fusions) %in% na.index,
      paste(coordA, coordB, sep = ","),
      NA
    )

    variant.g.fus = ifelse(!1:nrow(non_silent_fusions) %in% na.index,
      as.character(glue::glue('{coordA} ({non_silent_fusions$cytoA}), {coordB} ({non_silent_fusions$cytoB})')),
      NA_character_
    )

    non_silent_fusions$variant.g = variant.g.fus


    get_queries = c("exonA", "aminoA", "exonB", "aminoB")
    query_variables = base::mget(
      get_queries,
      as.environment(as.list(non_silent_fusions)),
      ifnotfound = rep_len(
        list(rep_len(NA_character_, NROW(non_silent_fusions))),
        NROW(get_queries)
      )
    )
    exonA_label = glue::glue('Exon {query_variables$exonA}')
    aminoA_label = ifelse(is.na(query_variables$aminoA), "", glue::glue(' (p.{query_variables$aminoA})'))
    exonB_label = glue::glue('Exon {query_variables$exonB}')
    aminoB_label = ifelse(is.na(query_variables$aminoB), "", glue::glue(' (p.{query_variables$aminoB})'))

    variant.p.parsed = as.character(
      glue::glue(
        '{exonA_label}',
        '{aminoA_label}', # Note space is encoded by aminoA_label
        '::',
        '{exonB_label}',
        '{aminoB_label}' # Note space is encoded by aminoB_label
      )
    )
    variant.p.parsed = trimws(variant.p.parsed)
    variant.p.parsed = gsub("[[:space:]]{2,}", "", variant.p.parsed, perl = TRUE)

    non_silent_fusions$variant.p = variant.p.parsed
    out <- non_silent_fusions[, .(
      gene = Hugo_Symbol,
      gene_summary = GENE_SUMMARY,
      role = Role,
      variant.g,
      variant.p,
      value = min_cn,
      vartype,
      type = "fusion",
      tier = tier,
      tier_description = tier_factor,
      variant_summary = VARIANT_SUMMARY,
      therapeutics = tx_string,
      resistances = rx_string,
      diagnoses = dx_string,
      prognoses = px_string,
      effect = MUTATION_EFFECT,
      effect_description = MUTATION_EFFECT_DESCRIPTION,
      fusion_genes,
      fusion_gene_coords,
      track = "variants",
      source = "oncokb_fusions"
    )]
  }

  return(out)
}

#'
#' Helper function to parse oncokb outputs
#' levels of evidence and assign tier.
#' Tiering is simply:
#' 1 = Clinically Actionable:
#' Therapeutic sensitivity, resistance,
#' diagnostic, or prognostic information is assigned
#' to the variant.
#' 2 = Clinically Significant:
#' Clinically significant variant entails that the variant
#' does not have therapeutic, diagnostic, or prognostic effect
#' validated at FDA level, but is oncogenic.
#' 3 = All others (VUS)
#'
#' @param oncokb data.table object holding oncokb outputs
#' @author Kevin Hadi
parse_oncokb_tier <- function(
    oncokb,
    tx_cols = c("LEVEL_1", "LEVEL_2"),
    rx_cols = c("LEVEL_R1"),
    dx_cols = c("LEVEL_Dx1"),
    px_cols = c("LEVEL_Px1")) {
  .concat_string <- function(oncokb, cols) {
    out_string <- lapply(base::subset(oncokb, select = cols), function(y) (strsplit(y, ",")))
    concat_out <- IRanges::CharacterList(Reduce(concat_vectors, out_string))
    concat_out <- S4Vectors::unique(concat_out)
    concat_out <- concat_out[!is.na(concat_out)]
    concat_out[S4Vectors::elementNROWS(concat_out) == 0] <- NA_character_
    concat_out <- stringi::stri_c_list(as.list(concat_out), sep = ",")
    return(concat_out)
  }
  is_actionable <- logical(NROW(oncokb))
  for (col in c(tx_cols, rx_cols, dx_cols, px_cols)) {
    oncokb[[col]] <- as.character(oncokb[[col]])
    is_actionable <- is_actionable | (!is.na(oncokb[[col]]) & base::nzchar(oncokb[[col]]))
  }
  oncokb$is_actionable <- is_actionable
  oncokb$is_oncogenic <- oncokb$ONCOGENIC %in% c("Likely Oncogenic", "Oncogenic")
  tier_factor <- ifelse(
    oncokb$is_actionable, "Clinically Actionable",
    ifelse(oncokb$is_oncogenic, "Clinically Significant", "VUS")
  ) %>% factor(c("Clinically Actionable", "Clinically Significant", "VUS"))

  oncokb$tier_factor <- tier_factor
  oncokb$tier <- as.integer(tier_factor)
  oncokb$tx_string <- .concat_string(oncokb, tx_cols)
  oncokb$rx_string <- .concat_string(oncokb, rx_cols)
  oncokb$dx_string <- .concat_string(oncokb, dx_cols)
  oncokb$px_string <- .concat_string(oncokb, px_cols)

  return(oncokb)
}

#' @title collect_oncokb
#' @description
#' Collects OncoKB mutation data from a specified file and processes it.
#'
#' @param oncokb_maf Path to the oncokb MAF file.
#' @param multiplicity Path to the multiplicity file.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed OncoKB mutation information.
collect_oncokb <- function(oncokb_maf, multiplicity = NA_character_, verbose = TRUE) {
  
  empty_output_oncokb = data.table(type = NA, source = "oncokb_maf")
  if (missing(oncokb_maf) || is.null(oncokb_maf) || !file.exists(oncokb_maf)) {
    if (verbose) message("OncoKB MAF file is missing or does not exist.")
    return(empty_output_oncokb)
  }

  # snpeff_ontology = readRDS(system.file("extdata", "data", "snpeff_ontology.rds", package = "Skilift"))
  oncokb <- data.table::fread(oncokb_maf)
  is_multiplicity_present = is.character(multiplicity) && NROW(multiplicity) == 1 && !is_loosely_na(multiplicity) && file.exists(multiplicity)
  is_oncokb_populated = NROW(oncokb) > 0
  
  is_multiplicity_populated = FALSE
  if (is_multiplicity_present) {
    multiplicity <- readRDS(multiplicity)
	is_multiplicity_populated = NROW(multiplicity) > 0
  }

  if (is_oncokb_populated && !is_multiplicity_populated) {
  	stop("Something's off - oncokb is populated with variants, but not multiplicity.")
  }
  
  if (is_oncokb_populated && is_multiplicity_populated) {
  	oncokb <- merge_oncokb_multiplicity(oncokb, multiplicity, overwrite = TRUE)
    oncokb = Skilift:::annotate_multihit(oncokb)
  }

  if (is_oncokb_populated) {
    ## oncokb$snpeff_ontology <- snpeff_ontology$short[match(oncokb$Consequence, snpeff_ontology$eff)]
	
    oncokb$short <- dplyr::case_when(
      grepl("frameshift", oncokb$Consequence) & grepl("fs$", oncokb$HGVSp) ~ "trunc",
      grepl("stop", oncokb$Consequence) & grepl("^p\\.", oncokb$HGVSp) ~ "trunc",
      grepl("lost", oncokb$Consequence) ~ "trunc",
      grepl("missense", oncokb$Consequence) & grepl("^p\\.", oncokb$HGVSp) ~ "missense",
      grepl("splice", oncokb$Consequence) ~ "splice",
      grepl("(5|3)_prime_UTR_variant", oncokb$Consequence) ~ "UTR",
      TRUE ~ NA_character_
    )

    oncokb <- parse_oncokb_tier(
      oncokb,
      tx_cols = c("LEVEL_1", "LEVEL_2"),
      rx_cols = c("LEVEL_R1"),
      dx_cols = c("LEVEL_Dx1"),
      px_cols = c("LEVEL_Px1")
    )

    # coerce_column_tuples <- list(
    #   c("altered_copies", "total_copies")
    # )
    # current_nms <- names(oncokb)
    # for (col in coerce_column_tuples) {
    #   # if (col[1] %in% current_nms && !col[2] %in% current_nms) {
    #   current_nms[current_nms == col[1]] <- col[2]
    #   # }
    # }
    # names(oncokb) <- current_nms
    return(oncokb[, .(
      gene = Hugo_Symbol,
      gene_summary = GENE_SUMMARY,
      role = Role,
      variant.g = paste(Chromosome, ":", Start_Position, "-", End_Position, " ", variant.g, sep = ""),
      variant.c = HGVSc,
      variant.p = HGVSp,
      annotation = Consequence,
      type = short,
      tier = tier,
      tier_description = tier_factor,
      variant_summary = VARIANT_SUMMARY,
      therapeutics = tx_string,
      resistances = rx_string,
      diagnoses = dx_string,
      prognoses = px_string,
      distance = NA_integer_,
      effect = MUTATION_EFFECT,
      effect_description = MUTATION_EFFECT_DESCRIPTION,
      major_count = major.count,
      minor_count = minor.count,
      major_snv_copies,
      minor_snv_copies,
      # total_copies, ## total_copies gets converted to estimated_altered_copies downstream
      altered_copies,
      segment_cn,
      ref,
      alt,
      VAF,
      vartype = "SNV",
      track = "variants",
      source = "oncokb_maf",
      is_multi_hit_per_gene
    )])
  }
  return(empty_output_oncokb)
}


#' @title oncotable
#' @description
#'
#' @param somatic_variant_annotations Path to annotated.bcf file
#' @param fusions Path to fusion.rds file
#' @param jabba_gg Path to jabba.simple.rds file
#' @param events Path to complex.rds file
#' @param signature_counts Path to signature_counts.txt file
#' @param karyograph Optional path to the karyograph.rds file
#' @param oncokb_snv Path to oncokb MAF file
#' @param oncokb_cna Path to oncokb CNA file
#' @param gencode_gr GRanges object with gencode annotations
#' @param amp.thresh SCNA amplification threshold to call an amp as a function of ploidy (4)
#' @param del.thresh SCNA deletion threshold for (het) del as a function of ploidy (by default cn = 1 will be called del, but this allows additoinal regions in high ploidy tumors to be considered het dels)
#' @param verbose logical flag
#' @export
oncotable <- function(
    pair,
    somatic_variant_annotations = NULL,
    fusions = NULL,
    jabba_gg = NULL,
    karyograph = NULL,
    events = NULL,
    signature_counts = NULL,
    multiplicity = NULL,
    oncokb_snv = NULL,
    oncokb_cna = NULL,
    oncokb_fusions = NULL,
    gencode,
    cytoband,
    verbose = TRUE,
    amp.thresh = 4,
    filter = "PASS",
    del.thresh = 0.5) {
  out <- data.table()

  if ("type" %in% names(mcols(gencode))) {
    pge <- gencode %Q% (type == "gene" & gene_type == "protein_coding")
  } else {
    pge <- gencode %Q% (gene_type == "protein_coding")
  }

  ## collect gene fusions
  # prefer fusions from oncokb
  if (!is.na(oncokb_fusions) && !is.null(oncokb_fusions) && file.exists(oncokb_fusions)) {
    out <- rbind(
      out,
      collect_oncokb_fusions(oncokb_fusions, pge, cytoband, verbose),
      fill = TRUE,
      use.names = TRUE
    )
  } else {
    out <- rbind(
      out,
      collect_gene_fusions(fusions, pge, verbose),
      fill = TRUE,
      use.names = TRUE
    )
  }

  ## collect complex events
  out <- rbind(
    out,
    collect_complex_events(events, verbose),
    fill = TRUE,
    use.names = TRUE
  )

  ## collect copy number
  # prefer oncokb cna
  if (!is.na(oncokb_cna) && !is.null(oncokb_cna) && file.exists(oncokb_cna)) {
    out <- rbind(
      out,
      collect_oncokb_cna(oncokb_cna, jabba_gg, pge, amp.thresh, del.thresh, karyograph, verbose),
      fill = TRUE,
      use.names = TRUE
    )
  } else {
    out <- rbind(
      out,
      collect_copy_number_jabba(jabba_gg, pge, amp.thresh, del.thresh, verbose, karyograph),
      fill = TRUE,
      use.names = TRUE
    )
  }

  ## collect gene mutations
  # prefer oncokb snv
  if (!is.na(oncokb_snv) && !is.null(oncokb_snv) && file.exists(oncokb_snv)) {
    out <- rbind(
      out,
      collect_oncokb(oncokb_snv, multiplicity, verbose),
      fill = TRUE,
      use.names = TRUE
    )
  } else {
    out <- rbind(
      out,
      collect_gene_mutations(somatic_variant_annotations, jabba_gg, filter, verbose),
      fill = TRUE,
      use.names = TRUE
    )
  }

  ## add gene locations
  gene_locations <- readRDS(system.file("extdata", "data", "gene_locations.rds", package = "Skilift"))

  # Merge gene locations and create location string
  if ("gene" %in% names(out)) {
    out[gene_locations,
      gene_location := paste0(i.seqnames, ":", i.start, "-", i.end),
      on = c(gene = "gene_name")
    ]
  }

  out$id <- pair

  # remove all rows for which data was not passed
  out <- out[!is.na(type)]

  # coerce all empty strings to NA without affecting columns of levels
  out <- out[, lapply(.SD, function(x) {
    if (is.factor(x)) x <- as.character(x) # Convert factors to characters
    ifelse(x == "", NA_character_, x) # Replace empty strings with NA
  })]

  if (verbose) message("done processing sample")
  return(out)
}

#' @name create_oncotable
#' @title create_oncotable
#' @description
#'
#' function to create oncotable for use with filtered_events_json
#'
#' @param cohort Cohort object containing sample information
#' @param amp_thresh_multiplier amp.thresh for oncotable is amp_thresh_multiplier*ploidy
#' @param gencode file to gencode annotations (uses v29lift37 by default)
#' @param outdir path to directory in which to write oncotable outputs
#' @param cores number of cores for parallel processing
#' @return None
#' @export
#' @author Shihab Dider, Joel Rosiene
create_oncotable <- function(
    cohort,
    amp_thresh_multiplier = 1.5,
    gencode = Skilift::get_default_gencode(),
    cytoband = system.file("extdata", "data", "cytoband.rds", package = "Skilift"),
    outdir,
    cores = 1) {
  if (!inherits(cohort, "Cohort")) {
    stop("Input must be a Cohort object")
  }

  if (nrow(cohort$inputs) == 0) {
    stop("No samples found in the cohort")
  }

  if (Sys.which("bcftools") == "") {
    stop("bcftools is not available on the system PATH. Try `module load htslib` first or install it.")
  } else {
    message("bcftools is available.")
  }

  if (gencode == "/gpfs/data/imielinskilab/DB/GENCODE/gencode.v19.annotation.gtf.nochr.rds") {
    message("using default gencode: /gpfs/data/imielinskilab/DB/GENCODE/gencode.v19.annotation.gtf.nochr.rds")
  }

  gencode <- process_gencode(gencode)
  cytoband <- process_cytoband(cytoband)

  if (amp_thresh_multiplier == 1.5) {
    message("using default amp_thres_multiplier: 1.5")
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  # Create a copy of the cohort to modify
  updated_cohort <- cohort

  # Add oncotable column if it doesn't exist
  if (!"oncotable" %in% names(updated_cohort$inputs)) {
    updated_cohort$inputs[, oncotable := NA_character_]
  }

  jabba_column = Skilift::DEFAULT_JABBA(object = cohort)
  
  results <- mclapply(seq_len(nrow(cohort$inputs)), function(i) {
    futile.logger::flog.threshold("ERROR")
    tryCatchLog(
      {
        row <- cohort$inputs[i]

        # Create output directory for this pair
        pair_outdir <- file.path(outdir, row$pair)
        if (!dir.exists(pair_outdir)) {
          dir.create(pair_outdir, recursive = TRUE)
        }

        # Get ploidy from jabba output, default to 2 if missing
        ploidy <- 2 # Default ploidy
        if (file.exists(row[[jabba_column]])) {
          ploidy_ggraph <- tryCatch(
            {
              process_jabba(row[[jabba_column]])
            },
            error = function(e) {
              msg <- sprintf("Error reading JaBbA file for %s: %s. Using default ploidy of 2.", row$pair, e$message)
              warning(msg)
            }
          )

          if (!is.null(ploidy_ggraph)) {
            ploidy <- ifelse(
              !is.null(ploidy_ggraph$meta$ploidy),
              ploidy_ggraph$meta$ploidy,
              ploidy_ggraph$ploidy
            )
          }
        }

        amp_thresh <- amp_thresh_multiplier * ploidy
        message(paste("Processing", row$pair, "using amp.thresh of", amp_thresh))


        # Run oncotable for this pair
        futile.logger::flog.threshold("ERROR")
        oncotable_result <- tryCatchLog(
          {
            oncotable(
              pair = row$pair,
              somatic_variant_annotations = row$somatic_variant_annotations,
              fusions = row$fusions,
              jabba_gg = row[[jabba_column]],
              karyograph = row$karyograph,
              events = row$events,
              signature_counts = row$signature_counts,
              multiplicity = row$multiplicity,
              oncokb_snv = row$oncokb_snv,
              oncokb_cna = row$oncokb_cna,
              oncokb_fusions = row$oncokb_fusions,
              gencode = gencode,
              cytoband = cytoband,
              verbose = TRUE,
              amp.thresh = amp_thresh,
              filter = "PASS",
              del.thresh = 0.5
            )
          },
          error = function(e) {
            print(sprintf("Error in oncotable for %s: %s", row$pair, e$message))
            NULL
          }
        )

        if (!is.null(oncotable_result)) {
          # Save successful results
          oncotable_path <- file.path(pair_outdir, "oncotable.rds")
          saveRDS(oncotable_result, oncotable_path)
          fwrite(oncotable_result, file.path(pair_outdir, "oncotable.txt"))
          return(list(index = i, path = oncotable_path))
        }
      },
      error = function(e) {
        print(sprintf("Unexpected error processing %s: %s", cohort$inputs[i]$pair, e$message))
        NULL
      }
    )
  }, mc.cores = cores, mc.preschedule = TRUE)


  # Update oncotable paths in the cohort
  results <- Filter(Negate(is.null), results)
  for (result in results) {
    updated_cohort$inputs[result$index, oncotable := result$path][]
  }
  message(sprintf("\nProcessing complete - results written to %s", outdir))

  return(updated_cohort)
}

  

#' @name create_filtered_events
#' @title create_filtered_events
#' @description
#' Create filtered events for a single sample
#'
#' @param pair patient id
#' @param oncotable oncotable task output
#' @param jabba_gg JaBbA output ggraph or complex
#' @param out_file path to write json
#' @param return_table TRUE/FALSE whether to return the data.table
#' @return data.table or NULL
#' @export
create_filtered_events <- function(
    pair,
    oncotable,
    jabba_gg,
    out_file,
    return_table = FALSE,
    cohort_type = "paired") {
  ot <- readRDS(oncotable)

  possible_drivers = empty_oncotable = structure(list(gene = character(0), gene_summary = character(0), 
    role = character(0), value = integer(0), vartype = character(0), 
    type = character(0), tier = integer(0), tier_description = character(0), 
    variant_summary = character(0), therapeutics = character(0), 
    resistances = logical(0), diagnoses = character(0), prognoses = character(0), 
    effect = character(0), effect_description = character(0), 
    fusion_genes = character(0), fusion_gene_coords = character(0), 
    track = character(0), source = character(0), variant.g = character(0), 
    variant.c = character(0), variant.p = character(0), annotation = character(0), 
    distance = logical(0), major_count = numeric(0), minor_count = numeric(0), 
    major_snv_copies = numeric(0), minor_snv_copies = numeric(0), 
    altered_copies = numeric(0), segment_cn = integer(0), ref = integer(0), 
    alt = integer(0), VAF = numeric(0), gene_location = character(0), 
    id = character(0)), row.names = c(NA, 0L), class = c("data.table", 
"data.frame"))

  if (NROW(ot) > 0) {
      # add a fusion_gene_coords column of NAs if no fusions
      if (!"fusion_gene_coords" %in% colnames(ot)) {
        ot[, fusion_genes := NA]
        ot[, fusion_gene_coords := NA]
      }
      # snvs <- ot[grepl("frameshift|missense|stop|disruptive", annotation, perl = TRUE)]
      snvs <- ot[ot$vartype == "SNV" & !is.na(ot$type)]
      if ("tier" %in% colnames(ot)) {
        snvs <- snvs[order(is.na(tier)), ]
      }
      if (NROW(snvs) > 0) {
        snvs[, is_unique_p := !is.na(variant.p) & !duplicated(cbind(gene, variant.p))]
        snvs[, is_unique_g := !duplicated(cbind(gene, variant.g))]
        # FIXME: variant.c should not be absent..
        ## remove_variant_c = FALSE
        ## if (is.null(snvs$variant.c)) {
        ##   snvs$variant.c = 1:NROW(snvs)
        ##   remove_variant_c = TRUE
        ## }
        snvs[, is_unique_c := !duplicated(cbind(gene, variant.c))]
        snvs <-  snvs[is_unique_p | (is_unique_g & is_unique_c)]
        snvs$is_unique_p = NULL
        snvs$is_unique_g = NULL
        snvs$is_unique_c = NULL
        ## if (remove_variant_c) snvs$variant.c = NA
        variant_concatenated = snvs[, paste(variant.p, "/", variant.c, sep = " ")]
        variant_concatenated = gsub(" / NA", "",
          gsub("NA / ", "",
            variant_concatenated,
            perl = TRUE
          ), perl = TRUE
        )
        snvs$variant.p = variant_concatenated
        ## snvs$variant.c = NULL
        rm("variant_concatenated")
        
      }

      homdels <- ot[ot$type == "homdel",][, vartype := "HOMDEL"][, type := "SCNA"]
      amps <- ot[ot$type == "amp",][, vartype := "AMP"][, type := "SCNA"]
      fusions <- ot[ot$type == "fusion",] ## fusion vartype is either fusion or outframe_fusion
      possible_drivers <- rbind(snvs, homdels, amps, fusions, fill = TRUE)
  }

  oncotable_col_to_filtered_events_col <- c(
    "id" = "id",
    "gene" = "gene",
    "fusion_genes" = "fusion_genes",
    "fusion_gene_coords" = "fusion_gene_coords",
    "value" = "fusion_cn",
    "vartype" = "vartype",
    "type" = "type",
    "variant.g" = "Variant_g",
    "variant.p" = "Variant",
    "altered_copies" = "estimated_altered_copies",
    "segment_cn" = "segment_cn",
    "ref" = "ref",
    "alt" = "alt",
    "VAF" = "VAF",
    "gene_location" = "Genome_Location",
    "tier" = "Tier",
    "role" = "role",
    "gene_summary" = "gene_summary",
    "variant_summary" = "variant_summary",
    "effect" = "effect",
    "effect_description" = "effect_description",
    "therapeutics" = "therapeutics",
    "resistances" = "resistances",
    "diagnoses" = "diagnoses",
    "prognoses" = "prognoses",
    "is_multi_hit_per_gene" = "is_multi_hit_per_gene"
  )
  filtered_events_columns <- names(possible_drivers)[names(possible_drivers) %in% names(oncotable_col_to_filtered_events_col)]

  res <- possible_drivers[, ..filtered_events_columns]
  intersected_columns <- intersect(filtered_events_columns, names(res))
  setnames(res, old = intersected_columns, new = oncotable_col_to_filtered_events_col[intersected_columns])

  res <- res %>% unique(., by = c("gene", "vartype", "Variant"))
  if (nrow(res) > 0) {
    res[, seqnames := tstrsplit(Genome_Location, ":", fixed = TRUE, keep = 1)]
    res[, start := tstrsplit(Genome_Location, "-", fixed = TRUE, keep = 1)]
    res[, start := tstrsplit(start, ":", fixed = TRUE, keep = 2)]
    res[, end := tstrsplit(Genome_Location, "-", fixed = TRUE, keep = 2)]
    res.mut <- res[vartype == "SNV"]
    if (nrow(res.mut) > 0) {
      # res.mut[, Variant := gsub("p.", "", Variant)]
      # res.mut[, vartype := "SNV"]
      # TODO:
      # truncating mutations are not always deletions
      # initial logic may be misleading calling all small mutations "SNV"
      # but we should encode this as something more robust
      # res.mut[type=="trunc", vartype := "DEL"]
      ## res.mut
      
      ## FIXME: Nothing seems to be necessary here at this point.
      NULL
    }
    res.fus = res[type == "fusion",] ## need to deal each class explicitly
    if (NROW(res.fus) > 0) {
      res.fus$gene = res.fus$fusion_genes
      is_inframe = res.fus$vartype == "fusion"
      is_outframe = res.fus$vartype == "outframe_fusion"
      ## If variant.p isn't found
      ## i.e. OncoKB isn't updated
      ## then only in frame or out of frame
      ## will be labeled
      variant.p.fus = base::get0(
        "Variant",
        as.environment(as.list(res.fus)),
        ifnotfound = rep_len("", NROW(res.fus))
      )
      fus_frame_label = ifelse(
        is_inframe,
        "In-Frame Fusion",
        ifelse(
          is_outframe,
          "Out-of-Frame Fusion",
          NA_character_
        )
      )
      if (any(is.na(fus_frame_label))) stop("A fusion was not labeled as in-frame or out-of-frame")
      variant_label = paste(
          fus_frame_label,
          variant.p.fus
      )
      variant_label = trimws(variant_label)
      variant_label = gsub("[[:space:]]{2,}", "", variant_label)
      res.fus$Variant = variant_label

      res.fus$estimated_altered_copies = res.fus$fusion_cn

      ## Filtering Fusions based on COSMIC and Duncavage:
      ## If the exact match is found -
      ## remove any reverse reciprocal fusion
      fusions_curated = Skilift:::get_curated_fusions()
      fusions_recurrent = fusions_curated$fusions_recurrent

      blacklist = fusions_curated$blacklist
      fusions_recurrent = fusions_recurrent[!fusions_recurrent %in% blacklist]

      fus_lst = data.table::tstrsplit(fusions_recurrent, "::")
      forwards = paste(fus_lst[[1]], "::", fus_lst[[2]], sep = "")
      reverse = paste(fus_lst[[2]], "::", fus_lst[[1]], sep = "")

      fgenes = gsub("@.*", "", res.fus$fusion_genes)
      fgenes = gsub("\\([0-9]+\\)", "", fgenes, perl = TRUE)

      is_in_forward = fgenes %in% forwards
      is_in_reverse = fgenes %in% reverse
      ## Some are present in both directions,
      ## but are not well-curated. So leave them in.
      is_both_f_and_r = fgenes %in% forwards[forwards %in% reverse]
      is_valid = (
          is_in_forward & !is_in_reverse
      ) | is_both_f_and_r
      ## Other is anything that's completely outside of the list
      ## in either direction.
      is_other = ! fgenes %in% c(forwards, reverse)

      res.fus = base::subset(res.fus, subset = is_valid | is_other)
      
    }
    res.cn.dt = res.cn <- res[(
      type == "SCNA" ## FIXME: redundant logic for now
      | vartype %in% c("AMP", "HOMDEL")
    )]
    if (nrow(res.cn) > 0) {
      # jab <- readRDS(jabba_gg)
      # res.cn.gr <- GRanges(res.cn)
      # res.cn.dt
      # res.cn.gr <- gr.val(res.cn.gr, jab$nodes$gr, c("cn", "cn.low", "cn.high"))
      # res.cn.dt <- as.data.table(res.cn.gr)
      # FIXME: value is actually min_cn is renamed fusion_cn!!
      # Also, we already picked out the corresponding segment cn to use (min_cn)
      res.cn.dt$min_cn = res.cn.dt$fusion_cn 
      res.cn.dt[, estimated_altered_copies := abs(min_cn - 2)]
      res.cn.dt[, segment_cn := min_cn]
      res.cn.dt[, Variant := vartype]
      res.cn.dt$fusion_cn = NA_real_

      # To fix very small CNAs that can appear due to forcing in reciprocal
      # junctions with very small gaps.
      ## res.cn.dt = res.cn.dt[!res.cn.dt$node_width_jab < 1000,]
      # remove redundant columns since already added to Variant
      
      res.cn.dt[, c("min_cn", "cn", "cn.high", "cn.low", "width", "strand") := NULL]
    }
    res.final <- rbind(res.mut, res.cn.dt, res.fus, fill = TRUE)
    res.final[, sample := pair]
    if (identical(cohort_type, "heme")) {
      res.final <- select_heme_events(res.final)
    }
    res.final$type = tools::toTitleCase(res.final$type)
    ### FIXME: REMOVING Y CHROMOSOME UNTIL WE UPDATE DRYCLEAN
    res.final <- res.final[res.final$seqnames != "Y"]
    ### FIXME ^^^: REMOVING Y CHROMOSOME UNTIL WE UPDATE DRYCLEAN
    write_json(res.final, out_file, pretty = TRUE)
    
    ## FIXME: Decide whether to either propagate return_table to lift_mvp in lift_heme somehow, or change this
    ## conditional logic to be based on cohort_type. Can also be left alone.
    ## Note that return_table for debugging purposes as well.
    if (return_table) {
      return(res.final)
    }
  }
  NULL
}


#' @name lift_filtered_events
#' @title lift_filtered_events
#' @description
#' Create filtered events for all samples in a cohort
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @return None
#' @export
lift_filtered_events <- function(cohort, output_data_dir, cores = 1, return_table = TRUE) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }
    
    # Validate required columns exist
    # required_cols <- c("pair", "oncotable", "jabba_gg")
	jabba_column = Skilift::DEFAULT_JABBA(object = cohort)
	required_cols <- c("pair", "oncotable", jabba_column)
    missing_cols <- required_cols[!required_cols %in% names(cohort$inputs)]
    if (length(missing_cols) > 0) {
        print("Missing required columns in cohort: ", paste(missing_cols, collapse = ", "))
    }
    
    cohort_type = cohort$type
	
	# Create a copy of the cohort to modify

	# Add oncotable column if it doesn't exist
	if (!"string_summary" %in% names(cohort$inputs)) {
		cohort$inputs[, string_summary := ""]
	}

    # Process each sample in parallel
    results = mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        row <- cohort$inputs[i,]
        pair_dir <- file.path(output_data_dir, row$pair)
        
        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }
        
        out_file <- file.path(pair_dir, "filtered.events.json")
        highlights_out_file <- file.path(pair_dir, "highlights.json")

        out = NULL
		string_summary = ""
        futile.logger::flog.threshold("ERROR")
        tryCatchLog({
            out <- create_filtered_events(
                pair = row$pair,
                oncotable = row$oncotable,
                jabba_gg = row[[jabba_column]],
                out_file = out_file,
                return_table = return_table,
                cohort_type = cohort_type
            )
            if (identical(cohort_type, "heme")) {
              create_heme_highlights(
                events_tbl = out,
                jabba_gg = row[[jabba_column]],
                out_file = highlights_out_file,
                tumor_type = row$tumor_type,
                cohorttuple = row
              )
            }
			string_summary = create_summary(
				events_tbl = out,
				cohort_type = cohort_type
			)
        }, error = function(e) {
            print(sprintf("Error processing %s: %s", row$pair, e$message))
            NULL
        })

        # return(out)
		return(
			list(
				pair = row$pair,
				string_summary = string_summary
			)
		)
    }, mc.cores = cores, mc.preschedule = TRUE)

	results = results[!sapply(results, is.null)]

	if (NROW(results) > 0) {
		# Transpose the list
		results = do.call(Map, c(f = c, results))
		cohort$inputs[results$pair, string_summary := results$string_summary]
	}

	return(cohort)
}

#' Select Heme events from Addy's hemedb
#' 
#' Coarse selection of Heme events for first pass filtering
#'
#' @export
select_heme_events <- function(
  filtered_events, 
  hemedb_path = Skilift:::HEMEDB()
) {
  hemedb = readRDS(hemedb_path) 

  # wtf = merge(events_tbl, hemedb_genes, by.x = "gene", by.y = "GENE")[5]
  is_small_mutation_heme <- (
    filtered_events$vartype == "SNV" &
      (
        filtered_events$gene %in% hemedb$GENE |
          filtered_events$Tier %in% c(1, 2)
      )
  )
  
  is_other_event = filtered_events$vartype != "SNV"
  # if (any(filtered_events$gene %in% c("FLT3", "DUX4", "KMT2A"))) {
  #   .NotYetImplemented()
  # }
  return(
    filtered_events[is_small_mutation_heme | is_other_event]
  )
}


#' Merge oncokb and multiplicity
#'
#' Finds overlaps between oncokb maf coordinates and multiplicity
#' Oncokb is in MAF coordinates which is based on the altered bases,
#' while multiplicity encodes VCF coordinates based on reference bases.
#'
#' @export
merge_oncokb_multiplicity <- function(
  oncokb, 
  multiplicity, 
  overwrite = FALSE,
  cols.keep = c(
    "annotation", "ref", "alt", "ref_denoised", "alt_denoised", "normal.ref", "normal.alt", "variant.c", "variant.p",
    "variant.g", "major.count", "minor.count", "major_snv_copies", "minor_snv_copies",
    "total_snv_copies", "total_copies", "VAF", "cn", "altered_copies"
  ),
  other.cols.keep = NULL
) {
  is_character_oncokb = is.character(oncokb)
  is_length_one = NROW(oncokb) == 1
  is_exists_oncokb = is_character_oncokb && is_length_one && file.exists(oncokb)
  is_oncokb_rds = is_exists_oncokb && grep("rds$", oncokb)
  is_oncokb_txt = is_exists_oncokb && grep("(c|t)sv$|maf$|txt$", oncokb)

  if (is_oncokb_txt) {
    oncokb <- data.table::fread(oncokb)
  } else if (is_oncokb_rds) {
    oncokb <- readRDS(oncokb)
  }

  is_character_mult = is.character(multiplicity)
  is_length_one = NROW(multiplicity) == 1
  is_exists_mult = is_character_mult && is_length_one && file.exists(multiplicity)
  is_mult_rds = is_exists_mult && grep("rds$", multiplicity)
  is_mult_txt = is_exists_mult && grep("(c|t)sv$|maf$|txt$", multiplicity)

  if (is_mult_rds) {
    multiplicity <- readRDS(multiplicity)
  } else if (is_mult_txt) {
    multiplicity = fread(multiplicity)
  }
  if (inherits(multiplicity, "GRanges")) {
    gr_multiplicity <- multiplicity
  } else {
    multiplicity_tmp = tryCatch(
      {
        dt2gr(multiplicity)
      }, error = function(e) tryCatch(
      {
        as(multiplicity, "GRanges")
      }, error = function(e) NULL)
    )
    if (is.null(multiplicity_tmp)) {
      stop("multiplicity must be coercible to GRanges")
    }
    gr_multiplicity = multiplicity_tmp
  }

  if (inherits(oncokb, "GRanges")) {
    gr_oncokb <- oncokb
    oncokb = gr2dt(oncokb)
  } else {
    oncokb_tmp = tryCatch(
      {
        dt2gr(oncokb)
      }, error = function(e) tryCatch(
      {
        as(oncokb, "GRanges")
      }, error = function(e) NULL)
    )
    if (is.null(oncokb_tmp)) {
      stop("oncokb must be coercible to GRanges")
    }
    gr_oncokb = oncokb_tmp
  }

  gr_oncokb$ALT <- gr_oncokb$Allele
  gr_oncokb$gene <- gr_oncokb$Hugo_Symbol

  mc <- S4Vectors::mcols(gr_multiplicity)
  checknormalcols <- c("normal.ref", "normal.alt")
  for (col in checknormalcols) {
    if (!col %in% names(mc)) {
      mc[[col]] <- NA_integer_
    }
  }
  S4Vectors::mcols(gr_multiplicity) <- mc
  is_oncokb_empty = NROW(gr_oncokb) == 0
  is_multiplicity_empty = NROW(gr_multiplicity) == 0

  # oncokb should always be a data.table at this point, regardless whether empty or not 
  # (columns are always there)
  columns_to_assign_na_in_oncokb = cols.keep[!cols.keep %in% names(oncokb)] 
  if (is_oncokb_empty || is_multiplicity_empty)  {
	for (col in columns_to_assign_na_in_oncokb) {
		oncokb[[col]] = rep_len(NA, NROW(oncokb))
	}
	return(oncokb) ## Should be returned as original empty data.table + multiplicity columns (na'd for downstream robustness)
  }

  ov <- gUtils::gr.findoverlaps(gr_oncokb, gr_multiplicity, by = c("gene", "ALT"), type = "equal")
  ovQuery <- data.table(query.id = integer(0), subject.id = integer(0))
  if (NROW(ov) > 0) {
    ovQuery <- gUtils::gr2dt(ov)[, .(query.id, subject.id)]
  }
  missingIds <- setdiff(1:NROW(gr_oncokb), ov$query.id)

  missingOvQuery <- data.table(query.id = integer(0), subject.id = integer(0))

  if (length(missingIds) > 0) {
    dt_oncokb <- gUtils::gr2dt(gr_oncokb[missingIds])
    invisible(dt_oncokb[, Reference_Allele_Fixed := ifelse(Reference_Allele == "-", "", Reference_Allele)])
    invisible({
      dt_oncokb[Variant_Type == "INS", end := start + nchar(Reference_Allele_Fixed)]
      dt_oncokb[Variant_Type == "DEL", start := end - nchar(Reference_Allele_Fixed)]
    })
    dt_oncokb$oid <- missingIds
    ovMissing <- gUtils::gr.findoverlaps(
      gUtils::dt2gr(dt_oncokb), gr_multiplicity,
      by = c("gene", "ALT"), type = "equal", qcol = c("oid")
    )
    if (NROW(ovMissing) > 0) {
      missingOvQuery <- gr2dt(ovMissing)[, .(query.id = oid, subject.id)]
    }
  }
  ovQuery <- rbind(ovQuery, missingOvQuery)

  # Logic for when there's exact match in coordinates but REF/ALT slight mismatch due to VCF -> maf parsing.
  missingIds = setdiff(1:NROW(gr_oncokb), ovQuery$query.id)
  missingOvQuery = data.table(query.id = integer(0), subject.id = integer(0))
  
  if (length(missingIds) > 0) {
    gr_oncokb_missing = gr_oncokb[missingIds]
    gr_oncokb_missing$oid = missingIds
    ovMissing = gUtils::gr.findoverlaps(
      gr_oncokb_missing, gr_multiplicity
     ,
     by = "gene",
     type = "equal",
      qcol = c("oid"),
     )
    if (NROW(ovMissing) > 0) {
      missingOvQuery = (
        gr2dt(ovMissing)
        [, .(query.id = oid, subject.id)]
        [!duplicated(query.id)]
      )
    }
  }
  ovQuery = rbind(ovQuery, missingOvQuery)

  # Logic for when there's inexact coordinate match due to VCF -> maf parsing.
  missingIds = setdiff(1:NROW(gr_oncokb), ovQuery$query.id)
  missingOvQuery = data.table(query.id = integer(0), subject.id = integer(0))

  if (length(missingIds) > 0) {
    gr_oncokb_missing = gr_oncokb[missingIds]
    gr_oncokb_missing$oid = missingIds
    ovMissing = gUtils::gr.findoverlaps(
      gr_oncokb_missing, gr_multiplicity
     ,
      by = "gene",
      qcol = c("oid")
    )
    if (NROW(ovMissing) > 0) {
      missingOvQuery = (
        gr2dt(ovMissing)
        [, .(query.id = oid, subject.id)]
        [!duplicated(query.id)]
      )
    }
  }
  ovQuery = rbind(ovQuery, missingOvQuery)

  # Logic for when there's inexact coordinate match due to VCF -> maf parsing.
  missingIds = setdiff(1:NROW(gr_oncokb), ovQuery$query.id)
  missingOvQuery = data.table(query.id = integer(0), subject.id = integer(0))

  if (length(missingIds) > 0) {
    gr_oncokb_missing = gr_oncokb[missingIds]
    gr_oncokb_missing$oid = missingIds
    ovMissing = gUtils::gr.findoverlaps(
      gr_oncokb_missing, gr_multiplicity
     ,
      ## by = "gene",
      qcol = c("oid")
    )
    if (NROW(ovMissing) > 0) {
      missingOvQuery = (
        gr2dt(ovMissing)
        [, .(query.id = oid, subject.id)]
        [!duplicated(query.id)]
      )
    }
  }
  ovQuery = rbind(ovQuery, missingOvQuery)

  # cols.keep <- c(
  #   "ref", "alt", "ref_denoised", "alt_denoised", "normal.ref", "normal.alt",
  #   "variant.g", "major.count", "minor.count", "major_snv_copies", "minor_snv_copies",
  #   "total_snv_copies", "total_copies", "VAF", "cn", "altered_copies"
  # )

  is_othercols_null = is.null(other.cols.keep)
  is_othercols_character = (
    is.character(other.cols.keep) 
    && NROW(other.cols.keep) > 0
  )
  is_othercols_valid = (
    is_othercols_character
    && all(nzchar(other.cols.keep))
    && !all(other.cols.keep %in% cols.keep)
  )
  if (!is_othercols_null && is_othercols_character && is_othercols_valid) {
    cols.keep = c(cols.keep, other.cols.keep[!other.cols.keep %in% cols.keep])
  }
  cols.keep <- cols.keep[cols.keep %in% names(mcols(gr_multiplicity))]
  subject <- base::subset(gUtils::gr2dt(gr_multiplicity), select = cols.keep)
  skey <- data.table::setkey(ovQuery, query.id)
  skey <- skey[list(1:NROW(oncokb))]

  if (!identical(overwrite, TRUE)) {
    oncokb_multiplicity <- cbind(oncokb, subject[skey$subject.id])
  } else {
    oncokb_multiplicity <- oncokb
    # mc_oncokb_multiplicity = S4Vectors::mcols(oncokb_multiplicity)
    subject_ord <- subject[skey$subject.id]
    for (multiplicity_col in cols.keep) {
      oncokb_multiplicity[[multiplicity_col]] <- NULL
      oncokb_multiplicity[[multiplicity_col]] <- subject_ord[[multiplicity_col]]
    }
    # S4Vectors::mcols(oncokb_multiplicity) = mc_oncokb_multiplicity
  }

  oncokb_multiplicity$multiplicity_id_match <- skey$subject.id
  return(oncokb_multiplicity)
}

