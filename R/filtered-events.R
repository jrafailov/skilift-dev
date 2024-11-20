#' @title process_gencode
#' @description
#'
#' Helper script to process gencode parameter
#'
#' @param gencode path to gencode file. Gencode file must be either rds or some format accepted by rtracklayer::import (e.g. GTF) with just a single entry for each gene (so gencode entries for each gene are collapse to a single range). The input could be .gtf or .rds with GRanges object, or a GRanges object i.e. resulting from importing the (appropriate) GENCODE .gtf via rtracklayer, note: this input is only used in CNA to gene mapping.
#' @return gencode_gr GRanges
#' @author Marcin Imielinski
process_gencode = function(gencode = NULL){
  if (is.null(gencode))
    stop('gencode file must be provided')
  else if (is.character(gencode))
  {
    if (grepl('.rds$', gencode))
      gencode = readRDS(gencode)
    else
      gencode = rtracklayer::import(gencode)
  }
  return(gencode)
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
    if (verbose) message('Fusions file is missing or does not exist.')
    return(data.table(type = NA, source = 'fusions'))
  }
  
  if (verbose) message('pulling fusions')
  fus <- readRDS(fusions)$meta
  
  if (nrow(fus) == 0) {
    if (verbose) message('No fusions found in the file.')
    return(data.table(type = NA, source = 'fusions'))
  }
  
  non_silent_fusions <- fus[silent == FALSE, ]
  unique_fusions <- non_silent_fusions[!duplicated(genes), ]
  unique_fusions[, vartype := ifelse(in.frame == TRUE, 'fusion', 'outframe_fusion')]
  
  split_genes <- function(genes) unlist(strsplit(genes, ','))
  gene_lengths <- function(genes) sapply(strsplit(genes, ','), length)
  
  fus <- unique_fusions[, .(
    gene = split_genes(genes),
    vartype = rep(vartype, gene_lengths(genes)),
    fusion_genes = rep(genes, gene_lengths(genes))
  )][, `:=`(track = 'variants', type = vartype, source = 'fusions')]
  
  get_gene_coords <- function(genes) {
    coords <- lapply(genes, function(gene) {
      gene_ranges <- pge[mcols(pge)$gene_name == gene]
      paste0(seqnames(gene_ranges), ":", start(gene_ranges), "-", end(gene_ranges))
    })
    paste(unlist(coords), collapse = ",")
  }
  
  fus[, fusion_gene_coords := unlist(lapply(strsplit(fusion_genes, ','), get_gene_coords))]
  
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
    if (verbose) message('Complex events file is missing or does not exist.')
    return(data.table(type = NA, source = 'complex'))
  }
  
  if (verbose) message('pulling complex events')
  sv <- readRDS(complex)$meta$events
  
  if (nrow(sv) == 0) {
    if (verbose) message('No complex events found in the file.')
    return(data.table(type = NA, source = 'complex'))
  }
  
  sv_summary <- sv[, .(value = .N), by = type]
  simple_sv_types <- c('del', 'dup', 'invdup', 'tra', 'inv')
  sv_summary[, track := ifelse(type %in% simple_sv_types, 'simple sv', 'complex sv')]
  sv_summary[, source := 'complex']
  
  return(sv_summary)
}

#' Get gene amplifications and deletions from jabba object
#'
#' Copied from github::mskilab-org/skitools.
#' This function takes in a jabba object and returns
#' amplifications and deletions
#'
#' @param jab character path to jabba file or gGraph object
get_gene_ampdels_from_jabba <- function(jab, pge, amp.thresh = 4, del.thresh = 0.5, nseg = NULL) {
    gg <- jab
    if (!inherits(gg, "gGraph")) {
        gg <- gG(jabba = jab)
    }
    gene_CN <- skitools::get_gene_copy_numbers(gg, gene_ranges = pge, nseg = nseg)
    gene_CN[, `:=`(type, NA_character_)]
    gene_CN[min_normalized_cn >= amp.thresh, `:=`(type, "amp")]
    gene_CN[min_cn > 1 & min_normalized_cn < del.thresh, `:=`(
        type,
        "del"
    )]
    gene_CN[min_cn == 1 & min_cn < ncn, `:=`(type, "hetdel")]
    gene_CN[min_cn == 0, `:=`(type, "homdel")]
    return(gene_CN[!is.na(type)])
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
  karyograph = NULL
) {
  if (is.null(jabba_rds) || !file.exists(jabba_rds)) {
    if (verbose) message('Jabba RDS file is missing or does not exist.')
    return(data.table(type = NA, source = 'jabba_rds'))
  }
  
  if (verbose) message('pulling jabba_rds to get SCNA and purity / ploidy')
  jab <- readRDS(jabba_rds)
  jabpurity <- ifelse(!is.null(jab$meta$purity), jab$meta$purity, jab$purity)  
  jabploidy <- ifelse(!is.null(jab$meta$ploidy), jab$meta$ploidy, jab$ploidy)
  result <- data.table(
    value = c(jabpurity, jabploidy),
    type = c('purity', 'ploidy'),
    track = 'pp'
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
    nseg = nseg
  )
  
  if (nrow(scna)) {
    scna[, track := 'variants'][, source := 'jabba_rds'][, vartype := 'scna']
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
  #TODO: update this function to use sigprofiler
  if (!is.null(signature_counts) && file.exists(signature_counts)) {
    if (verbose) message('pulling signature_counts')
    sig <- fread(signature_counts)
    sig <- sig[, .(value = num_events, type = Signature, etiology = Etiology, frac = frac.events, track = 'signature', source = 'signature_counts')]
    return(sig)
  } else {
    return(data.table(type = NA, source = 'signature_counts'))
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
  verbose = TRUE
) {
  if (is.null(annotated_bcf) || !file.exists(annotated_bcf)) {
    if (verbose) message('Annotated BCF file is missing or does not exist.')
    return(data.table(type = NA, source = 'annotated_bcf'))
  }

  if (verbose) message('pulling annotated_bcf using FILTER=', filter)
  
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
  if (verbose) message(length(bcf), ' variants pass filter')
  
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
    type = c('count', 'density'),
    track = 'tmb',
    source = 'annotated_bcf'
  )
  
  keepeff <- c('trunc', 'cnadel', 'cnadup', 'complexsv', 'splice', 'inframe_indel', 'fusion', 'missense', 'promoter', 'regulatory', 'mir')
  bcf <- bcf[bcf$short %in% keepeff]
  if (verbose) message(length(bcf), ' variants pass keepeff')
  
  if (length(bcf) == 0) {
    return(mut.density)
  }
  
  bcf$variant.g <- paste0(seqnames(bcf), ':', start(bcf), '-', end(bcf), ' ', bcf$REF, '>', bcf$ALT)
  vars <- gr2dt(bcf)[, .(gene, vartype, variant.g, variant.p, distance, annotation, type = short, track = 'variants', source = 'annotated_bcf')]
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
collect_oncokb_cna <- function(oncokb_cna, verbose = TRUE) {
  if (is.null(oncokb_cna) || !file.exists(oncokb_cna)) {
    if (verbose) message('OncoKB CNA file is missing or does not exist.')
    return(data.table(type = NA, source = 'oncokb_cna'))
  }

  oncokb_cna <- data.table::fread(oncokb_cna)
  
  if (NROW(oncokb_cna) > 0) {
    oncokb_cna = parse_oncokb_tier(
        oncokb_cna, 
        tx_cols = c("LEVEL_1", "LEVEL_2"), 
        rx_cols = c("LEVEL_R1"),
        dx_cols = c("LEVEL_Dx1"),
        px_cols = c("LEVEL_Px1")
    )
    return(oncokb_cna[, .(
            gene = HUGO_SYMBOL,
            role = Role,
            value = min_cn,
            type = ifelse(
              ALTERATION == "Amplification",
              "amp",
              ifelse(ALTERATION == "Deletion", "homdel", NA_character_)
            ),
            tier = tier,
            tier_description = tier_factor,
            therapeutics = tx_string,
            resistances = rx_string,
            diagnoses = dx_string,
            prognoses = px_string,
            effect = MUTATION_EFFECT,
            track = "scna",
            source = "oncokb_cna"
    )])
  }
  return(data.table(type = NA, source = 'oncokb_cna'))
}


#' @title collect_oncokb_fusions
#' @description
#' Collects OncoKB Fusion data from a specified file and processes it.
#'
#' @param oncokb_fusions Path to the oncokb fusions file.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed OncoKB Fusion information.
#' @author Kevin Hadi
collect_oncokb_fusions <- function(oncokb_fusions, pge, verbose = TRUE) {
  out = data.table(vartype = NA, source = 'oncokb_fusions')
  if (is.null(oncokb_fusions) || !file.exists(oncokb_fusions)) {
    if (verbose) message('OncoKB Fusions file is missing or does not exist.')
    return()
  }

  oncokb_fusions <- data.table::fread(oncokb_fusions)
  
  if (NROW(oncokb_fusions) > 0) {
    oncokb_fusions = parse_oncokb_tier(
        oncokb_fusions, 
        tx_cols = c("LEVEL_1", "LEVEL_2"), 
        rx_cols = c("LEVEL_R1"),
        dx_cols = c("LEVEL_Dx1"),
        px_cols = c("LEVEL_Px1")
    )
    
    non_silent_fusions <- oncokb_fusions[silent == FALSE, ] # already de-duped
    non_silent_fusions[, vartype := ifelse(in.frame == TRUE, 'fusion', 'outframe_fusion')]

    if (!NROW(non_silent_fusions) > 0) return(out)

    genes = strsplit(non_silent_fusions$FUSIONS, "-")
    genes_matrix = do.call(rbind, genes)
    ixA = match(genes_matrix[,1], pge$gene_name)
    ixB = match(genes_matrix[,2], pge$gene_name)
    coordA = gUtils::gr.string(pge[ixA])
    coordB = gUtils::gr.string(pge[ixB])
    non_silent_fusions$fusion_gene_coords = paste(coordA, coordB, sep = ",")
    out = non_silent_fusions[, .(
            gene = HUGO_SYMBOL,
            role = Role,
            value = min_cn,
            vartype,
            tier = tier,
            tier_description = tier_factor,
            therapeutics = tx_string,
            resistances = rx_string,
            diagnoses = dx_string,
            prognoses = px_string,
            effect = MUTATION_EFFECT,
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
parse_oncokb_tier = function(
    oncokb, 
    tx_cols = c("LEVEL_1", "LEVEL_2"), 
    rx_cols = c("LEVEL_R1"),
    dx_cols = c("LEVEL_Dx1"),
    px_cols = c("LEVEL_Px1")
) {
    .concat_string = function(oncokb, cols) {
        out_string = lapply(base::subset(oncokb, select = cols), function(y) (strsplit(y, ",")))
        concat_out = IRanges::CharacterList(Reduce(concat_vectors, out_string))
        concat_out = S4Vectors::unique(concat_out)
        concat_out = concat_out[!is.na(concat_out)]
        concat_out[S4Vectors::elementNROWS(concat_out) == 0] = NA_character_
        concat_out = stringi::stri_c_list(as.list(concat_out), sep = ",")
        return(concat_out)
    }
    is_actionable = logical(NROW(oncokb))
    for (col in c(tx_cols, rx_cols, dx_cols, px_cols)) {
        oncokb[[col]] = as.character(oncokb[[col]])
        is_actionable = is_actionable | (!is.na(oncokb[[col]]) & base::nzchar(oncokb[[col]]))
    }
    oncokb$is_actionable = is_actionable
    oncokb$is_oncogenic = oncokb$ONCOGENIC %in% c("Likely Oncogenic", "Oncogenic")
    tier_factor = ifelse(
        oncokb$is_actionable, "Clinically Actionable",
        ifelse(oncokb$is_oncogenic, "Clinically Significant", "VUS")
    ) %>% factor(c("Clinically Actionable", "Clinically Significant", "VUS"))

    oncokb$tier_factor = tier_factor
    oncokb$tier = as.integer(tier_factor)
    
    oncokb$tx_string = .concat_string(oncokb, tx_cols)
    oncokb$rx_string = .concat_string(oncokb, rx_cols)
    oncokb$dx_string = .concat_string(oncokb, dx_cols)
    oncokb$px_string = .concat_string(oncokb, px_cols)
    
    return(oncokb)
}

#' @title collect_oncokb
#' @description
#' Collects OncoKB mutation data from a specified file and processes it.
#'
#' @param oncokb_maf Path to the oncokb MAF file.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed OncoKB mutation information.
collect_oncokb <- function(oncokb_maf, verbose = TRUE) {
  if (is.null(oncokb_maf) || !file.exists(oncokb_maf)) {
    if (verbose) message('OncoKB MAF file is missing or does not exist.')
    return(data.table(type = NA, source = 'oncokb_maf'))
  }

  snpeff_ontology = readRDS(system.file("extdata", "data", "snpeff_ontology.rds", package = "Skilift"))
  oncokb <- data.table::fread(oncokb_maf)
  
  if (NROW(oncokb) > 0) {
    oncokb$snpeff_ontology <- snpeff_ontology$short[match(oncokb$Consequence, snpeff_ontology$eff)]
    oncokb = parse_oncokb_tier(
        oncokb, 
        tx_cols = c("LEVEL_1", "LEVEL_2"), 
        rx_cols = c("LEVEL_R1"),
        dx_cols = c("LEVEL_Dx1"),
        px_cols = c("LEVEL_Px1")
    )
    return(oncokb[, .(
            gene = Hugo_Symbol, 
            role = Role,
            variant.g = paste("g.",  Start_Position, "-", End_Position, sep = ""), 
            variant.c = HGVSc,
            variant.p = HGVSp,
            annotation = Consequence,
            type = snpeff_ontology,
            tier = tier,
            tier_description = tier_factor,
            therapeutics = tx_string,
            resistances = rx_string,
            diagnoses = dx_string,
            prognoses = px_string,
            distance = NA_integer_,
            effect = MUTATION_EFFECT,
            major_count = major.count, 
            minor_count = minor.count, 
            major_snv_copies, 
            minor_snv_copies,
            total_copies, 
            segment_cn,
            ref,
            alt,
            VAF,
            vartype = "SNV",
            track = "variants",
            source = "oncokb_maf"
    )])
  }
  return(data.table(type = NA, source = 'oncokb_maf'))
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
oncotable = function(
  pair,
  somatic_variant_annotations = NULL,
  fusions = NULL,
  jabba_gg = NULL,
  karyograph = NULL,
  events = NULL,
  signature_counts = NULL,
  oncokb_snv = NULL,
  oncokb_cna = NULL, 
  gencode,
  verbose = TRUE,
  amp.thresh = 4,
  filter = 'PASS',
  del.thresh = 0.5
) {
  out <- data.table()

  if ('type' %in% names(mcols(gencode))) {
    pge <- gencode %Q% (type == 'gene' & gene_type == 'protein_coding')
  } else {
    pge <- gencode %Q% (gene_type == 'protein_coding')
  }

  ## collect gene fusions
  out <- rbind(out, collect_gene_fusions(
    fusions,
    pge,
    verbose
  ), fill = TRUE, use.names = TRUE)

  ## collect complex events
  out <- rbind(
    out,
    collect_complex_events(events, verbose),
    fill = TRUE,
    use.names = TRUE
  )

  ## collect copy number / jabba
  out <- rbind(
    out,
    collect_copy_number_jabba(jabba_gg, pge, amp.thresh, del.thresh, verbose, karyograph),
    fill = TRUE,
    use.names = TRUE
  )

  ## collect oncokb_cna
  out <- rbind(
    out,
    collect_oncokb_cna(oncokb_cna, verbose),
    fill = TRUE,
    use.names = TRUE
  )

  ## collect gene mutations
  out <- rbind(
    out,
    collect_gene_mutations(somatic_variant_annotations, jabba_gg, filter, verbose),
    fill = TRUE,
    use.names = TRUE
  )

  ## collect oncokb
  out <- rbind(
    out,
    collect_oncokb(oncokb_snv, verbose),
    fill = TRUE,
    use.names = TRUE
  )

  ## add gene locations
  gene_locations = readRDS(system.file("extdata", "data", "gene_locations.rds", package = "Skilift"))

  # Merge gene locations and create location string
  if ("gene" %in% names(out)) {
    out[gene_locations, 
        gene_location := paste0(i.seqnames, ":", i.start, "-", i.end),
        on = c(gene = "gene_name")]
  }

  out$id = pair

  # remove all rows for which data was not passed 
  out <- out[!is.na(type)]

  if (verbose) message('done processing sample')
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
    gencode = "~/DB/GENCODE/gencode.v29lift37.annotation.nochr.rds",
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

    if (gencode == "~/DB/GENCODE/gencode.v29lift37.annotation.nochr.rds") {
        message("using default gencode: ~/DB/GENCODE/gencode.v29lift37.annotation.nochr.rds")
    }

    gencode <- process_gencode(gencode)

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
    results <- mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        tryCatch({
            row <- cohort$inputs[i]
            
            # Create output directory for this pair
            pair_outdir <- file.path(outdir, row$pair)
            if (!dir.exists(pair_outdir)) {
                dir.create(pair_outdir, recursive = TRUE)
            }
            
            # Get ploidy from jabba output, default to 2 if missing
            ploidy <- 2  # Default ploidy
            if (file.exists(row$jabba_gg)) {
                ploidy_ggraph <- tryCatch({
                    readRDS(row$jabba_gg)
                }, error = function(e) {
                    msg <- sprintf("Error reading JaBbA file for %s: %s. Using default ploidy of 2.", row$pair, e$message)
                    warning(msg)
                })
                
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
            oncotable_result <- tryCatch({
                oncotable(
                    pair = row$pair,
                    somatic_variant_annotations = row$somatic_variant_annotations,
                    fusions = row$fusions,
                    jabba_gg = row$jabba_gg,
                    karyograph = row$karyograph,
                    events = row$events,
                    signature_counts = row$signature_counts,
                    oncokb_snv = row$oncokb_snv,
                    oncokb_cna = row$oncokb_cna,
                    gencode = gencode,
                    verbose = TRUE,
                    amp.thresh = amp_thresh,
                    filter = "PASS",
                    del.thresh = 0.5
                )
            }, error = function(e) {
                msg <- sprintf("Error in oncotable for %s: %s", row$pair, e$message)
                warning(msg)
                return(NULL)
            })

            if (!is.null(oncotable_result)) {
                # Save successful results
                oncotable_path <- file.path(pair_outdir, "oncotable.rds")
                saveRDS(oncotable_result, oncotable_path)
                fwrite(oncotable_result, file.path(pair_outdir, "oncotable.txt"))
                return(list(index = i, path = oncotable_path))
            }

        }, error = function(e) {
            msg <- sprintf("Unexpected error processing %s: %s", cohort$inputs[i]$pair, e$message)
            warning(msg)
        })
    }, mc.cores = cores)

    # Update oncotable paths in the cohort
    results <- Filter(Negate(is.null), results)
    for (result in results) {
        updated_cohort$inputs[result$index, oncotable := result$path]
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
#' @param temp_fix TRUE/FALSE whether to apply temporary fix
#' @param return_table TRUE/FALSE whether to return the data.table
#' @return data.table or NULL
#' @export
create_filtered_events <- function(
    pair,
    oncotable,
    jabba_gg,
    out_file,
    temp_fix = FALSE,
    return_table = FALSE) {

    ot <- readRDS(oncotable)

    # add a fusion_gene_coords column of NAs if no fusions
    if (!"fusion_gene_coords" %in% colnames(ot)) {
        ot[, fusion_genes := NA]
        ot[, fusion_gene_coords := NA]
    }
    snvs <- ot[grepl("frameshift|missense|stop|disruptive", annotation, perl = TRUE)]
    snvs <- snvs[order(is.na(tier)), ]
    snvs <- snvs[!duplicated(variant.p), ]
    homdels <- ot[type == "homdel"][, vartype := "HOMDEL"][, type := "SCNA"]
  	amps <- ot[type == "amp"][, vartype := "AMP"][, type := "SCNA"]
    fusions <- ot[type == "fusion"]
    possible_drivers <- rbind(snvs, homdels, amps, fusions)
    filtered_events_columns <- c("gene", "fusion_genes", "id", "vartype", "type", "variant.g", "variant.p", "gene_location", "fusion_gene_coords")
    if ("tier" %in% colnames(possible_drivers)) {
        filtered_events_columns <- c(filtered_events_columns, "tier", "therapeutics", "resistances", "diagnoses", "prognoses")
    }
    if ("total_copies" %in% colnames(possible_drivers)) {
        filtered_events_columns <- c(filtered_events_columns, "total_copies")
    }
    
    res <- possible_drivers[, ..filtered_events_columns]
    oncotable_col_to_filtered_events_col <- c(
        "gene" = "gene",
        "fusion_genes" = "fusion_genes",
        "id" = "id",
        "vartype" = "vartype",
        "type" = "type",
        "variant.g" = "Variant_g",
        "variant.p" = "Variant",
        "gene_location" = "Genome_Location",
        "fusion_gene_coords" = "fusion_gene_coords",
        "tier" = "Tier",
        "therapeutics" = "therapeutics",
        "resistances" = "resistances",
        "diagnoses" = "diagnoses",
        "prognoses" = "prognoses",
        "total_copies" = "total_copies"
    )
    intersected_columns <- intersect(filtered_events_columns, names(res))
    setnames(res, old = intersected_columns, new = oncotable_col_to_filtered_events_col[intersected_columns])

  res <- res %>% unique(., by = c("gene", "vartype", "Variant"))
  browser()
    if (nrow(res) > 0) {
        res[, seqnames := tstrsplit(Genome_Location, ":", fixed = TRUE, keep = 1)]
        res[, start := tstrsplit(Genome_Location, "-", fixed = TRUE, keep = 1)]
        res[, start := tstrsplit(start, ":", fixed = TRUE, keep = 2)]
        res[, end := tstrsplit(Genome_Location, "-", fixed = TRUE, keep = 2)]
        res.mut <- res[!is.na(Variant), ]
        if (nrow(res.mut) > 0) {
            #res.mut[, Variant := gsub("p.", "", Variant)]
            res.mut[, vartype := "SNV"]
            res.mut[type=="trunc", vartype := "DEL"]
        }
        res.cn <- res[is.na(Variant) & !is.na(Genome_Location), ]
        if (nrow(res.cn) > 0) {
            jab <- readRDS(jabba_gg)
            res.cn.gr <- GRanges(res.cn)
            res.cn.gr <- gr.val(res.cn.gr, jab$nodes$gr, c("cn", "cn.low", "cn.high"), na.rm = T)
            res.cn.dt <- as.data.table(res.cn.gr)
            res.cn.dt[, total_copies := cn]
            res.cn.dt[, Variant := vartype]
            #res.cn.dt[!is.na(cn) & !is.na(cn.low) & !is.na(cn.high),
            #          Variant := paste0("Total CN:", round(cn, digits = 3), "; CN Minor:", round(cn.low, digits = 3), "; CN Major:", round(cn.high, digits = 3))]
            #res.cn.dt[!is.na(cn) & is.na(cn.low) & is.na(cn.high), Variant := paste0("Total CN:", round(cn, digits = 3))]
            if (temp_fix) {
                res.cn.dt <- res.cn.dt[!(type == "homdel" & cn != 0), ]
                res.cn.dt <- res.cn.dt[!(type == "amp" & cn <= 2), ]
            }
            res.cn.dt[, c("cn", "cn.high", "cn.low", "width", "strand") := NULL] # make null, already added to Variant
            res.final <- rbind(res.mut, res.cn.dt)
        } else {
            res.final <- res.mut
            res.final[, c("seqnames", "start", "end") := NULL]
        }
        write_json(res.final, out_file, pretty = TRUE)
        res.final[, sample := pair]
        if (return_table) {
            return(res.final)
        }
    }
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
lift_filtered_events <- function(cohort, output_data_dir, cores = 1) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }
    
    # Validate required columns exist
    required_cols <- c("pair", "oncotable", "jabba_gg")
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
        
        out_file <- file.path(pair_dir, "filtered.events.json")
        
        tryCatch({
            create_filtered_events(
                pair = row$pair,
                oncotable = row$oncotable,
                jabba_gg = row$jabba_gg,
                out_file = out_file,
                temp_fix = FALSE,
                return_table = FALSE
            )
        }, error = function(e) {
            warning(sprintf("Error processing %s: %s", row$pair, e$message))
        })
    }, mc.cores = cores)
    
    invisible(NULL)
}
