#' @importFrom gUtils dt2gr
#' @importFrom magrittr `%>%`
#' @useDynLib Skilift, .registration=TRUE
#' @import VariantAnnotation
#' @import skidb
#' @import Biostrings
#' @import skitools

internal_settings_path <- system.file("extdata", "test_data", "settings.json", package = "Skilift")

#' @name cov2arrowPGV
#' @description
#'
#' Prepares an scatter plot arrow file with coverage info for PGV (https://github.com/mskilab/pgv)
#'
#' @param cov input coverage data (GRanges)
#' @param field which field of the input data to use for the Y axis
#' @param output_file output file path.
#' @param ref the name of the reference to use. If not provided, then the default reference that is defined in the meta.js file will be loaded.
#' @param cov.color.field a field in the input GRanges object to use to determine the color of each point
#' @param overwrite (logical) by default, if the output path already exists, it will not be overwritten.
#' @param meta.js path to JSON file with metadata for PGV (should be located in "public/settings.json" inside the repository)
#' @param bin.width (integer) bin width for rebinning the coverage (default: 1e4)
#' @author Alon Shaiber
cov2arrowPGV <- function(cov,
                         field = "ratio",
                         output_file = "coverage.arrow",
                         ref = "hg19",
                         cov.color.field = NULL,
                         overwrite = TRUE,
                         meta.js = NULL,
                         ...) {
    outdir <- dirname(output_file)
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    if (!file.exists(output_file) | overwrite) {
        if (!requireNamespace("arrow", quietly = TRUE)) {
            stop('You must have the package "arrow" installed in order for this function to work. Please install it.')
        }

        message("Converting coverage format")
        dat <- gGnome:::cov2cov.js(cov,
            meta.js = meta.js, js.type = "PGV", field = field,
            ref = ref, cov.color.field = cov.color.field, ...
        )
        message("Done converting coverage format")

        if (!is.null(cov.color.field)) {
            dat[, color := gGnome:::color2numeric(get(cov.color.field))]
        } else {
            if (!is.null(meta.js)) {
                ref_meta <- gGnome:::get_ref_metadata_from_PGV_json(meta.js, ref)
                setkey(ref_meta, "chromosome")
                dat$color <- gGnome:::color2numeric(ref_meta[dat$seqnames]$color)
            } else {
                # no cov.color.field and no meta.js so set all colors to black
                dat$color <- 0
            }
        }

        outdt <- dat[, .(x = new.start, y = get(field), color)]

        # if there are any NAs for colors then set those to black
        outdt[is.na(color), color := 0]

        # remove NAs
        outdt <- outdt[!is.na(y)]

        # sort according to x values (that is what PGV expects)
        outdt <- outdt[order(x)]

        message("Writing arrow file (using write_feather)")
        arrow_table <- arrow::Table$create(outdt, schema = arrow::schema(x = arrow::float32(), y = arrow::float32(), color = arrow::float32()))
        arrow::write_feather(arrow_table, output_file, compression = "uncompressed")
    } else {
        message('arrow file, "', output_file, '" already exists.')
    }
    return(output_file)
}

#' @name grab.hets
#' @title grab.hets
#'
#' @description
#'
#' returns allele gtrack given sites.txt from het pileup
#'
#' @param agt.fname (character) path to sites.txt
#' @param min.frac (numeric) between 0 and 1, min frequency in normal to count as het site
#' @param max.frac (numeric) between 0 and 1, max frequency in normal to count as het site
#'
#' @return allele gTrack
grab.hets <- function(agt.fname = NULL,
                      min.frac = 0.2,
                      max.frac = 0.8) {
    if (is.null(agt.fname) || !file.exists(agt.fname)) {
        stop("agt.fname does not exist")
    }

    ## prepare and filter
    agt.dt <- fread(agt.fname)[alt.frac.n > min.frac & alt.frac.n < max.frac, ]
    ## add major and minor
    agt.dt[, which.major := ifelse(alt.count.t > ref.count.t, "alt", "ref")]
    agt.dt[, major.count := ifelse(which.major == "alt", alt.count.t, ref.count.t)]
    agt.dt[, minor.count := ifelse(which.major == "alt", ref.count.t, alt.count.t)]

    ## melt the data frame
    agt.melted <- rbind(
        agt.dt[, .(seqnames, start, end, count = major.count, allele = "major")],
        agt.dt[, .(seqnames, start, end, count = minor.count, allele = "minor")]
    )

    ## make GRanges
    agt.gr <- dt2gr(agt.melted[, .(seqnames, start, end, count, allele)])

    return(agt.gr)
}

#' @name subsample_hetsnps
#' @title subsample_hetsnps
#' @description subsets the hetsnps to masked unique sites and colors by major/minor allele
#'
#' @param het_pileups_wgs sites.txt from het pileup
#' @param mask rds file with masked regions
#' @param sample_size number of snps to randomly sample
#' @export
#' @author Jonathan Rafailov, Shihab Dider
subsample_hetsnps <- function(
    het_pileups_wgs,
    mask = NULL,
    sample_size = 100000) {
    if (is.null(het_pileups_wgs)) {
        stop("het_pileups_wgs does not exist")
    }
    if (is.null(mask)) {
        warning("mask does not exist, using default mask included with package.")
        maska_path <- system.file("extdata", "data", "maskA_re.rds", package = "Skilift")
        maska <- readRDS(maska_path)
    }

    hets.gr <- grab.hets(het_pileups_wgs)
    hets.gr <- gr.val(hets.gr, maska, "mask")
    hets.gr <- hets.gr %Q% (is.na(mask))
    hets.gr$mask <- NULL
    # lets call major blue and minor red
    hets.gr$col <- c("major" = "red", "minor" = "blue")[hets.gr$allele]

    # lets subset a random amount of SNPS so we're under 250k points
    unique.snps <- unique(gr2dt(hets.gr)[, .(seqnames, start, end)])
    n_snps <- nrow(unique.snps)
    message(paste(n_snps, "snps found"))
    if (!is.na(sample_size) && n_snps > sample_size) {
        message(paste("subsampling", sample_size, "points..."))
        snps.to.include <- unique.snps[sample(n_snps, sample_size)] %>% dt2gr()
        subset.hets.gr <- hets.gr %&% snps.to.include
    } else {
        snps.to.include <- unique.snps %>% dt2gr()
        subset.hets.gr <- hets.gr %&% snps.to.include
    }

    return(subset.hets.gr)
}


#' @description
#' Create coverage arrow plot JSON file.
#'
#' @param plot_metadata data.table with Plot metadata for a single plot, columns = (patient.id, source, x (contains path to data file or reference to data file object itself), ref, overwrite).
#' @param datadir Path to data directory with patient directories.
#' @param settings Path to settings.json file.
#'
#' @return NULL.
create_cov_arrow <- function(plot_metadata, datadir, settings = internal_settings_path, color_field = NULL) {
    cov_json_path <- file.path(
        datadir,
        plot_metadata$patient.id,
        plot_metadata$source
    )

    if (!("field" %in% names(plot_metadata))) {
        stop(warning("Please include a 'field' column which indicates the column name that contains the coverage data."))
    }

    # check if the input is a GRanges or a nested list
    is_list <- is(plot_metadata$x[[1]], "list")
    if (is_list) {
        plot_metadata$x <- plot_metadata$x[[1]]
    }

    is_granges <- is(plot_metadata$x[[1]], "GRanges")
    is_path <- is(plot_metadata$x[[1]], "character")
    if (is_granges) {
        # save granges to a temp file and assign path to plot_metadata$x
        print("Saving GRanges to temp file")
        temp_file <- tempfile(fileext = ".rds")
        saveRDS(plot_metadata$x[[1]], temp_file)
        plot_metadata$x <- temp_file
    } else if (is_path) {
        print("Using path to GRanges")
    } else {
        stop(warning("Please provide a GRanges object or a path to a GRanges object."))
    }

    if (!file.exists(cov_json_path) || plot_metadata$overwrite) {
        if (file.exists(plot_metadata$x)) {
            cov2arrowPGV(plot_metadata$x,
                field = plot_metadata$field,
                meta.js = settings,
                ref = plot_metadata$ref,
                output_file = cov_json_path,
                cov.color.field = color_field
            )
        } else {
            warning(
                paste0(
                    "Input coverage file does not exist for name: ",
                    plot_metadata$sample,
                    " so no coverage will be generated."
                )
            )
        }
    } else {
        message(cov_json_path, " already exists! Set overwrite = TRUE if you want to overwrite it.")
    }
}


#' @description
#' Add gGraph JSON file to PGV datafiles.json.
#'
#' @param plot_metadata data.table with Plot metadata, columns = (patient.id, source, x (contains path to data file or reference to data file object itself), ref, overwrite).
#' @param datadir Path to data directory with patient directories.
#' @param settings Path to settings.json file.
#'
#' @return NULL.
create_ggraph_json <- function(plot_metadata, datadir, settings = internal_settings_path) {
    ggraph_json_path <- file.path(
        datadir,
        plot_metadata$patient.id,
        plot_metadata$source
    )
    if (!file.exists(ggraph_json_path) || plot_metadata$overwrite) {
        if (is(plot_metadata$x[[1]], "gGraph")) {
            ggraph <- plot_metadata$x[[1]]
        } else {
            message(paste0("reading in ", plot_metadata$x))
            if (grepl(plot_metadata$x, pattern = ".rds")) {
                ggraph <- readRDS(plot_metadata$x[[1]])
            } else {
                message("Expected .rds ending for gGraph. Attempting to read anyway: ", plot_metadata$x[[1]])
                ggraph <- readRDS(plot_metadata$x[[1]])
            }
        }
        if (any(class(ggraph) == "gGraph")) {
            seq_lengths <- gGnome::parse.js.seqlengths(
                settings,
                js.type = "PGV",
                ref = plot_metadata$ref
            )
            # check for overlap in sequence names
            ggraph.reduced <- ggraph[seqnames %in% names(seq_lengths)]
            if (length(ggraph.reduced) == 0) {
                stop(sprintf(
                    'There is no overlap between the sequence names in the reference
                    used by PGV and the sequences in your gGraph. Here is an
                    example sequence from your gGraph: "%s". And here is an
                    example sequence from the reference used by gGnome.js: "%s"',
                    seqlevels(ggraph$nodes$gr)[1], names(seq_lengths)[1]
                ))
            }
            # add maxcn from plot_metadata if exists
            if ("max.cn" %in% colnames(plot_metadata)) {
                maxcn <- plot_metadata$max.cn
            } else {
                maxcn <- 100
            }
            # sedge.id or other field
            if ("col" %in% names(mcols(ggraph$nodes$gr))) {
                nfields <- "col"
            } else {
                nfields <- NULL
            }
            if ("annotation" %in% colnames(plot_metadata)) {
                annotations <- unlist(plot_metadata$annotation)
            } else {
                annotations <- NULL
            }
            ## if ("annotation" %in% colnames(plot_metadata)) {
            # probably check for other cid.field names?
            # field = 'sedge.id'
            gGnome::refresh(ggraph[seqnames %in% names(seq_lengths)])$json(
                filename = ggraph_json_path,
                verbose = TRUE,
                annotations = annotations,
                maxcn = maxcn,
                nfields = nfields
                # cid.field = field
            )
            ## } else {
            ##     gGnome::refresh(ggraph[seqnames %in% names(seq_lengths)])$json(
            ##                                                                   filename = ggraph_json_path,
            ##                                                                   verbose = TRUE,
            ##                                                                   maxcn = maxcn,
            ##                                                                   nfields = nfields
            ##                                                               )
            ## }
        } else {
            warning(plot_metadata$x, " rds read was not a gGraph")
        }
    } else {
        warning("file ", ggraph_json_path, "already exists. Set overwrite = TRUE if you want to overwrite it.")
    }
}

#' @description
#' Create allelic gGraph JSON file.
#'
#' @param plot_metadata data.table with Plot metadata, columns = (patient.id, source, x (contains path to data file or reference to data file object itself), ref, overwrite).
#' @param datadir Path to data directory with patient directories.
#' @param settings Path to settings.json file.
#'
#' @return NULL.
create_allelic_json <- function(plot_metadata, datadir, settings = internal_settings_path) {
    ggraph_json_path <- file.path(
        datadir,
        plot_metadata$patient.id,
        plot_metadata$source
    )
    if (!file.exists(ggraph_json_path) || plot_metadata$overwrite) {
        if (is(plot_metadata$x[[1]], "gGraph")) {
            ggraph <- plot_metadata$x[[1]]
        } else {
            message(paste0("reading in ", plot_metadata$x))
            if (grepl(plot_metadata$x, pattern = ".rds")) {
                ggraph <- readRDS(plot_metadata$x[[1]])
            } else {
                message("Expected .rds ending for gGraph. Attempting to read anyway: ", plot_metadata$x)
                ggraph <- readRDS(plot_metadata$x[[1]])
            }
        }
        if (any(class(ggraph) == "gGraph")) {
            seq_lengths <- gGnome::parse.js.seqlengths(
                settings,
                js.type = "PGV",
                ref = plot_metadata$ref
            )
            # check for overlap in sequence names
            ggraph.reduced <- ggraph[seqnames %in% names(seq_lengths)]
            if (length(ggraph.reduced) == 0) {
                stop(sprintf(
                    'There is no overlap between the sequence names in the reference
                    used by PGV and the sequences in your gGraph. Here is an
                    example sequence from your gGraph: "%s". And here is an
                    example sequence from the reference used by gGnome.js: "%s"',
                    seqlevels(ggraph$nodes$gr)[1], names(seq_lengths)[1]
                ))
            }
            # add maxcn from plot_metadata if exists
            if ("max.cn" %in% colnames(plot_metadata)) {
                maxcn <- plot_metadata$max.cn
            } else {
                maxcn <- 100
            }
            if ("annotation" %in% colnames(plot_metadata)) {
                annotations <- unlist(plot_metadata$annotation)
            } else {
                annotations <- NULL
            }


            # sedge.id or other field
            ## if ("annotation" %in% colnames(plot_metadata)) {
            # probably check for other cid.field names?
            # field = 'sedge.id'
            gGnome::refresh(ggraph[seqnames %in% names(seq_lengths)])$json(
                filename = ggraph_json_path,
                verbose = TRUE,
                annotations = annotations,
                maxcn = maxcn,
                nfields = "col",
                save = TRUE
                # offset = TRUE
                # cid.field = field
            )
            ## for(x in 1:length(gg.js$intervals)) {
            ##     col1 = gg.js$intervals[[x]]$metadata$color
            ##     if(gg.js$intervals[[x]]$y != 0) {
            ##         if(col1 == "#0000FF80") {
            ##             gg.js$intervals[[x]]$y = gg.js$intervals[[x]]$y + 0.1
            ##         } else if(col1 == "#FF000080") {
            ##             gg.js$intervals[[x]]$y = gg.js$intervals[[x]]$y - 0.1
            ##         }
            ##     }
            ## }

            ## } else {
            ##     gg.js = gGnome::refresh(ggraph[seqnames %in% names(seq_lengths)])$json(
            ##                                                                   verbose = TRUE,
            ##                                                                   maxcn = maxcn,
            ##                                                                   nfields = "col",
            ##                                                                   save = FALSE
            ##                                                                   )
            ##     ## filename = ggraph_json_path,
            ## }
            ## temporary fix for adding padding to allelic graphs for major and minor- should probably be implented into ggnome
            # major : "#0000FF80"
            # minor : "#FF000080"
        } else {
            warning(plot_metadata$x, " rds read was not a gGraph")
        }
    } else {
        warning("file ", ggraph_json_path, "already exists. Set overwrite = TRUE if you want to overwrite it.")
    }
}

#' @description
#' Create gWalk JSON file
#'
#' @param plot_metadata data.table with Plot metadata, columns = (patient.id, source, x (contains path to data file or reference to data file object itself), ref, overwrite).
#' @param datadir Path to data directory with patient directories.
#' @param settings Path to settings.json file.
#'
#' @return NULL.
create_gwalk_json <- function(plot_metadata, datadir, settings = internal_settings_path) {
    gwalk_json_path <- file.path(datadir, plot_metadata$patient.id, plot_metadata$source)
    if (!file.exists(gwalk_json_path) || plot_metadata$overwrite == TRUE) {
        if (is(plot_metadata$x[[1]], "gWalk")) {
            gwalk <- plot_metadata$x[[1]] %>% gGnome::refresh()
        } else {
            message(paste0("reading in ", plot_metadata$x))
            if (grepl(plot_metadata$x, pattern = ".rds")) {
                gwalk <- readRDS(plot_metadata$x) %>% gGnome::refresh()
            } else {
                message("Expected .rds ending for gWalk Attempting to read anyway: ", plot_metadata$x)
                gwalk <- readRDS(plot_metadata$x) %>% gGnome::refresh()
            }
        }

        if (gwalk$length == 0) {
            warning(sprintf("Zero walks in gWalk .rds file provided for sample %s!
                No walks json will be produced!", plot_metadata$sample))
            return(NA)
        }

        gwalk$json(
            filename = gwalk_json_path, verbose = TRUE,
            annotations = unlist(plot_metadat$annotation),
            include.graph = FALSE
        )
    } else {
        message(gwalk_json_path, "already exists! Set overwrite = TRUE if you want to overwrite it.")
    }
}

#' @description
#' Create mutations gGraph JSON file for SOMATIC MUTATIONS.
#'
#' @param plot_metadata data.table with Plot metadata, columns = (patient.id, source, x (contains path to data file or reference to data file object itself), ref, overwrite).
#' @param datadir Path to data directory with patient directories.
#' @param settings Path to settings.json file.
#'
#' @return NULL.
create_somatic_json <- function(plot_metadata, datadir, settings = internal_settings_path) {
    message("using function create_somatic_json")
    somatic_json_path <- file.path(
        datadir,
        plot_metadata$patient.id,
        plot_metadata$source
    )
    if (!file.exists(somatic_json_path) || plot_metadata$overwrite) {
        if (any(class(plot_metadata$x[[1]]) == "data.table")) {
            mutations.dt <- plot_metadata$x[[1]]
        } else {
            message(paste0("reading in ", plot_metadata$x))
            if (grepl(plot_metadata$x, pattern = ".rds")) {
                mutations.dt <- as.data.table(readRDS(plot_metadata$x[[1]]))
            } else {
                message("Expected .rds ending for mutations. Attempting to read anyway: ", plot_metadata$x)
                mutations.dt <- as.data.table(readRDS(plot_metadata$x))
            }
        }
        if (any(class(mutations.dt) == "data.table")) {
            seq_lengths <- gGnome::parse.js.seqlengths(
                settings,
                js.type = "PGV",
                ref = plot_metadata$ref
            )
            # check for overlap in sequence names
            mutations.reduced <- mutations.dt[seqnames %in% names(seq_lengths), ]
            if (length(mutations.reduced) == 0) {
                stop(sprintf(
                    'There is no overlap between the sequence names in the reference
                    used by PGV and the sequences in your mutations. Here is an
                    example sequence from your mutations: "%s". And here is an
                    example sequence from the reference used by gGnome.js: "%s"',
                    mutations$seqnames[1], names(seq_lengths)[1]
                ))
            }
            yfield <- plot_metadata$field[1]
            # coerce vaf col name to lower case to avoid col name mismatch for vaf/VAF column
            setnames(mutations.dt, old = "VAF", new = "vaf", skip_absent = TRUE)
            mutations.dt <- mutations.dt[!is.na(get(yfield)), ]
            mutations.dt[start == end, end := end + 1]
            mutations.dt[, strand := NULL]
            # create an empty mutation annotation string, then add attributes to it if they exist
            mut_ann <- ""
            if ("annotation" %in% colnames(mutations.dt)) {
                mut_ann <- paste0("Type: ", mutations.dt$annotation, "; ")
            }
            if ("gene" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Gene: ", mutations.dt$gene, "; ")
            }
            if ("variant.c" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Variant: ", mutations.dt$variant.c, "; ")
            }
            if ("variant.p" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Protein_variant: ", mutations.dt$variant.p, "; ")
            }
            if ("variant.g" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Genomic_variant: ", mutations.dt$variant.g, "; ")
            }
            if ("vaf" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "VAF: ", mutations.dt$vaf, "; ")
            }
            if ("alt" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Alt_count: ", mutations.dt$alt, "; ")
            }
            if ("ref" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Ref_count: ", mutations.dt$ref, "; ")
            }
            if ("normal.alt" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Normal_alt_count: ", mutations.dt$normal.alt, "; ")
            }
            if ("normal.ref" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Normal_ref_count: ", mutations.dt$normal.ref, "; ")
            }
            if ("FILTER" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Filter: ", mutations.dt$FILTER, "; ")
            }
            mutations.dt[, annotation := mut_ann]
            dt2json_mut(
                dt = mutations.dt,
                ref = plot_metadata$ref,
                settings = settings,
                meta_data = c("gene", "feature_type", "annotation", "REF", "ALT", "variant.c", "variant.p", "vaf", "transcript_type", "impact", "rank"),
                y_col = yfield,
                file_name = somatic_json_path
            )
        } else {
            warning(plot_metadata$x, " rds read was not mutations")
        }
    } else {
        warning("file ", somatic_json_path, "already exists. Set overwrite = TRUE if you want to overwrite it.")
    }
}

#' @description
#' Create mutations gGraph JSON file for GERMLINE MUTATIONS.
#' This function is unique from create_somatic_json in that it will subset only important mutations by MODIFIER code.
#'
#' @param plot_metadata data.table with Plot metadata, columns = (patient.id, source, x (contains path to data file or reference to data file object itself), ref, overwrite).
#' @param datadir Path to data directory with patient directories.
#' @param settings Path to settings.json file.
#'
#' @return NULL.
create_germline_json <- function(plot_metadata, datadir, settings = internal_settings_path) {
    message("using function create_germline_json")
    germline_json_path <- file.path(
        datadir,
        plot_metadata$patient.id,
        plot_metadata$source
    )
    if (!file.exists(germline_json_path) || plot_metadata$overwrite) {
        if (any(class(plot_metadata$x[[1]]) == "data.table")) {
            mutations.dt <- plot_metadata$x[[1]]
        } else {
            message(paste0("reading in ", plot_metadata$x))
            if (grepl(plot_metadata$x, pattern = ".rds")) {
                mutations.dt <- as.data.table(readRDS(plot_metadata$x))
            } else {
                message("Expected .rds ending for mutations. Attempting to read anyway: ", plot_metadata$x)
                mutations.dt <- as.data.table(readRDS(plot_metadata$x))
            }
        }
        if (any(class(mutations.dt) == "data.table")) {
            seq_lengths <- gGnome::parse.js.seqlengths(
                settings,
                js.type = "PGV",
                ref = plot_metadata$ref
            )
            # check for overlap in sequence names
            mutations.reduced <- mutations.dt[seqnames %in% names(seq_lengths), ]
            if (length(mutations.reduced) == 0) {
                stop(sprintf(
                    'There is no overlap between the sequence names in the reference
                    used by PGV and the sequences in your mutations. Here is an
                    example sequence from your mutations: "%s". And here is an
                    example sequence from the reference used by gGnome.js: "%s"',
                    mutations$seqnames[1], names(seq_lengths)[1]
                ))
            }
            yfield <- plot_metadata$field[1]
            # coerce vaf col name to lower case to avoid col name mismatch for vaf/VAF column
            setnames(mutations.dt, old = "VAF", new = "vaf", skip_absent = TRUE)
            mutations.dt <- mutations.dt[!is.na(get(yfield)), ]
            mutations.dt[start == end, end := end + 1]
            mutations.dt[, strand := NULL]
            # create an empty mutation annotation string, then add attributes to it if they exist
            mut_ann <- ""
            if ("annotation" %in% colnames(mutations.dt)) {
                mut_ann <- paste0("Type: ", mutations.dt$annotation, "; ")
            }
            if ("gene" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Gene: ", mutations.dt$gene, "; ")
            }
            if ("variant.c" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Variant: ", mutations.dt$variant.c, "; ")
            }
            if ("variant.p" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Protein_variant: ", mutations.dt$variant.p, "; ")
            }
            if ("variant.g" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Genomic_variant: ", mutations.dt$variant.g, "; ")
            }
            if ("vaf" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "VAF: ", mutations.dt$vaf, "; ")
            }
            if ("alt" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Alt_count: ", mutations.dt$alt, "; ")
            }
            if ("ref" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Ref_count: ", mutations.dt$ref, "; ")
            }
            if ("normal.alt" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Normal_alt_count: ", mutations.dt$normal.alt, "; ")
            }
            if ("normal.ref" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Normal_ref_count: ", mutations.dt$normal.ref, "; ")
            }
            if ("FILTER" %in% colnames(mutations.dt)) {
                mut_ann <- paste0(mut_ann, "Filter: ", mutations.dt$FILTER, "; ")
            }
            mutations.dt[, annotation := mut_ann]
            dt2json_mut(
                dt = mutations.dt,
                ref = plot_metadata$ref,
                settings = settings,
                meta_data = c("gene", "feature_type", "annotation", "REF", "ALT", "variant.c", "variant.p", "vaf", "transcript_type", "impact", "rank"),
                y_col = yfield,
                file_name = germline_json_path
            )
        } else {
            warning(plot_metadata$x, " rds read was not mutations")
        }
    } else {
        warning("file ", germline_json_path, "already exists. Set overwrite = TRUE if you want to overwrite it.")
    }
}


#' @description
#' Create ppfit gGraph JSON file.
#'
#' @param plot_metadata data.table with Plot metadata, columns = (patient.id, source, x (contains path to data file or reference to data file object itself), ref, overwrite).
#' @param datadir Path to data directory with patient directories.
#' @param settings Path to settings.json file.
#' @param cov_path Path to coverage file (necessary for segstats information, optional).
#'
#' @return NULL.
create_ppfit_genome_json <- function(plot_metadata, datadir, settings = internal_settings_path, cov_path = NULL) {
    ppfit_json_path <- file.path(
        datadir,
        plot_metadata$patient.id,
        plot_metadata$source
    )
    tryCatch(
        {
            if (!file.exists(ppfit_json_path) || plot_metadata$overwrite) {
                if (is(plot_metadata$x[[1]], "list")) {
                    ggraph <- plot_metadata$x[[1]]
                } else {
                    message(paste0("reading in ", plot_metadata$x))
                    if (grepl(plot_metadata$x, pattern = ".rds")) {
                        ggraph <- readRDS(plot_metadata$x[[1]])
                    } else {
                        message("Expected .rds ending for gGraph. Attempting to read anyway: ", plot_metadata$x[[1]])
                        ggraph <- readRDS(plot_metadata$x[[1]])
                    }
                }
                if (any(class(ggraph) == "gGraph")) {
                    seq_lengths <- gGnome::parse.js.seqlengths(
                        settings,
                        js.type = "PGV",
                        ref = plot_metadata$ref
                    )
                    colnames_check <- c(
                        "start_ix",
                        "end_ix",
                        "eslack_in",
                        "eslack_out",
                        "edges_in",
                        "edges_out",
                        "tile_id",
                        "snode_id",
                        "loose_left",
                        "loose_right",
                        "loose_cn_left",
                        "loose_cn_right",
                        "node_id",
                        "raw_mean",
                        "raw_var",
                        "nbins",
                        "nbins_tot",
                        "nbins_nafrac",
                        "wbins_nafrac",
                        "wbins_ok",
                        "mean",
                        "bad",
                        "max_na",
                        "loess_var",
                        "tau_sq_post",
                        "post_var",
                        "var",
                        "sd"
                    )

                    if (all(colnames_check %in% names(ggraph$nodes$dt))) {
                        ggraph2 <- ggraph
                        fields.keep <- c(
                            colnames_check,
                            "seqnames",
                            "start",
                            "end",
                            "strand",
                            "width",
                            "loose",
                            "index"
                        )
                    } else {
                        ## segstats information
                        segstats.dt <- create_ppfit_json(
                            balanced_gg = ggraph,
                            cov_path = cov_path,
                            path_obj = plot_metadata$x,
                            return_table = TRUE,
                            write_json = FALSE
                        )
                        segstats.gr <- GRanges(segstats.dt, seqlengths = seq_lengths) %>% trim()
                        ggraph2 <- gG(nodes = segstats.gr, edges = ggraph$edges$dt)
                        fields.keep <- names(segstats.dt) %>% grep("cn", ., invert = TRUE, value = TRUE)
                    }
                    ggraph2$set(y.field = "cn")
                    ## check for overlap in sequence names
                    ggraph.reduced <- ggraph2[seqnames %in% names(seq_lengths)]
                    if (length(ggraph.reduced) == 0) {
                        stop(sprintf(
                            'There is no overlap between the sequence names in the reference
                    used by PGV and the sequences in your gGraph. Here is an
                    example sequence from your gGraph: "%s". And here is an
                    example sequence from the reference used by gGnome.js: "%s"',
                            seqlevels(ggraph.reduced$nodes$gr)[1], names(seq_lengths)[1]
                        ))
                    }
                    # sedge.id or other field
                    if ("annotation" %in% colnames(plot_metadata)) {
                        # probably check for other cid.field names?
                        # field = 'sedge.id'
                        ## gGnome::refresh(ggraph[seqnames %in% names(seq_lengths)])$json(
                        gGnome::refresh(ggraph.reduced)$json(
                            filename = ppfit_json_path,
                            verbose = TRUE,
                            annotations = unlist(plot_metadata$annotation),
                            maxcn = 500,
                            nfields = fields.keep,
                            save = TRUE
                            # cid.field = field
                        )
                    } else {
                        ## gGnome::refresh(ggraph[seqnames %in% names(seq_lengths)])$json(
                        gGnome::refresh(ggraph.reduced)$json(
                            filename = ppfit_json_path,
                            verbose = TRUE,
                            maxcn = 500,
                            nfields = fields.keep,
                            save = TRUE
                        )
                    }
                } else {
                    warning(plot_metadata$x, " rds read was not a gGraph")
                }
            } else {
                warning("file ", ppfit_json_path, "already exists. Set overwrite = TRUE if you want to overwrite it.")
            }
        },
        error = function(e) {
            message("Error in creating ppfit plot for sample: ", plot_metadata$patient.id, "\n JaBbA:::segstats needs to be run or provided first")
            print(e)
        }
    )
}

#' @name create_somatic_json
#' @title create_somatic_json
#' @description
#'
#' @param somatic_snv_cn somatic copy number data table
#' @param out_file where the json should be written
#' @param pair patient.id to be used for pgvdb or case reports
#' @param pgv_settings the settings file used for pgvdb or case reports
#' @param return_table_pgv TRUE/FALSE whether to return the data.table to add to pgvdb
#' @param meta_keep meta data to write in the json
#' @param y_col column to be used for y axis values
#' @return data.table to add to pgv or NULL, depends on return_table_pgv
#' @export
#' @author Stanley Clarke
## create_somatic_json = function(somatic_snv_cn, out_file, pair, pgv_settings, return_table_pgv = FALSE, meta_keep = NULL, y_col = "est_cn_llrm", ref = "hg19") {
##     som.dt = readRDS(somatic_snv_cn)
##     som.dt = som.dt[!is.na(get(y_col)),]
##     som.dt[start == end, end := end +1]
##     som.dt[, strand := NULL]
##     som.dt[variant.p != "",annotation := paste0("Type: ", annotation, "; Gene: ", gene, "; Variant: ",variant.c, "; Protein_variant: ", variant.p, "; VAF: ",vaf)]
##     som.dt[variant.p == "",annotation := paste0("Type: ", annotation, "; Gene: ", gene, "; Variant: ",variant.c, "; VAF: ",vaf)]
##     dt2json_mut(dt = som.dt, ref = ref,settings = pgv_settings, meta_data = meta_keep, y_col = y_col,
##                 file_name = out_file)
##     if(return_table_pgv) {
##         dt.add = data.table(patient.id = pair, type = "genome",visible = TRUE, title = "Copy Number Mutations", source = "mutations.json")
##     }
## }

#' @name dt2json_mut
#' @title dt2json_mut
#' @description
#'
#' function to create a mutation json, used in create_somatic_json
#'
#' @param dt data.table with seqnames,start and end
#' @param patient.id patient id to be added to pgvdb or case reports
#' @param ref reference for pgv or case reports
#' @param file_name the file the json should be written to
#' @param meta_data extra information to add to the json
#' @param y_col column to be used for y axis values
#' @return NULL
#' @export
#' @author Stanley Clarke
dt2json_mut <- function(dt, patient.id, ref, settings, file_name = paste(getwd(), "test.json", sep = "/"), meta_data = NULL, y_col = NULL) {
    # create vector of seqlengths
    settings_data <- jsonlite::fromJSON(settings)
    chrom_lengths <- as.data.table(settings_data$coordinates$sets[[ref]])[, .(chromosome, startPoint, endPoint)]
    colnames(chrom_lengths) <- c("seqnames", "start", "end")

    if (nrow(chrom_lengths[grepl("chr", seqnames), ]) > 0) {
        chrom_lengths[!grepl("chr", seqnames), seqnames := paste0("chr", seqnames)] # weird fix because hg38_chr does not have chr on Y and M
    }
    # add y value specified
    if (is.null(y_col)) {
        dt$y_value <- 1
    } else {
        dt[, y_value := get(y_col)]
    }
    # convert to ggraph and create json
    if (nrow(chrom_lengths[grepl("chr", seqnames), ]) > 0) {
        gr1 <- dt2gr(dt[order(seqnames, start), ]) %>%
            sortSeqlevels() %>%
            gr.chr()
    } else {
        gr1 <- dt2gr(dt[order(seqnames, start), ]) %>%
            sortSeqlevels() %>%
            gr.nochr()
    }
    if (any(gr1@seqinfo@seqlengths > chrom_lengths[seqnames %in% names(seqlengths(gr1))]$end)) {
        stop(paste("the seqlengths of your granges has ranges that are not contained in the seqlengths of", ref))
    }
    jab <- gG(nodes = gr1)
    settings_y <- list(y_axis = list(
        title = "copy number",
        visible = TRUE
    ))
    node.dt <- gr2dt(jab$nodes$gr[, c("snode.id", "y_value", meta_data)])
    node.json <- node.dt[, .(chromosome = seqnames, startPoint = start, endPoint = end, iid = snode.id, title = snode.id, type = "interval", y = y_value, annotation = node.dt$annotation)]
    gg.js <- list(intervals = node.json, connections = data.table())
    gg.js <- c(list(settings = settings_y), gg.js)
    message(paste0("Writing json to ", file_name))
    jsonlite::write_json(gg.js, file_name,
        pretty = TRUE, auto_unbox = TRUE, digits = 4
    )
}

#' @name sigprofiler_decomposed_probs_json
#' @title sigprofiler_decomposed_probs_json
#' @description
#'
#' function to create decomposed probabilities json for case reports as adapted to sigprofiler outputs
#'
#' @param probs path to sigprofiler decomposed probabilities
#' @param is_indel Boolean to indicate if the input matrix is for indels
#' @param data_dir path to data directory of which each subdirectory contains jsons per patient
#' @param pairs vector of samples for which to create jsons, by default = NULL, indicating process all patients
#' @param cores number of cores to parallelize with
#' @return data.table or NULL
#' @export
#' @author Johnathan Rafailov
sigprofiler_decomposed_probs_json <- function(probs, is_indel, data_dir, pairs, cores = 1) {
    fread(probs) -> decomposed.probs
    ### lets change first column to something more predictable
    colnames(decomposed.probs)[1] <- "samples"
    ### seems like _somatic is systematically ended to the end of all sample names in the beginning
    decomposed.probs[, samples := gsub("_somatic", "", samples)]
    if (!is.null(pairs)) {
        decomposed.probs <- decomposed.probs[samples %in% pairs, ]
    }
    decomposed.probs.per.sample <- melt(decomposed.probs,
        measure.vars = c(3:ncol(decomposed.probs)),
        variable.name = "signature",
        value.name = "p"
    )
    samples <- unique(decomposed.probs.per.sample$samples)
    mclapply(samples, function(pair) {
        pair_data <- decomposed.probs.per.sample[samples == pair]
        if (is_indel) {
            samp_data <- data.frame(
                signature = pair_data[, signature],
                insdel = pair_data[[2]],
                p = pair_data[, p]
            )
            catalog_file_name <- "id_decomposed_prob.json"
        } else {
            samp_data <- data.frame(
                signature = pair_data[, signature],
                tnc = pair_data[[2]],
                p = pair_data[, p]
            )
            catalog_file_name <- "sbs_decomposed_prob.json"
        }
        write_json(samp_data, paste0(data_dir, "/", pair, "/", catalog_file_name))
    }, mc.cores = cores)
}

#' @title process_gencode
#' @description
#'
#' Helper script to process gencode parameter
#'
#' @param gencode path to gencode file. Gencode file must be either rds or some format accepted by rtracklayer::import (e.g. GTF)
#' @return gencode_gr GRanges
#' @author Marcin Imielinski
process_gencode <- function(gencode = NULL) {
    if (is.null(gencode)) {
        gencode <- skidb::read_gencode()
        print(gencode)
    } else if (is.character(gencode)) {
        if (grepl(".rds$", gencode)) {
            gencode <- readRDS(gencode)
        } else {
            gencode <- rtracklayer::import(gencode)
        }
    }
    return(gencode)
}

#' @name grok_bcf
#' @rdname
#' @title Reads and parses bcf via bcftools call
#' @param bcf path to bcf file
#' @param gr optional granges to query
#' @param bpath path to bcftools binary executable
#' @export
grok_bcf <- function(bcf, gr = NULL, bpath = "/nfs/sw/bcftools/bcftools-1.9/bin/bcftools", label = NA, filter = "PASS", snv = FALSE, indel = FALSE, het = FALSE, hom = FALSE, keep.modifier = TRUE, long = FALSE, oneliner = FALSE, verbose = FALSE) {
    cmd <- sprintf("%s view %s", bpath, bcf)

    if (is.na(label)) {
        label <- bcf
    }

    if (!is.null(gr)) {
        wins <- paste(gr.string(gr.stripstrand(gr)), collapse = ",")
        cmd <- paste(cmd, "-r", wins)
    }

    if (!is.null(filter) && !is.na(filter)) {
        cmd <- paste(cmd, sprintf('-i \'FILTER="%s"\'', filter))
    }

    if (het) {
        cmd <- paste(cmd, "-g het")
    }

    if (indel) {
        cmd <- paste(cmd, "-v indels")
    }

    if (snv) {
        cmd <- paste(cmd, "-v snps")
    }

    if (het) {
        cmd <- paste(cmd, "-g het")
    }

    if (hom) {
        cmd <- paste(cmd, "-g hom")
    }

    ## quick dunlist
    .dunlist <- function(x) {
        ## simplified dunlist to output integer listid and also listiid
        out <- data.table(listid = rep(1:length(x), elementNROWS(x)), V1 = unlist(x))[, listiid := 1:.N, by = listid]
        return(out)
    }

    if (verbose) {
        message(cmd)
    }

    p <- pipe(cmd)
    lines <- readLines(p)
    close(p)

    is.header <- grepl("^\\#", lines)
    header <- lines[is.header]
    contigs <- strsplit(gsub("^\\#\\#", "", grep("contig", header, value = TRUE)), ",")
    sl <- suppressWarnings(structure(
        names = gsub(".*ID=\\<", "", sapply(contigs, "[", 1)),
        as.numeric(gsub(">$", "", gsub(".*length=\\<", "", sapply(contigs, "[", 2))))
    ))

    other <- lines[!is.header]
    if (length(other)) {
        out <- fread(paste(c(other, ""), collapse = "\n"), sep = "\t", header = FALSE)
        sn <- unlist(strsplit(gsub("^\\#", "", header[length(header)]), "\\t"))
        sfields <- sn[-c(1:9)]
        setnames(out, sn)
        out[, seqnames := as.character(CHROM)]
        out[, start := POS]
        out[, end := POS]
        out[, listid := 1:.N] ## set listid to keep track of lists
        ## unpack bcf "format" + sample fields


        if (!is.null(out$FORMAT)) {
            fdat <- .dunlist(strsplit(out$FORMAT, ":"))
            setnames(fdat, 2, "field")
            #  out$FORMAT = NULL ### keep in for now for sanity checks
            for (sfield in sfields) ## can be more than one sample field
            {
                fdatm <- fdat %>% merge(.dunlist(strsplit(out[[sfield]], ":")), by = c("listid", "listiid")) ## merge on both listid and listiid
                fdatc <- dcast.data.table(copy(fdatm)[, field := paste(sfield, field, sep = "_")], listid ~ field, value.var = "V1")
                out <- merge(out, fdatc, by = "listid", all.x = TRUE) ## order of out should be maintained here since keyed by listid which (now) is an integer
            }
        }

        if (!is.null(out$INFO)) {
            ## unpack "info" field
            idat <- .dunlist(strsplit(out$INFO, ";"))
            idat <- cbind(idat, colsplit(idat$V1, pattern = "=", names = c("field", "value")))
            idatc <- dcast.data.table(idat, listid ~ field, value.var = "value")
            out$INFO <- NULL
            mcols <- setdiff(names(idatc), c("REF", "ALT"))
            out <- merge(out, idatc[, mcols, with = FALSE], by = "listid", all.x = TRUE) ##
        }

        out <- dt2gr(out, seqlengths = sl)
        out <- grok_vcf(out, keep.modifier = keep.modifier, long = long, oneliner = oneliner, verbose = verbose, label = label)
    } else {
        out <- GRanges(seqlengths = sl)
    }
    return(out)
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

#' Parse oncokb outputs and tier
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


#' @name oncotable
#' @title oncotable
#' @description
#'
#' THIS ONCOTABLE IS ADAPTED FOR SKILIFT
#' in that it will report only the putatively impacted transcript based on SnpEff outputs to remove duplicated SNV entries
#'
#' Takes as input (keyed) "tumors" (aka pairs) table which a metadata table with specific
#' columns pointing to paths corresponding to one or more of the following pipeline outputs:
#'
#' $annotated_bcf  Path to annotated.bcf file that is the primary output of SnpEff module from which TMB and basic mutation
#' descriptions are extracted along with their basic categories (these will comprising the core of the oncoplot are computed)
#'
#' $fusions  Path to fusion.rds file that is the primary output of the Fusions modjle, from which protein coding fusions will
#' be computed for
#'
#' $jabba_rds  Path to jabba.simple.rds output representing primary output of JaBbA module from which SCNA and overall
#' junction burden are computed
#'
#' $complex    Path to complex.rds gGnome cached object that is the primary output of Events module, from which simple
#' and complex event burdens are computed
#'
#' $signature_counts Path to signature_counts.txt that is the primary output of Signatures module from which SNV signature
#' counts are computed
#'
#' The function then outputs a melted data.table of "interesting" features that can be saved and/or immediately output
#' into oncoprint.  This data.table will at the very least have fields $id $type (event type), $track, and  $source
#' populated in addition to a few other data type specific columns.
#'
#' The $source column is the name of the column of tumors from which that data was extracted, and track is a grouping
#' variable that allows separation of the various data types.
#'
#' All the paths above can be NA or non existent, in which case a dummy row is inserted into the table so that downstream
#' applications know that data is missing for that sample.
#'
#' @param tumors keyed data.table i.e. keyed by unique tumor id with specific columns corresponding to  paths to pipeline outputs(see description)
#' @param gencode path to gencode with just a single entry for each gene (so gencode entries for each gene are collapse to a single range). The input could be .gtf or .rds with GRanges object, or a GRanges object i.e. resulting from importing the (appropriate) GENCODE .gtf via rtracklayer, note: this input is only used in CNA to gene mapping. If nothing is provided then 'http://mskilab.com/fishHook/hg19/gencode.v19.genes.gtf' is used by default.
#' @param amp.thresh SCNA amplification threshold to call an amp as a function of ploidy (4)
#' @param del.thresh SCNA deletion threshold for (het) del as a function of ploidy (by default cn = 1 will be called del, but this allows additoinal regions in high ploidy tumors to be considered het dels)
#' @param mc.cores number of cores for multithreading
#' @param verbose logical flag
#' @author Marcin Imielinski
#' @export
oncotable <- function(tumors, gencode = "http://mskilab.com/fishHook/hg19/gencode.v19.genes.gtf", verbose = TRUE, amp.thresh = 4, filter = "PASS", del.thresh = 0.5, mc.cores = 1) {
    gencode <- process_gencode(gencode)

    if ("type" %in% names(mcols(gencode))) {
        # This is a bit hacky. The hg38 object does not contain the "type" column so we check if it is there and only use it when it is present
        pge <- gencode %Q% (type == "gene" & gene_type == "protein_coding")
    } else {
        pge <- gencode %Q% (gene_type == "protein_coding")
    }

    .oncotable <- function(dat, x = dat[[key(dat)]][1], pge, verbose = TRUE, amp.thresh = 2, del.thresh = 0.5, filter = "PASS") {
        snpeff_ontology <- readRDS(system.file("extdata/data", "snpeff_ontology.rds", package = "Skilift"))
        out <- data.table()

        ## collect gene fusions
        if (!is.null(dat$fusions) && file.exists(dat[x, fusions])) {
            if (verbose) {
                message("pulling $fusions for ", x)
            }
            fus <- readRDS(dat[x, fusions])$meta
            if (nrow(fus)) {
                fus <- fus[silent == FALSE, ][!duplicated(genes), ]
                fus[, vartype := ifelse(in.frame == TRUE, "fusion", "outframe_fusion")] # annotate out of frame fusions
                fus <- fus[, .(
                    gene = strsplit(genes, ",") %>% unlist(),
                    vartype = rep(vartype, sapply(strsplit(genes, ","), length)),
                    fusion_genes = rep(genes, sapply(strsplit(genes, ","), length))
                )][, id := x][, track := "variants"][, type := vartype][, source := "fusions"]
                # get coordinates for fusion genes
                fus[, fusion_gene_coords := unlist(lapply(strsplit(fusion_genes, ","), function(genes) {
                    coords <- lapply(genes, function(gene) {
                        gene_ranges <- pge[mcols(pge)$gene_name == gene]
                        paste0(seqnames(gene_ranges), ":", start(gene_ranges), "-", end(gene_ranges))
                    })
                    paste(unlist(coords), collapse = ",")
                }))]
                out <- rbind(out, fus, fill = TRUE, use.names = TRUE)
            }
        } else { ## signal missing result
            out <- rbind(out, data.table(id = x, type = NA, source = "fusions"), fill = TRUE, use.names = TRUE)
        }

        ## collect complex events
        if (!is.null(dat$complex) && file.exists(dat[x, complex])) {
            if (verbose) {
                message("pulling $complex events for ", x)
            }
            sv <- readRDS(dat[x, complex])$meta$events
            if (nrow(sv)) {
                sv <- sv[, .(value = .N), by = type][, id := x][, track := ifelse(type %in% c("del", "dup", "invdup", "tra", "inv"), "simple sv", "complex sv")][, source := "complex"]
                out <- rbind(out, sv, fill = TRUE, use.names = TRUE)
            }
        } else {
            out <- rbind(out, data.table(id = x, type = NA, source = "complex"), fill = TRUE, use.names = TRUE)
        }

        ## collect copy number / jabba
        if (!is.null(dat$jabba_rds) && file.exists(dat[x, jabba_rds])) {
            if (verbose) {
                message("pulling $jabba_rds to get SCNA and purity / ploidy for ", x)
            }
            jab <- readRDS(dat[x, jabba_rds])
            jabpurity <- jab$meta$purity
            jabploidy <- jab$meta$ploidy
            if (is.null(jab$meta$purity)) jabpurity <- jab$purity
            if (is.null(jab$meta$ploidy)) jabploidy <- jab$ploidy
            out <- rbind(out,
                data.table(id = x, value = c(jabpurity, jabploidy), type = c("purity", "ploidy"), track = "pp"),
                fill = TRUE, use.names = TRUE
            )

            # get the ncn data from jabba
            nseg <- NULL

            ## Don't fail out if karyograph isn't found
            ## Just provide nseg is NULL, and get_gene_ampdels_from_jabba
            ## assumes ncn = 2 (pretty safe assumption)
            if (!is.null(dat$karyograph) && file.exists(dat[x, karyograph])) {
                nseg <- readRDS(dat$karyograph)$segstats[, c("ncn")]
            }

            scna <- Skilift::get_gene_ampdels_from_jabba(jab,
                amp.thresh = amp.thresh,
                del.thresh = del.thresh, pge = pge, nseg = nseg
            )

            if (nrow(scna)) {
                scna[, track := "variants"][, source := "jabba_rds"][, vartype := "scna"]
                out <- rbind(out,
                    scna[, .(id = x, value = min_cn, type, track, gene = gene_name)],
                    fill = TRUE, use.names = TRUE
                )
            }
        } else {
            out <- rbind(out, data.table(id = x, type = NA, source = "jabba_rds"), fill = TRUE, use.names = TRUE)
        }

        if (!is.null(dat$oncokb_cna) && file.exists(dat[x, oncokb_cna])) {
            oncokb_cna <- data.table::fread(dat[x, oncokb_cna])
            concat_out <- data.table(id = x, source = "oncokb_cna")

            if (NROW(oncokb) > 0) {
                oncokb_cna$snpeff_ontology <- snpeff_ontology$short[match(oncokb$Consequence, snpeff_ontology$eff)]
                # Classify variant based on levels of evidence
                # T1&2 for TX
                # T1 for DX & PX
                # TX/DX/PX Present = "Actionable"
                # Show drugs in separate column?
                # Oncogenic/Likely Oncogenic = "Relevant"
                # Rest = "VUS"
                # TODO: Tumor Type Specific annotations - filter with annotations.tsv from OncoKB
                oncokb_cna <- parse_oncokb_tier(
                    oncokb_cna,
                    tx_cols = c("LEVEL_1", "LEVEL_2"),
                    rx_cols = c("LEVEL_R1"),
                    dx_cols = c("LEVEL_Dx1"),
                    px_cols = c("LEVEL_Px1")
                )

                # scna[, .(id = x, value = min_cn, type, track, gene = gene_name)]
                concat_out <- oncokb_cna[, .(
                    id = x,
                    gene = Hugo_Symbol,
                    value = min_cn,
                    type = ifelse(ALTERATION == "Amplification", "amp", ifelse(ALTERATION == "Deletion", "homdel", NA_character_)),
                    tier = tier,
                    tier_description = tier_factor,
                    therapeutics = tx_string, # comes from parse_oncokb_tier
                    resistances = rx_string,
                    diagnoses = dx_string,
                    prognoses = px_string,
                    track = "scna"
                )]
            }
            out <- rbind(out, concat_out, fill = TRUE, use.names = TRUE)
        }


        ## collect signatures
        if (!is.null(dat$signature_counts) && file.exists(dat[x, signature_counts])) {
            if (verbose) {
                message("pulling $signature_counts for ", x)
            }
            sig <- fread(dat[x, signature_counts])
            sig <- sig[, .(id = x, value = num_events, type = Signature, etiology = Etiology, frac = frac.events, track = "signature", source = "signature_counts")]
            out <- rbind(out, sig, fill = TRUE, use.names = TRUE)
        } else {
            out <- rbind(out, data.table(id = x, type = NA, source = "signature_counts"), fill = TRUE, use.names = TRUE)
        }

        ## collect gene mutations
        keepeff <- c("trunc", "cnadel", "cnadup", "complexsv", "splice", "inframe_indel", "fusion", "missense", "promoter", "regulatory", "mir")

        if (!is.null(dat$annotated_bcf) && file.exists(dat[x, annotated_bcf])) {
            if (verbose) {
                message("pulling $annotated_bcf for ", x, " using FILTER=", filter)
            }
            local_bcftools_path <- Sys.which("bcftools")
            local_bcftools_path <- ifelse(local_bcftools_path == "", stop("bcftools not found in the system PATH. Please install or moudule load bcftools."), local_bcftools_path)
            message("bcftools found at: ", local_bcftools_path)
            bcf <- skitools::grok_bcf(dat[x, annotated_bcf], label = x, long = TRUE, filter = filter, bpath = local_bcftools_path)
            if (verbose) {
                message(length(bcf), " variants pass filter")
            }
            genome.size <- sum(seqlengths(bcf), na.rm = TRUE) / 1e6
            if (is.na(genome.size)) { ## something went wrong with vcf
                genome.size <- sum(seqlengths(gG(jabba = dat[x, jabba_rds])), na.rm = TRUE) / 1e6
            }
            nmut <- data.table(as.character(seqnames(bcf)), start(bcf), end(bcf), bcf$REF, bcf$ALT) %>%
                unique() %>%
                nrow()
            mut.density <- data.table(id = x, value = c(nmut, nmut / genome.size), type = c("count", "density"), track = "tmb", source = "annotated_bcf")
            out <- rbind(out, mut.density, fill = TRUE, use.names = TRUE)
            bcf <- bcf[bcf$short %in% keepeff]
            if (verbose) {
                message(length(bcf), " variants pass keepeff")
            }
            vars <- NULL
            if (length(bcf)) {
                bcf$variant.g <- paste0(seqnames(bcf), ":", start(bcf), "-", end(bcf), " ", bcf$REF, ">", bcf$ALT)
                vars <- gr2dt(bcf)[, .(id = x, gene, vartype, variant.g, variant.p, distance, annotation, type = short, track = "variants", source = "annotated_bcf")]
                setkey(vars, variant.g)
                vars <- vars[, .SD[1], by = variant.g]
            }
            out <- rbind(out, vars, fill = TRUE, use.names = TRUE)
        } else {
            out <- rbind(out, data.table(id = x, type = NA, source = "annotated_bcf"), fill = TRUE, use.names = TRUE)
        }

        if (!is.null(dat$oncokb_maf) && file.exists(dat[x, oncokb_maf])) {
            oncokb <- data.table::fread(dat[x, oncokb_maf])
            concat_out <- data.table(id = x, source = "oncokb_maf")

            if (NROW(oncokb) > 0) {
                oncokb$snpeff_ontology <- snpeff_ontology$short[match(oncokb$Consequence, snpeff_ontology$eff)]
                # Classify variant based on levels of evidence
                # T1&2 for TX
                # T1 for DX & PX
                # TX/DX/PX Present = "Actionable"
                # Show drugs in separate column?
                # Oncogenic/Likely Oncogenic = "Relevant"
                # Rest = "VUS"
                # TODO: Tumor Type Specific annotations - filter with annotations.tsv from OncoKB
                oncokb <- oncokb[oncokb$snpeff_ontology %in% keepeff]
                oncokb <- parse_oncokb_tier(
                    oncokb,
                    tx_cols = c("LEVEL_1", "LEVEL_2"),
                    rx_cols = c("LEVEL_R1"),
                    dx_cols = c("LEVEL_Dx1"),
                    px_cols = c("LEVEL_Px1")
                )
                concat_out <- oncokb[, .(
                    id = x,
                    gene = Hugo_Symbol,
                    variant.g = paste("g.", Start_Position, "-", End_Position, sep = ""),
                    variant.c = HGVSc,
                    variant.p = HGVSp,
                    annotation = Consequence,
                    type = snpeff_ontology,
                    tier = tier,
                    tier_description = tier_factor,
                    therapeutics = tx_string, # comes from parse_oncokb_tier
                    resistances = rx_string,
                    diagnoses = dx_string,
                    prognoses = px_string,
                    distance = NA_integer_,
                    major.count,
                    minor.count,
                    major_snv_copies,
                    minor_snv_copies,
                    total_copies,
                    VAF,
                    track = "variants"
                )]
            }
            out <- rbind(out, concat_out, fill = TRUE, use.names = TRUE)
        }

        if (verbose) {
            message("done ", x)
        }

        return(out)
    }

    if (is.null(key(tumors))) {
        if (is.null(tumors$id)) {
            stop("Input tumors table must be keyed or have column $id")
        } else {
            setkey(tumors, id)
        }
    }

    out <- mclapply(tumors[[key(tumors)]], .oncotable,
        dat = tumors, pge = pge, amp.thresh = amp.thresh, filter = filter, del.thresh = del.thresh, verbose = verbose, mc.cores = mc.cores
    )
    out <- rbindlist(out, fill = TRUE, use.names = TRUE)

    setnames(out, "id", key(tumors))
    return(out)
}


#' @name create_oncotable
#' @title create_oncotable
#' @description
#'
#' function to create oncotable for use with filtered_events_json
#'
#' @param pair patient id
#' @param annotated_bcf bcf with variant annotations (snpeff)
#' @param signature_counts output of Signatures task (optional)
#' @param jabba_simple jabba.simple.rds from JaBbA
#' @param fusions output of fusions task
#' @param events output of events task
#' @param gencode file to gencode annotations (uses v29lift37 by default)
#' @param amp_thresh_multiplier amp.thresh for oncotable is amp_thresh_multiplier*ploidy
#' @param outdir path to directory in which to write oncotable outputs
#' @return data.table or NULL
#' @export
#' @author Shihab Dider, Joel Rosiene
create_oncotable <- function(
    pair,
    annotated_bcf,
    signature_counts = NULL,
    fusions,
    jabba_simple,
    karyograph,
    events,
    gencode = NULL,
    amp_thresh_multiplier = NULL,
    outdir) {
    tumors <- data.table(
        id = pair,
        annotated_bcf = annotated_bcf,
        signature_counts = signature_counts,
        fusions = fusions,
        jabba_rds = jabba_simple,
        karyograph = karyograph,
        complex = events
    )
    setkey(tumors, id)

    if (system("which bcftools", intern = TRUE) == "") {
        stop("bcftools is not available on the system PATH. Try `module load htslib` first or install it.")
    } else {
        message("bcftools is available.")
    }

    if (is.null(gencode)) {
        warning("path to gencode was not passed, setting to default ~/DB/GENCODE/gencode.v29lift37.annotation.nochr.rds")
        gencode <- "~/DB/GENCODE/gencode.v29lift37.annotation.nochr.rds"
    }

    if (is.null(amp_thresh_multiplier)) {
        warning("amp_thres_multiplier was not passed, setting to default 1.5")
        amp_thresh_multiplier <- 1.5
    }

    ploidy <- ifelse(is.null(readRDS(jabba_simple)$ploidy), readRDS(jabba_simple)$meta$ploidy, readRDS(jabba_simple)$ploidy)

    amp_thresh <- amp_thresh_multiplier * ploidy
    message(paste("using amp.thresh of", amp_thresh))

    oncotable <- oncotable(
        tumors,
        gencode = gencode,
        filter = "PASS",
        verbose = TRUE,
        amp.thresh = amp_thresh
    )

    if (!dir.exists(outdir)) {
        dir.create(outdir, recursive = TRUE)
    }

    saveRDS(oncotable, paste0(outdir, "/", "oncotable.rds"))
    fwrite(oncotable, paste0(outdir, "/", "oncotable.txt"))
}

#' @name filtered_events_json
#' @title filtered_events_json
#' @description
#'
#' function to create filtered events json for case reports
#'
#' @param pair patient id to be added to pgvdb or case reports
#' @param oncotable oncotable task output
#' @param jabba_gg JaBbA output ggraph or complex
#' @param out_file path to write json
#' @param cgc_file path to cgc file to annotate drivers
#' @param oncokab_file path to oncokb gene tiers to annotate drivers
#' @param return_table TRUE/FALSE whether to return the data.table that is used for creating the json
#' @return data.table or NULL
#' @export
#' @author Stanley Clarke, Tanubrata Dey

filtered_events_json <- function(
    pair,
    oncotable,
    jabba_gg,
    out_file,
    cgc_file = "/gpfs/commons/groups/imielinski_lab/DB/COSMIC/v99_GRCh37/cancer_gene_census_fixed.csv",
    oncokb_file = NULL,
    temp_fix = FALSE,
    return_table = FALSE) {
    ## Driver CNA windows
    ## Load details from oncotable
    ot <- readRDS(oncotable)
    # add a fusion_gene_coords column of NAs if no fusions
    if (!"fusion_gene_coords" %in% colnames(ot)) {
        ot[, fusion_genes := NA]
        ot[, fusion_gene_coords := NA]
    }
    snvs <- ot[grepl("frameshift|missense|stop|disruptive", annotation)]
    snvs <- snvs[!duplicated(variant.p)]
    ## Note here probably have to crossreference these missense muts with hetdels
    # hetdel_snvs = snvs[gene %in% ot[type == 'hetdel',gene]]
    # possible_drivers = rbind(hetdel_snvs,homdels)
    homdels <- ot[type == "homdel"]
    amps <- ot[type == "amp"]
    fusions <- ot[type == "fusion"]
    jab <- readRDS(jabba_gg)
    possible_drivers <- rbind(snvs, homdels, amps, fusions)
    cgc <- fread(cgc_file)
    names(cgc) <- gsub(" ", ".", names(cgc))
    cgc$gene <- cgc$Gene.Symbol
    longlist <- merge.data.table(possible_drivers, cgc, by = "gene")
    res <- longlist[, .(gene, fusion_genes, id, vartype, type, variant.g, variant.p, Name, Genome.Location, fusion_gene_coords, Tier, Role.in.Cancer)]
    names(res) <- c("gene", "fusion_genes", "id", "vartype", "type", "Variant_g", "Variant", "Name", "Genome_Location", "fusion_gene_coords", "Tier", "Role_in_Cancer")
    if (!is.null(oncokb_file)) {
        oncokb <- fread(oncokb_file)
        names(oncokb) <- gsub(" ", "_", names(oncokb))
        oncokb$gene <- oncokb$Gene
        oncokb$drugs <- oncokb[, "Drugs_(for_therapeutic_implications_only)"]
        longlist <- merge(longlist, oncokb, by = "gene", all.x = TRUE)
        res <- longlist[, .(gene, fusion_genes, id, vartype, type, variant.g, variant.p, Name, Genome.Location, fusion_gene_coords, Tier, Level, Alterations, Cancer_Types, drugs, Role.in.Cancer)]
        names(res) <- c("gene", "fusion_genes", "id", "vartype", "type", "Variant_g", "Variant", "Name", "Genome_Location", "fusion_gene_coords", "Tier", "Oncokb_Tier", "Oncokb_Alterations", "Oncokb_Cancer_Types", "Oncokb_Drugs", "Role_in_Cancer")
    }
    res <- res %>% unique(., by = c("gene", "Variant"))
    if (nrow(res) > 0) {
        res[, seqnames := tstrsplit(Genome_Location, ":", fixed = TRUE, keep = 1)]
        res[, start := tstrsplit(Genome_Location, "-", fixed = TRUE, keep = 1)]
        res[, start := tstrsplit(start, ":", fixed = TRUE, keep = 2)]
        res[, end := tstrsplit(Genome_Location, "-", fixed = TRUE, keep = 2)]
        res.mut <- res[!is.na(Variant), ]
        if (nrow(res.mut) > 0) {
            res.mut[, Variant := gsub("p.", "", Variant)]
        }
        res.cn <- res[is.na(Variant), ]
        if (nrow(res.cn) > 0) {
            res.cn.gr <- GRanges(res.cn)
            res.cn.gr <- gr.val(res.cn.gr, jab$nodes$gr, c("cn", "cn.low", "cn.high"))
            res.cn.dt <- as.data.table(res.cn.gr)
            res.cn.dt[!is.na(cn) & !is.na(cn.low) & !is.na(cn.high), Variant := paste0("Total CN:", round(cn, digits = 3), "; CN Minor:", round(cn.low, digits = 3), "; CN Major:", round(cn.high, digits = 3))]
            res.cn.dt[!is.na(cn) & is.na(cn.low) & is.na(cn.high), Variant := paste0("Total CN:", round(cn, digits = 3))]
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
        message(paste0("Writing json to ", out_file))
        write_json(res.final, out_file, pretty = TRUE)
        res.final[, sample := pair]
        if (return_table) {
            return(res.final)
        }
    }
}

#' @name strelka2counts
#' @title strelka2counts
#' @description
#' takes in a strelka vcf and returns a data.table with total count and count of variants with normal vaf greater than 0
#'
#' @param vcf patient id to be added to pgvdb or case reports
#' @param seqnames_genome_width chromosomes to count variants in
#' @param type_return counts or dt, if counts returns counts for snv and snv count with a normal vaf greater than 0. If dt, returns the parsed vcf
#' @return data.table
#' @export
#' @author Stanley Clarke, Tanubrata Dey
strelka2counts <- function(vcf, seqnames_genome_width = c(1:22, "X", "Y"), type_return = "counts") {
    somatic.filtered.vcf <- read.delim(vcf, header = F, comment.char = "#", col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR")) %>% as.data.table()
    ## strip chr
    somatic.filtered.vcf[, CHROM := gsub("chr", "", CHROM)]
    ##
    sub.vcf <- somatic.filtered.vcf[CHROM %in% seqnames_genome_width, ]
    if (grepl("INDEL", sub.vcf$FORMAT)) {
        sub.vcf[, c("DP", "FDP", "SDP", "SUBDP", "AU", "CU", "GU", "TU", "INDEL") := tstrsplit(NORMAL, ":", fixed = TRUE)]
    } else if (grepl("DP2", sub.vcf$FORMAT)) {
        sub.vcf[, c("DP", "DP2", "TAR", "TIR", "TOR", "DP50", "FDP50", "SUBDP50", "BCN50") := tstrsplit(NORMAL, ":", fixed = TRUE)]
        sub.vcf[, ref_count := tstrsplit(TAR, ",", fixed = TRUE, keep = 1)]
        sub.vcf[, alt_count := tstrsplit(TIR, ",", fixed = TRUE, keep = 1)]
    } else {
        sub.vcf[, c("DP", "FDP", "SDP", "SUBDP", "AU", "CU", "GU", "TU") := tstrsplit(NORMAL, ":", fixed = TRUE)]

        # Define a function to perform conditional assignment
        conditional_assign <- function(ref, alt, ref_col, alt_col, ref_data, alt_data) {
            condition <- sub.vcf[REF == ref & ALT == alt]
            if (nrow(condition) > 0) {
                sub.vcf[REF == ref & ALT == alt, c(ref_col, alt_col) := .(tstrsplit(get(ref_data), ",", fixed = TRUE)[[1]], tstrsplit(get(alt_data), ",", fixed = TRUE)[[1]])]
            }
        }

        # Define the conditions and corresponding columns
        conditions <- list(
            list(ref = "A", alt = "T", ref_data = "AU", alt_data = "TU"),
            list(ref = "A", alt = "C", ref_data = "AU", alt_data = "CU"),
            list(ref = "A", alt = "G", ref_data = "AU", alt_data = "GU"),
            list(ref = "C", alt = "T", ref_data = "CU", alt_data = "TU"),
            list(ref = "C", alt = "A", ref_data = "CU", alt_data = "AU"),
            list(ref = "C", alt = "G", ref_data = "CU", alt_data = "GU"),
            list(ref = "G", alt = "T", ref_data = "GU", alt_data = "TU"),
            list(ref = "G", alt = "C", ref_data = "GU", alt_data = "CU"),
            list(ref = "G", alt = "A", ref_data = "GU", alt_data = "AU"),
            list(ref = "T", alt = "A", ref_data = "TU", alt_data = "AU"),
            list(ref = "T", alt = "C", ref_data = "TU", alt_data = "CU"),
            list(ref = "T", alt = "G", ref_data = "TU", alt_data = "GU")
        )

        # Apply the function to each condition for normal counts
        lapply(conditions, function(cond) {
            conditional_assign(cond$ref, cond$alt, "ref_count", "alt_count", cond$ref_data, cond$alt_data)
        })
    }
    # Continue with the rest of your code
    sub.vcf[, ref_count := as.numeric(ref_count)]
    sub.vcf[, alt_count := as.numeric(alt_count)]
    sub.vcf[, normal_vaf := alt_count / (ref_count + alt_count)]

    if (grepl("INDEL", sub.vcf$FORMAT)) {
        sub.vcf[, c("DP", "FDP", "SDP", "SUBDP", "AU", "CU", "GU", "TU", "INDEL") := tstrsplit(TUMOR, ":", fixed = TRUE)]
    } else if (grepl("DP2", sub.vcf$FORMAT)) {
        sub.vcf[, c("DP", "DP2", "TAR", "TIR", "TOR", "DP50", "FDP50", "SUBDP50", "BCN50") := tstrsplit(TUMOR, ":", fixed = TRUE)]
        sub.vcf[, ref_count_tumor := tstrsplit(TAR, ",", fixed = TRUE, keep = 1)]
        sub.vcf[, alt_count_tumor := tstrsplit(TIR, ",", fixed = TRUE, keep = 1)]
    } else {
        sub.vcf[, c("DP", "FDP", "SDP", "SUBDP", "AU", "CU", "GU", "TU") := tstrsplit(TUMOR, ":", fixed = TRUE)]
        # Apply the function to each condition for tumor counts
        lapply(conditions, function(cond) {
            conditional_assign(cond$ref, cond$alt, "ref_count_tumor", "alt_count_tumor", cond$ref_data, cond$alt_data)
        })
    }
    sub.vcf[, ref_count_tumor := as.numeric(ref_count_tumor)]
    sub.vcf[, alt_count_tumor := as.numeric(alt_count_tumor)]
    sub.vcf[, tumor_vaf := alt_count_tumor / (ref_count_tumor + alt_count_tumor)]

    if (type_return == "counts") {
        snv_count <- nrow(sub.vcf)
        snv_count_normal_vaf_greater0 <- nrow(sub.vcf[normal_vaf > 0, ])
        return(data.table(
            category = c("snv_count", "snv_count_normal_vaf_greater0"),
            counts = c(snv_count, snv_count_normal_vaf_greater0)
        ))
    } else if (type_return == "dt") {
        return(sub.vcf)
    }
}


#' @name parse_vcf_strelka2
#' @title parse_vcf_strelka2
#' @description
#' takes in a strelka vcf and returns a data.table format of the vcf file with required fields for strelka_qc function
#'
#' @param vcf patient id to be added to pgvdb or case reports
#' @param seqnames_genome_width chromosomes to count variants in
#' @return data.table
#' @export
#' @author Tanubrata Dey, Stanley Clarke
parse_vcf_strelka2 <- function(vcf, seqnames_genome_width = c(1:22, "X", "Y")) {
    somatic.filtered.vcf <- read.delim(vcf, header = F, comment.char = "#", col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR")) %>% as.data.table()
    ## strip chr
    somatic.filtered.vcf[, CHROM := gsub("chr", "", CHROM)]
    ## normal counts and VAF
    sub.vcf <- somatic.filtered.vcf[CHROM %in% seqnames_genome_width, ]
    sub.vcf[, somatic_EVS := as.numeric(gsub(".*SomaticEVS=([^;]+).*", "\\1", INFO))]
    sub.vcf[, MQ := as.numeric(gsub(".*MQ=([^;]+).*", "\\1", INFO))]
    sub.vcf[, c("N_DP", "N_FDP", "N_SDP", "N_SUBDP", "N_AU", "N_CU", "N_GU", "N_TU", "N_INDEL") := tstrsplit(NORMAL, ":", fixed = TRUE)]
    sub.vcf[REF == "A" & ALT == "T", c("ref_count_N", "alt_count_N") := c(tstrsplit(N_AU, ",", fixed = TRUE, keep = 1), tstrsplit(N_TU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "A" & ALT == "C", c("ref_count_N", "alt_count_N") := c(tstrsplit(N_AU, ",", fixed = TRUE, keep = 1), tstrsplit(N_CU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "A" & ALT == "G", c("ref_count_N", "alt_count_N") := c(tstrsplit(N_AU, ",", fixed = TRUE, keep = 1), tstrsplit(N_GU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "C" & ALT == "T", c("ref_count_N", "alt_count_N") := c(tstrsplit(N_CU, ",", fixed = TRUE, keep = 1), tstrsplit(N_TU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "C" & ALT == "A", c("ref_count_N", "alt_count_N") := c(tstrsplit(N_CU, ",", fixed = TRUE, keep = 1), tstrsplit(N_AU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "C" & ALT == "G", c("ref_count_N", "alt_count_N") := c(tstrsplit(N_CU, ",", fixed = TRUE, keep = 1), tstrsplit(N_GU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "G" & ALT == "T", c("ref_count_N", "alt_count_N") := c(tstrsplit(N_GU, ",", fixed = TRUE, keep = 1), tstrsplit(N_TU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "G" & ALT == "C", c("ref_count_N", "alt_count_N") := c(tstrsplit(N_GU, ",", fixed = TRUE, keep = 1), tstrsplit(N_CU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "G" & ALT == "A", c("ref_count_N", "alt_count_N") := c(tstrsplit(N_GU, ",", fixed = TRUE, keep = 1), tstrsplit(N_AU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "T" & ALT == "A", c("ref_count_N", "alt_count_N") := c(tstrsplit(N_TU, ",", fixed = TRUE, keep = 1), tstrsplit(N_AU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "T" & ALT == "C", c("ref_count_N", "alt_count_N") := c(tstrsplit(N_TU, ",", fixed = TRUE, keep = 1), tstrsplit(N_CU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "T" & ALT == "G", c("ref_count_N", "alt_count_N") := c(tstrsplit(N_TU, ",", fixed = TRUE, keep = 1), tstrsplit(N_GU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[, c("N_DP", "N_DP2", "N_TAR", "N_TIR", "N_TOR", "N_DP50", "N_FDP50", "N_SUBDP50", "N_BCN50") := tstrsplit(NORMAL, ":", fixed = TRUE)]
    sub.vcf[is.na(ref_count_N), ref_count_N := tstrsplit(N_TAR, ",", fixed = TRUE, keep = 1)]
    sub.vcf[is.na(alt_count_N), alt_count_N := tstrsplit(N_TIR, ",", fixed = TRUE, keep = 1)]
    sub.vcf[, ref_count_N := as.numeric(ref_count_N)]
    sub.vcf[, alt_count_N := as.numeric(alt_count_N)]
    sub.vcf[, VAF_N := alt_count_N / (ref_count_N + alt_count_N)]
    sub.vcf[, c("N_FDP", "N_SDP", "N_SUBDP", "N_AU", "N_CU", "N_GU", "N_TU", "N_INDEL", "N_DP2", "N_TAR", "N_TIR", "N_TOR", "N_DP50", "N_FDP50", "N_SUBDP50", "N_BCN50")] <- NULL
    ## tumor counts and VAF
    sub.vcf[, c("T_DP", "T_FDP", "T_SDP", "T_SUBDP", "T_AU", "T_CU", "T_GU", "T_TU", "T_INDEL") := tstrsplit(TUMOR, ":", fixed = TRUE)]
    sub.vcf[REF == "A" & ALT == "T", c("ref_count_T", "alt_count_T") := c(tstrsplit(T_AU, ",", fixed = TRUE, keep = 1), tstrsplit(T_TU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "A" & ALT == "C", c("ref_count_T", "alt_count_T") := c(tstrsplit(T_AU, ",", fixed = TRUE, keep = 1), tstrsplit(T_CU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "A" & ALT == "G", c("ref_count_T", "alt_count_T") := c(tstrsplit(T_AU, ",", fixed = TRUE, keep = 1), tstrsplit(T_GU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "C" & ALT == "T", c("ref_count_T", "alt_count_T") := c(tstrsplit(T_CU, ",", fixed = TRUE, keep = 1), tstrsplit(T_TU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "C" & ALT == "A", c("ref_count_T", "alt_count_T") := c(tstrsplit(T_CU, ",", fixed = TRUE, keep = 1), tstrsplit(T_AU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "C" & ALT == "G", c("ref_count_T", "alt_count_T") := c(tstrsplit(T_CU, ",", fixed = TRUE, keep = 1), tstrsplit(T_GU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "G" & ALT == "T", c("ref_count_T", "alt_count_T") := c(tstrsplit(T_GU, ",", fixed = TRUE, keep = 1), tstrsplit(T_TU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "G" & ALT == "C", c("ref_count_T", "alt_count_T") := c(tstrsplit(T_GU, ",", fixed = TRUE, keep = 1), tstrsplit(T_CU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "G" & ALT == "A", c("ref_count_T", "alt_count_T") := c(tstrsplit(T_GU, ",", fixed = TRUE, keep = 1), tstrsplit(T_AU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "T" & ALT == "A", c("ref_count_T", "alt_count_T") := c(tstrsplit(T_TU, ",", fixed = TRUE, keep = 1), tstrsplit(T_AU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "T" & ALT == "C", c("ref_count_T", "alt_count_T") := c(tstrsplit(T_TU, ",", fixed = TRUE, keep = 1), tstrsplit(T_CU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[REF == "T" & ALT == "G", c("ref_count_T", "alt_count_T") := c(tstrsplit(T_TU, ",", fixed = TRUE, keep = 1), tstrsplit(T_GU, ",", fixed = TRUE, keep = 1)), ]
    sub.vcf[, c("T_DP", "T_DP2", "T_TAR", "T_TIR", "T_TOR", "T_DP50", "T_FDP50", "T_SUBDP50", "T_BCN50") := tstrsplit(TUMOR, ":", fixed = TRUE)]
    sub.vcf[is.na(ref_count_T), ref_count_T := tstrsplit(T_TAR, ",", fixed = TRUE, keep = 1)]
    sub.vcf[is.na(alt_count_T), alt_count_T := tstrsplit(T_TIR, ",", fixed = TRUE, keep = 1)]
    sub.vcf[, ref_count_T := as.numeric(ref_count_T)]
    sub.vcf[, alt_count_T := as.numeric(alt_count_T)]
    sub.vcf[, VAF_T := alt_count_T / (ref_count_T + alt_count_T)]
    sub.vcf[, c("T_FDP", "T_SDP", "T_SUBDP", "T_AU", "T_CU", "T_GU", "T_TU", "T_INDEL", "T_DP2", "T_TAR", "T_TIR", "T_TOR", "T_DP50", "T_FDP50", "T_SUBDP50", "T_BCN50")] <- NULL
    return(sub.vcf)
}


#' @name strelka_qc
#' @title strelka_qc
#' @description
#'
#' function to create json for strelka qc plotting in case reports
#'
#' @param path to strelka2 vcf file to be used
#' @param chromosomes to select for. default: c(1:22,"X","Y")
#' @param outfile path to write json
#' @param write_json TRUE/FALSE whether to write the json
#' @param return_table TRUE/FALSE whether to return the data
#' @return NULL
#' @export
#' @author Tanubrata Dey
strelka_qc <- function(vcf, seqnames_genome_width = c(1:22, "X", "Y"), outfile, write_json = TRUE, return_table = TRUE) {
    sq <- parse_vcf_strelka2(vcf, seqnames_genome_width = seqnames_genome_width) %>% dplyr::select(CHROM, POS, REF, ALT, FILTER, T_DP, N_DP, alt_count_N, alt_count_T, MQ, VAF_T, somatic_EVS)
    names(sq) <- c("chromosome", "position", "reference", "alternate", "filter", "tumor_depth", "normal_depth", "normal_alt_counts", "tumor_alt_counts", "mapping_quality", "tumor_VAF", "somatic_EVS")
    consider_numeric <- c("tumor_depth", "normal_depth", "normal_alt_counts", "tumor_alt_counts", "mapping_quality", "tumor_vaf", "somatic_EVS")
    sq[, (consider_numeric) := lapply(.SD, as.numeric), .SDcols = consider_numeric]
    sq[, id := .I]
    if (write_json) {
        ## write the json
        message(paste0("Writing json to ", outfile))
        write_json(sq, outfile, pretty = TRUE)
        if (return_table) {
            return(sq)
        }
    } else {
        return(sq)
    }
}

#' @name sage_count
#' @title sage_count
#' @description
#' takes in a SAGE vcf and returns a data.table with total count and count of variants with normal vaf greater than 0
#'
#' @param vcf_path Path to the SAGE VCF file
#' @param genome Reference genome name (e.g., "hg19", "hg38")
#' @return data.table
#' @export
#' @author Shihab Dider, Tanubrata Dey
sage_count <- function(
    vcf_path,
    genome) {
    if (!file.exists(vcf_path)) {
        stop("VCF file does not exist.")
    }
    vcf <- readVcf(vcf_path, genome)

    # Filter for PASS variants
    pass_variants <- rowRanges(vcf)$FILTER == "PASS"
    vcf <- vcf[pass_variants, ]

    snv_count <- length(vcf)

    # Extract VAF_N from the genotype (geno) slot
    geno_data <- geno(vcf)
    normal <- colnames(geno_data$DP)[1]
    tumor <- colnames(geno_data$DP)[2]

    if (normal %in% colnames(geno_data$DP)) {
        if (!is.null(geno_data$AF)) {
            VAF_N <- as.numeric(geno_data$AF[, normal])
        } else if (!is.null(geno_data$TAR)) {
            alt_count_N <- as.numeric(geno_data$TAR[, normal, 1])
            total_read_depth_N <- as.numeric(geno_data$DP[, normal])
            VAF_N <- alt_count_N / total_read_depth_N
        } else if (!is.null(geno_data$AU)) {
            vcf_dt <- as.data.table(rowRanges(vcf))

            # ALT is stored as DNAStringSet, so we need to extract as a character
            get_alt_allele <- function(dna_string_set) {
                as.character(dna_string_set[[1]])
            }

            alt_alleles <- sapply(vcf_dt$ALT, get_alt_allele)

            # Calculate the number of alternate alleles for each variant
            alt_allele_counts <- geno_data$AU[, normal, 1] * (alt_alleles == "A") +
                geno_data$CU[, normal, 1] * (alt_alleles == "C") +
                geno_data$GU[, normal, 1] * (alt_alleles == "G") +
                geno_data$TU[, normal, 1] * (alt_alleles == "T")

            total_alleles <- geno_data$DP[, normal]

            vcf_dt[, normal_VAF := round(alt_allele_counts / total_alleles, 2)]
            VAF_N <- vcf_dt$normal_VAF
        } else {
            stop("No AF or TAR present in the VCF file.")
        }
        snv_count_normal_vaf_greater0 <- sum(VAF_N > 0, na.rm = TRUE)
        return(data.table(
            category = c("snv_count", "snv_count_normal_vaf_greater0"),
            counts = c(snv_count, snv_count_normal_vaf_greater0)
        ))
    } else {
        print("Tumor only run VCF provided, no VAF_N present.")
        return(data.table(
            category = c("snv_count", "snv_count_normal_vaf_greater0"),
            counts = c(snv_count, NA)
        ))
    }
}


#' @name sage_qc
#' @title sage_qc
#' @description
#'
#' function to create json for sage qc plotting in case reports
#'
#' @param vcf_path Path to sage VCF file to be used (paired T-N/ Tumor only)
#' @param genome Reference genome name (e.g., "hg19", "hg38")
#' @param outfile Path to write JSON
#' @param write_json TRUE/FALSE whether to write the JSON
#' @param return_table TRUE/FALSE whether to return the data
#' @return NULL
#' @export
#' @author Shihab Dider, Tanubrata Dey, Stanley Clarke
sage_qc <- function(
    vcf_path,
    genome,
    outfile,
    write_json = TRUE,
    return_table = TRUE) {
    vcf <- readVcf(vcf_path, genome)

    # Filter for PASS variants
    pass_variants <- rowRanges(vcf)$FILTER == "PASS"
    vcf <- vcf[pass_variants, ]

    # Extract necessary information from VCF object
    chrom <- as.character(seqnames(rowRanges(vcf)))
    pos <- start(rowRanges(vcf))
    ref <- as.character(ref(vcf))
    alt <- as.character(unlist(alt(vcf)))
    filter <- as.character(rowRanges(vcf)$FILTER)
    qual <- as.numeric(rowRanges(vcf)$QUAL)

    # Extract depth and allele count information from the genotype (geno) slot
    geno_data <- geno(vcf)
    normal <- colnames(geno_data$DP)[1]
    tumor <- colnames(geno_data$DP)[2]

    T_DP <- as.numeric(geno_data$DP[, tumor])
    alt_count_T <- sapply(geno_data$AD[, tumor], function(x) as.numeric(x[2])) # Extract the second element for alternate allele depth
    T_ABQ <- as.numeric(geno_data$ABQ[, tumor])
    VAF_T <- as.numeric(geno_data$AF[, tumor])

    # Check if normal sample data exists
    if (normal %in% colnames(geno_data$DP)) {
        N_DP <- as.numeric(geno_data$DP[, normal])
        alt_count_N <- sapply(geno_data$AD[, normal], function(x) as.numeric(x[2])) # Extract the second element for alternate allele depth
        VAF_N <- as.numeric(geno_data$AF[, normal])

        sq <- data.table(
            chromosome = chrom,
            position = pos,
            reference = ref,
            alternate = alt,
            filter = filter,
            mapping_quality = qual,
            tumor_depth = T_DP,
            normal_depth = N_DP,
            normal_alt_counts = alt_count_N,
            tumor_alt_counts = alt_count_T,
            tumor_abq = T_ABQ,
            tumor_vaf = VAF_T,
            normal_vaf = VAF_N
        )

        consider_numeric <- c("tumor_depth", "normal_depth", "normal_alt_counts", "tumor_alt_counts", "tumor_abq", "tumor_vaf", "normal_vaf")
    } else {
        sq <- data.table(
            chromosome = chrom,
            position = pos,
            reference = ref,
            alternate = alt,
            filter = filter,
            mapping_quality = qual,
            tumor_depth = T_DP,
            tumor_alt_counts = alt_count_T,
            tumor_abq = T_ABQ,
            tumor_vaf = VAF_T
        )

        consider_numeric <- c("tumor_depth", "tumor_alt_counts", "tumor_abq", "tumor_vaf")
    }

    sq[, (consider_numeric) := lapply(.SD, as.numeric), .SDcols = consider_numeric]
    sq[, id := .I]

    if (write_json) {
        message(paste0("Writing json to ", outfile))
        write_json(sq, outfile, pretty = TRUE)
        if (return_table) {
            return(sq)
        }
    } else {
        return(sq)
    }
}

#' @name dlrs
#' @title derivative log ratio spread
#' @description
#'
#' function to get the dlrs for coverage data
#' used in meta_data_json
#'
#' @param x foreground from dryclean
#' @return dlrs
#' @export
#' @author Joel Rosiene
dlrs <- function(x) {
    nx <- length(x)
    if (nx < 3) {
        stop("Vector length>2 needed for computation")
    }
    tmp <- embed(x, 2)
    diffs <- tmp[, 2] - tmp[, 1]
    dlrs <- IQR(diffs, na.rm = TRUE) / (sqrt(2) * 1.34)
    return(dlrs)
}


#' @name create_mutations_catalog_json
#' @title create_mutations_catalog_json
#' @description
#'
#' Builds the mutations catalog for case reports using the output of sigprofiler matrix generator
#'
#' @param sig_matrix_path Path to either the SBS or ID matrix
#' @param is_indel Boolean to indicate if the matrix is for indels
#' @param output_dir directory to write the jsons, will create directory for each patient id
#' @return
#' @export
#' @author Shihab Dider, Sukanya Panja
create_mutations_catalog_json <- function(
    sig_matrix_path,
    is_indel,
    output_dir) {
    sig_matrix <- fread(sig_matrix_path)
    sig_matrix_dt <- as.data.frame(sig_matrix)
    samples <- colnames(sig_matrix_dt)
    for (i in 2:(ncol(sig_matrix_dt))) { # skip first column as it has the tnc
        pair <- samples[i]

        if (is_indel) {
            samp_data <- data.frame(
                id = 1:nrow(sig_matrix_dt),
                insdel = sig_matrix_dt[, 1],
                mutations = sig_matrix_dt[, i]
            )
            catalog_file_name <- "id_mutation_catalog.json"
        } else {
            samp_data <- data.frame(
                id = 1:nrow(sig_matrix_dt),
                tnc = sig_matrix_dt[, 1],
                mutations = sig_matrix_dt[, i]
            )
            catalog_file_name <- "mutation_catalog.json"
        }
        pair_data <- list(pair = pair, data = samp_data)
        system(paste("mkdir -p", paste0(output_dir, "/", pair)))
        output_path <- paste0(output_dir, "/", catalog_file_name)
        write_json(
            pair_data,
            output_path,
            pretty = TRUE,
            auto_unbox = TRUE
        )
    }
}

#' @name meta_data_json
#' @title meta_data_json
#' @description
#'
#' creates the meta data summary json for case reports
#'
#' @param pair patient id to be added to pgvdb or case reports
#' @param outdir path to parent directory containing patient directories
#' @param coverage path dryclean coverage output
#' @param coverage_qc list of coverage qc metrics (tbd)
#' @param jabba_gg path to JaBbA output ggraph (non-integer-balanced)
#' @param complex path to events output ggraph
#' @param strelka2_vcf path to strelka vcf to get snv count (need either this or sage_vcf)
#' @param sage_vcf path to SAGE vcf to get snv count (sage preferred over strelka2)
#' @param tumor_type_final tumor type abbreviation of the sample
#' @param disease full length tumor type
#' @param primary_site primary site of tumor
#' @param inferred_sex sex of the patient
#' # @param karyograph JaBbA outputted karygraph
#' @param signatures_pair_name pair name that appears in the sigprofiler output (sometimes different from pair)
#' @param indel_sigprofiler path to Assignment_Solution_Activities.txt, output of sigprofiler
#' @param sbs_sigprofiler path to Assignment_Solution_Activities.txt, output of sigprofiler
#' @param sbs_deconstructSigs path to rds output from deconstructSigs
#' @param seqnames_loh chromosomes to be used to calculate LOH
#' @param seqnames_genome_width chromosomes to be used to calculate tmb
#' @param hrdetect path to HRDetect output
#' @param onenesstwoness path to onenesstwoness output
#' @param write_json TRUE/FALSE to write the json
#' @param overwrite TRUE/FALSE to overwrite the present json
#' @param return_table TRUE/FALSE to return the data.table output
#' @param make_dir TRUE/FALSE make the directory for the patient sample if it does not exists
#' @return data.table or NULL
#' @export
#' @author Stanley Clarke, Tanubrata Dey, Joel Rosiene, Shihab Dider

meta_data_json <- function(
    pair,
    outdir = "./",
    genome = "hg19",
    coverage = NULL,
    coverage_qc = NULL,
    het_pileups_wgs = NULL,
    jabba_gg = NULL,
    complex = NULL,
    strelka2_vcf = NULL,
    sage_vcf = NULL,
    tumor_type = NULL,
    disease = NULL,
    primary_site = NULL,
    inferred_sex = NULL,
    signatures_pair_name = NULL, # if different from pair
    indel_sigprofiler = NULL,
    sbs_sigprofiler = NULL,
    sbs_deconstructSigs = NULL,
    seqnames_loh = c(1:22),
    seqnames_genome_width = c(1:22, "X", "Y"),
    hrdetect = NULL,
    onenesstwoness = NULL,
    write_json = TRUE,
    overwrite = FALSE,
    return_table = FALSE,
    make_dir = FALSE) {
    out_file <- paste0(outdir, "/", pair, "/metadata.json")
    if (!overwrite && write_json == TRUE) {
        if (file.exists(out_file)) {
            print(paste0("Output already exists! - skipping sample ", pair))
            return(NA)
        }
    }

    ## get patient directory
    split_file_path <- strsplit(out_file, "/")[[1]]
    folder_path <- paste0(split_file_path[1:(length(split_file_path) - 1)], collapse = "/")
    ## check if directory exists
    if (!file.exists(folder_path)) {
        if (make_dir) {
            cmd <- paste0("mkdir -p ", folder_path)
            print(paste0("Making directory ", folder_path))
            system(cmd)
        } else {
            print(paste0("Folder does not exist; skipping sample ", pair, ". Use make_dir = TRUE to make directory"))
            return(NA)
        }
    }

    meta.dt <- data.table(pair = pair)

    if (!is.null(tumor_type)) {
        meta.dt[, tumor_type := tumor_type]
    }
    if (!is.null(disease)) {
        meta.dt[, disease := disease]
    }
    if (!is.null(primary_site)) {
        meta.dt[, primary_site := primary_site]
    }
    if (!is.null(inferred_sex)) {
        # check if inferred sex is one of "male" or "female"
        if (!inferred_sex == "male" & !inferred_sex == "female") {
            warning("inferred sex must be one of `male` or `female`, you passed: ")
            print(inferred_sex)
        }
        meta.dt[, inferred_sex := inferred_sex]
    } else if (!is.null(jabba_gg)) {
        gg <- readRDS(jabba_gg)
        ncn.x <- gg$nodes$dt[
            (seqnames == "X" | seqnames == "chrX"),
            weighted.mean(cn,
                w = end - start + 1,
                na.rm = TRUE
            )
        ]
        sex <- ifelse(ncn.x < 1.4, "male", "female")
        meta.dt[, inferred_sex := sex]
    } else if (!is.null(coverage)) {
        #' look at mean-normalized relative foreground
        ncn.x <- unique(gr2dt(coverage)[, foreground.chr := mean(foreground), by = seqnames][, .(seqnames, foreground.chr)])[seqnames %in% c("X", "chrX")]$foreground.chr
        sex <- ifelse(ncn.x < 0.7, "male", "female")
        meta.dt[, inferred_sex := sex]
    } else {
        warning("Could not infer sex. Pass it directly as `inferred_sex=[male|female]` or pass jabba_gg or coverage data to infer it from the copy-number/coverage.")
    }

    ## get derivative log ratio spread
    coverage_variance <- NULL
    if (!is.null(coverage)) {
        cov_var <- dlrs(readRDS(coverage)$foreground) / 100
        coverage_variance <- list(coverage_variance = cov_var)
    }

    ## coverage qc
    if (!is.null(coverage_qc)) {
        coverage_qc <- fread(coverage_qc)

        # Define a function to replace symbols with string equivalents
        replace_symbols <- function(x) {
            x <- gsub(" ", "_", x)
            x <- gsub("%", "percent", x)
            x <- gsub("\u2265", "greater_than_or_equal_to", x)
            # Add more replacements as needed
            return(x)
        }

        # Apply the function to column names
        setnames(coverage_qc, tolower(replace_symbols(names(coverage_qc))))

        if ("tumor_median_coverage" %in% names(coverage_qc)) {
            tumor_row <- coverage_qc[grepl(pair, sample_name)]
            meta.dt$tumor_median_coverage <- as.numeric(gsub("X", "", tumor_row$tumor_median_coverage))
        } else {
            tumor_row <- coverage_qc[grepl(pair, sample_name) & sample_class == "Tumor"]
            meta.dt$tumor_median_coverage <- as.numeric(gsub("X", "", tumor_row$median_cov))
        }

        if ("normal_median_coverage" %in% names(coverage_qc)) {
            normal_row <- coverage_qc[grepl(pair, sample_name)]
            meta.dt$normal_median_coverage <- as.numeric(gsub("X", "", normal_row$normal_median_coverage))
        } else {
            normal_row <- coverage_qc[grepl(pair, sample_name) & sample_class == "Normal"]
            meta.dt$normal_median_coverage <- as.numeric(gsub("X", "", normal_row$median_cov))
        }

        convert_to_numeric <- function(x) {
            if (is.character(x)) {
                is_percent <- grepl("%", x)
                x <- gsub("[^0-9.]", "", x)
                x <- as.numeric(x)
                if (is_percent) {
                    x <- x / 100
                }
            }
            return(x)
        }

        tumor_row[, (names(tumor_row)) := lapply(.SD, convert_to_numeric)]
        if ("patient.x" %in% names(tumor_row)) {
            coverage_qc <- subset(tumor_row, select = -c(median_cov, sample, sample_name, patient.x, patient.y, sample_class))
        } else {
            coverage_qc <- tumor_row
        }

        if (!is.null(coverage_variance)) {
            meta.dt$coverage_qc <- list(as.list(c(coverage_qc, coverage_variance)))
        } else {
            meta.dt$coverage_qc <- list(as.list(coverage_qc))
        }
    }

    ## snv counts
    if (!is.null(strelka2_vcf)) {
        snv.counts.dt <- strelka2counts(strelka2_vcf, seqnames_genome_width = seqnames_genome_width)
        meta.dt$snv_count <- snv.counts.dt[category == "snv_count", ]$counts
        meta.dt$snv_count_normal_vaf_greater0 <- snv.counts.dt[category == "snv_count_normal_vaf_greater0", ]$counts
    }

    if (is.null(sage_vcf) || sage_vcf == "") {
        warning("SAGE VCF not found as input, will only consider Strelka2 downstream...")
    } else {
        print("Found SAGE vcf, will use SAGE counts in meta file over Strelka2...")
        sage.snv.counts.dt <- sage_count(sage_vcf, genome = genome)
        meta.dt$snv_count <- sage.snv.counts.dt[category == "snv_count", ]$counts
        meta.dt$snv_count_normal_vaf_greater0 <- sage.snv.counts.dt[category == "snv_count_normal_vaf_greater0", ]$counts
    }

    sv_like_types_count <- NULL
    ## sv counts, purity, ploidy, loh, total genome length
    if (!is.null(jabba_gg)) {
        ## count := junctions + (loose ends / 2)
        gg <- readRDS(jabba_gg)
        if (nrow(gg$junctions$dt) > 0) {
            meta.dt$junction_count <- nrow(gg$junctions$dt[type != "REF", ])
        } else {
            warning("No junctions found in JaBbA graph, skipping SV counts...")
            meta.dt$junction_count <- 0
        }
        if (length(gg$loose) > 0) {
            meta.dt$loose_count <- nrow(as.data.table(gg$loose)[terminal == FALSE, ])
        } else {
            warning("No loose ends found in JaBbA graph, skipping loose counts...")
            meta.dt$loose_count <- 0
        }
        if (nrow(gg$edges$dt) > 0) {
            sv_like_types_count <- table(gg$edges$dt[type == "ALT" & class != "REF"]$class)
        } else {
            warning("No edges found in JaBbA graph, skipping SV-like types count...")
        }
        meta.dt[, sv_count := (junction_count + (loose_count / 2))]

        # get loh
        nodes.dt <- gg$nodes$dt
        nodes.dt[, seqnames := gsub("chr", "", seqnames)] # strip chr
        nodes.dt <- gg$nodes$dt[seqnames %in% seqnames_loh]
        totalseglen <- nodes.dt$width %>% sum()

        if ("cn.low" %in% names(nodes.dt)) {
            lohsegs <- nodes.dt[cn.low == 0, ] %>%
                .[cn.high > 0] %>%
                .$width %>%
                sum()
            loh_frc <- lohsegs / totalseglen
            meta.dt[, loh_fraction := loh_frc]
            meta.dt[, loh_seglen := lohsegs]
            meta.dt[, loh_total_genome := totalseglen]
        } else {
            meta.dt$loh_fraction <- "Not Allelic Jabba"
        }

        ## add purity and ploidy
        meta.dt$purity <- gg$meta$purity
        meta.dt$ploidy <- gg$meta$ploidy
        meta.dt$beta <- gg$meta$ploidy # karyograph not needed
        meta.dt$gamma <- 1 - gg$meta$purity # karyograph not needed
        ## add the total seqlengths by using the seqlengths in the jabba object
        nodes.gr <- gg$nodes$gr
        seqlengths.dt <- suppressWarnings(
            as.data.table(
                seqinfo(nodes.gr),
                keep.rownames = "seqnames"
            )
        ) # had to supressWarnings, says other arguments ignored
        seqlengths.dt[, seqnames := gsub("chr", "", seqnames)] # strip chr
        seqlengths.dt <- seqlengths.dt[seqnames %in% seqnames_genome_width, ]
        meta.dt$total_genome_length <- sum(seqlengths.dt$seqlengths)
    }

    # init qrp counts separately to pass to hrd since we want it in that tooltip
    qrp_counts <- data.table(qrpmin = 0, qrpmix = 0, qrppos = 0)

    # sv types count
    if (!is.null(complex) && !is.null(sv_like_types_count)) {
        complex <- readRDS(complex)
        if (!is.null(complex$meta$events$type)) {
            sv_types_count <- table(complex$meta$events$type)
            # update qrp counts
            if ("qrpmin" %in% names(sv_types_count)) {
                qrp_counts$qrpmin <- sv_types_count[["qrpmin"]]
            }
            if ("qrpmix" %in% names(sv_types_count)) {
                qrp_counts$qrpmix <- sv_types_count[["qrpmix"]]
            }
            if ("qrppos" %in% names(sv_types_count)) {
                qrp_counts$qrppos <- sv_types_count[["qrppos"]]
            }

            meta.dt$sv_types_count <- list(as.list(c(sv_types_count, sv_like_types_count)))
        }
    } else {
        warning("Complex events and JaBbA graph not found as inputs, skipping SV types count...")
    }

    ## Load beta/gamma for karyograph
    if (!is.null(coverage)) {
        rel2abs.cov <- skitools::rel2abs(readRDS(coverage),
            field = "foreground.X",
            purity = meta.dt$purity,
            ploidy = meta.dt$ploidy,
            return.params = T
        )
        meta.dt$cov_slope <- rel2abs.cov[1] %>% unname()
        meta.dt$cov_intercept <- rel2abs.cov[2] %>% unname()
    }

    if (!is.null(het_pileups_wgs)) {
        hets.read <- grab.hets(het_pileups_wgs) %>% gr2dt()
        hets.read <- dt2gr(hets.read[, .(count = sum(count)), by = c("seqnames", "start", "end")])
        rel2abs.hets <- skitools::rel2abs(hets.read,
            field = "count",
            purity = meta.dt$purity,
            ploidy = meta.dt$ploidy,
            return.params = T
        )
        meta.dt$hets_slope <- rel2abs.hets[1] %>% unname()
        meta.dt$hets_intercept <- rel2abs.hets[2] %>% unname()
    }

    ## add tmb
    if (("snv_count" %in% names(meta.dt)) & ("total_genome_length" %in% names(meta.dt))) {
        meta.dt[, tmb := (snv_count / (as.numeric(meta.dt$total_genome_length) / 1e6))]
        meta.dt[, tmb := round(tmb, digits = 3)]
    }

    ## add signatures
    meta.dt$deconstructsigs_sbs_fraction <- list()
    meta.dt$deletionInsertion <- list()
    meta.dt$sigprofiler_indel_fraction <- list()
    meta.dt$sigprofiler_indel_count <- list()
    meta.dt$sigprofiler_sbs_fraction <- list()
    meta.dt$sigprofiler_sbs_count <- list()
    meta.dt$signatures <- list()

    if (!is.null(signatures_pair_name)) {
        sig_sample_name <- signatures_pair_name
    } else {
        sig_sample_name <- pair
    }

    if (!is.null(sbs_deconstructSigs)) {
        signatures <- sbs_deconstructSigs2meta(sbs_file = sbs_deconstructSigs, sample = sig_sample_name)
        meta.dt$deconstructsigs_sbs_fraction <- list(as.list(signatures))
    }

    if (!is.null(indel_sigprofiler)) {
        deletionInsertion <- indels_sigprofiler2meta(indel_file = indel_sigprofiler, sample = sig_sample_name)
        meta.dt$deletionInsertion <- list(as.list(deletionInsertion[["indel_fraction"]]))
        meta.dt$sigprofiler_indel_fraction <- list(as.list(deletionInsertion[["indel_fraction"]]))
        meta.dt$sigprofiler_indel_count <- list(as.list(deletionInsertion[["indel_count"]]))
    }

    if (!is.null(sbs_sigprofiler)) {
        signatures <- sbs_sigprofiler2meta(sbs_file = sbs_sigprofiler, sample = sig_sample_name)
        meta.dt$sigprofiler_sbs_fraction <- list(as.list(signatures[["sbs_fraction"]]))
        meta.dt$sigprofiler_sbs_count <- list(as.list(signatures[["sbs_count"]]))
        meta.dt$signatures <- list(as.list(signatures[["sbs_fraction"]])) # Changed default to sigprofiler
    }

    if (!is.null(hrdetect)) {
        hrd <- readRDS(hrdetect)
        hrd_values <- data.table(
            dels_mh = hrd$indels_classification_table$del.mh.count,
            rs3 = hrd$exposures_rearr["RefSigR3", ],
            rs5 = hrd$exposures_rearr["RefSigR5", ],
            sbs3 = hrd$exposures_subs["SBS3", ],
            sbs8 = hrd$exposures_subs["SBS8", ],
            del_rep = hrd$indels_classification_table$del.rep.count
        )
        hrd_score <- hrd$hrdetect_output[1, "Probability"]
        meta.dt$hrd_score <- hrd_score
        meta.dt$hrd <- list(as.list(c(hrd_values, qrp_counts)))
    }

    if (!is.null(onenesstwoness)) {
        onetwo <- readRDS(onenesstwoness)
        meta.dt$b1_2 <- onetwo$ot_scores[, "SUM12"]
        meta.dt$b1 <- onetwo$ot_scores[, "BRCA1"]
        meta.dt$b2 <- onetwo$ot_scores[, "BRCA2"]
    }

    ## write the json
    if (write_json) {
        message(paste0("Writing json to ", out_file))
        write_json(
            meta.dt,
            out_file,
            pretty = TRUE,
            auto_unbox = TRUE
        )
        if (return_table) {
            return(meta.dt)
        }
    } else {
        return(meta.dt)
    }
}

#' @name create_distributions
#' @title create_distributions
#' @description
#'
#' Creates the background distributions for aggregated metadata KPIs across all patients
#'
#' @param case_reports_datadir directory containing all patient directories
#' @param filter_patients list() - subset of samples to filter on to create the distributions
#' @param cores number of cores to use for parallel processing
#' @param write_to_json boolean flag to indicate whether to write the files to JSON or not
#' @param common_dir path to the common directory in which to write the distribution JSONs (required if write_to_json is TRUE)
#' @param overwrite boolean flag to indicate whether to overwrite existing files (required if write_to_json is TRUE)
#' @return List of data tables containing distributions if write_to_json is FALSE, otherwise NULL
#' @export
create_distributions <- function(
    case_reports_datadir,
    filter_patients = NULL,
    cores = 1,
    write_to_json = FALSE,
    common_dir = NULL,
    overwrite = FALSE,
    haveSignatures = TRUE) {
    case_reports_data_folder <- paste0(case_reports_datadir, "/")

    files <- list.files(
        case_reports_data_folder,
        pattern = "metadata.json",
        recursive = TRUE,
        full.names = TRUE
    )

    patient_ids <- basename(dirname(files))

    meta.dt <- data.table(
        meta_json = files,
        patient_id = patient_ids
    )

    meta.dt <- meta.dt[file.exists(meta_json), ]

    if (!is.null(filter_patients)) {
        meta.dt <- meta.dt[patient_id %in% filter_patients, ]
    }

    jsons.lst <- mclapply(1:nrow(meta.dt), function(i) {
        meta.sub.dt <- meta.dt[i, ]
        read_meta_data_json(meta_json = meta.sub.dt$meta_json, patient_id = meta.sub.dt$patient_id)
    }, mc.cores = cores)
    jsons.dt <- rbindlist(jsons.lst, fill = TRUE)

    subset_and_setnames <- function(dt, cols, new_names) {
        if (all(cols %in% names(dt))) {
            return(dt[, ..cols] %>% setnames(., new_names))
        }
        return(NULL)
    }

    # add new metadata attributes here
    column_map <- list(
        snv.dt = list(
            cols = c("pair", "snv_count", "tumor_type"),
            new_names = c("pair", "value", "tumor_type")
        ),
        sv.dt = list(
            cols = c("pair", "sv_count", "tumor_type"),
            new_names = c("pair", "value", "tumor_type")
        ),
        loh.dt = list(
            cols = c("pair", "tumor_type", "loh_fraction", "loh_seglen", "loh_total_genome"),
            new_names = c("pair", "tumor_type", "value", "LOH_seg_len", "genome_width")
        ),
        ploidy.dt = list(
            cols = c("pair", "tumor_type", "ploidy", "purity"),
            new_names = c("pair", "tumor_type", "value", "purity")
        ),
        purity.dt = list(
            cols = c("pair", "tumor_type", "ploidy", "purity"),
            new_names = c("pair", "tumor_type", "ploidy", "value")
        ),
        cov_var.dt = list(
            cols = c("pair", "tumor_type", "dlrs"),
            new_names = c("pair", "tumor_type_mod", "value")
        ),
        tmb.dt = list(
            cols = c("pair", "tmb", "tumor_type"),
            new_names = c("pair", "value", "tumor_type")
        )
    )

    results <- lapply(names(column_map), function(name) {
        subset_and_setnames(jsons.dt, column_map[[name]]$cols, column_map[[name]]$new_names)
    })

    names(results) <- names(column_map)

    # signatures
    if (haveSignatures) {
        sigs.lst <- convert_signature_meta_json(jsons.dt = jsons.dt, cores = cores)

        sigs_names <- c(
            "deconstruct_sigs.dt",
            "sigprofilerassignment_indels.dt",
            "sigprofilerassignment_indels_counts.dt",
            "sigprofilerassignment_sbs_counts.dt",
            "sigprofilerassignment_sbs.dt"
        )

        sigs_results <- lapply(sigs_names, function(name) {
            if (name %in% names(sigs.lst)) {
                return(sigs.lst[[name]])
            }
            return(NULL)
        })

        json.lst <- c(results, sigs_results)

        names(json.lst) <- c(
            "snvCount",
            "svCount",
            "lohFraction",
            "ploidy",
            "purity",
            "coverageVariance",
            "tmb",
            "sbs",
            "sigprofiler_indel_fraction",
            "sigprofiler_indel_count",
            "sigprofiler_sbs_count",
            "sigprofiler_sbs_fraction"
        )
    } else {
        json.lst <- results
        names(json.lst) <- c(
            "snvCount",
            "svCount",
            "lohFraction",
            "ploidy",
            "purity",
            "coverageVariance",
            "tmb"
        )
    }

    if (write_to_json) {
        if (is.null(common_dir)) {
            stop("common_dir must be provided if write_to_json is TRUE")
        }
        write_distributions_to_json(json.lst, common_dir, cores, overwrite, haveSignatures)
        return(NULL)
    } else {
        return(json.lst)
    }
}

#' @name write_signature_jsons
#' @title write_signature_jsons
#' @description
#' Writes the different distribution JSONs for signatures from the list outputted by convert_signature_meta_json
#'
#' @param signatures_list list outputted from convert_sign
#' @param common_folder location of the signatures folder
#' @param cores cores for generating the signature distribution files
#' @return NULL
#' @export
#' @author Shihab Dider, Stanley Clarke
write_signature_jsons <- function(signatures_list, common_folder, cores = 1) {
    write_signature_json <- function(sig.dt, folder_path, sig) {
        sig.dt[, sig := sig]
        write_json(
            sig.dt,
            paste0(common_folder, folder_path, sig, ".json"),
            pretty = TRUE
        )
    }

    process_signatures <- function(sig_data, folder_path) {
        sig_add <- names(sig_data) %>% grep("pair|tumor_type", ., invert = TRUE, value = TRUE)
        empty.lst <- mclapply(sig_add, function(sig) {
            cols_keep <- c("pair", "tumor_type", sig)
            sig.dt <- sig_data[, ..cols_keep] %>% setnames(., c("pair", "tumor_type", "value"))
            write_signature_json(sig.dt, folder_path, sig)
            return(NULL)
        }, mc.cores = cores)
    }

    if ("sbs" %in% names(signatures_list)) {
        process_signatures(
            signatures_list[["sbs"]],
            "signatures/deconstructsigs_sbs_fraction/"
        )
    }
    if ("sigprofiler_indel_fraction" %in% names(signatures_list)) {
        process_signatures(
            signatures_list[["sigprofiler_indel_fraction"]],
            "signatures/sigprofiler_indel_fraction/"
        )
    }
    if ("sigprofiler_indel_count" %in% names(signatures_list)) {
        process_signatures(
            signatures_list[["sigprofiler_indel_count"]],
            "signatures/sigprofiler_indel_count/"
        )
    }
    if ("sigprofiler_sbs_count" %in% names(signatures_list)) {
        process_signatures(
            signatures_list[["sigprofiler_sbs_count"]],
            "signatures/sigprofiler_sbs_count/"
        )
    }
    if ("sigprofiler_sbs_fraction" %in% names(signatures_list)) {
        process_signatures(
            signatures_list[["sigprofiler_sbs_fraction"]],
            "signatures/sigprofiler_sbs_fraction/"
        )
        process_signatures(
            signatures_list[["sigprofiler_sbs_fraction"]],
            "signatures/sbs/"
        )
    }
}

#' @name write_distributions_to_json
#' @title write_distributions_to_json
#' @description
#'
#' Writes the distributions data table to JSON files
#'
#' @param distributions List of data tables containing distributions
#' @param common_dir path to the common directory in which to write the distribution JSONs
#' @param cores number of cores to use for parallel processing
#' @param overwrite boolean flag to indicate whether to overwrite existing files
#' @return NULL
#' @export
write_distributions_to_json <- function(
    distributions,
    common_dir,
    cores = 1,
    overwrite = FALSE,
    haveSignatures = TRUE) {
    common_folder <- paste0(common_dir, "/")

    if (!file.exists(common_folder)) {
        stop("common_folder does not exist. Make the directory first")
    }

    signature_folder_paths <- paste0(
        common_folder,
        c(
            "signatures/sbs/",
            "signatures/deconstructsigs_sbs_fraction/",
            "signatures/insertionDeletion/",
            "signatures/sigprofiler_indel_fraction/",
            "signatures/sigprofiler_indel_count/",
            "signatures/sigprofiler_sbs_count/",
            "signatures/sigprofiler_sbs_fraction/"
        )
    )

    if (!all(file.exists(signature_folder_paths))) {
        empty.lst <- mclapply(signature_folder_paths, function(folder_path) {
            cmd <- paste0("mkdir -p ", folder_path)
            print(paste0("Making directory ", folder_path))
            system(cmd)
            return(NULL)
        }, mc.cores = cores)
    }

    write_json_file <- function(data, file_path) {
        if (overwrite || !file.exists(file_path)) {
            write_json(data, file_path, pretty = TRUE)
        } else {
            message(paste0("File ", file_path, " already exists. Skipping."))
        }
    }

    message(paste0("writing jsons to ", common_folder))
    write_json_file(distributions$snvCount, paste0(common_folder, "/snvCount.json"))
    write_json_file(distributions$svCount, paste0(common_folder, "/svCount.json"))
    write_json_file(distributions$lohFraction, paste0(common_folder, "/lohFraction.json"))
    write_json_file(distributions$ploidy, paste0(common_folder, "/ploidy.json"))
    write_json_file(distributions$purity, paste0(common_folder, "/purity.json"))
    write_json_file(distributions$coverageVariance, paste0(common_folder, "/coverageVariance.json"))
    write_json_file(distributions$tmb, paste0(common_folder, "/tmb.json"))

    if (haveSignatures) {
        write_signature_jsons(
            signatures_list = distributions,
            common_folder = common_folder,
            cores = cores
        )
    }
}

#' @name load_distributions
#' @title load_distributions
#' @description
#'
#' function to read in distributions files, same format returned as create_distributions with write_jsons = FALSE
#'
#' @param common_dir path to a folder to write all 7 jsons
#' @param filter_pateints list of samples to filter on to read in the distributions
#' @return NULL
#' @export
#' @author Stanley Clarke
load_distributions <- function(
    common_dir,
    filter_patients = NULL) {
    common_dir <- paste0(common_dir, "/")
    files.lst <- paste0(
        common_dir,
        c(
            "coverageVariance.json",
            "lohFraction.json",
            "ploidy.json",
            "purity.json",
            "snvCount.json",
            "svCount.json",
            "tmb.json"
        )
    )
    files.dt <- data.table(
        name = c(
            "coverageVariance",
            "lohFraction",
            "ploidy",
            "purity",
            "snvCount",
            "svCount",
            "tmb"
        ),
        file = files.lst
    )

    json.lst <- lapply(setNames(nm = files.dt$name), function(name_i) {
        file_i <- files.dt[name == name_i, ]$file
        if (!is.null(filter_patients)) {
            json.dt <- as.data.table(fromJSON(file_i))
            json.sub.dt <- json.dt[pair %in% filter_patients, ]
            return(json.sub.dt)
        }
        return(as.data.table(fromJSON(file_i)))
    })
    return(json.lst)
}

#' @name add_patient_to_distributions
#' @title add_patient_to_distributions
#' @description
#'
#' Adds a new patient to the existing distributions and writes the updated distributions to JSON files.
#'
#' @param metadata_json_path Path to the new patient's metadata.json file.
#' @param common_dir Path to the common directory where the distributions are stored.
#' @param cores Number of cores to use for parallel processing.
#' @param overwrite Boolean flag to indicate whether to overwrite existing files.
#' @return Updated distributions.
#' @export
#' @author Shihab Dider
add_patient_to_distributions <- function(
    metadata_json_path,
    common_dir,
    cores = 1,
    overwrite = FALSE,
    write_to_json = FALSE) {
    distributions <- load_distributions(common_dir)

    new_patient_data <- fromJSON(metadata_json_path)

    # Extract patient ID from the metadata.json path
    patient_id <- basename(dirname(metadata_json_path))

    # Update distributions with the new patient's data
    update_distribution <- function(distribution, new_data, key) {
        if (!is.null(distribution)) {
            new_entry <- data.table(
                pair = patient_id,
                value = new_data[[key]],
                tumor_type = new_data$tumor_type
            )
            return(rbind(distribution, new_entry, fill = TRUE))
        }
        return(distribution)
    }

    distributions$snvCount <- update_distribution(distributions$snvCount, new_patient_data, "snv_count")
    distributions$svCount <- update_distribution(distributions$svCount, new_patient_data, "sv_count")
    distributions$lohFraction <- update_distribution(distributions$lohFraction, new_patient_data, "loh_fraction")
    distributions$ploidy <- update_distribution(distributions$ploidy, new_patient_data, "ploidy")
    distributions$purity <- update_distribution(distributions$purity, new_patient_data, "purity")
    distributions$coverageVariance <- update_distribution(distributions$coverageVariance, new_patient_data, "dlrs")
    distributions$tmb <- update_distribution(distributions$tmb, new_patient_data, "tmb")

    # Update signature distributions
    update_signature_distribution <- function(distribution, new_data, key_prefix) {
        if (!is.null(distribution)) {
            sig_keys <- names(new_data)[grep(paste0("^", key_prefix), names(new_data))]
            for (sig_key in sig_keys) {
                new_entry <- data.table(
                    pair = patient_id,
                    value = new_data[[sig_key]],
                    tumor_type = new_data$tumor_type
                )
                distribution <- rbind(distribution, new_entry, fill = TRUE)
            }
            return(distribution)
        }
        return(distribution)
    }

    distributions$sbs <- update_signature_distribution(
        distributions$sbs,
        new_patient_data,
        "sbs"
    )
    distributions$sigprofiler_indel_fraction <- update_signature_distribution(
        distributions$sigprofiler_indel_fraction,
        new_patient_data,
        "sigprofiler_indel_fraction"
    )
    distributions$sigprofiler_indel_count <- update_signature_distribution(
        distributions$sigprofiler_indel_count,
        new_patient_data,
        "sigprofiler_indel_count"
    )
    distributions$sigprofiler_sbs_count <- update_signature_distribution(
        distributions$sigprofiler_sbs_count,
        new_patient_data,
        "sigprofiler_sbs_count"
    )
    distributions$sigprofiler_sbs_fraction <- update_signature_distribution(
        distributions$sigprofiler_sbs_fraction,
        new_patient_data,
        "sigprofiler_sbs_fraction"
    )

    # Save the updated distributions
    if (write_to_json) {
        write_distributions_to_json(distributions, common_dir, cores, overwrite)
    } else {
        print("Returning updated distributions. If you want to save them, set write_to_json = TRUE.")
        return(distributions)
    }
}

#' @name bw_temp
#' @title bw_temp
#' @description
#'
#' function to create data.table for bigwig files for pgvdb
#'
#' @param patient_id patient.id to be added to pgvdb
#' @param order optional entry if you order plots with a column order
#' @param x path to a granges rds or a granges object
#' @param ref reference to use for pgvdb
#' @param chart_type defaultChartType for pgvdb, default is area, can also be scatterplot
#' @param visible TRUE/FALSE whether the plot is hidden or showing in pgv
#' @param title title of the plot in pgvdb
#' @param type bigwig, should not change for this plot type
#' @param field field to plot as the y value, default is foreground
#' @param overwrite TRUE/FALSE to overwrite an existing bw
#' @return NULL
#' @export
#' @author Stanley Clarke

bw_temp <- function(
    patient_id = NA,
    order = NA,
    x = list(NA),
    ref = NA,
    chart_type = "area",
    visible = TRUE,
    title = NA,
    type = "bigwig",
    field = "foreground",
    overwrite = FALSE) {
    dt1 <- data.table(
        patient.id = patient_id,
        visible = visible,
        x = x,
        type = type,
        field = field,
        ref = ref,
        title = title,
        order = order,
        defaultChartType = chart_type,
        overwrite = overwrite
    )
    return(dt1)
}


#' @name arrow_temp
#' @title arrow_temp
#' @description
#'
#' function to create data.table for arrow files for pgvdb
#'
#' @param patient_id patient.id to be added to pgvdb
#' @param order optional entry if you order plots with a column order
#' @param x path to a granges rds or a granges object
#' @param ref reference to use for pgvdb
#' @param chart_type defaultChartType for pgvdb, default is scatterplot, can also be area
#' @param visible TRUE/FALSE whether the plot is hidden or showing in pgv
#' @param title title of the plot in pgvdb
#' @param type scatterplot, should not change for this plot type
#' @param field field to plot as the y value, default is foreground
#' @param overwrite TRUE/FALSE to overwrite an existing arrow
#' @return NULL
#' @export
#' @author Stanley Clarke

arrow_temp <- function(
    patient_id = NA,
    order = NA,
    x = list(NA),
    ref = NA,
    source = "coverage.arrow",
    chart_type = "scatterplot",
    visible = TRUE,
    title = NA,
    type = "scatterplot",
    field = "foreground",
    overwrite = FALSE) {
    dt1 <- data.table(
        patient.id = patient_id,
        visible = visible,
        x = x,
        type = type,
        field = field,
        ref = ref,
        source = source,
        title = title,
        order = order,
        defaultChartType = chart_type,
        overwrite = overwrite
    )
    return(dt1)
}

#' @name genome_temp
#' @title genome_temp
#' @description
#'
#' function to create data.table for genome graphs for pgvdb
#'
#' @param patient_id patient.id to be added to pgvdb
#' @param order optional entry if you order plots with a column order
#' @param x path to a JaBbA ggraph object or object itself as a list
#' @param ref reference to use for pgvdb
#' @param visible TRUE/FALSE whether the plot is hidden or showing in pgv
#' @param title title of the plot in pgvdb
#' @param max.cn override max cn of 100
#' @param type genome, can be changed to allelic to render allelic graphs
#' @param annotation default is list of SVs, make null if no annotations present in object
#' @param overwrite TRUE/FALSE to overwrite an existing genome json
#' @return NULL
#' @export
#' @author Stanley Clarke

genome_temp <- function(
    patient_id = NA,
    order = NA,
    x = list(NA),
    ref = NA,
    source = "genome.json",
    type = "genome",
    visible = TRUE,
    title = NA,
    max.cn = NULL,
    annotation = list(c("bfb", "chromoplexy", "chromothripsis", "del", "dm", "cpxdm", "dup", "pyrgo", "rigma", "simple", "tic", "tyfonas")),
    overwrite = FALSE) {
    # use type = allelic to make a color a genome graph
    dt1 <- data.table(
        patient.id = patient_id,
        type = type,
        visible = visible,
        title = title,
        x = x,
        ref = ref,
        source = source,
        max.cn = max.cn,
        order = order,
        annotation = annotation,
        overwrite = overwrite
    )
    return(dt1)
}

#' @name walks_temp
#' @title walks_temp
#' @description
#'
#' function to create data.table for walks for pgvdb
#'
#' @param patient_id patient.id to be added to pgvdb
#' @param order optional entry if you order plots with a column order
#' @param x walks object as a list
#' @param ref reference to use for pgvdb
#' @param visible TRUE/FALSE whether the plot is hidden or showing in pgv
#' @param title title of the plot in pgvdb
#' @param type walk, do not change for this plot type
#' @param annotation default is list of SVs, make null if no annotations present in object
#' @param overwrite TRUE/FALSE to overwrite an existing genome json
#' @param tag optional argument, can be binset to override y spacing in pgv
#' @return NULL
#' @export
#' @author Stanley Clarke

walks_temp <- function(
    patient_id = NA,
    order = NA,
    x = list(NA),
    ref = NA,
    source = "walks.json",
    type = "walk",
    visible = TRUE,
    title = NA,
    tag = NA,
    overwrite = FALSE) {
    dt1 <- data.table(
        patient.id = patient_id,
        visible = visible,
        x = x,
        type = type,
        order = order,
        ref = ref,
        source = source,
        title = title,
        overwrite = overwrite
    )
    return(dt1)
}

#' @name mutations_temp
#' @title mutations_temp
#' @description
#'
#' function to create data.table for mutations for pgvdb
#'
#' @param patient_id patient.id to be added to pgvdb
#' @param order optional entry if you order plots with a column order
#' @param x mutations object as a list
#' @param ref reference to use for pgvdb
#' @param source name of json file to save as. For 'somatic', somatic.json is suggested whereas for 'germline', germline.json is suggested. by default, plot will save as mutations.json.
#' @param subtype enum of 'somatic' or 'germline' which has downstream effects on mutation generation code
#' @param visible TRUE/FALSE whether the plot is hidden or showing in pgv
#' @param title title of the plot in pgvdb
#' @param type mutations, do not change for this plot type
#' @param annotation default is list of SVs, make null if no annotations present in object
#' @param overwrite TRUE/FALSE to overwrite an existing genome json
#' @param tag optional argument, can be binset to override y spacing in pgv
#' @return NULL
#' @export
#' @author Stanley Clarke

mutations_temp <- function(
    patient_id = NA,
    order = NA,
    x = list(NA),
    ref = NA,
    source = "mutations.json",
    field,
    type = "mutations",
    subtype = "somatic",
    visible = TRUE,
    title = NA,
    tag = NA,
    overwrite = FALSE) {
    dt1 <- data.table(
        patient.id = patient_id,
        visible = visible,
        type = type,
        x = x,
        field = field,
        order = order,
        ref = ref,
        source = source,
        subtype = subtype,
        title = title,
        overwrite = overwrite
    )
    return(dt1)
}

#' @name ppfit_temp
#' @title ppfit_temp
#' @description
#'
#' function to create data.table for ppfit for pgvdb
#'
#' @param patient_id patient.id to be added to pgvdb
#' @param order optional entry if you order plots with a column order
#' @param x jabba.rds
#' @param ref reference to use for pgvdb
#' @param visible TRUE/FALSE whether the plot is hidden or showing in pgv
#' @param title title of the plot in pgvdb
#' @param type ppfit, do not change for this plot type
#' @param annotation default is NULL, currently won't work for ppfit
#' @param overwrite TRUE/FALSE to overwrite an existing genome json
#' @param tag optional argument, can be binset to override y spacing in pgv
#' @return NULL
#' @export
#' @author Stanley Clarke

ppfit_temp <- function(
    patient_id = NA,
    order = NA,
    x = list(NA),
    ref = NA,
    source = "ppfit.json",
    type = "ppfit", # use type = allelic to make a colored genome graph
    visible = TRUE,
    title = NA,
    annotation = NULL,
    overwrite = FALSE) {
    dt1 <- data.table(
        patient.id = patient_id,
        type = type,
        visible = visible,
        title = title,
        x = x,
        ref = ref,
        source = source,
        order = order,
        annotation = annotation,
        overwrite = overwrite
    )
    return(dt1)
}

#' @name create_ppfit_json
#' @title create_ppfit_json
#' @description
#'
#' function to create segmentation plots in case reports
#'
#' @param balanced_gg the balanced_gg ggraph or the path to it
#' @param cov path to associated coverage, if null, will try to pull path from pairs table
#' @param path_obj path to flow directory, if null, must supply coverage path
#' @param settings_json path to settings.json (from PGV)
#' @param ref reference name
#' @param max_na max.na used with JaBbA
#' @param out_file location to write json
#' @param write_json TRUE/FALSE whether to write the json
#' @param return_table TRUE/FALSE whether to return the data
#' @param overwrite TRUE/FALSE whether to overwrite existing file
#' @param cores cores for JaBbA:::segstats
#' @return NULL or segstats table
#' @export
#' @author Stanley Clarke, Tanubrata Dey

create_ppfit_json <- function(
    balanced_gg,
    cov_path = NULL,
    path_obj = NULL,
    settings_json,
    ref = NULL,
    max_na = NULL,
    out_file = NULL,
    write_json = TRUE,
    overwrite = FALSE,
    return_table = FALSE,
    cores = 1) {
    if (!is.null(out_file)) {
        if (!overwrite) {
            if (file.exists(out_file)) {
                print("Output already exists! - skipping sample")
                return(NA)
            }
        }
    }

    #' drcln_cov = readRDS(thisp$decomposed_cov[i])
    ## make this work with complex where the cov file was not an input and with jabba_gg
    ## x = path_obj %>% sniff %>% inputs %>% select(CovFile, maxna) #get coverage that was used for the jabba run
    if (!is.null(cov_path)) {
        cov <- readRDS(cov_path)
    } else if (!is.null(path_obj)) {
        inputs.dt <- path_obj %>%
            sniff() %>%
            inputs()
        if (!any(grepl("CovFile", names(inputs.dt))) & any(grepl("jabba", names(inputs.dt)))) {
            x <- inputs.dt$jabba %>%
                sniff() %>%
                inputs() %>%
                .[, .(CovFile, maxna)]
        } else if (!any(grepl("CovFile", names(inputs.dt))) & any(grepl("jab", names(inputs.dt)))) {
            x <- inputs.dt$jab %>%
                sniff() %>%
                inputs() %>%
                .[, .(CovFile, maxna)]
        } else {
            x <- path_obj %>%
                sniff() %>%
                inputs() %>%
                .[, .(CovFile, maxna)]
        }
        cov <- readRDS(x$CovFile)
        max_na <- as.numeric(x$maxna)
    } else {
        stop("Must supply either coverage path or flow directory path object")
    }

    if (is.null(ref)) {
        warning("ref was not passed, setting to default: hg19")
        ref <- "hg19"
    }

    if (is.null(max_na)) {
        warning("max_na was not passed/found, setting to default: 0.9")
        max_na <- 0.9
    }

    if ("ratio" %in% names(mcols(cov))) {
        message(paste0("Raw 'cov.rds' was used as input for JaBbA ", path_obj, ", will consider field as 'ratio''\n"))
        field <- "ratio"
    } else if ("foreground" %in% names(mcols(cov))) {
        message(paste0("Drycleaned 'drycleaned.cov.rds' was used as input for JaBbA ", path_obj, ", will consider field as 'foreground''\n"))
        field <- "foreground"
    }
    if (!(field %in% c("ratio", "foreground"))) {
        stop("Cov file is not clear. Ratio nor foreground in the the columns of the coverage file")
    }
    ## need to replace NaN with NA or JaBbA:::segstats breaks
    if (field == "ratio") {
        cov$ratio <- gsub("NaN", NA, cov$ratio) %>% as.numeric()
    } else if (field == "foreground") {
        cov$foreground <- gsub("NaN", NA, cov$foreground) %>% as.numeric()
    }
    if (is(balanced_gg, "gGraph")) {
        balanced_gg_gr <- balanced_gg$nodes$gr
    } else if (is(balanced_gg, "character")) {
        balanced_gg_gr <- readRDS(balanced_gg)$nodes$gr
    }

    segstats <- JaBbA:::segstats(balanced_gg_gr,
        cov,
        field = field,
        prior_weight = 1,
        max.chunk = 1e8,
        ## subsample = subsample,
        mc.cores = cores,
        verbose = FALSE,
        max.na = max_na,
        lp = FALSE
    )
    segstats.dt <- gr2dt(segstats)
    names(segstats.dt) <- gsub("\\.", "_", names(segstats.dt))
    if (write_json) {
        message("Cleaning up for writing ppfit to json")
        seq_lengths <- gGnome::parse.js.seqlengths(
            settings_json,
            js.type = "PGV",
            ref = ref
        )
        segstats.gr <- GRanges(segstats.dt, seqlengths = seq_lengths) %>% trim()
        ggraph <- readRDS(balanced_gg)
        ggraph2 <- gG(nodes = segstats.gr, edges = ggraph$edges$dt)
        fields.keep <- names(segstats.dt) %>% grep("cn", ., invert = TRUE, value = TRUE)
        ggraph2$set(y.field = "cn")
        ## check for overlap in sequence names
        ggraph.reduced <- ggraph2[seqnames %in% names(seq_lengths)]
        if (length(ggraph.reduced) == 0) {
            stop(sprintf(
                'There is no overlap between the sequence names in the reference
                    used by PGV and the sequences in your gGraph. Here is an
                    example sequence from your gGraph: "%s". And here is an
                    example sequence from the reference used by gGnome.js: "%s"',
                seqlevels(ggraph.reduced$nodes$gr)[1], names(seq_lengths)[1]
            ))
        }
        message(paste0("Writing json to ", out_file))
        gGnome::refresh(ggraph.reduced)$json(
            filename = out_file,
            verbose = TRUE,
            maxcn = 500,
            nfields = fields.keep,
            save = TRUE
        )
    }
    if (return_table) {
        return(segstats.dt)
    }
}


#' @name cov_abs
#' @title cov_abs
#' @description
#'
#' function to run rel2abs on coverage taking either a ggraph or purity & ploidy
#'
#' @param dryclean_cov dryclean coverage
#' @param jabba_gg optional jabba ggraph. If null, needs purity or ploidy
#' @param purity optional purity value. If null needs ggraph
#' @param ploidy optional ploidy value. If null needs ggraph
#' @param field column in granges to convert with rel2abs
#' @param new_col new column to add to add for the converted rel2abs
#' @return NULL or segstats table
#' @export
#' @author Stanley Clarke
cov2abs <- function(
    dryclean_cov,
    jabba_gg = NULL,
    purity = NULL,
    ploidy = NULL,
    field = "foreground",
    new_col = "foregroundabs") {
    cov_gr <- readRDS(dryclean_cov)
    if (!is.null(jabba_gg)) {
        gg <- readRDS(jabba_gg)
        purity <- gg$meta$purity
        ploidy <- gg$meta$ploidy
    }
    if (!is.null(purity) && !is.null(ploidy)) {
        purity <- purity
        ploidy <- ploidy
    }
    mcols(cov_gr)[new_col] <- rel2abs(
        gr = cov_gr,
        purity = purity,
        ploidy = ploidy,
        field = field
    )
    return(cov_gr)
}

#' @name cov2arrow_pgv
#' @title cov2arrow_pgv
#' @description
#'
#' function to create arrow with rel2abs from coverage file. Can take a ggraph or specified purity and ploidy to calculate rel2abs
#'
#' @param patient.id patient name to add to pgv
#' @param dryclean_cov dryclean coverage
#' @param ref reference for pgv
#' @param field column in granges to convert with rel2abs
#' @param jabba_gg optional jabba ggraph. If null, needs purity or ploidy
#' @param purity optional purity value. If null needs ggraph
#' @param ploidy optional ploidy value. If null needs ggraph
#' @param title title for ggraph
#' @param mask optional mask to use after rebinning
#' @param title title for ggraph
#' @param ref reference to use for pgvdb
#' @param title optional title for plot for pgvdb
#' @param seq.fix optional seqlengths vector to fix granges seqlengths
#' @param chart_type default chart type for plot in pgvdb. scatterplot or area
#' @param visible whether the plot is visible in pgv
#' @param field column to do rel2abs on
#' @param new_col new column after rel2abs
#' @param overwrite whether to overwrite the current bigwig
#' @param order optional order if using a column order in your pgvdb object to sort
#' @param binsize size to rebin coverages, default 1e4
#' @return NULL or segstats table
#' @export
#' @author Stanley Clarke
cov2arrow_pgv <- function(
    patient.id,
    dryclean_cov,
    jabba_gg = NULL,
    purity = NULL,
    ploidy = NULL,
    mask = NULL,
    ref,
    title = NA,
    seq.fix = NULL,
    chart_type = "scatterplot",
    visible = TRUE,
    field = "foreground",
    new_col = "foregroundabs",
    overwrite = FALSE,
    order = NA,
    binsize = 1e4) {
    if (!is.null(jabba_gg)) {
        cov_gr <- cov2abs(dryclean_cov, jabba_gg, field = field, new_col = new_col)
    } else if (!is.null(purity) && !is.null(ploidy)) {
        cov_gr <- cov2abs(dryclean_cov, purity = purity, ploidy = ploidy)
    } else {
        cov_gr <- readRDS(dryclean_cov)
    }
    if (!is.null(mask)) {
        cov_gr <- gr.val(cov_gr, mask, "mask")
        cov_gr <- cov_gr %Q% (is.na(mask))
        cov_gr$mask <- NULL
    }
    cov_gr2 <- rebin(cov_gr, binsize, field = field)
    if (length(cov_gr2$foregroundabs[cov_gr2$foregroundabs < 0]) > 0) {
        cov_gr2$foregroundabs[cov_gr2$foregroundabs < 0] <- 0
    }
    add.dt <- arrow_temp(patient_id = patient.id, ref = ref, field = field, x = list(cov_gr2), title = title, overwrite = overwrite, order = NA, chart_type = chart_type, visible = visible)
    return(add.dt)
}

#' @name cov2bw_pgv
#' @title cov2bw_pgv
#' @description
#'
#' function to create bigwig with rel2abs from coverage file. Can take a ggraph or specified purity and ploidy to calculate rel2abs
#'
#' @param patient.id patient name to add to pgv
#' @param dryclean_cov dryclean coverage
#' @param ref reference for pgv
#' @param field column in granges to convert with rel2abs
#' @param jabba_gg optional jabba ggraph. If null, needs purity or ploidy
#' @param purity optional purity value. If null needs ggraph
#' @param ploidy optional ploidy value. If null needs ggraph
#' @param mask mask to remove coverages from
#' @param title title for ggraph
#' @param ref reference to use for pgvdb
#' @param title optional title for plot for pgvdb
#' @param seq.fix optional seqlengths vector to fix granges seqlengths
#' @param chart_type default chart type for plot in pgvdb. scatterplot or area
#' @param visible whether the plot is visible in pgv
#' @param field column to do rel2abs on
#' @param new_col new column after rel2abs
#' @param overwrite whether to overwrite the current bigwig
#' @param order optional order if using a column order in your pgvdb object to sort
#' @param mask optional mask to use after rebinning
#' @return NULL or segstats table
#' @export
#' @author Stanley Clarke

## function to not rebin using higlass but mask
cov2bw_pgv <- function(
    patient.id,
    dryclean_cov,
    jabba_gg = NULL,
    purity = NULL,
    ploidy = NULL,
    mask = NULL,
    ref,
    title = NA,
    seq.fix = NULL,
    chart_type = "scatterplot",
    visible = TRUE,
    field = "foreground",
    new_col = "foregroundabs",
    overwrite = FALSE,
    order = NA) {
    if (!is.null(jabba_gg)) {
        cov_gr <- cov2abs(dryclean_cov, jabba_gg, field = field, new_col = new_col)
    } else if (!is.null(purity) && !is.null(ploidy)) {
        cov_gr <- cov2abs(dryclean_cov, purity = purity, ploidy = ploidy)
    } else {
        cov_gr <- readRDS(dryclean_cov)
    }
    ## cov_gr = cov2abs(dryclean_cov, jabba_gg, field = field, new_col = new_col)
    if (!is.null(mask)) {
        cov_gr <- gr.val(cov_gr, mask, "mask")
        cov_gr <- cov_gr %Q% (is.na(mask))
        cov_gr$mask <- NULL
    }
    if (length(cov_gr$foregroundabs[cov_gr$foregroundabs < 0]) > 0) {
        cov_gr$foregroundabs[cov_gr$foregroundabs < 0] <- 0
    }
    if (!is.null(seq.fix)) {
        ## fix seqlengths to specified seqlengths
        cov_gr <- GRanges(as.data.table(cov_gr), seqlengths = seq.fix) %>% trim()
    }
    add.dt <- bw_temp(
        patient_id = patient.id,
        ref = ref,
        field = new_col,
        x = list(cov_gr),
        title = title,
        overwrite = overwrite,
        order = NA,
        chart_type = chart_type,
        visible = visible
    )
    return(add.dt)
}

#' @name indels_sigprofiler2meta
#' @title indels_sigprofiler2meta
#' @description
#' function to generate metadata indel signature list from sigprofiler output
#'
#' @param indel_file path to Assignment_Solution_Activities.txt, output of sigprofiler
#' @param sample sample or pair name-needed to subset to the correct sample
#'
#' @return
#' @export
#' @author Sukanya Panja, Stanley Clarke
indels_sigprofiler2meta <- function(indel_file, sample) {
    ## indel signatures first
    indel_sig <- fread(indel_file)
    indel_sig[, pair := gsub("_somatic", "", Samples)]
    indel_sig[, Samples := NULL]
    indel_sig <- indel_sig[pair == sample, ]
    indel_sig_avg <- copy(indel_sig)
    indel_sig_avg[, pair := NULL]
    row_sum <- indel_sig_avg %>% rowSums()
    indels.dt <- melt.data.table(indel_sig_avg, measure.vars = names(indel_sig_avg)) %>% setnames(., c("signature", "value"))
    indels.dt[, avg_value := value / row_sum]
    indels.vect <- indels.dt$avg_value
    names(indels.vect) <- indels.dt$signature
    indels_counts.vect <- indels.dt$value
    names(indels_counts.vect) <- indels.dt$signature
    indels.lst <- list(indels.vect, indels_counts.vect)
    names(indels.lst) <- c("indel_fraction", "indel_count")
    return(indels.lst)
}


#' @name sbs_deconstructSigs2meta
#' @title sbs_deconstructSigs2meta
#' @description
#' generates list of SBS signatures from the output of deconstruct sigs to add to metadata
#'
#' @param sbs_file path to rds outpyut from deconstructSigs
#' @param sample sample or pair name-needed to subset to the correct sample
#'
#' @return
#' @export
#' @author Sukanya Panja, Stanley Clarke
sbs_deconstructSigs2meta <- function(sbs_file, sample) {
    sig_file <- readRDS(sbs_file)
    weights <- as.data.table(sig_file$weights)
    weights <- t(weights)
    sigs.dt <- as.data.table(weights, keep.rownames = "Signature") %>% setnames(., c("Signature", "weights"))
    sigs.dt[, Signature := gsub("Signature.", "SBS", Signature)]
    sigs.vect <- sigs.dt$weights
    names(sigs.vect) <- sigs.dt$Signature
    ## pair.vect = sample
    ## names(pair.vect) = "pair"
    ## sigs.vect = c(pair.vect,sigs.vect)
    ## signatures = list(signatures = as.list(sigs.vect))
    return(sigs.vect)
}


#' @name sbs_sigprofiler2meta
#' @title sbs_sigprofiler2meta
#' @description
#' generates list of SBS signature counts and percentages from the output of deconstruct sigs to add to metadata
#'
#' @param sbs_file path to Assignment_Solution_Activities.txt, output of sigprofiler
#' @param sample sample or pair name-needed to subset to the correct sample
#'
#' @return list of length 2, first being the signatures as a fraction and second being signatures as counts
#' @export
#' @author Stanley Clarke
sbs_sigprofiler2meta <- function(sbs_file, sample) {
    sig.dt <- fread(sbs_file)
    sig.dt[, pair := gsub("_somatic", "", Samples)]
    sig.dt <- sig.dt[pair == sample, ]
    sig.dt[, Samples := NULL]
    sig.dt_avg <- copy(sig.dt)
    sig.dt_avg[, pair := NULL]
    row_sum <- sig.dt_avg %>% rowSums()
    sigs.dt <- melt.data.table(sig.dt_avg, measure.vars = names(sig.dt_avg)) %>% setnames(., c("signature", "value"))
    sigs.dt[, avg_value := value / row_sum]
    sigs.vect <- sigs.dt$avg_value
    names(sigs.vect) <- sigs.dt$signature
    sigs_counts.vect <- sigs.dt$value
    names(sigs_counts.vect) <- sigs.dt$signature
    sigs.lst <- list(sigs.vect, sigs_counts.vect)
    names(sigs.lst) <- c("sbs_fraction", "sbs_count")
    return(sigs.lst)
}


#' @name read_meta_data_json
#' @title read_meta_data_json
#' @description
#' Reads in meta data jsons and gets signatures that are lists into the correct list format in the data.table. Used in create_distributions
#'
#' @param meta_json path to metadata.json for case reports
#' @param patient_id sample name for metadata_json
#'
#' @return data.table of the metadata.json for each sample
#' @export
#' @author Stanley Clarke
read_meta_data_json <- function(meta_json, patient_id) {
    json.dt <- as.data.table(jsonlite::read_json(meta_json, simplifyVector = TRUE)) ## have to change from this to work with signatures
    keep_cols <- names(json.dt) %>% grep("signatures|deletionInsertion|sigprofiler|deconstructsigs", ., value = TRUE, invert = TRUE)
    json.dt <- json.dt[, ..keep_cols]
    json.lst <- jsonlite::read_json(meta_json, simplifyVector = FALSE)
    if (!is.null(json.lst[[1]]$signatures)) {
        json.dt$signatures <- list(json.lst[[1]]$signatures)
    }
    if (!is.null(json.lst[[1]]$deletionInsertion)) {
        json.dt$deletionInsertion <- list(json.lst[[1]]$deletionInsertion)
    }
    if (!is.null(json.lst[[1]]$deconstructsigs_sbs_fraction)) {
        json.dt$deconstructsigs_sbs_fraction <- list(json.lst[[1]]$deconstructsigs_sbs_fraction)
    }
    if (!is.null(json.lst[[1]]$sigprofiler_indel_fraction)) {
        json.dt$sigprofiler_indel_fraction <- list(json.lst[[1]]$sigprofiler_indel_fraction)
    }
    if (!is.null(json.lst[[1]]$sigprofiler_indel_count)) {
        json.dt$sigprofiler_indel_count <- list(json.lst[[1]]$sigprofiler_indel_count)
    }
    if (!is.null(json.lst[[1]]$sigprofiler_sbs_fraction)) {
        json.dt$sigprofiler_sbs_fraction <- list(json.lst[[1]]$sigprofiler_sbs_fraction)
    }
    if (!is.null(json.lst[[1]]$sigprofiler_sbs_count)) {
        json.dt$sigprofiler_sbs_count <- list(json.lst[[1]]$sigprofiler_sbs_count)
    }
    return(json.dt)
}

#' @name convert_signature_meta_json
#' @title convert_signature_meta_json
#' @description
#' converts signatures from the list objects in the meta data to data.tables to write jsons for distributions. Used in create_distributions
#'
#' @param jsons.dt json dt after using read_meta_data_json
#' @param cores cores for generating each signature
#' @return list object with all of the signature distribution data.tables
#' @export
#' @author Stanley Clarke
convert_signature_meta_json <- function(jsons.dt, cores = 1) {
    ## signatures-sbs deconstruct sigs ## in column signatures but also in deconstructsigs_sbs_fraction
    if ("signatures" %in% names(jsons.dt)) {
        jsons.dt2 <- copy(jsons.dt)
        jsons.dt2[, length_signatures := sapply(signatures, function(x) length(unlist(x)))]
        jsons.dt2 <- jsons.dt2[length_signatures != 0, ]
        sigs.lst <- mclapply(1:nrow(jsons.dt2), function(x) {
            jsons.sub.dt2 <- jsons.dt2[x, ]
            sigs.dt <- jsons.sub.dt2$signatures[[1]] %>% as.data.table()
            sigs.dt[, pair := jsons.sub.dt2$pair]
            sigs.dt[, tumor_type := jsons.sub.dt2$tumor_type]
            return(sigs.dt)
        }, mc.cores = cores)
        deconstruct_sigs.dt <- rbindlist(sigs.lst, fill = TRUE)
    }
    ## signatures-indels sigprofiler assignment## in column deletionInsertion and sigprofiler_indel_fraction
    if ("deletionInsertion" %in% names(jsons.dt)) {
        jsons.dt2 <- copy(jsons.dt)
        jsons.dt2[, length_signatures := sapply(deletionInsertion, function(x) length(unlist(x)))]
        jsons.dt2 <- jsons.dt2[length_signatures != 0, ]
        sigs.lst <- mclapply(1:nrow(jsons.dt2), function(x) {
            jsons.sub.dt2 <- jsons.dt2[x, ]
            sigs.dt <- jsons.sub.dt2$deletionInsertion[[1]] %>% as.data.table()
            sigs.dt[, pair := jsons.sub.dt2$pair]
            sigs.dt[, tumor_type := jsons.sub.dt2$tumor_type]
            return(sigs.dt)
        }, mc.cores = cores)
        sigprofilerassignment_indels.dt <- rbindlist(sigs.lst, fill = TRUE)
    }
    ## signatures-indel counts sigprofiler assignment
    if ("sigprofiler_indel_count" %in% names(jsons.dt)) {
        jsons.dt2 <- copy(jsons.dt)
        jsons.dt2[, length_signatures := sapply(sigprofiler_indel_count, function(x) length(unlist(x)))]
        jsons.dt2 <- jsons.dt2[length_signatures != 0, ]
        sigs.lst <- mclapply(1:nrow(jsons.dt2), function(x) {
            jsons.sub.dt2 <- jsons.dt2[x, ]
            sigs.dt <- jsons.sub.dt2$sigprofiler_indel_count[[1]] %>% as.data.table()
            sigs.dt[, pair := jsons.sub.dt2$pair]
            sigs.dt[, tumor_type := jsons.sub.dt2$tumor_type]
            return(sigs.dt)
        }, mc.cores = cores)
        sigprofilerassignment_indels_counts.dt <- rbindlist(sigs.lst, fill = TRUE)
    } else {
        sigprofilerassignment_indels_counts.dt <- NULL
    }
    ## signatures-sigprofiler sbs fraction
    if ("sigprofiler_sbs_fraction" %in% names(jsons.dt)) {
        jsons.dt2 <- copy(jsons.dt)
        jsons.dt2[, length_signatures := sapply(sigprofiler_sbs_fraction, function(x) length(unlist(x)))]
        jsons.dt2 <- jsons.dt2[length_signatures != 0, ]
        sigs.lst <- mclapply(1:nrow(jsons.dt2), function(x) {
            jsons.sub.dt2 <- jsons.dt2[x, ]
            sigs.dt <- jsons.sub.dt2$sigprofiler_sbs_fraction[[1]] %>% as.data.table()
            sigs.dt[, pair := jsons.sub.dt2$pair]
            sigs.dt[, tumor_type := jsons.sub.dt2$tumor_type]
            return(sigs.dt)
        }, mc.cores = cores)
        sigprofilerassignment_sbs.dt <- rbindlist(sigs.lst, fill = TRUE)
    } else {
        sigprofilerassignment_sbs.dt <- NULL
    }

    ## signatures-sigprofiler sbs count
    if ("sigprofiler_sbs_count" %in% names(jsons.dt)) {
        jsons.dt2 <- copy(jsons.dt)
        jsons.dt2[, length_signatures := sapply(sigprofiler_sbs_count, function(x) length(unlist(x)))]
        jsons.dt2 <- jsons.dt2[length_signatures != 0, ]
        sigs.lst <- mclapply(1:nrow(jsons.dt2), function(x) {
            jsons.sub.dt2 <- jsons.dt2[x, ]
            sigs.dt <- jsons.sub.dt2$sigprofiler_sbs_count[[1]] %>% as.data.table()
            sigs.dt[, pair := jsons.sub.dt2$pair]
            sigs.dt[, tumor_type := jsons.sub.dt2$tumor_type]
            return(sigs.dt)
        }, mc.cores = cores)
        sigprofilerassignment_sbs_count.dt <- rbindlist(sigs.lst, fill = TRUE)
    } else {
        sigprofilerassignment_sbs_count.dt <- NULL
    }
    ## create list object to return
    list_return <- list(deconstruct_sigs.dt, sigprofilerassignment_indels.dt, sigprofilerassignment_indels_counts.dt, sigprofilerassignment_sbs.dt, sigprofilerassignment_sbs_count.dt)
    names(list_return) <- c("deconstruct_sigs.dt", "sigprofilerassignment_indels.dt", "sigprofilerassignment_indels_counts.dt", "sigprofilerassignment_sbs.dt", "sigprofilerassignment_sbs_counts.dt")
    return(list_return)
}
