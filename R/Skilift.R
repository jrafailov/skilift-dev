source("src/higlass_and_ggraph.R")  # import auxiliary functions for import datatables and granges bigwigs

#' Global Variables
empty_plot_default = data.table(
    patient.id = character(),
    sample = character(),
    type = character(),
    source = character(),
    title = character(),
    visible = logical(),
    figure = character(),
    tag = character(),
    server = character(),
    uuid = character(),
    defaultChartType = character()
)
default_plot_types = vector("list", base::NCOL(empty_plot_default))
names(default_plot_types) = names(empty_plot_default)
jx = 1:NCOL(empty_plot_default)
for (j in jx) {
  na = NA
  class(na) = class(empty_plot_default[[j]])
  default_plot_types[[j]] = na
}


#' @title Skilift Object
#'
#' @description
#' Class representing a PGV/Case-reports database. Contains metadata, plots, and methods
#' for interacting with the database and converting to and from JSON.
#'
#' @section Methods:
#' \code{initialize()} Initialize the Skilift object
#' \code{load_json()} Load data from JSON file into metadata and plots
#' \code{update_datafiles_json()} Update the JSON data files
#' \code{to_datatable()} Convert metadata and plots to a single data table
#' \code{add_plots()} Add new plots to the database
#' \code{remove_plots()} Remove plots from the database
#' \code{validate()} Validate metadata and plot data
#' \code{list_higlass_tilesets()} Upload file to higlass server
#' \code{upload_to_higlass()} Upload file to higlass server
#' \code{delete_from_higlass()} Upload file to higlass server
#' \code{init_pgv()} Download and launch PGV instance
#'
#' @export
Skilift <- R6Class("Skilift",
  private = list(
    #' @field datafiles_json_path (`character(1)`)
    #' Path to the datafiles.json
    datafiles_json_path = NULL
  ),
  public = list(
    #' @field metadata (`data.table`)\cr
    #' Data table containing patient metadata
    metadata = NULL,

    #' @field plots (`data.table`)\cr
    #' Data table containing plot metadata
    plots = NULL,

    #' @field datadir (`character(1)`).
    datadir = NULL,

    #' @field settings (`character(1)`).
    settings = NULL,

    #' @field higlass_metadata (`list`)\cr
    #' Attribute containing HiGlass metadata
    higlass_metadata = NULL,

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param public_dir (`character(1)`)\cr
    #'   Public or build directory path.
    #' @param datafiles_json_path (`character(1)`)\cr
    #'   JSON file path.
    #' @param datadir (`character(1)`)\cr
    #'   Data directory path.
    #' @param settings (`character(1)`)\cr
    #'   Settings object path.
    #' @param higlass_metadata (`list`, optional)\cr
    #'   HiGlass metadata containing endpoint, username, and password.
    initialize = function(
      public_dir = NULL,
      datafiles_json_path = NULL,
      datadir = NULL,
      settings = NULL,
      higlass_metadata = NULL
    ) {
      internal_settings_path = system.file("extdata", "test_data", "settings.json", package = "Skilift")
      if (!is.null(public_dir) && dir.exists(public_dir)) {
        # Check if the specific paths exist within the public_dir
        self$datadir <- file.path(public_dir, "data")
        self$settings <- file.path(public_dir, "settings.json")
        private$datafiles_json_path <- file.path(public_dir, "datafiles.json")
      } 

      # Overwrite the derived paths if provided as arguments
      self$datadir <- if (!is.null(datadir)) datadir else self$datadir
      self$settings <- if (!is.null(settings)) settings else internal_settings_path
      private$datafiles_json_path <- if (!is.null(datafiles_json_path)) datafiles_json_path else private$datafiles_json_path
      
      # Make sure all required paths are non-null 
      if (is.null(self$datadir) || is.null(self$settings) || is.null(private$datafiles_json_path)) {
        print(private$datafiles_json_path)
        stop("datadir, settings, and datafiles_json_path must be provided or derivable from public_dir")
      }

      # Make sure all required paths exist
      if (
        !dir.exists(self$datadir) 
        || !file.exists(self$settings) 
        # || !file.exists(private$datafiles_json_path) # unnecessary if empty json is being written later by default?
        ) {
        # stop("datadir, settings, and datafiles.json must all exist")
        stop("datadir and settings must both exist")
      }
      
      self$plots = Skilift:::empty_plot_default

      # Load the JSON data
      self$load_json(private$datafiles_json_path)
      
      # Handle higlass_metadata
      if (!is.null(higlass_metadata)) {
        self$higlass_metadata <- list(
          endpoint = ifelse(is.null(higlass_metadata$endpoint), "http://10.1.29.225:41800", higlass_metadata$endpoint),
          username = ifelse(is.null(higlass_metadata$username), "admin", higlass_metadata$username),
          password = ifelse(is.null(higlass_metadata$password), "higlass_test", higlass_metadata$password)
        )
      } else {
        self$higlass_metadata <- list(
          endpoint = "http://10.1.29.225:41800",
          username = "admin",
          password = "higlass_test"
        )
      }
    },

    #' @description
    #' Loads patient metadata and plot data from a JSON file.
    #'
    #' @param datafiles_json_path (`character(1)`)\cr
    #'   Datafiles.json file path.
    #' @return NULL.
    load_json = function(datafiles_json_path) {
      if (!file.exists(datafiles_json_path)) {
        # If the JSON file does not exist, create an empty file
        empty_json <- toJSON(list(), auto_unbox = TRUE)
        # Write the JSON file
        write(empty_json, file = datafiles_json_path)
      }

      json_data <- jsonlite::fromJSON(datafiles_json_path)

      if (length(json_data) == 0) {
        # If the JSON file is empty, return an empty data.table with the correct columns
        self$metadata <- data.table(patient.id = character(),
                                    ref = character(),
                                    description = list())
        setnames(self$metadata, old = names(self$metadata), new = c("patient.id", "ref", "description")) # Rename columns

        ## global variable defined above
        self$plots <- Skilift:::empty_plot_default
        self$validate()

        return(invisible())
      }

      # Extract metadata
      metadata_list <- lapply(names(json_data), function(patient_id) {
                                patient_data <- json_data[[patient_id]]
                                ref <- patient_data$reference
                                description <- list(patient_data$description) # Store description as a list

                                data.table(
                                           patient.id = patient_id,
                                           ref = ref,
                                           description = description
                                )
                                    })
      self$metadata <- rbindlist(metadata_list)
      setnames(self$metadata, old = names(self$metadata), new = c("patient.id", "ref", "description")) # Rename columns

      plots_list <- lapply(names(json_data), function(patient_id) {
                             patient_data <- json_data[[patient_id]]
                             patient_plots <- patient_data$plots

                             plot_data  <- copy(patient_plots)

                             data.table(
                                        patient.id = rep(patient_id, nrow(plot_data)),
                                        plot_data
                             )
                                    })

      self$plots <- rbindlist(plots_list, fill = TRUE)
      self$fill_plots_names()

      self$validate()
    },

    #' @description
    #' Update the JSON data files with current metadata and plot data.
    #' @return NULL.
    update_datafiles_json = function() {
      missing_data <- self$validate()
      if (!is.null(missing_data)) {
        stop(warning("Did not update datafiles because of malformed Skilift object."))
      }

      datafiles_json_path <- private$datafiles_json_path
      # Create a backup file with timestamp
      timestamp <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
      backup_file <- paste0(datafiles_json_path, ".", timestamp)
      file.copy(datafiles_json_path, backup_file)

      # Build the JSON structure
      json_data <- list()

      # Iterate over patient.ids
      for (patient_id in self$metadata$patient.id) {
        # Subset metadata for the current patient
        metadata <- self$metadata[self$metadata$patient.id == patient_id, ]
        description <- metadata$description[[1]]
        reference <- metadata$ref

        # Subset plots for the current patient
        patient_plots <- self$plots[self$plots$patient.id == patient_id, ]
        plots <- list()
        for (i in seq_len(nrow(patient_plots))) {
          plot_entry <- list(
            "type" = patient_plots$type[i],
            "visible" = patient_plots$visible[i]
          )

          # Get the column names of patient_plots
          cols <- names(patient_plots)

          # Loop over the column names
          for (col in cols) {
            # Use [[ ]] to access the column by name and check if it's NA
            if (col != "patient.id") {
              if (!is.na(patient_plots[[col]][i])) {
                # If it's not NA, add it to plot_entry
                plot_entry[[col]] <- patient_plots[[col]][i]
              }
            }
          }

          plots[[i]] <- plot_entry
        }

        # Create patient data entry
        patient_data <- list(
          "description" = description,
          "reference" = reference,
          "plots" = plots
        )

        # Add the patient data to the main json_data list
        json_data[[patient_id]] <- patient_data
      }

      # Convert JSON to string and write to file
      json_string <- jsonlite::toJSON(json_data, auto_unbox = TRUE, pretty = TRUE)
      write(json_string, file = datafiles_json_path)
    },

    #' @description
    #' Convert metadata and plots to a single data table.
    #'
    #' @param filter (`list()`)\cr
    #'  Filter to apply. Consists of two elements: the column to filter on and the string to filter with.
    #'
    #' @return (`data.table`).
    to_datatable = function(filter = NULL) {
      filtered_datatable <- data.table()

      if (!is.null(filter)) {
        column_name <- filter[[1]]
        filter_string <- filter[[2]]

        if (column_name %in% names(self$metadata)) {
          subset_metadata <- self$metadata[grepl(filter_string, self$metadata[[column_name]])]
          filtered_datatable <- merge(subset_metadata, self$plots, by = "patient.id")
        } else if (column_name %in% names(self$plots)) {
          subset_plots <- self$plots[grepl(filter_string, self$plots[[column_name]])]
          filtered_datatable <- merge(self$metadata, subset_plots, by = "patient.id")
        } else {
          stop(paste("Column", column_name, "not found in Skilift"))
        }
      } else {
        filtered_datatable <- merge(self$metadata, self$plots, by = "patient.id")
      }

      return(filtered_datatable)
    },

    #' @description
    #' Add new plots to the Skilift
    #'
    #' @param plots_to_add (`data.table`)\cr
    #'  Data table of plots to add. Must have a minimal set of columns: patient.id, x. 
    #'  x is either a list of the [server, uuid] or an object (GRanges, gWalk, gGraph, JSON), or a filepath
    #'  If adding a new patient, table must include a ref column.
    #' @param cores (`number(1)`)\cr
    #'  Number of cores to use for parallel execution
    #' @return NULL.
  add_plots = function(plots_to_add, cores=2) {
      new_plots <- data.table::setDT(plots_to_add) # Convert to data.table if required
      if (!("overwrite" %in% tolower(names(new_plots)))) {
        new_plots[, overwrite := FALSE]
      }
      # Check if required columns exist
      required_columns <- c("patient.id", "x")
      empty <- data.table::data.table(
        patient.id = character(),
        x = list(), # either a list of the [server, uuid] or an object (GRanges, gWalk, gGraph, JSON), or a filepath
        stringsAsFactors = FALSE
        )

      if (!all(required_columns %in% names(new_plots))) {
        warning("Required columns not found, creating an empty data.table (called 'empty') with required columns instead...")
        message("(All plots must have an 'x' column containing either list(list(server=, uuid=)) or an object (GRanges, gWalk, gGraph, JSON), or a filepath to a valid RDS data file)")

        return(empty)
      }

      if (any(!new_plots$patient.id %in% self$metadata$patient.id)) {
        if (!("ref" %in% names(new_plots)) | any(is.na(new_plots$ref)) | any(is.null(new_plots$ref))) {
          warning("You are trying to add a new patient without specifying its reference. Make sure all references are non-null")
          message("Creating an empty data.table (called 'empty') with required columns instead...")
          empty$ref <- character()
          return(empty)
        }
      }

      handle_granges <- function(plot) {
        if (!"type" %in% names(plot) || is.na(plot$type) || !(plot$type %in% c("scatterplot", "bigwig"))) {
          warning("Plot type must be specified by the user for GRanges objects and should be either type='scatterplot' or type='bigwig'. This plot will be skipped.")
          return(NULL)
        }
        if (plot$type == "scatterplot") {
          if (!"source" %in% names(plot) || is.na(plot$source)) {
            plot$source <- 'coverage.arrow'
          }
          if (is.null(plot$field)) {
            warning("The column name in GRanges containing scores for scatterplot was not specified, using 'foreground' as default.")
            plot$field <- "foreground" # Use default value
          }
        } else if (plot$type == "bigwig") {
          timestamp <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
          unique_filename <- paste0(timestamp, '.bw')
          plot$source  <- unique_filename
        }
        return(plot)
      }

      handle_gobject = function(gobject, plot) {
        # GRanges
        if (is(gobject, "GRanges")) {
          if (plot$type == "mutations") {
            plot$source = plot$source #' default to mutations.json but allow override via source for mutations
            # plot$x = plot$x # if this is to unlist the object - maybe do that explicitly?
          } else {
            plot <- handle_granges(plot)
          }
        # gGraph
        } else if (is(gobject, "gGraph")) {
          if (is.null(plot$type)) {
            warning("Plot type for gGraph has not been specified, using default type: 'genome'. If plotting an allelic graph, ppfit, or mutations plot, specify type as 'allelic', 'ppfit', or 'mutations' respectively.")
            plot$type <- 'genome'
          }

          if (!"source" %in% names(plot) || is.na(plot$source)) {
            plot$source <- switch(plot$type,
            "genome" = 'genome.json',
            "allelic" = 'allelic.json',
            "ppfit" = 'ppfit.json',
            "mutations" = 'mutations.json',
            'genome.json') # default case
          }
        # gWalk
        } else if (is(gobject, "gWalk")) {
          plot$type <- 'walk'
          if (!"source" %in% names(plot) || is.na(plot$source)) {
            plot$source <- 'walks.json'
          }
        }
        return(plot)
      }

      handle_gtrack = function(gtrack_obj, plot) {
        gtrack_data <- gtrack_obj@data[[1]]
        if (is.null(gtrack_data) || length(gtrack_data) == 0) {
          warning("No data found in gTrack gRanges This plot will be skipped.")
          return(NULL)
        }

        is_scatterplot = gtrack_obj$circles
        is_barplot = gtrack_obj$bars
        is_small_scatterplot = length(gtrack_data) <= 200000
        is_bigwig = !is_small_scatterplot || is_barplot
        is_ggraph = length(gtrack_obj@edges[[1]]) != 0

        if (is_scatterplot) {
          if (!"type" %in% names(plot) || is.na(plot$type) || !(plot$type %in% c("scatterplot", "bigwig"))) {
            plot$type = is_bigwig ? "bigwig" : "scatterplot"
            warning("Plot type not specified for gTrack scatterplot. Using default type based on number of ranges (<200k will produce an scatterplot arrow, otherwise a bigwig): ", plot$type, ". To override, specify type='scatterplot' or type='bigwig'.")
          }
          if (plot$type == "scatterplot") {
            if (!"source" %in% names(plot) || is.na(plot$source)) {
              plot$source <- 'coverage.arrow'
            }
            if (is.null(plot$field)) {
              warning("The column name in GRanges containing scores for scatterplot was not specified, using 'foreground' as default.")
              plot$field <- "foreground" # Use default value
            }
          } else if (plot$type == "bigwig") {
            timestamp <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
            unique_filename <- paste0(timestamp, '.bw')
            plot$source  <- unique_filename
          }
        }

        return(plot)
      }

      source_vector <- rep(NA, nrow(new_plots))
      type_vector <- rep(NA, nrow(new_plots))
      #x_vector <- rep(NA, nrow(new_plots))  # Initialize as a list
      x_vector <- vector("list", nrow(new_plots))  # Initialize as a list

      # Loop through each row of the plots_to_add table
      for (i in seq_len(nrow(new_plots))) {
        plot <- new_plots[i, ]

        if (is.null(plot)) {
          warning("Plot is NULL. Skipping...")
          next
        }

        # Check if patient dir exists, if not, create patient directory
        patient_dir <- file.path(self$datadir, plot$patient.id)
        if (!dir.exists(patient_dir)) {
          warning(paste("Patient directory does not exist. Creating a directory for new patient: ", plot$patient.id))
          dir.create(patient_dir)
        }
        
        is_plot_gtrack <- FALSE
        is_plot_gobject <- FALSE
        is_plot_higlass <- FALSE
        is_plot_json <- FALSE
        is_plot_rds <- FALSE
        if (is(plot$x, "list")) {
          is_plot_gtrack <- is(plot$x[[1]], "gTrack")
          is_plot_gobject <- is(plot$x[[1]], "GRanges") || is(plot$x[[1]], "gGraph") || is(plot$x[[1]], "gWalk")
          is_plot_higlass <- length(plot$x) == 2 && all(sapply(plot$x, is.character))
        } else {
          is_plot_json <- file.exists(plot$x) && tools::file_ext(plot$x) == 'json'
          is_plot_rds <- file.exists(plot$x) && tools::file_ext(plot$x) == 'rds'
        }

        # gTrack
        if (is_plot_gtrack) {
          plot <- handle_gtrack(plot$x[[1]], plot)
          plot$x <- list(plot$x[[1]]@data)
        # GRanges, gGraph, gWalk
        } else if (is_plot_gobject) {
          plot <- handle_gobject(plot$x[[1]], plot)
        # higlass
        } else if (is_plot_higlass) {
          plot$type <- 'bigwig'
          plot$server <- plot$x[[1]][[1]]
          plot$uuid  <- plot$x[[1]][[2]]
          new_plots[i, "server"] <- plot$server
          new_plots[i, "uuid"] <- plot$uuid
        # json
        } else if (is_plot_json) {
          if (!"type" %in% names(plot) || is.na(plot$type)) {
            warning("Type must be specified by the user for JSON files. This plot will be skipped.")
            next
          } else {
            if (!"source" %in% names(plot) || is.na(plot$source)) {
              plot$source  <- basename(plot$x)
            }
          }
        # rds
        } else if (is_plot_rds) {
          rds_object <- readRDS(plot$x)
          if (is(rds_object, "gTrack")) {
            plot <- handle_gtrack(rds_object, plot)
          } else if (is(rds_object, "GRanges") || is(rds_object, "gGraph") || is(rds_object, "gWalk")) {
            plot <- handle_gobject(rds_object, plot)
          }
        } else {
          stop("Invalid file/object. File must be a valid GRanges, gGraph, gWalk, gTrack, JSON, or .rds file. Please check your extension and/or filetype and make sure the file/object exists.")
        }


        # Check if source already exists, if so, increment
        if (!is.null(plot$source)) {
          source_full_path <- file.path(self$datadir, plot$patient.id, plot$source)
          if (plot$overwrite && file.exists(source_full_path)) {
            warning("Existing source for plot will be overwritten.")
          } else if (!plot$overwrite && file.exists(source_full_path)) {
            base_name <- tools::file_path_sans_ext(plot$source)
            ext <- tools::file_ext(plot$source)
            counter <- 2
            while (file.exists(file.path(self$datadir, plot$patient.id, paste0(base_name, counter, ".", ext)))) {
              counter <- counter + 1
            }
            plot$source <- paste0(base_name, counter, ".", ext)
          }
          if (!is.null(plot$source)) {
            source_vector[i] <- plot$source
          }
        }
        if (!is.null(plot$type)) {
          type_vector[i] <- plot$type
        }
        if (!is.null(plot$x)) {
          xval <- plot$x
          while (is.list(xval)) xval = xval[[1]]
          x_vector[[i]] <- xval
        }
      } # Break the loop into three pieces to parallize the plot creation
      new_plots$source <- source_vector
      new_plots$type <- type_vector
      ## Need to wrap in list due to odd data.table behavior
      ## If new_plots is nrow == 1, directly assigning will unlist x
      ## This will work regardless of list length
      if(class(x_vector[[1]]) != "character")
          new_plots$x <- list(x_vector) # 
      
      print('new_plots')
      print(new_plots)

      # Define the function to create plot files
      create_plot_file <- function(plot, plot_index) {
        tryCatch({
          if (!is.null(plot$source)) {
            plot_file <- file.path(self$datadir, plot$patient.id, plot$source)

            if (plot$type == "bigwig") {

              # field contains the column name in GRanges that corresponds to bigwig scores
              if (is.null(plot$field)) {
                warning("The column name in GRanges containing scores for bigwig was not specified, using 'foreground' as default")
                score_col_name <- "foreground" # Use default value
              } else {
                score_col_name <- plot$field
              }

              # Find ref chrom lengths with matching patient.id
              settings_data <- jsonlite::fromJSON(self$settings)
              chrom_lengths <- as.data.table(settings_data$coordinates$sets[[plot$ref]])[,.(chromosome,startPoint,endPoint)]
              colnames(chrom_lengths) = c("seqnames","start","end")
              chrom_lengths[!grepl("chr",seqnames), seqnames := paste0("chr",seqnames)] # weird fix because hg38_chr does not have chr on Y and M

              if (is(plot$x[[1]], "GRanges")) {
                bigwig_grange <- plot$x[[1]]
              } else if (endsWith(plot$x, ".rds")) {
                bigwig_grange <- readRDS(plot$x)
              } else {
                stop(message("Invalid file/object. File must be a valid GRanges or .rds file. Please check your extension and/or filetype."))
              }

              gr2bw(gr = bigwig_grange,
                    output_filepath = plot_file,
                    score_col_name = score_col_name,
                    chrom_lengths = chrom_lengths
              )

              # Call upload_to_higlass function with parameters
              uuid <- self$upload_to_higlass(patient.id = plot$patient.id,
                                             datafile = plot_file,
                                             filetype = "bigwig",
                                             datatype = "vector",
                                             name = plot$source,
                                             coordSystem = plot$ref)

              # Bigwig should be deleted after uploading
              if (plot$overwrite && !is.null(uuid)) {
                warning("Overwrite = TRUE: removing bigwig file after successful upload to higlass")
                file.remove(plot_file)
              }
              return(uuid)
            } else if (grepl("\\.json$", plot$x)) {
                file.copy(plot$x, plot_file)
            } else {
                # Use the plot type to determine the conversion function to use
                switch(plot$type,
                "genome" = create_ggraph_json(plot, self$datadir, self$settings),
                "allelic" = create_allelic_json(plot, self$datadir, self$settings),
                "ppfit" = create_ppfit_genome_json(plot, self$datadir, self$settings),
                "scatterplot" = create_cov_arrow(plot, self$datadir, self$settings),
                "walk" = create_gwalk_json(plot, self$datadir, self$settings),
                "mutations" = if(plot$subtype == "somatic"){
                                create_somatic_json(plot, self$datadir, self$settings)
                              } else if (plot$subtype == "germline") {
                                create_germline_json(plot, self$datadir, self$settings)
                              }
              )
            }
          }
        }, error = function(e) {
            message("Error in creating plot file: ", e$message)
            print(e)
            traceback()
          })
      }

      #debug(create_plot_file)
      if (!any(is.null(plot$source))) {
                                        # Use mclapply to create the plot files in parallel
        new_plots <- parallel::mclapply(
          seq_len(nrow(new_plots)), 
          function(i) {
            plot <- new_plots[i, ]
            is_x_empty = is.list(plot$x) && base::lengths(plot$x) == 0
            source_file_path = file.path(self$datadir, plot$patient.id, plot$source)
            # don't do anything if plot already exists
            if (is_x_empty && file.exists(source_file_path)) {
              return(plot)
            }
            if (plot$type == "bigwig") {
              uuid <- create_plot_file(plot)
              if (!(is.na(uuid) || is.null(uuid))) {
                plot$server <- gsub("/$", "", self$higlass_metadata$endpoint)
                plot$uuid <- uuid
                plot$source <- NULL
              } else {
                message("No uuid for plot, skipping adding to plots table")
              }
            } else {
              create_plot_file(plot)
            }
            return(plot)  # Return the modified plot
        }, mc.cores = cores)

        # Assign the updated new_plots object
        new_plots <- rbindlist(new_plots, fill=TRUE)
      }
      # change allelic & mutations to genome for type to render in pgv
      new_plots[type %in% c("allelic", "mutations", "ppfit"), type := "genome"]
      # Last piece of the loop
      for (i in seq_len(nrow(new_plots))) {
        plot <- new_plots[i, ]
        if (!is.null(plot$source) && !is.na(plot$source)) {
          source_full_path  <- file.path(self$datadir, plot$patient.id, plot$source)
          if (!file.exists(source_full_path)) {
            warning("File does not exist: ", source_full_path, " skipping adding to Skilift")
            next
          }
        }

        common_columns <- intersect(names(plot), names(self$plots))
        extended_data <- plot[, ..common_columns]

        is_duplicate <- duplicated(rbind(self$plots, extended_data, fill = TRUE))

        if (any(is_duplicate)) {
          warning("Row already in plots, skipping (plot file will still be overwritten if overwrite flag was set): ")
          print(plot)
        } else {
          self$plots <- rbind(self$plots, extended_data, fill = TRUE)
        }

        if (plot$patient.id %in% self$metadata$patient.id) {
          # If patient.id exists, update description and ref
          if ("description" %in% colnames(plot)) {
            self$metadata[self$metadata$patient.id == plot$patient.id, "description"] <- plot$description
          }
          if ("ref" %in% colnames(plot) & !is.null(plot$ref) & !is.na(plot$ref) ) {
            self$metadata[self$metadata$patient.id == plot$patient.id, "ref"] <- plot$ref
          }
        } else {
          if ("description" %in% colnames(plot)) {
            self$metadata <- rbind(self$metadata, data.frame(patient.id = plot$patient.id, description = plot$description, ref = plot$ref), fill = TRUE)
          } else {
            self$metadata <- rbind(self$metadata, data.frame(patient.id = plot$patient.id, ref = plot$ref), fill = TRUE)
          }
        }
      }

      
      plot_colnames = names(self$plots)
      default_plot_types = Skilift:::default_plot_types # global variable in namespace
      jx = NROW(default_plot_types)
      for (j in jx) {
        nm = names(default_plot_types[j])
        if (!nm %in% plot_colnames) {
          # message("adding missing plot column name: ", nm)
          self$plots[[nm]] = default_plot_types[[j]]
        }

          
      }

      self$update_datafiles_json()
    },

    fill_plots_names = function() {
      plot_colnames = names(self$plots)
      default_plot_types = Skilift:::default_plot_types # global variable in namespace
      jx = seq_len(NROW(default_plot_types))
      for (j in jx) {
        nm = names(default_plot_types[j])
        if (!nm %in% plot_colnames) {
          # message("adding missing plot column name: ", nm)
          self$plots[[nm]] = default_plot_types[[j]]
        }          
      }
    }
    ,

    add_plots__DEV = function(plots_to_add, cores=2) {
      new_plots <- data.table::setDT(plots_to_add) # Convert to data.table if required
      if (!("overwrite" %in% tolower(names(new_plots)))) {
        new_plots[, overwrite := FALSE]
      }
      # Check if required columns exist
      required_columns <- c("patient.id", "x")
      empty <- data.table::data.table(
        patient.id = character(),
        x = list(), # either a list of the [server, uuid] or an object (GRanges, gWalk, gGraph, JSON), or a filepath
        stringsAsFactors = FALSE
        )
      
      # intializing empty defaults
      if (is.null(new_plots$source)) new_plots$source = NA_character_
      if (is.null(new_plots$type)) new_plots$type = NA_character_
      ## Need to wrap in list due to odd data.table behavior
      ## If new_plots is nrow == 1, directly assigning will unlist x
      ## This will work regardless of list length
      if (is.null(new_plots$x)) new_plots$x = list(x_vector)
      new_plots$is_valid = TRUE ## assume valid until proven not 


      if (!all(required_columns %in% names(new_plots))) {
        warning("Required columns not found, creating an empty data.table (called 'empty') with required columns instead...")
        message("(All plots must have an 'x' column containing either list(list(server=, uuid=)) or an object (GRanges, gWalk, gGraph, JSON), or a filepath to a valid RDS data file)")

        return(empty)
      }

      if (any(!new_plots$patient.id %in% self$metadata$patient.id)) {
        if (!("ref" %in% names(new_plots)) | any(is.na(new_plots$ref)) | any(is.null(new_plots$ref))) {
          warning("You are trying to add a new patient without specifying its reference. Make sure all references are non-null")
          message("Creating an empty data.table (called 'empty') with required columns instead...")
          empty$ref <- character()
          return(empty)
        }
      }

      handle_granges <- function(plot) {
        if (!"type" %in% names(plot) || is.na(plot$type) || !(plot$type %in% c("scatterplot", "bigwig"))) {
          warning("Plot type must be specified by the user for GRanges objects and should be either type='scatterplot' or type='bigwig'. This plot will be skipped.")
          return(NULL)
        }
        if (plot$type == "scatterplot") {
          if (!"source" %in% names(plot) || is.na(plot$source)) {
            plot$source <- 'coverage.arrow'
          }
          if (is.null(plot$field)) {
            warning("The column name in GRanges containing scores for scatterplot was not specified, using 'foreground' as default.")
            plot$field <- "foreground" # Use default value
          }
        } else if (plot$type == "bigwig") {
          timestamp <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
          unique_filename <- paste0(timestamp, '.bw')
          plot$source  <- unique_filename
        }
        return(plot)
      }

      handle_gobject = function(gobject, plot) {
        # GRanges
        if (is(gobject, "GRanges")) {
          if (plot$type == "mutations") {
            plot$source = plot$source #' default to mutations.json but allow override via source for mutations
            # plot$x = plot$x # if this is to unlist the object - maybe do that explicitly?
          } else {
            plot <- handle_granges(plot)
          }
        # gGraph
        } else if (is(gobject, "gGraph")) {
          if (is.null(plot$type)) {
            warning("Plot type for gGraph has not been specified, using default type: 'genome'. If plotting an allelic graph, ppfit, or mutations plot, specify type as 'allelic', 'ppfit', or 'mutations' respectively.")
            plot$type <- 'genome'
          }

          if (!"source" %in% names(plot) || is.na(plot$source)) {
            plot$source <- switch(plot$type,
            "genome" = 'genome.json',
            "allelic" = 'allelic.json',
            "ppfit" = 'ppfit.json',
            "mutations" = 'mutations.json',
            'genome.json') # default case
          }
        # gWalk
        } else if (is(gobject, "gWalk")) {
          plot$type <- 'walk'
          if (!"source" %in% names(plot) || is.na(plot$source)) {
            plot$source <- 'walks.json'
          }
        }
        return(plot)
      }

      handle_gtrack = function(gtrack_obj, plot) {
        gtrack_data <- gtrack_obj@data[[1]]
        if (is.null(gtrack_data) || length(gtrack_data) == 0) {
          warning("No data found in gTrack gRanges This plot will be skipped.")
          return(NULL)
        }

        is_scatterplot = gtrack_obj$circles
        is_barplot = gtrack_obj$bars
        is_small_scatterplot = length(gtrack_data) <= 200000
        is_bigwig = !is_small_scatterplot || is_barplot
        is_ggraph = length(gtrack_obj@edges[[1]]) != 0

        if (is_scatterplot) {
          if (!"type" %in% names(plot) || is.na(plot$type) || !(plot$type %in% c("scatterplot", "bigwig"))) {
            plot$type = is_bigwig ? "bigwig" : "scatterplot"
            warning("Plot type not specified for gTrack scatterplot. Using default type based on number of ranges (<200k will produce an scatterplot arrow, otherwise a bigwig): ", plot$type, ". To override, specify type='scatterplot' or type='bigwig'.")
          }
          if (plot$type == "scatterplot") {
            if (!"source" %in% names(plot) || is.na(plot$source)) {
              plot$source <- 'coverage.arrow'
            }
            if (is.null(plot$field)) {
              warning("The column name in GRanges containing scores for scatterplot was not specified, using 'foreground' as default.")
              plot$field <- "foreground" # Use default value
            }
          } else if (plot$type == "bigwig") {
            timestamp <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
            unique_filename <- paste0(timestamp, '.bw')
            plot$source  <- unique_filename
          }
        }

        return(plot)
      }

      # Define the function to create plot files
      create_plot_file <- function(plot, plot_index) {
        tryCatch({
          if (!is.null(plot$source)) {
            plot_file <- file.path(self$datadir, plot$patient.id, plot$source)

            if (plot$type == "bigwig") {

              # field contains the column name in GRanges that corresponds to bigwig scores
              if (is.null(plot$field)) {
                warning("The column name in GRanges containing scores for bigwig was not specified, using 'foreground' as default")
                score_col_name <- "foreground" # Use default value
              } else {
                score_col_name <- plot$field
              }

              # Find ref chrom lengths with matching patient.id
              settings_data <- jsonlite::fromJSON(self$settings)
              chrom_lengths <- as.data.table(settings_data$coordinates$sets[[plot$ref]])[,.(chromosome,startPoint,endPoint)]
              colnames(chrom_lengths) = c("seqnames","start","end")
              chrom_lengths[!grepl("chr",seqnames), seqnames := paste0("chr",seqnames)] # weird fix because hg38_chr does not have chr on Y and M

              if (is(plot$x[[1]], "GRanges")) {
                bigwig_grange <- plot$x[[1]]
              } else if (endsWith(plot$x, ".rds")) {
                bigwig_grange <- readRDS(plot$x)
              } else {
                stop(message("Invalid file/object. File must be a valid GRanges or .rds file. Please check your extension and/or filetype."))
              }

              gr2bw(gr = bigwig_grange,
                    output_filepath = plot_file,
                    score_col_name = score_col_name,
                    chrom_lengths = chrom_lengths
              )

              # Call upload_to_higlass function with parameters
              uuid <- self$upload_to_higlass(patient.id = plot$patient.id,
                                             datafile = plot_file,
                                             filetype = "bigwig",
                                             datatype = "vector",
                                             name = plot$source,
                                             coordSystem = plot$ref)

              # Bigwig should be deleted after uploading
              if (plot$overwrite && !is.null(uuid)) {
                warning("Overwrite = TRUE: removing bigwig file after successful upload to higlass")
                file.remove(plot_file)
              }
              return(uuid)
            } else if (grepl("\\.json$", plot$x)) {
                file.copy(plot$x, plot_file)
            } else {
                # Use the plot type to determine the conversion function to use
                switch(plot$type,
                "genome" = create_ggraph_json(plot, self$datadir, self$settings),
                "allelic" = create_allelic_json(plot, self$datadir, self$settings),
                "ppfit" = create_ppfit_genome_json(plot, self$datadir, self$settings),
                "scatterplot" = create_cov_arrow(plot, self$datadir, self$settings),
                "walk" = create_gwalk_json(plot, self$datadir, self$settings),
                "mutations" = if(plot$subtype == "somatic"){
                                create_somatic_json(plot, self$datadir, self$settings)
                              } else if (plot$subtype == "germline") {
                                create_germline_json(plot, self$datadir, self$settings)
                              }
              )
            }
          }
        }, error = function(e) {
            message("Error in creating plot file: ", e$message)
            print(e)
            traceback()
          })
      }


      ix = seq_len(NROW(new_plots)) 
      iterate_plot_table = function(i, new_plots) {
        tryCatch({
          plot <- data.table::copy(new_plots[i,])
          patient_dir <- file.path(self$datadir, plot$patient.id)
          is_x_empty = is.list(plot$x) && base::lengths(plot$x) == 0
          source_full_path = NA_character_
          if (!is.na(plot$source))
            source_full_path = file.path(self$datadir, plot$patient.id, plot$source)
          is_file_to_be_directly_uploaded = is_x_empty && file.exists(source_full_path)

          if (!dir.exists(patient_dir)) {
            warning(paste("Patient directory does not exist. Creating a directory for new patient: ", plot$patient.id))
            dir.create(patient_dir)
          }

          is_plot_gtrack <- FALSE
          is_plot_gobject <- FALSE
          is_plot_higlass <- FALSE
          is_plot_json <- FALSE
          is_plot_rds <- FALSE

          if (is(plot$x, "list")) {
            is_plot_gtrack <- is(plot$x[[1]], "gTrack")
            is_plot_gobject <- is(plot$x[[1]], "GRanges") || is(plot$x[[1]], "gGraph") || is(plot$x[[1]], "gWalk")
            is_plot_higlass <- length(plot$x) == 2 && all(sapply(plot$x, is.character))
          } else if (NROW(plot$x) > 0) {
            is_plot_json <- file.exists(plot$x) && tools::file_ext(plot$x) == 'json'
            is_plot_rds <- file.exists(plot$x) && tools::file_ext(plot$x) == 'rds'
          }

          # gTrack
          if (is_plot_gtrack) {
            plot <- handle_gtrack(plot$x[[1]], plot)
            plot$x <- list(plot$x[[1]]@data)
          # GRanges, gGraph, gWalk
          } else if (is_plot_gobject) {
            plot <- handle_gobject(plot$x[[1]], plot)
          # higlass
          } else if (is_plot_higlass) {
            plot$type <- 'bigwig'
            plot$server <- plot$x[[1]][[1]]
            plot$uuid  <- plot$x[[1]][[2]]
          # json
          } else if (is_plot_json) {
            if (!"type" %in% names(plot) || is.na(plot$type)) {
              warning("Type must be specified by the user for JSON files. This plot will be skipped.")
              next
            } else {
              if (!"source" %in% names(plot) || is.na(plot$source)) {
                plot$source  <- basename(plot$x)
              }
            }
          # rds
          } else if (is_plot_rds) {
            rds_object <- readRDS(plot$x)
            if (is(rds_object, "gTrack")) {
              plot <- handle_gtrack(rds_object, plot)
            } else if (is(rds_object, "GRanges") || is(rds_object, "gGraph") || is(rds_object, "gWalk")) {
              plot <- handle_gobject(rds_object, plot)
            }
          } else if (is_file_to_be_directly_uploaded) {
            return(plot)
          } else {
            plot$is_valid = FALSE
            offending = paste("Row", i, plot$patient.id, plot$source)
            warning(offending, ": Invalid file/object. File must be a valid GRanges, gGraph, gWalk, gTrack, JSON, or .rds file. Please check your extension and/or filetype and make sure the file/object exists.")
            return(plot)
          }


          # Check if source already exists, if so, increment
          if (!is.na(plot$source)) {
            if (plot$overwrite && file.exists(source_full_path)) {
              warning("Existing source for plot will be overwritten.")
            } else if (!plot$overwrite && file.exists(source_full_path)) {
              base_name <- tools::file_path_sans_ext(plot$source)
              ext <- tools::file_ext(plot$source)
              counter <- 2
              while (file.exists(file.path(self$datadir, plot$patient.id, paste0(base_name, counter, ".", ext)))) {
                counter <- counter + 1
              }
              plot$source <- paste0(base_name, counter, ".", ext)
            }
          }

          ## create plot file
          if (plot$type == "bigwig") {
            uuid <- create_plot_file(plot)
            if (!(is.na(uuid) || is.null(uuid))) {
              plot$server <- gsub("/$", "", self$higlass_metadata$endpoint)
              plot$uuid <- uuid
              plot$source <- NULL
            } else {
              message("No uuid for plot, skipping adding to plots table")
            }
          } else {
            create_plot_file(plot)
          }
        }, error = function(e) {
          message("Error in plot iteration: ", e$message)
          print(e)
          print(plot)
          traceback()
        })
      }

      lst_new_plots = mclapply(
        X = seq_len(NROW(new_plots)), 
        FUN = iterate_plot_table, 
        new_plots = new_plots,
        mc.cores = cores
      )

      new_plots = data.table::rbindlist(lst_new_plots,fill = TRUE)
      new_plots[type %in% c("allelic", "mutations", "ppfit"), type := "genome"]
      
      is_na_source = is.na(new_plots$source)
      source_full_path  <- ifelse(
        !is_na_source,
        file.path(self$datadir, new_plots$patient.id, new_plots$source),
        NA_character_
      )

      is_source_path_present = file.exists(source_full_path)
      non_existent_paths = source_full_path[!file.exists(source_full_path) & !is_na_source]
      if (NROW(non_existent_paths) > 0) {
        for (path in non_existent_paths) warning("File does not exist: ", path, " skipping adding to Skilift")
      }
      new_plots = new_plots[is_source_path_present, ]

      common_columns <- intersect(names(new_plots), names(self$plots))
      extended_data <- new_plots[, ..common_columns]

      newplots_to_add = rbind(self$plots, extended_data, fill = TRUE)

      is_duplicate_row = duplicated(newplots_to_add)
      if (any(is_duplicate_row)) {
        dup_rows = newplots_to_add[is_duplicate_row,]
        for (i in seq_len(NROW(dup_rows))) {
          warning("Row already in plots, skipping (plot file will still be overwritten if overwrite flag was set): ")
          print(dup_rows[i,])
        }
      }

      self$plots <- newplots_to_add[!is_duplicate_row,]

      if (any(new_plots$patient.id %in% self$metadata$patient.id)) {
        match_ix = base::match(new_plots$patient.id, self$metadata$patient.id)
        is_matched = NROW(match_ix) > 0
        if (is_matched && !is.null(new_plots$description) && any(!is.na(new_plots$description))) {
          self$metadata[match_ix,"description"] = ifelse(
            !is.na(new_plots$description),
            new_plots$description,
            self$metadata[match_ix,"description"]
          )
        }
        if (is_matched && !is.null(new_plots$ref) && any(!is.na(new_plots$ref))) {
          self$metadata[match_ix,"ref"] = ifelse(
            !is.na(new_plots$ref),
            new_plots$ref,
            self$metadata[match_ix,"ref"]
          )
        }
      } 

      is_plot_new_patient = !new_plots$patient.id %in% self$metadata$patient.id
      if (any(is_plot_new_patient)) {
        new_plots_for_metadata = new_plots[is_plot_new_patient,]
        new_metadata = data.frame(
          patient.id = new_plots_for_metadata$patient.id, 
          description = NA_character_, 
          ref = new_plots_for_metadata$ref
        )
        if ("description" %in% colnames(new_plots)) {
          new_metadata$description = new_plots_for_metadata$description
        }
        self$metadata <- rbind(self$metadata, new_metadata, fill = TRUE)
      }
      
      
      self$fill_plots_names()

      self$update_datafiles_json()
    }
    ,

    #' @description
    #' Remove plots from the Skilift
    #'
    #' @param plots_to_remove (`data.table`)\cr 
    #'   Data table of plots to remove
    #' @param delete (`logical(1)`)\cr 
    #'   Delete plot files if TRUE.
    #'
    #' @return NULL
    remove_plots = function(plots_to_remove, delete = FALSE) {
      plots_to_remove <- data.table::setDT(plots_to_remove) # Convert to data.table if required

      # Determine if only patient IDs are provided
      is_remove_patients <- length(names(plots_to_remove)) == 1 && names(plots_to_remove) == "patient.id"
      is_remove_server <- "server" %in% names(plots_to_remove) && any(!is.na(plots_to_remove$server)) && any(!is.na(plots_to_remove$uuid))

      # Check if required columns exist
      required_columns <- if (is_remove_patients) "patient.id" else c("patient.id", "source")
      server_required_columns <- c("patient.id", "server", "uuid")
      if (!(all(required_columns %in% names(plots_to_remove)) || all(server_required_columns %in% names(plots_to_remove)))) {
        return("Error: Required columns (patient.id OR patient.id, source OR patient.id, server, uuid) not found in plots_to_remove data.table.")
      }

      # If delete is TRUE, remove plot files and patient directories if empty
      if (delete) {
        # Loop through each row of the plots_to_remove table
        print("Deleting plots from data directory...")
        for (i in seq_len(nrow(plots_to_remove[!is.na("source")]))) {
          plot <- plots_to_remove[i, ]

          # Select plots/patients to remove
          if (is_remove_patients) {
            patient_plots <- self$plots[patient.id == plot$patient.id, .(patient.id, source)]
          } else {
            patient_plots <- self$plots[patient.id == plot$patient.id & source %in% plot$source, .(patient.id, source)]
          }

          # Loop through each plot file and remove it
          for (j in seq_len(nrow(patient_plots))) {
            plot_file <- file.path(self$datadir, patient_plots$patient.id[j], patient_plots$source[j])
            if (file.exists(plot_file)) {
              file.remove(plot_file)
            }
          }

          # Check if patient directories are empty and remove if so
          unique_patients <- unique(patient_plots$patient.id)
          for (patient in unique_patients) {
            patient_dir <- file.path(self$datadir, patient)
            if (dir.exists(patient_dir) && length(dir(paste0(patient_dir, "/*"))) == 0) {
              message(patient_dir, "has no plots, deleting...")
              unlink(patient_dir)
              self$metadata <- self$metadata[patient.id != patient]
            }
          }
        }
      }

      # Remove plots from Skilift$plots
      print("Removing plots from Skilift")
      initial_rows <- nrow(self$plots)

      if (is_remove_patients) {
        self$plots <- self$plots[!self$plots$patient.id %in% plots_to_remove$patient.id | is.na(self$plots$patient.id) | is.na(plots_to_remove$patient.id), ]
      } else {
        if (!is.null(plots_to_remove$source)) {
          self$plots <- self$plots[
                                   !(self$plots$patient.id %in% plots_to_remove$patient.id & !is.na(self$plots$patient.id) & !is.na(plots_to_remove$patient.id) &
                                     self$plots$source %in% plots_to_remove$source & !is.na(self$plots$source) & !is.na(plots_to_remove$source))
                                   ]
        }
        if (is_remove_server) {
          self$plots <- self$plots[
                                   !(self$plots$patient.id %in% plots_to_remove$patient.id & !is.na(self$plots$patient.id) & !is.na(plots_to_remove$patient.id) &
                                     self$plots$server %in% plots_to_remove$server & !is.na(self$plots$server) & !is.na(plots_to_remove$server) &
                                     self$plots$uuid %in% plots_to_remove$uuid & !is.na(self$plots$uuid) & !is.na(plots_to_remove$uuid))
                                   ]
        }
      }

      final_rows <- nrow(self$plots)
      if (initial_rows == final_rows) {
        warning("No rows were removed. Are you sure the row exists in plots?")
      }
      self$update_datafiles_json()

    },

    #' @description
    #' Validate metadata and plot data.
    #'
    #' @return NULL.
    validate = function() {

      plots_empty = nrow(self$plots) == 0
      metadata_empty = nrow(self$metadata) == 0

      if (!(plots_empty || metadata_empty)) {
        # Check if there are any duplicate columns in self$plots
        if (any(duplicated(colnames(self$plots)))) {
          warning("Duplicate columns found in plots table. Removing duplicate...")
          print(colnames(self$plots))
          self$plots[, which(duplicated(names(self$plots))) := NULL]
        }

        # Check if there are any duplicate columns in self$metadata
        if (any(duplicated(colnames(self$metadata)))) {
          warning("Duplicate columns found in plots table. Removing duplicate...")
          print(colnames(self$metadata))
          self$metadata[, which(duplicated(names(self$metadata))) := NULL]
        }

        # Check if all patients have at least one plot, otherwise remove patient from metadata
        patients_without_plots <- self$metadata[!patient.id %in% unique(self$plots$patient.id), patient.id]
        if (length(patients_without_plots) > 0) {
          self$metadata <- self$metadata[!patient.id %in% patients_without_plots]
        }

        # Check if all source files exist in proper directory and server is not null
        # missing_files <- self$plots[(is.na(server) | is.null(server)) & !file.exists(file.path(self$datadir, patient.id, source)), .(patient.id, source)]
        missing_files <- self$plots[is.na(server) & !file.exists(file.path(self$datadir, patient.id, source)), .(patient.id, source)]
        

        # missing_servers <- self$plots[(is.na(source) | is.null(source)) & (is.na(server) | is.null(server) & (is.na(uuid) | is.null(uuid)))]
        missing_servers <- self$plots[is.na(source) & is.na(server) & is.na(uuid)]

        # All plots must have a patient.id and either a source or server and uuid
        missing_values <- self$plots[
          is.na(patient.id) |
          (
            is.na(source) & 
            is.na(server) &
            is.na(uuid)
          ),
          ,.SDcols = c("patient.id", "source", "server", "uuid")
        ]

        # Construct error message and table with missing source and missing values
        error_message <- ""
        missing_data <- data.table()
        if (nrow(missing_files) > 0) {
          error_message <- paste(error_message, "Missing Files:\n")
          error_message <- paste(error_message, paste(missing_files$patient.id, missing_files$source, "\n", sep = "\t"), collapse = "\n")
          missing_data <- rbind(missing_data, missing_files, fill=TRUE)
        }
        if (nrow(missing_servers) > 0) {
          error_message <- paste(error_message, "Missing Servers (or uuids for the servers):\n")
          error_message <- paste(error_message, paste(missing_servers$patient.id, missing_servers$server, missing_servers$uuid, "\n", sep = "\t"), collapse = "\n")
          missing_data <- rbind(missing_data, missing_servers, fill=TRUE)
        }
        if (nrow(missing_values) > 0) {
          error_message <- paste(error_message, "Missing Values:\n")
          error_message <- paste(error_message, paste(missing_values$patient.id, missing_values$source, "\n", sep = "\t"), collapse = "\n")
          missing_data <- rbind(missing_data, missing_values, fill=TRUE)
        }

        # Coerce visible to be boolean
        self$plots$visible = as.logical(self$plots$visible)

        # Return error message if there are any missing files or values
        if (error_message != "") {
          warning(error_message)
          print("Returning data.table with the invalid rows...")
          return(missing_data)
        }
      } else {
        if (plots_empty) {
          warning("plots table is empty")
        }
        if (metadata_empty) {
          warning("metadata table is empty")
        }
      }
    },


    #' @description 
    #' List tilesets on higlass server
    #'
    #' @return tilesets
    #'  Data.table containing list of tilesets
    list_higlass_tilesets = function(endpoint = self$higlass_metadata$endpoint,
                                    username = self$higlass_metadata$username, 
                                    password = self$higlass_metadata$password) {

      endpoint <- gsub("/$", "", endpoint)
      url <- paste0(endpoint, "/api/v1/tilesets/")
      response <- httr::GET(url,
                            authenticate(username, password, "basic"),
                            encode = "multipart"
      )
      # Convert the JSON response to a data.table
      if (response$status_code == 200) {
        json_content <- httr::content(response, as = "text", encoding="UTF-8")
        parsed_json <- jsonlite::fromJSON(json_content)
        count  <- parsed_json$count
      }

      url <- paste0(endpoint, "/api/v1/tilesets/?limit=", count)
      response <- httr::GET(url,
                            authenticate(username, password, "basic"),
                            encode = "multipart"
      )

      # Convert the JSON response to a data.table
      if (response$status_code == 200) {
        json_content <- httr::content(response, as = "text", encoding="UTF-8")
        parsed_json <- jsonlite::fromJSON(json_content)
        results <- parsed_json$results

        dt_tilesets <- as.data.table(results)
        return(dt_tilesets)
      } else {
        stop("Failed to retrieve tilesets from the server.")
      }
    },


    #' @description
    #' Upload file to higlass server
    #' 
    #' @return uuid
    #'  The uuid returned by the response from higlass
    upload_to_higlass = function(endpoint = self$higlass_metadata$endpoint, 
                                 patient.id = "TEST_HIGLASS",
                                 datafile, 
                                 filetype, 
                                 datatype, 
                                 coordSystem, 
                                 name, 
                                 username = self$higlass_metadata$username, 
                                 password = self$higlass_metadata$password) {
      # Define the API endpoint
      print(paste("Uploading datafile to higlass:", datafile))

      # Convert datafile argument to a multipart object
      datafile <- httr::upload_file(datafile)

      # Package other data as key-value pairs in a list
      body <- list(
                   'datafile' = datafile,
                   'filetype' = filetype,
                   'datatype' = datatype,
                   'coordSystem' = coordSystem,
                   'name' = name
      )

      endpoint <- gsub("/$", "", endpoint)
      url <- paste0(endpoint, "/api/v1/tilesets/")

      # Create the response object
      response <- httr::POST(
                             url,
                             authenticate(username, password, "basic"),
                             body = body,
                             encode = "multipart"
      )

      # Parse the response
      response_content <- httr::content(response, "parsed")

      # Store the UUID
      uuid <- response_content$uuid
      filetype  <- response_content$filetype

      print(paste("UUID:", uuid))
      print(paste("filetype:", filetype))
      return(uuid)
    },

    #' @description
    #' Remove files in higlass server
    #'
    #' @return httr:response
    delete_from_higlass = function(
                                   patient.id = "TEST_ADD",
                                   endpoint = self$higlass_metadata$endpoint,
                                   username = self$higlass_metadata$username, 
                                   password = self$higlass_metadata$password,
                                   uuids,
                                   cores = 2
                                   ) {

      endpoint <- gsub("/$", "", endpoint)
      for (uuid in uuids) {
        # Define the API endpoint
        url <- paste0(endpoint, "/api/v1/tilesets/", uuid, "/")

        # Create the response object
        response <- DELETE(url, authenticate(username, password, "basic"))

        response_content <- httr::content(response, "parsed")
        print(response_content)

        remove_higlass <- data.table(patient.id = patient.id, server = endpoint, uuid = uuid)
        self$remove_plots(remove_higlass)
      }
    },

    #' @description
    #' Download and instantiate a PGV instance with symlinked data
    #'
    #' @param pgv_dir (`character(1)`)\cr 
    #'   Directory where the pgv instance will be installed
    #' @param build (`logical(1)`)\cr 
    #'   Flag to indicate whether to build pgv or launch local instance
    #'
    #' @return NULL
    init_pgv = function(pgv_dir, build = FALSE) {
      init_script_path  <- system.file("src", "init_pgv.sh", package="Skilift")
      # Check if node is installed
      cmd <- paste("which", "node")
      is_installed <- length(system(cmd, intern = TRUE)) != 0
      if (!is_installed) {
        stop(paste(program_name, "is not installed. Please install it before proceeding."))
      } else {
        system(paste("bash", init_script_path, private$datafiles_json_path, self$datadir, self$settings, pgv_dir, build))
      }
    }
  )
)
