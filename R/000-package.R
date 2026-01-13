"_package"

## Global Variables

#' Paired reads factor 
#' 
#' Multiply paired reads by this number for coverage
#' 
#' To get binned coverage into base level coverage,
#' multiply by paired reads factor (2)
PAIRED_READS_FACTOR = 2

#' Read length
#' 
#' Multiply paired reads by this number for read length
#' 
#' To get binned coverage into read length level coverage,
#' multiply by read length (151)
READ_LENGTH = 151

#' Snpeff Protein Coding Annotations
#' 
#' Nuff said.
#' 
#' 
snpeff_protein_coding_annotations = c("frameshift_variant", "stop_lost", "start_lost", "stop_gained", 
"chromosome", "exon_loss_variant", "feature_ablation", "duplication", 
"splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", 
"inframe_insertion", "disruptive_inframe_insertion", "inframe_deletion", 
"disruptive_inframe_deletion", "coding_sequence_variant", "missense_variant", 
"protein_protein_contact", "structural_interaction_variant", 
"rare_amino_acid_variant")


#' QC Flags Tresholds
#'
#' Parse QC flags into strings
#'
#' QC metrics from picard need to be parsed based on coverage, insert size
#' total number of reads, duplicate rate.
#' The strings should be parsable into a form
#' digested in gOS and shown as a single "PASS"/Checkmark", "Warning", or "Fail". 
#' The actual metrics should show up on hover.
qc_flag_thresholds = list(
    list("FAIL", "greater_than_or_equal_to_50x", `<`, 0.99, "Fraction of genome covered at 50X","99%"),
    list("WARN", "purity", `<`, 0.2, "Purity", "20%"),
    list("WARN", "insert_size", `<`, 300, "Insert Size", "300 bp"),
    list("WARN", "percent_duplication", `>`, 0.3, "Duplicate percent", "30%"),
    list("FAIL", "fraction_of_reads_aligned", `<`, 0.9, "Percent of reads aligned", "90%"),
    list("FAIL", "conpair_concordance_metric", `<`, 0.9, "Tumor/Normal SNP concordance", "90%")
)


#' Echtvar annotated vcf fields
#' 
#' Fields from echtvar annotated vcf to keep
echtvar_vcf_fields = c("dbNSFP_gene", "dbNSFP_transcript_id", "dbNSFP_VEP_canonical", 
"dbNSFP_aaref", "dbNSFP_aaalt", "dbNSFP_AlphaMissense_score", 
"dbNSFP_AlphaMissense_pred", "dbNSFP_REVEL_score", "dbNSFP_CADD_phred", 
"dbNSFP_SIFT_pred", "dbNSFP_Polyphen2_HVAR_pred", "dbNSFP_ESM1b_pred", 
"dbNSFP_ESM1b_score", "dbNSFP_RegeneronME_ALL_AF", "dbNSFP_gnomAD211_exomes_controls_AF", 
"dbNSFP_gnomAD211_exomes_controls_POPMAX_AF", "dbNSFP_gnomAD41_joint_AF", 
"dbNSFP_gnomAD41_joint_POPMAX_AF", "dbNSFP_clinvar_clnsig", "dbNSFP_clinvar_trait", 
"dbNSFP_clinvar_review", "dbNSFP_clinvar_hgvs", "clinvar_CLNSIG", 
"clinvar_CLNDN", "clinvar_CLNREVSTAT", "clinvar_ALLELEID", "clinvar_GENEINFO", 
"clinvar_CLNVC", "clinvar_MC", "clinvar_CLNHGVS", "clinvar_ONC", 
"clinvar_ONCDN", "clinvar_ONCCONF", "clinvar_ONCREVSTAT", "clinvar_SCI", 
"clinvar_SCIDN", "clinvar_SCIREVSTAT", "civic_GN", "civic_VT", 
"civic_CSQ")

template_filtered_events_class_icon = list(
    "id" = "alphamissense",
    "title" = "AlphaMissense",
    "dataIndex" = "AlphaMissense",
    "viewType" = "class-icon"
)

template_filtered_events_formatted_number = list(
    "id" =  "altCounts",
    "title" = "components.filtered-events-panel.altCounts",
    "dataIndex" = "altCounts",
    "type" = "numeric",
    "viewType" = "formatted-number",
    "width" = 100L,
    "sortable" = TRUE,
    "rendererProps" = list(
        "format" = ","
    )
)


default_echtvar_json_fields = list(
    list(
        "id" = "alphamissense",
        "title" = "AlphaMissense",
        "dataIndex" = "AlphaMissense",
        "viewType" = "class-icon"
    )
    ,
    list(
        "id" = "sift",
        "title" = "SIFT",
        "dataIndex" = "SIFT",
        "viewType" = "class-icon"
    )
    ,
    list(
        "id" = "polyphen2",
        "title" = "Polyphen2",
        "dataIndex" = "Polyphen2_HVAR",
        "viewType" = "class-icon"
    )
    ,
    list(
        "id" = "clinvar",
        "title" = "Clinvar",
        "dataIndex" = "Clinvar",
        "viewType" = "class-icon"
    )
)

# {
#     "id": "altCounts",
#     "title": "components.filtered-events-panel.altCounts",
#     "dataIndex": "altCounts",
#     "type": "numeric",
#     "viewType": "formatted-number",
#     "width": 100,
#     "sortable": true,
#     "rendererProps": {
#     "format": ","
#     }
# },

# {
#     "id": "complex_events.rDelDup",
#     "title": "rDelDup Events Count",
#     "type": "numeric",
#     "group": "complex-events",
#     "groupTitle": "Complex Events",
#     "groupOrder": 6,
#     "kpiPlot": true,
#     "format": ",",
#     "scale": "linear"
# }

template_metadata = list(
    "id" = "complex_events.rDelDup",
    "title" = "rDelDup Events Count",
    "type" = "numeric",
    "group" = "complex-events",
    "groupTitle" = "Complex Events",
    "groupOrder" = 6,
    "kpiPlot" = TRUE,
    "format" = ",",
    "scale" = "linear"
)
