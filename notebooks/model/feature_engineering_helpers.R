#' Parse CD-HIT cluster file and create mapping dataframe
#'
#' @param clust_file Path to CD-HIT .clstr file
#' @return Dataframe with Protein_ID and Cluster_CDHit columns
parse_cdhit_clusters <- function(clust_file) {
    lines <- readLines(clust_file)
    protein_names <- c()
    cluster_ids <- c()
    current_cluster <- NULL

    for (line in lines) {
        if (grepl("^>Cluster", line)) {
            current_cluster <- strsplit(line, " ")[[1]][2]
        } else {
            if (!is.null(current_cluster)) {
                protein_name <- sub(".*, >", "", line)
                protein_name <- sub("\\s.*", "", protein_name)
                protein_names <- c(protein_names, protein_name)
                cluster_ids <- c(cluster_ids, current_cluster)
            }
        }
    }

    df_clusters <- data.frame(
        Protein_ID = protein_names,
        Cluster_CDHit = cluster_ids
    ) %>%
        arrange(Protein_ID) %>%
        mutate(Protein_ID = str_replace_all(Protein_ID, "\\...", ""))

    return(df_clusters)
}

#' Load protein sequences from FASTA file
#'
#' @param fasta_file Path to FASTA file
#' @return Dataframe with Protein_ID and sequence columns
load_protein_sequences <- function(fasta_file) {
    sequences <- readAAStringSet(fasta_file)
    seq_ids <- sapply(strsplit(names(sequences), "[ |]"), `[`, 1)

    seq_df <- data.frame(
        Protein_ID = seq_ids,
        sequence = as.character(sequences),
        stringsAsFactors = FALSE
    )

    return(seq_df)
}

#' Add f_ prefix to all feature columns (excluding metadata)
#'
#' @param df Dataframe with mixed feature and metadata columns
#' @param metadata_cols Vector of column names that should NOT get f_ prefix
#' @return Dataframe with renamed feature columns
add_feature_prefix <- function(df, metadata_cols) {
    # Get feature column names (all columns except metadata)
    feature_cols <- setdiff(colnames(df), metadata_cols)

    cat("Adding f_ prefix to", length(feature_cols), "feature columns\n")
    cat("Preserving", length(metadata_cols), "metadata columns\n")

    # Add f_ prefix to all feature columns
    new_colnames <- colnames(df)
    for (col in feature_cols) {
        new_colnames[new_colnames == col] <- paste0('f_', col)
    }
    colnames(df) <- new_colnames

    # Verify
    cat("\nFirst 10 feature columns:\n")
    print(head(grep('^f_', colnames(df), value = TRUE), 10))

    cat("\nMetadata columns:\n")
    print(colnames(df)[!grepl('^f_', colnames(df))])

    return(df)
}

#' Assess feature variance and generate report
#'
#' @param df Dataframe with feature columns (starting with f_)
#' @param output_file Path to save variance report TSV
#' @param plot Whether to generate and display variance plot
#' @return Dataframe with variance report
assess_feature_variance <- function(df, output_file = NULL, plot = TRUE) {
    feature_cols <- grep('^f_', colnames(df), value = TRUE)

    var_values <- sapply(df[feature_cols], function(x) {
        x_clean <- suppressWarnings(as.numeric(as.character(x)))
        if(all(is.na(x_clean))) return(NA)
        return(var(x_clean, na.rm = TRUE))
    })

    var_report <- data.frame(
        Feature = names(var_values),
        Variance = var_values
    ) %>% arrange(desc(Variance))

    # Save report if output file specified
    if (!is.null(output_file)) {
        write.table(var_report, output_file, sep = '\t', quote = FALSE, row.names = FALSE)
        cat("Variance report saved to:", output_file, "\n")
    }

    # Print statistics
    cat("\nVariance Statistics:\n")
    cat("  Mean:", mean(var_values, na.rm = TRUE), "\n")
    cat("  Median:", median(var_values, na.rm = TRUE), "\n")
    cat("  Features with var=0:", sum(var_values == 0, na.rm = TRUE), "\n")
    cat("  Features with var<0.01:", sum(var_values < 0.01, na.rm = TRUE), "\n")

    # Visualize if requested
    if (plot) {
        p <- ggplot(var_report, aes(x = reorder(Feature, Variance), y = Variance)) +
            geom_bar(stat = "identity") +
            coord_flip() +
            labs(title = paste0("Feature Variance (n=", nrow(df), ")"),
                 x = "Feature", y = "Variance") +
            theme_minimal()
        print(p)
    }

    cat("\nNOTE: This is for documentation only. Variance filtering will be applied during modeling.\n")

    return(var_report)
}

#' Report and remove rows with missing values in feature columns
#'
#' @param df Dataframe with feature columns (starting with f_)
#' @return Dataframe with complete cases only
handle_missing_values <- function(df) {
    feature_cols <- grep('^f_', colnames(df), value = TRUE)

    # Count NAs
    na_per_sample <- rowSums(is.na(df[feature_cols]))
    na_per_feature <- colSums(is.na(df[feature_cols]))

    cat("Samples with NAs:", sum(na_per_sample > 0), "\n")

    if (any(na_per_feature > 0)) {
        cat("\nFeatures with NAs:\n")
        print(na_per_feature[na_per_feature > 0])
    }

    # Remove rows with NAs
    df_no_na <- df %>%
        filter(complete.cases(select(., all_of(feature_cols))))

    cat("\nAfter removing rows with NAs:", nrow(df_no_na), "samples\n")
    cat("Removed:", nrow(df) - nrow(df_no_na), "samples\n")

    return(df_no_na)
}

#' Reorder columns (metadata first, then features) and save to file
#'
#' @param df Dataframe with mixed columns
#' @param output_file Path to save output TSV
#' @param metadata_cols Vector of metadata column names (in desired order)
#' @return Reordered dataframe (also saved to file)
save_feature_table <- function(df, output_file, metadata_cols) {
    feature_cols <- grep('^f_', colnames(df), value = TRUE)

    # Reorder: metadata first (except Protein_ID which becomes rowname), then features
    df_out <- df %>%
        column_to_rownames(var = 'Protein_ID') %>%
        select(all_of(c(setdiff(metadata_cols, 'Protein_ID'), feature_cols)))

    # Save
    write.table(df_out, output_file, sep = '\t', quote = FALSE, row.names = TRUE)

    cat("\n=== FILE SAVED ===\n")
    cat("Output file:", output_file, "\n")
    cat("\n=== SUMMARY ===\n")
    cat("  Samples:", nrow(df_out), "\n")
    cat("  Features:", length(feature_cols), "\n")
    cat("  Metadata columns:", length(metadata_cols) - 1, "\n")
    cat("  Total columns:", ncol(df_out), "\n")

    return(df_out)
}
