library(jsonlite)
library(digest)

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
debug <- "--debug" %in% args
args <- setdiff(args, "--debug")

if (length(args) < 1) stop("Usage: Rscript reverse_sync.R <notebook_path> [--update] [--force] [--debug]")

notebook_path <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
do_update <- "--update" %in% args
force_update <- "--force" %in% args

git_root <- normalizePath(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), winslash = "/", mustWork = TRUE)
notebook_rel <- sub(paste0("^", git_root, "/?"), "", notebook_path)

cat("🔄 Reverse syncing from:", notebook_rel, "\n")

notebook <- fromJSON(notebook_path, simplifyVector = FALSE)

# Group notebook cells by block_path
cells_by_block <- split(
  notebook$cells,
  sapply(notebook$cells, function(cell) cell$metadata$block_path %||% NA)
)

updated_blocks <- list()

for (block_path_rel in names(cells_by_block)) {
  if (is.na(block_path_rel)) next

  block_path <- file.path(git_root, block_path_rel)
  if (!file.exists(block_path)) {
    cat("⚠️  Block not found on disk:", block_path_rel, "\n")
    next
  }

  block <- fromJSON(block_path, simplifyVector = FALSE)
  is_multi_cell <- !is.null(block$cells)
  notebook_cells <- cells_by_block[[block_path_rel]]

  if (debug) cat("📦 Block:", block_path_rel, "| Notebook cells:", length(notebook_cells), "\n")

  replaced_indexes <- c()
  block_updated <- FALSE

  if (is_multi_cell) {
    block_cells <- block$cells

    for (nb_cell in notebook_cells) {
      nb_source_hash <- nb_cell$metadata$source_hash %||% NULL
      if (is.null(nb_source_hash)) {
        if (debug) cat("⚠️  Skipping notebook cell with no source_hash\n")
        next
      }

      # Find matching cell in the block
      match_index <- NULL
      for (i in seq_along(block_cells)) {
        original_source <- paste0(block_cells[[i]]$source, collapse = "")
        original_hash <- digest(original_source, algo = "sha1")
        if (original_hash == nb_source_hash) {
          match_index <- i
          break
        }
      }

      if (is.null(match_index)) {
        cat("⚠️  No matching cell found in block for source_hash:", nb_source_hash, "\n")
        next
      }

      # Compare current notebook cell source to block version
      new_source_digest <- digest(paste0(nb_cell$source, collapse = ""), algo = "sha1")
      current_block_digest <- digest(paste0(block_cells[[match_index]]$source, collapse = ""), algo = "sha1")

      if (new_source_digest != current_block_digest) {
        block_cells[[match_index]] <- nb_cell
        replaced_indexes <- c(replaced_indexes, match_index)
        block_updated <- TRUE
        if (debug) cat(sprintf("🔁 Replaced cell %d in %s\n", match_index, block_path_rel))
      } else {
        if (debug) cat(sprintf("✅ Skipped unchanged cell %d in %s\n", match_index, block_path_rel))
      }
    }

    if (block_updated) {
      if (do_update || force_update) {
        write_json(list(cells = block_cells), block_path, auto_unbox = TRUE, pretty = TRUE)
        cat("✅ Updated block:", block_path_rel, "\n")
      } else {
        cat("📝 Detected changes in:", block_path_rel, "(dry run)\n")
      }

      updated_blocks <- c(updated_blocks, block_path_rel)
      if (debug && length(replaced_indexes) > 0) {
        cat("🔍 Replaced cell(s):", paste(replaced_indexes, collapse = ", "), "\n")
      }
    }

  } else {
    # Single-cell block
    nb_cell <- notebook_cells[[1]]
    nb_source_hash <- nb_cell$metadata$source_hash %||% NULL
    if (is.null(nb_source_hash)) {
      if (debug) cat("⚠️  Skipping notebook cell with no source_hash (single-cell block)\n")
      next
    }

    block_source_digest <- digest(paste0(block$source, collapse = ""), algo = "sha1")
    nb_source_digest <- digest(paste0(nb_cell$source, collapse = ""), algo = "sha1")

    if (nb_source_digest != block_source_digest) {
      if (do_update || force_update) {
        write_json(nb_cell, block_path, auto_unbox = TRUE, pretty = TRUE)
        cat("✅ Updated block:", block_path_rel, "\n")
      } else {
        cat("📝 Change detected in:", block_path_rel, "(dry run)\n")
      }
      updated_blocks <- c(updated_blocks, block_path_rel)
    } else {
      if (debug) cat("✅ Skipped unchanged single-cell block:", block_path_rel, "\n")
    }
  }
}

# Summary
cat("\n🎉 Reverse sync complete.\n")
if (length(updated_blocks) > 0) {
  cat("🧱 Updated blocks:\n")
  for (b in updated_blocks) cat(" -", b, "\n")
} else {
  cat("ℹ️ No blocks updated.\n")
}
if (!do_update && !force_update) {
  cat("💡 Use --update to apply changes.\n")
}
