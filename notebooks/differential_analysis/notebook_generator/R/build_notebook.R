# ----- build_notebook.R -----
library(jsonlite)
library(yaml)
library(digest)

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) > 0) args[1] else "configs/default.yaml"
if (!file.exists(config_file)) stop("❌ Config file not found: ", config_file)
cat("📄 Using config file:", config_file, "\n")

config <- yaml::read_yaml(config_file)
block_root <- config$block_root %||% "blocks/"
git_root <- normalizePath(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), winslash = "/", mustWork = TRUE)

load_block_with_metadata <- function(block_path) {
  if (!file.exists(block_path)) return(NULL)
  obj <- tryCatch({ fromJSON(block_path, simplifyVector = FALSE) }, error = function(e) {
    warning("⚠️ Failed to parse JSON: ", block_path)
    return(NULL)
  })

  inject_metadata <- function(cell) {
    cell$metadata <- cell$metadata %||% list()
    cell_source <- paste0(cell$source, collapse = "")
    cell$metadata$source_hash <- digest(cell_source, algo = "sha1")
    abs_path <- normalizePath(block_path, winslash = "/", mustWork = TRUE)
    rel_path <- sub(paste0("^", git_root, "/?"), "", abs_path)
    cell$metadata$block_path <- rel_path
    if (cell$cell_type == "code") {
      cell$execution_count <- cell$execution_count %||% NULL
      cell$outputs <- cell$outputs %||% list()
    }
    return(cell)
  }

  if (!is.null(obj$cells)) {
    return(lapply(obj$cells, inject_metadata))
  } else {
    if (is.null(obj$cell_type)) {
      warning("⚠️ Missing 'cell_type' in block:", block_path)
      return(NULL)
    }
    return(list(inject_metadata(obj)))
  }
}

build_notebook <- function(project) {
  cat("\n📘 Building notebook for:", project$name, "\n")
  cat("Expected blocks:\n")

  found_cells <- list()
  missing_blocks <- list()

  for (rel_path in project$blocks) {
    block_path <- file.path(block_root, rel_path)
    cat(" - ", block_path)
    if (file.exists(block_path)) {
      cat(" ✅\n")
      cells <- load_block_with_metadata(block_path)
      found_cells <- c(found_cells, cells)
    } else {
      cat(" ❌ MISSING\n")
      missing_blocks[[length(missing_blocks) + 1]] <- block_path
    }
  }

  if (length(found_cells) == 0) {
    cat("❗ No blocks found. Skipping notebook creation for:", project$name, "\n")
    return()
  }

  notebook <- list(
    cells = found_cells,
    metadata = list(
      kernelspec = list(display_name = "R", language = "R", name = "ir"),
      language_info = list(
        name = "R",
        codemirror_mode = "r",
        file_extension = ".r",
        mimetype = "text/x-r-source",
        pygments_lexer = "r",
        version = "4.2.2"
      )
    ),
    nbformat = 4,
    nbformat_minor = 2
  )

  block_root_abs <- normalizePath(block_root, winslash = "/", mustWork = TRUE)
  notebook_dir <- normalizePath(file.path(block_root_abs, "..", "notebooks"), winslash = "/", mustWork = FALSE)
  dir.create(notebook_dir, showWarnings = FALSE, recursive = TRUE)
  output_path <- file.path(notebook_dir, paste0(project$name, ".ipynb"))
  write_json(notebook, output_path, auto_unbox = TRUE, pretty = TRUE)

  cat("\n✔️ Notebook written to:", output_path, "\n")
  cat("✅ Cells loaded:", length(found_cells), "\n")
  cat("📦 Blocks processed:", length(project$blocks), "\n")
  cat("❌ Blocks missing:", length(missing_blocks), "\n")
  if (length(missing_blocks) > 0) {
    cat("Missing block files:\n")
    for (miss in missing_blocks) cat(" - ", miss, "\n")
  }
}

for (project in config$projects) {
  build_notebook(project)
}
