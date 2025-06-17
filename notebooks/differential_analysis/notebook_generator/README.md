# 📘 Notebook Generator + Reverse Sync

This project allows you to **compose project-specific Jupyter notebooks from reusable content blocks**, and safely **sync edits made in notebooks back into those blocks**.

---

## 📂 Directory Structure

```
notebooks/
  ProjectA.ipynb         # Generated notebook
blocks/
  _shared/
    00_intro.ipynb       # A reusable content block
configs/
  my_config.yaml         # Describes how to build notebooks
R/
  build_notebook.R       # Assembles notebooks from blocks
  reverse_sync.R         # Syncs edited notebooks back to blocks
run/
  build.sh
  sync.sh
```

---

## 🛠️ 1. Build a Notebook

To generate notebooks from config-defined blocks:

```bash
./run/build.sh configs/my_config.yaml
```

- Loads block paths defined in `my_config.yaml`
- Assembles notebook cells in order
- Injects metadata (`source_hash`, `block_path`) for reverse sync

You can re-run this any time — only missing or changed blocks will be rebuilt.

---

## 🔄 2. Reverse Sync Notebook Changes

After editing a notebook (e.g. `notebooks/ProjectA.ipynb`), run:

```bash
./run/sync.sh notebooks/ProjectA.ipynb --update
```

This will:

- Scan each cell for `metadata.block_path` and `metadata.source_hash`
- Locate and load the corresponding block file
- Identify the original cell within that block
- Update the block **only if the cell’s content has changed**

### Optional flags:

- `--debug` – Print detailed output for cell matching and changes
- `--force` – Force updates even if block and notebook both changed

---

## 🧠 How Reverse Sync Works

Each cell in a generated notebook contains metadata like:

```json
"metadata": {
  "block_path": "blocks/_shared/00_intro.ipynb",
  "source_hash": "abc123..."
}
```

This allows reverse sync to:

- Identify which block the cell came from
- Match it back to the correct original block cell
- Replace it **only if the content was modified**

No guessing, no unnecessary overwrites.

---

## ✅ Best Practices

- Always build notebooks using `build.sh`
- Use `--update` to apply edits made in notebooks
- Keep block files small and composable
- Avoid manual editing of block files unless intentional
