# elden-ring
**Effelsberg Large-scale Data Exploration with Nextflow for Robust Identification of New Globular cluster pulsars.**

![elden-ring-transformed](https://github.com/user-attachments/assets/3e1a1c35-d055-4266-9ff7-5380ab1d463f)

## Flowchart
![mermaid-diagram-2025-04-24-115529](https://github.com/user-attachments/assets/92b86bb0-0b40-4050-8526-5f380136cf62)

## 📦 Quick Start

Follow these steps to get up and running in minutes.

---

### 1. Clone & make your config

```bash
git clone https://github.com/erc-compact/elden-ring.git
cd elden-ring
```

or 

```bash
nextflow pull erc-compact/elden-ring
```

# Copy the example into your working config
cp example/params.config.example params.config

# Edit `params.config`:
# • Tweak any global parameters (e.g. singularity cacheDir, threads)
```

---

### 2. Prepare your input list

```bash
# Copy the sample CSV, or create your own with the same fields:
cp examples/inputfile.txt .

# By default, the pipeline looks for `inputfile.txt` in the elden-ring folder unless you override it:
# In nextflow.config: params.files_list = 'path/to/inputfile.txt'
```

---

### 3. (Optional) Tune pipeline parameters

Open **`params.config`** and adjust any of:
* **File copy**
You can switch this to true to copy files from other clusters. passwordless login required.

* **RFI filtering**

  ```groovy
  params.generateRfiFilter.run_rfi_filter = true
  params.filtool.run_filtool         = true
  ```
* **Search/fold settings** (e.g. DM range, nbins, etc..)
* **Parfile folding**

  ```groovy
  params.parfold.parfile_path = 'path/to/my.par'
  ```

You can also drop JSON templates into `templates/` for the candidate‐sifting step.

---

### 4. Choose & run an entry workflow

| Workflow              | Description                                                                               |
| --------------------- | ----------------------------------------------------------------------------------------- |
| `full`                | Run end-to-end: RFI mitigation → search & fold → ML classification                        |
| `generate_rfi_filter` | Generate per‐file SK/kurtosis plots (no filtool)                                          |
| `rfi_clean`           | Build SK‐based RFI mask **and** run `filtool` (requires `run_rfi_filter = true`)          |
| `run_search_fold`     | Run search + fold on already‐cleaned data                                                 |
| `fold_par`            | Fold using pre‐computed `.par` files (set `parfold.parfile_path`)                         |
| `candypolice`         | Re‐fold “T1”/“T2” candidates from an existing `candyjar.csv` across all beams & bands     |


Replace `<workflow>` below with any of the above names.

```bash
nextflow run erc-compact/elden-ring \
  -c params.config \
  --files_list inputfile.txt \
  -profile <profile> \
  --entry  <workflow> \
  -resume
```

* **`-profile <profile>`**
  One of the profiles defined in your `nextflow.config` (`local`, `slurm`, `hercules`, etc.).
  If you need a new profile, copy the profile example file from `conf/profiles` and use 
* **`-c <nextflow.config>`**
  Your local config overrides (you can chain `-c` multiple times).
* **`--entry <workflow>`**
  Selects which sub‐workflow to run (default: `full`).
* **`-resume`**
  Re‐uses any existing `work/` cache—skips stages that haven’t changed.

---

### 5. Inspect outputs & reports

* Pipeline outputs land in `basedir` or in `output_dir` of certain workflows as configured in the config file.
* Nextflow’s HTML **report**, **trace**, and **timeline** are generated automatically:

  ```bash
  nextflow run … --with-report --with-trace --with-timeline
  ```
---

❓ **Need help?**
– Open an issue at [https://github.com/erc-compact/elden-ring/issues](https://github.com/erc-compact/elden-ring/issues)
- mail me at fkareem[at]mpifr-bonn.mpg.de

Happy processing! 🚀
