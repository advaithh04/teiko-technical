# Immune Cell Population Analysis — Loblaw Bio Clinical Trial

A data pipeline and interactive dashboard for analyzing how the drug candidate **miraclib** affects immune cell populations across a clinical trial cohort.

---

## Dashboard

**Live dashboard:** [https://teiko-technical.streamlit.app](https://teiko-technical.streamlit.app)

> To deploy: go to [share.streamlit.io](https://share.streamlit.io), connect GitHub repo `advaithh04/teiko-technical`, set main file to `dashboard.py`, and update this link with the generated URL.

To run locally: `make dashboard`, then open `http://localhost:8501`

---

## Quickstart (GitHub Codespaces or local)

```bash
make setup      # install all Python dependencies
make pipeline   # initialize database, load data, run full analysis → outputs/
make dashboard  # start the interactive Streamlit dashboard
```

No arguments, no manual steps. Run the three commands in order.
**Requirements:** Python 3.10+. All dependencies are installed by `make setup`.

---

## Project Structure

```
teiko_technical/
├── cell-count.csv           # raw input data
├── load_data.py             # Part 1  — schema creation + CSV → SQLite
├── analysis.py              # Parts 2–4 — analysis, saves all outputs
├── dashboard.py             # interactive Streamlit dashboard (3 tabs)
├── requirements.txt
├── Makefile
└── outputs/
    ├── part2_summary_table.csv       # relative frequencies, all samples
    ├── part3_boxplot.png             # boxplot: responders vs non-responders
    ├── part3_stats_results.csv       # Mann-Whitney U results per population
    ├── part4_baseline_samples.csv    # melanoma / miraclib / PBMC / time=0
    ├── part4_project_counts.csv      # samples per project
    ├── part4_response_counts.csv     # responders / non-responders count
    ├── part4_gender_counts.csv       # males / females count
    └── part4_bcell_average.txt       # avg B cell count (Part 4 final answer)
```

---

## Code Design

**`load_data.py`** is self-contained. Running `python load_data.py` from the repo root creates `cell_counts.db` and ingests every row of `cell-count.csv`. No CLI arguments required.

**`analysis.py`** organizes each part (2, 3, 4) as separate functions with a `main()` guard. It queries the database, writes output files to `outputs/`, and prints a summary to stdout. The function structure means `dashboard.py` can reuse query logic without re-running the script.

**`dashboard.py`** queries SQLite directly at runtime using `@st.cache_data` to avoid redundant queries. It does not depend on pre-computed CSVs, so it always reflects the current database state. Three tabs map to the three analysis parts, with Plotly charts and interactive sidebar filters.

This separation keeps each file focused: `load_data.py` owns ingest, `analysis.py` owns reproducible static output, and `dashboard.py` owns interactivity.

---

## Database Schema

```
projects    (project_id PK)
subjects    (subject_id PK, project_id FK, condition, age, sex, treatment, response)
samples     (sample_id PK, subject_id FK, sample_type, time_from_treatment_start)
cell_counts (id PK, sample_id FK, population TEXT, count INTEGER)
```

### Rationale

A normalized four-table design eliminates redundancy and makes each concept independently queryable.

- **`projects`** — one row per project. Cross-project queries join cleanly; adding a new project is one INSERT with no schema change.
- **`subjects`** — patient demographics and clinical attributes stored once. `treatment` and `response` are patient-level facts and should not be duplicated across every sample or time-point.
- **`samples`** — one row per biological sample collection event. Multiple time-point samples per patient are separate rows, keeping the table normalized.
- **`cell_counts`** — stored in **long (tidy) format**: one row per (sample, population) pair rather than one column per population. This is the most important design decision (see below).

### Why long format for `cell_counts`?

| Concern | Wide (one column per population) | Long (one row per population) |
|---|---|---|
| Add a new cell type | `ALTER TABLE` + re-ingest | INSERT new rows, no schema change |
| Query one population | Parse all columns | `WHERE population = 'b_cell'` |
| Aggregate across populations | Unpivot in application code | Native `GROUP BY` |
| Sparse data | Wastes space on NULLs | Stores only what exists |

### Scalability

| Scale | How the schema handles it |
|---|---|
| Hundreds of projects | `project_id` FK keeps project metadata centralized; new projects never touch other tables |
| Thousands of samples | `cell_counts` grows linearly (n_samples × n_populations); indexed on `sample_id` and `population` |
| New analytic requirements | Long format supports any population-level `GROUP BY` without schema changes; new metadata fields are new columns on `subjects` or `samples` |
| New sample types / time-points | `sample_type` and `time_from_treatment_start` on `samples` allow fine-grained filtering without separate tables |

---

## Analysis Results

### Part 2 — Data Overview
52,500 rows (10,500 samples × 5 populations). Each row reports the relative frequency (%) of one cell population within one sample.

### Part 3 — Statistical Analysis
Filter: melanoma patients, miraclib treatment, PBMC samples only.
Test: **Mann-Whitney U** (two-sided, α = 0.05) — non-parametric, robust to non-normal distributions common in small clinical cohorts.

| Population | Mean % (Resp.) | Mean % (Non-Resp.) | p-value | Significant |
|---|---|---|---|---|
| B Cell | 9.80 | 10.00 | 0.056 | No |
| CD8 T Cell | 24.88 | 24.94 | 0.639 | No |
| **CD4 T Cell** | **30.54** | **29.90** | **0.013** | **Yes** |
| NK Cell | 14.84 | 15.07 | 0.121 | No |
| Monocyte | 19.94 | 20.08 | 0.163 | No |

**CD4 T cell frequency is significantly higher in responders (p = 0.013)**, supporting its candidacy as a predictive biomarker for miraclib response in melanoma.

### Part 4 — Subset Analysis
Melanoma PBMC miraclib samples at baseline (time = 0):

| Metric | Value |
|---|---|
| Total baseline samples | 656 |
| Samples — prj1 | 384 |
| Samples — prj3 | 272 |
| Responders (subjects) | 331 |
| Non-responders (subjects) | 325 |
| Male subjects | 344 |
| Female subjects | 312 |
| **Avg B cells — melanoma males, all types, time=0, responders** | **10206.15** |
