# Immune Cell Population Analysis — Loblaw Bio

This project analyzes immune cell population data from a clinical trial of miraclib. It covers loading the data into a relational database, computing cell population frequencies, running statistical comparisons between responders and non-responders, and exploring baseline subsets — all wrapped in an interactive dashboard.

---

## Live Dashboard

[https://teiko-bio-analysis.streamlit.app](https://teiko-bio-analysis.streamlit.app)

Or run it locally with `make dashboard` and go to `http://localhost:8501`.

---

## How to Run

Tested on Python 3.10+. Works out of the box in GitHub Codespaces.

```bash
make setup       # installs dependencies from requirements.txt
make pipeline    # builds the database, loads the CSV, runs all analysis
make dashboard   # starts the Streamlit dashboard
```

Run them in that order. No arguments needed anywhere.

---

## File Structure

```
├── cell-count.csv        # the raw input data
├── load_data.py          # sets up the SQLite schema and loads the CSV
├── analysis.py           # runs parts 2–4, saves tables and plots to outputs/
├── dashboard.py          # the Streamlit dashboard
├── requirements.txt
├── Makefile
└── outputs/
    ├── part2_summary_table.csv
    ├── part3_boxplot.png
    ├── part3_stats_results.csv
    ├── part4_baseline_samples.csv
    ├── part4_project_counts.csv
    ├── part4_response_counts.csv
    ├── part4_gender_counts.csv
    └── part4_bcell_average.txt
```

---

## Why I Structured It This Way

I kept three separate scripts instead of one big file because each piece has a different job. `load_data.py` only cares about getting data into the database correctly. `analysis.py` runs all the actual analysis and saves outputs you can share or reference later — PNGs, CSVs, printed results. The dashboard reads directly from the database at runtime so it always reflects the current state, not a stale cached file.

The dashboard auto-initializes the database if it doesn't find one, which is what makes it work on Streamlit Cloud without needing to run the pipeline manually first.

---

## Database Schema

```
projects    (project_id PK)
subjects    (subject_id PK, project_id FK, condition, age, sex, treatment, response)
samples     (sample_id PK, subject_id FK, sample_type, time_from_treatment_start)
cell_counts (id PK, sample_id FK, population TEXT, count INTEGER)
```

### Design decisions

I went with a normalized four-table design. The main thing I wanted to avoid was repeating patient-level information (like treatment and response) across every sample row — those are facts about the subject, not the sample, so they live in the subjects table.

For `cell_counts` I stored the data in long format — one row per (sample, population) pair — instead of five separate columns. This felt like the right call because:

- You can add a new cell population without touching the schema
- Filtering and aggregating by population is just a `WHERE` clause
- The wide format would need unpivoting in code every time you want to compare across populations

I added indexes on `sample_id`, `population`, `condition`, and `treatment` since those are the columns everything filters on.

### How it scales

If the trial grows to hundreds of projects and thousands of samples, the structure holds up. New projects are just new rows in `projects`. New samples are new rows in `samples` and `cell_counts`. The long format in `cell_counts` means the table grows linearly with samples × populations, and the indexes keep aggregation queries fast.

If you needed to add new types of analytics — say, comparing across sample types or adding new metadata fields — you'd just add columns to `subjects` or `samples` rather than restructuring anything.

---

## Results Summary

### Part 2
Built a summary table with relative frequencies for all 10,500 samples across 5 populations (52,500 rows total). Each row has sample id, total count, population name, raw count, and percentage.

### Part 3
Filtered to melanoma patients on miraclib, PBMC samples only. Compared relative frequencies between responders and non-responders using the Mann-Whitney U test (two-sided, α = 0.05). I chose Mann-Whitney because the sample sizes per group aren't huge and I didn't want to assume normality.

**CD4 T cell was the only significant result (p = 0.013)** — responders had slightly higher CD4 T cell frequencies on average. The other four populations didn't reach significance.

### Part 4
Filtered to melanoma PBMC samples at baseline (time = 0) on miraclib:

- 656 total samples — 384 from prj1, 272 from prj3
- 331 responders, 325 non-responders
- 344 male subjects, 312 female subjects
- Average B cell count for melanoma males (all sample/treatment types, time = 0, responders): **10206.15**
