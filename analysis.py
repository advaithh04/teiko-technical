"""
Runs Parts 2–4 of the analysis and saves outputs to the outputs/ directory.
Requires load_data.py to have been run first (cell_counts.db must exist).
"""

import os
import sqlite3

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

DB_PATH = "cell_counts.db"
OUT_DIR = "outputs"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def get_conn():
    if not os.path.exists(DB_PATH):
        raise FileNotFoundError(f"{DB_PATH} not found — run load_data.py first")
    return sqlite3.connect(DB_PATH)


# ---------------------------------------------------------------------------
# Part 2: Summary table
# ---------------------------------------------------------------------------

def compute_summary_table(conn) -> pd.DataFrame:
    sql = """
        WITH totals AS (
            SELECT sample_id, SUM(count) AS total_count
            FROM cell_counts
            GROUP BY sample_id
        )
        SELECT
            cc.sample_id  AS sample,
            t.total_count,
            cc.population,
            cc.count,
            ROUND(100.0 * cc.count / t.total_count, 4) AS percentage
        FROM cell_counts cc
        JOIN totals t ON cc.sample_id = t.sample_id
        ORDER BY cc.sample_id, cc.population
    """
    return pd.read_sql_query(sql, conn)


# ---------------------------------------------------------------------------
# Part 3: Statistical analysis — melanoma / miraclib / PBMC
# ---------------------------------------------------------------------------

def compute_part3_data(conn) -> pd.DataFrame:
    sql = """
        WITH totals AS (
            SELECT sample_id, SUM(count) AS total_count
            FROM cell_counts
            GROUP BY sample_id
        )
        SELECT
            cc.sample_id,
            sub.response,
            sub.subject_id,
            cc.population,
            cc.count,
            t.total_count,
            ROUND(100.0 * cc.count / t.total_count, 4) AS percentage
        FROM cell_counts cc
        JOIN samples  s   ON cc.sample_id  = s.sample_id
        JOIN subjects sub ON s.subject_id  = sub.subject_id
        JOIN totals   t   ON cc.sample_id  = t.sample_id
        WHERE sub.condition   = 'melanoma'
          AND sub.treatment   = 'miraclib'
          AND s.sample_type   = 'PBMC'
    """
    return pd.read_sql_query(sql, conn)


def run_statistics(df: pd.DataFrame) -> pd.DataFrame:
    results = []
    for pop in POPULATIONS:
        pop_df = df[df["population"] == pop]
        resp = pop_df[pop_df["response"] == "yes"]["percentage"].values
        non_resp = pop_df[pop_df["response"] == "no"]["percentage"].values
        if len(resp) < 2 or len(non_resp) < 2:
            continue
        stat, pval = stats.mannwhitneyu(resp, non_resp, alternative="two-sided")
        results.append(
            {
                "population": pop,
                "n_responders": len(resp),
                "n_non_responders": len(non_resp),
                "mean_pct_responders": round(resp.mean(), 4),
                "mean_pct_non_responders": round(non_resp.mean(), 4),
                "median_pct_responders": round(float(pd.Series(resp).median()), 4),
                "median_pct_non_responders": round(float(pd.Series(non_resp).median()), 4),
                "mann_whitney_u": round(stat, 4),
                "p_value": round(pval, 6),
                "significant": pval < 0.05,
            }
        )
    return pd.DataFrame(results)


def make_boxplot(df: pd.DataFrame, out_path: str):
    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(18, 6), sharey=False)
    fig.suptitle(
        "Melanoma / Miraclib / PBMC: Cell Population Frequencies\nResponders vs Non-Responders",
        fontsize=13,
        fontweight="bold",
    )

    colors = {"yes": "#4c9be8", "no": "#e8744c"}
    labels = {"yes": "Responder", "no": "Non-Responder"}

    for ax, pop in zip(axes, POPULATIONS):
        pop_df = df[df["population"] == pop]
        groups = [
            pop_df[pop_df["response"] == "yes"]["percentage"].values,
            pop_df[pop_df["response"] == "no"]["percentage"].values,
        ]
        bp = ax.boxplot(
            groups,
            patch_artist=True,
            widths=0.5,
            medianprops=dict(color="black", linewidth=2),
        )
        for patch, key in zip(bp["boxes"], ["yes", "no"]):
            patch.set_facecolor(colors[key])
            patch.set_alpha(0.75)

        # overlay individual data points
        for i, (grp, key) in enumerate(zip(groups, ["yes", "no"]), start=1):
            ax.scatter(
                [i] * len(grp),
                grp,
                color=colors[key],
                edgecolors="black",
                zorder=5,
                s=30,
                alpha=0.8,
            )

        ax.set_title(pop.replace("_", " ").title(), fontsize=10)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Resp.", "Non-Resp."], fontsize=8)
        ax.set_ylabel("Relative Frequency (%)" if pop == POPULATIONS[0] else "")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Boxplot saved → {out_path}")


# ---------------------------------------------------------------------------
# Part 4: Subset analysis
# ---------------------------------------------------------------------------

def run_part4(conn):
    # Baseline melanoma PBMC miraclib samples
    baseline_sql = """
        SELECT
            s.sample_id,
            sub.subject_id,
            sub.project_id,
            sub.response,
            sub.sex,
            s.sample_type,
            s.time_from_treatment_start
        FROM samples  s
        JOIN subjects sub ON s.subject_id = sub.subject_id
        WHERE sub.condition               = 'melanoma'
          AND s.sample_type               = 'PBMC'
          AND s.time_from_treatment_start = 0
          AND sub.treatment               = 'miraclib'
        ORDER BY sub.project_id, s.sample_id
    """
    baseline = pd.read_sql_query(baseline_sql, conn)

    # Samples per project
    proj_counts = baseline.groupby("project_id")["sample_id"].count().reset_index()
    proj_counts.columns = ["project", "sample_count"]

    # Distinct subjects → responder / non-responder counts
    subj_unique = baseline.drop_duplicates(subset="subject_id")
    response_counts = subj_unique["response"].value_counts().reset_index()
    response_counts.columns = ["response", "subject_count"]

    # Distinct subjects → sex counts
    gender_counts = subj_unique["sex"].value_counts().reset_index()
    gender_counts.columns = ["sex", "subject_count"]

    # Average B cells: melanoma males, all sample/treatment types, time=0, responders
    bcell_sql = """
        SELECT ROUND(AVG(cc.count), 2) AS avg_b_cells
        FROM cell_counts cc
        JOIN samples  s   ON cc.sample_id = s.sample_id
        JOIN subjects sub ON s.subject_id = sub.subject_id
        WHERE sub.condition               = 'melanoma'
          AND sub.sex                     = 'M'
          AND sub.response                = 'yes'
          AND s.time_from_treatment_start = 0
          AND cc.population               = 'b_cell'
    """
    avg_bcell = pd.read_sql_query(bcell_sql, conn)["avg_b_cells"].iloc[0]

    return baseline, proj_counts, response_counts, gender_counts, avg_bcell


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    conn = get_conn()

    print("=== Part 2: Summary Table ===")
    summary = compute_summary_table(conn)
    summary.to_csv(f"{OUT_DIR}/part2_summary_table.csv", index=False)
    print(f"  Rows: {len(summary)}  →  saved to {OUT_DIR}/part2_summary_table.csv")
    print(summary.head(10).to_string(index=False))

    print("\n=== Part 3: Statistical Analysis ===")
    p3 = compute_part3_data(conn)
    stats_df = run_statistics(p3)
    stats_df.to_csv(f"{OUT_DIR}/part3_stats_results.csv", index=False)
    make_boxplot(p3, f"{OUT_DIR}/part3_boxplot.png")
    print("\nStatistical Results (Mann-Whitney U, two-sided):")
    print(stats_df.to_string(index=False))
    sig = stats_df[stats_df["significant"]]["population"].tolist()
    print(f"\nSignificant populations (p < 0.05): {sig if sig else 'None'}")

    print("\n=== Part 4: Subset Analysis ===")
    baseline, proj_counts, response_counts, gender_counts, avg_bcell = run_part4(conn)
    baseline.to_csv(f"{OUT_DIR}/part4_baseline_samples.csv", index=False)
    proj_counts.to_csv(f"{OUT_DIR}/part4_project_counts.csv", index=False)
    response_counts.to_csv(f"{OUT_DIR}/part4_response_counts.csv", index=False)
    gender_counts.to_csv(f"{OUT_DIR}/part4_gender_counts.csv", index=False)

    print(f"\nTotal baseline samples (melanoma/PBMC/miraclib/time=0): {len(baseline)}")
    print("\nSamples per project:")
    print(proj_counts.to_string(index=False))
    print("\nSubjects by response:")
    print(response_counts.to_string(index=False))
    print("\nSubjects by sex:")
    print(gender_counts.to_string(index=False))
    print(f"\nAverage B cells (melanoma males, all types, time=0, responders): {avg_bcell:.2f}")

    with open(f"{OUT_DIR}/part4_bcell_average.txt", "w") as f:
        f.write(f"Average B cells (melanoma males, all sample/treatment types, time=0, responders): {avg_bcell:.2f}\n")

    conn.close()
    print("\nAll outputs saved to outputs/")


if __name__ == "__main__":
    main()
