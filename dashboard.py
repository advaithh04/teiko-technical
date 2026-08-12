"""
Interactive Streamlit dashboard for immune cell population analysis.
Run with: streamlit run dashboard.py
"""

import os
import sqlite3

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats
import streamlit as st

DB_PATH = "cell_counts.db"
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
POP_LABELS = {p: p.replace("_", " ").title() for p in POPULATIONS}


# ---------------------------------------------------------------------------
# DB helpers (cached)
# ---------------------------------------------------------------------------

@st.cache_data
def load_summary_table():
    conn = sqlite3.connect(DB_PATH)
    sql = """
        WITH totals AS (
            SELECT sample_id, SUM(count) AS total_count
            FROM cell_counts GROUP BY sample_id
        )
        SELECT
            cc.sample_id  AS sample,
            sub.project_id,
            sub.condition,
            sub.treatment,
            sub.response,
            sub.sex,
            s.sample_type,
            s.time_from_treatment_start,
            t.total_count,
            cc.population,
            cc.count,
            ROUND(100.0 * cc.count / t.total_count, 4) AS percentage
        FROM cell_counts cc
        JOIN samples  s   ON cc.sample_id = s.sample_id
        JOIN subjects sub ON s.subject_id = sub.subject_id
        JOIN totals   t   ON cc.sample_id = t.sample_id
        ORDER BY cc.sample_id, cc.population
    """
    df = pd.read_sql_query(sql, conn)
    conn.close()
    return df


@st.cache_data
def load_part3_data():
    conn = sqlite3.connect(DB_PATH)
    sql = """
        WITH totals AS (
            SELECT sample_id, SUM(count) AS total_count
            FROM cell_counts GROUP BY sample_id
        )
        SELECT
            cc.sample_id,
            sub.subject_id,
            sub.response,
            cc.population,
            cc.count,
            t.total_count,
            ROUND(100.0 * cc.count / t.total_count, 4) AS percentage
        FROM cell_counts cc
        JOIN samples  s   ON cc.sample_id = s.sample_id
        JOIN subjects sub ON s.subject_id = sub.subject_id
        JOIN totals   t   ON cc.sample_id = t.sample_id
        WHERE sub.condition = 'melanoma'
          AND sub.treatment = 'miraclib'
          AND s.sample_type = 'PBMC'
    """
    df = pd.read_sql_query(sql, conn)
    conn.close()
    return df


@st.cache_data
def load_part4_data():
    conn = sqlite3.connect(DB_PATH)

    baseline_sql = """
        SELECT s.sample_id, sub.subject_id, sub.project_id,
               sub.response, sub.sex, s.sample_type, s.time_from_treatment_start
        FROM samples  s
        JOIN subjects sub ON s.subject_id = sub.subject_id
        WHERE sub.condition = 'melanoma'
          AND s.sample_type = 'PBMC'
          AND s.time_from_treatment_start = 0
          AND sub.treatment = 'miraclib'
        ORDER BY sub.project_id, s.sample_id
    """
    baseline = pd.read_sql_query(baseline_sql, conn)

    bcell_sql = """
        SELECT ROUND(AVG(cc.count), 2) AS avg_b_cells
        FROM cell_counts cc
        JOIN samples  s   ON cc.sample_id = s.sample_id
        JOIN subjects sub ON s.subject_id = sub.subject_id
        WHERE sub.condition = 'melanoma'
          AND sub.sex       = 'M'
          AND sub.response  = 'yes'
          AND s.time_from_treatment_start = 0
          AND cc.population = 'b_cell'
    """
    avg_bcell = pd.read_sql_query(bcell_sql, conn)["avg_b_cells"].iloc[0]
    conn.close()
    return baseline, float(avg_bcell)


def run_stats(df):
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
                "Population": POP_LABELS[pop],
                "N Responders": len(resp),
                "N Non-Responders": len(non_resp),
                "Mean % (Resp.)": round(resp.mean(), 2),
                "Mean % (Non-Resp.)": round(non_resp.mean(), 2),
                "Median % (Resp.)": round(float(pd.Series(resp).median()), 2),
                "Median % (Non-Resp.)": round(float(pd.Series(non_resp).median()), 2),
                "Mann-Whitney U": round(stat, 2),
                "p-value": round(pval, 6),
                "Significant": "✅" if pval < 0.05 else "❌",
            }
        )
    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Pages
# ---------------------------------------------------------------------------

def page_overview():
    st.header("Part 2 — Data Overview: Cell Population Frequencies")
    st.markdown(
        "Relative frequency of each immune cell population per sample. "
        "Use the sidebar filters to explore subsets of the data."
    )

    df = load_summary_table()

    with st.sidebar:
        st.subheader("Filters")
        projects = ["All"] + sorted(df["project_id"].unique().tolist())
        sel_proj = st.selectbox("Project", projects)

        conditions = ["All"] + sorted(df["condition"].unique().tolist())
        sel_cond = st.selectbox("Condition / Indication", conditions)

        treatments = ["All"] + sorted(df["treatment"].unique().tolist())
        sel_treat = st.selectbox("Treatment", treatments)

        sample_types = ["All"] + sorted(df["sample_type"].unique().tolist())
        sel_stype = st.selectbox("Sample Type", sample_types)

        responses = ["All"] + sorted(df["response"].dropna().unique().tolist())
        sel_resp = st.selectbox("Response", responses)

    filtered = df.copy()
    if sel_proj != "All":
        filtered = filtered[filtered["project_id"] == sel_proj]
    if sel_cond != "All":
        filtered = filtered[filtered["condition"] == sel_cond]
    if sel_treat != "All":
        filtered = filtered[filtered["treatment"] == sel_treat]
    if sel_stype != "All":
        filtered = filtered[filtered["sample_type"] == sel_stype]
    if sel_resp != "All":
        filtered = filtered[filtered["response"] == sel_resp]

    display = filtered[["sample", "total_count", "population", "count", "percentage"]].copy()
    display["population"] = display["population"].map(POP_LABELS)

    st.metric("Samples shown", filtered["sample"].nunique())
    st.dataframe(display, use_container_width=True, height=400)

    st.subheader("Distribution of Relative Frequencies by Population")
    fig = px.box(
        filtered,
        x="population",
        y="percentage",
        color="population",
        labels={"population": "Cell Population", "percentage": "Relative Frequency (%)"},
        category_orders={"population": POPULATIONS},
    )
    fig.update_xaxes(
        ticktext=[POP_LABELS[p] for p in POPULATIONS],
        tickvals=POPULATIONS,
    )
    fig.update_layout(showlegend=False)
    st.plotly_chart(fig, use_container_width=True)

    if st.checkbox("Show stacked bar chart (avg % per condition / treatment)"):
        pivot = (
            filtered.groupby(["condition", "treatment", "population"])["percentage"]
            .mean()
            .reset_index()
        )
        fig2 = px.bar(
            pivot,
            x="treatment",
            y="percentage",
            color="population",
            facet_col="condition",
            barmode="stack",
            labels={"percentage": "Mean Relative Frequency (%)", "treatment": "Treatment"},
            title="Average Cell Population Frequencies by Condition & Treatment",
        )
        st.plotly_chart(fig2, use_container_width=True)


def page_stats():
    st.header("Part 3 — Statistical Analysis: Responders vs Non-Responders")
    st.markdown(
        "Melanoma patients treated with **miraclib**, PBMC samples only. "
        "Comparing cell population relative frequencies between responders and non-responders "
        "using the **Mann-Whitney U test** (two-sided, α = 0.05)."
    )

    df = load_part3_data()
    stats_df = run_stats(df)

    # Boxplots
    fig = make_subplots(
        rows=1,
        cols=len(POPULATIONS),
        subplot_titles=[POP_LABELS[p] for p in POPULATIONS],
        shared_yaxes=False,
    )

    colors = {"yes": "#2E86AB", "no": "#E84855"}
    resp_labels = {"yes": "Responder", "no": "Non-Responder"}

    for col_idx, pop in enumerate(POPULATIONS, start=1):
        pop_df = df[df["population"] == pop]
        for resp_key in ["yes", "no"]:
            grp = pop_df[pop_df["response"] == resp_key]
            fig.add_trace(
                go.Box(
                    y=grp["percentage"],
                    name=resp_labels[resp_key],
                    marker_color=colors[resp_key],
                    showlegend=(col_idx == 1),
                    boxpoints="all",
                    jitter=0.3,
                    pointpos=-1.8,
                    legendgroup=resp_key,
                ),
                row=1,
                col=col_idx,
            )

    fig.update_layout(
        height=500,
        title_text="Cell Population Frequencies: Responders vs Non-Responders (Melanoma / Miraclib / PBMC)",
        boxmode="group",
    )
    for i in range(1, len(POPULATIONS) + 1):
        fig.update_yaxes(title_text="Relative Frequency (%)" if i == 1 else "", row=1, col=i)

    st.plotly_chart(fig, use_container_width=True)

    st.subheader("Statistical Results")
    st.dataframe(stats_df, use_container_width=True)

    sig_pops = stats_df[stats_df["Significant"] == "✅"]["Population"].tolist()
    if sig_pops:
        st.success(
            f"**Significant difference (p < 0.05):** {', '.join(sig_pops)}\n\n"
            "These populations show a statistically significant difference in relative frequency "
            "between responders and non-responders, supporting their potential as predictive biomarkers."
        )
    else:
        st.info("No population reached statistical significance (p < 0.05) with the current sample size.")

    st.subheader("Raw Data")
    st.dataframe(
        df[["sample_id", "response", "population", "percentage"]].rename(
            columns={"sample_id": "Sample", "response": "Response",
                     "population": "Population", "percentage": "% Frequency"}
        ),
        use_container_width=True,
        height=300,
    )


def page_subset():
    st.header("Part 4 — Subset Analysis: Baseline Melanoma / Miraclib / PBMC")
    st.markdown(
        "All melanoma PBMC samples at **baseline (time = 0)** from patients treated with **miraclib**."
    )

    baseline, avg_bcell = load_part4_data()

    col1, col2, col3 = st.columns(3)
    col1.metric("Total Baseline Samples", len(baseline))
    col2.metric("Unique Subjects", baseline["subject_id"].nunique())
    col3.metric("Projects", baseline["project_id"].nunique())

    st.subheader("Samples per Project")
    proj_counts = (
        baseline.groupby("project_id")["sample_id"].count().reset_index()
        .rename(columns={"project_id": "Project", "sample_id": "Sample Count"})
    )
    fig_proj = px.bar(
        proj_counts, x="Project", y="Sample Count", color="Project",
        title="Number of Baseline Samples per Project",
    )
    st.plotly_chart(fig_proj, use_container_width=True)
    st.dataframe(proj_counts, use_container_width=True)

    subj_unique = baseline.drop_duplicates(subset="subject_id")

    col_r, col_g = st.columns(2)

    with col_r:
        st.subheader("Subjects by Response")
        resp_counts = subj_unique["response"].value_counts().reset_index()
        resp_counts.columns = ["Response", "Subject Count"]
        fig_resp = px.pie(
            resp_counts, names="Response", values="Subject Count",
            color="Response",
            color_discrete_map={"yes": "#2E86AB", "no": "#E84855"},
            title="Responders vs Non-Responders",
        )
        st.plotly_chart(fig_resp, use_container_width=True)
        st.dataframe(resp_counts, use_container_width=True)

    with col_g:
        st.subheader("Subjects by Sex")
        sex_counts = subj_unique["sex"].value_counts().reset_index()
        sex_counts.columns = ["Sex", "Subject Count"]
        fig_sex = px.pie(
            sex_counts, names="Sex", values="Subject Count",
            color="Sex",
            color_discrete_map={"M": "#6baed6", "F": "#fd8d3c"},
            title="Males vs Females",
        )
        st.plotly_chart(fig_sex, use_container_width=True)
        st.dataframe(sex_counts, use_container_width=True)

    st.subheader("Average B Cells — Melanoma Males, All Sample/Treatment Types, Time = 0, Responders")
    st.metric(
        label="Average B Cell Count (Melanoma, Male, Responder, Time=0, All types)",
        value=f"{avg_bcell:.2f}",
    )
    st.markdown(
        "_Computed across all sample types and treatment types for melanoma male subjects "
        "who responded, at time = 0._"
    )

    st.subheader("Full Baseline Sample Table")
    st.dataframe(baseline, use_container_width=True, height=350)


# ---------------------------------------------------------------------------
# App entry point
# ---------------------------------------------------------------------------

def main():
    st.set_page_config(
        page_title="Immune Cell Population Analysis",
        page_icon="🧬",
        layout="wide",
    )

    if not os.path.exists(DB_PATH):
        with st.spinner("Initializing database — this runs once and takes ~10 seconds..."):
            import load_data
            conn = sqlite3.connect(DB_PATH)
            load_data.create_schema(conn)
            load_data.load_csv(conn, load_data.CSV_PATH)
            conn.close()
        st.cache_data.clear()

    st.title("🧬 Immune Cell Population Analysis Dashboard")
    st.caption("Loblaw Bio Clinical Trial — Miraclib Immune Profiling")

    tabs = st.tabs(["📊 Data Overview", "📈 Statistical Analysis", "🔬 Subset Analysis"])

    with tabs[0]:
        page_overview()

    with tabs[1]:
        page_stats()

    with tabs[2]:
        page_subset()


main()
