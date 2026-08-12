import sqlite3
import pandas as pd
import os

DB_PATH = "cell_counts.db"
CSV_PATH = "cell-count.csv"


def create_schema(conn):
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS projects (
            project_id TEXT PRIMARY KEY
        );

        CREATE TABLE IF NOT EXISTS subjects (
            subject_id  TEXT PRIMARY KEY,
            project_id  TEXT NOT NULL REFERENCES projects(project_id),
            condition   TEXT,
            age         INTEGER,
            sex         TEXT,
            treatment   TEXT,
            response    TEXT
        );

        CREATE TABLE IF NOT EXISTS samples (
            sample_id                 TEXT PRIMARY KEY,
            subject_id                TEXT NOT NULL REFERENCES subjects(subject_id),
            sample_type               TEXT,
            time_from_treatment_start INTEGER
        );

        CREATE TABLE IF NOT EXISTS cell_counts (
            id         INTEGER PRIMARY KEY AUTOINCREMENT,
            sample_id  TEXT    NOT NULL REFERENCES samples(sample_id),
            population TEXT    NOT NULL,
            count      INTEGER NOT NULL
        );

        CREATE INDEX IF NOT EXISTS idx_cc_sample    ON cell_counts(sample_id);
        CREATE INDEX IF NOT EXISTS idx_cc_pop       ON cell_counts(population);
        CREATE INDEX IF NOT EXISTS idx_samples_subj ON samples(subject_id);
        CREATE INDEX IF NOT EXISTS idx_subj_proj    ON subjects(project_id);
        CREATE INDEX IF NOT EXISTS idx_subj_cond    ON subjects(condition);
        CREATE INDEX IF NOT EXISTS idx_subj_treat   ON subjects(treatment);
    """)
    conn.commit()


def load_csv(conn, csv_path):
    df = pd.read_csv(csv_path)
    cursor = conn.cursor()

    for proj in df["project"].unique():
        cursor.execute("INSERT OR IGNORE INTO projects VALUES (?)", (proj,))

    subj_df = df[["subject", "project", "condition", "age", "sex", "treatment", "response"]].drop_duplicates(
        subset=["subject"]
    )
    for _, r in subj_df.iterrows():
        cursor.execute(
            "INSERT OR IGNORE INTO subjects VALUES (?,?,?,?,?,?,?)",
            (r["subject"], r["project"], r["condition"], int(r["age"]), r["sex"], r["treatment"], r["response"]),
        )

    samp_df = df[["sample", "subject", "sample_type", "time_from_treatment_start"]]
    for _, r in samp_df.iterrows():
        cursor.execute(
            "INSERT OR IGNORE INTO samples VALUES (?,?,?,?)",
            (r["sample"], r["subject"], r["sample_type"], int(r["time_from_treatment_start"])),
        )

    populations = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
    rows = [
        (r["sample"], pop, int(r[pop]))
        for _, r in df.iterrows()
        for pop in populations
    ]
    cursor.executemany(
        "INSERT INTO cell_counts (sample_id, population, count) VALUES (?,?,?)", rows
    )

    conn.commit()
    print(f"Loaded {len(df)} samples → {len(rows)} cell count records")


if __name__ == "__main__":
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
    conn = sqlite3.connect(DB_PATH)
    create_schema(conn)
    load_csv(conn, CSV_PATH)
    conn.close()
    print(f"Database ready: {DB_PATH}")
