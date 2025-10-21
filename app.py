# ---- REMOTE EDGES + LOCAL (TMP) IDMAP ----
import os, tempfile, pathlib, requests
import duckdb
import pandas as pd
import streamlit as st

# URLs from environment (Render -> Environment)
EDGES_URL = os.environ["EDGES_PARQUET_URL"]
IDMAP_URL = os.environ["IDMAP_PARQUET_URL"]

# Normalized location buckets
LOC_CAT_EXPR = """
CASE
  WHEN m.location IS NULL OR TRIM(m.location) = '' THEN 'Unknown'
  WHEN lower(m.location) LIKE '%mitochond%' THEN 'Mitochondrion'
  WHEN lower(m.location) LIKE '%endoplasmic reticulum%' THEN 'Endoplasmic Reticulum'
  WHEN lower(m.location) LIKE '%golgi%' THEN 'Golgi apparatus'
  WHEN lower(m.location) LIKE '%plasma membrane%' OR lower(m.location) LIKE '%cell membrane%' OR lower(m.location) LIKE '%cell surface%' THEN 'Plasma membrane'
  WHEN lower(m.location) LIKE '%lysosom%' THEN 'Lysosome'
  WHEN lower(m.location) LIKE '%endosom%' THEN 'Endosome'
  WHEN lower(m.location) LIKE '%peroxisom%' THEN 'Peroxisome'
  WHEN lower(m.location) LIKE '%cytosol%' OR lower(m.location) LIKE '%cytoplasm%' THEN 'Cytoplasm'
  WHEN lower(m.location) LIKE '%nucleus%' OR lower(m.location) LIKE '%nucleoplasm%' OR lower(m.location) LIKE '%nucleolus%' OR lower(m.location) LIKE '%chromosom%' OR lower(m.location) LIKE '%chromatin%' THEN 'Nucleus'
  WHEN lower(m.location) LIKE '%cytoskeleton%' OR lower(m.location) LIKE '%microtubule%' OR lower(m.location) LIKE '%actin%' OR lower(m.location) LIKE '%intermediate filament%' THEN 'Cytoskeleton'
  WHEN lower(m.location) LIKE '%extracellular%' OR lower(m.location) LIKE '%secreted%' THEN 'Extracellular'
  WHEN lower(m.location) LIKE '%ribosom%' THEN 'Ribosome'
  WHEN lower(m.location) LIKE '%centrosom%' THEN 'Centrosome'
  WHEN lower(m.location) LIKE '%membrane%' THEN 'Membrane'
  ELSE 'Other'
END
"""

@st.cache_resource(show_spinner=True)
def _download_to_tmp(url: str, filename: str) -> str:
    """Download a small file (idmap.parquet) with proper headers to a tmp path; reuse if already present."""
    tmp_dir = pathlib.Path(tempfile.gettempdir()) / "ppi_explorer"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    dst = tmp_dir / filename
    if dst.exists() and dst.stat().st_size > 0:
        return str(dst)
    headers = {
        "User-Agent": "Mozilla/5.0 (DuckDB httpfs workaround)",
        "Accept": "application/octet-stream",
    }
    with requests.get(url, headers=headers, stream=True, timeout=300) as r:
        r.raise_for_status()
        with open(dst, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
    return str(dst)

@st.cache_resource
def get_con():
    con = duckdb.connect()              # in-memory DB
    con.execute("INSTALL httpfs;")
    con.execute("LOAD httpfs;")
    con.execute("SET enable_http_metadata_cache = true;")
    con.execute("SET http_user_agent = 'Mozilla/5.0';")  # important for GitHub

    # 1) Read EDGES directly over HTTP
    edges_url_sql = EDGES_URL.replace("'", "''")
    con.execute(f"CREATE OR REPLACE VIEW edges AS SELECT * FROM parquet_scan('{edges_url_sql}');")

    # 2) Download the small IDMAP to /tmp, then read locally (avoids GitHub 403 with httpfs)
    idmap_local = _download_to_tmp(IDMAP_URL, "idmap.parquet")
    idmap_local_sql = idmap_local.replace("'", "''")
    con.execute(f"CREATE OR REPLACE VIEW idmap AS SELECT * FROM parquet_scan('{idmap_local_sql}');")
    return con

# -------------- rest of your existing code below (search, fetch_interactions, UI) --------------

def get_candidates(con, query_text, limit=50):
    if query_text and " " not in query_text and len(query_text) <= 12:
        return con.execute(f"""
            WITH id_hits AS (
                SELECT id,
                       COALESCE(NULLIF(TRIM(protein_name), ''), id) AS pname,
                       {LOC_CAT_EXPR} AS loc_cat,
                       1 AS rank
                FROM idmap m WHERE id = ?
                UNION ALL
                SELECT id,
                       COALESCE(NULLIF(TRIM(protein_name), ''), id) AS pname,
                       {LOC_CAT_EXPR} AS loc_cat,
                       2 AS rank
                FROM idmap m WHERE id ILIKE ? || '%'
            ),
            name_hits AS (
                SELECT id,
                       COALESCE(NULLIF(TRIM(protein_name), ''), id) AS pname,
                       {LOC_CAT_EXPR} AS loc_cat,
                       3 AS rank
                FROM idmap m WHERE protein_name ILIKE '%' || ? || '%'
            )
            SELECT id, pname, loc_cat FROM (
                SELECT * FROM id_hits
                UNION ALL
                SELECT * FROM name_hits
            )
            GROUP BY id, pname, loc_cat, rank
            ORDER BY rank, pname
            LIMIT ?
        """, [query_text, query_text, query_text, limit]).df()
    return con.execute(f"""
        SELECT id,
               COALESCE(NULLIF(TRIM(protein_name), ''), id) AS pname,
               {LOC_CAT_EXPR} AS loc_cat
        FROM idmap m
        WHERE protein_name ILIKE '%' || ? || '%'
        ORDER BY pname
        LIMIT ?
    """, [query_text, limit]).df()

def get_display_name(con, uniprot_id):
    row = con.execute(
        "SELECT COALESCE(NULLIF(TRIM(protein_name), ''), ?) FROM idmap WHERE id = ?",
        [uniprot_id, uniprot_id]
    ).fetchone()
    return row[0] if row else uniprot_id

def get_partner_location_categories(con, uniprot_id, min_score=0.0, strength=None):
    df = con.execute(f"""
    WITH nbrs AS (
      SELECT CASE WHEN e.src = ? THEN e.dst ELSE e.src END AS partner,
             e.score, e.strength
      FROM edges e
      WHERE e.src = ? OR e.dst = ?
    )
    SELECT DISTINCT {LOC_CAT_EXPR} AS loc_cat
    FROM nbrs n
    LEFT JOIN idmap m ON m.id = n.partner
    WHERE n.score >= ?
      AND (? IS NULL OR n.strength = ?)
    ORDER BY loc_cat;
    """, [uniprot_id, uniprot_id, uniprot_id, min_score, strength, strength]).df()
    return df["loc_cat"].tolist()

def fetch_interactions(con, uniprot_id, min_score=0.0, strength=None, loc_cats=None, limit=5000):
    cat_filter = ""
    cat_params = []
    if loc_cats:
        placeholders = ",".join(["?"] * len(loc_cats))
        cat_filter = f"AND {LOC_CAT_EXPR} IN ({placeholders})"
        cat_params = loc_cats

    df = con.execute(f"""
        WITH nbrs AS (
            SELECT
              CASE WHEN e.src = ? THEN e.dst ELSE e.src END AS partner,
              e.score, e.strength
            FROM edges e
            WHERE e.src = ? OR e.dst = ?
        )
        SELECT
          n.partner AS partner_id,
          COALESCE(NULLIF(TRIM(m.protein_name), ''), n.partner) AS partner_protein_name,
          {LOC_CAT_EXPR} AS partner_location_category,
          n.score, n.strength
        FROM nbrs n
        LEFT JOIN idmap m ON m.id = n.partner
        WHERE n.score >= ?
          AND (? IS NULL OR n.strength = ?)
          {cat_filter}
        ORDER BY n.score DESC
        LIMIT ?
    """, [uniprot_id, uniprot_id, uniprot_id, min_score, strength, strength] + cat_params + [limit]).df()
    return df[["partner_id", "partner_protein_name", "partner_location_category", "score", "strength"]]

# ---------- UI ----------
st.set_page_config(page_title="PPI Explorer", layout="wide")
st.title("Protein–Protein Interaction Explorer")

con = get_con()

with st.sidebar:
    st.header("Search")
    q = st.text_input("Type a UniProt ID or protein name", placeholder="e.g., Q5T5U3 or kinase")
    candidates = get_candidates(con, q) if q else pd.DataFrame()
    selected = None
    if not candidates.empty:
        labels = [f"{r['id']} — {r['pname']} [{r['loc_cat']}]" for _, r in candidates.iterrows()]
        ids = candidates["id"].tolist()
        choice = st.selectbox("Pick a protein", options=labels)
        selected = ids[labels.index(choice)]
    else:
        st.caption("Start typing to search IDs/names.")

    st.divider()
    st.header("Filters")
    strength = st.selectbox("Strength", ["any","strong","moderate","weak"], index=0)
    strength = None if strength == "any" else strength
    min_score = st.slider("Min probability", 0.0, 1.0, 0.0, 0.01)
    topk = st.number_input("Max rows", min_value=100, max_value=100000, value=5000, step=100)

    loc_cat_choices = []
    if selected:
        opts = get_partner_location_categories(con, selected, min_score=min_score, strength=strength)
        if opts:
            loc_cat_choices = st.multiselect("Partner location (category)", options=opts, default=[])

if selected:
    selected_name = get_display_name(con, selected)
    st.subheader(f"Interactions for **{selected} — {selected_name}**")
    df = fetch_interactions(con, selected, min_score=min_score, strength=strength,
                            loc_cats=loc_cat_choices if loc_cat_choices else None,
                            limit=int(topk))
    if df.empty:
        st.info("No interactions found with current filters.")
    else:
        st.dataframe(df, height=600, width="stretch")
        st.download_button("Download CSV",
                           df.to_csv(index=False).encode("utf-8"),
                           file_name=f"{selected}_interactions.csv",
                           mime="text/csv")
else:
    st.info("Search for a UniProt ID or protein name on the left to begin.")
