#!/bin/sh
set -eu

if [ "$#" -gt 0 ]; then
  exec radio-transfer "$@"
fi

DATA_DIR="${DATA_DIR:-/app/data/raw}"
PORT="${PORT:-8000}"
METADATA="${METADATA:-/app/config/clinical.csv}"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"
OUTPUT_DIR="${OUTPUT_DIR:-/app/results/radio-transfer-${RUN_ID}}"
BASELINE="$DATA_DIR/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz"
POST_TREATMENT="$DATA_DIR/GSE248378_Durva_Post_FPKMs.txt.gz"

if [ -e "$OUTPUT_DIR" ] && [ -n "$(find "$OUTPUT_DIR" -mindepth 1 -print -quit)" ]; then
  echo "Output directory must be empty to preserve previous reports: $OUTPUT_DIR" >&2
  exit 2
fi
mkdir -p "$DATA_DIR" "$OUTPUT_DIR"

if [ ! -f "$BASELINE" ] && [ ! -f "$POST_TREATMENT" ]; then
  radio-transfer download --dataset all --output "$DATA_DIR"
else
  if [ ! -f "$BASELINE" ]; then
    radio-transfer download --dataset baseline --output "$DATA_DIR"
  fi
  if [ ! -f "$POST_TREATMENT" ]; then
    radio-transfer download --dataset post-treatment --output "$DATA_DIR"
  fi
fi

radio-transfer explore \
  --baseline-expression "$BASELINE" \
  --post-treatment-expression "$POST_TREATMENT" \
  --signatures /app/config/rss.json /app/config/immune.json \
  --output "$OUTPUT_DIR/exploration"

ANALYSIS_STATUS="metadata_missing"
if [ -f "$METADATA" ]; then
  if radio-transfer analyze \
    --expression "$BASELINE" \
    --metadata "$METADATA" \
    --signatures /app/config/rss.json /app/config/immune.json \
    --output "$OUTPUT_DIR/analysis" >"$OUTPUT_DIR/analysis.log" 2>&1; then
    ANALYSIS_STATUS="complete"
    cat "$OUTPUT_DIR/analysis.log"
  else
    ANALYSIS_STATUS="failed"
    cat "$OUTPUT_DIR/analysis.log" >&2
  fi
fi

python - "$OUTPUT_DIR" "$PORT" "$ANALYSIS_STATUS" "$METADATA" <<'PY'
from html import escape
from pathlib import Path
import sys

output = Path(sys.argv[1])
port, status, metadata = sys.argv[2:]
if status == "complete":
    analysis = '<a class="card" href="analysis/analysis_report.html"><strong>Outcome analysis</strong><span>Open fixed-signature MPR comparisons and reproducibility files</span></a>'
    note = "The statistical analysis completed using the supplied metadata. Review the report's small-cohort limitations before interpreting estimates."
elif status == "failed":
    analysis = '<a class="card" href="analysis.log"><strong>Analysis log</strong><span>Open the reason the outcome analysis could not run</span></a>'
    note = "The data explorer was generated, but the outcome analysis could not complete. Check the log and resolve any metadata, identifier, or signature-gene issues."
else:
    analysis = '<article class="card disabled"><strong>Outcome analysis is waiting for verified labels</strong><span>Add a verified clinical CSV at the configured metadata path, then rerun the container.</span></article>'
    note = f"No clinical metadata file was found at <code>{escape(metadata)}</code>. The explorer does not infer treatment arms or outcomes from sample names."

html = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="description" content="Project dashboard for public NSCLC radiotherapy RNA-seq data exploration and fixed-signature analysis.">
<title>Radiotherapy signature project</title>
<style>
:root{{--ink:#173047;--muted:#607186;--paper:#f3f6f8;--line:#dce5ec;--teal:#187d78;--white:#fff}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font:16px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif}}
header{{padding:50px max(22px,calc((100vw - 900px)/2));background:linear-gradient(120deg,#10263d,#1d5364);color:#fff}}header p{{color:#d8e8ef;max-width:720px;margin:12px 0 0}}
main{{max-width:900px;margin:30px auto;padding:0 20px}}h1{{font-size:clamp(30px,5vw,46px);line-height:1.08;margin:8px 0}}h2{{font-size:20px;margin:0 0 14px}}.eyebrow{{text-transform:uppercase;letter-spacing:.12em;font-size:11px;font-weight:750;color:#8ed0c9;margin:0}}
section{{background:#fff;border:1px solid var(--line);border-radius:15px;padding:24px;margin-bottom:16px}}.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(245px,1fr));gap:13px}}.card{{display:flex;flex-direction:column;gap:5px;padding:18px;background:#f8fafb;border:1px solid var(--line);border-radius:11px;color:var(--ink);text-decoration:none}}a.card:hover{{border-color:var(--teal);box-shadow:0 3px 14px #17304712}}.card strong{{font-size:18px}}.card span,.muted{{color:var(--muted);font-size:14px}}.disabled{{background:#fff8eb}}.note{{padding:14px 16px;border-left:4px solid #d78b45;background:#fff8eb;border-radius:5px}}code{{overflow-wrap:anywhere}}footer{{color:var(--muted);padding:5px 0 30px}}
</style></head><body>
<header><p class="eyebrow">Cross-cancer radiotherapy signature transfer</p><h1>Project dashboard</h1>
<p>Explore the public NSCLC RNA-seq matrices, inspect their parameters and figures, and review the fixed-signature response analysis when verified clinical labels are available.</p></header>
<main><section><h2>Open a project view</h2><div class="cards">
<a class="card" href="exploration/index.html"><strong>Radiotherapy RNA-seq data explorer</strong><span>Baseline and post-treatment figures, data dictionary, and complete expression profiles</span></a>
{analysis}</div></section><section><h2>Analysis status</h2><p class="note">{note}</p><p class="muted">The primary comparison tests the locked breast-cancer radiosensitivity signature within the durvalumab + SBRT arm. The immune signature is secondary. Scores are association measures, not a clinical prediction model.</p></section>
<footer>Dashboard served on port {escape(port)}. Downloaded expression matrices and generated files are stored in the mounted data and results folders.</footer></main></body></html>'''
(output / "index.html").write_text(html)
PY

printf '\nDashboard: http://localhost:%s/\n' "$PORT"
printf 'Report files: %s\n' "$OUTPUT_DIR"
exec python -m http.server "$PORT" --bind 0.0.0.0 --directory "$OUTPUT_DIR"
