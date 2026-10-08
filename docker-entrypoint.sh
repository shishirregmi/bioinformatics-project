#!/bin/sh
set -eu

if [ "$#" -gt 0 ]; then
  exec radio-transfer "$@"
fi

DATA_DIR="${DATA_DIR:-/app/data/raw}"
PORT="${PORT:-8000}"
METADATA="${METADATA:-/app/data/clinical.csv}"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"
OUTPUT_DIR="${OUTPUT_DIR:-/app/results/data-exploration-${RUN_ID}}"
BASELINE="$DATA_DIR/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz"
POST_TREATMENT="$DATA_DIR/GSE248378_Durva_Post_FPKMs.txt.gz"

mkdir -p "$DATA_DIR" "$(dirname "$OUTPUT_DIR")"

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

set -- \
  --baseline-expression "$BASELINE" \
  --post-treatment-expression "$POST_TREATMENT"
if [ -f "$METADATA" ]; then
  set -- "$@" --metadata "$METADATA"
fi
set -- "$@" --signatures /app/config/rss.json /app/config/immune.json --output "$OUTPUT_DIR"
radio-transfer explore "$@"

printf '\nDashboard: http://localhost:%s/\n' "$PORT"
printf 'Report files: %s\n' "$OUTPUT_DIR"
exec python -m http.server "$PORT" --bind 0.0.0.0 --directory "$OUTPUT_DIR"
