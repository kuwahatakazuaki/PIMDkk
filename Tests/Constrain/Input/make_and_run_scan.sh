#!/usr/bin/env bash
set -euo pipefail
# set +x

# Edit these values to change the scan points.
# DISTANCE_START=0.50
# DISTANCE_END=2.60
# DISTANCE_STEP=0.05
DISTANCE_START=2.70
DISTANCE_END=4.00
DISTANCE_STEP=0.2

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
parent_dir="$(cd "$script_dir/.." && pwd)"
template_input="$script_dir/input.inp"
template_runner="$script_dir/run_pimd.sh"
scan_log="$parent_dir/scan.log"

if [[ ! -f "$template_input" ]]; then
  echo "ERROR: template input is missing: $template_input" >&2
  exit 1
fi

if [[ ! -f "$template_runner" ]]; then
  echo "ERROR: runner script is missing: $template_runner" >&2
  exit 1
fi

last_run_index="$(
  find "$parent_dir" -maxdepth 1 -type d -name 'run[0-9]*' -print |
    sed 's|.*/run||' |
    awk '/^[0-9]+$/ { if ($1 > max) max = $1 } END { print max + 0 }'
)"

run_index=$((last_run_index + 1))

printf "scan start=%s distance_start=%s distance_end=%s distance_step=%s first_run=run%d\n" \
  "$(date '+%Y-%m-%d %H:%M:%S')" \
  "$DISTANCE_START" \
  "$DISTANCE_END" \
  "$DISTANCE_STEP" \
  "$run_index" >> "$scan_log"

while IFS= read -r distance; do
  run_dir="$parent_dir/run${run_index}"
  mkdir "$run_dir"
  mkdir "$run_dir/Result"

  awk -v distance="$distance" '
    /^\$cons_val/ {
      print
      getline
      printf "%sd0\n", distance
      next
    }

    /^\$Coords/ {
      in_coords = 1
      coord_count = 0
      print
      next
    }

    in_coords && /^\$end Coords/ {
      in_coords = 0
      print
      next
    }

    in_coords && $1 !~ /^\$/ && NF > 0 {
      coord_count++
      if (coord_count == 2) {
        printf "%-4s %.2f   0.0   0.0\n", $1, distance + 0.0
        next
      }
    }

    { print }
  ' "$template_input" > "$run_dir/input.inp"

  cp "$template_runner" "$run_dir/run_pimd.sh"
  chmod u+x "$run_dir/run_pimd.sh"

  {
    printf "run%d cons_val=%sd0 start=%s\n" "$run_index" "$distance" "$(date '+%Y-%m-%d %H:%M:%S')"
    (
      cd "$run_dir"
      bash run_pimd.sh
    )
    printf "run%d cons_val=%sd0 end=%s\n" "$run_index" "$distance" "$(date '+%Y-%m-%d %H:%M:%S')"
  } >> "$scan_log" 2>&1

  run_index=$((run_index + 1))
done < <(
  awk \
    -v start="$DISTANCE_START" \
    -v end="$DISTANCE_END" \
    -v step="$DISTANCE_STEP" \
    'BEGIN {
      if (step <= 0) {
        exit 1
      }
      for (distance = start; distance <= end + step / 2; distance += step) {
        printf "%.2f\n", distance
      }
    }'
)

printf "scan end=%s last_run=run%d\n" \
  "$(date '+%Y-%m-%d %H:%M:%S')" \
  "$((run_index - 1))" >> "$scan_log"
