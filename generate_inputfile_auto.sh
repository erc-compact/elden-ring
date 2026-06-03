#!/bin/bash
# Automatically generate inputfile.txt by running readfile via the presto
# Singularity image to extract cluster, RA, Dec, and UTC from each FITS file.
#
# File structure assumed:
#   search<beam_id>/<project>/<date>/<scan_no>/<receiver>/<SOURCE>_<FREQ>_<YYYYMMDD>-<HH:MM:SS>.fits
#
# Files sharing the same UTC timestamp across beams belong to the same pointing.
# readfile is run once per unique UTC (i.e. one representative file per pointing)
# to avoid redundant container launches. Results are cached for safe re-runs.
#
# Usage:
#   bash generate_inputfile_auto.sh --sif /path/to/presto.sif [--cdm "47.0 101.0"] [--output FILE] DATA_DIR [DATA_DIR...]

set -uo pipefail

usage() {
  cat << 'USAGE'
Automatically generate inputfile.txt using readfile (via Presto Singularity image).

USAGE:
  generate_inputfile_auto.sh --sif /path/to/presto.sif [--cdm "47.0 101.0"] [--output FILE] DATA_DIR [DATA_DIR...]

DESCRIPTION:
  Scans DATA_DIR(s) for .fits/.fil/.sf/.rf files, runs readfile inside the
  provided Presto Singularity image on each file, and extracts:
    - cluster   : Source Name field from readfile output
    - ra        : RA J2000 field from readfile output
    - dec       : Dec J2000 field from readfile output
    - utc_start : parsed from the filename (YYYYMMDD-HH:MM:SS -> DD-MM-YYYYTHH:MM:SS)
    - beam_id   : extracted from search<N> in the directory path
    - pointing  : auto-incremented per unique UTC timestamp (same UTC = same pointing)

  readfile is run once per file. Results are cached so re-running is fast
  and safe if the script is interrupted.

OPTIONS:
  --sif FILE     Path to the Presto Singularity .sif image (required).
                 Example: /hercules/scratch/fkareem/singularity_img/presto4.sif
                 Can also be a docker URI: docker://vishnubk/presto5:latest

  --cdm LIST     Space or comma-separated CDM values (optional).
                 If not provided, the script attempts to extract CDM from the
                 filename (pattern: _cdm_[0-9]+[.0-9]*).
                 If neither is available, 0.0 is written with a warning.

  --output FILE  Output file name (default: inputfile.txt).

  --bind PATHS   Extra bind paths passed to singularity/apptainer -B flag.
                 Example: --bind "/data,/scratch"

  --avoid LIST   Space or comma-separated list of source/cluster names to skip.
                 Matched against the Source Name field from readfile output.
                 Example: --avoid "CAL,NGC1234" or --avoid "CAL NGC1234"

  -h, --help     Show this help message and exit.

OUTPUT FORMAT:
  pointing,cluster,beam_name,beam_id,utc_start,ra,dec,fits_files,cdm

NOTES:
  - File structure: search<N>/<project>/<date>/<scan>/<rcvr>/<SOURCE>_<FREQ>_<YYYYMMDD>-<HH:MM:SS>.fits
  - beam_id/beam_name come from search<N> in the path; filename patterns
    (cfbf, baseband, beam, Band) are checked first as a fallback.
  - Pointing is assigned per unique UTC — all beams with the same UTC get
    the same pointing index.
  - Supported file extensions: .fil, .fits, .sf, .rf

EXAMPLES:
  bash generate_inputfile_auto.sh \
    --sif /fpra/timing/01/fazal/singularity_img/presto4.sif \
    --cdm "360.0" \
    /bEDD/EDD_pipeline_data/production/pipeline_data/search{3..7}/107-23/*/*/P170mm/

USAGE
}

SIF=""
CDM_LIST=""
OUTPUT="inputfile.txt"
EXTRA_BIND=""
AVOID_LIST=""
DATA_DIRS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sif)    SIF="${2:-}";        shift 2 ;;
    --cdm)    CDM_LIST="${2:-}";   shift 2 ;;
    --output) OUTPUT="${2:-}";     shift 2 ;;
    --bind)   EXTRA_BIND="${2:-}"; shift 2 ;;
    --avoid)  AVOID_LIST="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) DATA_DIRS+=("$1"); shift ;;
  esac
done

if [[ -z "${SIF}" || ${#DATA_DIRS[@]} -eq 0 ]]; then
  echo "ERROR: --sif and at least one DATA_DIR are required."
  usage
  exit 1
fi

# Resolve singularity/apptainer binary
if command -v apptainer &>/dev/null; then
  SNG_BIN="apptainer"
elif command -v singularity &>/dev/null; then
  SNG_BIN="singularity"
else
  echo "ERROR: neither apptainer nor singularity found in PATH."
  exit 1
fi

CDM_LIST_CLEAN="$(echo "${CDM_LIST}" | tr ',' ' ' | xargs)"
if [[ -n "${CDM_LIST_CLEAN}" ]]; then
  read -r -a CDMS <<< "${CDM_LIST_CLEAN}"
else
  CDMS=()
fi

AVOID_CLEAN="$(echo "${AVOID_LIST}" | tr ',' ' ' | xargs)"
if [[ -n "${AVOID_CLEAN}" ]]; then
  read -r -a AVOID <<< "${AVOID_CLEAN}"
else
  AVOID=()
fi

# ---------------------------------------------------------------------------
# extract_beam: directory path (search<N>) takes priority; filename patterns
# (cfbf, baseband, beam, Band) are a fallback for other data layouts.
# Returns: "<beam_id> <beam_name>"
# ---------------------------------------------------------------------------
extract_beam() {
  local filename="$1"
  local filepath="$2"
  if [[ "${filepath}" =~ /search([0-9]+)/ ]]; then
    echo "${BASH_REMATCH[1]} search${BASH_REMATCH[1]}"
  elif [[ "${filename}" =~ baseband([0-9]+) ]]; then
    local id="${BASH_REMATCH[1]}"
    printf "%s %s\n" "${id}" "$(printf "cfbf%05d" "${id}")"
  elif [[ "${filename}" =~ cfbf([0-9]+) ]]; then
    echo "${BASH_REMATCH[1]} cfbf${BASH_REMATCH[1]}"
  elif [[ "${filename}" =~ beam([0-9]+) ]]; then
    echo "${BASH_REMATCH[1]} beam${BASH_REMATCH[1]}"
  elif [[ "${filename}" =~ Band([0-9]+) ]]; then
    echo "${BASH_REMATCH[1]} Band${BASH_REMATCH[1]}"
  else
    echo "0 beam0"
  fi
}

# ---------------------------------------------------------------------------
# extract_cdm: look for _cdm_<value> in filename
# ---------------------------------------------------------------------------
extract_cdm() {
  local filename="$1"
  if [[ "${filename}" =~ _cdm_([0-9]+\.?[0-9]*) ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    echo ""
  fi
}

# ---------------------------------------------------------------------------
# extract_utc_from_filename: YYYYMMDD-HH:MM:SS -> DD-MM-YYYYTHH:MM:SS
# ---------------------------------------------------------------------------
extract_utc_from_filename() {
  local filename="$1"
  if [[ "${filename}" =~ ([0-9]{8})-([0-9]{2}:[0-9]{2}:[0-9]{2}) ]]; then
    local d="${BASH_REMATCH[1]}"
    local t="${BASH_REMATCH[2]}"
    echo "${d:6:2}-${d:4:2}-${d:0:4}T${t}"
  else
    echo ""
  fi
}

# ---------------------------------------------------------------------------
# run_readfile: one singularity exec per file, binding the storage root
# ---------------------------------------------------------------------------
run_readfile() {
  local filepath="$1"
  local bind_arg="-B $(dirname "${filepath}")"
  [[ -n "${EXTRA_BIND}" ]] && bind_arg="${bind_arg} -B ${EXTRA_BIND}"
  "${SNG_BIN}" exec ${bind_arg} "${SIF}" readfile "${filepath}" 2>/dev/null
}

parse_cluster() { grep -m1 'Source Name' <<< "$1" | awk -F'=' '{print $2}' | xargs; }
parse_ra()      { grep -m1 'RA J2000 '  <<< "$1" | awk -F'=' '{print $2}' | xargs; }
parse_dec()     { grep -m1 'Dec J2000 ' <<< "$1" | awk -F'=' '{print $2}' | xargs; }

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
CACHE_FILE="${OUTPUT%.txt}.readfile_cache.tsv"

echo "SIF image  : ${SIF}"
echo "Output     : ${OUTPUT}"
echo "Cache      : ${CACHE_FILE}"
echo "CDM list   : ${CDM_LIST_CLEAN:-<from filename>}"
echo "Avoid      : ${AVOID_CLEAN:-<none>}"
echo "Data dirs  : ${#DATA_DIRS[@]} total"
echo ""

mapfile -t ALL_FILES < <(find "${DATA_DIRS[@]}" -type f \( -name "*.fil" -o -name "*.fits" -o -name "*.sf" -o -name "*.rf" \) | sort)

if [[ ${#ALL_FILES[@]} -eq 0 ]]; then
  echo "WARNING: no matching files found."
  exit 0
fi

echo "Found ${#ALL_FILES[@]} files."

# Load existing cache (filepath -> cluster<TAB>ra<TAB>dec)
declare -A RF_CACHE
if [[ -f "${CACHE_FILE}" ]]; then
  echo "Loading cache from ${CACHE_FILE}"
  while IFS=$'\t' read -r fp cluster ra dec; do
    RF_CACHE["${fp}"]="${cluster}	${ra}	${dec}"
  done < "${CACHE_FILE}"
fi

# Run readfile on any file not already cached
echo "Running readfile on uncached files..."
for filepath in "${ALL_FILES[@]}"; do
  if [[ -n "${RF_CACHE["${filepath}"]+x}" ]]; then
    continue
  fi
  filename="$(basename "${filepath}")"
  echo "  readfile: ${filename}"
  rf_out="$(run_readfile "${filepath}")"
  cluster="$(parse_cluster "${rf_out}")"
  ra="$(parse_ra "${rf_out}")"
  dec="$(parse_dec "${rf_out}")"
  if [[ -z "${cluster}" || -z "${ra}" || -z "${dec}" ]]; then
    echo "  WARNING: readfile failed or missing fields for ${filename} — will skip."
    continue
  fi
  RF_CACHE["${filepath}"]="${cluster}	${ra}	${dec}"
  printf '%s\t%s\t%s\t%s\n' "${filepath}" "${cluster}" "${ra}" "${dec}" >> "${CACHE_FILE}"
done

# Build pointing index: unique UTC -> incrementing integer
# Same UTC across different beams = same pointing
declare -A UTC_TO_POINTING
pointing_counter=0
for filepath in "${ALL_FILES[@]}"; do
  filename="$(basename "${filepath}")"
  utc_raw="$(extract_utc_from_filename "${filename}")"
  [[ -z "${utc_raw}" ]] && continue
  if [[ -z "${UTC_TO_POINTING["${utc_raw}"]+x}" ]]; then
    UTC_TO_POINTING["${utc_raw}"]="${pointing_counter}"
    (( pointing_counter++ )) || true
  fi
done

echo ""
echo "Found ${pointing_counter} unique pointings (UTC timestamps)."
echo "Writing ${OUTPUT}..."

echo "pointing,cluster,beam_name,beam_id,utc_start,ra,dec,fits_files,cdm" > "${OUTPUT}"

skipped=0
for filepath in "${ALL_FILES[@]}"; do
  filename="$(basename "${filepath}")"

  if [[ -z "${RF_CACHE["${filepath}"]+x}" ]]; then
    (( skipped++ )) || true
    continue
  fi
  IFS=$'\t' read -r cluster ra dec <<< "${RF_CACHE["${filepath}"]}"

  # Skip sources in the avoid list
  for avoid_name in "${AVOID[@]}"; do
    if [[ "${cluster}" == "${avoid_name}" ]]; then
      echo "  [avoid] ${filename} (source: ${cluster})"
      (( skipped++ )) || true
      continue 2
    fi
  done

  utc="$(extract_utc_from_filename "${filename}")"
  if [[ -z "${utc}" ]]; then
    echo "  WARNING: could not parse UTC from ${filename} — skipping."
    (( skipped++ )) || true
    continue
  fi

  pointing="${UTC_TO_POINTING["${utc}"]}"
  read -r beam_id beam_name < <(extract_beam "${filename}" "${filepath}")

  if [[ ${#CDMS[@]} -gt 0 ]]; then
    for cdm in "${CDMS[@]}"; do
      echo "${pointing},${cluster},${beam_name},${beam_id},${utc},${ra},${dec},${filepath},${cdm}" >> "${OUTPUT}"
    done
  else
    cdm_from_file="$(extract_cdm "${filename}")"
    if [[ -n "${cdm_from_file}" ]]; then
      echo "${pointing},${cluster},${beam_name},${beam_id},${utc},${ra},${dec},${filepath},${cdm_from_file}" >> "${OUTPUT}"
    else
      echo "  WARNING: no CDM found for ${filename}, writing 0.0"
      echo "${pointing},${cluster},${beam_name},${beam_id},${utc},${ra},${dec},${filepath},0.0" >> "${OUTPUT}"
    fi
  fi
done

lines=$(( $(wc -l < "${OUTPUT}") - 1 ))
echo ""
echo "Generated ${OUTPUT} with ${lines} entries (${skipped} skipped)"
