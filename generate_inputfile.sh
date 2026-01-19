#!/bin/bash
# Generate inputfile.txt or dada_files.csv from data directories.
#
# Usage:
#   bash generate_inputfile.sh --cluster CLUSTER --ra RA --dec DEC --utc UTC --cdm "47.0 101.0" /path/to/data
#   bash generate_inputfile.sh --dada --cluster CLUSTER --ra RA --dec DEC --utc UTC --cdm "47.0 101.0" /path/to/dada_dir1 /path/to/dada_dir2

set -uo pipefail

usage() {
  cat << 'USAGE'
Generate inputfile.txt or dada_files.csv from data directories.

USAGE:
  generate_inputfile.sh --cluster CLUSTER --ra RA --dec DEC --utc UTC [--cdm "47.0 101.0"] DATA_DIR [DATA_DIR...]
  generate_inputfile.sh --dada --cluster CLUSTER --ra RA --dec DEC --utc UTC [--cdm "47.0 101.0"] DADA_DIR [DADA_DIR...]

DESCRIPTION:
  This script scans data directories for filterbank (.fil), FITS (.fits),
  sigproc (.sf), or raw (.rf) files and generates a CSV file with metadata
  for processing. It extracts beam information from filenames and can optionally
  extract CDM (Cold Dispersion Measure) values from filenames if not provided.

OPTIONS:
  --dada           Generate dada_files.csv using directory/*dada globs instead
                   of searching for individual files.

  --cluster NAME   Cluster name to write in CSV (required).
                   Example: NGC1904, NGC6401

  --ra RA          Right ascension in HH:MM:SS.SS format (required).
                   Example: 18:24:32.89, 05:24:10.09

  --dec DEC        Declination in DD:MM:SS.S format (required).
                   Example: -24:52:11.4, +15:30:45.2

  --utc UTC        UTC start time in ISO 8601 format (required).
                   Example: 2025-05-13T01:41:02, 2025-11-07T04:21:56

  --cdm LIST       Space or comma separated CDM values (optional).
                   Example: "37.0 101.0" or "62.3"
                   If not provided, the script will attempt to extract CDM
                   values from filenames (e.g., NGC6401_..._cdm_101.0.sf).
                   If CDM is in filename, that value will be used.

  --output FILE    Output file name (optional).
                   Default: inputfile.txt (fits mode) or dada_files.csv (dada mode)

  -h, --help       Show this help message and exit.

EXAMPLES:
  # Generate inputfile.txt with CDM from command line
  bash generate_inputfile.sh \
    --cluster NGC1904 \
    --ra 05:24:10.09 \
    --dec -24:31:23.40 \
    --utc 2025-11-07T04:21:56 \
    --cdm "62.3" \
    /path/to/data/M79_01/

  # Generate inputfile.txt with CDM extracted from filenames
  bash generate_inputfile.sh \
    --cluster NGC6401 \
    --ra 17:38:36.83 \
    --dec -23:54:33.9 \
    --utc 2025-09-14T16:01:33 \
    /path/to/data/

  # Generate dada_files.csv for DADA mode
  bash generate_inputfile.sh --dada \
    --cluster NGC1904 \
    --ra 05:24:10.09 \
    --dec -24:31:23.40 \
    --utc 2025-11-07T04:21:56 \
    --cdm "47.0 101.0" \
    /path/to/dada_dir1 /path/to/dada_dir2

NOTES:
  - All entries will have pointing=0 (single pointing mode)
  - Supported file extensions: .fil, .fits, .sf, .rf
  - Beam information is extracted from filenames matching patterns:
    cfbf[0-9]+, baseband[0-9]+, beam[0-9]+, or Band[0-9]+
  - CDM extraction looks for pattern: _cdm_[0-9]+[.0-9]*

OUTPUT FORMAT:
  inputfile.txt (fits mode):
    pointing,cluster,beam_name,beam_id,utc_start,ra,dec,fits_files,cdm

  dada_files.csv (dada mode):
    pointing,dada_files,cluster,beam_name,beam_id,utc_start,ra,dec,cdm_list

USAGE
}

MODE="fits"
OUTPUT=""
CLUSTER=""
RA=""
DEC=""
UTC=""
CDM_LIST=""
DATA_DIRS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dada) MODE="dada"; shift ;;
    --cluster) CLUSTER="${2:-}"; shift 2 ;;
    --ra) RA="${2:-}"; shift 2 ;;
    --dec) DEC="${2:-}"; shift 2 ;;
    --utc) UTC="${2:-}"; shift 2 ;;
    --cdm) CDM_LIST="${2:-}"; shift 2 ;;
    --output) OUTPUT="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) DATA_DIRS+=("$1"); shift ;;
  esac
done

if [[ -z "${CLUSTER}" || -z "${RA}" || -z "${DEC}" || -z "${UTC}" || ${#DATA_DIRS[@]} -eq 0 ]]; then
  echo "ERROR: missing required arguments."
  usage
  exit 1
fi

CDM_LIST_CLEAN="$(echo "${CDM_LIST}" | tr ',' ' ' | xargs)"
if [[ -n "${CDM_LIST_CLEAN}" ]]; then
  read -r -a CDMS <<< "${CDM_LIST_CLEAN}"
else
  CDMS=()
fi

if [[ -z "${OUTPUT}" ]]; then
  if [[ "${MODE}" == "dada" ]]; then
    OUTPUT="dada_files.csv"
  else
    OUTPUT="inputfile.txt"
  fi
fi

echo "Mode: ${MODE}"
echo "Output: ${OUTPUT}"
echo "Cluster: ${CLUSTER}"
echo "RA/DEC: ${RA} ${DEC}"
echo "UTC: ${UTC}"
echo "CDM list: ${CDM_LIST_CLEAN}"
echo "Data dirs (${#DATA_DIRS[@]} total):"
for d in "${DATA_DIRS[@]}"; do
  echo "  - $d"
done
echo ""

extract_beam() {
  local text="$1"
  if [[ "${text}" =~ baseband([0-9]+) ]]; then
    local id="${BASH_REMATCH[1]}"
    printf "%s %s\n" "${id}" "$(printf "cfbf%05d" "${id}")"
  elif [[ "${text}" =~ cfbf([0-9]+) ]]; then
    echo "${BASH_REMATCH[1]} cfbf${BASH_REMATCH[1]}"
  elif [[ "${text}" =~ beam([0-9]+) ]]; then
    echo "${BASH_REMATCH[1]} beam${BASH_REMATCH[1]}"
  elif [[ "${text}" =~ Band([0-9]+) ]]; then
    echo "${BASH_REMATCH[1]} Band${BASH_REMATCH[1]}"
  else
    echo "0 beam0"
  fi
}

extract_cdm() {
  local filename="$1"
  # Match patterns like: cdm_101.0 or cdm_47.0
  if [[ "${filename}" =~ _cdm_([0-9]+\.?[0-9]*) ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    echo ""
  fi
}

pointing=0

if [[ "${MODE}" == "dada" ]]; then
  echo "pointing,dada_files,cluster,beam_name,beam_id,utc_start,ra,dec,cdm_list" > "${OUTPUT}"
  for dir in "${DATA_DIRS[@]}"; do
    if [[ ! -d "${dir}" ]]; then
      echo "WARNING: ${dir} is not a directory; skipping"
      continue
    fi
    dir="${dir%/}"
    read -r beam_id beam_name < <(extract_beam "${dir}")
    dada_glob="${dir}/*dada"
    echo "${pointing},${dada_glob},${CLUSTER},${beam_name},${beam_id},${UTC},${RA},${DEC},${CDM_LIST_CLEAN}" >> "${OUTPUT}"
  done
else
  echo "pointing,cluster,beam_name,beam_id,utc_start,ra,dec,fits_files,cdm" > "${OUTPUT}"
  mapfile -t files < <(find "${DATA_DIRS[@]}" -type f \( -name "*.fil" -o -name "*.fits" -o -name "*.sf" -o -name "*.rf" \) | sort)
  if [[ ${#files[@]} -eq 0 ]]; then
    echo "WARNING: no matching files found."
  fi
  for filepath in "${files[@]}"; do
    filename="$(basename "${filepath}")"
    read -r beam_id beam_name < <(extract_beam "${filename}")

    # If CDM list provided, use it; otherwise extract from filename
    if [[ ${#CDMS[@]} -gt 0 ]]; then
      for cdm in "${CDMS[@]}"; do
        echo "${pointing},${CLUSTER},${beam_name},${beam_id},${UTC},${RA},${DEC},${filepath},${cdm}" >> "${OUTPUT}"
      done
    else
      cdm_from_file=$(extract_cdm "${filename}")
      if [[ -n "${cdm_from_file}" ]]; then
        echo "${pointing},${CLUSTER},${beam_name},${beam_id},${UTC},${RA},${DEC},${filepath},${cdm_from_file}" >> "${OUTPUT}"
      else
        echo "WARNING: No CDM provided and couldn't extract from filename: ${filename}"
        echo "${pointing},${CLUSTER},${beam_name},${beam_id},${UTC},${RA},${DEC},${filepath},0.0" >> "${OUTPUT}"
      fi
    fi
  done
fi

lines=$(( $(wc -l < "${OUTPUT}") - 1 ))
echo "Generated ${OUTPUT} with ${lines} entries"
