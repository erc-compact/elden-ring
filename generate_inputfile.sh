#!/bin/bash
# Generate inputfile.txt or dada_files.csv from data directories.
#
# Usage:
#   bash generate_inputfile.sh --cluster CLUSTER --ra RA --dec DEC --utc UTC --cdm "47.0 101.0" /path/to/data
#   bash generate_inputfile.sh --dada --cluster CLUSTER --ra RA --dec DEC --utc UTC --cdm "47.0 101.0" /path/to/dada_dir1 /path/to/dada_dir2

set -euo pipefail

usage() {
  cat << 'USAGE'
Usage:
  generate_inputfile.sh --cluster CLUSTER --ra RA --dec DEC --utc UTC --cdm "47.0 101.0" DATA_DIR [DATA_DIR...]
  generate_inputfile.sh --dada --cluster CLUSTER --ra RA --dec DEC --utc UTC --cdm "47.0 101.0" DADA_DIR [DADA_DIR...]

Options:
  --dada           Generate dada_files.csv using directory/*dada globs.
  --cluster NAME   Cluster name to write in CSV.
  --ra RA          Right ascension (e.g. 18:24:32.89).
  --dec DEC        Declination (e.g. -24:52:11.4).
  --utc UTC        UTC start time (e.g. 2025-05-13T01:41:02).
  --cdm LIST       Space or comma separated CDM values (e.g. "37.0 101.0").
  --output FILE    Output file name (optional).
  -h, --help       Show this help.
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

if [[ -z "${CLUSTER}" || -z "${RA}" || -z "${DEC}" || -z "${UTC}" || -z "${CDM_LIST}" || ${#DATA_DIRS[@]} -eq 0 ]]; then
  echo "ERROR: missing required arguments."
  usage
  exit 1
fi

CDM_LIST_CLEAN="$(echo "${CDM_LIST}" | tr ',' ' ' | xargs)"
if [[ -z "${CDM_LIST_CLEAN}" ]]; then
  echo "ERROR: CDM list is empty."
  exit 1
fi
read -r -a CDMS <<< "${CDM_LIST_CLEAN}"

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
    ((pointing++))
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
    for cdm in "${CDMS[@]}"; do
      echo "${pointing},${CLUSTER},${beam_name},${beam_id},${UTC},${RA},${DEC},${filepath},${cdm}" >> "${OUTPUT}"
      ((pointing++))
    done
  done
fi

lines=$(( $(wc -l < "${OUTPUT}") - 1 ))
echo "Generated ${OUTPUT} with ${lines} entries"
