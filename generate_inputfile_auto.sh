#!/bin/bash
# Automatically generate inputfile.txt by running readfile via the presto
# Singularity image on each FITS/filterbank file to extract cluster, RA, Dec,
# and UTC from the file itself.
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
    - cluster   : Source Name field
    - ra        : RA J2000 field
    - dec       : Dec J2000 field
    - utc_start : parsed from the filename  (YYYYMMDD-HH:MM:SS -> DD-MM-YYYYTHH:MM:SS)
  Beam info is extracted from the filename using the same patterns as
  generate_inputfile.sh.

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
                 Default: attempts to bind the data directory automatically.

  -h, --help     Show this help message and exit.

OUTPUT FORMAT:
  pointing,cluster,beam_name,beam_id,utc_start,ra,dec,fits_files,cdm

NOTES:
  - UTC is extracted from the filename. Expected pattern: YYYYMMDD-HH:MM:SS
    Output format: DD-MM-YYYYTHH:MM:SS
  - All entries will have pointing=0
  - Supported file extensions: .fil, .fits, .sf, .rf
  - Beam patterns recognised: cfbf[0-9]+, baseband[0-9]+, beam[0-9]+, Band[0-9]+

EXAMPLES:
  bash generate_inputfile_auto.sh \
    --sif /hercules/scratch/fkareem/singularity_img/presto4.sif \
    --cdm "101.0" \
    /bEDD/EDD_pipeline_data/production/pipeline_data/search1/107-23/2026-05-25/

  bash generate_inputfile_auto.sh \
    --sif docker://vishnubk/presto5:latest \
    /bEDD/EDD_pipeline_data/production/pipeline_data/search1/107-23/2026-05-25/

USAGE
}

SIF=""
CDM_LIST=""
OUTPUT="inputfile.txt"
EXTRA_BIND=""
DATA_DIRS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sif)    SIF="${2:-}";        shift 2 ;;
    --cdm)    CDM_LIST="${2:-}";   shift 2 ;;
    --output) OUTPUT="${2:-}";     shift 2 ;;
    --bind)   EXTRA_BIND="${2:-}"; shift 2 ;;
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

# ---------------------------------------------------------------------------
# extract_beam: identical logic to generate_inputfile.sh
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# extract_cdm: identical logic to generate_inputfile.sh
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
# extract_utc_from_filename: parse YYYYMMDD-HH:MM:SS from filename and
# reformat to DD-MM-YYYYTHH:MM:SS
# ---------------------------------------------------------------------------
extract_utc_from_filename() {
  local filename="$1"
  # Match YYYYMMDD-HH:MM:SS  (colons may be escaped as \: in the actual path
  # but basename will have the real colon)
  if [[ "${filename}" =~ ([0-9]{8})-([0-9]{2}:[0-9]{2}:[0-9]{2}) ]]; then
    local date_compact="${BASH_REMATCH[1]}"   # 20260525
    local time_part="${BASH_REMATCH[2]}"      # 00:54:49
    local yyyy="${date_compact:0:4}"
    local mm="${date_compact:4:2}"
    local dd="${date_compact:6:2}"
    echo "${dd}-${mm}-${yyyy}T${time_part}"
  else
    echo ""
  fi
}

# ---------------------------------------------------------------------------
# run_readfile: run readfile inside singularity on a single file and echo
# the output. Bind the file's parent directory automatically.
# ---------------------------------------------------------------------------
run_readfile() {
  local filepath="$1"
  local filedir
  filedir="$(dirname "${filepath}")"

  local bind_arg="-B ${filedir}"
  if [[ -n "${EXTRA_BIND}" ]]; then
    bind_arg="${bind_arg} -B ${EXTRA_BIND}"
  fi

  "${SNG_BIN}" exec ${bind_arg} "${SIF}" readfile "${filepath}" 2>/dev/null
}

# ---------------------------------------------------------------------------
# parse_readfile_output: extract Source Name, RA J2000, Dec J2000
# ---------------------------------------------------------------------------
parse_cluster() { grep -m1 'Source Name' <<< "$1" | awk -F'=' '{print $2}' | xargs; }
parse_ra()      { grep -m1 'RA J2000 '  <<< "$1" | awk -F'=' '{print $2}' | xargs; }
parse_dec()     { grep -m1 'Dec J2000 ' <<< "$1" | awk -F'=' '{print $2}' | xargs; }

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
echo "SIF image : ${SIF}"
echo "Output    : ${OUTPUT}"
echo "CDM list  : ${CDM_LIST_CLEAN:-<from filename>}"
echo "Data dirs (${#DATA_DIRS[@]} total):"
for d in "${DATA_DIRS[@]}"; do echo "  - $d"; done
echo ""

echo "pointing,cluster,beam_name,beam_id,utc_start,ra,dec,fits_files,cdm" > "${OUTPUT}"

pointing=0
mapfile -t files < <(find "${DATA_DIRS[@]}" -type f \( -name "*.fil" -o -name "*.fits" -o -name "*.sf" -o -name "*.rf" \) | sort)

if [[ ${#files[@]} -eq 0 ]]; then
  echo "WARNING: no matching files found."
  exit 0
fi

for filepath in "${files[@]}"; do
  filename="$(basename "${filepath}")"
  echo "Processing: ${filename}"

  # --- readfile ---
  rf_output="$(run_readfile "${filepath}")"

  cluster="$(parse_cluster "${rf_output}")"
  ra="$(parse_ra "${rf_output}")"
  dec="$(parse_dec "${rf_output}")"

  if [[ -z "${cluster}" || -z "${ra}" || -z "${dec}" ]]; then
    echo "  WARNING: readfile failed or missing fields for ${filename} — skipping."
    continue
  fi

  # --- UTC from filename ---
  utc="$(extract_utc_from_filename "${filename}")"
  if [[ -z "${utc}" ]]; then
    echo "  WARNING: could not parse UTC from filename: ${filename} — skipping."
    continue
  fi

  # --- beam ---
  read -r beam_id beam_name < <(extract_beam "${filename}")

  # --- CDM ---
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
echo "Generated ${OUTPUT} with ${lines} entries"
