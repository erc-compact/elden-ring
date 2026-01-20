# Input Formats

This page describes the input file formats required by ELDEN-RING.

## Standard Filterbank Input (`inputfile.txt`)

The primary input format for filterbank files.

### Format

CSV file with the following columns:

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `pointing` | Integer | Yes | Pointing number (0-indexed) |
| `cluster` | String | Yes | Target cluster/source name |
| `beam_name` | String | Yes | Beam identifier (e.g., `cfbf00003`) |
| `beam_id` | Integer | Yes | Numeric beam ID |
| `utc_start` | String | Yes | Observation start time (ISO 8601) |
| `ra` | String | Yes | Right Ascension (HH:MM:SS.ss) |
| `dec` | String | Yes | Declination (DD:MM:SS.s) |
| `fits_files` | String | Yes | Path to filterbank file |
| `cdm` | Float | Yes | Coherent DM value |

### Example

```csv
pointing,cluster,beam_name,beam_id,utc_start,ra,dec,fits_files,cdm
0,NGC6401,cfbf00003,3,2025-09-14T16:01:33.9,17:38:36.53,-23:54:34.6,/data/NGC6401/beam3_cdm47.fil,47.0
0,NGC6401,cfbf00003,3,2025-09-14T16:01:33.9,17:38:36.53,-23:54:34.6,/data/NGC6401/beam3_cdm101.fil,101.0
0,NGC6401,cfbf00004,4,2025-09-14T16:01:33.9,17:38:36.53,-23:54:34.6,/data/NGC6401/beam4_cdm47.fil,47.0
1,NGC6401,cfbf00003,3,2025-09-14T17:01:33.9,17:38:36.53,-23:54:34.6,/data/NGC6401/beam3_p2.fil,47.0
```

### Notes

- Multiple rows with the same `beam_name` but different `cdm` values will be stacked when `stack_by_cdm = true`
- Different `pointing` values represent separate observations of the same target
- File paths can be absolute or relative to the working directory

### Generating with Helper Script

```bash
bash generate_inputfile.sh \
    --cluster NGC6401 \
    --ra 17:38:36.83 \
    --dec -23:54:33.9 \
    --utc 2025-09-14T16:01:33 \
    --cdm "47.0 101.0" \
    --beam-pattern "cfbf*" \
    /path/to/filterbank/data/
```

---

## DADA Baseband Input (`dada_input.csv`)

Input format for DADA baseband files (used with `run_dada_search`).

### Format

CSV file with the following columns:

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `pointing` | Integer | Yes | Pointing number |
| `dada_files` | String | Yes | Glob pattern for DADA files |
| `cluster` | String | Yes | Target cluster/source name |
| `beam_name` | String | Yes | Beam identifier |
| `beam_id` | Integer | Yes | Numeric beam ID |
| `utc_start` | String | Yes | Observation start time |
| `ra` | String | Yes | Right Ascension |
| `dec` | String | Yes | Declination |
| `cdm_list` | String | Yes | Space-separated coherent DM values |

### Example

```csv
pointing,dada_files,cluster,beam_name,beam_id,utc_start,ra,dec,cdm_list
0,/data/NGC6401/beam3/*dada,NGC6401,cfbf00003,3,2025-09-14T16:01:33.9,17:38:36.53,-23:54:34.6,47.0 101.0
0,/data/NGC6401/beam4/*dada,NGC6401,cfbf00004,4,2025-09-14T16:01:33.9,17:38:36.53,-23:54:34.6,47.0 101.0
```

### Notes

- `dada_files` supports glob patterns to match multiple DADA files
- `cdm_list` specifies all coherent DM values for digifits conversion
- Each DM value in `cdm_list` produces a separate FITS file

---

## Known Pulsars File (`known_pulsars.csv`)

Optional file listing known pulsars for candidate comparison.

### Format

| Column | Type | Description |
|--------|------|-------------|
| `name` | String | Pulsar name (e.g., J1748-2446A) |
| `ra` | String | Right Ascension |
| `dec` | String | Declination |
| `period` | Float | Period in seconds |
| `dm` | Float | Dispersion measure |

### Example

```csv
name,ra,dec,period,dm
J1748-2446A,17:48:04.85,-24:46:48.8,0.01165,242.0
J1748-2446B,17:48:05.12,-24:46:47.1,0.00354,243.1
B1937+21,19:39:38.56,+21:34:59.1,0.00156,71.04
```

---

## RFI Filter File (`rfi_filters.txt`)

Optional file defining frequency ranges to mask.

### Format

Text file with frequency ranges (MHz):

```
# Comment lines start with #
# Format: start_freq end_freq

# Mobile phone bands
925 960
1805 1880

# WiFi
2400 2500

# Satellite
1525 1559
```

### Notes

- Frequencies in MHz
- One range per line
- Comments start with `#`

---

## Birdie File (`birdies.txt`)

Known RFI frequencies (birdies) to flag.

### Format

Text file with frequencies:

```
# Known RFI frequencies (Hz)
50.0
60.0
100.0
```

---

## Ephemeris File (`.par`)

Pulsar ephemeris file for `fold_par` workflow.

### Format

Standard TEMPO/TEMPO2 parameter file:

```
PSRJ           J1748-2446A
RAJ            17:48:04.85
DECJ           -24:46:48.8
F0             85.8424
F1             -1.234e-15
DM             242.0
PEPOCH         58000.0
```

---

## File Validation

Before running the pipeline, validate your input files:

```bash
nextflow run elden.nf -entry validate_inputs \
    -c params.config \
    --files_list inputfile.txt
```

This checks:
- CSV format validity
- Required columns present
- File paths exist
- Coordinate format validity
- Numeric values in valid ranges

---

## Common Issues

### Invalid CSV Format

**Error**: `Error parsing CSV file`

**Solution**: Ensure:
- No trailing commas
- All columns present
- Consistent column count per row
- No special characters in values

### Missing Files

**Error**: `File not found: /path/to/file.fil`

**Solution**:
- Check file paths are absolute
- Verify file permissions
- Use `ls` to confirm files exist

### Coordinate Format

**Error**: `Invalid coordinate format`

**Solution**:
- RA format: `HH:MM:SS.ss` (e.g., `17:38:36.53`)
- Dec format: `DD:MM:SS.s` or `+DD:MM:SS.s` (e.g., `-23:54:34.6`)
