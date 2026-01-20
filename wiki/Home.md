# ELDEN-RING Pipeline

**Effelsberg Large-scale Data Exploration with Nextflow for Robust Identification of New Globular cluster pulsars**

ELDEN-RING is a GPU-accelerated Nextflow pipeline for detecting pulsar candidates in large-scale radio astronomy observational data.

## Key Features

- **RFI Mitigation**: Automated Radio Frequency Interference detection and filtering using spectral kurtosis
- **GPU-Accelerated Search**: Fast periodicity searches using peasoup with NVIDIA CUDA
- **Candidate Processing**: Folding with PulsarX and intelligent filtering
- **Machine Learning Classification**: PICS-based candidate scoring and Alpha-Beta-Gamma ranking
- **Multi-beam Support**: Parallel processing of multiple telescope beams
- **Advanced Capabilities**: Coherent dedispersion, filterbank stacking, segmented searches, and resume support

## Quick Links

| Page | Description |
|------|-------------|
| [Getting Started](Getting-Started) | Installation, requirements, and quick start guide |
| [Configuration](Configuration) | Detailed parameter reference and cluster profiles |
| [Workflows](Workflows) | Available workflows and their purposes |
| [Input Formats](Input-Formats) | CSV input file specifications |
| [Output Structure](Output-Structure) | Understanding pipeline outputs |
| [Troubleshooting](Troubleshooting) | Common issues and solutions |

## Pipeline Overview

```
Input Files (filterbank/DADA)
         ↓
    RFI Filtering
         ↓
  Beam Stacking (optional)
         ↓
    Segmentation
         ↓
  GPU Periodicity Search (peasoup)
         ↓
   Candidate Parsing
         ↓
  PulsarX Folding
         ↓
  ML Classification (PICS)
         ↓
  CandyJar Tarballs
```

## Supported Telescopes

- Effelsberg 100-m Radio Telescope
- MeerKAT

## Requirements

| Component | Version |
|-----------|---------|
| Nextflow | 21.10.0+ |
| Singularity/Apptainer | 3.0+ |
| NVIDIA GPU | CUDA-enabled (for peasoup) |

## Citation

If you use ELDEN-RING in your research, please cite the appropriate publications.

## Support

For issues and feature requests, please use the GitHub Issues tracker.
