# makeIndexFqs.py

A Python script for extracting index sequences from FASTQ read headers and writing them to separate index FASTQ files. Supports two different FASTQ format modes: **3lvl** (ScaleRNA v1) and **QS** (QuantumScale).

## Overview

This script processes FASTQ files and extracts index sequences from the read headers to create separate index FASTQ files (I1 and I2). It's designed to work with ScaleBio sequencing data and supports different barcode formats used in various kit versions.

## Usage

```bash
python makeIndexFqs.py <read1_fastq> [options]
```

### Basic Command

```bash
python makeIndexFqs.py sample_R1_001.fastq.gz --outDir ./fastqDir --mode 3lvl
```

### Command Line Options

| Option | Description | Default | Required |
|--------|-------------|---------|----------|
| `--outDir` | Output directory | `.` | No |
| `--mode` | FASTQ format mode: `3lvl` or `QS` | `3lvl` | No |
| `--no-index2` | Don't extract Index2 read | False | No |

## Modes

### 3lvl Mode (ScaleRNA v1)

**Format**: Standard Illumina FASTQ headers with index sequences in the attributes field.

**Header Example**:
```
@VH02171:15:2227VLHNX:1:1101:19144:1000 1:N:0:GCTCTCGCCT+TCGGATTCGG
```

**Index Extraction**:
- **Index1 (I1 or i7)**: Extracted from the 4th colon-separated field before the `+`
- **Index2 (I2 or i5)**: Extracted from the 4th colon-separated field after the `+`

**Usage**:
```bash
python makeIndexFqs.py sample_R1_001.fastq.gz --mode 3lvl --outDir ./fastqDir
```

### QS Mode (QuantumScale)

**Format**: QuantumScale format with partial index1 in read name and full indices in attributes.

**Header Example**:
```
@LH00659:241:22T7WLLT4:1:1101:42065:1140:CTGTCCTAATGGGGTTACCGAAGA 1:N:0:TNCAGACA+GTTCGATA
```

**Index Extraction**:
- **Index1 (I1 or i7)**: Concatenation of attribute index1 + read name index (cell barcode)
- **Index2 (I2 or i5)**: Extracted from attributes after the `+` (PCR barcode)

**Usage**:
```bash
python makeIndexFqs.py sample_R1_001.fastq.gz --mode QS --outDir ./fastqDir
```

## Output Files

The script generates compressed FASTQ files:

- **I1 file**: Contains Index1 sequences (e.g., `sample_I1_001.fastq.gz`)
- **I2 file**: Contains Index2 sequences (e.g., `sample_I2_001.fastq.gz`) - only if `--no-index2` is not used

## Examples

### ScaleRNA v1 Processing
```bash
# Process ScaleRNA v1 data
python makeIndexFqs.py experiment_R1_001.fastq.gz --mode 3lvl --outDir ./fastqDir
```

### QuantumScale Processing
```bash
# Process QuantumScale data
python makeIndexFqs.py quantum_R1_001.fastq.gz --mode QS --outDir ./fastqDir
```

## Notes

- The script automatically handles gzip compression/decompression
- Quality scores are set to a default value (ASCII 37, offset 33)
- Index sequences are validated for length consistency across all reads
- It's important to keep the `I1/I2` files together in the same directory as `R1/R2` fastqs
