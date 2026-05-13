# 1b - running OrthoFinder

## Installation

### 1. Install Miniconda (if not already installed)
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### 2. Create and activate the conda environment
```bash
conda create -n orthofinder_env -c bioconda -c conda-forge \
    orthofinder \
    mafft \
    iqtree \
    diamond
conda activate orthofinder_env
```

---

## Running OrthoFinder

Place one `.fasta` file per species in an input directory, then run:

```bash
orthofinder \
    -f <input_dir> \
    -t 16 \
    -M msa \
    -A mafft \
    -T iqtree \
    -I 1.5 \
    -S diamond_ultra_sens \
    -o <output_dir>
```

---

## Parameters

| Flag | Value | Description |
|------|-------|-------------|
| `-t` | 16 | Threads |
| `-M` | msa | MSA-based orthogroup inference |
| `-A` | mafft | Alignment tool |
| `-T` | iqtree | Tree inference tool |
| `-I` | 1.5 | MCL inflation parameter |
| `-S` | diamond_ultra_sens | Sequence search sensitivity |

---