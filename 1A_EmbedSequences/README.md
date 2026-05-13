
## 1A - Running ESM

### 1. Install Miniconda (if not already installed)
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### 2. Create and activate the conda environment
```bash
conda create -n esm2 -c pytorch -c conda-forge python=3.9 pytorch
conda activate esm2
pip install fair-esm
```

---

## Running ESM2 Embedding Extraction

```bash
python extract.py esm2_t30_150M_UR50D \
    <input.fasta> \
    <output_dir> \
    --repr_layers 10 15 20 25 30 \
    --include mean \
    --toks_per_batch 4096
```

`extract.py` is provided by the [fair-esm](https://github.com/facebookresearch/esm) 
package and found in `esm/scripts/`.

---

## Parameters

| Flag | Value | Description |
|------|-------|-------------|
| model | esm2_t30_150M_UR50D | ESM2 150M parameter model (30 layers) |
| `--repr_layers` | 10 15 20 25 30 | Transformer layers from which to extract representations |
| `--include` | mean | Pooling strategy (mean across residues) |
| `--toks_per_batch` | 4096 | Tokens per batch (reduce if running out of memory) |

---

## Output

For each input FASTA, produces one `.pt` (PyTorch tensor) file per protein 
containing the mean-pooled embedding at each specified layer. Output files 
are organized by species name matching the input FASTA filename.

---