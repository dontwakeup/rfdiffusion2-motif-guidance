# Empirical Characterization of Classifier-Free Guidance Strength in RFdiffusion2 Motif Scaffolding

## Overview

This repository implements a controlled empirical study of **classifier-free guidance strength** in RFdiffusion2 motif scaffolding.

We systematically evaluate how varying:

```
inference.classifier_free_guidance_scale
```

affects:

1. Motif core fidelity
    
2. Motif boundary continuity
    
3. Structural diversity
    
4. Failure mode frequency
    

The goal is not to introduce a new design method, but to **quantify the empirical conditioning behavior of RFdiffusion2 under a reproducible sweep protocol**.

This is a focused 3–4 day experimental sprint emphasizing determinism, falsifiability, and measurement discipline.

---

## Research Questions

1. Does increasing classifier-free guidance monotonically improve motif core RMSD?
    
2. Is there a region of diminishing returns (saturation)?
    
3. Does stronger guidance reduce structural diversity?
    
4. Does conditioning preserve motif geometry while distorting its immediate structural boundary?
    

Importantly:

> If no saturation region is observed and the trade-off is strictly monotonic, that is considered an equally valid and informative result.

---

## Motif Definition

**Status: To be committed before full sweep execution**

```
PDB ID: <TBD>
Chain: <TBD>
Motif residues: <TBD>
Length: <TBD>
```

Selection criteria:

- Mixed secondary structure (e.g., helix–loop or beta–loop)
    
- Moderate solvent exposure
    
- Not located at extreme N/C terminus
    
- Avoids trivial α-helix-only motif
    

The motif will be fixed prior to running the sweep and recorded in `run_manifest.json`.

Results are motif-context specific.

---

## Experimental Design

### Swept Parameter

We sweep the RFdiffusion2 parameter:

```
inference.classifier_free_guidance_scale
```

This controls classifier-free guidance strength during diffusion sampling.

### Sweep Levels (Pre-Committed, Log-Spaced)

```
C1 = 0.5
C2 = 1.0
C3 = 2.0
C4 = 4.0
C5 = 8.0
```

Logarithmic spacing is used due to nonlinear response behavior typical in diffusion guidance mechanisms.

Each level generates ~80–100 designs.

Total designs: ~400–500.

Sweep levels are fixed before execution and will not be adjusted post hoc.

---

## Evaluation Metrics

All metrics are computed automatically and deterministically.

---

### 1️⃣ Motif Core Fidelity

Procedure:

- Extract motif backbone atoms (N, CA, C)
    
- Kabsch alignment to reference motif
    
- Compute backbone RMSD over motif residues
    

Metric:

```
motif_core_rmsd
```

---

### 2️⃣ Motif Boundary Continuity

Definition:

- Up to ±5 residues flanking motif on each side
    
- Automatically truncated near termini
    
- Excluded if exceeding scaffold boundaries
    

Rationale:

Five residues chosen as compromise between capturing local geometric continuity and avoiding contamination from distant scaffold regions.

After motif alignment:

- Compute RMSD over boundary residues
    

Metric:

```
motif_boundary_rmsd
```

---

### 3️⃣ Structural Diversity

Per conditioning level:

- Pairwise global backbone RMSD
    
- Clustering with fixed RMSD threshold (e.g., 3 Å)
    
- Record:
    
    - `cluster_count`
        
    - `largest_cluster_fraction`
        
    - `mean_pairwise_rmsd`
        

These quantify diversity collapse under stronger guidance.

---

### 4️⃣ Steric Clash Metric

Definition:

- Count heavy-atom pairs within 2.2 Å
    
- Computed between motif atoms and scaffold atoms only
    
- Hydrogen atoms excluded
    

Rationale:

2.2 Å threshold detects severe steric clashes only. Moderate steric strain (e.g., 2.7–3.0 Å) is not captured.

Metric:

```
clash_count
```

---

### 5️⃣ Failure Mode Taxonomy

Each design is labeled if exhibiting:

- Motif distortion (core RMSD above threshold)
    
- Boundary collapse
    
- Global fold ignoring motif geometry
    
- Excessive steric clashes
    
- Degenerate/repetitive scaffold structures
    

Metric:

```
failure_label
```

Failure rates are reported per conditioning level.

---

## Reproducibility Protocol

Each sweep execution produces:

```
runs/<timestamp>/
    ├── raw_outputs/
    ├── run_manifest.json
    └── config.yaml
```

`run_manifest.json` includes:

- RFdiffusion2 commit hash
    
- inference.classifier_free_guidance_scale value
    
- Random seeds
    
- Number of designs generated
    
- Timestamp
    
- Hardware metadata
    
- Motif definition (PDB, chain, residues)
    

All analysis scripts are deterministic.

No manual filtering prior to scoring.

---

## Repository Structure

```
inputs/        Reference PDB + motif definition
runs/          Raw RFdiffusion2 outputs (not version controlled)
analysis/      scores.csv + clustering results
scripts/       Sweep + scoring + analysis scripts
docs/          Generated plots
MEMO.md        Interpretation and observations
```

Primary output:

```
analysis/scores.csv
```

Columns:

- conditioning_level
    
- design_id
    
- motif_core_rmsd
    
- motif_boundary_rmsd
    
- clash_count
    
- cluster_id
    
- mean_pairwise_rmsd
    
- failure_label
    

---

## Workflow

1. Execute conditioning sweep:
    

```
bash scripts/run_sweep.sh
```

2. Score designs:
    

```
python scripts/score_designs.py
```

3. Analyze diversity:
    

```
python scripts/analyze_diversity.py
```

4. Generate plots:
    

```
python scripts/plot_results.py
```

5. Document findings in `MEMO.md`.
    

---

## Interpretation Criteria

Evidence for saturation would include:

- Plateau in motif_core_rmsd improvement
    
- Continued reduction in diversity beyond plateau
    
- Increased boundary instability at high guidance
    

If RMSD and diversity trade off smoothly without plateau, the result will be reported as a monotonic conditioning response.

---

## Limitations

- Single motif and protein context
    
- No wet-lab validation
    
- No Rosetta refinement
    
- Diversity measured via backbone RMSD only
    
- Moderate sample size (~500 designs)
    

This is a small-scale empirical characterization study, not a general theory of diffusion conditioning.

---

## Why This Matters

Classifier-free guidance strength in RFdiffusion2 is often tuned heuristically.

This repository provides:

- Pre-committed sweep design
    
- Deterministic evaluation pipeline
    
- Boundary continuity analysis
    
- Failure mode characterization
    
- Explicit null-result framing
    

The emphasis is disciplined empirical analysis rather than performance demonstration.
