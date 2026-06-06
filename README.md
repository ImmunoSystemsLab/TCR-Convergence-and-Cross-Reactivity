# TCR Convergence and Cross-Reactivity

## Overview

This project explores T-cell receptor (TCR) repertoire organization across multiple antigen epitopes using publicly available TCR sequence datasets.

The analysis focuses on two key questions:

1. Do epitope-specific TCR repertoires exhibit convergent receptor architectures?
2. Are identical TCR clonotypes shared across different epitopes, suggesting potential cross-reactivity?

The project combines repertoire characterization, diversity analysis, motif discovery, clonotype mapping, and machine learning-based epitope prediction.

---

## Objectives

- Characterize TCR repertoires associated with known epitopes.
- Identify conserved α-chain features within dominant repertoires.
- Detect shared clonotypes across epitopes.
- Quantify repertoire diversity and overlap.
- Evaluate epitope prediction performance using TCR sequence features.

---

## Dataset

Analyzed epitopes:

| Epitope | Number of TCRs |
|----------|----------|
| ELAGIGILTV | 2289 |
| IMDQVPFSV | 56 |
| KTWGQYWQV | 16 |
| AAGIGILTV | 13 |
| YLEPGPVTA | 11 |

Total mapped TCRs: **2385**

---

## Analysis Pipeline

### 1. TCR–Epitope Mapping

Construction of complete TCR-to-epitope associations.

Outputs:

- Epitope-specific clonotype sets
- Shared clonotype identification
- Repertoire overlap assessment

---

### 2. Conserved α-Chain Analysis

Investigation of α-chain repertoire architecture.

Metrics:

- TRAV usage
- TRAJ usage
- CDR3α motifs
- Shannon diversity

Key observation:

- ELAGIGILTV-associated repertoires displayed a highly constrained α-chain architecture despite extensive β-chain diversity.

---

### 3. Cross-Reactivity Analysis

Identification of clonotypes observed in more than one epitope-specific repertoire.

Examples:

| Shared TCR | Associated Epitopes |
|------------|--------------------|
| CAISEVGVGQPQHF | ELAGIGILTV, AAGIGILTV |
| CASSLSFGTEAFF | ELAGIGILTV, AAGIGILTV |
| CASSWSFGTEAFF | ELAGIGILTV, AAGIGILTV |
| CAWSETGLGMGGWQF | ELAGIGILTV, AAGIGILTV |
| CAWSETGLGTGELFF | ELAGIGILTV, AAGIGILTV |

These shared clonotypes are consistent with potential TCR cross-reactivity.

---

### 4. Diversity Analysis

Metrics calculated:

- Shannon Diversity Index
- Simpson Diversity Index

Purpose:

- Quantify repertoire complexity
- Compare diversity between epitopes
- Evaluate α-chain versus β-chain variability

---

### 5. Repertoire Network Analysis

Network-based comparison of epitope repertoires using shared motif distributions.

Outputs:

- Motif overlap matrix
- Epitope similarity relationships
- Potential cross-reactive clusters

---

### 6. Epitope Prediction

Machine learning model trained to predict epitope specificity from TCR repertoire features.

Evaluation:

- Precision
- Recall
- F1-score
- Overall accuracy

Observation:

- High performance for dominant repertoires.
- Prediction of rare epitopes remained challenging because of severe class imbalance and repertoire overlap.

---

## Main Findings

### Conserved Repertoire Architecture

ELAGIGILTV-specific repertoires exhibited:

- Strong α-chain convergence
- Reduced α-chain diversity
- Extensive β-chain diversity

This suggests convergent immune solutions for antigen recognition.

---

### Shared Clonotypes Across Epitopes

Multiple clonotypes were observed in more than one epitope-associated repertoire.

This finding is consistent with:

- TCR promiscuity
- Repertoire overlap
- Potential cross-reactivity

---

### Prediction Challenges

Despite complete TCR–epitope mapping:

- Repertoire overlap persisted.
- Rare epitope classes remained underrepresented.
- Accurate epitope prediction was limited by class imbalance.

---

## Repository Structure

```text
.
├── analysis.py
├── cross_reactivity.py
├── canonical_repertoire.py
├── network.py
├── TCR_specificity_prediction.py
├── alpha_chain_conserved_signature.py
├── ELAGIGILTV.py
├── AAGIGILTV.py
├── IMDQVPFSV.py
├── KTWGQYWQV.py
├── YLEPGPVTA.py
└── README.md
```

---

## Requirements

```bash
pip install pandas numpy scipy scikit-learn networkx matplotlib seaborn
```

---

## Running the Analysis

Example:

```bash
python analysis.py
```

Cross-reactivity:

```bash
python cross_reactivity.py
```

Conserved α-chain analysis:

```bash
python alpha_chain_conserved_signature.py
```

Epitope prediction:

```bash
python TCR_specificity_prediction.py
```

---

## Conclusion

TCR antigen recognition is shaped by both convergent α-chain architectures and a degree of clonotype-level cross-reactivity, indicating that immune repertoires reuse common structural solutions rather than exploring sequence space randomly.
