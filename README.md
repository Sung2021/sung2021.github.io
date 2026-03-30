# sung2021.github.io

**Personal portfolio site of Sung Rye Park, Computational Biologist**  
Live: https://sung2021.github.io

---

## About

Computational Biologist with 7+ years of experience in multi-omics integration and translational biomarker analysis.  
Pharma collaborations: BMS, Takeda, Japan Tobacco  
Currently at Battelle — applying AI/ML to analytical pipelines and scalable omics workflows.

**Domains:** Immuno-oncology · Metabolic disease (NAFLD) · Neurodegeneration (AD, TBI)  
**Expertise:** Single-cell · Multi-omics · Machine Learning · Pipeline development

---

## Portfolio Sections

| Section | Description |
|---|---|
| **About** | Background, domains, expertise, CV / GitHub / LinkedIn |
| **Featured Projects** | TLE3 multi-omics (NI 2024), STING Agonist (BMS), Virtual Patient Cohort, Cell Image + RNA-seq |
| **Selected Publications** | 5 peer-reviewed papers (Nature Immunology, Cell Reports Medicine, Cell, Cell Reports, AJP) |
| **Analysis Workflows** | Linked pages: scRNA-seq, Bulk RNA-seq, ChIP/ATAC, NK/ML, General Pipelines |
| **Experience** | Battelle (2025–), Dana Farber (2022–2024), Hackensack Meridian (2021–2022), UMich (2017–2020) |
| **GitHub Repositories** | Selected repos: multi-omics, RAG tools, deconvolution, RNA-seq |
| **Contact** | Email, GitHub, LinkedIn |

---

## Stack

- **Source:** `index.Rmd` → rendered to `index.html` (R Markdown, `rmdformats::flatly`)
- **Workflow pages:** `samplePages/`, `etc/` — same flatly theme
- **Styling:** Inline CSS + Bootstrap 4 grid + Font Awesome 6.5

---

## Repository Structure

```
index.Rmd / index.html     # Main portfolio
samplePages/               # Analysis workflow pages (scRNAseq, RNAseq, ChIP, TCR...)
etc/                       # Additional workflow and reference pages
info/resume/               # CV PDF
sampledata/                # Example datasets
png/                       # Figures
```

> Developer notes / page structure details → see [INDEX_STRUCTURE.md](INDEX_STRUCTURE.md)
