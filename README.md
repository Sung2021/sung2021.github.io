# Personal Portfolio Website — sung2021.github.io

Live: https://sung2021.github.io

## Overview

- Source: `index.Rmd` → rendered to `index.html` via `rmarkdown::render()`
- Theme: `flatly` (Bootstrap 4), floating TOC
- External dependency: Font Awesome 6.5 (CDN)
- All CSS is inlined in the `<style>` block inside `index.Rmd` — no external stylesheet

---

## CSS Classes

| Class | Usage | Border/Style |
|---|---|---|
| `.card-comp` | Base card — white box with shadow | shadow only |
| `.card-project` | Featured Projects override | pink `#e8a0b0` |
| `.pub-card` | Publications card | grey `#d0d5dd` |
| `.workflow-card` | Analysis Workflows card | blue `#a8c4e0` |
| `.about-right` | About section box | shadow only |
| `.pub-journal` | Journal / workflow name | bold, `1.0em` |
| `.pub-year` | Year accent | green `#18bc9c` |
| `.pub-title` | Paper / workflow description | `0.80em`, grey |
| `.pub-badge` | Author role / links | `0.72em`, italic |
| `.card-comp-2` | ML card icon color | blue `#2c7be5` |
| `.card-comp-3` | Translational card icon color | orange `#e67e22` |

---

## Section Structure

### 1. About `{#about}`
- Single `.about-right` box (shadow, rounded)
- Inner grey box: Domains + Expertise
- Buttons: CV (primary), GitHub, LinkedIn

### 2. Featured Projects `{#projects}`
- Pandoc fenced div `::: {.row}` + `::: {.col-md-6}` — 2열 그리드
- 4 cards: `<div class="card-comp card-project">`
- Projects: TLE3 multi-omics (NI 2024), STING Agonist (BMS), Virtual Patient Cohort (Battelle), Cell Image + RNA-seq MOFA2

### 3. Selected Publications `{#publications}`
- `<div class="row">` + `<div class="col-md-4">` — 3열 그리드 (Bootstrap HTML)
- 5 cards: `.pub-card` with `.pub-journal`, `.pub-year`, `.pub-title`, `.pub-badge`
- All 5 papers linked (PMC or DOI)

### 4. Analysis Workflows `{#workflows}`
- 5 `.workflow-card` in 3열 그리드
- Links to live pages: scRNAseq.html, RNAseq.html, ChIP_seq.html, General_workflow2.html

### 5. Experience `{#experience}`
- Bullet list: Battelle (2025–), Dana Farber (2022–2024), Hackensack (2021–2022), UMich (2017–2020)

### 6. GitHub Repositories `{#repos}`
- R chunk: `knitr::kable()` + `kableExtra::kable_styling()`
- Manual `data.frame` — 4 repos: `tle3-multiomics-NI2024`, `pubmed-local-rag`, `lung-tme-deconv-profiler`, `rnaseq_local_chat`

### 7. Contact `{#contact}`
- Email, GitHub, LinkedIn

---

## Design Principles

- **가독성 우선**: 복잡한 구조 지양, 내용이 명확하게 보이는 레이아웃
- 그리드: Pandoc fenced div (`.row`, `.col-md-*`) 또는 Bootstrap HTML div
- 헤딩: `h1 1.75em`, `h2 1.4em`, `h3 1.1em` (flatly 기본값보다 한 단계 축소)
- 아이콘: Font Awesome 6.5 (`<i class="fa-solid fa-...">` / `<i class="fa fa-...">`)
- 태그 뱃지: 인라인 코드 `` `텍스트` ``

---

## Files

| Path | 설명 |
|---|---|
| `index.Rmd` | 메인 포트폴리오 소스 |
| `index.html` | knit 결과물 (git push 배포 대상) |
| `info/resume/SungryePark_resume.pdf` | CV 다운로드 링크 대상 |
| `samplePages/scRNAseq.html` | scRNA-seq workflow (Seurat v5, Harmony) |
| `samplePages/RNAseq.html` | Bulk RNA-seq workflow (DESeq2, GSEA, ssGSEA) |
| `samplePages/ChIP_seq.html` | Genomics workflow (MACS2, peak analysis) |
| `etc/General_workflow2.html` | Omics general pipeline overview |

---

## How to update

```r
# 1. Edit index.Rmd
# 2. Knit
rmarkdown::render("index.Rmd")

# 3. Push
git add index.Rmd index.html
git commit -m "Update portfolio"
git push origin main
```
