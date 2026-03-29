# Personal Portfolio Website — index.Rmd

## Overview

- R Markdown (`index.Rmd`) → rendered to `index.html` via `rmarkdown::render()`
- Theme: `flatly` (Bootstrap 4), floating TOC
- External dependency: Font Awesome 6.5 (CDN)
- No separate `style.css` — all CSS is inlined in `<style>` block at the top of `index.Rmd`

---

## CSS Classes (defined in `<style>` block)

| Class | Usage |
|---|---|
| `.card-comp` | Project & Competency cards — white box with shadow |
| `.card-comp-2` | ML card icon color (blue) |
| `.card-comp-3` | Translational Medicine card icon color (orange) |
| `.about-right` | About section box — white, rounded, shadow |
| `.pub-card` | Publication card — white box with shadow |
| `.pub-journal` | Journal name — bold, dark |
| `.pub-year` | Year — green accent |
| `.pub-title` | Paper title — small, grey |
| `.pub-badge` | Author role (e.g. First author) — small italic |

---

## Section Structure

### 1. About `{#about}`

```
<div class="about-right">
  큰 타이틀 (font-size 1.5em)
  서브타이틀 (italic)
  본문 2단락
  내부 박스 (배경 #f4f6f8): Domains + Expertise
  버튼: CV / GitHub / LinkedIn
</div>
```

### 2. Featured Projects `{#projects}`

- Pandoc fenced div `::: {.row}` + `::: {.col-md-6}` 2열 그리드
- 각 카드: `<div class="card-comp">` 로 감싸기
- 내용: 제목(`###`), bullet 3~4개, 태그 (인라인 코드 `` ` ``)

### 3. Selected Publications `{#publications}`

- `<div class="row">` + `<div class="col-md-4">` 3열 그리드 (순수 HTML)
- 각 카드: `<div class="pub-card">`
- 내부 구조:
  ```html
  <div class="pub-journal">저널명 <span class="pub-year">연도</span></div>
  <div class="pub-title">논문 제목</div>
  <div class="pub-badge">기여 역할 (선택)</div>
  ```

### 4. Experience `{#experience}`

- 단순 Markdown bullet list

### 5. Contact `{#contact}`

- 단순 Markdown bullet list (Email, GitHub, LinkedIn)

---

## Design Principles

- **가독성 우선**: 복잡한 구조 지양, 내용이 명확하게 보이는 레이아웃
- 그리드: Pandoc fenced div (`.row`, `.col-md-*`) 또는 Bootstrap HTML div 직접 사용
- 박스: `.card-comp` 또는 `about-right` 클래스, 또는 인라인 `style` 속성으로 경량 박스
- 아이콘: Font Awesome (`<i class="fa-solid fa-...">`)
- 태그 뱃지: 인라인 코드 `` `텍스트` ``

---

## Files

| Path | 설명 |
|---|---|
| `index.Rmd` | 메인 포트폴리오 소스 |
| `index.html` | knit 결과물 (git push 대상) |
| `info/resume/SungryePark_resume.pdf` | CV 다운로드 링크 대상 |
