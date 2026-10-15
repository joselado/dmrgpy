---
name: docs-sync
description: Update the user guide website (docs/guide/*.md, built with MkDocs) and the development notes docs/development.md after a change to dmrgpy, keeping the guide's pages building under strict mkdocs. Use this whenever a feature, method, submode, backend behaviour or dispatch rule has just landed or changed, whenever the user asks to document something or to update the user guide, and whenever a doc page is suspected of having drifted from the code. The standard failures are a guide page missing from the mkdocs.yml nav (so it is in neither the menu nor the PDF).
---

# Keeping the documentation in step

Two documents:

- The user guide, `docs/guide/*.md`, is physics-facing. What calculation each
  method performs, the formula behind it, the arguments that matter and what to
  look for in the result. A new model, a new `Many_Body_Chain` method, a new
  dynamical-correlator submode, a new post-processing tool all land here. It is
  a website (<https://joselado.github.io/dmrgpy/>) with one page per topic, and
  these pages are its only source: the PDF is generated from them, and
  `docs/user_guide.md` is only a pointer.
- `docs/development.md` is the development document, architecture-facing. Backends, dispatch,
  directory layout, the conventions a later change has to respect. A change
  that only moves numbers does not belong here; a change to which backend
  answers which call does.

Decide which document the change belongs in before writing anything. Many
changes belong in the user guide alone.

## The user guide pages

Find the page the change belongs to (the `nav` in `mkdocs.yml` lists them by
topic) and edit it there, under the `##` heading of its topic. A topic that
needs a page of its own gets a new file and an entry in `nav`, or it is in
neither the menu nor the PDF. Rules that keep the pages building:

- Math is `$...$` inline and `$$...$$` displayed. Inside a table row a `|` in
  math is written `\vert`, `\lvert` or `\rvert`, or the table splits.
- A reference to another topic is a link by page name,
  `[Ground state](ground-state.md)` or `[...](ground-state.md#energy-fluctuation)`
  for one heading, never "section N".
- A backend restriction that deserves a box is a `!!! note "Backends"`
  admonition with its body indented four spaces; the PDF turns it into a quote.

Build both before declaring it done:

```bash
pip install mkdocs==1.6.1 mkdocs-material==9.7.7   # once
mkdocs build --site-dir /tmp/site                  # strict: fails on a broken link or anchor
python docs/build_guide_pdf.py /tmp/guide.pdf      # pandoc + xelatex
```

Pushing to `master` rebuilds and deploys the site and the PDF
(`.github/workflows/guide.yml`).

## The development notes

`docs/development.md` is a single Markdown file, with no other format to keep
in step. Add to the section the change belongs to rather than appending at the
end (`grep -n '^## ' docs/development.md` lists them).

## What a page should say

The user guide's style is fixed (see `CLAUDE.md`, "Documentation"): brief,
physics-first, for physicists who do not care how the code works inside.
The "Ground state" page, `docs/guide/ground-state.md`, is the reference; read
it before writing and match its length per topic and its tone. No history, no internals,
no timings, every snippet run before it goes in.

The house shape, and the reason the guide reads the way it does: open with one
sentence saying what we will now see, give the physics before the code, then
read the result aloud afterwards, saying what to look for in the arrays rather
than promising a figure. Tie the quantity to what measures it in a clause, be
specific about approximations, and point at the example under `examples/` that
runs it. An argument list is a catalogue, one bullet per argument, what it buys
you and when to reach for it.

Where the change moved numbers rather than behaviour, say so in the same terms
the audit records use: which calls, how large, and that results from before are
not comparable. That goes in `docs/numerical_changes.md`, never in the user
guide, which stays brief and physics-first: a short physics explanation per
capability, then the code, then the few keywords that matter, in simple English.
