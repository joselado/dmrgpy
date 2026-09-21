---
name: docs-sync
description: Update docs/user_guide.{md,tex} and docs/documentation.{md,tex} together after a change to dmrgpy, keeping the Markdown and LaTeX versions of each in step and verifying the .tex still compiles under pdflatex. Use this whenever a feature, method, submode, backend behaviour or dispatch rule has just landed or changed, whenever the user asks to document something or to update the user guide, and whenever a doc section is suspected of having drifted from the code. Forgetting the .tex half is the standard failure here, and it is invisible until someone rebuilds the PDF.
---

# Keeping the documentation in step

Four files, two documents, two formats each:

- `docs/user_guide.{md,tex}` is physics-facing. What calculation each method
  performs, the formula behind it, the arguments that matter and what to look
  for in the result. A new model, a new `Many_Body_Chain` method, a new
  dynamical-correlator submode, a new post-processing tool all land here.
- `docs/documentation.{md,tex}` is architecture-facing. Backends, dispatch,
  directory layout, the conventions a later change has to respect. A change
  that only moves numbers does not belong here; a change to which backend
  answers which call does.

Decide which document the change belongs in before writing anything. Many
changes belong in the user guide alone.

## The two formats are one document

The `.md` and the `.tex` carry the same content in the same order under the same
section names, and they drift silently because the `.md` is the one people read
in the repository and the `.tex` is only noticed when the PDF is rebuilt.

Write the `.md` first, then port the same text into the `.tex`, matching the
section it belongs to rather than appending at the end. Cheap checks that the
two are still in step:

```bash
grep -n '^## ' docs/user_guide.md                    # sections, in order
grep -n '\\section' docs/user_guide.tex              # the same list, same order
grep -c '<what you just added>' docs/user_guide.md docs/user_guide.tex
```

Match the starred forms too (`\subsection*{...}`), or a section that is present
reads as missing.

Translating as you go: a fenced code block becomes `lstlisting`, inline code
becomes `\texttt{}`, and mathematics is written properly in both, `$t_2=0.2$` in
the Markdown and the same in the LaTeX rather than left as plain text.
Underscores, ampersands, percent signs and carets in prose have to be escaped in
the `.tex`, which is the usual cause of a build that was fine yesterday.

## Compile before declaring it done

```bash
cd docs && pdflatex -interaction=nonstopmode <whichever you edited>.tex
```

Run it twice when the table of contents or a reference moved. It has to come out
with no errors, and ideally with no overfull-hbox warnings, which in practice
come from a long code line or an unbroken identifier in prose. Read the log
rather than trusting the exit status.

## What a section should say

The house shape, and the reason the guide reads the way it does: open with one
sentence saying what we will now see, give the physics before the code, then
read the result aloud afterwards, saying what to look for in the arrays rather
than promising a figure. Tie the quantity to what measures it in a clause, be
specific about approximations, and point at the example under `examples/` that
runs it. An argument list is a catalogue, one bullet per argument, what it buys
you and when to reach for it.

Where the change moved numbers rather than behaviour, say so in the same terms
the audit records use: which calls, how large, and that results from before are
not comparable. `docs/user_guide.md`'s closing section is where that lives.
