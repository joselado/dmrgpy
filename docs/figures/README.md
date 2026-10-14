# Figures of the user guide

Each figure of the user guide has its own folder here, `<id>/`, with two
scripts. `compute.py` runs the calculation with dmrgpy and writes
`data.npz`. `plot.py` reads only `data.npz` and writes
`docs/guide/img/<id>.png` (light theme) and `<id>-dark.png` (dark theme).
`plot.py` does not import dmrgpy, so a figure can be replotted anywhere.
The plan of all figures, with the conventions, is `plan.md`.

## Shared files

- `style.py`: the matplotlib style. A `plot.py` defines `draw(theme)` and
  ends with `style.main(draw, "<id>")`, which writes both PNG files.
- `common.py`: `common.start(params)` prints the parameters,
  `common.save(__file__, params, **arrays)` writes `data.npz` with the
  parameters, the git commit, the host and the wall time as JSON under
  `params`, and `common.load(path)` reads it back.
- `run.sh <id>`: runs `compute.py` with one BLAS thread, `PYTHONPATH` on
  this checkout, at most three jobs per machine (each pinned to its own
  P-core), and a memory cap of 8 GB. It writes `<id>/compute.log`, and
  `<id>/DONE` or `<id>/FAILED` at the end.
- `d22.sh`: runs the jobs on the workstation D22 from here.
- `ground-state_convergence/`: figure 1, the template for a new figure.

## Recompute and replot here

```bash
bash docs/figures/run.sh <id>          # compute, writes <id>/data.npz
cd docs/figures/<id> && python3 plot.py
```

## Run on D22

```bash
bash docs/figures/d22.sh push          # copy the code (never results, never *.so)
bash docs/figures/d22.sh run <id> ...  # start the jobs, detached; they queue on 3 slots
bash docs/figures/d22.sh status        # RUNNING, DONE or FAILED per figure
bash docs/figures/d22.sh pull <id> ... # fetch data.npz and compute.log
bash docs/figures/d22.sh plot <id> ... # plot here
```

On D22 the jobs use CPUs 0, 2 and 8, one thread of three different
P-cores, since CPUs 4 to 7 are often busy with other work. `D22_CPUS`
changes this.

## Add a figure to a page

```markdown
![One sentence that reads the physics off the plot.](img/<id>.png#only-light)
![One sentence that reads the physics off the plot.](img/<id>-dark.png#only-dark)
```

The website shows one image per theme, and the PDF
(`docs/build_guide_pdf.py`) keeps the light one. Check both with
`mkdocs build` and `python docs/build_guide_pdf.py /tmp/guide.pdf`.
