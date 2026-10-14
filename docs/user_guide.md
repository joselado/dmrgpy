# DMRGPY user guide

The user guide is a website: <https://joselado.github.io/dmrgpy/>. It is
also available as one PDF, <https://joselado.github.io/dmrgpy/user_guide.pdf>.

The source of the guide is the pages in `docs/guide/`, one page per topic,
in the order given by `nav` in `mkdocs.yml`. To build it locally:

```bash
pip install mkdocs==1.6.1 mkdocs-material==9.7.7
mkdocs serve                       # website at http://127.0.0.1:8000
python docs/build_guide_pdf.py     # docs/user_guide.pdf, needs pandoc and xelatex
```
