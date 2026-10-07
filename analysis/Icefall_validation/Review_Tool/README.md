# Icefall Model Review Tool

Lokales Review-Tool fuer die Modellvalidierung anhand der Bilder in `analysis/Icefall_validation/UIDs`.

## Start

1. Daten neu bauen:
   `python build_review_data.py`

2. Server starten:
   `python serve_review_tool.py`

3. Browser oeffnen:
   `http://127.0.0.1:8765/Review_Tool/review_tool.html`

## Ausgabe

- `review_data.json`: automatisch generierte Faelle und Modellreihen.
- `validation_results.json`: gespeicherte Entscheidungen.
- `validation_results.csv`: tabellarische Entscheidungen fuer Excel/R.

Ein Fall ist jeweils eine UID an einem Datum, wenn in diesem UID-Ordner mindestens ein Bild mit diesem Datum liegt. Videos vom selben Datum werden zusaetzlich angezeigt.

## Public version (English)

Only the review tool code and a sanitized export of the review decisions are published:
`build_review_data.py` (builds the review cases), `review_tool.html` (review interface),
`serve_review_tool.py` (local server that stores the decisions) and
`review_decisions_public.csv` (83 decisions of the final review of 27 Sep 2026: case, route,
date, decision, model state). Free-text comments, condition reports and image file names are
removed. The photographs and condition reports were posted by ice climbers in a private chat
group and on a website and are not published (third-party material, privacy); without them
the tool cannot be run.
