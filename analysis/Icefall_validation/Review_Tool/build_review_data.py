from __future__ import annotations

import csv
import json
import re
import statistics
from pathlib import Path
from urllib.parse import quote


HERE = Path(__file__).resolve().parent
VALIDATION_ROOT = HERE.parent
REPO = VALIDATION_ROOT.parent.parent
UID_ROOT = VALIDATION_ROOT / "UIDs"
META_CSV = REPO / "data" / "Koordinaten_Wasserfaelle" / "eisklettern_links_entries_diff.csv"
MODEL_DIRS = [
    REPO / "data" / "ModelRuns",
    REPO / "ALT" / "adj_model_legacy" / "data" / "plots" / "ModelRuns",
]
OUT_JSON = HERE / "review_data.json"
RESULTS_JSON = HERE / "validation_results.json"
RESULTS_CSV = HERE / "validation_results.csv"

DATE_PREFIX_RE = re.compile(r"^(\d{2})-(\d{2})-(\d{4})_")
UID_FOLDER_RE = re.compile(r"^(\d+)\s+_\s+(.+)$")
IMAGE_EXTS = {".jpg", ".jpeg", ".png", ".webp"}
VIDEO_EXTS = {".mp4", ".mov", ".m4v"}
NORTH_TYROL_BBOX = {
    "lat_min": 46.7,
    "lon_min": 10.1,
    "lat_max": 47.7,
    "lon_max": 12.2,
}


def parse_float(value: object) -> float | None:
    text = str(value or "").strip().replace(",", ".")
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def load_uid_coordinates() -> dict[int, dict[str, float | str]]:
    if not META_CSV.exists():
        return {}

    out: dict[int, dict[str, float | str]] = {}
    with META_CSV.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=";")
        for row in reader:
            try:
                uid = int(str(row.get("uid", "")).strip())
            except ValueError:
                continue
            lat = parse_float(row.get("latitude"))
            lon = parse_float(row.get("longitude"))
            if lat is None or lon is None:
                continue
            out[uid] = {
                "latitude": lat,
                "longitude": lon,
                "name": row.get("name", ""),
            }
    return out


def is_in_north_tyrol_bbox(uid: int, coordinates: dict[int, dict[str, float | str]]) -> bool:
    coord = coordinates.get(uid)
    if not coord:
        return False
    lat = coord.get("latitude")
    lon = coord.get("longitude")
    if not isinstance(lat, (float, int)) or not isinstance(lon, (float, int)):
        return False
    return (
        NORTH_TYROL_BBOX["lat_min"] <= float(lat) <= NORTH_TYROL_BBOX["lat_max"]
        and NORTH_TYROL_BBOX["lon_min"] <= float(lon) <= NORTH_TYROL_BBOX["lon_max"]
    )


def iso_from_dd_mm_yyyy(date_text: str) -> str:
    try:
        dd, mm, yyyy = date_text.split("-")
        return f"{yyyy}-{mm}-{dd}"
    except ValueError:
        return ""


def dd_mm_yyyy_from_iso(date_iso: str) -> str:
    try:
        yyyy, mm, dd = date_iso.split("-")
        return f"{dd}-{mm}-{yyyy}"
    except ValueError:
        return date_iso


def web_path(path: Path) -> str:
    rel = path.relative_to(VALIDATION_ROOT)
    return "../" + "/".join(quote(part, safe="") for part in rel.parts)


def parse_uid_folder(folder: Path) -> tuple[int | None, str]:
    match = UID_FOLDER_RE.match(folder.name)
    if not match:
        return None, folder.name
    return int(match.group(1)), match.group(2)


def parse_conditions(txt_path: Path) -> dict[str, dict[str, object]]:
    if not txt_path.exists():
        return {}

    text = txt_path.read_text(encoding="utf-8-sig")
    lines = text.splitlines()
    in_summary = False
    current: dict[str, object] | None = None
    out: dict[str, dict[str, object]] = {}

    for line in lines:
        stripped = line.strip()
        if stripped == "Zusammenfassung nach Datum":
            in_summary = True
            continue
        if stripped == "Originalnachrichten":
            if current and current.get("date"):
                out[str(current["date"])] = current
            break
        if not in_summary:
            continue
        if re.match(r"^\d{2}-\d{2}-\d{4}$", stripped):
            if current and current.get("date"):
                out[str(current["date"])] = current
            current = {"date": stripped, "condition": "", "images_listed": [], "videos_listed": []}
            continue
        if current is None:
            continue
        if stripped.startswith("Zustand:"):
            current["condition"] = stripped.split(":", 1)[1].strip()
        elif stripped.startswith("- "):
            item = stripped[2:].strip()
            ext = Path(item).suffix.lower()
            if ext in IMAGE_EXTS:
                current["images_listed"].append(item)
            elif ext in VIDEO_EXTS:
                current["videos_listed"].append(item)

    return out


def collect_cases() -> tuple[list[dict[str, object]], set[int], dict[str, int]]:
    cases: list[dict[str, object]] = []
    uids: set[int] = set()
    coordinates = load_uid_coordinates()
    filter_summary = {
        "uids_outside_north_tyrol_bbox": 0,
        "uids_missing_coordinates": 0,
        "cases_outside_north_tyrol_bbox": 0,
        "cases_missing_coordinates": 0,
    }

    for folder in sorted(p for p in UID_ROOT.iterdir() if p.is_dir()):
        uid, name = parse_uid_folder(folder)
        if uid is None:
            continue
        if uid not in coordinates:
            filter_summary["uids_missing_coordinates"] += 1
            media_case_count = 0
        elif not is_in_north_tyrol_bbox(uid, coordinates):
            filter_summary["uids_outside_north_tyrol_bbox"] += 1
            media_case_count = 0
        else:
            media_case_count = -1
        uids.add(uid)
        conditions = parse_conditions(folder / "verhaeltnisse.txt")

        media_by_date: dict[str, dict[str, list[dict[str, str]]]] = {}
        for path in sorted(p for p in folder.iterdir() if p.is_file()):
            match = DATE_PREFIX_RE.match(path.name)
            if not match:
                continue
            date = f"{match.group(1)}-{match.group(2)}-{match.group(3)}"
            ext = path.suffix.lower()
            entry = {
                "file": path.name,
                "path": web_path(path),
                "source": "website" if "_website_" in path.name else ("whatsapp" if "_whatsapp" in path.name else ""),
            }
            bucket = media_by_date.setdefault(date, {"images": [], "videos": []})
            if ext in IMAGE_EXTS:
                bucket["images"].append(entry)
            elif ext in VIDEO_EXTS:
                bucket["videos"].append(entry)

        if media_case_count == 0:
            for media in media_by_date.values():
                if media["images"]:
                    if uid not in coordinates:
                        filter_summary["cases_missing_coordinates"] += 1
                    else:
                        filter_summary["cases_outside_north_tyrol_bbox"] += 1
            uids.discard(uid)
            continue

        for date, media in sorted(media_by_date.items(), key=lambda item: iso_from_dd_mm_yyyy(item[0])):
            if not media["images"]:
                continue
            condition = conditions.get(date, {})
            sources = sorted({m["source"] for m in media["images"] + media["videos"] if m.get("source")})
            date_iso = iso_from_dd_mm_yyyy(date)
            coord = coordinates.get(uid, {})
            cases.append(
                {
                    "id": f"{uid:03d}_{date_iso}",
                    "uid": f"{uid:03d}",
                    "uid_int": uid,
                    "name": name,
                    "latitude": coord.get("latitude"),
                    "longitude": coord.get("longitude"),
                    "date": date,
                    "date_iso": date_iso,
                    "sources": sources,
                    "condition": condition.get("condition", ""),
                    "images": media["images"],
                    "videos": media["videos"],
                    "folder": folder.name,
                    "txt_path": web_path(folder / "verhaeltnisse.txt") if (folder / "verhaeltnisse.txt").exists() else "",
                }
            )

    return cases, uids, filter_summary


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def read_model_series(uid: int) -> dict[str, object] | None:
    path = next((model_dir / f"model_uid{uid}.csv" for model_dir in MODEL_DIRS if (model_dir / f"model_uid{uid}.csv").exists()), None)
    if path is None:
        return None

    daily: dict[str, dict[str, list[float] | bool]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            date = (row.get("date") or row.get("time") or "")[:10]
            if not re.match(r"^\d{4}-\d{2}-\d{2}$", date):
                continue
            bucket = daily.setdefault(
                date,
                {"thickness": [], "climbability": [], "tlz": [], "score_h": [], "forecast": False},
            )
            for col, key in [
                ("thickness_m", "thickness"),
                ("climbability", "climbability"),
                ("TLz", "tlz"),
                ("score_h", "score_h"),
            ]:
                raw = str(row.get(col, "")).strip()
                if not raw or raw.upper() in {"NA", "NAN", "NULL"}:
                    continue
                try:
                    bucket[key].append(float(raw))  # type: ignore[index, union-attr]
                except ValueError:
                    pass
            if str(row.get("is_forecast", "")).strip().upper() == "TRUE":
                bucket["forecast"] = True

    series: list[dict[str, object]] = []
    for date, bucket in sorted(daily.items()):
        thickness = bucket["thickness"]  # type: ignore[assignment]
        if not thickness:
            continue
        clim = bucket["climbability"]  # type: ignore[assignment]
        tlz = bucket["tlz"]  # type: ignore[assignment]
        score_h = bucket["score_h"]  # type: ignore[assignment]
        series.append(
            {
                "date": date,
                "date_display": dd_mm_yyyy_from_iso(date),
                "thickness_m_median": round(median(thickness) or 0, 4),
                "thickness_m_max": round(max(thickness), 4),
                "thickness_m_min": round(min(thickness), 4),
                "climbability_median": round(median(clim), 4) if clim else None,
                "TLz_median_C": round(median(tlz), 2) if tlz else None,
                "score_h_median": round(median(score_h), 4) if score_h else None,
                "forecast": bool(bucket["forecast"]),
            }
        )

    if not series:
        return None

    by_date = {row["date"]: row for row in series}
    return {
        "uid": f"{uid:03d}",
        "file": str(path),
        "first_date": series[0]["date"],
        "last_date": series[-1]["date"],
        "series": series,
        "by_date": by_date,
    }


def build_data() -> dict[str, object]:
    cases, uids, filter_summary = collect_cases()
    models: dict[str, object] = {}
    for uid in sorted(uids):
        model = read_model_series(uid)
        if model is not None:
            models[f"{uid:03d}"] = model

    for case in cases:
        model = models.get(str(case["uid"]))
        if not model:
            case["model_status"] = "kein Modelllauf"
            continue
        by_date = model.get("by_date", {})  # type: ignore[union-attr]
        if case["date_iso"] not in by_date:
            case["model_status"] = "kein Modellwert am Tag"
            continue
        case["model_status"] = "ok"
        case["model_day"] = by_date[case["date_iso"]]

    return {
        "version": 1,
        "created_by": "build_review_data.py",
        "notes": {
            "model_value": "Tageswert = Median der 10-min-thickness_m-Werte; max/min werden zusaetzlich angezeigt.",
            "scope": "Ein Review-Fall pro UID und Datum, wenn mindestens ein Bild im UID-Ordner liegt und die UID in der North-Tyrol-BBox liegt.",
            "north_tyrol_bbox": NORTH_TYROL_BBOX,
        },
        "paths": {
            "uid_root": str(UID_ROOT),
            "metadata_csv": str(META_CSV),
            "model_dirs": [str(path) for path in MODEL_DIRS],
        },
        "summary": {
            "cases": len(cases),
            "uids_with_cases": len(uids),
            "uids_with_model": len(models),
            "cases_with_model_day": sum(1 for c in cases if c.get("model_status") == "ok"),
            "cases_without_model_run": sum(1 for c in cases if c.get("model_status") == "kein Modelllauf"),
            "cases_without_model_day": sum(1 for c in cases if c.get("model_status") == "kein Modellwert am Tag"),
            **filter_summary,
        },
        "cases": cases,
        "models": models,
    }


def ensure_empty_results() -> None:
    if not RESULTS_JSON.exists():
        RESULTS_JSON.write_text("{}", encoding="utf-8")
    if not RESULTS_CSV.exists():
        with RESULTS_CSV.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.writer(handle, delimiter=";")
            writer.writerow(
                [
                    "case_id",
                    "uid",
                    "name",
                    "date",
                    "decision",
                    "comment",
                    "model_status",
                    "model_thickness_m_median",
                    "updated_at",
                ]
            )


def main() -> None:
    data = build_data()
    OUT_JSON.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    ensure_empty_results()
    print(json.dumps(data["summary"], ensure_ascii=False, indent=2))
    print(f"Wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
