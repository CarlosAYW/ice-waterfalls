from __future__ import annotations

import argparse
import csv
import json
import socket
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse


HERE = Path(__file__).resolve().parent
VALIDATION_ROOT = HERE.parent
RESULTS_JSON = HERE / "validation_results.json"
RESULTS_CSV = HERE / "validation_results.csv"


CSV_FIELDS = [
    "case_id",
    "uid",
    "name",
    "date",
    "date_iso",
    "decision",
    "comment",
    "model_status",
    "model_thickness_m_median",
    "model_thickness_m_max",
    "model_climbability_median",
    "condition",
    "image_files",
    "video_files",
    "updated_at",
]


def write_results(payload: dict[str, object]) -> None:
    results = payload.get("results", payload)
    if not isinstance(results, dict):
        raise ValueError("Expected JSON object with a 'results' object")

    RESULTS_JSON.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")

    with RESULTS_CSV.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, delimiter=";")
        writer.writeheader()
        for key in sorted(results):
            row = results[key]
            if not isinstance(row, dict):
                continue
            writer.writerow({field: row.get(field, "") for field in CSV_FIELDS})


class ReviewHandler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(VALIDATION_ROOT), **kwargs)

    def log_message(self, format: str, *args) -> None:  # noqa: A002
        print(format % args)

    def do_GET(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        if parsed.path == "/":
            self.send_response(302)
            self.send_header("Location", "/Review_Tool/review_tool.html")
            self.end_headers()
            return
        if parsed.path == "/api/results":
            body = RESULTS_JSON.read_text(encoding="utf-8") if RESULTS_JSON.exists() else "{}"
            self.send_response(200)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.end_headers()
            self.wfile.write(body.encode("utf-8"))
            return
        super().do_GET()

    def do_POST(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        if parsed.path != "/api/save":
            self.send_error(404)
            return

        try:
            length = int(self.headers.get("Content-Length", "0"))
            raw = self.rfile.read(length)
            payload = json.loads(raw.decode("utf-8"))
            write_results(payload)
        except Exception as exc:  # pragma: no cover - local utility
            self.send_response(400)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.end_headers()
            self.wfile.write(json.dumps({"ok": False, "error": str(exc)}, ensure_ascii=False).encode("utf-8"))
            return

        self.send_response(200)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.end_headers()
        self.wfile.write(json.dumps({"ok": True}, ensure_ascii=False).encode("utf-8"))


def port_is_free(port: int) -> bool:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.settimeout(0.2)
        return sock.connect_ex(("127.0.0.1", port)) != 0


def find_port(start: int) -> int:
    for port in range(start, start + 100):
        if port_is_free(port):
            return port
    raise RuntimeError(f"No free port found from {start} to {start + 99}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()

    port = find_port(args.port)
    server = ThreadingHTTPServer(("127.0.0.1", port), ReviewHandler)
    print(f"Review tool running: http://127.0.0.1:{port}/Review_Tool/review_tool.html", flush=True)
    print(f"Saving to: {RESULTS_JSON}", flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
