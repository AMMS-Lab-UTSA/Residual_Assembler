"""Synchronous real-browser request check; owns and always stops its server PID."""

import argparse
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
import urllib.error
import urllib.request

from playwright.sync_api import expect, sync_playwright


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    args.out.mkdir(parents=True, exist_ok=False)
    with socket.socket() as listener:
        listener.bind(("127.0.0.1", 0))
        port = listener.getsockname()[1]
    command = [sys.executable, "-m", "streamlit", "run", str(root / "scripts/app.py"),
               "--server.address=127.0.0.1", "--server.port=" + str(port),
               "--server.headless=true", "--browser.gatherUsageStats=false"]
    report = {"url": "http://127.0.0.1:%s" % port, "command": command}
    with (args.out / "server.log").open("w") as log:
        server = subprocess.Popen(command, cwd=root, env=os.environ.copy(), stdout=log, stderr=subprocess.STDOUT)
        report["server_pid"] = server.pid
        try:
            with sync_playwright() as playwright:
                browser = playwright.chromium.launch(headless=True)
                page = browser.new_page(viewport={"width": 1440, "height": 1100})
                for attempt in range(100):
                    if server.poll() is not None:
                        raise RuntimeError("Streamlit exited before readiness")
                    try:
                        with urllib.request.urlopen(report["url"] + "/_stcore/health", timeout=1) as response:
                            if response.status == 200:
                                break
                    except (urllib.error.URLError, TimeoutError):
                        page.wait_for_timeout(100)
                else:
                    raise RuntimeError("Streamlit readiness timeout")
                page.goto(report["url"])
                expect(page.get_by_role("tab", name="Sensitivity Request", exact=True)).to_have_attribute("aria-selected", "true")
                expect(page.get_by_role("heading", name="Residual Sensitivity Solver", exact=True)).to_be_visible()
                for name in ("OTI_UMAT.obj", "Analysis.inp", "Analysis.odb", "sensitivity_request.json"):
                    page.get_by_label(name + " path", exact=True).fill(str((args.inputs / name).resolve()))
                page.get_by_label("Output directory", exact=True).fill(str((args.out / "results").resolve()))
                page.get_by_role("button", name="Run sensitivity request", exact=True).click()
                expect(page.get_by_text("Executed: 4 scalar results. Independent validation: not run.", exact=True)).to_be_visible(timeout=90000)
                filenames = ["sensitivity_results.json", "sensitivity_tables.csv", "run_report.txt"]
                for filename in filenames:
                    button = page.get_by_role("button", name=filename, exact=True)
                    expect(button).to_be_visible()
                    with page.expect_download() as downloading:
                        button.click()
                    download = downloading.value
                    assert download.suggested_filename == filename
                    download.save_as(args.out / ("download_" + filename))
                page.screenshot(path=str(args.out / "desktop.png"), full_page=True)
                page.set_viewport_size({"width": 390, "height": 844})
                expect(page.get_by_role("heading", name="Residual Sensitivity Solver", exact=True)).to_be_visible()
                page.screenshot(path=str(args.out / "mobile.png"), full_page=True)
                report["horizontal_overflow"] = page.evaluate("document.documentElement.scrollWidth > window.innerWidth")
                assert not report["horizontal_overflow"]
                assert page.locator('[data-testid="stException"]').count() == 0
                report.update(passed=True, downloads=filenames, viewports=[[1440, 1100], [390, 844]])
                browser.close()
        finally:
            server.terminate()
            try:
                server.wait(timeout=15)
            except subprocess.TimeoutExpired:
                server.kill()
                server.wait(timeout=5)
            report["server_stopped"] = server.poll() is not None
            (args.out / "browser_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()