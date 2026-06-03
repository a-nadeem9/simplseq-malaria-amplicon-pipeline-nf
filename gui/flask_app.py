"""Flask browser app for SIMPLseq-nf App."""

from __future__ import annotations

import csv
import datetime as dt
import os
import re
import signal
import subprocess
import sys
import threading
import time
import webbrowser
import zipfile
from io import BytesIO
from pathlib import Path
from typing import Any

from flask import Flask, abort, jsonify, render_template, request, send_file, send_from_directory, url_for

from simplseq import __version__
from simplseq.job_state import read_json, write_json
from simplseq.pathutils import user_path
from simplseq.progress import read_events
from simplseq.resources import human_bytes
from simplseq.runner import (
    STAGES,
    check_environment,
    local_runtime_env,
    progress_summary,
    project_root,
    results_manifest,
)
from simplseq.samplesheet import SAMPLE_FIELDS, FastqPair, FastqScan, scan_fastqs, write_samples_csv


RUN_PROCESSES: dict[str, subprocess.Popen[str]] = {}
RUN_LOCK = threading.Lock()
DOWNLOAD_SLUG_RE = re.compile(r"[^a-z0-9]+")
RUN_NAME_RE = re.compile(r"[^A-Za-z0-9._-]+")
ANSI_RE = re.compile(r"\x1B\][^\x07]*(?:\x07|\x1B\\)|\x1B\[[0-?]*[ -/]*[@-~]|\x1B[@-Z\\-_]")
CORE_RESULT_LABELS = {
    "ASV count table",
    "Mapped ASV table",
    "ASV to CIGAR map",
    "CIGAR count table",
}
REPORT_LABEL = "Run summary"
BUNDLE_RESULT_LABELS = CORE_RESULT_LABELS | {REPORT_LABEL, "Input FASTQ MD5s"}


def resolve_app_path(root: Path, value: str | os.PathLike[str] | None, default: str | Path) -> Path:
    raw = str(value if value not in {None, ""} else default).strip()
    path = user_path(raw)
    if not path.is_absolute():
        path = root / path
    return path.expanduser().resolve()


def rel_or_abs(root: Path, path: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return str(path)


def bool_payload(data: dict[str, Any], key: str, default: bool = False) -> bool:
    value = data.get(key, default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.lower() in {"1", "true", "yes", "on"}
    return bool(value)


def int_payload(data: dict[str, Any], key: str, default: int = 0) -> int:
    try:
        return int(data.get(key, default) or default)
    except (TypeError, ValueError):
        return default


def dinemites_model_settings(data: dict[str, Any]) -> dict[str, int | float | str | bool]:
    def require_int(key: str, default: int, minimum: int) -> int:
        try:
            value = int(data.get(key, default))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{key} must be at least {minimum}.") from exc
        if value < minimum:
            raise ValueError(f"{key} must be at least {minimum}.")
        return value

    def require_float_range(key: str, default: float, minimum: float, maximum: float) -> float:
        try:
            value = float(data.get(key, default))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{key} must be between {minimum} and {maximum}.") from exc
        if value < minimum or value > maximum:
            raise ValueError(f"{key} must be between {minimum} and {maximum}.")
        return value

    try:
        n_lags = int(data.get("n_lags", 3))
    except (TypeError, ValueError) as exc:
        raise ValueError("n_lags must be at least 1.") from exc
    if n_lags < 1:
        raise ValueError("n_lags must be at least 1.")

    try:
        min_abundance_pct = float(data.get("min_abundance_pct", 0.3))
    except (TypeError, ValueError) as exc:
        raise ValueError("min_abundance_pct must be between 0 and 100.") from exc
    if min_abundance_pct < 0 or min_abundance_pct > 100:
        raise ValueError("min_abundance_pct must be between 0 and 100.")
    abundance_denominator = str(data.get("abundance_denominator", "locus")).strip().lower()
    if abundance_denominator not in {"locus", "sample"}:
        raise ValueError("abundance_denominator must be either locus or sample.")

    seed = require_int("seed", 1, 1)
    refresh = require_int("refresh", 100, 0)
    bayesian_lag_days = require_int("bayesian_lag_days", 30, 1)
    bayesian_chains = require_int("bayesian_chains", 1, 1)
    bayesian_parallel_chains = require_int("bayesian_parallel_chains", 1, 1)
    bayesian_iter_warmup = require_int("bayesian_iter_warmup", 500, 1)
    bayesian_iter_sampling = require_int("bayesian_iter_sampling", 500, 1)
    bayesian_adapt_delta = require_float_range("bayesian_adapt_delta", 0.99, 0.000001, 0.999999)
    bayesian_drop_out = bool_payload(data, "bayesian_drop_out", False)

    common_settings: dict[str, int | float | str | bool] = {
        "min_abundance_pct": min_abundance_pct,
        "abundance_denominator": abundance_denominator,
        "seed": seed,
        "refresh": refresh,
        "bayesian_lag_days": bayesian_lag_days,
        "bayesian_chains": bayesian_chains,
        "bayesian_parallel_chains": bayesian_parallel_chains,
        "bayesian_iter_warmup": bayesian_iter_warmup,
        "bayesian_iter_sampling": bayesian_iter_sampling,
        "bayesian_adapt_delta": bayesian_adapt_delta,
        "bayesian_drop_out": bayesian_drop_out,
    }

    no_day_cutoff = bool_payload(data, "no_day_cutoff", False)
    raw_t_lag = str(data.get("t_lag", "Inf")).strip()
    if no_day_cutoff or raw_t_lag.lower() in {"", "inf", "infinity", "none"}:
        return {
            "n_lags": n_lags,
            "t_lag": "Inf",
            **common_settings,
        }

    try:
        numeric_t_lag = float(raw_t_lag)
    except ValueError as exc:
        raise ValueError("t_lag must be a non-negative number or Inf.") from exc
    if numeric_t_lag < 0:
        raise ValueError("t_lag must be a non-negative number or Inf.")

    return {
        "n_lags": n_lags,
        "t_lag": str(int(numeric_t_lag) if numeric_t_lag.is_integer() else numeric_t_lag),
        **common_settings,
    }


def dcifer_settings(data: dict[str, Any]) -> dict[str, int | float | str]:
    def require_int(key: str, default: int, minimum: int) -> int:
        try:
            value = int(data.get(key, default))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{key} must be at least {minimum}.") from exc
        if value < minimum:
            raise ValueError(f"{key} must be at least {minimum}.")
        return value

    try:
        min_abundance_pct = float(data.get("min_abundance_pct", 0.3))
    except (TypeError, ValueError) as exc:
        raise ValueError("min_abundance_pct must be between 0 and 100.") from exc
    if min_abundance_pct < 0 or min_abundance_pct > 100:
        raise ValueError("min_abundance_pct must be between 0 and 100.")

    abundance_denominator = str(data.get("abundance_denominator", "locus")).strip().lower()
    if abundance_denominator not in {"locus", "sample"}:
        raise ValueError("abundance_denominator must be either locus or sample.")

    try:
        alpha = float(data.get("alpha", 0.05))
    except (TypeError, ValueError) as exc:
        raise ValueError("alpha must be greater than 0 and less than 1.") from exc
    if alpha <= 0 or alpha >= 1:
        raise ValueError("alpha must be greater than 0 and less than 1.")

    afreq_mode = str(data.get("afreq_mode", "current_run")).strip().lower()
    if afreq_mode != "current_run":
        raise ValueError("afreq_mode currently supports only current_run.")

    return {
        "min_abundance_pct": min_abundance_pct,
        "abundance_denominator": abundance_denominator,
        "coi_lrank": require_int("coi_lrank", 2, 1),
        "ibd_grid_nr": require_int("ibd_grid_nr", 1000, 1),
        "alpha": alpha,
        "afreq_mode": afreq_mode,
    }


def slugify(label: str) -> str:
    slug = DOWNLOAD_SLUG_RE.sub("-", label.lower()).strip("-")
    return slug or "file"


def default_run_name(now: dt.datetime | None = None) -> str:
    stamp = (now or dt.datetime.now()).strftime("%Y-%m-%d_%H%M%S")
    return f"SIMPLseq_{stamp}"


def safe_run_name(value: str | None, *, now: dt.datetime | None = None) -> str:
    raw = str(value or "").strip() or default_run_name(now)
    name = RUN_NAME_RE.sub("_", raw).strip("._-")
    return name or default_run_name(now)


def allocate_run_outdir(
    parent: Path,
    run_name: str | None,
    *,
    now: dt.datetime | None = None,
    reuse_existing: bool = False,
) -> Path:
    parent.mkdir(parents=True, exist_ok=True)
    name = safe_run_name(run_name, now=now)
    candidate = parent / name
    if reuse_existing:
        candidate.mkdir(parents=True, exist_ok=True)
        return candidate.resolve()
    for index in range(1, 1000):
        unique = candidate if index == 1 else parent / f"{name}_{index:02d}"
        try:
            unique.mkdir(parents=True, exist_ok=False)
            return unique.resolve()
        except FileExistsError:
            continue
    raise RuntimeError(f"Could not allocate a unique run folder under {parent}")


def json_error(message: str, status: int = 400, **extra: Any):
    payload = {"ok": False, "error": message}
    payload.update(extra)
    return jsonify(payload), status


def sample_pair_json(pair: FastqPair, root: Path) -> dict[str, str]:
    return {
        "sample_id": pair.sample_id,
        "participant_id": pair.participant_id,
        "collection_date": pair.collection_date,
        "replicate": pair.replicate,
        "sample_type": pair.sample_type,
        "fastq_1": rel_or_abs(root, pair.fastq_1),
        "fastq_2": rel_or_abs(root, pair.fastq_2),
    }


def scan_json(scan: FastqScan, root: Path, *, preview_limit: int = 100) -> dict[str, Any]:
    missing_pairs = len(scan.missing_r2) + len(scan.orphan_r2)
    return {
        "fastq_dir": str(scan.fastq_dir),
        "pair_count": len(scan.pairs),
        "md5_files": scan.md5_files,
        "total_fastq_bytes": scan.total_fastq_bytes,
        "total_fastq_size": human_bytes(scan.total_fastq_bytes),
        "missing_pairs": missing_pairs,
        "missing_r2": scan.missing_r2[:100],
        "orphan_r2": scan.orphan_r2[:100],
        "duplicate_sample_ids": scan.duplicate_sample_ids,
        "preview": [sample_pair_json(pair, root) for pair in scan.pairs[:preview_limit]],
    }


def read_samples_preview(path: Path, limit: int = 100) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        rows = []
        for index, row in enumerate(csv.DictReader(handle)):
            if index >= limit:
                break
            rows.append({field: row.get(field, "") for field in SAMPLE_FIELDS})
        return rows


def file_tail(path: Path, max_bytes: int) -> tuple[str, bool]:
    if not path.exists():
        return "", False
    size = path.stat().st_size
    with path.open("rb") as handle:
        if size > max_bytes:
            handle.seek(size - max_bytes)
            data = handle.read()
            return data.decode("utf-8", errors="replace"), True
        return handle.read().decode("utf-8", errors="replace"), False


def clean_log_text(text: str) -> str:
    cleaned = ANSI_RE.sub("", text)
    return "".join(ch for ch in cleaned if ch in {"\n", "\r", "\t"} or ord(ch) >= 32).replace("\r", "\n")


def run_log_tail(outdir: Path, max_bytes: int) -> tuple[str, bool]:
    text, truncated = file_tail(outdir / "technical_log.txt", max_bytes)
    if text.strip():
        return text, truncated

    parts = []
    for path in (outdir / "logs" / "flask-run.stdout.log", outdir / "logs" / "flask-run.stderr.log"):
        log_text, log_truncated = file_tail(path, max_bytes)
        if log_text.strip():
            parts.append(log_text.rstrip())
        truncated = truncated or log_truncated
    return "\n".join(parts), truncated


def safe_exists(path: Path) -> bool:
    try:
        return path.exists()
    except OSError:
        return False


def safe_is_dir(path: Path) -> bool:
    try:
        return path.is_dir()
    except OSError:
        return False


def safe_is_file(path: Path) -> bool:
    try:
        return path.is_file()
    except OSError:
        return False


def safe_resolve(path: Path) -> Path | None:
    try:
        return path.expanduser().resolve()
    except (OSError, RuntimeError):
        return None


def is_wsl() -> bool:
    try:
        proc_version = Path("/proc/version").read_text(encoding="utf-8", errors="ignore").lower()
    except OSError:
        return False
    return "microsoft" in proc_version or "wsl" in proc_version


def path_style() -> str:
    if os.name == "nt":
        return "windows"
    if is_wsl():
        return "wsl"
    return "posix"


def is_macos() -> bool:
    return sys.platform == "darwin"


def wsl_to_windows_path(path: Path) -> str:
    text = str(path)
    match = re.match(r"^/mnt/([a-zA-Z])(?:/(.*))?$", text)
    if match:
        drive = match.group(1).upper()
        rest = (match.group(2) or "").replace("/", "\\")
        return f"{drive}:\\{rest}" if rest else f"{drive}:\\"
    try:
        completed = subprocess.run(
            ["wslpath", "-w", text],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
            timeout=3,
        )
    except (OSError, subprocess.TimeoutExpired):
        return text
    return completed.stdout.strip() or text


def windows_to_wsl_path(value: str) -> str:
    text = value.strip().replace("\\", "/")
    match = re.match(r"^([A-Za-z]):/?(.*)$", text)
    if not match:
        return value.strip()
    drive = match.group(1).lower()
    rest = match.group(2).strip("/")
    return f"/mnt/{drive}/{rest}" if rest else f"/mnt/{drive}"


def select_windows_folder_dialog(initial: Path | None = None) -> dict[str, Any]:
    env = os.environ.copy()
    if initial:
        env["SIMPLSEQ_PICKER_INITIAL"] = wsl_to_windows_path(initial)
    script = r"""
Add-Type -AssemblyName System.Windows.Forms
Add-Type -AssemblyName System.Drawing
[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new()
[System.Windows.Forms.Application]::EnableVisualStyles()

$owner = New-Object System.Windows.Forms.Form
$owner.Text = "SIMPLseq folder picker"
$owner.Size = New-Object System.Drawing.Size(1, 1)
$owner.StartPosition = "CenterScreen"
$owner.ShowInTaskbar = $false
$owner.TopMost = $true
$owner.Opacity = 0
$owner.Show()
$owner.Activate()

$dialog = New-Object System.Windows.Forms.FolderBrowserDialog
$dialog.Description = "Select the folder containing FASTQ files"
$dialog.ShowNewFolderButton = $false
if ($env:SIMPLSEQ_PICKER_INITIAL -and (Test-Path -LiteralPath $env:SIMPLSEQ_PICKER_INITIAL)) {
  $dialog.SelectedPath = $env:SIMPLSEQ_PICKER_INITIAL
}
$result = $dialog.ShowDialog($owner)
$owner.Dispose()
if ($result -eq [System.Windows.Forms.DialogResult]::OK) {
  Write-Output $dialog.SelectedPath
}
"""
    try:
        completed = subprocess.run(
            ["powershell.exe", "-NoLogo", "-NoProfile", "-Sta", "-ExecutionPolicy", "Bypass", "-Command", script],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=300,
            env=env,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"ok": False, "error": f"Folder picker could not open: {exc}"}
    selected = completed.stdout.strip().splitlines()[-1].strip() if completed.stdout.strip() else ""
    if completed.returncode != 0 and not selected:
        detail = " ".join(completed.stderr.split())[:220] if completed.stderr else "PowerShell folder picker failed"
        return {"ok": False, "error": detail}
    if not selected:
        return {"ok": True, "selected": False}
    return {"ok": True, "selected": True, "path": windows_to_wsl_path(selected), "windows_path": selected}


def select_macos_folder_dialog(initial: Path | None = None) -> dict[str, Any]:
    initial_path = initial.as_posix() if initial and safe_exists(initial) else ""
    script = r"""
on run argv
  set promptText to "Select the folder containing FASTQ files"
  if (count of argv) > 0 and item 1 of argv is not "" then
    try
      set initialAlias to POSIX file (item 1 of argv) as alias
      set chosenFolder to choose folder with prompt promptText default location initialAlias
    on error
      set chosenFolder to choose folder with prompt promptText
    end try
  else
    set chosenFolder to choose folder with prompt promptText
  end if
  return POSIX path of chosenFolder
end run
"""
    try:
        completed = subprocess.run(
            ["osascript", "-e", script, "--", initial_path],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=300,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"ok": False, "error": f"Folder picker could not open: {exc}"}
    selected = completed.stdout.strip().splitlines()[-1].strip() if completed.stdout.strip() else ""
    if completed.returncode != 0 and not selected:
        detail = " ".join(completed.stderr.split())
        if "User canceled" in detail or "(-128)" in detail:
            return {"ok": True, "selected": False}
        return {"ok": False, "error": detail[:220] or "macOS folder picker failed"}
    if not selected:
        return {"ok": True, "selected": False}
    path = selected.rstrip("/") or "/"
    return {"ok": True, "selected": True, "path": path}


def select_folder_dialog(initial: Path | None = None) -> dict[str, Any]:
    if is_macos():
        return select_macos_folder_dialog(initial)
    if is_wsl():
        return select_windows_folder_dialog(initial)
    return {"ok": False, "error": "Native folder picker is only available in WSL or macOS for this build."}


def common_paths(workspace_root: Path, app_root: Path) -> list[dict[str, str]]:
    paths: list[tuple[str, Path]] = [
        ("Current folder", workspace_root),
        ("Data in current folder", workspace_root / "data"),
        ("Home", Path.home()),
    ]
    if safe_exists(workspace_root / "test-data"):
        paths.append(("Test data", workspace_root / "test-data"))
    desktop = Path.home() / "Desktop"
    if safe_exists(desktop):
        paths.append(("Desktop", desktop))
    windows_home = Path("/mnt/c/Users") / Path.home().name
    for label, path in [
        ("Windows home", windows_home),
        ("Windows Desktop", windows_home / "Desktop"),
        ("Windows Downloads", windows_home / "Downloads"),
        ("Windows Documents", windows_home / "Documents"),
    ]:
        if safe_exists(path):
            paths.append((label, path))
    for mount in [Path("/mnt/c"), Path("/mnt/d")]:
        if safe_exists(mount):
            paths.append((str(mount), mount))
    windows_users = Path("/mnt/c/Users")
    if safe_exists(windows_users):
        try:
            candidates = sorted(windows_users.iterdir(), key=lambda item: item.name.lower())
        except OSError:
            candidates = []
        for candidate in candidates:
            desktop_dir = candidate / "Desktop"
            if safe_exists(desktop_dir):
                paths.append((f"{candidate.name} Desktop", desktop_dir))
    seen: set[str] = set()
    result = []
    for label, path in paths:
        resolved_path = safe_resolve(path)
        if resolved_path is None:
            continue
        resolved = str(resolved_path)
        if resolved in seen:
            continue
        seen.add(resolved)
        result.append({"label": label, "path": resolved})
    return result


def fastq_count(path: Path) -> int:
    count = 0
    try:
        for item in path.iterdir():
            if item.name.endswith((".fastq.gz", ".fq.gz")) and safe_is_file(item):
                count += 1
    except OSError:
        return 0
    return count


def browse_payload(path: Path) -> dict[str, Any]:
    if not safe_exists(path):
        return {"path": str(path), "exists": False, "directories": [], "parent": str(path.parent)}
    if not safe_is_dir(path):
        return {"path": str(path), "exists": True, "is_dir": False, "directories": [], "parent": str(path.parent)}
    directories = []
    try:
        children = sorted((item for item in path.iterdir() if safe_is_dir(item)), key=lambda item: item.name.lower())
    except OSError:
        children = []
    for child in children[:250]:
        directories.append(
            {
                "name": child.name,
                "path": str(child),
                "fastq_files": fastq_count(child),
            }
        )
    current_fastq_files = fastq_count(path)
    scan = scan_fastqs(path) if current_fastq_files else FastqScan(path, [], [], [], 0, 0, [])
    return {
        "path": str(path),
        "exists": True,
        "is_dir": True,
        "parent": str(path.parent),
        "directories": directories,
        "fastq_files": current_fastq_files,
        "pair_count": len(scan.pairs),
        "missing_pairs": len(scan.missing_r2) + len(scan.orphan_r2),
    }


def result_files_with_downloads(root: Path, outdir: Path) -> dict[str, Any]:
    manifest = results_manifest(outdir)
    files = []
    for item in manifest["files"]:
        label = str(item["label"])
        slug = slugify(label)
        exists = bool(item["exists"])
        size_bytes = int(item["size_bytes"])
        files.append(
            {
                "label": label,
                "slug": slug,
                "path": item["path"],
                "relative_path": rel_or_abs(root, Path(str(item["path"]))),
                "exists": exists,
                "size_bytes": size_bytes,
                "size": human_bytes(size_bytes),
                "status": "ready" if exists and size_bytes else "missing",
                "download_url": url_for("download_result", file_key=slug, out=str(outdir)) if exists else "",
                "view_url": url_for("download_result", file_key=slug, out=str(outdir), inline=1)
                if exists and label == REPORT_LABEL
                else "",
            }
        )
    manifest["files"] = files
    manifest["report"] = next((item for item in files if item["label"] == REPORT_LABEL), None)
    manifest["core_files"] = [item for item in files if item["label"] in CORE_RESULT_LABELS]
    manifest["support_files"] = [
        item for item in files if item["label"] != REPORT_LABEL and item["label"] not in CORE_RESULT_LABELS
    ]
    manifest["ready_counts"] = {
        "core": sum(1 for item in manifest["core_files"] if item["status"] == "ready"),
        "support": sum(1 for item in manifest["support_files"] if item["status"] == "ready"),
    }
    manifest["bundle_ready"] = any(
        item["label"] in BUNDLE_RESULT_LABELS and item["status"] == "ready" for item in files
    )
    manifest["bundle_url"] = url_for("download_bundle", out=str(outdir)) if manifest["bundle_ready"] else ""
    return manifest


def bundle_result_paths(outdir: Path) -> list[tuple[Path, str]]:
    outdir = outdir.resolve()
    bundled: list[tuple[Path, str]] = []
    for item in results_manifest(outdir)["files"]:
        label = str(item["label"])
        if label not in BUNDLE_RESULT_LABELS:
            continue
        path = Path(str(item["path"])).resolve()
        try:
            arcname = path.relative_to(outdir).as_posix()
        except ValueError:
            continue
        if path.exists() and path.is_file() and path.stat().st_size > 0:
            bundled.append((path, arcname))
    return bundled


def active_state(outdir: Path) -> bool:
    state = read_json(outdir / "run_state.json")
    return state.get("status") in {"starting", "running"}


def process_active(outdir: Path) -> bool:
    key = str(outdir)
    with RUN_LOCK:
        process = RUN_PROCESSES.get(key)
        if process is None:
            return False
        if process.poll() is None:
            return True
        RUN_PROCESSES.pop(key, None)
    return False


def active_process_outdir() -> Path | None:
    with RUN_LOCK:
        for key, process in list(RUN_PROCESSES.items()):
            if process.poll() is None:
                return Path(key)
            RUN_PROCESSES.pop(key, None)
    return None


def any_process_active() -> bool:
    return active_process_outdir() is not None


def terminate_process_tree(process: subprocess.Popen[str], timeout: float = 5.0) -> None:
    if process.poll() is not None:
        return
    try:
        if os.name == "nt":
            process.terminate()
        else:
            os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=timeout)
        return
    except subprocess.TimeoutExpired:
        try:
            if os.name == "nt":
                process.kill()
            else:
                os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            return
        process.wait(timeout=timeout)


def stop_tracked_process(outdir: Path) -> bool:
    key = str(outdir)
    with RUN_LOCK:
        process = RUN_PROCESSES.get(key)
        if process is None:
            return False
        if process.poll() is not None:
            RUN_PROCESSES.pop(key, None)
            return False

    terminate_process_tree(process)

    with RUN_LOCK:
        RUN_PROCESSES.pop(key, None)

    state_file = outdir / "run_state.json"
    state = read_json(state_file)
    state.update(
        {
            "status": "stopped",
            "completed_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "outdir": str(outdir),
            "detail": "Run stopped by user from the SIMPLseq browser app.",
        }
    )
    write_json(state_file, state)
    return True



def payload_resume_enabled(data: dict[str, Any]) -> bool:
    if "resume" in data:
        return bool_payload(data, "resume", False)
    if "clean" in data:
        return not bool_payload(data, "clean", True)
    return False


def headless_run_command(samples: Path, outdir: Path, data: dict[str, Any]) -> list[str]:
    command = [
        sys.executable,
        "-m",
        "simplseq",
        "run-headless",
        "--samples",
        str(samples),
        "--out",
        str(outdir),
    ]
    cpus = int_payload(data, "cpus", 0)
    memory = str(data.get("memory", "")).strip()
    if cpus:
        command.extend(["--cpus", str(cpus)])
    if memory:
        command.extend(["--memory", memory])
    if not payload_resume_enabled(data):
        command.append("--no-resume")
    if bool_payload(data, "dry_run", False):
        command.append("--dry-run")
    return command


def start_run_process(root: Path, samples: Path, outdir: Path, data: dict[str, Any]) -> subprocess.Popen[str]:
    outdir.mkdir(parents=True, exist_ok=True)
    logs_dir = outdir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    command = headless_run_command(samples, outdir, data)
    env = local_runtime_env(root)
    env["SIMPLSEQ_PROJECT_ROOT"] = str(root)
    env["PYTHONPATH"] = str(root / "src") + os.pathsep + env.get("PYTHONPATH", "")
    stdout_handle = (logs_dir / "flask-run.stdout.log").open("w", encoding="utf-8")
    stderr_handle = (logs_dir / "flask-run.stderr.log").open("w", encoding="utf-8")
    kwargs: dict[str, Any] = {}
    if os.name == "nt":
        kwargs["creationflags"] = getattr(subprocess, "CREATE_NO_WINDOW", 0)
    else:
        kwargs["start_new_session"] = True
    process = subprocess.Popen(
        command,
        cwd=root,
        env=env,
        text=True,
        stdout=stdout_handle,
        stderr=stderr_handle,
        **kwargs,
    )
    stdout_handle.close()
    stderr_handle.close()
    with RUN_LOCK:
        RUN_PROCESSES[str(outdir)] = process
    return process


def create_app(root: Path | None = None, workspace_root: Path | None = None) -> Flask:
    app_root = (root or project_root()).resolve()
    workspace = (workspace_root or Path.cwd()).expanduser().resolve()
    gui_root = Path(__file__).resolve().parent
    app = Flask(
        __name__,
        static_folder=str(gui_root / "static"),
        template_folder=str(gui_root / "templates"),
    )
    app.config["SIMPLSEQ_PROJECT_ROOT"] = app_root
    app.config["SIMPLSEQ_WORKSPACE_ROOT"] = workspace
    app.config["TEMPLATES_AUTO_RELOAD"] = True

    @app.get("/")
    def index():
        static_paths = [gui_root / "static" / "css" / "app.css", gui_root / "static" / "js" / "app.js"]
        asset_version = str(max((path.stat().st_mtime_ns for path in static_paths if path.exists()), default=0))
        return render_template(
            "index.html",
            workspace_root=str(workspace),
            asset_version=asset_version,
        )

    @app.get("/assets/<path:filename>")
    def assets(filename: str):
        assets_dir = (app_root / "assets").resolve()
        requested = (assets_dir / filename).resolve()
        try:
            requested.relative_to(assets_dir)
        except ValueError:
            abort(404)
        return send_from_directory(assets_dir, filename)

    @app.get("/api/health")
    def api_health():
        return jsonify(
            {
                "ok": True,
                "app": "SIMPLseq-nf App",
                "version": f"v{__version__}-dev" if __version__ == "0.1.0" else f"v{__version__}",
                "app_root": str(app_root),
                "workspace_root": str(workspace),
                "path_style": path_style(),
                "common_paths": common_paths(workspace, app_root),
            }
        )

    @app.get("/api/browse")
    def api_browse():
        path = resolve_app_path(workspace, request.args.get("path"), workspace)
        return jsonify({"ok": True, **browse_payload(path)})

    @app.post("/api/select-folder")
    def api_select_folder():
        data = request.get_json(silent=True) or {}
        initial = resolve_app_path(workspace, data.get("initial"), workspace)
        return jsonify(select_folder_dialog(initial))

    @app.post("/api/scan")
    def api_scan():
        data = request.get_json(silent=True) or {}
        fastq_dir = resolve_app_path(workspace, data.get("fastq_dir"), "data")
        samples_out = resolve_app_path(workspace, data.get("samples_out"), "samples.csv")
        include_pool = bool_payload(data, "include_pool_in_sample_id", False)
        absolute = bool_payload(data, "absolute_paths", True)
        write_samples = bool_payload(data, "write_samples", True)

        scan = scan_fastqs(fastq_dir, include_pool_in_sample_id=include_pool)
        written = False
        duplicates: list[str] = scan.duplicate_sample_ids
        count = 0
        if write_samples and not duplicates:
            count, duplicates = write_samples_csv(
                fastq_dir,
                samples_out,
                include_pool_in_sample_id=include_pool,
                absolute=absolute,
            )
            written = not duplicates
        response = scan_json(scan, workspace)
        response.update(
            {
                "ok": True,
                "samples_out": str(samples_out),
                "samples_relative": rel_or_abs(workspace, samples_out),
                "samples_written": written,
                "sample_rows_written": count,
                "sample_preview": read_samples_preview(samples_out),
            }
        )
        return jsonify(response)

    @app.post("/api/check")
    def api_check():
        data = request.get_json(silent=True) or {}
        samples = data.get("samples")
        samples_path = resolve_app_path(workspace, samples, "samples.csv") if samples else None
        outdir = resolve_app_path(workspace, data.get("outdir"), "results")
        rows = check_environment(app_root, samples_path, outdir=outdir)
        failed = sum(1 for row in rows if row.get("status") not in {"ok", "warn"})
        return jsonify({"ok": failed == 0, "failed": failed, "checks": rows})

    @app.post("/api/run")
    def api_run():
        data = request.get_json(silent=True) or {}
        samples = resolve_app_path(workspace, data.get("samples"), "samples.csv")
        output_parent = resolve_app_path(workspace, data.get("outdir"), "results")
        resume = payload_resume_enabled(data)
        outdir = allocate_run_outdir(output_parent, str(data.get("run_name", "")), reuse_existing=resume)
        dry_run = bool_payload(data, "dry_run", False)
        if not samples.exists() and not dry_run:
            return json_error(f"Sample sheet not found: {samples}", 400)
        if any_process_active():
            return json_error("A SIMPLseq run is already active.", 409, outdir=str(outdir))
        try:
            process = start_run_process(app_root, samples, outdir, data)
        except Exception as exc:  # pragma: no cover - reported to browser
            return json_error(str(exc), 500)
        return jsonify(
            {
                "ok": True,
                "pid": process.pid,
                "outdir": str(outdir),
                "output_parent": str(output_parent),
                "samples": str(samples),
                "dry_run": dry_run,
                "status_url": url_for("api_status", out=str(outdir)),
            }
        )

    @app.post("/api/stop-run")
    def api_stop_run():
        data = request.get_json(silent=True) or {}
        outdir = resolve_app_path(workspace, data.get("outdir") or data.get("out"), "results")
        stopped = stop_tracked_process(outdir)
        if not stopped:
            return json_error("No active SIMPLseq run is attached to this app session.", 409, outdir=str(outdir))
        state = read_json(outdir / "run_state.json")
        summary = progress_summary(outdir)
        summary = dict(summary)
        summary["status"] = "stopped"
        return jsonify(
            {
                "ok": True,
                "stopped": True,
                "active": False,
                "outdir": str(outdir),
                "state": state,
                "summary": summary,
            }
        )

    @app.get("/api/active-run")
    def api_active_run():
        outdir = active_process_outdir()
        if outdir is None:
            return jsonify({"ok": True, "active": False})
        return jsonify(
            {
                "ok": True,
                "active": True,
                "outdir": str(outdir),
                "status_url": url_for("api_status", out=str(outdir)),
            }
        )

    @app.get("/api/status")
    def api_status():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        state = read_json(outdir / "run_state.json")
        summary = progress_summary(outdir)
        tracked_active = process_active(outdir)
        if tracked_active:
            state = dict(state)
            summary = dict(summary)
            status = str(state.get("status") or summary.get("status") or "pending")
            if status == "pending":
                status = "starting"
            state.setdefault("status", status)
            summary["status"] = status
            if str(summary.get("current_stage", "")) not in STAGES:
                summary["current_stage"] = "pending"
        if active_state(outdir) and not tracked_active:
            state = dict(state)
            state["status"] = "stale"
            state.setdefault("detail", "The previous run state is not attached to this app session.")
            summary = dict(summary)
            summary["status"] = "stale"
        return jsonify(
            {
                "ok": True,
                "outdir": str(outdir),
                "state": state,
                "summary": summary,
                "active": tracked_active,
            }
        )

    @app.get("/api/progress")
    def api_progress():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        events = read_events(outdir / "progress.jsonl")
        return jsonify(
            {
                "ok": True,
                "outdir": str(outdir),
                "events": events,
                "summary": progress_summary(outdir),
            }
        )

    @app.get("/api/results")
    def api_results():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        manifest = result_files_with_downloads(workspace, outdir)
        return jsonify({"ok": True, **manifest})

    @app.get("/api/logs")
    def api_logs():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        max_bytes = max(1000, min(int_payload(request.args, "max_bytes", 50000), 250000))
        log_text, truncated = run_log_tail(outdir, max_bytes)
        return jsonify(
            {
                "ok": True,
                "outdir": str(outdir),
                "path": str(outdir / "technical_log.txt"),
                "text": clean_log_text(log_text),
                "truncated": truncated,
            }
        )

    @app.get("/download/<file_key>")
    def download_result(file_key: str):
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        manifest = results_manifest(outdir)
        for item in manifest["files"]:
            label = str(item["label"])
            if slugify(label) != file_key:
                continue
            path = Path(str(item["path"])).resolve()
            try:
                path.relative_to(outdir.resolve())
            except ValueError:
                abort(404)
            if not path.exists() or not path.is_file():
                abort(404)
            inline = request.args.get("inline") == "1" and label == REPORT_LABEL
            return send_file(path, as_attachment=not inline, download_name=path.name)
        abort(404)

    @app.get("/download-bundle")
    def download_bundle():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        bundle_paths = bundle_result_paths(outdir)
        if not bundle_paths:
            abort(404)
        archive = BytesIO()
        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as zip_handle:
            for path, arcname in bundle_paths:
                zip_handle.write(path, arcname)
        archive.seek(0)
        bundle_name = f"{outdir.name or 'simplseq-results'}-output-bundle.zip"
        return send_file(archive, as_attachment=True, download_name=bundle_name, mimetype="application/zip")

    # ------------------------------------------------------------------
    # DINEMITES analysis endpoints
    # ------------------------------------------------------------------

    DINEMITES_THREADS: dict[str, threading.Thread] = {}
    DINEMITES_LOCK = threading.Lock()

    DINEMITES_FILE_KEYS = {
        "allele_probabilities": "dinemites_allele_probabilities.tsv",
        "allele_key": "dinemites_allele_key.tsv",
        "molfoi": "dinemites_molfoi.tsv",
        "new_infections": "dinemites_new_infections.tsv",
    }

    def dinemites_outdir(outdir: Path) -> Path:
        return outdir / "dinemites"

    def dinemites_plots_dir(outdir: Path) -> Path:
        return dinemites_outdir(outdir) / "dinemites_plots"

    def dinemites_state_path(outdir: Path) -> Path:
        return dinemites_outdir(outdir) / "dinemites_state.json"

    def dinemites_row_value(row: dict[str, str], *keys: str) -> str:
        for key in keys:
            value = row.get(key, "")
            if value not in {"", None}:
                return value
        return ""

    def dinemites_tsv_rows(path: Path, limit: int = 200) -> list[dict[str, str]]:
        if not safe_is_file(path):
            return []
        rows: list[dict[str, str]] = []
        try:
            with path.open(newline="", encoding="utf-8-sig", errors="replace") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for index, row in enumerate(reader):
                    if index >= limit:
                        break
                    rows.append({
                        str(key or "").lstrip("\ufeff"): str(value or "")
                        for key, value in row.items()
                    })
        except Exception:
            return []
        return rows

    def dinemites_plot_subject(filename: str) -> str:
        stem = Path(filename).stem
        return stem[len("subject_") :] if stem.startswith("subject_") else stem

    def dinemites_plot_entries(outdir: Path) -> list[dict[str, Any]]:
        plots_dir = dinemites_plots_dir(outdir)
        if not safe_is_dir(plots_dir):
            return []
        entries: list[dict[str, Any]] = []
        try:
            paths = sorted(plots_dir.glob("*.png"), key=lambda item: item.name.lower())
        except OSError:
            return []
        for path in paths:
            if not safe_is_file(path):
                continue
            size_bytes = path.stat().st_size
            filename = path.name
            entries.append({
                "filename": filename,
                "subject": dinemites_plot_subject(filename),
                "exists": True,
                "size_bytes": size_bytes,
                "size": human_bytes(size_bytes),
                "view_url": url_for(
                    "download_dinemites_plot",
                    filename=filename,
                    out=str(outdir),
                    inline=1,
                ),
                "download_url": url_for(
                    "download_dinemites_plot",
                    filename=filename,
                    out=str(outdir),
                ),
            })
        return entries

    def dinemites_process_active(outdir: Path) -> bool:
        key = str(outdir)
        with DINEMITES_LOCK:
            thread = DINEMITES_THREADS.get(key)
            if thread is None:
                return False
            if thread.is_alive():
                return True
            DINEMITES_THREADS.pop(key, None)
        return False

    def run_dinemites_worker(
        root: Path,
        outdir: Path,
        samples: Path,
        model_type: str,
        model_settings: dict[str, int | float | str | bool],
    ) -> None:
        dm_dir = dinemites_outdir(outdir)
        dm_dir.mkdir(parents=True, exist_ok=True)
        state_file = dinemites_state_path(outdir)
        cigar_path = outdir / "run_dada2" / "seqtab_cigar.tsv"
        dm_input = dm_dir / "dinemites_input.tsv"

        def state_payload(status: str, **extra: Any) -> dict[str, Any]:
            payload = {
                "status": status,
                "model": model_type,
                **model_settings,
                "outdir": str(dm_dir),
            }
            payload.update(extra)
            return payload

        write_json(state_file, state_payload(
            "running",
            started_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        ))

        log_path = dm_dir / "dinemites.log"
        env = os.environ.copy()
        env["SIMPLSEQ_PROJECT_ROOT"] = str(root)

        try:
            # Step 1: Convert SIMPLseq output to DINEMITES input
            cmd_convert = [
                "Rscript",
                str(root / "workflow" / "scripts" / "simplseq_to_dinemites.R"),
                "--cigar", str(cigar_path),
                "--samples", str(samples),
                "--out", str(dm_input),
                "--min_abundance_pct", str(model_settings["min_abundance_pct"]),
                "--abundance_denominator", str(model_settings["abundance_denominator"]),
            ]
            kwargs: dict[str, Any] = {}
            if os.name == "nt":
                kwargs["creationflags"] = getattr(subprocess, "CREATE_NO_WINDOW", 0)
            with log_path.open("w", encoding="utf-8") as log_handle:
                result = subprocess.run(
                    cmd_convert,
                    cwd=root,
                    env=env,
                    text=True,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    check=False,
                    timeout=600,
                    **kwargs,
                )
            if result.returncode != 0:
                write_json(state_file, state_payload(
                    "failed",
                    detail="simplseq_to_dinemites.R failed. Check dinemites.log.",
                    completed_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                ))
                return

            # Step 2: Run DINEMITES analysis
            cmd_run = [
                "Rscript",
                str(root / "workflow" / "scripts" / "run_dinemites.R"),
                "--input", str(dm_input),
                "--model", model_type,
                "--outdir", str(dm_dir),
                "--n_lags", str(model_settings["n_lags"]),
                "--t_lag", str(model_settings["t_lag"]),
                "--seed", str(model_settings["seed"]),
                "--refresh", str(model_settings["refresh"]),
                "--bayesian_lag_days", str(model_settings["bayesian_lag_days"]),
                "--bayesian_chains", str(model_settings["bayesian_chains"]),
                "--bayesian_parallel_chains", str(model_settings["bayesian_parallel_chains"]),
                "--bayesian_iter_warmup", str(model_settings["bayesian_iter_warmup"]),
                "--bayesian_iter_sampling", str(model_settings["bayesian_iter_sampling"]),
                "--bayesian_adapt_delta", str(model_settings["bayesian_adapt_delta"]),
                "--bayesian_drop_out", str(model_settings["bayesian_drop_out"]).lower(),
            ]
            with log_path.open("a", encoding="utf-8") as log_handle:
                result = subprocess.run(
                    cmd_run,
                    cwd=root,
                    env=env,
                    text=True,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    check=False,
                    timeout=7200,
                    **kwargs,
                )
            if result.returncode != 0:
                write_json(state_file, state_payload(
                    "failed",
                    detail="run_dinemites.R failed. Check dinemites.log.",
                    completed_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                ))
                return

            write_json(state_file, state_payload(
                "complete",
                completed_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            ))

        except Exception as exc:
            write_json(state_file, state_payload(
                "failed",
                detail=str(exc),
                completed_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            ))
        finally:
            with DINEMITES_LOCK:
                DINEMITES_THREADS.pop(str(outdir), None)

    @app.post("/api/dinemites/run")
    def api_dinemites_run():
        data = request.get_json(silent=True) or {}
        model_type = str(data.get("model_type", "simple")).strip()
        if model_type not in {"simple", "clustering", "bayesian"}:
            return json_error(f"Invalid DINEMITES model type: {model_type}", 400)
        try:
            model_settings = dinemites_model_settings(data)
        except ValueError as exc:
            return json_error(str(exc), 400)

        outdir = resolve_app_path(workspace, data.get("outdir"), "results")
        samples = resolve_app_path(workspace, data.get("samples"), "samples.csv")

        cigar_path = outdir / "run_dada2" / "seqtab_cigar.tsv"
        if not cigar_path.exists():
            return json_error(
                "Pipeline results not found. Run the main SIMPLseq pipeline first "
                "(seqtab_cigar.tsv is required).",
                400,
            )
        if not samples.exists():
            return json_error(f"Sample sheet not found: {samples}", 400)

        if dinemites_process_active(outdir):
            return json_error("A DINEMITES analysis is already running for this output folder.", 409)

        thread = threading.Thread(
            target=run_dinemites_worker,
            args=(app_root, outdir, samples, model_type, model_settings),
            daemon=True,
        )
        with DINEMITES_LOCK:
            DINEMITES_THREADS[str(outdir)] = thread
        thread.start()

        return jsonify({
            "ok": True,
            "status": "running",
            "model": model_type,
            **model_settings,
            "outdir": str(outdir),
        })

    @app.get("/api/dinemites/status")
    def api_dinemites_status():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        state = read_json(dinemites_state_path(outdir))
        active = dinemites_process_active(outdir)
        status = state.get("status", "idle")
        if active and status not in {"running"}:
            status = "running"
        return jsonify({
            "ok": True,
            "status": status,
            "active": active,
            "state": state,
            "outdir": str(outdir),
        })

    @app.get("/api/dinemites/results")
    def api_dinemites_results():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        dm_dir = dinemites_outdir(outdir)
        state = read_json(dinemites_state_path(outdir))

        files: dict[str, dict[str, Any]] = {}
        for key, filename in DINEMITES_FILE_KEYS.items():
            path = dm_dir / filename
            exists = safe_is_file(path)
            size_bytes = path.stat().st_size if exists else 0
            files[key] = {
                "exists": exists,
                "size_bytes": size_bytes,
                "size": human_bytes(size_bytes),
                "download_url": url_for("download_dinemites", file_key=key, out=str(outdir)) if exists else "",
            }

        plots = dinemites_plot_entries(outdir)
        has_outputs = any(item["exists"] for item in files.values()) or bool(plots)
        if state.get("status") in {None, "", "idle"} and has_outputs:
            state = {
                **state,
                "status": "complete",
                "detail": "Existing DINEMITES outputs found.",
                "outdir": str(dm_dir),
            }

        # Read summary from results files if available
        summary: dict[str, Any] = {"new_infections": "--", "molfoi": "--", "subjects": "--"}
        new_inf_path = dm_dir / "dinemites_new_infections.tsv"
        molfoi_path = dm_dir / "dinemites_molfoi.tsv"
        subjects_data: list[dict[str, str]] = []
        molfoi_by_subject: dict[str, str] = {}
        allele_key_rows: list[dict[str, str]] = []

        for row_item in dinemites_tsv_rows(dm_dir / "dinemites_allele_key.tsv", limit=500):
            allele_key_rows.append({
                "short_allele_id": dinemites_row_value(row_item, "short_allele_id", "short_id", "allele_id"),
                "locus": dinemites_row_value(row_item, "locus"),
                "allele": dinemites_row_value(row_item, "allele"),
            })

        if safe_is_file(molfoi_path):
            try:
                with molfoi_path.open(newline="", encoding="utf-8-sig", errors="replace") as fh:
                    reader = csv.DictReader(fh, delimiter="\t")
                    rows = list(reader)
                    vals = []
                    for r in rows:
                        subject = dinemites_row_value(r, "subject", "participant_id")
                        molfoi_value = dinemites_row_value(r, "molFOI", "molfoi", "mol_foi")
                        if subject:
                            molfoi_by_subject[subject] = molfoi_value
                        try:
                            vals.append(float(molfoi_value))
                        except (ValueError, TypeError):
                            pass
                    if vals:
                        summary["molfoi"] = round(sum(vals) / len(vals), 3)
                    if rows:
                        summary["subjects"] = len({dinemites_row_value(r, "subject", "participant_id") for r in rows})
            except Exception:
                pass

        if safe_is_file(new_inf_path):
            try:
                with new_inf_path.open(newline="", encoding="utf-8-sig", errors="replace") as fh:
                    reader = csv.DictReader(fh, delimiter="\t")
                    rows = list(reader)
                    new_infection_values = []
                    summary["subjects"] = len({dinemites_row_value(r, "subject", "participant_id") for r in rows})
                    for r in rows:
                        new_value = dinemites_row_value(r, "new_infections", "n_new")
                        try:
                            new_infection_values.append(float(new_value))
                        except (ValueError, TypeError):
                            pass
                    for r in rows[:200]:
                        subject = dinemites_row_value(r, "subject", "participant_id")
                        new_value = dinemites_row_value(r, "new_infections", "n_new")
                        subjects_data.append({
                            "subject": subject,
                            "new_infections": new_value,
                            "molfoi": molfoi_by_subject.get(subject, dinemites_row_value(r, "molFOI", "molfoi", "mol_foi")),
                            "time_points": dinemites_row_value(r, "time_points", "n_timepoints"),
                        })
                    if new_infection_values:
                        total_new = sum(new_infection_values)
                        summary["new_infections"] = int(total_new) if total_new.is_integer() else round(total_new, 3)
            except Exception:
                pass

        if not subjects_data and molfoi_by_subject:
            for subject, molfoi_value in list(molfoi_by_subject.items())[:200]:
                subjects_data.append({
                    "subject": subject,
                    "new_infections": "",
                    "molfoi": molfoi_value,
                    "time_points": "",
                })

        return jsonify({
            "ok": True,
            "state": state,
            "files": files,
            "plots": plots,
            "summary": summary,
            "subjects": subjects_data,
            "allele_key": allele_key_rows,
            "outdir": str(outdir),
        })

    @app.get("/download/dinemites/<file_key>")
    def download_dinemites(file_key: str):
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        filename = DINEMITES_FILE_KEYS.get(file_key)
        if not filename:
            abort(404)
        path = (dinemites_outdir(outdir) / filename).resolve()
        try:
            path.relative_to(dinemites_outdir(outdir).resolve())
        except ValueError:
            abort(404)
        if not path.exists() or not path.is_file():
            abort(404)
        return send_file(path, as_attachment=True, download_name=path.name)

    @app.get("/download/dinemites-plot/<path:filename>")
    def download_dinemites_plot(filename: str):
        if Path(filename).name != filename or "\\" in filename:
            abort(404)
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        plots_dir = dinemites_plots_dir(outdir).resolve()
        path = (plots_dir / filename).resolve()
        try:
            path.relative_to(plots_dir)
        except ValueError:
            abort(404)
        if path.suffix.lower() != ".png" or not path.exists() or not path.is_file():
            abort(404)
        inline = request.args.get("inline") == "1"
        return send_file(path, as_attachment=not inline, download_name=path.name, mimetype="image/png")

    # ------------------------------------------------------------------
    # dcifer analysis endpoints
    # ------------------------------------------------------------------

    DCIFER_THREADS: dict[str, threading.Thread] = {}
    DCIFER_LOCK = threading.Lock()

    DCIFER_FILE_KEYS = {
        "input": "dcifer_input_long.tsv",
        "filter_summary": "dcifer_filter_summary.tsv",
        "replicate_summary": "dcifer_replicate_summary.tsv",
        "coi": "dcifer_coi.tsv",
        "allele_frequencies": "dcifer_allele_frequencies.tsv",
        "pairwise_relatedness": "dcifer_pairwise_relatedness.tsv",
        "relatedness_matrix": "dcifer_relatedness_matrix.tsv",
        "pvalue_matrix": "dcifer_pvalue_matrix.tsv",
        "summary": "dcifer_summary.json",
    }

    def dcifer_outdir(outdir: Path) -> Path:
        return outdir / "dcifer"

    def dcifer_plots_dir(outdir: Path) -> Path:
        return dcifer_outdir(outdir) / "dcifer_plots"

    def dcifer_state_path(outdir: Path) -> Path:
        return dcifer_outdir(outdir) / "dcifer_state.json"

    def dcifer_plot_entries(outdir: Path) -> list[dict[str, Any]]:
        plots_dir = dcifer_plots_dir(outdir)
        if not safe_is_dir(plots_dir):
            return []
        entries: list[dict[str, Any]] = []
        try:
            paths = sorted(plots_dir.glob("*.png"), key=lambda item: item.name.lower())
        except OSError:
            return []
        for path in paths:
            if not safe_is_file(path):
                continue
            size_bytes = path.stat().st_size
            filename = path.name
            entries.append({
                "filename": filename,
                "title": Path(filename).stem.replace("_", " "),
                "exists": True,
                "size_bytes": size_bytes,
                "size": human_bytes(size_bytes),
                "view_url": url_for(
                    "download_dcifer_plot",
                    filename=filename,
                    out=str(outdir),
                    inline=1,
                ),
                "download_url": url_for(
                    "download_dcifer_plot",
                    filename=filename,
                    out=str(outdir),
                ),
            })
        return entries

    def dcifer_matrix_payload(path: Path, value_label: str, limit: int = 80) -> dict[str, Any]:
        if not safe_is_file(path):
            return {"labels": [], "rows": [], "value_label": value_label, "truncated": False}
        try:
            with path.open(newline="", encoding="utf-8-sig", errors="replace") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                fieldnames = reader.fieldnames or []
                row_label_field = str(fieldnames[0] or "sample_id") if fieldnames else "sample_id"
                all_labels = fieldnames[1:]
                labels = all_labels[:limit]
                rows: list[dict[str, Any]] = []
                row_count_seen = 0
                truncated = len(all_labels) > limit
                for index, row_item in enumerate(reader):
                    row_count_seen = index + 1
                    if index >= limit:
                        truncated = True
                        break
                    values: list[float | None] = []
                    for label in labels:
                        raw_value = (row_item.get(label) or "").strip()
                        if raw_value.upper() in {"", "NA", "NAN"}:
                            values.append(None)
                            continue
                        try:
                            values.append(float(raw_value))
                        except ValueError:
                            values.append(None)
                    row_label = dinemites_row_value(row_item, "sample_id", row_label_field)
                    rows.append({"sample_id": row_label, "values": values})
                return {
                    "labels": labels,
                    "rows": rows,
                    "value_label": value_label,
                    "total_columns": len(all_labels),
                    "total_rows": row_count_seen,
                    "truncated": truncated,
                }
        except Exception:
            return {"labels": [], "rows": [], "value_label": value_label, "truncated": False}

    def dcifer_process_active(outdir: Path) -> bool:
        key = str(outdir)
        with DCIFER_LOCK:
            thread = DCIFER_THREADS.get(key)
            if thread is None:
                return False
            if thread.is_alive():
                return True
            DCIFER_THREADS.pop(key, None)
        return False

    def run_dcifer_worker(
        root: Path,
        outdir: Path,
        samples: Path,
        settings: dict[str, int | float | str],
    ) -> None:
        dc_dir = dcifer_outdir(outdir)
        dc_dir.mkdir(parents=True, exist_ok=True)
        state_file = dcifer_state_path(outdir)
        cigar_path = outdir / "run_dada2" / "seqtab_cigar.tsv"
        dc_input = dc_dir / "dcifer_input_long.tsv"
        log_path = dc_dir / "dcifer.log"

        def state_payload(status: str, **extra: Any) -> dict[str, Any]:
            payload = {
                "status": status,
                **settings,
                "outdir": str(dc_dir),
            }
            payload.update(extra)
            return payload

        write_json(state_file, state_payload(
            "running",
            started_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        ))

        env = os.environ.copy()
        env["SIMPLSEQ_PROJECT_ROOT"] = str(root)
        kwargs: dict[str, Any] = {}
        if os.name == "nt":
            kwargs["creationflags"] = getattr(subprocess, "CREATE_NO_WINDOW", 0)

        try:
            cmd_convert = [
                "Rscript",
                str(root / "workflow" / "scripts" / "simplseq_to_dcifer.R"),
                "--cigar", str(cigar_path),
                "--samples", str(samples),
                "--out", str(dc_input),
                "--filter_summary", str(dc_dir / "dcifer_filter_summary.tsv"),
                "--replicate_summary", str(dc_dir / "dcifer_replicate_summary.tsv"),
                "--min_abundance_pct", str(settings["min_abundance_pct"]),
                "--abundance_denominator", str(settings["abundance_denominator"]),
            ]
            with log_path.open("w", encoding="utf-8") as log_handle:
                result = subprocess.run(
                    cmd_convert,
                    cwd=root,
                    env=env,
                    text=True,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    check=False,
                    timeout=600,
                    **kwargs,
                )
            if result.returncode != 0:
                write_json(state_file, state_payload(
                    "failed",
                    detail="simplseq_to_dcifer.R failed. Check dcifer.log.",
                    completed_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                ))
                return

            cmd_run = [
                "Rscript",
                str(root / "workflow" / "scripts" / "run_dcifer.R"),
                "--input", str(dc_input),
                "--outdir", str(dc_dir),
                "--coi_lrank", str(settings["coi_lrank"]),
                "--ibd_grid_nr", str(settings["ibd_grid_nr"]),
                "--alpha", str(settings["alpha"]),
                "--afreq_mode", str(settings["afreq_mode"]),
            ]
            with log_path.open("a", encoding="utf-8") as log_handle:
                result = subprocess.run(
                    cmd_run,
                    cwd=root,
                    env=env,
                    text=True,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    check=False,
                    timeout=7200,
                    **kwargs,
                )
            if result.returncode != 0:
                write_json(state_file, state_payload(
                    "failed",
                    detail="run_dcifer.R failed. Check dcifer.log.",
                    completed_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                ))
                return

            write_json(state_file, state_payload(
                "complete",
                completed_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            ))
        except Exception as exc:
            write_json(state_file, state_payload(
                "failed",
                detail=str(exc),
                completed_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            ))
        finally:
            with DCIFER_LOCK:
                DCIFER_THREADS.pop(str(outdir), None)

    @app.post("/api/dcifer/run")
    def api_dcifer_run():
        data = request.get_json(silent=True) or {}
        try:
            settings = dcifer_settings(data)
        except ValueError as exc:
            return json_error(str(exc), 400)

        outdir = resolve_app_path(workspace, data.get("outdir"), "results")
        samples = resolve_app_path(workspace, data.get("samples"), "samples.csv")
        cigar_path = outdir / "run_dada2" / "seqtab_cigar.tsv"
        if not cigar_path.exists():
            return json_error(
                f"Pipeline results not found. seqtab_cigar.tsv is required: {cigar_path}",
                400,
            )
        if not samples.exists():
            return json_error(f"Sample sheet not found: {samples}", 400)
        if dcifer_process_active(outdir):
            return json_error("A dcifer analysis is already running for this output folder.", 409)

        thread = threading.Thread(
            target=run_dcifer_worker,
            args=(app_root, outdir, samples, settings),
            daemon=True,
        )
        with DCIFER_LOCK:
            DCIFER_THREADS[str(outdir)] = thread
        thread.start()

        return jsonify({
            "ok": True,
            "status": "running",
            **settings,
            "outdir": str(outdir),
        })

    @app.get("/api/dcifer/status")
    def api_dcifer_status():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        state = read_json(dcifer_state_path(outdir))
        active = dcifer_process_active(outdir)
        status = state.get("status", "idle")
        if active and status != "running":
            status = "running"
        return jsonify({
            "ok": True,
            "status": status,
            "active": active,
            "state": state,
            "outdir": str(outdir),
        })

    @app.get("/api/dcifer/results")
    def api_dcifer_results():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        dc_dir = dcifer_outdir(outdir)
        state = read_json(dcifer_state_path(outdir))

        files: dict[str, dict[str, Any]] = {}
        for key, filename in DCIFER_FILE_KEYS.items():
            path = dc_dir / filename
            exists = safe_is_file(path)
            size_bytes = path.stat().st_size if exists else 0
            files[key] = {
                "exists": exists,
                "size_bytes": size_bytes,
                "size": human_bytes(size_bytes),
                "download_url": url_for("download_dcifer", file_key=key, out=str(outdir)) if exists else "",
            }

        plots = dcifer_plot_entries(outdir)
        has_outputs = any(item["exists"] for item in files.values()) or bool(plots)
        if state.get("status") in {None, "", "idle"} and has_outputs:
            state = {
                **state,
                "status": "complete",
                "detail": "Existing dcifer outputs found.",
                "outdir": str(dc_dir),
            }

        summary: dict[str, Any] = {
            "samples": "--",
            "pairs": "--",
            "max_relatedness": "--",
            "raw_p_le_alpha": "--",
        }
        summary_json = read_json(dc_dir / "dcifer_summary.json")
        for key in ("samples", "pairs", "max_relatedness", "raw_p_le_alpha"):
            if key in summary_json:
                summary[key] = summary_json[key]

        coi_rows = dinemites_tsv_rows(dc_dir / "dcifer_coi.tsv", limit=10000)
        if coi_rows:
            summary["samples"] = len({
                dinemites_row_value(row_item, "sample_id")
                for row_item in coi_rows
                if dinemites_row_value(row_item, "sample_id")
            })

        pair_rows = dinemites_tsv_rows(dc_dir / "dcifer_pairwise_relatedness.tsv", limit=10000)
        pairs: list[dict[str, str]] = []
        estimates: list[float] = []
        raw_p_count = 0
        alpha = float(state.get("alpha", summary_json.get("alpha", 0.05)) or 0.05)
        for index, row_item in enumerate(pair_rows):
            if index < 200:
                pairs.append({
                    "sample_a": dinemites_row_value(row_item, "sample_a"),
                    "sample_b": dinemites_row_value(row_item, "sample_b"),
                    "estimate": dinemites_row_value(row_item, "estimate"),
                    "p_value": dinemites_row_value(row_item, "p_value"),
                    "ci_lower": dinemites_row_value(row_item, "ci_lower", "CI_lower"),
                    "ci_upper": dinemites_row_value(row_item, "ci_upper", "CI_upper"),
                    "comparison_type": dinemites_row_value(row_item, "comparison_type"),
                })
            try:
                estimates.append(float(dinemites_row_value(row_item, "estimate")))
            except (TypeError, ValueError):
                pass
            raw_flag = dinemites_row_value(row_item, "raw_p_le_alpha").lower()
            if raw_flag in {"true", "t", "1", "yes"}:
                raw_p_count += 1
            elif not raw_flag:
                try:
                    if float(dinemites_row_value(row_item, "p_value")) <= alpha:
                        raw_p_count += 1
                except (TypeError, ValueError):
                    pass
        if pair_rows:
            summary["pairs"] = len(pair_rows)
            summary["raw_p_le_alpha"] = raw_p_count
        if estimates:
            summary["max_relatedness"] = round(max(estimates), 6)

        return jsonify({
            "ok": True,
            "state": state,
            "files": files,
            "plots": plots,
            "matrices": {
                "relatedness": dcifer_matrix_payload(dc_dir / "dcifer_relatedness_matrix.tsv", "Relatedness"),
                "pvalue": dcifer_matrix_payload(dc_dir / "dcifer_pvalue_matrix.tsv", "p-value"),
            },
            "summary": summary,
            "pairs": pairs,
            "outdir": str(outdir),
        })

    @app.get("/download/dcifer/<file_key>")
    def download_dcifer(file_key: str):
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        filename = DCIFER_FILE_KEYS.get(file_key)
        if not filename:
            abort(404)
        dc_dir = dcifer_outdir(outdir).resolve()
        path = (dc_dir / filename).resolve()
        try:
            path.relative_to(dc_dir)
        except ValueError:
            abort(404)
        if not path.exists() or not path.is_file():
            abort(404)
        return send_file(path, as_attachment=True, download_name=path.name)

    @app.get("/download/dcifer-plot/<path:filename>")
    def download_dcifer_plot(filename: str):
        if Path(filename).name != filename or "\\" in filename:
            abort(404)
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        plots_dir = dcifer_plots_dir(outdir).resolve()
        path = (plots_dir / filename).resolve()
        try:
            path.relative_to(plots_dir)
        except ValueError:
            abort(404)
        if path.suffix.lower() != ".png" or not path.exists() or not path.is_file():
            abort(404)
        inline = request.args.get("inline") == "1"
        return send_file(path, as_attachment=not inline, download_name=path.name, mimetype="image/png")

    return app


def open_browser_later(url: str) -> None:
    time.sleep(1.0)
    if is_wsl():
        try:
            subprocess.Popen(["cmd.exe", "/c", "start", "", url], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return
        except Exception:
            pass
    try:
        webbrowser.open(url)
    except Exception:
        pass


def run_server(root: Path | None = None, host: str = "127.0.0.1", port: int = 8501, open_browser: bool = True) -> int:
    if host not in {"127.0.0.1", "localhost", "::1"}:
        raise ValueError("SIMPLseq-nf App only serves the browser UI on a loopback host")
    app = create_app(root, Path.cwd())
    url = f"http://{host}:{port}"
    if open_browser:
        threading.Thread(target=open_browser_later, args=(url,), daemon=True).start()
    app.run(host=host, port=port, debug=False, use_reloader=False, threaded=True)
    return 0


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Run the SIMPLseq-nf App Flask GUI.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8501)
    parser.add_argument("--no-browser", action="store_true")
    args = parser.parse_args()
    return run_server(host=args.host, port=args.port, open_browser=not args.no_browser)


if __name__ == "__main__":
    raise SystemExit(main())
