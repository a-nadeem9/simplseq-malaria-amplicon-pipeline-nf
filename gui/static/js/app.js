const STAGES = [
  "prepare_inputs",
  "dada2",
  "prepare_stage2",
  "asv_mapping",
  "prepare_stage3",
  "cigar_check",
  "asv_to_cigar",
  "report"
];

const STAGE_LABELS = {
  prepare_inputs: "Checking sample sheet",
  dada2: "Running DADA2",
  prepare_stage2: "Cleaning ASV table",
  asv_mapping: "Mapping ASVs",
  prepare_stage3: "Preparing CIGAR inputs",
  cigar_check: "Checking CIGAR inputs",
  asv_to_cigar: "Converting ASVs to CIGAR",
  report: "Writing report"
};

const STORE_KEY = "simplseq.flask.settings";
const EMPTY_LOG_TEXT = "> Waiting for run.\n> Live progress appears here while SIMPLseq is active.";
let pollTimer = null;
let browseParent = "";
let scanInFlight = false;
let logInFlight = false;
let latestStatusPayload = null;
let lastRunStatus = "";
let completedRedirectKey = "";
let displayedProgressPercent = 0;
let progressAnimationFrame = null;
let activeOutdir = "";
let pathStyle = "";
let followLog = true;

function $(selector) {
  return document.querySelector(selector);
}

function text(node, value) {
  if (!node) return;
  node.textContent = value == null ? "" : String(value);
}

function terminalLineClass(line) {
  const trimmed = String(line || "").trim();
  if (trimmed.startsWith("[OK]")) return "log-ok";
  if (trimmed.startsWith("[WARN]") || trimmed.startsWith("WARN:")) return "log-warn";
  if (trimmed.startsWith("[ERROR]") || trimmed.startsWith("ERROR") || trimmed.includes("failed")) return "log-error";
  if (trimmed.startsWith(">")) return "log-muted";
  return "";
}

function renderTerminalLog(node, value) {
  if (!node) return;
  node.replaceChildren();
  const lines = String(value == null ? "" : value).split("\n");
  lines.forEach((line, index) => {
    const span = document.createElement("span");
    const lineClass = terminalLineClass(line);
    if (lineClass) span.className = lineClass;
    span.textContent = line;
    node.appendChild(span);
    if (index < lines.length - 1) node.appendChild(document.createTextNode("\n"));
  });
}

function scrollLogToBottom(node) {
  if (!node) return;
  requestAnimationFrame(() => {
    node.scrollTop = node.scrollHeight;
  });
}

function getSettings() {
  try {
    return JSON.parse(localStorage.getItem(STORE_KEY) || "{}");
  } catch (_error) {
    return {};
  }
}

function normalizeDrivePath(value) {
  const raw = String(value || "").trim();
  if (pathStyle !== "wsl") return raw;
  const match = raw.match(/^([A-Za-z]):[\\/]?(.*)$/);
  if (!match) return raw;
  const drive = match[1].toLowerCase();
  const rest = match[2].replaceAll("\\", "/");
  return `/mnt/${drive}/${rest}`;
}

function normalizePathInput(selector) {
  const node = $(selector);
  if (node && node.value) node.value = normalizeDrivePath(node.value);
}

function normalizePathInputs() {
  ["#fastq-dir", "#samples-out", "#run-samples", "#outdir", "#results-outdir", "#browse-path"].forEach(normalizePathInput);
  syncPathTitles();
}

function syncPathTitles() {
  ["#fastq-dir", "#samples-out", "#run-samples", "#outdir", "#results-outdir"].forEach((selector) => {
    const node = $(selector);
    if (node) node.title = node.value || "";
  });
}

function saveSettings() {
  normalizePathInputs();
  const settings = {
    fastqDir: $("#fastq-dir").value,
    samplesOut: $("#samples-out").value,
    runSamples: $("#run-samples").value,
    outdir: $("#outdir").value,
    runName: $("#run-name").value,
    resultsOutdir: $("#results-outdir").value,
    resumeRun: $("#resume-run").checked,
    dryRun: $("#dry-run").checked,
    cpus: $("#cpus").value,
    memory: $("#memory").value,
    dinemitesEnabled: $("#dinemites-enable") ? $("#dinemites-enable").checked : false,
    dinemitesModel: $("#dinemites-model") ? $("#dinemites-model").value : "simple",
    dinemitesNLags: $("#dinemites-n-lags") ? $("#dinemites-n-lags").value : "3",
    dinemitesTLag: $("#dinemites-t-lag") ? $("#dinemites-t-lag").value : "90",
    dinemitesMinAbundancePct: $("#dinemites-min-abundance-pct") ? $("#dinemites-min-abundance-pct").value : "0.3",
    dinemitesAbundanceDenominator: $("#dinemites-abundance-denominator") ? $("#dinemites-abundance-denominator").value : "locus",
    dinemitesNoDayCutoff: $("#dinemites-no-day-cutoff") ? $("#dinemites-no-day-cutoff").checked : true,
    dinemitesSeed: $("#dinemites-seed") ? $("#dinemites-seed").value : "1",
    dinemitesRefresh: $("#dinemites-refresh-interval") ? $("#dinemites-refresh-interval").value : "100",
    dinemitesBayesianLagDays: $("#dinemites-bayesian-lag-days") ? $("#dinemites-bayesian-lag-days").value : "30",
    dinemitesBayesianChains: $("#dinemites-bayesian-chains") ? $("#dinemites-bayesian-chains").value : "1",
    dinemitesBayesianParallelChains: $("#dinemites-bayesian-parallel-chains") ? $("#dinemites-bayesian-parallel-chains").value : "1",
    dinemitesBayesianWarmup: $("#dinemites-bayesian-warmup") ? $("#dinemites-bayesian-warmup").value : "500",
    dinemitesBayesianSampling: $("#dinemites-bayesian-sampling") ? $("#dinemites-bayesian-sampling").value : "500",
    dinemitesBayesianAdaptDelta: $("#dinemites-bayesian-adapt-delta") ? $("#dinemites-bayesian-adapt-delta").value : "0.99",
    dinemitesBayesianDropOut: $("#dinemites-bayesian-drop-out") ? $("#dinemites-bayesian-drop-out").checked : false,
    dciferEnabled: $("#dcifer-enable") ? $("#dcifer-enable").checked : false,
    dciferMinAbundancePct: $("#dcifer-min-abundance-pct") ? $("#dcifer-min-abundance-pct").value : "0.3",
    dciferAbundanceDenominator: $("#dcifer-abundance-denominator") ? $("#dcifer-abundance-denominator").value : "locus",
    dciferCoiLrank: $("#dcifer-coi-lrank") ? $("#dcifer-coi-lrank").value : "2",
    dciferIbdGridNr: $("#dcifer-ibd-grid-nr") ? $("#dcifer-ibd-grid-nr").value : "1000",
    dciferAlpha: $("#dcifer-alpha") ? $("#dcifer-alpha").value : "0.05"
  };
  localStorage.setItem(STORE_KEY, JSON.stringify(settings));
}

function restoreSettings() {
  const settings = getSettings();
  if (settings.fastqDir) $("#fastq-dir").value = settings.fastqDir;
  if (settings.samplesOut) $("#samples-out").value = settings.samplesOut;
  if (settings.runSamples) $("#run-samples").value = settings.runSamples;
  if (settings.outdir) $("#outdir").value = settings.outdir;
  if (settings.runName) $("#run-name").value = settings.runName;
  if (settings.resultsOutdir) $("#results-outdir").value = settings.resultsOutdir;
  if (typeof settings.resumeRun === "boolean") $("#resume-run").checked = settings.resumeRun;
  if (typeof settings.dryRun === "boolean") $("#dry-run").checked = settings.dryRun;
  if (settings.cpus != null) $("#cpus").value = settings.cpus;
  if (settings.memory != null) $("#memory").value = settings.memory;
  if (typeof settings.dinemitesEnabled === "boolean" && $("#dinemites-enable")) $("#dinemites-enable").checked = settings.dinemitesEnabled;
  if (settings.dinemitesModel && $("#dinemites-model")) $("#dinemites-model").value = settings.dinemitesModel;
  if (settings.dinemitesNLags != null && $("#dinemites-n-lags")) $("#dinemites-n-lags").value = settings.dinemitesNLags;
  if (settings.dinemitesTLag != null && $("#dinemites-t-lag")) $("#dinemites-t-lag").value = settings.dinemitesTLag;
  if (settings.dinemitesMinAbundancePct != null && $("#dinemites-min-abundance-pct")) {
    $("#dinemites-min-abundance-pct").value = settings.dinemitesMinAbundancePct;
  }
  if (settings.dinemitesAbundanceDenominator && $("#dinemites-abundance-denominator")) {
    $("#dinemites-abundance-denominator").value = settings.dinemitesAbundanceDenominator;
  }
  if (typeof settings.dinemitesNoDayCutoff === "boolean" && $("#dinemites-no-day-cutoff")) {
    $("#dinemites-no-day-cutoff").checked = settings.dinemitesNoDayCutoff;
  }
  if (settings.dinemitesSeed != null && $("#dinemites-seed")) $("#dinemites-seed").value = settings.dinemitesSeed;
  if (settings.dinemitesRefresh != null && $("#dinemites-refresh-interval")) $("#dinemites-refresh-interval").value = settings.dinemitesRefresh;
  if (settings.dinemitesBayesianLagDays != null && $("#dinemites-bayesian-lag-days")) $("#dinemites-bayesian-lag-days").value = settings.dinemitesBayesianLagDays;
  if (settings.dinemitesBayesianChains != null && $("#dinemites-bayesian-chains")) $("#dinemites-bayesian-chains").value = settings.dinemitesBayesianChains;
  if (settings.dinemitesBayesianParallelChains != null && $("#dinemites-bayesian-parallel-chains")) {
    $("#dinemites-bayesian-parallel-chains").value = settings.dinemitesBayesianParallelChains;
  }
  if (settings.dinemitesBayesianWarmup != null && $("#dinemites-bayesian-warmup")) $("#dinemites-bayesian-warmup").value = settings.dinemitesBayesianWarmup;
  if (settings.dinemitesBayesianSampling != null && $("#dinemites-bayesian-sampling")) $("#dinemites-bayesian-sampling").value = settings.dinemitesBayesianSampling;
  if (settings.dinemitesBayesianAdaptDelta != null && $("#dinemites-bayesian-adapt-delta")) {
    $("#dinemites-bayesian-adapt-delta").value = settings.dinemitesBayesianAdaptDelta;
  }
  if (typeof settings.dinemitesBayesianDropOut === "boolean" && $("#dinemites-bayesian-drop-out")) {
    $("#dinemites-bayesian-drop-out").checked = settings.dinemitesBayesianDropOut;
  }
  if (typeof settings.dciferEnabled === "boolean" && $("#dcifer-enable")) $("#dcifer-enable").checked = settings.dciferEnabled;
  if (settings.dciferMinAbundancePct != null && $("#dcifer-min-abundance-pct")) {
    $("#dcifer-min-abundance-pct").value = settings.dciferMinAbundancePct;
  }
  if (settings.dciferAbundanceDenominator && $("#dcifer-abundance-denominator")) {
    $("#dcifer-abundance-denominator").value = settings.dciferAbundanceDenominator;
  }
  if (settings.dciferCoiLrank != null && $("#dcifer-coi-lrank")) $("#dcifer-coi-lrank").value = settings.dciferCoiLrank;
  if (settings.dciferIbdGridNr != null && $("#dcifer-ibd-grid-nr")) $("#dcifer-ibd-grid-nr").value = settings.dciferIbdGridNr;
  if (settings.dciferAlpha != null && $("#dcifer-alpha")) $("#dcifer-alpha").value = settings.dciferAlpha;
  normalizePathInputs();
}

async function fetchJson(url, options = {}, settings = {}) {
  const {allowAppError = false} = settings;
  const response = await fetch(url, {
    headers: {"Content-Type": "application/json"},
    ...options
  });
  const payload = await response.json().catch(() => ({}));
  if (!response.ok || (!allowAppError && payload.ok === false)) {
    throw new Error(payload.error || `Request failed: ${response.status}`);
  }
  return payload;
}

function postJson(url, body, settings = {}) {
  return fetchJson(url, {
    method: "POST",
    body: JSON.stringify(body)
  }, settings);
}

function delay(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

async function withButtonFeedback(button, busyText, action) {
  if (!button) return action();
  const oldText = button.textContent;
  button.disabled = true;
  text(button, busyText);
  try {
    await Promise.all([action(), delay(900)]);
  } finally {
    text(button, oldText);
    button.disabled = false;
  }
}

function bindPlotWheelScrolling() {
  document.addEventListener("wheel", (event) => {
    const target = event.target instanceof Element ? event.target : null;
    const plotCard = target?.closest(".dinemites-plot-card, .dcifer-heatmap-card");
    if (!plotCard || event.defaultPrevented || event.shiftKey) return;
    if (Math.abs(event.deltaY) <= Math.abs(event.deltaX)) return;

    const scrollRoot = document.scrollingElement || document.documentElement;
    const maxScrollTop = scrollRoot.scrollHeight - window.innerHeight;
    const scrollingDown = event.deltaY > 0;
    const canScrollPage = scrollingDown
      ? window.scrollY < maxScrollTop - 1
      : window.scrollY > 1;
    if (!canScrollPage) return;

    const deltaScale = event.deltaMode === 1
      ? 24
      : event.deltaMode === 2
        ? window.innerHeight
        : 1;
    event.preventDefault();
    window.scrollBy({top: event.deltaY * deltaScale, left: 0, behavior: "auto"});
  }, {passive: false});
}

function setPill(node, label, status) {
  if (!node) return;
  node.classList.remove("ok", "warn", "bad");
  if (status) node.classList.add(status);
  text(node, label);
}

function payloadStatus(payload) {
  const state = payload?.state || {};
  const summary = payload?.summary || {};
  const status = state.status || summary.status || "pending";
  if (payload?.active && status === "pending") return "starting";
  return status;
}

function isActiveStatus(status) {
  return status === "starting" || status === "running";
}

function checkStatusReady() {
  return ($("#check-status")?.textContent || "").trim().toLowerCase() === "ready";
}

function stageClass(status) {
  if (status === "complete") return "complete";
  if (status === "started" || status === "running") return "running";
  if (status === "failed" || status === "error") return "failed";
  return "pending";
}

function stageStatusLabel(status) {
  const klass = stageClass(status);
  if (klass === "running") return "running";
  if (klass === "complete") return "complete";
  if (klass === "failed") return "failed";
  return "waiting";
}

function setProgressDisplay(value) {
  const percent = Math.max(0, Math.min(100, Math.round(Number(value) || 0)));
  const progressFill = $("#progress-fill");
  const globalFill = $("#global-progress-fill");
  if (progressFill) progressFill.style.width = `${percent}%`;
  if (globalFill) globalFill.style.width = `${percent}%`;
  text($("#progress-percent"), `${percent}%`);
  text($("#pipeline-percent"), `${percent}%`);
}

function animateProgressTo(percent) {
  const target = Math.max(0, Math.min(100, Math.round(Number(percent) || 0)));
  const start = displayedProgressPercent;
  if (progressAnimationFrame) {
    cancelAnimationFrame(progressAnimationFrame);
    progressAnimationFrame = null;
  }
  if (Math.abs(target - start) < 1) {
    displayedProgressPercent = target;
    setProgressDisplay(target);
    return;
  }
  const duration = Math.min(1400, Math.max(650, Math.abs(target - start) * 42));
  const startedAt = performance.now();
  const tick = (now) => {
    const elapsed = Math.min(1, (now - startedAt) / duration);
    const eased = 1 - Math.pow(1 - elapsed, 3);
    setProgressDisplay(start + (target - start) * eased);
    if (elapsed < 1) {
      progressAnimationFrame = requestAnimationFrame(tick);
      return;
    }
    displayedProgressPercent = target;
    progressAnimationFrame = null;
    setProgressDisplay(target);
  };
  progressAnimationFrame = requestAnimationFrame(tick);
}

function compactTerminalLog(raw, statusPayload) {
  const value = String(raw || "").replace(/\r\n/g, "\n").trimEnd();
  const status = payloadStatus(statusPayload);
  if (!value) {
    return isActiveStatus(status)
      ? "> Starting SIMPLseq.\n> Live Nextflow progress will appear here."
      : EMPTY_LOG_TEXT;
  }

  const lines = value.split("\n");
  const nextflowStart = lines.findIndex((line) => {
    const trimmed = line.trim();
    return (
      trimmed.includes("N E X T F L O W") ||
      trimmed.startsWith("WARN:") ||
      trimmed.startsWith("Launching `") ||
      trimmed.startsWith("executor >") ||
      /^\[[0-9a-f]{2}\//i.test(trimmed)
    );
  });
  if (nextflowStart >= 0) {
    return lines.slice(nextflowStart).slice(-220).join("\n").trimEnd();
  }

  const filtered = lines.filter((line) => {
    const trimmed = line.trim();
    if (!trimmed) return false;
    if (trimmed.startsWith("[SIMPLseq/App]")) return false;
    if (trimmed.startsWith("[OK]")) return false;
    if (trimmed.startsWith("[WARN] Large/high-depth datasets")) return false;
    if (trimmed.startsWith("[SIMPLseq] preparing")) return false;
    if (trimmed.startsWith("[SIMPLseq] wrote runtime versions")) return false;
    if (trimmed.startsWith("[SIMPLseq] wrote input FASTQ MD5 table")) return false;
    return true;
  });
  const cleaned = filtered.length ? filtered : lines.slice(-40);
  return cleaned.slice(-180).join("\n").trimEnd() || EMPTY_LOG_TEXT;
}

function setStep(name, active) {
  const node = document.querySelector(`.step[data-step="${name}"]`);
  if (node) node.classList.toggle("is-active", active);
}

function selectTab(name) {
  document.querySelectorAll(".tab").forEach((tab) => {
    tab.classList.toggle("is-active", tab.dataset.tab === name);
  });
  document.querySelectorAll(".tab-panel").forEach((panel) => {
    panel.classList.toggle("is-active", panel.id === `tab-${name}`);
  });
  if (name === "dinemites") {
    loadDinemitesResults();
  } else if (name === "dcifer") {
    loadDciferResults();
  }
}

function setFolderMessage(message, className = "") {
  const node = $("#folder-message");
  if (!node) return;
  node.className = `field-help ${className}`.trim();
  text(node, message);
}

function row(cells) {
  const tr = document.createElement("tr");
  cells.forEach((value) => {
    const td = document.createElement("td");
    text(td, value);
    tr.appendChild(td);
  });
  return tr;
}

function emptyRow(colspan, message) {
  const tr = document.createElement("tr");
  const td = document.createElement("td");
  td.colSpan = colspan;
  td.className = "empty";
  text(td, message);
  tr.appendChild(td);
  return tr;
}

function displayMissing(value) {
  return value === null || value === undefined || value === "" ? "--" : value;
}

function formatNumber(value, digits = 2) {
  if (value === null || value === undefined || value === "") return "--";
  const number = Number(value);
  if (!Number.isFinite(number)) return String(value);
  if (Number.isInteger(number)) return String(number);
  return number.toFixed(digits).replace(/\.?0+$/, "");
}

function formatPValue(value) {
  if (value === null || value === undefined || value === "") return "--";
  const number = Number(value);
  if (!Number.isFinite(number)) return String(value);
  if (number > 0 && number < 0.001) return number.toExponential(2);
  return formatNumber(number, 3);
}

function updatePlotJump(gallery, countNode, button, count, emptyText, singularText, pluralText) {
  if (countNode) {
    text(countNode, count ? `${count} ${count === 1 ? singularText : pluralText} available` : emptyText);
  }
  if (!button) return;
  button.disabled = !count;
  button.onclick = count && gallery
    ? () => gallery.scrollIntoView({ behavior: "smooth", block: "start" })
    : null;
}

function renderScanPreview(items) {
  const tbody = $("#scan-preview");
  tbody.replaceChildren();
  if (!items.length) {
    tbody.appendChild(emptyRow(4, "No paired FASTQs found."));
    return;
  }
  items.forEach((item) => {
    tbody.appendChild(row([
      item.sample_id,
      item.participant_id,
      item.collection_date,
      item.replicate
    ]));
  });
}

function renderSamplePreview(items) {
  const tbody = $("#sample-preview");
  tbody.replaceChildren();
  if (!items.length) {
    tbody.appendChild(emptyRow(5, "No sample rows available."));
    return;
  }
  items.forEach((item) => {
    tbody.appendChild(row([
      item.sample_id,
      item.fastq_1,
      item.fastq_2,
      item.participant_id,
      item.collection_date
    ]));
  });
}

function renderWarnings(scan) {
  const box = $("#scan-warnings");
  box.replaceChildren();
  const notices = [];
  const pairCount = Number(scan.pair_count || 0);
  const totalBytes = Number(scan.total_fastq_bytes || 0);
  if (pairCount >= 50 || totalBytes >= 5 * 1024 * 1024 * 1024) {
    notices.push({
      className: "notice",
      text: "Large dataset detected. Full runs can require much more memory than a typical laptop. Run a small test first if this is a new setup."
    });
  }
  if (scan.duplicate_sample_ids && scan.duplicate_sample_ids.length) {
    notices.push({
      className: "notice bad",
      text: `Duplicate sample IDs: ${scan.duplicate_sample_ids.slice(0, 8).join(", ")}`
    });
  }
  if (scan.missing_r2 && scan.missing_r2.length) {
    notices.push({className: "notice", text: `Missing R2 files for ${scan.missing_r2.length} R1 files.`});
  }
  if (scan.orphan_r2 && scan.orphan_r2.length) {
    notices.push({className: "notice", text: `Found ${scan.orphan_r2.length} R2 files without R1.`});
  }
  notices.forEach((notice) => {
    const div = document.createElement("div");
    div.className = notice.className;
    text(div, notice.text);
    box.appendChild(div);
  });
}

async function scanFastqs() {
  if (scanInFlight) return;
  scanInFlight = true;
  saveSettings();
  setPill($("#scan-status"), "Scanning", "warn");
  $("#scan-button").disabled = true;
  try {
    const payload = await postJson("/api/scan", {
      fastq_dir: $("#fastq-dir").value,
      samples_out: $("#samples-out").value,
      include_pool_in_sample_id: false,
      absolute_paths: true,
      write_samples: true
    });
    text($("#metric-pairs"), payload.pair_count);
    text($("#metric-size"), payload.total_fastq_size);
    text($("#metric-md5"), payload.md5_files);
    text($("#metric-missing"), payload.missing_pairs);
    text($("#sample-sheet-title"), payload.samples_relative || payload.samples_out);
    renderScanPreview(payload.preview || []);
    renderSamplePreview(payload.sample_preview || []);
    renderWarnings(payload);
    $("#run-samples").value = $("#samples-out").value;
    if (payload.samples_written) {
      setPill($("#scan-status"), `${payload.pair_count} pairs`, payload.pair_count ? "ok" : "warn");
      setPill($("#sample-sheet-status"), `Wrote ${payload.sample_rows_written} rows`, payload.pair_count ? "ok" : "warn");
      setStep("review", true);
    } else if (payload.duplicate_sample_ids && payload.duplicate_sample_ids.length) {
      setPill($("#scan-status"), "Duplicates", "bad");
      setPill($("#sample-sheet-status"), "Not written", "bad");
    } else {
      setPill($("#scan-status"), "No pairs", "warn");
      setPill($("#sample-sheet-status"), "Empty", "warn");
    }
    saveSettings();
  } catch (error) {
    setPill($("#scan-status"), "Scan failed", "bad");
    renderWarnings({duplicate_sample_ids: [error.message]});
  } finally {
    scanInFlight = false;
    $("#scan-button").disabled = false;
  }
}

function renderChecks(checks) {
  const tbody = $("#check-table");
  tbody.replaceChildren();
  if (!checks.length) {
    tbody.appendChild(emptyRow(3, "No checks returned."));
    return;
  }
  checks.forEach((item) => {
    const tr = document.createElement("tr");
    const name = document.createElement("td");
    const status = document.createElement("td");
    const detail = document.createElement("td");
    const badge = document.createElement("span");
    badge.className = `check-badge ${item.status === "ok" ? "ok" : item.status === "warn" ? "warn" : "bad"}`;
    text(name, item.name);
    text(badge, item.status);
    text(detail, item.detail);
    status.appendChild(badge);
    tr.append(name, status, detail);
    tbody.appendChild(tr);
  });
}

async function runCheck() {
  saveSettings();
  setPill($("#check-status"), "Checking", "warn");
  $("#check-button").disabled = true;
  try {
    const payload = await postJson("/api/check", {samples: $("#run-samples").value, outdir: $("#outdir").value}, {allowAppError: true});
    renderChecks(payload.checks || []);
    await loadLog({silent: true});
    if (payload.failed) {
      setPill($("#check-status"), `${payload.failed} need attention`, "bad");
    } else {
      setPill($("#check-status"), "Ready", "ok");
      setStep("runtime", true);
      resetRunDisplay();
    }
  } catch (error) {
    setPill($("#check-status"), "Check failed", "bad");
    renderChecks([{name: "Runtime check", status: "missing", detail: error.message}]);
  } finally {
    $("#check-button").disabled = false;
  }
}

function runPayload() {
  return {
    samples: $("#run-samples").value,
    outdir: $("#outdir").value,
    run_name: $("#run-name").value,
    resume: $("#resume-run").checked,
    dry_run: $("#dry-run").checked,
    cpus: Number($("#cpus").value || 0),
    memory: $("#memory").value
  };
}

async function startRun() {
  saveSettings();
  $("#run-button").disabled = true;
  $("#run-message").className = "inline-message";
  text($("#run-message"), $("#dry-run").checked ? "Starting preview..." : "Starting run...");
  try {
    const payload = await postJson("/api/run", runPayload());
    activeOutdir = payload.outdir;
    $("#results-outdir").value = payload.outdir;
    $("#run-name").value = "";
    saveSettings();
    text($("#run-message"), payload.dry_run ? `Preview ready in ${payload.outdir}.` : `Run started in ${payload.outdir}.`);
    $("#run-message").classList.add("ok");
    if (payload.dry_run) {
      setPill($("#run-state-pill"), "Preview", "warn");
      renderStages([], {status: "dry_run"}, {status: "dry_run"});
      await loadLog({forceScroll: true, statusPayload: {active: false, state: {status: "dry_run"}}});
    } else {
      showLaunchState(payload.outdir);
    }
    setStep("collect", true);
    startPolling();
    setTimeout(refreshAllRunState, 1200);
  } catch (error) {
    text($("#run-message"), error.message);
    $("#run-message").classList.add("bad");
  } finally {
    $("#run-button").disabled = false;
  }
}

async function stopRun() {
  if (!activeOutdir) return;
  const confirmed = window.confirm("Stop the active SIMPLseq run? The output folder will stay on disk so you can resume later.");
  if (!confirmed) return;
  const button = $("#stop-button");
  const oldLabel = button ? button.textContent : "";
  if (button) {
    button.disabled = true;
    text(button, "Stopping...");
  }
  $("#run-message").className = "inline-message";
  text($("#run-message"), "Stopping run...");
  try {
    const payload = await postJson("/api/stop-run", {outdir: currentRunOutdir()});
    latestStatusPayload = payload;
    lastRunStatus = "stopped";
    renderStatus(payload);
    renderStages(payload.events || [], payload.summary || {status: "stopped"}, payload.state || {status: "stopped"});
    await loadLog({forceScroll: true, statusPayload: payload});
    text($("#run-message"), `Run stopped. To continue this output folder later, keep the same run name and enable resume.`);
    $("#run-message").classList.add("ok");
    await refreshAllRunState();
  } catch (error) {
    text($("#run-message"), error.message);
    $("#run-message").classList.add("bad");
  } finally {
    if (button) {
      text(button, oldLabel || "Stop run");
    }
  }
}

function showLaunchState(outdir) {
  followLog = true;
  lastRunStatus = "starting";
  displayedProgressPercent = 0;
  setProgressDisplay(0);
  setPill($("#run-state-pill"), "Running", "ok");
  renderStages([], {status: "starting", current_stage: "pending"}, {status: "starting"});
  renderTerminalLog(
    $("#technical-log"),
    `> Run started.\n> Output folder: ${outdir}\n> Waiting for first Nextflow update.`
  );
  scrollLogToBottom($("#technical-log"));
}

function openFolderModal() {
  const modal = $("#folder-modal");
  if (!modal) return;
  modal.classList.add("is-open");
  modal.setAttribute("aria-hidden", "false");
}

function closeFolderModal() {
  const modal = $("#folder-modal");
  if (!modal) return;
  modal.classList.remove("is-open");
  modal.setAttribute("aria-hidden", "true");
}

async function openFallbackFolderBrowser() {
  $("#browse-path").value = $("#fastq-dir").value || $("#browse-path").value || ".";
  openFolderModal();
  await loadBrowse($("#browse-path").value);
}

async function chooseFastqFolder() {
  const button = $("#browse-button");
  const scanButton = $("#scan-button");
  const oldLabel = button ? button.textContent : "";
  if (button) {
    button.disabled = true;
    text(button, "Opening...");
  }
  if (scanButton && !scanInFlight) scanButton.disabled = true;
  setFolderMessage("Opening folder picker...");
  const pickerHint = setTimeout(() => {
    if (button) text(button, "Waiting...");
    setFolderMessage("Native folder picker is open. If you do not see it, check behind this browser or on the taskbar/Dock.", "ok");
  }, 1500);
  try {
    const payload = await postJson("/api/select-folder", {initial: $("#fastq-dir").value});
    if (payload.selected && payload.path) {
      $("#fastq-dir").value = payload.path;
      $("#browse-path").value = payload.path;
      saveSettings();
      setFolderMessage("Folder selected. Click Scan folder when ready.", "ok");
      return;
    }
    setFolderMessage("Folder selection was cancelled.");
  } catch (error) {
    setFolderMessage("Native folder picker was not available. Opening the manual browser instead.", "warn");
    await openFallbackFolderBrowser();
  } finally {
    clearTimeout(pickerHint);
    if (button) {
      button.disabled = false;
      text(button, oldLabel || "Choose folder");
    }
    if (scanButton && !scanInFlight) scanButton.disabled = false;
  }
}

async function chooseOutputFolder() {
  const button = $("#choose-outdir-button");
  const oldLabel = button ? button.textContent : "";
  if (button) {
    button.disabled = true;
    text(button, "Opening...");
  }
  text($("#run-message"), "Opening output folder picker...");
  $("#run-message").className = "inline-message";
  try {
    const payload = await postJson("/api/select-folder", {initial: $("#outdir").value});
    if (payload.selected && payload.path) {
      $("#outdir").value = payload.path;
      activeOutdir = "";
      saveSettings();
      text($("#run-message"), "Output parent folder selected.");
      $("#run-message").classList.add("ok");
      return;
    }
    text($("#run-message"), "Output folder selection was cancelled.");
  } catch (error) {
    text($("#run-message"), "Folder picker was not available. Type the output folder path manually.");
    $("#run-message").classList.add("bad");
  } finally {
    if (button) {
      button.disabled = false;
      text(button, oldLabel || "Choose folder");
    }
  }
}

function renderStages(events, summary, state) {
  const statusByStage = {};
  const messageByStage = {};
  const staleState = state?.status === "stale" || state?.status === "stopped";
  const runStatus = state?.status || summary.status || "pending";
  const hasStageEvents = events.some((event) => STAGES.includes(event.stage));
  STAGES.forEach((stage) => {
    statusByStage[stage] = "pending";
  });
  events.forEach((event) => {
    if (!STAGES.includes(event.stage)) return;
    if (staleState && ["started", "running"].includes(event.status)) return;
    statusByStage[event.stage] = event.status;
    if (event.message) messageByStage[event.stage] = event.message;
  });
  if (state && state.status === "failed") {
    const current = summary.current_stage;
    if (STAGES.includes(current) && statusByStage[current] !== "complete") {
      statusByStage[current] = "failed";
    }
  }
  if (isActiveStatus(runStatus) && !hasStageEvents && STAGES.length) {
    statusByStage[STAGES[0]] = "running";
    messageByStage[STAGES[0]] = "Launching workflow";
  }

  const list = $("#stage-list");
  list.replaceChildren();
  STAGES.forEach((stage) => {
    const li = document.createElement("li");
    const status = statusByStage[stage];
    li.className = stageClass(status);
    const dot = document.createElement("span");
    dot.className = "stage-dot";
    const label = document.createElement("span");
    label.className = "stage-label";
    label.title = messageByStage[stage] || "";
    text(label, STAGE_LABELS[stage] || stage);
    const stateNode = document.createElement("span");
    stateNode.className = "stage-status";
    text(stateNode, stageStatusLabel(status));
    li.append(dot, label, stateNode);
    list.appendChild(li);
  });

  const completed = summary.completed_stages || 0;
  const total = summary.total_stages || STAGES.length;
  const percent = runStatus === "dry_run" || runStatus === "complete" ? 100 : Math.round((completed / Math.max(total, 1)) * 100);
  const indeterminate = isActiveStatus(runStatus) && percent === 0;
  const progressFill = $("#progress-fill");
  const progressBar = progressFill?.parentElement;
  const globalFill = $("#global-progress-fill");
  const globalBar = globalFill?.parentElement;
  [progressBar, globalBar].forEach((bar) => {
    if (!bar) return;
    bar.classList.remove("is-running", "is-complete", "is-failed", "is-stopped");
    if (isActiveStatus(runStatus)) bar.classList.add("is-running");
    if (runStatus === "complete" || runStatus === "dry_run") bar.classList.add("is-complete");
    if (runStatus === "failed") bar.classList.add("is-failed");
    if (runStatus === "stopped") bar.classList.add("is-stopped");
    bar.classList.toggle("is-indeterminate", indeterminate);
  });
  animateProgressTo(percent);
  let pipelineLabel = checkStatusReady() ? "Ready" : "Idle";
  let pipelineClass = checkStatusReady() ? "is-ready" : "is-idle";
  let showPipelinePercent = false;
  if (isActiveStatus(runStatus)) {
    pipelineLabel = "Running";
    pipelineClass = "is-running";
    showPipelinePercent = true;
  }
  if (runStatus === "complete") pipelineLabel = "Complete";
  if (runStatus === "complete") {
    pipelineClass = "is-complete";
    showPipelinePercent = true;
  }
  if (runStatus === "dry_run") {
    pipelineLabel = "Preview ready";
    pipelineClass = "is-preview";
  }
  if (runStatus === "failed") {
    pipelineLabel = "Failed";
    pipelineClass = "is-failed";
    showPipelinePercent = true;
  }
  if (runStatus === "stopped") {
    pipelineLabel = "Stopped";
    pipelineClass = "is-unknown";
    showPipelinePercent = true;
  }
  if (!["pending", "starting", "running", "complete", "dry_run", "failed", "stopped"].includes(runStatus)) {
    pipelineLabel = "Check status";
    pipelineClass = "is-unknown";
  }
  text($("#pipeline-label"), pipelineLabel);

  const pipelineStatus = $(".pipeline-status");
  if (pipelineStatus) {
    pipelineStatus.classList.remove("is-idle", "is-ready", "is-running", "is-complete", "is-preview", "is-failed", "is-unknown", "has-percent");
    pipelineStatus.classList.add(pipelineClass);
    if (showPipelinePercent) pipelineStatus.classList.add("has-percent");
  }
  const current = STAGES.includes(summary.current_stage) ? summary.current_stage : "";
  let currentLabel = current ? (STAGE_LABELS[current] || current) : "No active run";
  if (isActiveStatus(runStatus) && !current) currentLabel = "Launching workflow";
  if (runStatus === "dry_run") currentLabel = "Preview ready";
  if (runStatus === "complete") currentLabel = "Run complete";
  if (runStatus === "failed") currentLabel = "Run failed";
  if (runStatus === "stopped") currentLabel = "Run stopped";
  if (runStatus === "stale") currentLabel = "Check run status";
  text($("#current-stage"), currentLabel);
}

function renderStatus(payload) {
  const status = payloadStatus(payload);
  let pillStatus = "";
  if (status === "complete" || status === "dry_run") pillStatus = "ok";
  if (status === "failed") pillStatus = "bad";
  if (status === "stopped") pillStatus = "warn";
  if (status === "running" || status === "starting") pillStatus = "ok";
  $("#stop-button").hidden = !Boolean(payload.active);
  $("#stop-button").disabled = !Boolean(payload.active);
  setPill($("#run-state-pill"), status, pillStatus);
  $("#run-button").disabled = Boolean(payload.active && !$("#dry-run").checked);
}

function resetRunDisplay() {
  latestStatusPayload = null;
  lastRunStatus = "";
  completedRedirectKey = "";
  if (progressAnimationFrame) {
    cancelAnimationFrame(progressAnimationFrame);
    progressAnimationFrame = null;
  }
  displayedProgressPercent = 0;
  setProgressDisplay(0);
  renderStatus({active: false, state: {status: "pending"}, summary: {status: "pending"}});
  renderStages([], {status: "pending", completed_stages: 0, total_stages: STAGES.length}, {status: "pending"});
  renderTerminalLog($("#technical-log"), EMPTY_LOG_TEXT);
}

async function fetchActiveRun() {
  return fetchJson("/api/active-run");
}

async function refreshStatus() {
  const out = encodeURIComponent(currentRunOutdir());
  const payload = await fetchJson(`/api/status?out=${out}`);
  renderStatus(payload);
  return payload;
}

async function refreshProgress(statusPayload = null) {
  const out = encodeURIComponent(currentRunOutdir());
  const payload = await fetchJson(`/api/progress?out=${out}`);
  const state = statusPayload ? statusPayload.state : {};
  renderStages(payload.events || [], payload.summary || {}, state || {});
  return payload;
}

async function refreshResults() {
  const out = $("#results-outdir").value || $("#outdir").value || "results";
  const payload = await fetchJson(`/api/results?out=${encodeURIComponent(out)}`);
  updateBundleButton(payload);
  const list = $("#results-list");
  list.replaceChildren();
  if (!payload.files || !payload.files.length) {
    const div = document.createElement("div");
    div.className = "notice";
    text(div, "No result manifest entries found.");
    list.appendChild(div);
    return;
  }
  renderResultsDashboard(payload, list);
}

function updateBundleButton(payload) {
  const button = $("#download-bundle");
  if (!button) return;
  button.disabled = !payload.bundle_ready || !payload.bundle_url;
  button.dataset.url = payload.bundle_url || "";
  button.title = button.disabled ? "Output bundle is available after result files are ready." : "Download report and core result tables as a zip.";
}

function resultAction(url, label, primary = false) {
  const link = document.createElement("a");
  link.href = url;
  if (label.startsWith("Open")) link.target = "_blank";
  const button = document.createElement("button");
  if (primary) button.className = "primary";
  text(button, label);
  link.appendChild(button);
  return link;
}

function resultPill(item) {
  const pill = document.createElement("span");
  const ready = item && item.status === "ready";
  pill.className = `status-pill ${ready ? "ok" : "warn"}`;
  text(pill, ready ? "ready" : "pending");
  return pill;
}

function disabledResultAction(label) {
  const button = document.createElement("button");
  button.disabled = true;
  text(button, label);
  return button;
}

function resultSummaryCard(label, value, detail, item = null) {
  const card = document.createElement("div");
  card.className = "results-summary-card";
  const badge = document.createElement("div");
  badge.className = "results-summary-badge";
  text(badge, value);
  const meta = document.createElement("div");
  const title = document.createElement("strong");
  text(title, label);
  const description = document.createElement("span");
  text(description, detail);
  meta.append(title, description);
  card.append(badge, meta);
  if (item) card.appendChild(resultPill(item));
  return card;
}

function resultTableRow(item) {
  const row = document.createElement("tr");
  if (item.status !== "ready") row.className = "is-pending";
  const fileCell = document.createElement("td");
  fileCell.className = "result-file-cell";
  const title = document.createElement("strong");
  text(title, item.label);
  const detail = document.createElement("span");
  text(detail, item.relative_path);
  fileCell.append(title, detail);

  const sizeCell = document.createElement("td");
  text(sizeCell, item.status === "ready" ? item.size : "--");

  const actionsCell = document.createElement("td");
  const actions = document.createElement("div");
  actions.className = "result-actions";
  if (item.view_url) actions.appendChild(resultAction(item.view_url, "Open"));
  if (item.download_url) actions.appendChild(resultAction(item.download_url, "Download"));
  if (!item.view_url && !item.download_url) actions.appendChild(resultPill(item));
  actionsCell.appendChild(actions);

  row.append(fileCell, sizeCell, actionsCell);
  return row;
}

function resultTable(files, emptyText) {
  const section = document.createElement("section");
  section.className = "results-table-section";
  const wrap = document.createElement("div");
  wrap.className = "results-table-wrap";
  const table = document.createElement("table");
  table.className = "results-table";
  const thead = document.createElement("thead");
  thead.innerHTML = "<tr><th>Table</th><th>Size</th><th></th></tr>";
  const tbody = document.createElement("tbody");
  if (!files.length) {
    const row = document.createElement("tr");
    const cell = document.createElement("td");
    cell.colSpan = 3;
    cell.className = "empty";
    text(cell, emptyText);
    row.appendChild(cell);
    tbody.appendChild(row);
  } else {
    files.forEach((item) => tbody.appendChild(resultTableRow(item)));
  }
  table.append(thead, tbody);
  wrap.appendChild(table);
  section.appendChild(wrap);
  return section;
}

function renderResultsDashboard(payload, list) {
  const report = payload.report;
  const coreFiles = payload.core_files || [];

  const dashboard = document.createElement("div");
  dashboard.className = "results-dashboard";

  const overview = document.createElement("section");
  overview.className = "results-overview";
  const reportBlock = document.createElement("div");
  reportBlock.className = "results-report";
  const reportText = document.createElement("div");
  const reportTitle = document.createElement("h3");
  text(reportTitle, "Run report");
  const reportDetail = document.createElement("span");
  text(reportDetail, report ? report.relative_path : "reports/run_summary.html");
  reportText.append(reportTitle, reportDetail);
  const reportActions = document.createElement("div");
  reportActions.className = "result-actions";
  if (report && report.view_url) reportActions.appendChild(resultAction(report.view_url, "Open report"));
  if (report && report.download_url) reportActions.appendChild(resultAction(report.download_url, "Download"));
  if (!report || report.status !== "ready") reportActions.appendChild(resultPill(report));
  reportBlock.append(reportText, reportActions);

  const count = document.createElement("div");
  count.className = "results-count";
  const countValue = document.createElement("strong");
  text(countValue, `${payload.ready_counts?.core || 0}/${coreFiles.length}`);
  const countLabel = document.createElement("span");
  text(countLabel, "tables ready");
  count.append(countValue, countLabel);
  overview.append(reportBlock, count);

  const tableHeader = document.createElement("div");
  tableHeader.className = "results-section-heading";
  const heading = document.createElement("h3");
  text(heading, "Download tables");
  const subheading = document.createElement("span");
  text(subheading, "Primary SIMPLseq output tables");
  tableHeader.append(heading, subheading);

  dashboard.append(overview, tableHeader, resultTable(coreFiles, "Core result tables are not listed yet."));

  list.appendChild(dashboard);
}

async function refreshAllRunState() {
  if (!activeOutdir) {
    resetRunDisplay();
    return {active: false};
  }
  const previousStatus = lastRunStatus;
  try {
    const status = await refreshStatus();
    latestStatusPayload = status;
    try {
      await refreshProgress(status);
    } catch (_error) {
      // Keep polling even if one progress read races the file writer.
    }
    try {
      await refreshResults();
    } catch (_error) {
      // Results may not exist yet while the pipeline is still active.
    }
    await loadLog({silent: true, statusPayload: status});
    const currentStatus = payloadStatus(status);
    const wasRunning = isActiveStatus(previousStatus);
    if (isActiveStatus(currentStatus) && !pollTimer) {
      startPolling();
    }
    if ((currentStatus === "complete" || currentStatus === "dry_run") && wasRunning) {
      const state = status.state || {};
      const key = `${currentRunOutdir()}:${state.completed_at || currentStatus}`;
      if (completedRedirectKey !== key) {
        completedRedirectKey = key;
        $("#results-outdir").value = currentRunOutdir();
        saveSettings();
        await refreshResults();
        selectTab("results");
        saveSettings();
      }
    }
    lastRunStatus = currentStatus;
    const active = Boolean(status.active);
    if (!active && !isActiveStatus(currentStatus) && pollTimer) {
      clearInterval(pollTimer);
      pollTimer = null;
    }
    return status;
  } catch (_error) {
    if (!isActiveStatus(lastRunStatus) && pollTimer) {
      clearInterval(pollTimer);
      pollTimer = null;
    }
    return null;
  }
}

function startPolling() {
  if (pollTimer) return;
  pollTimer = setInterval(refreshAllRunState, 3000);
}

async function loadLog() {
  const options = arguments[0] || {};
  if (logInFlight) return;
  const logNode = $("#technical-log");
  if (!activeOutdir) {
    renderTerminalLog(logNode, EMPTY_LOG_TEXT);
    return;
  }
  logInFlight = true;
  const out = encodeURIComponent(currentRunOutdir());
  const statusPayload = options.statusPayload || latestStatusPayload;
  const shouldStick =
    options.forceScroll ||
    followLog ||
    Boolean(statusPayload?.active);
  try {
    const status = payloadStatus(statusPayload);
    const hasRunState = Boolean(statusPayload?.state?.status);
    const canFetchLog =
      statusPayload?.active ||
      hasRunState ||
      isActiveStatus(status) ||
      status === "complete" ||
      status === "failed" ||
      status === "dry_run";
    if (!canFetchLog) {
      renderTerminalLog(logNode, EMPTY_LOG_TEXT);
      return;
    }
    const payload = await fetchJson(`/api/logs?out=${out}&max_bytes=120000`);
    renderTerminalLog(logNode, compactTerminalLog(payload.text, statusPayload));
    if (shouldStick && logNode) {
      scrollLogToBottom(logNode);
    }
  } catch (error) {
    if (!options.silent) {
      renderTerminalLog(logNode, error.message);
    }
  } finally {
    logInFlight = false;
  }
}

function renderCommonPaths(paths) {
  const box = $("#common-paths");
  if (!box) return;
  box.replaceChildren();
  paths.forEach((item) => {
    const button = document.createElement("button");
    text(button, item.label);
    button.title = item.path;
    button.addEventListener("click", () => {
      $("#fastq-dir").value = item.path;
      $("#browse-path").value = item.path;
      saveSettings();
      loadBrowse(item.path);
    });
    box.appendChild(button);
  });
}

function resetInstalledPath(selector, appRoot, fallback) {
  const node = $(selector);
  if (!node || !appRoot || !node.value) return;
  const normalizedValue = node.value.replaceAll("\\", "/").toLowerCase();
  const normalizedRoot = appRoot.replaceAll("\\", "/").toLowerCase();
  if (normalizedValue.startsWith(normalizedRoot)) {
    node.value = fallback;
  }
}

function currentRunOutdir() {
  return activeOutdir || $("#outdir").value || "results";
}

async function loadHealth() {
  const payload = await fetchJson("/api/health");
  pathStyle = payload.path_style || pathStyle;
  renderCommonPaths(payload.common_paths || []);
  if (payload.workspace_root && $("#workspace-root")) {
    text($("#workspace-root"), payload.workspace_root);
  }
  resetInstalledPath("#fastq-dir", payload.app_root, "data");
  resetInstalledPath("#samples-out", payload.app_root, "samples.csv");
  resetInstalledPath("#run-samples", payload.app_root, "samples.csv");
  resetInstalledPath("#outdir", payload.app_root, "results");
  resetInstalledPath("#results-outdir", payload.app_root, "results");
  if (!$("#browse-path").value) {
    $("#browse-path").value = $("#fastq-dir").value || payload.workspace_root || ".";
  }
  saveSettings();
}

async function loadBrowse(path) {
  const payload = await fetchJson(`/api/browse?path=${encodeURIComponent(path || $("#browse-path").value || ".")}`);
  $("#browse-path").value = payload.path;
  browseParent = payload.parent || "";
  const list = $("#browse-list");
  list.replaceChildren();
  const current = document.createElement("div");
  current.className = "browse-row";
  const meta = document.createElement("div");
  const title = document.createElement("strong");
  text(title, "Use this folder");
  const detail = document.createElement("span");
  text(detail, `${payload.path} | ${payload.pair_count || 0} pairs`);
  meta.append(title, detail);
  const action = document.createElement("button");
  text(action, "Select");
  action.addEventListener("click", () => {
    $("#fastq-dir").value = payload.path;
    saveSettings();
    setFolderMessage("Folder selected. Click Scan folder when ready.", "ok");
    closeFolderModal();
  });
  current.append(meta, action);
  list.appendChild(current);

  if (!payload.exists || payload.is_dir === false) {
    const div = document.createElement("div");
    div.className = "notice bad";
    text(div, "Folder is not available.");
    list.appendChild(div);
    return;
  }

  (payload.directories || []).forEach((item) => {
    const div = document.createElement("div");
    div.className = "browse-row";
    const info = document.createElement("div");
    const name = document.createElement("strong");
    text(name, item.name);
    const details = document.createElement("span");
    text(details, `${item.path} | ${item.fastq_files} FASTQ files`);
    info.append(name, details);
    const button = document.createElement("button");
    text(button, "Open");
  button.addEventListener("click", () => {
      $("#browse-path").value = item.path;
      loadBrowse(item.path);
    });
    div.append(info, button);
    list.appendChild(div);
  });
}

function bindEvents() {
  document.querySelectorAll(".tab").forEach((tab) => {
    tab.addEventListener("click", () => {
      selectTab(tab.dataset.tab);
      saveSettings();
    });
  });
  document.querySelectorAll("input").forEach((input) => {
    input.addEventListener("change", saveSettings);
  });
  $("#scan-button").addEventListener("click", scanFastqs);
  $("#check-button").addEventListener("click", runCheck);
  $("#run-button").addEventListener("click", startRun);
  $("#stop-button").addEventListener("click", stopRun);
  $("#refresh-status").addEventListener("click", () => withButtonFeedback($("#refresh-status"), "Checking...", refreshAllRunState));
  $("#refresh-results").addEventListener("click", () => {
    withButtonFeedback($("#refresh-results"), "Checking...", async () => {
      if (!$("#results-outdir").value) $("#results-outdir").value = currentRunOutdir();
      saveSettings();
      await refreshResults();
    });
  });
  $("#download-bundle").addEventListener("click", () => {
    const url = $("#download-bundle").dataset.url;
    if (url) window.location.href = url;
  });
  $("#load-log").addEventListener("click", () => withButtonFeedback($("#load-log"), "Checking...", loadLog));
  const logNode = $("#technical-log");
  if (logNode) {
    logNode.addEventListener("scroll", () => {
      followLog = logNode.scrollTop + logNode.clientHeight >= logNode.scrollHeight - 36;
    });
  }
  $("#browse-button").addEventListener("click", chooseFastqFolder);
  $("#choose-outdir-button").addEventListener("click", chooseOutputFolder);
  $("#browse-parent").addEventListener("click", () => {
    if (browseParent) loadBrowse(browseParent);
  });
  $("#browse-path").addEventListener("change", () => loadBrowse($("#browse-path").value));
  const folderClose = $(".folder-modal__close");
  if (folderClose) folderClose.addEventListener("click", closeFolderModal);
  const folderScrim = $(".folder-modal__scrim");
  if (folderScrim) folderScrim.addEventListener("click", closeFolderModal);
  $("#dry-run").addEventListener("change", () => {
    text($("#run-button"), $("#dry-run").checked ? "Preview command" : "Start run");
    saveSettings();
  });
  $("#outdir").addEventListener("change", () => {
    activeOutdir = "";
    saveSettings();
    refreshAllRunState();
  });
  $("#samples-out").addEventListener("change", () => {
    $("#run-samples").value = $("#samples-out").value;
    saveSettings();
  });
}


// ---------------------------------------------------------------------------
// DINEMITES analysis
// ---------------------------------------------------------------------------

let dinemitesPollTimer = null;
let dinemitesPlotItems = [];
let dinemitesPlotIndex = 0;
let dinemitesRunSummary = {};
let dinemitesSubjectRows = [];

function dinemitesOutdir() {
  return activeOutdir || $("#results-outdir").value || $("#outdir").value || "results";
}

function updateDinemitesRunButton() {
  const btn = $("#dinemites-run");
  if (!btn) return;
  const hasOutdir = !!(activeOutdir || $("#results-outdir").value);
  btn.disabled = !hasOutdir;
  btn.title = hasOutdir ? "" : "Run the main SIMPLseq pipeline first.";
}

function dinemitesTLagValue() {
  const noCutoff = $("#dinemites-no-day-cutoff")?.checked;
  if (noCutoff) return "Inf";
  return $("#dinemites-t-lag")?.value || "Inf";
}

function updateDinemitesModelSettingsVisibility() {
  const model = $("#dinemites-model")?.value || "simple";
  const simpleSettings = $("#dinemites-simple-settings");
  const bayesianSettings = $("#dinemites-bayesian-settings");
  if (simpleSettings) simpleSettings.hidden = model !== "simple";
  if (bayesianSettings) bayesianSettings.hidden = model !== "bayesian";
}

function handleDinemitesDayCutoffToggle() {
  const noCutoff = $("#dinemites-no-day-cutoff");
  const tLag = $("#dinemites-t-lag");
  if (tLag && noCutoff) {
    tLag.disabled = noCutoff.checked;
  }
  saveSettings();
}

function handleDinemitesToggle() {
  const enabled = $("#dinemites-enable").checked;
  const controls = $("#dinemites-controls");
  const pill = $("#dinemites-status");
  if (controls) controls.hidden = !enabled;
  updateDinemitesModelSettingsVisibility();
  if (enabled) {
    setPill(pill, "Ready", "ok");
    updateDinemitesRunButton();
  } else {
    setPill(pill, "Disabled", "");
  }
  saveSettings();
}

async function runDinemites() {
  const btn = $("#dinemites-run");
  const msg = $("#dinemites-message");
  btn.disabled = true;
  msg.className = "inline-message";
  text(msg, "Starting DINEMITES analysis...");
  try {
    const payload = await postJson("/api/dinemites/run", {
      model_type: $("#dinemites-model").value,
      outdir: dinemitesOutdir(),
      samples: $("#run-samples").value,
      n_lags: Number($("#dinemites-n-lags").value || 3),
      t_lag: dinemitesTLagValue(),
      min_abundance_pct: Number($("#dinemites-min-abundance-pct").value || 0.3),
      abundance_denominator: $("#dinemites-abundance-denominator").value,
      no_day_cutoff: $("#dinemites-no-day-cutoff").checked,
      seed: Number($("#dinemites-seed").value || 1),
      refresh: Number($("#dinemites-refresh-interval").value || 100),
      bayesian_lag_days: Number($("#dinemites-bayesian-lag-days").value || 30),
      bayesian_chains: Number($("#dinemites-bayesian-chains").value || 1),
      bayesian_parallel_chains: Number($("#dinemites-bayesian-parallel-chains").value || 1),
      bayesian_iter_warmup: Number($("#dinemites-bayesian-warmup").value || 500),
      bayesian_iter_sampling: Number($("#dinemites-bayesian-sampling").value || 500),
      bayesian_adapt_delta: Number($("#dinemites-bayesian-adapt-delta").value || 0.99),
      bayesian_drop_out: $("#dinemites-bayesian-drop-out").checked
    });
    text(msg, "DINEMITES analysis started.");
    msg.classList.add("ok");
    setPill($("#dinemites-status"), "Running", "warn");
    startDinemitesPolling();
  } catch (error) {
    text(msg, error.message);
    msg.classList.add("bad");
    setPill($("#dinemites-status"), "Failed", "bad");
    btn.disabled = false;
  }
}

function startDinemitesPolling() {
  if (dinemitesPollTimer) return;
  dinemitesPollTimer = setInterval(pollDinemitesStatus, 3000);
}

function stopDinemitesPolling() {
  if (dinemitesPollTimer) {
    clearInterval(dinemitesPollTimer);
    dinemitesPollTimer = null;
  }
}

async function pollDinemitesStatus() {
  const out = encodeURIComponent(dinemitesOutdir());
  try {
    const payload = await fetchJson(`/api/dinemites/status?out=${out}`);
    const status = payload.status || "idle";
    if (status === "running") {
      setPill($("#dinemites-status"), "Running", "warn");
    } else if (status === "complete") {
      setPill($("#dinemites-status"), "Complete", "ok");
      stopDinemitesPolling();
      await loadDinemitesResults();
      updateDinemitesRunButton();
    } else if (status === "failed") {
      setPill($("#dinemites-status"), "Failed", "bad");
      const detail = payload.state?.detail || "DINEMITES analysis failed.";
      const msg = $("#dinemites-message");
      msg.className = "inline-message bad";
      text(msg, detail);
      stopDinemitesPolling();
      updateDinemitesRunButton();
    } else {
      stopDinemitesPolling();
      updateDinemitesRunButton();
    }
  } catch (_error) {
    // Silently retry on transient errors.
  }
}

function renderDinemitesPlots(plots) {
  const gallery = $("#dinemites-plot-gallery");
  if (!gallery) return;
  const items = Array.isArray(plots) ? plots.filter((plot) => plot && plot.exists && plot.view_url) : [];
  const previousFilename = dinemitesPlotItems[dinemitesPlotIndex]?.filename || $("#dm-plot-selector")?.value || "";
  dinemitesPlotItems = items;
  dinemitesPlotIndex = Math.max(0, items.findIndex((plot) => plot.filename === previousFilename));
  if (dinemitesPlotIndex < 0) dinemitesPlotIndex = 0;

  updatePlotJump(
    gallery,
    $("#dm-plot-count"),
    null,
    items.length,
    "No DINEMITES plots available yet.",
    "DINEMITES plot",
    "DINEMITES plots"
  );

  updateDinemitesPlotBrowser();
  if (!items.length) {
    gallery.replaceChildren();
    gallery.hidden = true;
    updateDinemitesSelectedMetrics();
    return;
  }

  gallery.hidden = false;
  renderSelectedDinemitesPlot();
}

function dinemitesPlotLabel(plot) {
  return plot?.subject || plot?.filename || "Subject";
}

function updateDinemitesSelectedMetrics() {
  const plot = dinemitesPlotItems[dinemitesPlotIndex];
  const selectedSubject = dinemitesPlotLabel(plot);
  const subjectRow = dinemitesSubjectRows.find((item) => {
    return String(item?.subject || "").trim() === String(selectedSubject || "").trim();
  });
  const metricSource = subjectRow || dinemitesRunSummary || {};
  text($("#dm-new-infections"), formatNumber(metricSource.new_infections, 3));
  text($("#dm-molfoi"), formatNumber(metricSource.molfoi, 3));
}

function updateDinemitesPlotBrowser() {
  const selector = $("#dm-plot-selector");
  const prev = $("#dm-prev-plot");
  const next = $("#dm-next-plot");
  const hasPlots = dinemitesPlotItems.length > 0;

  if (selector) {
    selector.replaceChildren();
    if (!hasPlots) {
      const option = document.createElement("option");
      option.value = "";
      text(option, "No plots");
      selector.appendChild(option);
    } else {
      dinemitesPlotItems.forEach((plot, index) => {
        const option = document.createElement("option");
        option.value = plot.filename || String(index);
        text(option, dinemitesPlotLabel(plot));
        selector.appendChild(option);
      });
      selector.selectedIndex = dinemitesPlotIndex;
    }
    selector.disabled = !hasPlots;
  }

  if (prev) prev.disabled = !hasPlots || dinemitesPlotIndex <= 0;
  if (next) next.disabled = !hasPlots || dinemitesPlotIndex >= dinemitesPlotItems.length - 1;
}

function setDinemitesPlotIndex(index) {
  if (!dinemitesPlotItems.length) return;
  dinemitesPlotIndex = Math.min(dinemitesPlotItems.length - 1, Math.max(0, index));
  updateDinemitesPlotBrowser();
  renderSelectedDinemitesPlot();
}

function renderSelectedDinemitesPlot() {
  const gallery = $("#dinemites-plot-gallery");
  if (!gallery) return;
  gallery.replaceChildren();
  const plot = dinemitesPlotItems[dinemitesPlotIndex];
  if (!plot) {
    gallery.hidden = true;
    updateDinemitesSelectedMetrics();
    return;
  }
  gallery.hidden = false;
  updateDinemitesSelectedMetrics();

  const figure = document.createElement("figure");
  figure.className = "dinemites-plot-card";

  const img = document.createElement("img");
  img.src = plot.view_url;
  img.alt = `DINEMITES plot for ${dinemitesPlotLabel(plot)}`;
  img.loading = "lazy";
  figure.appendChild(img);

  const caption = document.createElement("figcaption");
  const title = document.createElement("span");
  text(title, dinemitesPlotLabel(plot));
  caption.appendChild(title);

  if (plot.download_url) {
    const button = document.createElement("button");
    button.type = "button";
    text(button, "Download plot");
    button.addEventListener("click", () => {
      window.location.href = plot.download_url;
    });
    caption.appendChild(button);
  }

  figure.appendChild(caption);
  gallery.appendChild(figure);
}

function renderDinemitesAlleleKey(rows) {
  const tbody = $("#dinemites-allele-key");
  if (!tbody) return;
  tbody.replaceChildren();

  const items = Array.isArray(rows) ? rows.filter((item) => item && item.short_allele_id) : [];
  if (!items.length) {
    tbody.appendChild(emptyRow(3, "No allele key available."));
    return;
  }

  items.forEach((item) => {
    const tr = row([
      item.short_allele_id,
      item.locus,
      item.allele
    ]);
    tr.title = item.allele || "";
    tr.addEventListener("click", async () => {
      if (!item.allele || !navigator.clipboard) return;
      try {
        await navigator.clipboard.writeText(item.allele);
      } catch (_error) {
        // Clipboard access can be unavailable in some browser contexts.
      }
    });
    tbody.appendChild(tr);
  });
}

async function loadDinemitesResults() {
  const out = encodeURIComponent(dinemitesOutdir());
  try {
    const payload = await fetchJson(`/api/dinemites/results?out=${out}`);
    const state = payload.state || {};
    const status = state.status || "idle";
    const resultsPanel = $("#dinemites-results");

    if (status === "complete") {
      if (resultsPanel) resultsPanel.hidden = false;
      setPill($("#dinemites-results-status"), "Complete", "ok");

      const summary = payload.summary || {};
      const subjects = payload.subjects || [];
      dinemitesRunSummary = summary;
      dinemitesSubjectRows = Array.isArray(subjects) ? subjects : [];
      text($("#dm-subjects"), formatNumber(summary.subjects, 0));
      text($("#dm-model"), state.model || "--");
      renderDinemitesPlots(payload.plots || []);
      renderDinemitesAlleleKey(payload.allele_key || []);

      // Render per-subject table
      const tbody = $("#dinemites-table");
      tbody.replaceChildren();
      if (!subjects.length) {
        tbody.appendChild(emptyRow(4, "No per-subject data available."));
      } else {
        subjects.forEach((item) => {
          tbody.appendChild(row([
            displayMissing(item.subject),
            formatNumber(item.new_infections, 2),
            formatNumber(item.molfoi, 2),
            displayMissing(item.time_points)
          ]));
        });
      }

      // Enable download buttons
      const files = payload.files || {};
      enableDinemitesDownload("#dm-dl-probabilities", files.allele_probabilities);
      enableDinemitesDownload("#dm-dl-allele-key", files.allele_key);
      enableDinemitesDownload("#dm-dl-molfoi", files.molfoi);
      enableDinemitesDownload("#dm-dl-new-infections", files.new_infections);
    } else if (status === "running") {
      setPill($("#dinemites-status"), "Running", "warn");
      startDinemitesPolling();
    } else if (status === "failed") {
      if (resultsPanel) resultsPanel.hidden = true;
      dinemitesRunSummary = {};
      dinemitesSubjectRows = [];
      renderDinemitesPlots([]);
      renderDinemitesAlleleKey([]);
      setPill($("#dinemites-status"), "Failed", "bad");
    }
  } catch (_error) {
    // Results may not exist yet.
    dinemitesRunSummary = {};
    dinemitesSubjectRows = [];
    renderDinemitesPlots([]);
  }
}

async function refreshDinemitesResults() {
  const msg = $("#dinemites-message");
  if (msg) {
    msg.className = "inline-message";
    text(msg, "Loading DINEMITES results...");
  }
  await loadDinemitesResults();
  if (msg) {
    text(msg, "DINEMITES results loaded.");
    msg.classList.add("ok");
  }
}

function enableDinemitesDownload(selector, fileInfo) {
  const btn = $(selector);
  if (!btn) return;
  if (fileInfo && fileInfo.exists && fileInfo.download_url) {
    btn.disabled = false;
    btn.onclick = () => { window.location.href = fileInfo.download_url; };
  } else {
    btn.disabled = true;
    btn.onclick = null;
  }
}

function bindDinemitesEvents() {
  const toggle = $("#dinemites-enable");
  if (toggle) toggle.addEventListener("change", handleDinemitesToggle);
  const runBtn = $("#dinemites-run");
  if (runBtn) runBtn.addEventListener("click", runDinemites);
  const refreshBtn = $("#dinemites-refresh");
  if (refreshBtn) refreshBtn.addEventListener("click", refreshDinemitesResults);
  const plotSelector = $("#dm-plot-selector");
  if (plotSelector) {
    plotSelector.addEventListener("change", () => {
      setDinemitesPlotIndex(plotSelector.selectedIndex);
    });
  }
  const prevPlot = $("#dm-prev-plot");
  if (prevPlot) prevPlot.addEventListener("click", () => setDinemitesPlotIndex(dinemitesPlotIndex - 1));
  const nextPlot = $("#dm-next-plot");
  if (nextPlot) nextPlot.addEventListener("click", () => setDinemitesPlotIndex(dinemitesPlotIndex + 1));
  const modelSelect = $("#dinemites-model");
  if (modelSelect) {
    modelSelect.addEventListener("change", () => {
      updateDinemitesModelSettingsVisibility();
      saveSettings();
    });
  }
  const nLagsInput = $("#dinemites-n-lags");
  if (nLagsInput) nLagsInput.addEventListener("change", saveSettings);
  const tLagInput = $("#dinemites-t-lag");
  if (tLagInput) tLagInput.addEventListener("change", saveSettings);
  const noDayCutoff = $("#dinemites-no-day-cutoff");
  if (noDayCutoff) noDayCutoff.addEventListener("change", handleDinemitesDayCutoffToggle);
  const abundanceInput = $("#dinemites-min-abundance-pct");
  if (abundanceInput) abundanceInput.addEventListener("change", saveSettings);
  const abundanceDenominator = $("#dinemites-abundance-denominator");
  if (abundanceDenominator) abundanceDenominator.addEventListener("change", saveSettings);
  const seedInput = $("#dinemites-seed");
  if (seedInput) seedInput.addEventListener("change", saveSettings);
  const refreshInput = $("#dinemites-refresh-interval");
  if (refreshInput) refreshInput.addEventListener("change", saveSettings);
  [
    "#dinemites-bayesian-lag-days",
    "#dinemites-bayesian-chains",
    "#dinemites-bayesian-parallel-chains",
    "#dinemites-bayesian-warmup",
    "#dinemites-bayesian-sampling",
    "#dinemites-bayesian-adapt-delta",
    "#dinemites-bayesian-drop-out"
  ].forEach((selector) => {
    const control = $(selector);
    if (control) control.addEventListener("change", saveSettings);
  });
}

// ---------------------------------------------------------------------------
// dcifer analysis
// ---------------------------------------------------------------------------

let dciferPollTimer = null;

function dciferOutdir() {
  return activeOutdir || $("#results-outdir").value || $("#outdir").value || "results";
}

function updateDciferRunButton() {
  const btn = $("#dcifer-run");
  if (!btn) return;
  const hasOutdir = !!(activeOutdir || $("#results-outdir").value);
  btn.disabled = !hasOutdir;
  btn.title = hasOutdir ? "" : "Run the main SIMPLseq pipeline first.";
}

function handleDciferToggle() {
  const toggle = $("#dcifer-enable");
  if (!toggle) return;
  const enabled = toggle.checked;
  const controls = $("#dcifer-controls");
  const pill = $("#dcifer-status");
  if (controls) controls.hidden = !enabled;
  if (enabled) {
    setPill(pill, "Ready", "ok");
    updateDciferRunButton();
  } else {
    setPill(pill, "Disabled", "");
  }
  saveSettings();
}

async function runDcifer() {
  const btn = $("#dcifer-run");
  const msg = $("#dcifer-message");
  if (btn) btn.disabled = true;
  if (msg) {
    msg.className = "inline-message";
    text(msg, "Starting dcifer analysis...");
  }
  try {
    await postJson("/api/dcifer/run", {
      outdir: dciferOutdir(),
      samples: $("#run-samples").value,
      min_abundance_pct: Number($("#dcifer-min-abundance-pct").value || 0.3),
      abundance_denominator: $("#dcifer-abundance-denominator").value,
      coi_lrank: Number($("#dcifer-coi-lrank").value || 2),
      ibd_grid_nr: Number($("#dcifer-ibd-grid-nr").value || 1000),
      alpha: Number($("#dcifer-alpha").value || 0.05),
      afreq_mode: "current_run"
    });
    if (msg) {
      text(msg, "dcifer analysis started.");
      msg.classList.add("ok");
    }
    setPill($("#dcifer-status"), "Running", "warn");
    startDciferPolling();
  } catch (error) {
    if (msg) {
      text(msg, error.message);
      msg.classList.add("bad");
    }
    setPill($("#dcifer-status"), "Failed", "bad");
    updateDciferRunButton();
  }
}

function startDciferPolling() {
  if (dciferPollTimer) return;
  dciferPollTimer = setInterval(pollDciferStatus, 3000);
}

function stopDciferPolling() {
  if (dciferPollTimer) {
    clearInterval(dciferPollTimer);
    dciferPollTimer = null;
  }
}

async function pollDciferStatus() {
  const out = encodeURIComponent(dciferOutdir());
  try {
    const payload = await fetchJson(`/api/dcifer/status?out=${out}`);
    const status = payload.status || "idle";
    if (status === "running") {
      setPill($("#dcifer-status"), "Running", "warn");
    } else if (status === "complete") {
      setPill($("#dcifer-status"), "Complete", "ok");
      stopDciferPolling();
      await loadDciferResults();
      updateDciferRunButton();
    } else if (status === "failed") {
      setPill($("#dcifer-status"), "Failed", "bad");
      const detail = payload.state?.detail || "dcifer analysis failed.";
      const msg = $("#dcifer-message");
      if (msg) {
        msg.className = "inline-message bad";
        text(msg, detail);
      }
      stopDciferPolling();
      updateDciferRunButton();
    } else {
      stopDciferPolling();
      updateDciferRunButton();
    }
  } catch (_error) {
    // Leave polling active for transient local server hiccups.
  }
}

function matrixHasValues(matrix) {
  return Boolean(
    matrix &&
    Array.isArray(matrix.labels) &&
    matrix.labels.length &&
    Array.isArray(matrix.rows) &&
    matrix.rows.length
  );
}

function colorFromStops(stops, value) {
  const clamped = Math.max(0, Math.min(1, value));
  const scaled = clamped * (stops.length - 1);
  const index = Math.min(stops.length - 2, Math.floor(scaled));
  const local = scaled - index;
  const left = stops[index];
  const right = stops[index + 1];
  const mix = left.map((channel, offset) => Math.round(channel + (right[offset] - channel) * local));
  return `rgb(${mix[0]}, ${mix[1]}, ${mix[2]})`;
}

function dciferHeatmapColor(kind, value) {
  if (value === null || value === undefined || Number.isNaN(Number(value))) return "#dedede";
  const number = Number(value);
  if (kind === "pvalue") {
    const intensity = 1 - Math.max(0, Math.min(1, number / 0.5));
    return colorFromStops([
      [247, 247, 252],
      [218, 218, 235],
      [158, 154, 200],
      [106, 81, 163],
      [63, 0, 125]
    ], intensity);
  }
  return colorFromStops([
    [247, 251, 255],
    [198, 219, 239],
    [107, 174, 214],
    [33, 113, 181],
    [8, 48, 107]
  ], Math.max(0, Math.min(1, number)));
}

function matrixValueText(kind, value) {
  if (value === null || value === undefined || Number.isNaN(Number(value))) return "NA";
  return kind === "pvalue" ? formatPValue(value) : formatNumber(value, 3);
}

function dciferPlotForKind(plots, kind) {
  const needle = kind === "pvalue" ? "pvalue" : "relatedness";
  return (Array.isArray(plots) ? plots : []).find((plot) => {
    const name = `${plot?.filename || ""} ${plot?.title || ""}`.toLowerCase();
    return name.includes(needle);
  });
}

function dciferSampleDisplayMap(labels, rows) {
  const ordered = [];
  const seen = new Set();
  const addSample = (value) => {
    const sampleId = String(value || "").trim();
    if (!sampleId || seen.has(sampleId)) return;
    seen.add(sampleId);
    ordered.push(sampleId);
  };
  labels.forEach(addSample);
  rows.forEach((rowItem) => addSample(rowItem?.sample_id));
  const digits = Math.max(2, String(ordered.length || 1).length);
  const entries = ordered.map((sampleId, index) => ({
    sampleId,
    displayId: `S${String(index + 1).padStart(digits, "0")}`,
  }));
  return {
    entries,
    bySample: new Map(entries.map((entry) => [entry.sampleId, entry.displayId])),
  };
}

function dciferDisplayId(labelMap, sampleId, fallbackIndex) {
  const full = String(sampleId || "").trim();
  if (labelMap.bySample.has(full)) return labelMap.bySample.get(full);
  const digits = Math.max(2, String(Math.max(labelMap.entries.length, fallbackIndex + 1)).length);
  return `S${String(fallbackIndex + 1).padStart(digits, "0")}`;
}

function renderDciferHeatmapKey(entries) {
  const details = document.createElement("details");
  details.className = "dcifer-heatmap-key";

  const summary = document.createElement("summary");
  text(summary, "Sample key");
  details.appendChild(summary);

  const list = document.createElement("div");
  list.className = "dcifer-heatmap-key-list";
  entries.forEach((entry) => {
    const item = document.createElement("div");
    item.className = "dcifer-heatmap-key-item";

    const displayId = document.createElement("strong");
    text(displayId, entry.displayId);
    item.appendChild(displayId);

    const sample = document.createElement("span");
    text(sample, entry.sampleId);
    sample.title = entry.sampleId;
    item.appendChild(sample);

    list.appendChild(item);
  });
  details.appendChild(list);
  return details;
}

function renderDciferHeatmapCard(kind, matrix, plotInfo) {
  const namespace = "http://www.w3.org/2000/svg";
  const labels = matrix.labels || [];
  const rows = matrix.rows || [];
  const labelMap = dciferSampleDisplayMap(labels, rows);
  const dimension = Math.max(labels.length, rows.length);
  const cell = dimension > 30 ? 16 : dimension > 18 ? 20 : 28;
  const left = dimension > 30 ? 56 : dimension > 18 ? 58 : 64;
  const top = dimension > 30 ? 56 : 62;
  const labelStep = Math.max(1, Math.ceil(dimension / 26));
  const labelFont = dimension > 30 ? 8 : dimension > 18 ? 9 : 10;
  const width = left + labels.length * cell + 18;
  const height = top + rows.length * cell + 30;
  const titleText = kind === "pvalue"
    ? "p-value heatmap"
    : "Relatedness heatmap";

  const figure = document.createElement("figure");
  figure.className = "dcifer-heatmap-card";

  const svg = document.createElementNS(namespace, "svg");
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  svg.setAttribute("width", String(width));
  svg.setAttribute("height", String(height));
  svg.setAttribute("role", "img");
  svg.setAttribute("aria-label", `${titleText}; hover cells for sample pair values`);

  const title = document.createElementNS(namespace, "text");
  title.setAttribute("x", "12");
  title.setAttribute("y", "24");
  title.setAttribute("fill", "#181817");
  title.setAttribute("font-family", "Inter, Arial, sans-serif");
  title.setAttribute("font-size", "16");
  title.setAttribute("font-weight", "800");
  title.textContent = titleText;
  svg.appendChild(title);

  const subtitle = document.createElementNS(namespace, "text");
  subtitle.setAttribute("x", "12");
  subtitle.setAttribute("y", "42");
  subtitle.setAttribute("fill", "#5e625b");
  subtitle.setAttribute("font-family", "Inter, Arial, sans-serif");
  subtitle.setAttribute("font-size", "10");
  subtitle.setAttribute("font-weight", "600");
  subtitle.textContent = "Axis IDs map to full sample names in the key";
  svg.appendChild(subtitle);

  labels.forEach((label, index) => {
    if (index % labelStep !== 0 && index !== labels.length - 1) return;
    const x = left + index * cell + cell / 2;
    const xLabel = document.createElementNS(namespace, "text");
    xLabel.setAttribute("x", String(x));
    xLabel.setAttribute("y", String(top - 8));
    xLabel.setAttribute("fill", "#494b44");
    xLabel.setAttribute("font-family", "Inter, Arial, sans-serif");
    xLabel.setAttribute("font-size", String(labelFont));
    xLabel.setAttribute("font-weight", "700");
    xLabel.setAttribute("text-anchor", "middle");
    xLabel.textContent = dciferDisplayId(labelMap, label, index);
    svg.appendChild(xLabel);
  });

  const tooltip = document.createElement("div");
  tooltip.className = "dcifer-heatmap-tooltip";
  text(tooltip, "Hover a cell to inspect sample pair values.");

  rows.forEach((rowItem, rowIndex) => {
    const rowLabel = String(rowItem.sample_id || "");
    const rowDisplayId = dciferDisplayId(labelMap, rowLabel, rowIndex);
    const y = top + rowIndex * cell;
    if (rowIndex % labelStep === 0 || rowIndex === rows.length - 1) {
      const yLabel = document.createElementNS(namespace, "text");
      yLabel.setAttribute("x", String(left - 8));
      yLabel.setAttribute("y", String(y + cell * 0.68));
      yLabel.setAttribute("fill", "#494b44");
      yLabel.setAttribute("font-family", "Inter, Arial, sans-serif");
      yLabel.setAttribute("font-size", String(labelFont));
      yLabel.setAttribute("font-weight", "700");
      yLabel.setAttribute("text-anchor", "end");
      yLabel.textContent = rowDisplayId;
      svg.appendChild(yLabel);
    }

    labels.forEach((columnLabel, columnIndex) => {
      const columnFullLabel = String(columnLabel || "");
      const columnDisplayId = dciferDisplayId(labelMap, columnFullLabel, columnIndex);
      const value = Array.isArray(rowItem.values) ? rowItem.values[columnIndex] : null;
      const rect = document.createElementNS(namespace, "rect");
      rect.classList.add("dcifer-heatmap-cell");
      rect.setAttribute("x", String(left + columnIndex * cell));
      rect.setAttribute("y", String(y));
      rect.setAttribute("width", String(cell));
      rect.setAttribute("height", String(cell));
      rect.setAttribute("fill", dciferHeatmapColor(kind, value));

      const hoverText = `${rowDisplayId} (${rowLabel}) vs ${columnDisplayId} (${columnFullLabel}): ${matrix.value_label || titleText} = ${matrixValueText(kind, value)}`;
      const nativeTitle = document.createElementNS(namespace, "title");
      nativeTitle.textContent = hoverText;
      rect.appendChild(nativeTitle);
      rect.setAttribute("tabindex", "0");
      rect.addEventListener("mouseenter", () => text(tooltip, hoverText));
      rect.addEventListener("focus", () => text(tooltip, hoverText));
      svg.appendChild(rect);
    });
  });

  figure.appendChild(svg);

  const caption = document.createElement("figcaption");
  const captionTitle = document.createElement("span");
  const truncatedNote = matrix.truncated ? ` (previewing ${labels.length} x ${rows.length})` : "";
  text(captionTitle, `${titleText}${truncatedNote}`);
  caption.appendChild(captionTitle);

  if (plotInfo?.download_url) {
    const button = document.createElement("button");
    button.type = "button";
    text(button, "Download original PNG");
    button.addEventListener("click", () => {
      window.location.href = plotInfo.download_url;
    });
    caption.appendChild(button);
  }

  figure.appendChild(caption);
  figure.appendChild(renderDciferHeatmapKey(labelMap.entries));
  figure.appendChild(tooltip);
  return figure;
}

function renderDciferPlots(plots, matrices = {}) {
  const gallery = $("#dcifer-plot-gallery");
  if (!gallery) return;
  gallery.replaceChildren();
  const items = Array.isArray(plots) ? plots.filter((plot) => plot && plot.exists && plot.view_url) : [];
  const matrixItems = [
    ["relatedness", matrices?.relatedness],
    ["pvalue", matrices?.pvalue],
  ].filter(([, matrix]) => matrixHasValues(matrix));
  const visibleCount = matrixItems.length || items.length;
  updatePlotJump(
    gallery,
    $("#dcifer-plot-count"),
    $("#dcifer-view-plots"),
    visibleCount,
    "No dcifer heatmaps available yet.",
    "dcifer heatmap",
    "dcifer heatmaps"
  );
  if (matrixItems.length) {
    gallery.hidden = false;
    const grid = document.createElement("div");
    grid.className = "dcifer-interactive-grid";
    matrixItems.forEach(([kind, matrix]) => {
      grid.appendChild(renderDciferHeatmapCard(kind, matrix, dciferPlotForKind(items, kind)));
    });
    gallery.appendChild(grid);
    return;
  }

  if (!items.length) {
    gallery.hidden = true;
    return;
  }
  gallery.hidden = false;
  items.forEach((plot) => {
    const figure = document.createElement("figure");
    figure.className = "dcifer-plot-card";

    const img = document.createElement("img");
    img.src = plot.view_url;
    img.alt = plot.title || plot.filename || "dcifer plot";
    img.loading = "lazy";
    figure.appendChild(img);

    const caption = document.createElement("figcaption");
    const title = document.createElement("span");
    text(title, plot.title || plot.filename || "dcifer plot");
    caption.appendChild(title);

    if (plot.download_url) {
      const button = document.createElement("button");
      button.type = "button";
      text(button, "Download original PNG");
      button.addEventListener("click", () => {
        window.location.href = plot.download_url;
      });
      caption.appendChild(button);
    }

    figure.appendChild(caption);
    gallery.appendChild(figure);
  });
}

function renderDciferPairs(pairs) {
  const tbody = $("#dcifer-pairs-table");
  if (!tbody) return;
  tbody.replaceChildren();
  const rows = Array.isArray(pairs) ? pairs : [];
  if (!rows.length) {
    tbody.appendChild(emptyRow(5, "No pairwise relatedness rows available."));
    return;
  }
  rows.forEach((item) => {
    tbody.appendChild(row([
      displayMissing(item.sample_a),
      displayMissing(item.sample_b),
      formatNumber(item.estimate, 3),
      formatPValue(item.p_value),
      displayMissing(item.comparison_type)
    ]));
  });
}

async function loadDciferResults() {
  const out = encodeURIComponent(dciferOutdir());
  try {
    const payload = await fetchJson(`/api/dcifer/results?out=${out}`);
    const state = payload.state || {};
    const status = state.status || "idle";
    const resultsPanel = $("#dcifer-results");

    if (status === "complete") {
      if (resultsPanel) resultsPanel.hidden = false;
      setPill($("#dcifer-results-status"), "Complete", "ok");
      const summary = payload.summary || {};
      text($("#dcifer-samples"), formatNumber(summary.samples, 0));
      text($("#dcifer-pairs"), formatNumber(summary.pairs, 0));
      text($("#dcifer-max-relatedness"), formatNumber(summary.max_relatedness, 3));
      text($("#dcifer-raw-p-le-alpha"), formatNumber(summary.raw_p_le_alpha, 0));
      renderDciferPlots(payload.plots || [], payload.matrices || {});
      renderDciferPairs(payload.pairs || []);

      const files = payload.files || {};
      enableDciferDownload("#dcifer-dl-pairs", files.pairwise_relatedness);
      enableDciferDownload("#dcifer-dl-coi", files.coi);
      enableDciferDownload("#dcifer-dl-input", files.input);
      enableDciferDownload("#dcifer-dl-matrix", files.relatedness_matrix);
    } else if (status === "running") {
      setPill($("#dcifer-status"), "Running", "warn");
      startDciferPolling();
    } else if (status === "failed") {
      if (resultsPanel) resultsPanel.hidden = true;
      renderDciferPairs([]);
      renderDciferPlots([]);
      setPill($("#dcifer-status"), "Failed", "bad");
    }
  } catch (_error) {
    // Results may not exist yet.
    renderDciferPlots([]);
  }
}

async function refreshDciferResults() {
  const msg = $("#dcifer-message");
  if (msg) {
    msg.className = "inline-message";
    text(msg, "Loading dcifer results...");
  }
  await loadDciferResults();
  if (msg) {
    text(msg, "dcifer results loaded.");
    msg.classList.add("ok");
  }
}

function enableDciferDownload(selector, fileInfo) {
  const btn = $(selector);
  if (!btn) return;
  if (fileInfo && fileInfo.exists && fileInfo.download_url) {
    btn.disabled = false;
    btn.onclick = () => { window.location.href = fileInfo.download_url; };
  } else {
    btn.disabled = true;
    btn.onclick = null;
  }
}

function bindDciferEvents() {
  const toggle = $("#dcifer-enable");
  if (toggle) toggle.addEventListener("change", handleDciferToggle);
  const runBtn = $("#dcifer-run");
  if (runBtn) runBtn.addEventListener("click", runDcifer);
  const refreshBtn = $("#dcifer-refresh");
  if (refreshBtn) refreshBtn.addEventListener("click", refreshDciferResults);
  [
    "#dcifer-min-abundance-pct",
    "#dcifer-abundance-denominator",
    "#dcifer-coi-lrank",
    "#dcifer-ibd-grid-nr",
    "#dcifer-alpha"
  ].forEach((selector) => {
    const control = $(selector);
    if (control) control.addEventListener("change", saveSettings);
  });
}

// ---------------------------------------------------------------------------
// Init
// ---------------------------------------------------------------------------

async function init() {
  restoreSettings();
  bindEvents();
  bindPlotWheelScrolling();
  bindDinemitesEvents();
  bindDciferEvents();
  selectTab("inputs");
  text($("#run-button"), $("#dry-run").checked ? "Preview command" : "Start run");
  try {
    await loadHealth();
  } catch (_error) {
    renderCommonPaths([]);
  }
  resetRunDisplay();
  const activeRun = await fetchActiveRun();
  if (activeRun?.active && activeRun.outdir) {
    activeOutdir = activeRun.outdir;
    $("#results-outdir").value = activeRun.outdir;
    saveSettings();
  }
  const status = await refreshAllRunState();
  if (status?.active) {
    selectTab("run");
    startPolling();
  }
  // Restore DINEMITES toggle state and check for existing results
  handleDinemitesDayCutoffToggle();
  updateDinemitesModelSettingsVisibility();
  handleDinemitesToggle();
  loadDinemitesResults();
  handleDciferToggle();
  loadDciferResults();
}

init();
