(function () {
  "use strict";
  const engine = MPQS.createEngine();
  const form = document.getElementById("factor-form");
  const composite = document.getElementById("composite");
  const threads = document.getElementById("threads");
  const ecmPretesting = document.getElementById("ecm-pretesting");
  const threadError = document.getElementById("thread-error");
  const factors = document.getElementById("factors");
  const factorHelp = document.getElementById("factor-help");
  const defaultFactorHelp = factorHelp.textContent;
  const factorize = document.getElementById("factorize");
  const clear = document.getElementById("clear");
  const digitCount = document.getElementById("digit-count");
  const inputError = document.getElementById("input-error");
  const status = document.getElementById("status");
  const statusArea = document.querySelector(".status-area");
  const elapsed = document.getElementById("elapsed");
  const relationSummary = document.getElementById("relation-summary");
  const relationProgress = document.getElementById("relations-progress");
  const relationPercentage = document.getElementById("relation-percentage");
  const relationCounts = { fullRelations: 0, smoothRelations: 0, combinedRelations: 0,
    partialsFound: 0, estimatedRelations: 0 };
  let activeWorker = null, clock = null, startedAt = 0;

  function setStatus(message, state = "ready") {
    status.textContent = message;
    statusArea.dataset.state = state;
  }

  function updateMetrics(values = {}) {
    for (const key of Object.keys(relationCounts)) {
      if (values[key] !== undefined) relationCounts[key] = values[key];
    }
    const { fullRelations, smoothRelations, combinedRelations, partialsFound, estimatedRelations } = relationCounts;
    relationSummary.textContent = `${fullRelations.toLocaleString()} relations ` +
      `(${smoothRelations.toLocaleString()} full + ${combinedRelations.toLocaleString()} combined from ` +
      `${partialsFound.toLocaleString()} partial), need ${estimatedRelations ? estimatedRelations.toLocaleString() : "—"}`;
    const percent = estimatedRelations > 0 ? Math.min(100, 100 * fullRelations / estimatedRelations) : 0;
    const label = `${Math.floor(percent * 10) / 10}%`;
    relationProgress.value = percent;
    relationProgress.textContent = label;
    relationPercentage.textContent = estimatedRelations > 0 ? label : "—";
    relationProgress.setAttribute("aria-valuetext", estimatedRelations > 0 ?
      `${fullRelations.toLocaleString()} usable relations out of an estimated ${estimatedRelations.toLocaleString()} (${label})` :
      "No relation estimate available");
  }

  function resetMetrics() {
    updateMetrics({ fullRelations: 0, smoothRelations: 0, combinedRelations: 0,
      partialsFound: 0, estimatedRelations: 0 });
    relationSummary.scrollLeft = 0;
  }

  function validate() {
    const value = composite.value.trim();
    const result = engine.validateInput(value);
    digitCount.textContent = `${value.length} / 110 digits`;
    inputError.textContent = value && !result.valid ? result.message : "";
    composite.setAttribute("aria-invalid", String(Boolean(value && !result.valid)));
    const validThreads = MPQS.validateThreads(threads.value);
    threadError.textContent = validThreads ? "" : `Enter a whole number from 1 to ${MPQS.MAX_THREADS}.`;
    threads.setAttribute("aria-invalid", String(!validThreads));
    factorize.disabled = !result.valid || !validThreads || activeWorker !== null;
    return { ...result, valid: result.valid && validThreads };
  }

  function stopWorker() {
    if (activeWorker) activeWorker.terminate();
    activeWorker = null;
    if (clock) clearInterval(clock);
    clock = null;
    composite.readOnly = false;
    threads.disabled = false;
    ecmPretesting.disabled = false;
    form.setAttribute("aria-busy", "false");
    validate();
  }

  function updateElapsed() {
    const seconds = (performance.now() - startedAt) / 1000;
    elapsed.textContent = seconds < 60 ? `${seconds.toFixed(1)} s` :
      `${Math.floor(seconds / 60)} min ${Math.floor(seconds % 60)} s`;
  }

  function formatFactors(list) {
    const runs = [];
    for (const value of list) {
      const last = runs[runs.length - 1];
      if (last && last.value === value) last.count++;
      else runs.push({ value, count: 1 });
    }
    return runs.map(({ value, count }) => count === 1 ? value : `${value}^${count}`).join(" × ");
  }

  function resetFactors() {
    factors.value = "";
    factors.scrollLeft = 0;
    factorHelp.textContent = defaultFactorHelp;
  }

  function begin(event) {
    event.preventDefault();
    const validation = validate();
    if (!validation.valid || activeWorker) return;
    resetFactors();
    resetMetrics();
    startedAt = performance.now();
    updateElapsed();
    setStatus("Trial division through 1000…", "working");
    try {
      composite.readOnly = true;
      threads.disabled = true;
      ecmPretesting.disabled = true;
      form.setAttribute("aria-busy", "true");
      factorize.disabled = true;
      clock = setInterval(updateElapsed, 250);
      const worker = MPQS.startFactorization(validation.text, Number(threads.value), data => {
        if (activeWorker !== worker) return;
        if (data.type === "progress") {
          const progress = data.progress;
          if (progress.factors) {
            factors.value = formatFactors(progress.factors);
            factorHelp.textContent = progress.factorSource ?
              "Partial result: factors found so far. The remaining cofactor is unresolved." :
              "Partial result: trial-division factors. The remaining cofactor is unresolved.";
          }
          updateMetrics(progress);
          const detail = progress.ecmTarget !== undefined ?
            ` · ${progress.digits}-digit cofactor · ${progress.ecmThreads} thread${progress.ecmThreads === 1 ? "" : "s"}` +
            ` · ${progress.ecmCurves.toLocaleString()} / ${progress.ecmTotalCurves.toLocaleString()} curves` +
            ` · t${progress.ecmDepth.toFixed(2)} / t${progress.ecmTarget.toFixed(2)}` +
            (progress.ecmB1 ? ` · B1 ${progress.ecmB1.toLocaleString()} · B2 ${progress.ecmB2.toLocaleString()}` : "") :
            progress.squfofIterations !== undefined ?
            ` · ${progress.digits}-digit cofactor · ${progress.squfofThreads} thread${progress.squfofThreads === 1 ? "" : "s"}` +
            ` · ${progress.squfofMultipliers.toLocaleString()} multipliers · ${progress.squfofIterations.toLocaleString()} iterations` :
            progress.rhoBudgetMs !== undefined ?
            ` · ${progress.digits}-digit cofactor · ${progress.rhoThreads || 1} thread${(progress.rhoThreads || 1) === 1 ? "" : "s"}` +
            ` · ${progress.rhoAttempt} attempt${progress.rhoAttempt === 1 ? "" : "s"}` +
            ` · ${progress.rhoEvaluations.toLocaleString()} evaluations · ` +
            `${Math.min(100, Math.floor(100 * progress.rhoElapsedMs / Math.max(0.001, progress.rhoBudgetMs)))}% of allowance` :
            progress.matrixColumns !== undefined ?
            ` · ${progress.matrixRows.toLocaleString()} rows × ${progress.matrixColumns.toLocaleString()} columns` +
            (progress.matrixThreads !== undefined ?
              ` · ${progress.matrixThreads} thread${progress.matrixThreads === 1 ? "" : "s"}` : "") +
            (progress.matrixIteration !== undefined ?
              ` · iteration ${progress.matrixIteration.toLocaleString()} · ` +
              `${Math.min(100, Math.floor(100 * progress.matrixDimensionsSolved / Math.max(1, progress.matrixDimensions)))}% of dimensions` +
              (progress.matrixAttempt > 1 ? ` · attempt ${progress.matrixAttempt}` : "") : "") :
            progress.polynomials !== undefined ?
            ` · ${progress.digits}-digit cofactor · ${progress.polynomials.toLocaleString()} polynomials` +
            ` · ${progress.threads || 1} thread${(progress.threads || 1) === 1 ? "" : "s"}` : "";
          setStatus(`${progress.stage}${detail}…`, "working");
        } else if (data.type === "result") {
          const result = data.result;
          updateElapsed();
          stopWorker();
          updateMetrics(result.stats);
          factorHelp.textContent = defaultFactorHelp;
          if (result.kind === "factored") {
            factors.value = formatFactors(result.factors);
            setStatus(`Complete. ${result.factors.length} prime or probable-prime factors, including repetitions. Product verified.`, "complete");
          } else {
            factors.value = "";
            setStatus(result.message, "ready");
          }
        } else if (data.type === "error") {
          updateElapsed();
          stopWorker();
          setStatus(`Unable to complete factorization: ${data.message}`, "error");
        }
      }, { ecmPretesting: ecmPretesting.checked });
      activeWorker = worker;
    } catch (error) {
      stopWorker();
      setStatus(`Unable to start the calculation: ${error.message}`, "error");
    }
  }

  form.addEventListener("submit", begin);
  threads.addEventListener("input", validate);
  composite.addEventListener("input", () => {
    stopWorker();
    resetFactors();
    resetMetrics();
    elapsed.textContent = "";
    setStatus("Ready for a composite number.");
  });
  composite.addEventListener("keydown", event => {
    if (event.key === "Enter") {
      event.preventDefault();
      if (!factorize.disabled) form.requestSubmit();
    }
  });
  clear.addEventListener("click", () => {
    stopWorker();
    composite.value = "";
    resetFactors();
    composite.scrollLeft = 0;
    for (const viewport of document.querySelectorAll(".field-viewport")) viewport.scrollLeft = 0;
    elapsed.textContent = "";
    resetMetrics();
    setStatus("Ready for a composite number.");
    validate();
    composite.focus();
  });
  window.addEventListener("pagehide", stopWorker);
  validate();
})();
