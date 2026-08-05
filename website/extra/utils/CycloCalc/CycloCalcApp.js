(function () {
  "use strict";

  const Core = window.CyclotomicCore;
  const form = document.getElementById("calculator-form");
  const inputA = document.getElementById("input-a");
  const inputB = document.getElementById("input-b");
  const inputN = document.getElementById("input-n");
  const mirroredN = document.getElementById("mirrored-n");
  const normalizePerfectPowers = document.getElementById("normalize-perfect-powers");
  const signToggle = document.getElementById("sign-toggle");
  const canvas = document.getElementById("sign-canvas");
  const calculateButton = document.getElementById("calculate-button");
  const clearButton = document.getElementById("clear-button");
  const results = document.getElementById("results");

  let sign = "+";
  let signHovered = false;
  let signFocused = false;

  function integerText(value) {
    return /^\d+$/.test(value.trim());
  }

  function parsedInputs() {
    if (!integerText(inputA.value) || !integerText(inputB.value) || !integerText(inputN.value)) return null;
    const a = BigInt(inputA.value.trim());
    const b = BigInt(inputB.value.trim());
    const exponent = BigInt(inputN.value.trim());
    if (a < 2n || b < 1n || exponent < 1n) return null;
    return { a, b, exponent };
  }

  function updateValidity() {
    calculateButton.disabled = parsedInputs() === null;
  }

  function syncMirroredExponent() {
    mirroredN.textContent = inputN.value;
  }

  function roundedRectangle(context, x, y, width, height, radius) {
    const r = Math.min(radius, width / 2, height / 2);
    context.beginPath();
    context.moveTo(x + r, y);
    context.arcTo(x + width, y, x + width, y + height, r);
    context.arcTo(x + width, y + height, x, y + height, r);
    context.arcTo(x, y + height, x, y, r);
    context.arcTo(x, y, x + width, y, r);
    context.closePath();
  }

  function drawSign() {
    const rectangle = canvas.getBoundingClientRect();
    const scale = window.devicePixelRatio || 1;
    canvas.width = Math.max(1, Math.round(rectangle.width * scale));
    canvas.height = Math.max(1, Math.round(rectangle.height * scale));
    const context = canvas.getContext("2d");
    context.setTransform(scale, 0, 0, scale, 0, 0);
    context.clearRect(0, 0, rectangle.width, rectangle.height);

    const active = signHovered || signFocused;
    const gradient = context.createLinearGradient(0, 0, rectangle.width, rectangle.height);
    gradient.addColorStop(0, active ? "#0a968d" : "#087f78");
    gradient.addColorStop(1, active ? "#08746e" : "#075f5b");
    roundedRectangle(context, 1, 1, rectangle.width - 2, rectangle.height - 2, 15);
    context.fillStyle = gradient;
    context.fill();
    context.strokeStyle = active ? "#e2bd6e" : "rgba(255,255,255,.25)";
    context.lineWidth = active ? 2 : 1;
    context.stroke();

    context.fillStyle = active ? "#fff7e2" : "#ffffff";
    context.font = "600 31px Georgia, serif";
    context.textAlign = "center";
    context.textBaseline = "middle";
    context.fillText(sign === "+" ? "+" : "−", rectangle.width / 2, rectangle.height / 2 - 1);
  }

  function toggleSign() {
    sign = sign === "-" ? "+" : "-";
    signToggle.setAttribute("aria-label", `Change sign; current sign is ${sign === "+" ? "plus" : "minus"}`);
    signToggle.setAttribute("aria-pressed", sign === "+" ? "true" : "false");
    drawSign();
  }

  function clearAll() {
    inputA.value = "";
    inputB.value = "1";
    inputN.value = "";
    syncMirroredExponent();
    normalizePerfectPowers.checked = true;
    sign = "+";
    signToggle.setAttribute("aria-label", "Change sign; current sign is plus");
    signToggle.setAttribute("aria-pressed", "true");
    results.replaceChildren();
    results.hidden = true;
    updateValidity();
    drawSign();
    inputA.focus();
  }

  function n(value) {
    return value.toString();
  }

  function softBreakInteger(value) {
    const text = n(value);
    const signPrefix = text.startsWith("-") ? "−" : "";
    const digits = text.startsWith("-") ? text.slice(1) : text;
    return signPrefix + digits.replace(/(\d{5})(?=\d)/g, "$1<wbr>");
  }

  function phiHtml(index, a, b, simple) {
    const readable = simple && b === 1n
      ? `Phi ${n(index)} of ${n(a)}`
      : `Phi ${n(index)} of ${n(a)}, ${n(b)}`;
    const argumentsHtml = simple && b === 1n
      ? softBreakInteger(a)
      : `${softBreakInteger(a)},<wbr> ${softBreakInteger(b)}`;
    return `<span class="wrappable-phi" aria-label="${readable}">Φ<sub>${softBreakInteger(index)}</sub>( ${argumentsHtml} )<wbr></span>`;
  }

  function multiplyText(values) {
    return values.map(n).join(" × ");
  }

  function texProduct(values) {
    return values.map(n).join(" \\cdot ");
  }

  function classifiedNumber(value) {
    const classification = Core.classifyPrimality(value, 10n);
    if (classification.kind === "unit") return `<span class="verification-number">${n(value)}</span>`;
    return `<span class="verification-number ${classification.kind}" title="${classification.label}">${n(value)}</span>`;
  }

  function verificationEquation(values, productValue, originalValue) {
    const operators = `<span class="verification-operator" aria-hidden="true">·</span>`;
    const left = values.map(classifiedNumber).join(operators);
    const spoken = `${values.map(n).join(" times ")} equals ${n(productValue)} equals ${n(originalValue)}`;
    return `<div class="math-block verification-equation" aria-label="${spoken}">
      <div>${left}<span class="verification-operator">=</span>${classifiedNumber(productValue)}<span class="verification-operator">=</span>${classifiedNumber(originalValue)}</div>
    </div>`;
  }

  function phiTex(index, a, b, simple) {
    return simple && b === 1n ? `\\Phi_{${n(index)}}(${n(a)})` : `\\Phi_{${n(index)}}(${n(a)},${n(b)})`;
  }

  function displayMath(tex) {
    return `<div class="math-block">\\[${tex}\\]</div>`;
  }

  function inlineMath(tex) {
    return `\\(${tex}\\)`;
  }

  function perfectPowerCard(relation, simple, contextLabel) {
    const q = relation.index;
    const factorNames = relation.indices.map((index) => phiTex(index, relation.aRoot, relation.bRoot, simple));
    const left = phiTex(q, Core.pow(relation.aRoot, relation.exponent), Core.pow(relation.bRoot, relation.exponent), simple);
    const rule = `\\Phi_q(X^m,Y^m)=\\prod_{\\substack{e\\mid m\\\\\\gcd(q,m/e)=1}}\\Phi_{qe}(X,Y)`;
    const substitution = `${left}=${factorNames.join("\\cdot")}`;
    const numeric = `${n(Core.cyclotomic(q, Core.pow(relation.aRoot, relation.exponent), Core.pow(relation.bRoot, relation.exponent)))}=${texProduct(relation.factors)}`;
    return `<article class="identity-card">
      <p class="identity-label">Perfect-power composition${contextLabel ? ` · ${contextLabel}` : ""}</p>
      <h3>A common ${n(relation.exponent)}th power refines this factor</h3>
      <p class="identity-rule">Rule (shown before use):</p>
      ${displayMath(rule)}
      <p>Here ${inlineMath(`X=${n(relation.aRoot)}`)}, ${inlineMath(simple ? "Y=1" : `Y=${n(relation.bRoot)}`)}, ${inlineMath(`m=${n(relation.exponent)}`)}, and ${inlineMath(`q=${n(q)}`)}.</p>
      <div class="identity-use"><p><strong>Apply it:</strong></p>${displayMath(`${substitution}=${numeric}`)}</div>
    </article>`;
  }

  function balancedCard(relation, a, b, simple, contextLabel) {
    const value = Core.cyclotomic(3n, a, b);
    return `<article class="identity-card">
      <p class="identity-label">Cubic algebraic split${contextLabel ? ` · ${contextLabel}` : ""}</p>
      <h3>The square conditions produce two exact factors</h3>
      <p class="identity-rule">If ${inlineMath("Y=s^2")} and ${inlineMath("X+Y=3r^2")}, then:</p>
      ${displayMath("\\Phi_3(X,Y)=(3r^2-3rs+s^2)(3r^2+3rs+s^2)")}
      <div class="identity-use"><p><strong>Check and apply it:</strong> ${inlineMath(`s=${n(relation.s)}`)}, ${inlineMath(`r=${n(relation.r)}`)}.</p>
      ${displayMath(`${phiTex(3n, a, b, simple)}=${n(value)}=${n(relation.factors[0])}\\cdot${n(relation.factors[1])}`)}</div>
    </article>`;
  }

  function pellCard(relation, a, b, simple, contextLabel) {
    const value = Core.cyclotomic(3n, a, b);
    return `<article class="identity-card">
      <p class="identity-label">Pell-type cubic split${contextLabel ? ` · ${contextLabel}` : ""}</p>
      <h3>A square discriminant splits ${inlineMath("\\Phi_3")}</h3>
      <p class="identity-rule">If ${inlineMath("3X(X+Y)=k^2")}, then:</p>
      ${displayMath("\\Phi_3(X,Y)=(2X+Y-k)(2X+Y+k)")}
      <div class="identity-use"><p><strong>Check and apply it:</strong> ${inlineMath(`k=${n(relation.k)}`)}.</p>
      ${displayMath(`${phiTex(3n, a, b, simple)}=${n(value)}=${n(relation.factors[0])}\\cdot${n(relation.factors[1])}`)}</div>
    </article>`;
  }

  function shiftThreeCard(relation, simple) {
    const value = Core.cyclotomic(3n, relation.a, relation.b);
    return `<article class="identity-card">
      <p class="identity-label">Phi shift identity</p>
      <h3>The cubic minus factor is also a cubic plus factor</h3>
      <p class="identity-rule">Rule:</p>
      ${displayMath("\\Phi_3(X,Y)=X^2+XY+Y^2=\\Phi_6(X+Y,Y)")}
      <div class="identity-use"><p><strong>Apply it:</strong></p>
      ${displayMath(`${phiTex(3n, relation.a, relation.b, simple)}=${phiTex(6n, relation.shiftedA, relation.shiftedB, simple)}=${n(value)}`)}</div>
    </article>`;
  }

  function shiftSixCard(relation, simple) {
    const value = Core.cyclotomic(6n, relation.a, relation.b);
    const symmetry = relation.usedSymmetry ? " First use symmetry in the two arguments." : "";
    return `<article class="identity-card">
      <p class="identity-label">Phi shift identity</p>
      <h3>The cubic plus factor becomes a cubic minus factor</h3>
      <p class="identity-rule">Rule: ${inlineMath("\\Phi_6(X,Y)=X^2-XY+Y^2=\\Phi_3(X-Y,Y)")}.${symmetry}</p>
      ${displayMath("\\Phi_6(X,Y)=\\Phi_3(X-Y,Y)")}
      <div class="identity-use"><p><strong>Apply it:</strong></p>
      ${displayMath(`${phiTex(6n, relation.a, relation.b, simple)}=${phiTex(3n, relation.shiftedA, relation.shiftedB, simple)}=${n(value)}`)}</div>
    </article>`;
  }

  function nestedCard(relation, index, a, b, simple) {
    const representation = `${n(a)}=${phiTex(relation.r, relation.c, relation.d, simple)}`;
    const basePower = simple ? "" : `,\\quad ${n(b)}=${n(relation.d)}^{${n(relation.r - 1n)}}`;
    const rule = relation.form === "odd"
      ? "r\\equiv-1\\pmod{2t}\\Longrightarrow\\Phi_{2t}(C,D)\\mid\\Phi_t(\\Phi_r(C,D),D^{r-1})"
      : "r\\equiv-1\\pmod{2t}\\Longrightarrow\\Phi_t(C,D)\\mid\\Phi_{2t}(\\Phi_r(C,D),D^{r-1})";
    const whole = Core.cyclotomic(index, a, b);
    return `<article class="identity-card">
      <p class="identity-label">Nested phi composition</p>
      <h3>A cyclotomic input forces a divisor</h3>
      <p class="identity-rule">For odd primes ${inlineMath("r,t")}, the applicable composition rule is:</p>
      ${displayMath(rule)}
      <div class="identity-use"><p><strong>Recognize the input:</strong> ${inlineMath(`${representation}${basePower}`)}, with ${inlineMath(`r=${n(relation.r)}\\equiv-1\\pmod{${n(2n * relation.t)}}`)}.</p>
      ${displayMath(`${phiTex(index, a, b, simple)}=${n(whole)}=${n(relation.factor)}\\cdot${n(relation.cofactor)}`)}</div>
    </article>`;
  }

  function analysisCards(analyses, a, b, simple, contextLabel) {
    return analyses.map((analysis) => {
      if (analysis.type === "perfect-power") return perfectPowerCard(analysis, simple, contextLabel);
      if (analysis.type === "balanced") return balancedCard(analysis, a, b, simple, contextLabel);
      if (analysis.type === "pell") return pellCard(analysis, a, b, simple, contextLabel);
      return "";
    }).join("");
  }

  function relationshipCards(factor, calculation, simple) {
    let html = "";
    for (const relation of factor.relationships) {
      if (relation.type === "perfect-power") {
        html += perfectPowerCard(relation, simple, "inside the standard factor");
      } else if (relation.type === "balanced") {
        html += balancedCard(relation, calculation.reducedA, calculation.reducedB, simple, "directly applicable");
      } else if (relation.type === "pell") {
        html += pellCard(relation, calculation.reducedA, calculation.reducedB, simple, "directly applicable");
      } else if (relation.type === "shift-3-to-6") {
        html += shiftThreeCard(relation, simple);
      } else if (relation.type === "shift-6-to-3") {
        html += shiftSixCard(relation, simple);
        html += analysisCards(relation.analyses, relation.shiftedA, relation.shiftedB, simple, "after the shift");
      } else if (relation.type === "nested") {
        html += nestedCard(relation, factor.index, calculation.reducedA, calculation.reducedB, simple);
      }
    }
    return html;
  }

  function isGenuineSplit(values, expectedProduct) {
    if (!Array.isArray(values) || values.length < 2) return false;
    const nonUnits = values.filter((value) => value !== 1n && value !== -1n);
    return nonUnits.length >= 2 && Core.product(values) === expectedProduct;
  }

  function splitFromAnalyses(analyses, expectedProduct) {
    for (const analysis of analyses) {
      if (["perfect-power", "balanced", "pell"].includes(analysis.type) && isGenuineSplit(analysis.factors, expectedProduct)) {
        return analysis.factors;
      }
    }
    return null;
  }

  function nonStandardSplit(factor) {
    for (const relation of factor.relationships) {
      if (["perfect-power", "balanced", "pell"].includes(relation.type) && isGenuineSplit(relation.factors, factor.reducedValue)) {
        return relation.factors;
      }
      if (relation.type === "nested") {
        const factors = [relation.factor, relation.cofactor];
        if (isGenuineSplit(factors, factor.reducedValue)) return factors;
      }
      if (relation.type === "shift-6-to-3") {
        const shiftedSplit = splitFromAnalyses(relation.analyses, factor.reducedValue);
        if (shiftedSplit !== null) return shiftedSplit;
      }
    }
    return null;
  }

  function refinedVerification(factors) {
    const values = [];
    let foundExtraFactors = false;
    for (const factor of factors) {
      const split = nonStandardSplit(factor);
      if (split === null) {
        values.push(factor.value);
        continue;
      }
      foundExtraFactors = true;
      if (factor.scale !== 1n) values.push(factor.scale);
      values.push(...split);
    }
    return { values, foundExtraFactors };
  }

  function typesetResults() {
    if (window.MathJax && typeof window.MathJax.typesetPromise === "function") {
      window.MathJax.typesetClear([results]);
      window.MathJax.typesetPromise([results]).catch(function () {});
    }
  }

  function perfectPowerNormalization(parsed) {
    if (!normalizePerfectPowers.checked) return null;
    const power = Core.commonPerfectPower(parsed.a, parsed.b);
    if (power === null) return null;
    return {
      originalA: parsed.a,
      originalB: parsed.b,
      originalExponent: parsed.exponent,
      powerExponent: power.exponent,
      workingA: power.aRoot,
      workingB: power.bRoot,
      workingExponent: power.exponent * parsed.exponent
    };
  }

  function expressionTex(a, b, exponent, operator) {
    return b === 1n
      ? `${n(a)}^{${n(exponent)}}${operator}1`
      : `${n(a)}^{${n(exponent)}}${operator}${n(b)}^{${n(exponent)}}`;
  }

  function exactResultExpression(source, operator) {
    const spokenOperator = operator === "+" ? "plus" : "minus";
    const spoken = `${n(source.a)} to the power ${n(source.exponent)} ${spokenOperator} ${n(source.b)} to the power ${n(source.exponent)}`;
    const rightSide = source.b === 1n ? "1" : `${n(source.b)}^{${n(source.exponent)}}`;
    const spacedExpression = `${n(source.a)}^{${n(source.exponent)}}\\;${operator}\\;${rightSide}`;
    return `<div class="result-expression" tabindex="0" aria-label="${spoken}"><span class="result-expression-content">${inlineMath(spacedExpression)}</span></div>`;
  }

  function normalizationCard(normalization, operator) {
    const m = normalization.powerExponent;
    const originalN = normalization.originalExponent;
    const workingN = normalization.workingExponent;
    const originalExpression = expressionTex(normalization.originalA, normalization.originalB, originalN, operator);
    const workingExpression = expressionTex(normalization.workingA, normalization.workingB, workingN, operator);
    const rewritten = `\\left(${n(normalization.workingA)}^{${n(m)}}\\right)^{${n(originalN)}}${operator}\\left(${n(normalization.workingB)}^{${n(m)}}\\right)^{${n(originalN)}}`;
    return `<article class="result-card normalization-card">
      <h2>1. Perfect-power normalization</h2>
      <p>When both bases have the same perfect-power exponent, multiply that exponent into ${inlineMath("N")}:</p>
      ${displayMath("(X^m)^N\\mathbin{\\pm}(Y^m)^N=X^{mN}\\mathbin{\\pm}Y^{mN}")}
      <p><strong>Recognize the supplied bases:</strong></p>
      ${displayMath(`${n(normalization.originalA)}=${n(normalization.workingA)}^{${n(m)}},\\qquad${n(normalization.originalB)}=${n(normalization.workingB)}^{${n(m)}}`)}
      <p><strong>Apply the normalization:</strong></p>
      ${displayMath(`${originalExpression}=${rewritten}=${workingExpression}`)}
      <div class="working-values" aria-label="Normalized working values">
        <span><small>Working A</small><strong>${n(normalization.workingA)}</strong></span>
        <span><small>Working B</small><strong>${n(normalization.workingB)}</strong></span>
        <span><small>Working N</small><strong>${n(workingN)}</strong></span>
      </div>
    </article>`;
  }

  function render(calculation, source, normalization) {
    const simple = calculation.b === 1n;
    const operator = calculation.sign === "+" ? "+" : "-";
    const workingExpression = expressionTex(calculation.a, calculation.b, calculation.exponent, operator);
    const exactExpanded = source.b === 1n
      ? `${n(calculation.leftPower)} ${operator} 1 = ${n(calculation.value)}`
      : `${n(calculation.leftPower)} ${operator} ${n(calculation.rightPower)} = ${n(calculation.value)}`;
    const standardStep = normalization === null ? "1" : "2";
    const identityStep = normalization === null ? "2" : "3";

    let html = `<article class="result-card result-hero">
      <p class="result-kicker">Exact result</p>
      ${exactResultExpression(source, operator)}
      <p class="exact-value">${exactExpanded.replace(/-/g, "−")}</p>
    </article>`;

    if (normalization !== null) html += normalizationCard(normalization, operator);

    const indexCondition = calculation.sign === "-" ? "d\\mid N" : "d\\mid 2N,\\ d\\nmid N";
    const generalRule = calculation.sign === "-"
      ? (simple ? "A^N-1=\\prod_{d\\mid N}\\Phi_d(A)" : "A^N-B^N=\\prod_{d\\mid N}\\Phi_d(A,B)")
      : (simple ? "A^N+1=\\prod_{\\substack{d\\mid2N\\\\d\\nmid N}}\\Phi_d(A)" : "A^N+B^N=\\prod_{\\substack{d\\mid2N\\\\d\\nmid N}}\\Phi_d(A,B)");
    const symbolicFactors = calculation.factors.map((factor) => phiTex(factor.index, calculation.a, calculation.b, simple));
    const numericFactors = calculation.factors.map((factor) => factor.value);

    html += `<article class="result-card">
      <h2>${standardStep}. Standard cyclotomic decomposition</h2>
      <p>The index set is ${inlineMath(indexCondition)}. First use the standard rule:</p>
      ${displayMath(generalRule)}
      <p><strong>Substitute the inputs and evaluate every factor:</strong></p>
      ${displayMath(`${workingExpression}=${symbolicFactors.join("\\cdot")}=${texProduct(numericFactors)}=${n(calculation.value)}`)}
      <ul class="factor-list">${calculation.factors.map((factor) => {
        const scale = !simple && calculation.common > 1n
          ? `<small> = ${softBreakInteger(factor.scale)} × ${softBreakInteger(factor.reducedValue)} from the coprime core</small>` : "";
        return `<li class="factor-row"><span class="factor-name">${phiHtml(factor.index, calculation.a, calculation.b, simple)}</span><span class="factor-value">${softBreakInteger(factor.value)}${scale}</span></li>`;
      }).join("")}</ul>
    </article>`;

    const identityHtml = calculation.factors.map((factor) => relationshipCards(factor, calculation, simple)).join("");
    if (identityHtml !== "") {
      html += `<article class="result-card">
        <h2>${identityStep}. Extra relationships</h2>
        <div class="identity-stack">${identityHtml}</div>
      </article>`;
    }

    const productValue = Core.product(numericFactors);
    const refined = refinedVerification(calculation.factors);
    const refinedHtml = refined.foundExtraFactors
      ? `<div class="refined-check">
          <p><strong>Now substitute the extra factors found by the applicable identities:</strong></p>
          ${verificationEquation(refined.values, Core.product(refined.values), calculation.value)}
        </div>`
      : "";
    html += `<article class="result-card verification">
      <h2><span class="verification-mark">✓</span>Verification</h2>
      <div class="primality-legend" aria-label="Primality color legend">
        <span><i class="legend-swatch prime-swatch"></i>Prime or probable prime</span>
        <span><i class="legend-swatch composite-swatch"></i>Composite</span>
      </div>
      <p>The product of the displayed standard factors reconstructs the original integer:</p>
      ${verificationEquation(numericFactors, productValue, calculation.value)}
      ${refinedHtml}
      <p>These are cyclotomic factors, not a claim of complete prime factorization.</p>
    </article>`;

    results.innerHTML = html;
    results.hidden = false;
    typesetResults();
    results.scrollIntoView({ behavior: "smooth", block: "start" });
  }

  function showError(message) {
    results.innerHTML = `<article class="result-card error-card"><h2>Calculation stopped</h2><p>${message}</p></article>`;
    results.hidden = false;
  }

  function assumptionViolations(parsed) {
    const violations = [];
    if (parsed.a <= parsed.b) {
      violations.push(`<li><strong>A &gt; B &gt; 0</strong> was violated: A = ${n(parsed.a)} and B = ${n(parsed.b)}.</li>`);
    }
    const common = Core.gcd(parsed.a, parsed.b);
    if (common !== 1n) {
      violations.push(`<li><strong>gcd(A, B) = 1</strong> was violated: gcd(${n(parsed.a)}, ${n(parsed.b)}) = ${n(common)}.</li>`);
    }
    return violations;
  }

  function gcdExtractionHelp(parsed) {
    const common = Core.gcd(parsed.a, parsed.b);
    if (common === 1n) return "";
    const reducedA = parsed.a / common;
    const reducedB = parsed.b / common;
    const operator = sign === "+" ? "+" : "-";
    return `<div class="gcd-help">
      <h3>Extract the common factor</h3>
      <p>Let ${inlineMath("g=\\gcd(A,B)")}. At every homogeneous cyclotomic factor, use:</p>
      ${displayMath("\\Phi_n(gA,gB)=g^{\\varphi(n)}\\Phi_n(A,B)")}
      <p><strong>Substitute the supplied bases:</strong> ${inlineMath(`g=\\gcd(${n(parsed.a)},${n(parsed.b)})=${n(common)}`)}.</p>
      ${displayMath(`\\Phi_n(${n(parsed.a)},${n(parsed.b)})=${n(common)}^{\\varphi(n)}\\Phi_n(${n(reducedA)},${n(reducedB)})`)}
      <p>The corresponding normalization of the complete expression is:</p>
      ${displayMath(`${n(parsed.a)}^{${n(parsed.exponent)}}${operator}${n(parsed.b)}^{${n(parsed.exponent)}}=${n(common)}^{${n(parsed.exponent)}}\\left(${n(reducedA)}^{${n(parsed.exponent)}}${operator}${n(reducedB)}^{${n(parsed.exponent)}}\\right)`)}
    </div>`;
  }

  function showAssumptionViolations(violations, parsed) {
    results.innerHTML = `<article class="result-card error-card">
      <h2>Assumption violation</h2>
      <p>The calculation was not performed.</p>
      <ul class="assumption-list">${violations.join("")}</ul>
      ${gcdExtractionHelp(parsed)}
    </article>`;
    results.hidden = false;
    typesetResults();
    results.scrollIntoView({ behavior: "smooth", block: "start" });
  }

  form.addEventListener("submit", function (event) {
    event.preventDefault();
    const parsed = parsedInputs();
    if (parsed === null) {
      clearAll();
      return;
    }
    const violations = assumptionViolations(parsed);
    if (violations.length > 0) {
      showAssumptionViolations(violations, parsed);
      return;
    }
    calculateButton.disabled = true;
    calculateButton.textContent = "Calculating…";
    window.setTimeout(function () {
      try {
        const normalization = perfectPowerNormalization(parsed);
        const working = normalization === null
          ? parsed
          : { a: normalization.workingA, b: normalization.workingB, exponent: normalization.workingExponent };
        render(Core.calculate(working.a, working.b, working.exponent, sign), parsed, normalization);
      } catch (error) {
        showError(error instanceof Error ? error.message : "The calculation could not be completed.");
      } finally {
        calculateButton.textContent = "Calculate";
        updateValidity();
      }
    }, 0);
  });

  for (const input of [inputA, inputB, inputN]) input.addEventListener("input", updateValidity);
  inputN.addEventListener("input", syncMirroredExponent);
  signToggle.addEventListener("click", toggleSign);
  signToggle.addEventListener("mouseenter", function () { signHovered = true; drawSign(); });
  signToggle.addEventListener("mouseleave", function () { signHovered = false; drawSign(); });
  signToggle.addEventListener("focus", function () { signFocused = true; drawSign(); });
  signToggle.addEventListener("blur", function () { signFocused = false; drawSign(); });
  clearButton.addEventListener("click", clearAll);
  window.addEventListener("resize", drawSign);

  updateValidity();
  syncMirroredExponent();
  window.requestAnimationFrame(drawSign);
})();
