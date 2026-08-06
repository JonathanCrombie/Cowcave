(function () {
  'use strict';

  const M = globalThis.AurifeuilleanMath;
  const form = document.getElementById('calculatorForm');
  const inputA = document.getElementById('inputA');
  const inputB = document.getElementById('inputB');
  const inputN = document.getElementById('inputN');
  const calculateButton = document.getElementById('calculateButton');
  const clearButton = document.getElementById('clearButton');
  const canvas = document.getElementById('signToggle');
  const status = document.getElementById('status');
  const results = document.getElementById('results');
  const inputs = [inputA, inputB, inputN];
  let selectedSign = '+';
  let signHovered = false;
  let calculating = false;

  function validIntegerText(value, minimum) {
    if (!/^[0-9]+$/.test(value)) return false;
    try { return BigInt(value) >= minimum; }
    catch (_) { return false; }
  }

  function inputIsValid(input) {
    const minimum = input === inputB ? 1n : 2n;
    return validIntegerText(input.value.trim(), minimum);
  }

  function refreshValidity(showInvalid) {
    let allValid = true;
    for (const input of inputs) {
      const valid = inputIsValid(input);
      allValid = allValid && valid;
      input.classList.toggle('invalid', showInvalid && !valid && input.value.trim() !== '');
      input.setAttribute('aria-invalid', showInvalid && !valid ? 'true' : 'false');
    }
    calculateButton.disabled = calculating || !allValid;
    return allValid;
  }

  function clearResults() {
    results.replaceChildren();
    status.textContent = '';
    status.className = 'status';
  }

  function clearAll() {
    form.reset();
    inputA.value = '';
    inputB.value = '1';
    inputN.value = '';
    selectedSign = '+';
    canvas.setAttribute('aria-label', 'Plus selected. Click to change to minus.');
    canvas.setAttribute('aria-pressed', 'false');
    calculating = false;
    clearResults();
    drawSignToggle();
    refreshValidity(false);
    inputA.focus();
  }

  function drawSignToggle() {
    const rect = canvas.getBoundingClientRect();
    const ratio = window.devicePixelRatio || 1;
    canvas.width = Math.max(1, Math.round(rect.width * ratio));
    canvas.height = Math.max(1, Math.round(rect.height * ratio));
    const ctx = canvas.getContext('2d');
    ctx.scale(ratio, ratio);
    const width = rect.width;
    const height = rect.height;
    const radius = 14;

    ctx.clearRect(0, 0, width, height);
    ctx.beginPath();
    ctx.roundRect(1, 1, width - 2, height - 2, radius);
    const fill = ctx.createLinearGradient(0, 0, width, height);
    fill.addColorStop(0, signHovered ? '#244a52' : '#192c38');
    fill.addColorStop(1, signHovered ? '#59452a' : '#302a25');
    ctx.fillStyle = fill;
    ctx.fill();
    ctx.strokeStyle = signHovered ? '#62d6dc' : 'rgba(255,255,255,.18)';
    ctx.lineWidth = signHovered ? 2 : 1;
    ctx.stroke();

    ctx.fillStyle = selectedSign === '+' ? '#ffd58a' : '#62d6dc';
    ctx.font = '600 38px Georgia, serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.shadowColor = selectedSign === '+' ? 'rgba(255,189,89,.45)' : 'rgba(98,214,220,.45)';
    ctx.shadowBlur = signHovered ? 15 : 7;
    ctx.fillText(selectedSign === '+' ? '+' : '−', width / 2, height / 2 - 1);
  }

  function toggleSign() {
    selectedSign = selectedSign === '+' ? '-' : '+';
    canvas.setAttribute('aria-label', `${selectedSign === '+' ? 'Plus' : 'Minus'} selected. Click to change to ${selectedSign === '+' ? 'minus' : 'plus'}.`);
    canvas.setAttribute('aria-pressed', selectedSign === '-' ? 'true' : 'false');
    clearResults();
    drawSignToggle();
  }

  function escapeHtml(value) {
    return String(value).replace(/[&<>"]/g, character => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[character]));
  }

  function integer(value) {
    return `<div class="integer">${escapeHtml(value.toString())}</div>`;
  }

  function factorizationLatex(factors) {
    if (factors.length === 0) return '1';
    return factors.map(item => item.count === 1n ? item.prime.toString() : `${item.prime}^{${item.count}}`).join('\\cdot ');
  }

  function polynomialLatex(coefficients, degree, variable) {
    const terms = [];
    let index = 0n;
    while (index <= degree) {
      const coefficient = coefficients.get(index.toString()) ?? 0n;
      const power = degree - index;
      if (coefficient !== 0n) {
        const negative = coefficient < 0n;
        const magnitude = negative ? -coefficient : coefficient;
        let body = '';
        if (magnitude !== 1n || power === 0n) body += magnitude.toString();
        if (power > 0n) body += variable;
        if (power > 1n) body += `^{${power}}`;
        terms.push({ negative, body });
      }
      index += 1n;
    }
    if (terms.length === 0) return '0';
    return terms.map((term, position) => {
      if (position === 0) return `${term.negative ? '-' : ''}${term.body}`;
      return `${term.negative ? '-' : '+'}${term.body}`;
    }).join('');
  }

  function card(step, title, body) {
    return `<article class="result-card"><header><span class="step-number">${step}</span><h2>${title}</h2></header><div class="card-body">${body}</div></article>`;
  }

  function equation(latex) {
    return `<div class="equation">\\[${latex}\\]</div>`;
  }

  function factorBoxes(lLabel, lValue, mLabel, mValue) {
    return `<div class="factor-grid"><div class="factor-box"><strong>${lLabel}</strong>${integer(lValue)}</div><div class="factor-box"><strong>${mLabel}</strong>${integer(mValue)}</div></div>`;
  }

  function renderEligibility(data, oneBase) {
    const { a, b, exponent, sign, target, decomposition } = data;
    const expression = oneBase ? `${a}^{${exponent}}${sign}1` : `${a}^{${exponent}}${sign}${b}^{${exponent}}`;
    let body = `<p>The requested integer is ${equation(`${expression}=${target}`)}</p>`;
    body += `<div class="integer-label">Exact value</div>${integer(target)}`;

    if (oneBase) {
      body += equation(`${a}=${decomposition.squareRootPart}^{2}\\cdot ${decomposition.squareFreePart}`);
      body += `<p>The square-free part of the base is \\(r=${decomposition.squareFreePart}\\).</p>`;
    } else {
      body += equation(`A=${factorizationLatex(data.factorsA)},\\qquad B=${factorizationLatex(data.factorsB)}`);
      body += equation(`AB=${a * b}=${decomposition.squareRootPart}^{2}\\cdot ${decomposition.squareFreePart}`);
      body += `<p>The square class of \\(A/B\\), equivalently the square-free part of \\(AB\\), is \\(r=${decomposition.squareFreePart}\\).</p>`;
    }

    if (data.decomposition.squareFreePart === 1n) {
      body += `<div class="notice">The square-free part is 1, so this representation has no non-trivial Lucas–Aurifeuillean split.</div>`;
      return card(1, 'Eligibility and square-free decomposition', body);
    }

    const quotientText = data.h === null ? `${exponent}/${decomposition.squareFreePart}\\notin\\mathbb Z` : `${exponent}/${decomposition.squareFreePart}=${data.h}`;
    body += equation(`h=\\frac{N}{r}=${quotientText}`);
    body += `<p>The required sign for \\(r=${decomposition.squareFreePart}\\) is <strong>${escapeHtml(data.expectedSign)}</strong>: minus when \\(r\\equiv1\\pmod 4\\), plus otherwise.</p>`;

    if (!data.eligible) {
      let reason;
      if (!data.eligibleExponent) reason = data.h === null ? 'N is not divisible by r.' : 'N/r is not odd.';
      else reason = `The selected sign is ${sign}, but this square class requires ${data.expectedSign}.`;
      body += `<div class="notice">No Aurifeuillean \\(L,M\\) split occurs here. ${escapeHtml(reason)}</div>`;
    } else {
      body += `<div class="notice success-notice">Eligible: \\(N=rh\\) with odd \\(h=${data.h}\\), and the selected sign is the required one.</div>`;
    }
    return card(1, 'Eligibility and square-free decomposition', body);
  }

  function renderLucas(data) {
    const { lucas, decomposition } = data;
    const cPolynomial = polynomialLatex(lucas.gamma, lucas.degree, 'x');
    const dPolynomial = polynomialLatex(lucas.delta, lucas.degree - 1n, 'x');
    let body = `<p>Algorithm L uses \\(r'=${lucas.nPrime}\\) and degree \\(\\varphi(r')/2=${lucas.degree}\\). Its exact integer polynomials are</p>`;
    body += equation(`C_{${decomposition.squareFreePart}}(x)=${cPolynomial}`);
    body += equation(`D_{${decomposition.squareFreePart}}(x)=${dPolynomial}`);
    body += equation(`F_{${decomposition.squareFreePart}}(x)=C_{${decomposition.squareFreePart}}(x)^2-${decomposition.squareFreePart}xD_{${decomposition.squareFreePart}}(x)^2`);
    return card(2, 'Lucas polynomials', body);
  }

  function renderRaw(data, oneBase) {
    const { raw, decomposition } = data;
    let body;
    if (oneBase) {
      body = `<p>Set \\(x=A^h\\). The square root is integral because \\(A=u^2r\\) and \\(h\\) is odd.</p>`;
      body += equation(`x=${data.a}^{${raw.h}}=${raw.x}`);
      body += equation(`\\sqrt{rx}=\\sqrt{${decomposition.squareFreePart}\\cdot ${raw.x}}=${raw.radical}`);
      body += equation(`C_r(x)=${raw.cValue},\\qquad D_r(x)=${raw.dValue}`);
      body += equation(`L=${raw.cValue}-${raw.radical}\\cdot${raw.dValue}=${raw.L}`);
      body += equation(`M=${raw.cValue}+${raw.radical}\\cdot${raw.dValue}=${raw.M}`);
    } else {
      body = `<p>Homogenize at \\(X=A^h\\) and \\(Y=B^h\\), then clear the rational square root exactly.</p>`;
      body += equation(`X=${data.a}^{${raw.h}}=${raw.x},\\qquad Y=${data.b}^{${raw.h}}=${raw.y}`);
      body += equation(`T_h=\\sqrt{r(AB)^h}=${raw.radical}`);
      body += equation(`\\widehat C_{r,h}(A,B)=${raw.cValue},\\qquad \\widehat D_{r,h}(A,B)=${raw.dValue}`);
      body += equation(`L=${raw.cValue}-${raw.radical}\\cdot${raw.dValue}=${raw.L}`);
      body += equation(`M=${raw.cValue}+${raw.radical}\\cdot${raw.dValue}=${raw.M}`);
    }
    body += factorBoxes('Aurifeuillean L', raw.L, 'Aurifeuillean M', raw.M);
    body += equation(`L\\cdot M=${data.rawProduct}`);
    body += `<div class="integer-label">Remaining cofactor in the requested integer</div>${integer(data.rawCofactor)}`;
    return card(3, 'Aurifeuillean L and M factors', body);
  }

  function renderPrimitive(data, oneBase) {
    const rows = data.primitive.rows.map(row => {
      const action = row.mu === 1n ? 'multiply' : 'divide';
      const swap = row.chi === 1n ? 'keep L/M' : 'swap L/M';
      return `<tr><td>${row.divisor}</td><td>${row.raw.exponent}</td><td>${row.mu}</td><td>${row.chi}</td><td>${action}; ${swap}</td><td class="numeric">${row.raw.L}</td><td class="numeric">${row.raw.M}</td></tr>`;
    }).join('');

    let body = `<p>Möbius inversion removes lower-exponent components. The Jacobi symbol determines whether each conjugate pair keeps or swaps its \\(L,M\\) labels.</p>`;
    body += `<div class="table-wrap"><table><thead><tr><th>d</th><th>N/d</th><th>μ(d)</th><th>(r|d)</th><th>Action</th><th>Raw L</th><th>Raw M</th></tr></thead><tbody>${rows}</tbody></table></div>`;
    body += `<details><summary>Exact numerator and denominator products</summary><div>`;
    body += equation(`L_{prim}=\\frac{${data.primitive.lNumerator}}{${data.primitive.lDenominator}}=${data.primitive.L}`);
    body += equation(`M_{prim}=\\frac{${data.primitive.mNumerator}}{${data.primitive.mDenominator}}=${data.primitive.M}`);
    body += `</div></details>`;
    body += factorBoxes('Primitive L', data.primitive.L, 'Primitive M', data.primitive.M);
    const phiText = oneBase ? `\\Phi_{${data.cyclotomicIndex}}(${data.a})` : `\\Phi_{${data.cyclotomicIndex}}(${data.a},${data.b})`;
    body += equation(`L_{prim}M_{prim}=${phiText}=${data.primitiveProduct}`);
    body += `<div class="integer-label">Remaining cofactor in the requested integer</div>${integer(data.primitiveCofactor)}`;
    return card(4, 'Primitive L and M factors', body);
  }

  function renderResults(data) {
    const oneBase = data.b === 1n;
    let html = renderEligibility(data, oneBase);
    if (data.eligible) {
      html += renderLucas(data);
      html += renderRaw(data, oneBase);
      html += renderPrimitive(data, oneBase);
    }
    results.innerHTML = html;
    if (window.MathJax?.typesetPromise) window.MathJax.typesetPromise([results]).catch(() => {});
  }

  async function calculate() {
    if (!refreshValidity(true)) {
      clearAll();
      return;
    }
    calculating = true;
    refreshValidity(false);
    clearResults();
    status.textContent = 'Calculating exact factors…';
    await new Promise(resolve => window.setTimeout(resolve, 20));

    try {
      const a = BigInt(inputA.value.trim());
      const b = BigInt(inputB.value.trim());
      const exponent = BigInt(inputN.value.trim());
      const data = M.analyze(a, b, exponent, selectedSign);
      renderResults(data);
      status.textContent = data.eligible ? 'Exact calculation complete.' : 'Eligibility check complete.';
    } catch (error) {
      results.replaceChildren();
      status.textContent = `Calculation could not be completed: ${error.message}`;
      status.className = 'status error';
    } finally {
      calculating = false;
      refreshValidity(false);
    }
  }

  form.addEventListener('submit', event => { event.preventDefault(); calculate(); });
  clearButton.addEventListener('click', clearAll);
  for (const input of inputs) {
    input.addEventListener('input', () => {
      input.value = input.value.replace(/[^0-9]/g, '');
      clearResults();
      refreshValidity(false);
    });
    input.addEventListener('blur', () => refreshValidity(true));
  }

  canvas.addEventListener('click', toggleSign);
  canvas.addEventListener('keydown', event => {
    if (event.key === 'Enter' || event.key === ' ') {
      event.preventDefault();
      toggleSign();
    }
  });
  canvas.addEventListener('pointerenter', () => { signHovered = true; drawSignToggle(); });
  canvas.addEventListener('pointerleave', () => { signHovered = false; drawSignToggle(); });
  window.addEventListener('resize', drawSignToggle);

  drawSignToggle();
  refreshValidity(false);
})();
