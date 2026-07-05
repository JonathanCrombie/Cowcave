const valueA = document.querySelector("#value-a");
const valueB = document.querySelector("#value-b");
const valueGcd = document.querySelector("#value-gcd");
const valueM = document.querySelector("#value-m");
const valueN = document.querySelector("#value-n");
const clearButton = document.querySelector("#clear-button");
const computeButton = document.querySelector("#compute-button");

const integerPattern = /^[+-]?\d+$/;

function parseInteger(text) {
  const trimmed = text.trim();

  if (!integerPattern.test(trimmed)) {
    return null;
  }

  return BigInt(trimmed);
}

function absolute(value) {
  return value < 0n ? -value : value;
}

function extendedGcd(a, b) {
  if (a === 0n && b === 0n) {
    return { gcd: 0n, m: 0n, n: 0n };
  }

  let oldR = absolute(a);
  let r = absolute(b);
  let oldS = 1n;
  let s = 0n;
  let oldT = 0n;
  let t = 1n;

  while (r !== 0n) {
    const quotient = oldR / r;

    const nextR = oldR - quotient * r;
    oldR = r;
    r = nextR;

    const nextS = oldS - quotient * s;
    oldS = s;
    s = nextS;

    const nextT = oldT - quotient * t;
    oldT = t;
    t = nextT;
  }

  return {
    gcd: oldR,
    m: a < 0n ? -oldS : oldS,
    n: b < 0n ? -oldT : oldT,
  };
}

function clearOutputs() {
  valueGcd.value = "";
  valueM.value = "";
  valueN.value = "";
}

function updateComputeButton() {
  const a = parseInteger(valueA.value);
  const b = parseInteger(valueB.value);

  computeButton.disabled = a === null || b === null;
  clearOutputs();
}

function compute() {
  const a = parseInteger(valueA.value);
  const b = parseInteger(valueB.value);

  if (a === null || b === null) {
    return;
  }

  const result = extendedGcd(a, b);

  valueGcd.value = result.gcd.toString();
  valueM.value = result.m.toString();
  valueN.value = result.n.toString();
}

function clearAll() {
  valueA.value = "";
  valueB.value = "";
  clearOutputs();
  computeButton.disabled = true;
  valueA.focus();
}

valueA.addEventListener("input", updateComputeButton);
valueB.addEventListener("input", updateComputeButton);
clearButton.addEventListener("click", clearAll);
computeButton.addEventListener("click", compute);
