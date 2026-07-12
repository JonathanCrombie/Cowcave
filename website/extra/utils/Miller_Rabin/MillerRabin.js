const numberInput = document.getElementById("numberInput");
const iterationsSlider = document.getElementById("iterationsSlider");
const iterationsOutput = document.getElementById("iterationsOutput");
const currentIterationOutput = document.getElementById("currentIterationOutput");
const resultOutput = document.getElementById("resultOutput");
const testButton = document.getElementById("testButton");
const form = document.getElementById("millerRabinForm");

function isNumericInput(value) {
  return /^\s*\d+\s*$/.test(value);
}

function updateButtonState() {
  testButton.disabled = !isNumericInput(numberInput.value);
}

function modularPow(base, exponent, modulus) {
  let result = 1n;
  let currentBase = base % modulus;
  let currentExponent = exponent;

  while (currentExponent > 0n) {
    if (currentExponent % 2n === 1n) {
      result = (result * currentBase) % modulus;
    }

    currentExponent = currentExponent / 2n;
    currentBase = (currentBase * currentBase) % modulus;
  }

  return result;
}

function isPrimeByTrialDivision(n) {
  for (let divisor = 2n; divisor * divisor <= n; divisor = divisor + 1n) {
    if (n % divisor === 0n) {
      return false;
    }
  }

  return true;
}

function getWitness(iteration, n) {
  const range = n - 3n;

  return 2n + ((iteration - 1n) % range);
}

function waitForPaint() {
  return new Promise((resolve) => {
    requestAnimationFrame(resolve);
  });
}

async function runMillerRabin(n, maxIterations) {
  if (n < 2n) {
    currentIterationOutput.value = "0";
    return false;
  }

  if (n === 2n || n === 3n) {
    currentIterationOutput.value = "0";
    return true;
  }

  if (n % 2n === 0n) {
    currentIterationOutput.value = "0";
    return false;
  }

  let d = n - 1n;
  let r = 0n;

  while (d % 2n === 0n) {
    d = d / 2n;
    r = r + 1n;
  }

  for (let iteration = 1n; iteration <= maxIterations; iteration = iteration + 1n) {
    currentIterationOutput.value = String(iteration);
    await waitForPaint();

    const a = getWitness(iteration, n);
    let x = modularPow(a, d, n);

    if (x === 1n || x === n - 1n) {
      continue;
    }

    let mayBePrime = false;
    for (let round = 1n; round < r; round = round + 1n) {
      x = (x * x) % n;

      if (x === n - 1n) {
        mayBePrime = true;
        break;
      }
    }

    if (!mayBePrime) {
      return false;
    }
  }

  return true;
}

iterationsSlider.addEventListener("input", () => {
  iterationsOutput.value = iterationsSlider.value;
});

numberInput.addEventListener("input", () => {
  updateButtonState();
  currentIterationOutput.value = "0";
  resultOutput.value = "";
});

form.addEventListener("submit", async (event) => {
  event.preventDefault();

  if (!isNumericInput(numberInput.value)) {
    return;
  }

  currentIterationOutput.value = "0";
  resultOutput.value = "";
  testButton.disabled = true;

  const n = BigInt(numberInput.value.trim());

  if (n < 2n) {
    updateButtonState();
    return;
  }

  if (n >= 2n && n < 100n) {
    resultOutput.value = isPrimeByTrialDivision(n) ? "Prime" : "Composite";
    updateButtonState();
    return;
  }

  const maxIterations = BigInt(iterationsOutput.value);
  const isProbablyPrime = await runMillerRabin(n, maxIterations);

  resultOutput.value = isProbablyPrime ? "Probably Prime" : "Composite";
  updateButtonState();
});

updateButtonState();
