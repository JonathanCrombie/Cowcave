const baseInput = document.getElementById("base");
const exponentInput = document.getElementById("exponent");
const modulusInput = document.getElementById("modulus");
const resultInput = document.getElementById("result");
const message = document.getElementById("message");
const clearButton = document.getElementById("clear");
const computeButton = document.getElementById("compute");
const form = document.getElementById("modular-form");

const integerPattern = /^[+-]?\d+$/;

function isIntegerText(value) {
  return integerPattern.test(value.trim());
}

function normalizeMod(value, modulus) {
  return ((value % modulus) + modulus) % modulus;
}

function binaryExponentiation(base, exponent, modulus) {
  let result = 1n;
  let currentBase = normalizeMod(base, modulus);
  let currentExponent = exponent;

  while (currentExponent > 0n) {
    if (currentExponent % 2n === 1n) {
      result = (result * currentBase) % modulus;
    }

    currentBase = (currentBase * currentBase) % modulus;
    currentExponent = currentExponent / 2n;
  }

  return normalizeMod(result, modulus);
}

function validateForm() {
  const hasValidBase = isIntegerText(baseInput.value);
  const hasValidExponent = isIntegerText(exponentInput.value);
  const hasValidModulus = isIntegerText(modulusInput.value);

  computeButton.disabled = !(hasValidBase && hasValidExponent && hasValidModulus);
  message.textContent = "";
}

function computeResult() {
  const base = BigInt(baseInput.value.trim());
  const exponent = BigInt(exponentInput.value.trim());
  const modulus = BigInt(modulusInput.value.trim());

  if (exponent < 0n) {
    resultInput.value = "";
    message.textContent = "Exponent must be a non-negative integer.";
    return;
  }

  if (modulus === 0n) {
    resultInput.value = "";
    message.textContent = "Modulus must not be zero.";
    return;
  }

  resultInput.value = binaryExponentiation(base, exponent, modulus).toString();
  message.textContent = "";
}

[baseInput, exponentInput, modulusInput].forEach((input) => {
  input.addEventListener("input", validateForm);
  input.addEventListener("keydown", (event) => {
    if (event.key === "Enter") {
      event.preventDefault();

      if (!computeButton.disabled) {
        computeResult();
      }
    }
  });
});

form.addEventListener("submit", (event) => {
  event.preventDefault();
  computeResult();
});

clearButton.addEventListener("click", () => {
  baseInput.value = "";
  exponentInput.value = "";
  modulusInput.value = "";
  resultInput.value = "";
  message.textContent = "";
  validateForm();
  baseInput.focus();
});

validateForm();
