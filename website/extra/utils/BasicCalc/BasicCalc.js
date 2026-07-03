const expressionEl = document.querySelector("#expression");
const resultEl = document.querySelector("#result");
const keysEl = document.querySelector(".keys");

let currentValue = "0";
let storedValue = null;
let pendingOperator = null;
let shouldResetDisplay = false;

const ERROR_VALUE = "Error";
const REMAINDER_SEPARATOR = " rem ";

const operatorLabels = {
  "/": "÷",
  "*": "×",
  "-": "−",
  "+": "+",
  "**": "^",
};

function formatBigInt(value) {
  return value.toLocaleString("en-US");
}

function parseDisplayValue() {
  const valueText = currentValue.includes(REMAINDER_SEPARATOR)
    ? currentValue.split(REMAINDER_SEPARATOR)[0]
    : currentValue;

  try {
    return BigInt(valueText);
  } catch {
    return null;
  }
}

function updateDisplay() {
  resultEl.value = currentValue;
  resultEl.textContent = currentValue;

  if (storedValue === null || pendingOperator === null) {
    expressionEl.textContent = "";
    return;
  }

  expressionEl.textContent = `${formatBigInt(storedValue)} ${operatorLabels[pendingOperator]}`;
}

function clearCalculator() {
  currentValue = "0";
  storedValue = null;
  pendingOperator = null;
  shouldResetDisplay = false;
  updateDisplay();
}

function currentDisplayHasRemainder() {
  return currentValue.includes(REMAINDER_SEPARATOR);
}

function deleteDigit() {
  if (
    shouldResetDisplay ||
    currentValue === ERROR_VALUE ||
    currentDisplayHasRemainder()
  ) {
    currentValue = "0";
    shouldResetDisplay = false;
  } else if (currentValue.length === 1 || (currentValue.startsWith("-") && currentValue.length === 2)) {
    currentValue = "0";
  } else {
    currentValue = currentValue.slice(0, -1);
  }

  updateDisplay();
}

function appendDigit(digit) {
  if (
    currentValue === ERROR_VALUE ||
    shouldResetDisplay ||
    currentDisplayHasRemainder()
  ) {
    currentValue = digit;
    shouldResetDisplay = false;
    updateDisplay();
    return;
  }

  if (currentValue === "0") {
    currentValue = digit;
  } else {
    currentValue += digit;
  }

  updateDisplay();
}

function normalizePastedNumber(rawValue) {
  const compactValue = rawValue.trim().replace(/[,\s]/g, "");

  if (!/^[+-]?\d+$/.test(compactValue)) {
    return null;
  }

  const signedValue = compactValue.startsWith("+") ? compactValue.slice(1) : compactValue;
  return BigInt(signedValue).toString();
}

function pasteNumber(rawValue) {
  const pastedValue = normalizePastedNumber(rawValue);

  if (pastedValue === null) {
    return false;
  }

  currentValue = pastedValue;
  shouldResetDisplay = false;
  updateDisplay();
  return true;
}

function formatResult(value) {
  return value.toString();
}

function calculationError() {
  return {
    value: null,
    displayValue: ERROR_VALUE,
  };
}

function calculate(firstValue, operator, secondValue) {
  switch (operator) {
    case "+": {
      const value = firstValue + secondValue;
      return {
        value,
        displayValue: formatResult(value),
      };
    }
    case "-": {
      const value = firstValue - secondValue;
      return {
        value,
        displayValue: formatResult(value),
      };
    }
    case "*": {
      const value = firstValue * secondValue;
      return {
        value,
        displayValue: formatResult(value),
      };
    }
    case "/": {
      if (secondValue === 0n) {
        return calculationError();
      }

      const remainder = firstValue % secondValue;
      const quotient = (firstValue - remainder) / secondValue;

      return {
        value: quotient,
        displayValue:
          remainder === 0n ? quotient.toString() : `${quotient}${REMAINDER_SEPARATOR}${remainder}`,
      };
    }
    case "**": {
      if (secondValue < 0n) {
        return calculationError();
      }

      const value = firstValue ** secondValue;
      return {
        value,
        displayValue: formatResult(value),
      };
    }
    default:
      return {
        value: secondValue,
        displayValue: formatResult(secondValue),
      };
  }
}

function chooseOperator(operator) {
  if (currentValue === ERROR_VALUE || currentDisplayHasRemainder()) {
    clearCalculator();
    return;
  }

  const inputValue = parseDisplayValue();
  if (inputValue === null) {
    currentValue = ERROR_VALUE;
    shouldResetDisplay = true;
    updateDisplay();
    return;
  }

  if (storedValue !== null && pendingOperator !== null && !shouldResetDisplay) {
    const result = calculate(storedValue, pendingOperator, inputValue);
    currentValue = result.displayValue;
    storedValue = result.value;

    if (storedValue === null) {
      pendingOperator = null;
      shouldResetDisplay = true;
      updateDisplay();
      return;
    }
  } else {
    storedValue = inputValue;
  }

  pendingOperator = operator;
  shouldResetDisplay = true;
  updateDisplay();
}

function runCalculation() {
  if (currentDisplayHasRemainder()) {
    clearCalculator();
    return;
  }

  if (storedValue === null || pendingOperator === null || currentValue === ERROR_VALUE) {
    return;
  }

  const secondValue = parseDisplayValue();
  if (secondValue === null) {
    currentValue = ERROR_VALUE;
    storedValue = null;
    pendingOperator = null;
    shouldResetDisplay = true;
    updateDisplay();
    return;
  }

  const result = calculate(storedValue, pendingOperator, secondValue);

  expressionEl.textContent = `${formatBigInt(storedValue)} ${operatorLabels[pendingOperator]} ${formatBigInt(secondValue)} =`;
  currentValue = result.displayValue;
  storedValue = null;
  pendingOperator = null;
  shouldResetDisplay = true;

  if (currentValue === ERROR_VALUE) {
    shouldResetDisplay = true;
  }

  resultEl.value = currentValue;
  resultEl.textContent = currentValue;
}

function handleAction(action) {
  if (action === "clear") {
    clearCalculator();
  } else if (action === "delete") {
    deleteDigit();
  } else if (action === "calculate") {
    runCalculation();
  }
}

keysEl.addEventListener("click", (event) => {
  const button = event.target.closest("button");

  if (!button) {
    return;
  }

  if (button.dataset.value) {
    if (button.classList.contains("key-operator")) {
      chooseOperator(button.dataset.value);
    } else {
      appendDigit(button.dataset.value);
    }
  } else if (button.dataset.action) {
    handleAction(button.dataset.action);
  }
});

window.addEventListener("keydown", (event) => {
  const key = event.key;

  if (/^\d$/.test(key)) {
    appendDigit(key);
  } else if (["+", "-", "*", "/"].includes(key)) {
    chooseOperator(key);
  } else if (key === "^") {
    chooseOperator("**");
  } else if (key === "Enter" || key === "=") {
    event.preventDefault();
    runCalculation();
  } else if (key === "Backspace") {
    deleteDigit();
  } else if (key === "Escape") {
    clearCalculator();
  }
});

window.addEventListener("paste", (event) => {
  const pastedText = event.clipboardData?.getData("text") ?? "";

  if (pasteNumber(pastedText)) {
    event.preventDefault();
  }
});

updateDisplay();
