const cells = Array.from(document.querySelectorAll(".cell"));
const statusText = document.querySelector("#status");
const newRoundButton = document.querySelector("#new-round");
const resetScoreButton = document.querySelector("#reset-score");
const scoreX = document.querySelector("#score-x");
const scoreO = document.querySelector("#score-o");
const scoreDraw = document.querySelector("#score-draw");

const canvasButtons = new Map();
const scoreCards = new Map();

function roundedRect(context, x, y, width, height, radius) {
  const safeRadius = Math.min(radius, width / 2, height / 2);

  context.beginPath();
  context.roundRect(x, y, width, height, safeRadius);
}

function getButtonPalette(button) {
  if (button.classList.contains("reset-button") || button.classList.contains("new-round-button")) {
    return {
      top: "#26a69a",
      middle: "#0f766e",
      bottom: "#0b5d57",
      edge: "#074641",
      glow: "rgba(183, 255, 245, 0.8)",
      shadow: "rgba(5, 46, 44, 0.34)",
    };
  }

  if (button.classList.contains("winning")) {
    return {
      top: "#ecfff9",
      middle: "#c9f4e7",
      bottom: "#92d8c5",
      edge: "#3c9f88",
      glow: "rgba(255, 255, 255, 0.94)",
      shadow: "rgba(15, 118, 110, 0.28)",
    };
  }

  return {
    top: "#ffffff",
    middle: "#f4f7fb",
    bottom: "#dce3ee",
    edge: button.matches(":hover:not(:disabled)") ? "#0f766e" : "#aeb9c9",
    glow: "rgba(255, 255, 255, 0.96)",
    shadow: "rgba(30, 41, 59, 0.24)",
  };
}

function drawCanvasButton(button) {
  const canvas = canvasButtons.get(button);

  if (!canvas) {
    return;
  }

  const rect = button.getBoundingClientRect();

  if (!rect.width || !rect.height) {
    return;
  }

  const pixelRatio = Math.min(window.devicePixelRatio || 1, 2);
  const width = Math.round(rect.width * pixelRatio);
  const height = Math.round(rect.height * pixelRatio);

  if (canvas.width !== width || canvas.height !== height) {
    canvas.width = width;
    canvas.height = height;
  }

  const context = canvas.getContext("2d");
  context.setTransform(pixelRatio, 0, 0, pixelRatio, 0, 0);
  context.clearRect(0, 0, rect.width, rect.height);

  const palette = getButtonPalette(button);
  const isCell = button.classList.contains("cell");
  const depth = button.matches(":active:not(:disabled)") ? 3 : isCell ? 8 : 7;
  const inset = 2;
  const faceHeight = rect.height - depth - inset;
  const radius = isCell ? 12 : 10;

  context.save();
  context.shadowColor = palette.shadow;
  context.shadowBlur = isCell ? 9 : 7;
  context.shadowOffsetY = 4;
  roundedRect(context, inset, inset + depth, rect.width - inset * 2, faceHeight, radius);
  context.fillStyle = palette.edge;
  context.fill();
  context.restore();

  const sideGradient = context.createLinearGradient(0, faceHeight, 0, rect.height);
  sideGradient.addColorStop(0, palette.bottom);
  sideGradient.addColorStop(1, palette.edge);
  roundedRect(context, inset, inset + depth, rect.width - inset * 2, faceHeight, radius);
  context.fillStyle = sideGradient;
  context.fill();

  const faceGradient = context.createLinearGradient(0, inset, 0, faceHeight + inset);
  faceGradient.addColorStop(0, palette.top);
  faceGradient.addColorStop(0.52, palette.middle);
  faceGradient.addColorStop(1, palette.bottom);
  roundedRect(context, inset, inset, rect.width - inset * 2, faceHeight, radius);
  context.fillStyle = faceGradient;
  context.fill();

  context.lineWidth = isCell ? 1.5 : 1;
  context.strokeStyle = palette.edge;
  context.stroke();

  context.save();
  roundedRect(context, inset + 1.5, inset + 1.5, rect.width - inset * 2 - 3, faceHeight - 3, radius - 2);
  context.clip();
  const shine = context.createLinearGradient(0, inset, 0, faceHeight * 0.58);
  shine.addColorStop(0, palette.glow);
  shine.addColorStop(0.22, "rgba(255, 255, 255, 0.22)");
  shine.addColorStop(1, "rgba(255, 255, 255, 0)");
  context.fillStyle = shine;
  context.fillRect(inset, inset, rect.width - inset * 2, faceHeight * 0.6);
  context.restore();

  context.beginPath();
  context.moveTo(radius, faceHeight - 1);
  context.lineTo(rect.width - radius, faceHeight - 1);
  context.strokeStyle = "rgba(15, 23, 42, 0.14)";
  context.lineWidth = 1;
  context.stroke();
}

function addCanvasSurface(button) {
  let canvas = button.querySelector(".button-canvas");
  let label = button.querySelector(".button-label");

  if (!canvas || !label) {
    const buttonText = button.textContent;
    canvas = document.createElement("canvas");
    label = document.createElement("span");
    canvas.className = "button-canvas";
    canvas.setAttribute("aria-hidden", "true");
    label.className = "button-label";
    label.textContent = buttonText;
    button.replaceChildren(canvas, label);
  }

  button.dataset.canvasButton = "";
  canvasButtons.set(button, canvas);

  ["pointerenter", "pointerleave", "pointerdown", "pointerup", "focus", "blur"].forEach((eventName) => {
    button.addEventListener(eventName, () => requestAnimationFrame(() => drawCanvasButton(button)));
  });
}

document.querySelectorAll("button").forEach(addCanvasSurface);

function getScoreCardPalette(card) {
  if (card.classList.contains("score-card-x")) {
    return { accent: "#d94660", edge: "#9f2740", tint: "#fff0f3" };
  }

  if (card.classList.contains("score-card-o")) {
    return { accent: "#2563eb", edge: "#1746aa", tint: "#edf4ff" };
  }

  return { accent: "#0f766e", edge: "#0a554f", tint: "#eafaf7" };
}

function drawScoreCard(card) {
  const canvas = scoreCards.get(card);

  if (!canvas) {
    return;
  }

  const rect = card.getBoundingClientRect();

  if (!rect.width || !rect.height) {
    return;
  }

  const pixelRatio = Math.min(window.devicePixelRatio || 1, 2);
  const width = Math.round(rect.width * pixelRatio);
  const height = Math.round(rect.height * pixelRatio);

  if (canvas.width !== width || canvas.height !== height) {
    canvas.width = width;
    canvas.height = height;
  }

  const context = canvas.getContext("2d");
  const palette = getScoreCardPalette(card);
  const inset = 2;
  const depth = 9;
  const radius = 10;
  const faceHeight = rect.height - depth - inset;

  context.setTransform(pixelRatio, 0, 0, pixelRatio, 0, 0);
  context.clearRect(0, 0, rect.width, rect.height);

  context.save();
  context.shadowColor = "rgba(30, 41, 59, 0.26)";
  context.shadowBlur = 9;
  context.shadowOffsetY = 5;
  roundedRect(context, inset, inset + depth, rect.width - inset * 2, faceHeight, radius);
  context.fillStyle = palette.edge;
  context.fill();
  context.restore();

  const edgeGradient = context.createLinearGradient(0, faceHeight, 0, rect.height);
  edgeGradient.addColorStop(0, palette.accent);
  edgeGradient.addColorStop(1, palette.edge);
  roundedRect(context, inset, inset + depth, rect.width - inset * 2, faceHeight, radius);
  context.fillStyle = edgeGradient;
  context.fill();

  const faceGradient = context.createLinearGradient(0, inset, 0, faceHeight + inset);
  faceGradient.addColorStop(0, "#ffffff");
  faceGradient.addColorStop(0.58, "#fafcff");
  faceGradient.addColorStop(1, palette.tint);
  roundedRect(context, inset, inset, rect.width - inset * 2, faceHeight, radius);
  context.fillStyle = faceGradient;
  context.fill();
  context.strokeStyle = palette.accent;
  context.lineWidth = 1.35;
  context.stroke();

  context.save();
  roundedRect(context, inset + 2, inset + 2, rect.width - inset * 2 - 4, faceHeight - 4, radius - 2);
  context.clip();
  const highlight = context.createLinearGradient(0, inset, 0, faceHeight * 0.55);
  highlight.addColorStop(0, "rgba(255, 255, 255, 0.98)");
  highlight.addColorStop(1, "rgba(255, 255, 255, 0)");
  context.fillStyle = highlight;
  context.fillRect(inset, inset, rect.width - inset * 2, faceHeight * 0.55);
  context.restore();

  context.beginPath();
  context.moveTo(radius, faceHeight - 0.5);
  context.lineTo(rect.width - radius, faceHeight - 0.5);
  context.strokeStyle = "rgba(15, 23, 42, 0.16)";
  context.lineWidth = 1;
  context.stroke();
}

document.querySelectorAll(".score-card").forEach((card) => {
  scoreCards.set(card, card.querySelector(".score-card-canvas"));
});

const resizeObserver = new ResizeObserver((entries) => {
  entries.forEach(({ target }) => {
    if (canvasButtons.has(target)) {
      drawCanvasButton(target);
    } else {
      drawScoreCard(target);
    }
  });
});

canvasButtons.forEach((canvas, button) => {
  resizeObserver.observe(button);
  drawCanvasButton(button);
});

scoreCards.forEach((canvas, card) => {
  resizeObserver.observe(card);
  drawScoreCard(card);
});

function setCellMark(cell, mark) {
  cell.querySelector(".button-label").textContent = mark;
  drawCanvasButton(cell);
}

const winningCombos = [
  [0, 1, 2],
  [3, 4, 5],
  [6, 7, 8],
  [0, 3, 6],
  [1, 4, 7],
  [2, 5, 8],
  [0, 4, 8],
  [2, 4, 6],
];

let board = Array(9).fill("");
let currentPlayer = "X";
let roundOver = false;
let scores = {
  X: 0,
  O: 0,
  draw: 0,
};

function updateStatus(message) {
  statusText.textContent = message;
}

function updateScores() {
  scoreX.textContent = scores.X;
  scoreO.textContent = scores.O;
  scoreDraw.textContent = scores.draw;
}

function getWinner() {
  for (const combo of winningCombos) {
    const [a, b, c] = combo;

    if (board[a] && board[a] === board[b] && board[a] === board[c]) {
      return { player: board[a], combo };
    }
  }

  return null;
}

function finishRound(winner) {
  roundOver = true;

  if (winner) {
    scores[winner.player] += 1;
    winner.combo.forEach((index) => {
      cells[index].classList.add("winning");
      drawCanvasButton(cells[index]);
    });
    updateStatus(`Player ${winner.player} wins!`);
  } else {
    scores.draw += 1;
    updateStatus("It's a draw!");
  }

  cells.forEach((cell) => {
    cell.disabled = true;
  });
  updateScores();
}

function playTurn(index) {
  if (roundOver || board[index]) {
    return;
  }

  board[index] = currentPlayer;
  setCellMark(cells[index], currentPlayer);
  cells[index].classList.add(currentPlayer.toLowerCase());
  cells[index].disabled = true;

  const winner = getWinner();

  if (winner) {
    finishRound(winner);
    return;
  }

  if (board.every(Boolean)) {
    finishRound(null);
    return;
  }

  currentPlayer = currentPlayer === "X" ? "O" : "X";
  updateStatus(`Player ${currentPlayer}'s turn`);
}

function startNewRound() {
  board = Array(9).fill("");
  currentPlayer = "X";
  roundOver = false;
  updateStatus("Player X's turn");

  cells.forEach((cell) => {
    setCellMark(cell, "");
    cell.disabled = false;
    cell.classList.remove("x", "o", "winning");
    drawCanvasButton(cell);
  });
}

function resetScore() {
  scores = {
    X: 0,
    O: 0,
    draw: 0,
  };
  updateScores();
  startNewRound();
}

cells.forEach((cell, index) => {
  cell.addEventListener("click", () => playTurn(index));
});

newRoundButton.addEventListener("click", startNewRound);
resetScoreButton.addEventListener("click", resetScore);

updateScores();
