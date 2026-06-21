const cells = Array.from(document.querySelectorAll(".cell"));
const statusText = document.querySelector("#status");
const newRoundButton = document.querySelector("#new-round");
const resetScoreButton = document.querySelector("#reset-score");
const scoreX = document.querySelector("#score-x");
const scoreO = document.querySelector("#score-o");
const scoreDraw = document.querySelector("#score-draw");

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
    winner.combo.forEach((index) => cells[index].classList.add("winning"));
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
  cells[index].textContent = currentPlayer;
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
    cell.textContent = "";
    cell.disabled = false;
    cell.className = "cell";
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
