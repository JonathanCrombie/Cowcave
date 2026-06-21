const boardElement = document.querySelector("#board");
const statusText = document.querySelector("#statusText");
const redCount = document.querySelector("#redCount");
const blackCount = document.querySelector("#blackCount");
const turnLabel = document.querySelector("#turnLabel");
const turnDot = document.querySelector("#turnDot");
const captureLabel = document.querySelector("#captureLabel");
const moveCount = document.querySelector("#moveCount");
const moveHistory = document.querySelector("#moveHistory");
const resetButton = document.querySelector("#resetButton");
const undoButton = document.querySelector("#undoButton");

const BOARD_SIZE = 8;
const PLAYERS = {
  red: {
    name: "Red",
    direction: -1,
    kingRow: 0,
  },
  black: {
    name: "Black",
    direction: 1,
    kingRow: 7,
  },
};

let board = createInitialBoard();
let currentPlayer = "red";
let selectedSquare = null;
let legalMoves = [];
let forcedPiece = null;
let history = [];
let moveLog = [];
let moveNumber = 0;
let gameOver = false;

function createInitialBoard() {
  const nextBoard = Array.from({ length: BOARD_SIZE }, () => Array(BOARD_SIZE).fill(null));

  for (let row = 0; row < BOARD_SIZE; row += 1) {
    for (let col = 0; col < BOARD_SIZE; col += 1) {
      if (!isPlayable(row, col)) {
        continue;
      }

      if (row < 3) {
        nextBoard[row][col] = { player: "black", king: false };
      }

      if (row > 4) {
        nextBoard[row][col] = { player: "red", king: false };
      }
    }
  }

  return nextBoard;
}

function render() {
  boardElement.innerHTML = "";
  const allMoves = getAllMoves(currentPlayer);
  const mustCapture = allMoves.some((move) => move.captures.length > 0);

  for (let row = 0; row < BOARD_SIZE; row += 1) {
    for (let col = 0; col < BOARD_SIZE; col += 1) {
      const square = document.createElement("button");
      const piece = board[row][col];
      const playable = isPlayable(row, col);
      const selected = selectedSquare && selectedSquare.row === row && selectedSquare.col === col;
      const move = legalMoves.find((item) => item.to.row === row && item.to.col === col);

      square.className = `square ${playable ? "dark" : "light"}`;
      square.type = "button";
      square.dataset.row = row;
      square.dataset.col = col;
      square.disabled = !playable || gameOver;
      square.setAttribute("aria-label", getSquareLabel(row, col, piece));

      if (selected) {
        square.classList.add("selected");
      }

      if (move) {
        square.classList.add(move.captures.length > 0 ? "capture" : "legal");
      }

      if (piece) {
        const pieceElement = document.createElement("span");
        pieceElement.className = `piece ${piece.player}${piece.king ? " king" : ""}`;
        pieceElement.setAttribute("aria-hidden", "true");
        square.appendChild(pieceElement);
      }

      square.addEventListener("click", () => handleSquareClick(row, col));
      boardElement.appendChild(square);
    }
  }

  updatePanel(allMoves, mustCapture);
}

function updatePanel(allMoves, mustCapture) {
  const redPieces = countPieces("red");
  const blackPieces = countPieces("black");
  redCount.textContent = redPieces;
  blackCount.textContent = blackPieces;
  turnLabel.textContent = PLAYERS[currentPlayer].name;
  turnDot.classList.toggle("red", currentPlayer === "red");
  turnDot.classList.toggle("black", currentPlayer === "black");
  moveCount.textContent = moveNumber;
  undoButton.disabled = history.length === 0;

  if (gameOver) {
    const winner = redPieces === 0 ? "Black" : blackPieces === 0 ? "Red" : PLAYERS[getOpponent(currentPlayer)].name;
    statusText.textContent = `${winner} wins`;
    captureLabel.textContent = "Game over";
    return;
  }

  if (allMoves.length === 0) {
    statusText.textContent = `${PLAYERS[getOpponent(currentPlayer)].name} wins`;
    captureLabel.textContent = "No legal moves";
    gameOver = true;
    return;
  }

  if (forcedPiece) {
    statusText.textContent = `${PLAYERS[currentPlayer].name} continues`;
    captureLabel.textContent = "Multi-jump";
    return;
  }

  statusText.textContent = `${PLAYERS[currentPlayer].name} to move`;
  captureLabel.textContent = mustCapture ? "Capture required" : "Quiet board";
}

function handleSquareClick(row, col) {
  if (gameOver || !isPlayable(row, col)) {
    return;
  }

  const clickedPiece = board[row][col];
  const chosenMove = legalMoves.find((move) => move.to.row === row && move.to.col === col);

  if (chosenMove) {
    applyMove(chosenMove);
    return;
  }

  if (clickedPiece && clickedPiece.player === currentPlayer) {
    selectPiece(row, col);
    return;
  }

  clearSelection();
  render();
}

function selectPiece(row, col) {
  if (forcedPiece && (forcedPiece.row !== row || forcedPiece.col !== col)) {
    return;
  }

  const pieceMoves = getMovesForPiece(row, col);
  const allMoves = getAllMoves(currentPlayer);
  const captureRequired = allMoves.some((move) => move.captures.length > 0);
  legalMoves = captureRequired ? pieceMoves.filter((move) => move.captures.length > 0) : pieceMoves;

  if (legalMoves.length === 0) {
    clearSelection();
    render();
    return;
  }

  selectedSquare = { row, col };
  render();
}

function applyMove(move) {
  saveSnapshot();

  const piece = board[move.from.row][move.from.col];
  board[move.from.row][move.from.col] = null;
  board[move.to.row][move.to.col] = piece;

  move.captures.forEach((capture) => {
    board[capture.row][capture.col] = null;
  });

  const becameKing = !piece.king && move.to.row === PLAYERS[piece.player].kingRow;
  if (becameKing) {
    piece.king = true;
  }

  appendMove(move, becameKing);

  const moreCaptures = !becameKing ? getMovesForPiece(move.to.row, move.to.col).filter((nextMove) => nextMove.captures.length > 0) : [];
  if (move.captures.length > 0 && moreCaptures.length > 0) {
    forcedPiece = { row: move.to.row, col: move.to.col };
    selectedSquare = forcedPiece;
    legalMoves = moreCaptures;
    render();
    return;
  }

  currentPlayer = getOpponent(currentPlayer);
  forcedPiece = null;
  clearSelection();
  render();
}

function appendMove(move, becameKing) {
  const notation = `${formatSquare(move.from)}-${move.captures.length ? "x" : ">"}${formatSquare(move.to)}${becameKing ? "=K" : ""}`;
  const lastMove = moveLog[moveLog.length - 1];

  if (move.captures.length > 0 && forcedPiece && lastMove && lastMove.player === currentPlayer && lastMove.turn === moveNumber) {
    lastMove.notation = `${lastMove.notation} ${notation}`;
  } else {
    moveNumber += 1;
    moveLog.push({
      player: currentPlayer,
      turn: moveNumber,
      notation,
    });
  }

  renderMoveHistory();
}

function renderMoveHistory() {
  moveHistory.innerHTML = "";

  moveLog.slice(-12).forEach((move) => {
    const item = document.createElement("li");
    item.textContent = `${PLAYERS[move.player].name}: ${move.notation}`;
    moveHistory.appendChild(item);
  });

  moveHistory.scrollTop = moveHistory.scrollHeight;
}

function getAllMoves(player) {
  const moves = [];

  for (let row = 0; row < BOARD_SIZE; row += 1) {
    for (let col = 0; col < BOARD_SIZE; col += 1) {
      const piece = board[row][col];

      if (piece && piece.player === player) {
        moves.push(...getMovesForPiece(row, col));
      }
    }
  }

  const captures = moves.filter((move) => move.captures.length > 0);
  return captures.length > 0 ? captures : moves;
}

function getMovesForPiece(row, col) {
  const piece = board[row][col];

  if (!piece) {
    return [];
  }

  const moves = [];
  const directions = getDirections(piece);

  directions.forEach(([rowDirection, colDirection]) => {
    const nextRow = row + rowDirection;
    const nextCol = col + colDirection;
    const jumpRow = row + rowDirection * 2;
    const jumpCol = col + colDirection * 2;

    if (!isInside(nextRow, nextCol)) {
      return;
    }

    const adjacent = board[nextRow][nextCol];

    if (!adjacent) {
      moves.push({
        from: { row, col },
        to: { row: nextRow, col: nextCol },
        captures: [],
      });
      return;
    }

    if (adjacent.player !== piece.player && isInside(jumpRow, jumpCol) && !board[jumpRow][jumpCol]) {
      moves.push({
        from: { row, col },
        to: { row: jumpRow, col: jumpCol },
        captures: [{ row: nextRow, col: nextCol }],
      });
    }
  });

  return moves;
}

function getDirections(piece) {
  const forward = PLAYERS[piece.player].direction;

  if (piece.king) {
    return [
      [-1, -1],
      [-1, 1],
      [1, -1],
      [1, 1],
    ];
  }

  return [
    [forward, -1],
    [forward, 1],
  ];
}

function saveSnapshot() {
  history.push({
    board: cloneBoard(board),
    currentPlayer,
    selectedSquare: selectedSquare ? { ...selectedSquare } : null,
    legalMoves: cloneMoves(legalMoves),
    forcedPiece: forcedPiece ? { ...forcedPiece } : null,
    moveLog: moveLog.map((move) => ({ ...move })),
    moveNumber,
    gameOver,
  });
}

function undoMove() {
  const snapshot = history.pop();

  if (!snapshot) {
    return;
  }

  board = snapshot.board;
  currentPlayer = snapshot.currentPlayer;
  selectedSquare = snapshot.selectedSquare;
  legalMoves = snapshot.legalMoves;
  forcedPiece = snapshot.forcedPiece;
  moveLog = snapshot.moveLog;
  moveNumber = snapshot.moveNumber;
  gameOver = snapshot.gameOver;
  renderMoveHistory();
  render();
}

function resetGame() {
  board = createInitialBoard();
  currentPlayer = "red";
  selectedSquare = null;
  legalMoves = [];
  forcedPiece = null;
  history = [];
  moveLog = [];
  moveNumber = 0;
  gameOver = false;
  renderMoveHistory();
  render();
}

function clearSelection() {
  selectedSquare = null;
  legalMoves = [];
}

function countPieces(player) {
  return board.flat().filter((piece) => piece && piece.player === player).length;
}

function cloneBoard(sourceBoard) {
  return sourceBoard.map((row) => row.map((piece) => (piece ? { ...piece } : null)));
}

function cloneMoves(moves) {
  return moves.map((move) => ({
    from: { ...move.from },
    to: { ...move.to },
    captures: move.captures.map((capture) => ({ ...capture })),
  }));
}

function getOpponent(player) {
  return player === "red" ? "black" : "red";
}

function isPlayable(row, col) {
  return (row + col) % 2 === 1;
}

function isInside(row, col) {
  return row >= 0 && row < BOARD_SIZE && col >= 0 && col < BOARD_SIZE;
}

function formatSquare(square) {
  const file = String.fromCharCode(65 + square.col);
  const rank = BOARD_SIZE - square.row;
  return `${file}${rank}`;
}

function getSquareLabel(row, col, piece) {
  const square = formatSquare({ row, col });

  if (!piece) {
    return `${square} empty`;
  }

  return `${square} ${PLAYERS[piece.player].name} ${piece.king ? "king" : "piece"}`;
}

resetButton.addEventListener("click", resetGame);
undoButton.addEventListener("click", undoMove);

renderMoveHistory();
render();
