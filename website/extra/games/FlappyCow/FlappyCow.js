"use strict";

const LOGICAL_WIDTH = 1000;
const LOGICAL_HEIGHT = 620;
const GROUND_Y = 564;
const COW_X = 195;
const COW_WIDTH = 96;
const COW_HEIGHT = 66;
const OBSTACLE_WIDTH = 112;
const STORAGE_KEY = "cownoise-flappy-cow-progress-v1";

// Every gap and obstacle position is deliberately fixed. No randomness is used.
const LEVELS = [
  {
    name: "Sunny Stroll",
    blurb: "Wide gaps & a gentle breeze",
    intro: "Guide Bessie through 5 haystack gaps. This first pasture is nice and roomy.",
    speed: 132,
    gravity: 1120,
    flap: -405,
    gapHeight: 242,
    spacing: 435,
    gapCenters: [300, 255, 335, 275, 315],
  },
  {
    name: "Breezy Meadow",
    blurb: "A little quicker, a little tighter",
    intro: "The breeze is picking up. Stay light on the flap through 8 fixed haystacks.",
    speed: 144,
    gravity: 1160,
    flap: -410,
    gapHeight: 224,
    spacing: 410,
    gapCenters: [260, 330, 235, 355, 285, 215, 325, 270],
  },
  {
    name: "Twisty Trail",
    blurb: "Alternating gaps test your rhythm",
    intro: "Bessie's route gets twisty. Watch the next gap and settle into the rhythm.",
    speed: 156,
    gravity: 1190,
    flap: -415,
    gapHeight: 207,
    spacing: 380,
    gapCenters: [335, 225, 355, 245, 310, 205, 345, 260, 320, 220, 295],
  },
  {
    name: "Golden Gauntlet",
    blurb: "Faster stacks with less room",
    intro: "Fourteen golden gates stand between Bessie and the far pasture. Keep steady.",
    speed: 168,
    gravity: 1220,
    flap: -420,
    gapHeight: 190,
    spacing: 355,
    gapCenters: [235, 345, 255, 365, 280, 210, 335, 245, 355, 275, 220, 325, 260, 345],
  },
  {
    name: "Barnyard Legend",
    blurb: "The ultimate Bessie challenge",
    intro: "This is Bessie's biggest flight: 17 quick, narrow gaps. You know the route.",
    speed: 180,
    gravity: 1250,
    flap: -425,
    gapHeight: 176,
    spacing: 335,
    gapCenters: [340, 215, 360, 250, 315, 205, 350, 235, 330, 270, 205, 345, 245, 365, 280, 225, 325],
  },
];

const els = {
  stage: document.querySelector("#gameStage"),
  bessie: document.querySelector("#bessie"),
  obstacles: document.querySelector("#obstacleLayer"),
  particles: document.querySelector("#particleLayer"),
  finishLine: document.querySelector("#finishLine"),
  overlay: document.querySelector("#gameOverlay"),
  overlayEyebrow: document.querySelector("#overlayEyebrow"),
  overlayTitle: document.querySelector("#overlayTitle"),
  overlayCopy: document.querySelector("#overlayCopy"),
  controlHints: document.querySelector("#controlHints"),
  startButton: document.querySelector("#startButton"),
  startButtonLabel: document.querySelector("#startButton span:first-child"),
  secondaryButton: document.querySelector("#secondaryButton"),
  pauseButton: document.querySelector("#pauseButton"),
  hudLevel: document.querySelector("#hudLevel"),
  hudScore: document.querySelector("#hudScore"),
  progressFill: document.querySelector("#progressFill"),
  levelGrid: document.querySelector("#levelGrid"),
  resetProgressButton: document.querySelector("#resetProgressButton"),
  toast: document.querySelector("#toast"),
  gameStatus: document.querySelector("#gameStatus"),
};

let unlockedLevel = readProgress();
let selectedLevel = Math.min(unlockedLevel - 1, LEVELS.length - 1);
let phase = "ready";
let cowY = 270;
let cowVelocity = 0;
let distance = 0;
let passedCount = 0;
let obstacleModels = [];
let lastFrame = performance.now();
let animationFrame = 0;
let transitionTimer = 0;
let toastTimer = 0;
let mooCounter = 0;

function readProgress() {
  try {
    const saved = Number.parseInt(localStorage.getItem(STORAGE_KEY) || "1", 10);
    return Number.isFinite(saved) ? Math.max(1, Math.min(LEVELS.length, saved)) : 1;
  } catch {
    return 1;
  }
}

function saveProgress() {
  try {
    localStorage.setItem(STORAGE_KEY, String(unlockedLevel));
  } catch {
    // The game remains fully playable if a privacy setting blocks local storage.
  }
}

function pct(value, total) {
  return `${(value / total) * 100}%`;
}

function currentLevel() {
  return LEVELS[selectedLevel];
}

function renderLevelButtons() {
  els.levelGrid.replaceChildren();

  LEVELS.forEach((level, index) => {
    const isUnlocked = index < unlockedLevel;
    const isSelected = index === selectedLevel;
    const button = document.createElement("button");
    button.type = "button";
    button.className = `level-button${isSelected ? " is-selected" : ""}`;
    button.disabled = !isUnlocked;
    button.dataset.level = String(index);
    button.setAttribute(
      "aria-label",
      isUnlocked
        ? `Play level ${index + 1}, ${level.name}, ${level.gapCenters.length} obstacles`
        : `Level ${index + 1}, ${level.name}, locked`,
    );
    button.innerHTML = `
      <span class="level-topline">
        <span class="level-number">${index + 1}</span>
        <span class="lock-icon">${isUnlocked ? (isSelected ? "Selected" : "Open") : "🔒 Locked"}</span>
      </span>
      <h3>${level.name}</h3>
      <p>${level.blurb}</p>
      <span class="level-meta">${level.gapCenters.length} haystacks</span>
    `;
    button.addEventListener("click", () => selectLevel(index));
    els.levelGrid.append(button);
  });
}

function selectLevel(index) {
  if (index >= unlockedLevel || index < 0 || index >= LEVELS.length) return;
  clearTimeout(transitionTimer);
  selectedLevel = index;
  prepareLevel();
  renderLevelButtons();
  els.stage.scrollIntoView({ behavior: "smooth", block: "center" });
  els.stage.focus({ preventScroll: true });
}

function createObstacleElements() {
  els.obstacles.replaceChildren();
  const level = currentLevel();

  obstacleModels = level.gapCenters.map((center, index) => {
    const baseX = 810 + index * level.spacing;
    const gapTop = center - level.gapHeight / 2;
    const gapBottom = center + level.gapHeight / 2;
    const node = document.createElement("div");
    const top = document.createElement("div");
    const bottom = document.createElement("div");

    node.className = "obstacle";
    top.className = "hay-column top";
    bottom.className = "hay-column bottom";
    top.style.height = pct(gapTop, LOGICAL_HEIGHT);
    bottom.style.height = pct(GROUND_Y - gapBottom, LOGICAL_HEIGHT);
    node.append(top, bottom);
    els.obstacles.append(node);

    return {
      baseX,
      gapTop,
      gapBottom,
      node,
      top,
      bottom,
      passed: false,
    };
  });
}

function prepareLevel() {
  const level = currentLevel();
  phase = "ready";
  cowY = level.gapCenters[0] - COW_HEIGHT / 2;
  cowVelocity = 0;
  distance = 0;
  passedCount = 0;
  mooCounter = 0;

  createObstacleElements();
  updateHud();
  renderWorld();

  els.bessie.classList.add("is-waiting");
  els.bessie.classList.remove("is-crashed");
  els.overlay.classList.remove("is-hidden");
  els.overlayEyebrow.textContent = `Level ${selectedLevel + 1} · ${level.name}`;
  els.overlayTitle.textContent = selectedLevel === 0 ? "Ready to flap?" : "Bessie is ready";
  els.overlayCopy.textContent = level.intro;
  els.controlHints.classList.remove("is-hidden");
  els.startButton.classList.remove("is-hidden");
  els.secondaryButton.classList.add("is-hidden");
  els.startButtonLabel.textContent = `Play Level ${selectedLevel + 1}`;
  els.pauseButton.querySelector("span").textContent = "Ⅱ";
  els.pauseButton.setAttribute("aria-label", "Pause game");
  setStatus(`Level ${selectedLevel + 1}, ${level.name}, ready.`);
}

function beginFlight(withFlap = true) {
  if (phase !== "ready" && phase !== "crashed") return;

  if (phase === "crashed") {
    prepareLevel();
  }

  phase = "playing";
  els.overlay.classList.add("is-hidden");
  els.bessie.classList.remove("is-waiting", "is-crashed");
  lastFrame = performance.now();
  setStatus(`Level ${selectedLevel + 1} started.`);

  if (withFlap) flap();
}

function flap() {
  if (phase !== "playing") return;
  cowVelocity = currentLevel().flap;
  mooCounter += 1;
  if (mooCounter % 3 === 1) createMooParticle();
}

function createMooParticle() {
  const particle = document.createElement("span");
  particle.className = "particle";
  particle.textContent = mooCounter % 6 === 1 ? "Moo!" : "flap!";
  particle.style.left = pct(COW_X + 8, LOGICAL_WIDTH);
  particle.style.top = pct(cowY + 8, LOGICAL_HEIGHT);
  els.particles.append(particle);
  window.setTimeout(() => particle.remove(), 700);
}

function pauseGame() {
  if (phase === "playing") {
    phase = "paused";
    els.overlay.classList.remove("is-hidden");
    els.overlayEyebrow.textContent = `Level ${selectedLevel + 1} · Paused`;
    els.overlayTitle.textContent = "Taking a graze break";
    els.overlayCopy.textContent = "Bessie is hovering patiently. Continue whenever you're ready.";
    els.controlHints.classList.add("is-hidden");
    els.startButton.classList.remove("is-hidden");
    els.startButtonLabel.textContent = "Resume flight";
    els.secondaryButton.classList.add("is-hidden");
    els.pauseButton.querySelector("span").textContent = "▶";
    els.pauseButton.setAttribute("aria-label", "Resume game");
    setStatus("Game paused.");
  } else if (phase === "paused") {
    phase = "playing";
    els.overlay.classList.add("is-hidden");
    els.pauseButton.querySelector("span").textContent = "Ⅱ";
    els.pauseButton.setAttribute("aria-label", "Pause game");
    lastFrame = performance.now();
    setStatus("Flight resumed.");
  }
}

function crash() {
  if (phase !== "playing") return;
  phase = "crashed";
  els.bessie.classList.add("is-crashed");
  els.overlay.classList.remove("is-hidden");
  els.overlayEyebrow.textContent = `Level ${selectedLevel + 1} · ${passedCount}/${currentLevel().gapCenters.length} cleared`;
  els.overlayTitle.textContent = "Hay happens!";
  els.overlayCopy.textContent =
    "Bessie is okay—and the course never changes. Use what you learned and try that route again.";
  els.controlHints.classList.add("is-hidden");
  els.startButton.classList.remove("is-hidden");
  els.startButtonLabel.textContent = "Try again";
  els.secondaryButton.classList.remove("is-hidden");
  els.secondaryButton.textContent = "Choose a level";
  setStatus(`Bessie bumped a haystack after clearing ${passedCount} obstacles.`);
}

function completeLevel() {
  if (phase !== "playing") return;
  phase = "transition";

  const levelNumber = selectedLevel + 1;
  const hasNextLevel = selectedLevel < LEVELS.length - 1;
  let newlyUnlocked = false;

  if (hasNextLevel && unlockedLevel < levelNumber + 1) {
    unlockedLevel = levelNumber + 1;
    saveProgress();
    newlyUnlocked = true;
    const message = `🔓 Level ${unlockedLevel} unlocked: ${LEVELS[unlockedLevel - 1].name}!`;
    showToast(message);
    console.log(message);
  }

  renderLevelButtons();
  els.overlay.classList.remove("is-hidden");
  els.controlHints.classList.add("is-hidden");
  els.secondaryButton.classList.add("is-hidden");
  els.startButton.classList.add("is-hidden");

  if (hasNextLevel) {
    els.overlayEyebrow.textContent = `Level ${levelNumber} complete`;
    els.overlayTitle.textContent = "Udderly brilliant!";
    els.overlayCopy.textContent = newlyUnlocked
      ? `Level ${levelNumber + 1} is unlocked. Moving Bessie to ${LEVELS[selectedLevel + 1].name}…`
      : `Great flight. Moving Bessie to ${LEVELS[selectedLevel + 1].name}…`;
    setStatus(`Level ${levelNumber} complete. Level ${levelNumber + 1} unlocked.`);

    transitionTimer = window.setTimeout(() => {
      selectedLevel += 1;
      prepareLevel();
      renderLevelButtons();
    }, 1900);
  } else {
    els.overlayEyebrow.textContent = "All 5 levels complete";
    els.overlayTitle.textContent = "Barnyard legend!";
    els.overlayCopy.textContent =
      "Bessie conquered every fixed haystack course. Pick any unlocked pasture and fly it again.";
    els.startButton.classList.remove("is-hidden");
    els.startButtonLabel.textContent = "Fly Level 5 again";
    els.secondaryButton.classList.remove("is-hidden");
    els.secondaryButton.textContent = "Choose a level";
    phase = "complete";
    setStatus("Congratulations! Bessie completed all five levels.");
  }
}

function showToast(message) {
  clearTimeout(toastTimer);
  els.toast.textContent = message;
  els.toast.classList.add("is-visible");
  toastTimer = window.setTimeout(() => els.toast.classList.remove("is-visible"), 4300);
}

function setStatus(message) {
  els.gameStatus.textContent = message;
}

function updateHud() {
  const total = currentLevel().gapCenters.length;
  els.hudLevel.textContent = String(selectedLevel + 1);
  els.hudScore.textContent = `${passedCount} / ${total}`;
  els.progressFill.style.width = `${(passedCount / total) * 100}%`;
}

function update(dt) {
  const level = currentLevel();
  cowVelocity += level.gravity * dt;
  cowY += cowVelocity * dt;
  distance += level.speed * dt;

  // Render before measuring so collision stays accurate at every responsive aspect ratio.
  renderWorld();

  const stageRect = els.stage.getBoundingClientRect();
  const cowRect = els.bessie.getBoundingClientRect();
  const cowInsetX = cowRect.width * 0.16;
  const cowInsetTop = cowRect.height * 0.18;
  const cowInsetBottom = cowRect.height * 0.14;
  const cowLeft = cowRect.left + cowInsetX;
  const cowRight = cowRect.right - cowInsetX;
  const cowTop = cowRect.top + cowInsetTop;
  const cowBottom = cowRect.bottom - cowInsetBottom;
  const groundTop = stageRect.top + (GROUND_Y / LOGICAL_HEIGHT) * stageRect.height;

  if (cowTop <= stageRect.top || cowBottom >= groundTop) {
    crash();
    return;
  }

  for (const obstacle of obstacleModels) {
    const obstacleRect = obstacle.node.getBoundingClientRect();
    const topRect = obstacle.top.getBoundingClientRect();
    const bottomRect = obstacle.bottom.getBoundingClientRect();
    const obstacleInset = obstacleRect.width * 0.08;
    const obstacleLeft = obstacleRect.left + obstacleInset;
    const obstacleRight = obstacleRect.right - obstacleInset;
    const overlapsX = cowRight > obstacleLeft && cowLeft < obstacleRight;
    const hitsHay = cowTop < topRect.bottom || cowBottom > bottomRect.top;

    if (overlapsX && hitsHay) {
      crash();
      return;
    }

    if (!obstacle.passed && obstacleRight < cowLeft) {
      obstacle.passed = true;
      passedCount += 1;
      updateHud();
      setStatus(`${passedCount} of ${level.gapCenters.length} haystacks cleared.`);
    }
  }

  const lastObstacle = obstacleModels[obstacleModels.length - 1];
  if (
    lastObstacle?.passed &&
    lastObstacle.node.getBoundingClientRect().right < stageRect.left - 30
  ) {
    completeLevel();
  }
}

function renderWorld() {
  els.bessie.style.top = pct(cowY, LOGICAL_HEIGHT);
  const rotation = Math.max(-22, Math.min(62, cowVelocity * 0.075));
  els.bessie.style.transform = `rotate(${phase === "ready" ? -4 : rotation}deg)`;

  obstacleModels.forEach((obstacle) => {
    obstacle.node.style.left = pct(obstacle.baseX - distance, LOGICAL_WIDTH);
  });

  const lastObstacle = obstacleModels[obstacleModels.length - 1];
  if (lastObstacle) {
    const finishX = lastObstacle.baseX + 245 - distance;
    els.finishLine.style.display =
      finishX > -80 && finishX < LOGICAL_WIDTH + 80 ? "block" : "none";
    els.finishLine.style.left = pct(finishX, LOGICAL_WIDTH);
  }

  const parallax = -distance * 0.018;
  els.stage.style.backgroundPositionX = `calc(50% + ${parallax}px)`;
}

function gameLoop(now) {
  const dt = Math.min((now - lastFrame) / 1000, 0.033);
  lastFrame = now;

  if (phase === "playing") {
    update(dt);
  }

  animationFrame = requestAnimationFrame(gameLoop);
}

function primaryAction() {
  if (phase === "paused") {
    pauseGame();
  } else if (phase === "complete") {
    prepareLevel();
    beginFlight(true);
  } else {
    beginFlight(true);
  }
}

els.startButton.addEventListener("click", (event) => {
  event.stopPropagation();
  primaryAction();
});

els.secondaryButton.addEventListener("click", (event) => {
  event.stopPropagation();
  clearTimeout(transitionTimer);
  if (phase === "complete") prepareLevel();
  document.querySelector(".below-stage").scrollIntoView({ behavior: "smooth", block: "start" });
});

els.pauseButton.addEventListener("click", (event) => {
  event.stopPropagation();
  pauseGame();
});

els.stage.addEventListener("pointerdown", (event) => {
  if (event.target.closest("button")) return;
  if (phase === "playing") {
    flap();
  } else if (phase === "ready" || phase === "crashed") {
    beginFlight(true);
  }
});

window.addEventListener("keydown", (event) => {
  if (event.code === "Space" || event.code === "ArrowUp") {
    event.preventDefault();
    if (phase === "playing") {
      flap();
    } else if (phase === "ready" || phase === "crashed") {
      beginFlight(true);
    }
  } else if (event.code === "KeyP" || event.code === "Escape") {
    event.preventDefault();
    pauseGame();
  }
});

els.resetProgressButton.addEventListener("click", () => {
  clearTimeout(transitionTimer);
  unlockedLevel = 1;
  selectedLevel = 0;
  saveProgress();
  prepareLevel();
  renderLevelButtons();
  showToast("Progress reset. Level 1 is ready for Bessie!");
  console.log("Progress reset. Only Level 1 is unlocked.");
});

document.addEventListener("visibilitychange", () => {
  if (document.hidden && phase === "playing") pauseGame();
});

renderLevelButtons();
prepareLevel();
cancelAnimationFrame(animationFrame);
animationFrame = requestAnimationFrame(gameLoop);
