/* Plain JavaScript MPQS/SIQS, using the public-domain msieve implementations
 * by Jason Papadopoulos as a guide. Integer arithmetic uses BigInt; bounded
 * sieve residues and packed GF(2) vectors use exact 32-bit/Number arithmetic. */
(function (root) {
  "use strict";
  const createMatrixAlgebra = typeof module !== "undefined" && module.exports ?
    require("./matrix.js").createMatrixAlgebra : root.MPQSMatrix.createMatrixAlgebra;
  const createECM = typeof module !== "undefined" && module.exports ?
    require("./ecm.js").createECM : root.MPQSECM.createECM;

  // Keep the engine self-contained so the same source runs in a Blob Worker
  // from either file:// or a web server, and in Node for regression tests.
  function createEngine() {
    "use strict";
    const MAX_DIGITS = 110;
    const MR_ROUNDS = 30;
    const RHO_MIN_DIGITS = 40, RHO_MAX_MS = 30000, RHO_MAX_EVALUATIONS = 20000000;
    const SQUFOF_LIMIT = 10n ** 24n; // At most 24 decimal digits.
    const matrixAlgebra = createMatrixAlgebra();
    const ecm = createECM();
    const abs = n => n < 0n ? -n : n;
    const mod = (n, m) => ((n % m) + m) % m;

    function validateInput(value) {
      const text = String(value).trim();
      if (!text) return { valid: false, message: "Enter a composite number." };
      if (!/^[0-9]+$/.test(text)) {
        return { valid: false, message: "Use decimal digits only (0–9)." };
      }
      if (text.length > MAX_DIGITS) {
        return { valid: false, message: "Enter no more than 110 digits." };
      }
      return { valid: true, text };
    }

    function gcd(a, b) {
      a = abs(a);
      b = abs(b);
      while (b !== 0n) [a, b] = [b, a % b];
      return a;
    }

    function powMod(a, e, n) {
      let result = 1n;
      a = mod(a, n);
      while (e > 0n) {
        if (e & 1n) result = result * a % n;
        a = a * a % n;
        e >>= 1n;
      }
      return result;
    }

    function inverse(a, n) {
      let b = n, x = 1n, y = 0n;
      a = mod(a, n);
      while (b !== 0n) {
        const q = a / b;
        [a, b] = [b, a - q * b];
        [x, y] = [y, x - q * y];
      }
      if (a !== 1n) throw new Error("A modular inverse does not exist.");
      return mod(x, n);
    }

    function isqrt(n) {
      if (n < 0n) throw new RangeError("Negative square root.");
      if (n < 2n) return n;
      let x = 1n << BigInt(Math.ceil(n.toString(2).length / 2));
      for (;;) {
        const next = (x + n / x) >> 1n;
        if (next >= x) return x;
        x = next;
      }
    }

    function integerRoot(n, k) {
      if (k === 2) return isqrt(n);
      const exponent = BigInt(k);
      let x = 1n << BigInt(Math.ceil(n.toString(2).length / k));
      for (;;) {
        const next = ((exponent - 1n) * x + n / (x ** (exponent - 1n))) / exponent;
        if (next >= x) return x;
        x = next;
      }
    }

    // Rejection sampling, rather than reduction modulo the range, gives each
    // possible base the same probability. All 30 rounds use bases in [2,n-2].
    function randomBelow(limit) {
      if (limit <= 0n) throw new RangeError("Empty random range.");
      if (limit === 1n) return 0n;
      const bits = (limit - 1n).toString(2).length;
      const words = new Uint32Array(Math.ceil(bits / 32));
      const highBits = bits - (words.length - 1) * 32;
      for (;;) {
        globalThis.crypto.getRandomValues(words);
        if (highBits < 32) words[0] &= (2 ** highBits) - 1;
        let value = 0n;
        for (const word of words) value = (value << 32n) | BigInt(word);
        if (value < limit) return value;
      }
    }

    function primesThrough(limit) {
      const composite = new Uint8Array(limit + 1);
      const primes = [];
      for (let p = 2; p <= limit; p++) {
        if (composite[p]) continue;
        primes.push(p);
        for (let j = p * p; j <= limit; j += p) composite[j] = 1;
      }
      return primes;
    }

    const trialPrimes = primesThrough(1000);

    function isProbablePrime(n) {
      if (n < 2n) return false;
      for (const p of trialPrimes) {
        const prime = BigInt(p);
        if (n === prime) return true;
        if (n % prime === 0n) return false;
      }
      let d = n - 1n, s = 0;
      while ((d & 1n) === 0n) { d >>= 1n; s++; }
      for (let round = 0; round < MR_ROUNDS; round++) {
        const a = 2n + randomBelow(n - 3n);
        let x = powMod(a, d, n);
        if (x === 1n || x === n - 1n) continue;
        let passed = false;
        for (let r = 1; r < s; r++) {
          x = x * x % n;
          if (x === n - 1n) { passed = true; break; }
          if (x === 1n) return false;
        }
        if (!passed) return false;
      }
      return true;
    }

    function brentRhoBudget(n, threads = 1) {
      if (!Number.isInteger(threads) || threads < 1 || threads > 32) throw new Error("Invalid Brent-Rho thread count.");
      const digits = n.toString().length;
      if (digits < RHO_MIN_DIGITS) return { digits, estimatedTotalMs: 0, budgetMs: 0 };
      // A rough size model, not a prediction of a particular input's factors.
      // Historical local records: ~4 s at 52 digits on one worker and ~838 s
      // at 86 digits on eight. Convert the latter to a worker-time proxy,
      // interpolate logarithmically, and divide by the selected worker count.
      const estimatedTotalMs = 4000 * Math.exp(Math.log(6700 / 4) * (digits - 52) / 34) / threads;
      return { digits, estimatedTotalMs, budgetMs: Math.min(RHO_MAX_MS, estimatedTotalMs * 0.01) };
    }

    // BrentRho/BrentRho.c: x <- x^2+c, doubling cycle lengths, batches of
    // 128 differences, and replay if the batched GCD equals n. The controller
    // bounds every walk, including the replay; n itself is never a divisor.
    function* brentAttemptSteps(n, seed, c, controller) {
      let y = seed % n, r = 1;
      const step = value => {
        controller.evaluations++;
        return (value * value + c) % n;
      };
      let nextReport = -Infinity;
      const expired = () => {
        if (controller.evaluations >= controller.limit) return true;
        const current = controller.now();
        if (controller.onProgress && current >= nextReport) {
          controller.onProgress();
          nextReport = current + 200;
        }
        return current >= controller.deadline;
      };
      for (;;) {
        if (controller.shouldYield?.()) yield controller;
        if (expired()) return null;
        const x = y;
        for (let i = 0; i < r;) {
          if (controller.shouldYield?.()) yield controller;
          if (expired()) return null;
          const count = Math.min(128, r - i, controller.limit - controller.evaluations);
          for (let j = 0; j < count; j++) y = step(y);
          i += count;
        }
        for (let k = 0; k < r;) {
          if (controller.shouldYield?.()) yield controller;
          if (expired()) return null;
          let saved = y, product = 1n;
          const count = Math.min(128, r - k, controller.limit - controller.evaluations);
          for (let i = 0; i < count; i++) {
            y = step(y);
            product = product * abs(x - y) % n;
          }
          let divisor = gcd(product, n);
          controller.gcds++;
          k += count;
          if (divisor === 1n) continue;
          if (divisor === n) {
            controller.replays++;
            for (let i = 0; i < count; i++) {
              if (controller.shouldYield?.()) yield controller;
              if (expired()) return null;
              saved = step(saved);
              divisor = gcd(abs(x - saved), n);
              controller.gcds++;
              if (divisor > 1n) break;
            }
          }
          return divisor > 1n && divisor < n ? divisor : null;
        }
        r *= 2;
      }
    }

    function brentAttempt(n, seed, c, controller) {
      const steps = brentAttemptSteps(n, seed, c, controller);
      let step = steps.next();
      while (!step.done) step = steps.next();
      return step.value;
    }

    function* brentRhoSteps(n, budgetMs, options = {}) {
      const now = options.now || (() => performance.now());
      const maxEvaluations = Math.min(RHO_MAX_EVALUATIONS, options.maxEvaluations ?? RHO_MAX_EVALUATIONS);
      const maxAttempts = options.attemptIds?.length ?? options.maxAttempts ?? 8;
      if (!Number.isSafeInteger(maxAttempts) || maxAttempts < 1 || maxAttempts > 32) {
        throw new Error("Invalid Brent-Rho budget.");
      }
      const attemptIds = options.attemptIds ?? Array.from({ length: maxAttempts }, (_, i) => i);
      if (!Number.isFinite(budgetMs) || budgetMs < 0 ||
          !Number.isSafeInteger(maxEvaluations) || maxEvaluations < 0 ||
          !Number.isSafeInteger(maxAttempts) || maxAttempts < 1 || maxAttempts > 32 ||
          new Set(attemptIds).size !== maxAttempts ||
          attemptIds.some(id => !Number.isInteger(id) || id < 0 || id >= 32)) {
        throw new Error("Invalid Brent-Rho budget.");
      }
      const result = { divisor: null, attempts: 0, evaluations: 0, gcds: 0, replays: 0, elapsedMs: 0, threads: 1 };
      if (n.toString().length < RHO_MIN_DIGITS || budgetMs === 0 || maxEvaluations === 0) return result;
      const start = now(), deadline = start + Math.min(RHO_MAX_MS, budgetMs);
      let lastReport = start;
      for (let attempt = 0; attempt < maxAttempts; attempt++) {
        const current = now();
        if (current >= deadline || result.evaluations >= maxEvaluations) break;
        const remainingAttempts = maxAttempts - attempt;
        // Each new polynomial gets a share of the remaining allowance.
        const controller = { now, evaluations: result.evaluations, gcds: 0, replays: 0,
          limit: result.evaluations + Math.ceil((maxEvaluations - result.evaluations) / remainingAttempts),
          deadline: current + (deadline - current) / remainingAttempts,
          shouldYield: options.shouldYield, baseGcds: result.gcds, baseReplays: result.replays };
        result.attempts++;
        controller.attempts = result.attempts;
        controller.onProgress = () => options.onProgress?.({ ...result,
          evaluations: controller.evaluations, elapsedMs: Math.max(0, now() - start) });
        const id = attemptIds[attempt];
        const divisor = yield* brentAttemptSteps(n, BigInt(id + 2), BigInt(id + 1), controller);
        result.evaluations = controller.evaluations;
        result.gcds += controller.gcds;
        result.replays += controller.replays;
        result.elapsedMs = Math.max(0, now() - start);
        if (divisor !== null) {
          if (n % divisor !== 0n) throw new Error("Brent-Rho produced an invalid divisor.");
          result.divisor = divisor;
          break;
        }
        if (now() - lastReport >= 200) {
          options.onProgress?.(result);
          lastReport = now();
        }
      }
      result.elapsedMs = Math.max(0, now() - start);
      return result;
    }

    function brentRho(n, budgetMs, options = {}) {
      const steps = brentRhoSteps(n, budgetMs, options);
      let step = steps.next();
      while (!step.done) step = steps.next();
      return step.value;
    }

    // Preserve a walk between message-sized batches. Only time spent inside
    // advance() consumes the quota; waiting for other workers is not charged.
    function createRhoJob(n, budgetMs, options = {}) {
      const clock = options.clock || (() => performance.now());
      let spent = 0, started = 0, active = false, checks = 0, finished = null;
      const now = () => spent + (active ? Math.max(0, clock() - started) : 0);
      const steps = brentRhoSteps(n, budgetMs, { ...options, now,
        shouldYield: () => ++checks > 64 || (checks > 1 && clock() - started >= 4) });
      return {
        advance() {
          if (finished) return finished;
          checks = 0;
          started = clock(); active = true;
          let step;
          try { step = steps.next(); }
          finally { spent += Math.max(0, clock() - started); active = false; }
          const state = step.done ? step.value : {
            divisor: null, attempts: step.value.attempts, evaluations: step.value.evaluations,
            gcds: step.value.baseGcds + step.value.gcds,
            replays: step.value.baseReplays + step.value.replays
          };
          const result = { ...state, elapsedMs: spent, done: step.done };
          if (step.done) finished = result;
          return result;
        }
      };
    }

    let nextRhoJobId = 0;
    async function brentRhoParallel(n, budgetMs, ports, options = {}) {
      if (!ports.length) return brentRho(n, budgetMs, options);
      if (!Number.isFinite(budgetMs) || budgetMs < 0 || ports.length > 31) {
        throw new Error("Invalid Brent-Rho budget.");
      }
      const threads = ports.length + 1;
      budgetMs = Math.min(RHO_MAX_MS, budgetMs);
      const maxEvaluations = Math.min(RHO_MAX_EVALUATIONS, options.maxEvaluations ?? RHO_MAX_EVALUATIONS);
      if (!Number.isSafeInteger(maxEvaluations) || maxEvaluations < 0) throw new Error("Invalid Brent-Rho budget.");
      const attemptCount = Math.max(8, threads);
      const workerBudgetsMs = Array.from({ length: threads }, (_, i) =>
        budgetMs * (i + 1) / threads - budgetMs * i / threads);
      const workerEvaluationLimits = Array.from({ length: threads }, (_, i) =>
        Math.floor(maxEvaluations / threads) + (i < maxEvaluations % threads ? 1 : 0));
      const workerAttemptIds = Array.from({ length: threads }, (_, i) =>
        Array.from({ length: attemptCount }, (_, id) => id).filter(id => id % threads === i));
      const states = Array.from({ length: threads }, () => ({
        divisor: null, attempts: 0, evaluations: 0, gcds: 0, replays: 0, elapsedMs: 0, done: false
      }));
      const started = performance.now();
      function snapshot() {
        return { divisor: states.find(s => s.divisor !== null)?.divisor ?? null,
          attempts: states.reduce((sum, s) => sum + s.attempts, 0),
          evaluations: states.reduce((sum, s) => sum + s.evaluations, 0),
          gcds: states.reduce((sum, s) => sum + s.gcds, 0),
          replays: states.reduce((sum, s) => sum + s.replays, 0),
          elapsedMs: states.reduce((sum, s) => sum + s.elapsedMs, 0),
          wallMs: performance.now() - started, threads,
          workerBudgetsMs, workerEvaluationLimits, workerAttemptIds,
          workerEvaluations: states.map(s => s.evaluations),
          workerElapsedMs: states.map(s => s.elapsedMs) };
      }
      if (n.toString().length < RHO_MIN_DIGITS || budgetMs === 0 || maxEvaluations === 0) return snapshot();
      const jobId = ++nextRhoJobId, pending = new Map();
      let nextRequestId = 0;
      const previousHandlers = ports.map(port => port.onmessage);
      function request(index, operation, extra = {}) {
        return new Promise((resolve, reject) => {
          const requestId = ++nextRequestId;
          pending.set(requestId, { index, resolve, reject });
          try { ports[index].postMessage({ type: "rho", jobId, requestId, operation, ...extra }); }
          catch (error) { pending.delete(requestId); reject(error); }
        });
      }
      function accept(index, state) {
        const previous = states[index];
        if (!state || !Number.isSafeInteger(state.evaluations) ||
            state.evaluations < previous.evaluations || state.evaluations > workerEvaluationLimits[index] ||
            !Number.isSafeInteger(state.attempts) || state.attempts < previous.attempts ||
            state.attempts > workerAttemptIds[index].length ||
            !Number.isSafeInteger(state.gcds) || state.gcds < previous.gcds ||
            !Number.isSafeInteger(state.replays) || state.replays < previous.replays ||
            !Number.isFinite(state.elapsedMs) || state.elapsedMs < previous.elapsedMs ||
            typeof state.done !== "boolean" ||
            (state.divisor !== null && (typeof state.divisor !== "bigint" ||
              state.divisor <= 1n || state.divisor >= n || n % state.divisor !== 0n))) {
          throw new Error("Invalid Brent-Rho worker result.");
        }
        states[index] = state;
      }
      try {
        ports.forEach((port, index) => {
          port.onmessage = ({ data }) => {
            if (data.jobId !== jobId || (data.type !== "rho-result" && data.type !== "rho-error")) return;
            const call = pending.get(data.requestId);
            if (!call || call.index !== index) return;
            pending.delete(data.requestId);
            if (data.type === "rho-error") call.reject(new Error(data.message));
            else call.resolve(data.result);
          };
        });
        const local = createRhoJob(n, workerBudgetsMs[0], {
          maxEvaluations: workerEvaluationLimits[0], attemptIds: workerAttemptIds[0]
        });
        await Promise.all(ports.map((_, i) => request(i, "init", { n, budgetMs: workerBudgetsMs[i + 1],
          options: { maxEvaluations: workerEvaluationLimits[i + 1], attemptIds: workerAttemptIds[i + 1] } })));
        let lastReport = -Infinity;
        for (;;) {
          // One short batch per worker in flight. A successful batch stops the
          // next wave; existing batches finish before helpers are released.
          const replies = ports.map((_, i) => states[i + 1].done ? Promise.resolve() :
            request(i, "step").then(state => accept(i + 1, state)));
          replies.push(Promise.resolve().then(() => {
            if (!states[0].done) accept(0, local.advance());
          }));
          await Promise.all(replies);
          const result = snapshot();
          if (performance.now() - lastReport >= 200 || result.divisor !== null || states.every(s => s.done)) {
            options.onProgress?.(result);
            lastReport = performance.now();
          }
          if (result.divisor !== null || states.every(s => s.done)) return result;
        }
      } finally {
        await Promise.allSettled(ports.map((_, i) => request(i, "release")));
        ports.forEach((port, i) => { port.onmessage = previousHandlers[i]; });
      }
    }

    // Tonelli–Shanks for an odd prime p and a quadratic residue a.
    function sqrtMod(a, p) {
      a = mod(a, p);
      if (a === 0n || p === 2n) return a;
      if (powMod(a, (p - 1n) / 2n, p) !== 1n) return null;
      if (p % 4n === 3n) return powMod(a, (p + 1n) / 4n, p);
      let d = p - 1n, s = 0;
      while ((d & 1n) === 0n) { d >>= 1n; s++; }
      let z = 2n;
      while (powMod(z, (p - 1n) / 2n, p) !== p - 1n) z++;
      let c = powMod(z, d, p), x = powMod(a, (d + 1n) / 2n, p);
      let t = powMod(a, d, p), m = s;
      while (t !== 1n) {
        let v = t, i = 0;
        do { v = v * v % p; i++; } while (v !== 1n && i < m);
        if (i === m) throw new Error("Invalid modular square root.");
        const b = powMod(c, 1n << BigInt(m - i - 1), p);
        x = x * b % p;
        c = b * b % p;
        t = t * c % p;
        m = i;
      }
      return x;
    }

    function log2(n) {
      const shift = Math.max(0, n.toString(2).length - 53);
      return Math.log2(Number(n >> BigInt(shift))) + shift;
    }

    function parameters(digits) {
      const bound = digits <= 14 ? 300 : digits <= 18 ? 1000 :
        digits <= 24 ? 3000 : digits <= 30 ? 10000 :
        digits <= 36 ? 25000 : digits <= 45 ? 60000 : 120000;
      return { bound, halfInterval: digits <= 18 ? 8192 : digits <= 30 ? 16384 : 32768 };
    }

    function mergeExponents(a, b) {
      const counts = new Map(a);
      for (const [column, count] of b) counts.set(column, (counts.get(column) || 0) + count);
      return Array.from(counts);
    }

    function emptyRelationCounts() {
      // fullRelations retains its original API meaning: all usable matrix
      // rows. smoothRelations separates the directly smooth subset.
      return { fullRelations: 0, smoothRelations: 0, combinedRelations: 0,
        partialRelations: 0, partialsFound: 0, estimatedRelations: 0, maxLargePrimes: 1 };
    }

    function beginSieve(n, report, stats) {
      const counts = emptyRelationCounts();
      Object.assign(stats, counts);
      report({ stage: "Building the factor base", digits: n.toString().length, ...counts });
    }

    function tinyMpqs(n, report, stats) {
      const { bound, halfInterval: M } = parameters(n.toString().length);
      const length = 2 * M + 1;
      const fb = [];
      beginSieve(n, report, stats);
      for (const p of primesThrough(bound)) {
        const big = BigInt(p);
        if (n % big === 0n) return big;
        const squareRoot = sqrtMod(n, big);
        if (squareRoot !== null) fb.push({ p, big, root: squareRoot, weight: Math.log2(p) });
      }
      const largePrimeBound = BigInt(bound) ** 2n;
      const partials = new Map();
      const maxPartials = 50000;
      let rows = [], pivots = [];
      let q = isqrt(isqrt(2n * n) / BigInt(M));
      if (q < 3n) q = 3n;
      if (q % 2n === 0n) q++;
      let localPolynomials = 0, localRelations = 0;
      let smoothRelations = 0, combinedRelations = 0, partialsFound = 0;

      function progress() {
        const counts = { fullRelations: rows.length, smoothRelations, combinedRelations,
          partialRelations: partials.size, partialsFound, estimatedRelations: fb.length + 1 + 16,
          maxLargePrimes: 1 };
        Object.assign(stats, counts);
        report({ stage: "Sieving polynomials", digits: n.toString().length,
          polynomials: localPolynomials, relations: localRelations, ...counts });
      }
      progress();

      function extract(dependency) {
        let X = 1n, Y = 1n;
        const powers = new Uint32Array(fb.length + 1);
        let index = 0;
        while (dependency !== 0n) {
          if (dependency & 1n) {
            const row = rows[index];
            X = X * row.x % n;
            Y = Y * row.square % n;
            for (const [column, count] of row.exponents) powers[column] += count;
          }
          dependency >>= 1n;
          index++;
        }
        for (let column = 0; column < powers.length; column++) {
          if (powers[column] & 1) throw new Error("Invalid MPQS dependency.");
          if (column > 0 && powers[column]) {
            Y = Y * powMod(fb[column - 1].big, BigInt(powers[column] / 2), n) % n;
          }
        }
        // Always take GCDs with the actual number being factored and reject
        // trivial dependencies. A failed dependency is not a factoring failure.
        for (const difference of [X - Y, X + Y]) {
          const divisor = gcd(difference, n);
          if (divisor > 1n && divisor < n) return divisor;
        }
        return null;
      }

      function addRelation(row) {
        const index = rows.length;
        rows.push(row);
        localRelations++;
        stats.relations++;
        let parity = row.parity, dependency = 1n << BigInt(index);
        while (parity !== 0n) {
          const column = parity.toString(2).length - 1;
          if (!pivots[column]) {
            pivots[column] = { parity, dependency };
            return null;
          }
          parity ^= pivots[column].parity;
          dependency ^= pivots[column].dependency;
        }
        return extract(dependency);
      }

      for (;;) {
        let root;
        for (;;) {
          if (n % q === 0n && q < n) return q;
          if (isProbablePrime(q)) {
            root = sqrtMod(n, q);
            if (root !== null) break;
          }
          q += 2n;
        }
        const A = q * q;
        // Hensel lifting: B² = n (mod q²), hence C is an exact integer.
        let B = root + q * mod(((n - root * root) / q) * inverse(2n * root, q), q);
        if (B > A / 2n) B = A - B;
        const C = (B * B - n) / A;
        const evaluate = x => A * x * x + 2n * B * x + C;
        const scores = new Float32Array(length);
        const polynomialRoots = [];

        function mark(residue, modulus, weight) {
          // These values are sieve indices, bounded by the interval size.
          const start = ((residue + M) % modulus + modulus) % modulus;
          for (let index = start; index < length; index += modulus) scores[index] += weight;
        }

        for (const prime of fb) {
          const { p, big, weight } = prime;
          let roots;
          if (p === 2) {
            roots = [Number(mod(-C, 2n))];
          } else if (A % big === 0n) {
            roots = [Number(mod(-C * inverse(2n * B, big), big))];
          } else {
            const invA = inverse(A, big);
            roots = [Number(mod((prime.root - B) * invA, big)),
              Number(mod((-prime.root - B) * invA, big))];
          }
          polynomialRoots.push(roots);
          for (const r of roots) mark(r, p, weight);

          // Include prime powers in the logarithmic sieve. BigInt is used for
          // polynomial arithmetic; only small residues/indices become Number.
          let modulus = p;
          while (modulus * p <= length && roots.length) {
            const nextModulus = modulus * p;
            const lifted = [];
            for (const r of roots) {
              if (p === 2) {
                for (const candidate of [r, r + modulus]) {
                  if (evaluate(BigInt(candidate)) % BigInt(nextModulus) === 0n) lifted.push(candidate);
                }
              } else {
                const quotient = evaluate(BigInt(r)) / BigInt(modulus);
                const derivative = 2n * (A * BigInt(r) + B);
                const step = mod(-quotient * inverse(derivative, big), big);
                lifted.push(r + Number(step) * modulus);
              }
            }
            roots = lifted;
            modulus = nextModulus;
            for (const r of roots) mark(r, modulus, weight);
          }
        }

        localPolynomials++;
        stats.polynomials++;
        const threshold = Math.max(0, log2(abs(C) || 1n) - log2(largePrimeBound) - 6);
        for (let index = 0; index < length; index++) {
          if (scores[index] < threshold) continue;
          const offset = index - M, x = BigInt(offset);
          let remainder = evaluate(x);
          if (remainder === 0n) continue; // A nonsquare n cannot produce Q(x)=0.
          const exponents = [];
          let parity = 0n;
          if (remainder < 0n) {
            remainder = -remainder;
            exponents.push([0, 1]);
            parity = 1n;
          }
          for (let j = 0; j < fb.length && remainder !== 1n; j++) {
            const { p, big } = fb[j];
            const residue = ((offset % p) + p) % p;
            if (!polynomialRoots[j].includes(residue)) continue;
            let count = 0;
            while (remainder % big === 0n) { remainder /= big; count++; }
            if (count) {
              exponents.push([j + 1, count]);
              if (count & 1) parity |= 1n << BigInt(j + 1);
            }
          }
          if (remainder > largePrimeBound) continue;
          let row = { x: mod(A * x + B, n), square: q % n, exponents, parity };
          if (remainder !== 1n) {
            partialsFound++;
            const previous = partials.get(remainder);
            if (!previous) {
              // Bound memory even on inputs that need a very long search.
              if (partials.size >= maxPartials) partials.delete(partials.keys().next().value);
              partials.set(remainder, row);
              continue;
            }
            partials.delete(remainder);
            combinedRelations++;
            row = {
              x: previous.x * row.x % n,
              square: previous.square * row.square % n * remainder % n,
              exponents: mergeExponents(previous.exponents, row.exponents),
              parity: previous.parity ^ row.parity
            };
          } else smoothRelations++;
          const divisor = addRelation(row);
          if (divisor !== null) { progress(); return divisor; }
          if (rows.length >= fb.length + 256) {
            rows = [];
            pivots = [];
            smoothRelations = combinedRelations = 0;
          }
        }
        progress();
        q += 2n;
      }
    }

    // All moduli in these helpers are below 2^26, so products of two residues
    // are exact Numbers (<2^52). The cofactor and polynomial remain BigInts.
    function modSmall(a, p) {
      a %= p;
      return a < 0 ? a + p : a;
    }

    function inverseSmall(a, p) {
      let b = p, x = 1, y = 0;
      while (b) {
        const q = Math.floor(a / b);
        [a, b] = [b, a - q * b];
        [x, y] = [y, x - q * y];
      }
      if (a !== 1) throw new Error("Invalid small modular inverse.");
      return modSmall(x, p);
    }

    function siqsParameters(n) {
      // Interpolate the factor-base sizes from msieve's parameter table.
      // Bound the browser matrix/cache footprint for very large inputs.
      const table = [[64, 100], [128, 450], [183, 2000], [200, 3000],
        [212, 5400], [233, 10000], [249, 27000], [266, 50000],
        [283, 55000], [298, 60000], [315, 80000], [332, 100000],
        [348, 140000], [363, 210000], [379, 300000]];
      const bits = n.toString(2).length;
      let size = 100000;
      for (let i = 1; i < table.length; i++) {
        if (bits <= table[i][0]) {
          const [lo, a] = table[i - 1], [hi, b] = table[i];
          size = Math.round(a + (b - a) * (bits - lo) / (hi - lo));
          break;
        }
      }
      return { size: Math.max(100, Math.min(bits < 240 ? 12000 : 100000, size)),
        halfInterval: bits < 140 ? 16384 : bits < 190 ? 32768 : bits < 240 ? 65536 : 262144,
        largeMultiplier: bits < 190 ? 40 : 80, doubleLargePrimes: bits >= 282,
        largeSieve: bits >= 240 };
    }

    function chooseMultiplier(n) {
      // Modified Knuth–Schroeppel score, following msieve/mpqs/mpqs.c.
      const multipliers = [1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19,
        21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38, 39, 41, 42, 43,
        46, 47, 51, 53, 55, 57, 58, 59, 61, 62, 65, 66, 67, 69, 70, 71, 73];
      const scores = multipliers.map(k => {
        const residue = Number(n % 8n) * k % 8;
        return Math.log(k) / 2 - Math.LN2 * (residue === 1 ? 2 : residue === 5 ? 1 : residue & 1 ? 0.5 : 0);
      });
      for (const p of primesThrough(2000).slice(1, 300)) {
        const residues = new Uint8Array(p);
        for (let i = 1; i < p; i++) residues[i * i % p] = 1;
        const r = Number(n % BigInt(p)), contribution = Math.log(p) / (p - 1);
        for (let j = 0; j < multipliers.length; j++) {
          const a = r * multipliers[j] % p;
          scores[j] -= contribution * (a === 0 ? 1 : residues[a] ? 2 : 0);
        }
      }
      return multipliers[scores.indexOf(Math.min(...scores))];
    }

    function gcdSmall(a, b) {
      while (b) { const r = a % b; a = b; b = r; }
      return a;
    }

    // BigInt port of msieve-1.36/common/smallfact/squfof.c. The C routine
    // restricts k*n to 62 bits; BigInt keeps all scaled forms exact, including
    // inputs up to 24 digits that exceed Number's exact integer range.
    // Each job races disjoint multipliers, retaining forms between messages.
    function createSqufofJob(n, options = {}) {
      if (typeof n !== "bigint" || n < 4n || n >= SQUFOF_LIMIT) {
        throw new Error("SQUFOF requires a composite with at most 24 digits.");
      }
      const workerIndex = options.workerIndex ?? 0, workerCount = options.workerCount ?? 1;
      if (!Number.isInteger(workerCount) || workerCount < 1 || workerCount > 32 ||
          !Number.isInteger(workerIndex) || workerIndex < 0 || workerIndex >= workerCount) {
        throw new Error("Invalid SQUFOF worker assignment.");
      }
      const clock = options.clock || (() => performance.now());
      const result = { divisor: null, iterations: 0, forwardIterations: 0, reverseIterations: 0,
        multipliers: 0, failedMultipliers: 0, elapsedMs: 0, done: false };
      const base = [1155n, 385n, 231n, 165n, 105n, 77n, 55n, 35n, 33n, 21n, 15n, 11n, 7n, 5n, 3n, 1n];
      function* multiplierSequence() {
        yield* base;
        // Additional square-free multipliers if the original sixteen fail.
        // A new prime coprime to 3*5*7*11 makes every assignment distinct.
        for (let p = 13n; ; p += 2n) {
          let prime = true;
          for (let d = 3n; d * d <= p; d += 2n) {
            if (p % d === 0n) { prime = false; break; }
          }
          if (prime) for (const k of base) yield k * p;
        }
      }
      const sequence = multiplierSequence();
      let sequenceIndex = 0;
      function nextMultiplier() {
        for (;;) {
          const k = sequence.next().value;
          if (sequenceIndex++ % workerCount === workerIndex) return k;
        }
      }
      function accept(candidate) {
        const divisor = gcd(candidate, n);
        if (divisor > 1n && divisor < n && n % divisor === 0n) {
          result.divisor = divisor;
          result.done = true;
          return true;
        }
        return false;
      }
      const states = Array.from({ length: Math.ceil(16 / workerCount) }, () => null);
      let cursor = 0, initialized = false;
      function initialize(k) {
        result.multipliers++;
        if (accept(k)) return null;
        const scaled = n * k, s = isqrt(scaled), q = scaled - s * s;
        if (q === 0n) { accept(s); return null; }
        return { k, scaled, s, cutoff: isqrt(2n * s), saved: new Set(),
          p: s, q, previousQ: 1n, even: true, reverse: false, batch: 0 };
      }
      function save(state, q) {
        if (q < state.cutoff * 2n * state.k) {
          const reduced = q / gcd(q, 2n * state.k);
          if (reduced < state.cutoff) {
            if (state.saved.size >= 50 && !state.saved.has(reduced)) return false;
            state.saved.add(reduced);
          }
        }
        return true;
      }
      function advanceForm(state) {
        const quotient = (state.s + state.p) / state.q;
        const p = quotient * state.q - state.p;
        const q = state.previousQ + (state.p - p) * quotient;
        result.iterations++;
        if (state.reverse) {
          result.reverseIterations++;
          if (p === state.p) { accept(state.q); return false; }
        } else {
          result.forwardIterations++;
          if (!save(state, state.q)) return false;
          if (state.even && q > 0n) {
            const low = q & 15n;
            if (low === 0n || low === 1n || low === 4n || low === 9n) {
              const root = isqrt(q);
              if (root * root === q && !state.saved.has(root)) {
                if (root === 1n) return false;
                // Inverse square-root form; the reverse traversal shares the
                // same work/progress accounting and yields as the forward one.
                state.p = p + root * ((state.s - p) / root);
                state.previousQ = root;
                const numerator = state.scaled - state.p * state.p;
                if (numerator <= 0n || numerator % root !== 0n) return false;
                state.q = numerator / root;
                state.reverse = true;
                return true;
              }
            }
          }
          state.even = !state.even;
        }
        if (q <= 0n) return false;
        state.previousQ = state.q; state.q = q; state.p = p;
        // A complete forward period without a proper square cannot help.
        // Retire it instead of cycling forever through already queued forms.
        if (!state.reverse && p === state.s && q === state.scaled - state.s * state.s) return false;
        return true;
      }
      return {
        advance() {
          if (result.done) return { ...result };
          const started = clock();
          try {
            if (!initialized) {
              initialized = true;
              if (n % 2n === 0n) accept(2n);
              else {
                const root = isqrt(n);
                if (root * root === n) accept(root);
                else if (isProbablePrime(n)) result.done = true;
              }
            }
            // Bounded messages keep cancellation/worker handoffs responsive.
            // No total time cutoff: this is the full small-composite path.
            for (let work = 0; !result.done && work < 4096; work++) {
              if (work % 32 === 0 && clock() - started >= 4) break;
              let state = states[cursor];
              if (!state) {
                state = states[cursor] = initialize(nextMultiplier());
                if (result.done) break;
                if (!state) {
                  result.failedMultipliers++;
                  cursor = (cursor + 1) % states.length;
                  continue;
                }
              }
              if (!advanceForm(state)) {
                if (!result.done) result.failedMultipliers++;
                states[cursor] = null;
                cursor = (cursor + 1) % states.length;
              } else if (++state.batch === 600) { // C: 300 even/odd pairs.
                state.batch = 0;
                cursor = (cursor + 1) % states.length;
              }
            }
          } finally { result.elapsedMs += Math.max(0, clock() - started); }
          return { ...result };
        }
      };
    }

    function squfof(n, onProgress = () => {}) {
      const job = createSqufofJob(n);
      let state, lastReport = -Infinity;
      do {
        state = job.advance();
        if (state.done || performance.now() - lastReport >= 200) {
          onProgress({ ...state, threads: 1 });
          lastReport = performance.now();
        }
      } while (!state.done);
      return { ...state, threads: 1 };
    }

    let nextSqufofJobId = 0;
    async function squfofParallel(n, onProgress, ports) {
      if (!ports.length) return squfof(n, onProgress);
      const threads = ports.length + 1, jobId = ++nextSqufofJobId;
      const local = createSqufofJob(n, { workerIndex: 0, workerCount: threads });
      const states = Array.from({ length: threads }, () => ({
        divisor: null, iterations: 0, forwardIterations: 0, reverseIterations: 0,
        multipliers: 0, failedMultipliers: 0, elapsedMs: 0, done: false
      }));
      const previous = ports.map(port => port.onmessage), pending = new Map();
      let requestId = 0;
      function request(index, operation, extra = {}) {
        return new Promise((resolve, reject) => {
          const id = ++requestId;
          pending.set(id, { index, resolve, reject });
          try { ports[index].postMessage({ type: "squfof", jobId, requestId: id, operation, ...extra }); }
          catch (error) { pending.delete(id); reject(error); }
        });
      }
      const counters = ["iterations", "forwardIterations", "reverseIterations", "multipliers", "failedMultipliers"];
      function accept(index, state) {
        if (!state || counters.some(key => !Number.isSafeInteger(state[key]) || state[key] < states[index][key]) ||
            state.iterations !== state.forwardIterations + state.reverseIterations ||
            state.failedMultipliers > state.multipliers ||
            !Number.isFinite(state.elapsedMs) || state.elapsedMs < states[index].elapsedMs ||
            typeof state.done !== "boolean" ||
            (state.divisor !== null && (typeof state.divisor !== "bigint" ||
              state.divisor <= 1n || state.divisor >= n || n % state.divisor !== 0n))) {
          throw new Error("Invalid SQUFOF worker result.");
        }
        states[index] = state;
      }
      try {
        ports.forEach((port, index) => {
          port.onmessage = ({ data }) => {
            if (data.jobId !== jobId || (data.type !== "squfof-result" && data.type !== "squfof-error")) return;
            const call = pending.get(data.requestId);
            if (!call || call.index !== index) return;
            pending.delete(data.requestId);
            if (data.type === "squfof-error") call.reject(new Error(data.message));
            else call.resolve(data.result);
          };
        });
        await Promise.all(ports.map((_, i) => request(i, "init", {
          n, options: { workerIndex: i + 1, workerCount: threads }
        })));
        let lastReport = -Infinity;
        for (;;) {
          const replies = ports.map((_, i) => states[i + 1].done ? Promise.resolve() :
            request(i, "step").then(state => accept(i + 1, state)));
          replies.push(Promise.resolve().then(() => { if (!states[0].done) accept(0, local.advance()); }));
          await Promise.all(replies);
          const result = { divisor: states.find(s => s.divisor !== null)?.divisor ?? null, threads,
            elapsedMs: states.reduce((sum, s) => sum + s.elapsedMs, 0),
            workerIterations: states.map(s => s.iterations) };
          for (const key of counters) result[key] = states.reduce((sum, s) => sum + s[key], 0);
          result.done = result.divisor !== null || states.every(s => s.done);
          if (result.done || performance.now() - lastReport >= 200) {
            onProgress(result);
            lastReport = performance.now();
          }
          if (result.done) return result;
        }
      } finally {
        await Promise.allSettled(ports.map((_, i) => request(i, "release")));
        ports.forEach((port, i) => { port.onmessage = previous[i]; });
      }
    }

    function squfofSmall(n) {
      // A bounded Number version of msieve's racing SQUFOF. Scaled inputs
      // stay below 2^52, so every integer product/difference is exact.
      // The returned divisor is additionally checked with BigInt by callers.
      if (!Number.isSafeInteger(n) || n < 4 || n >= 2 ** 52) return 0;
      const multipliers = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155];
      const states = [];
      for (const k of multipliers) {
        const scaled = n * k;
        if (scaled >= 2 ** 52) break;
        const s = Math.floor(Math.sqrt(scaled)), q = scaled - s * s;
        if (!q) { const g = gcdSmall(s, n); if (g > 1 && g < n) return g; continue; }
        states.push({ k, scaled, s, cutoff: Math.floor(Math.sqrt(2 * s)), saved: [], q0: 1, p1: s, q1: q });
      }
      for (let work = 0; work < 20000 && states.length;) {
        for (let j = states.length - 1; j >= 0 && work < 20000; j--) {
          const state = states[j], { k, scaled, s, cutoff, saved } = state;
          let { q0, p1, q1 } = state, p0 = 0, square = 0, failed = false;
          function save(q) {
            if (q < cutoff * 2 * k) {
              const reduced = q / gcdSmall(q, 2 * k);
              if (reduced < cutoff) saved.push(reduced);
            }
            if (saved.length >= 50) failed = true;
          }
          for (let i = 0; i < 300; i++, work++) {
            let quotient = Math.floor((s + p1) / q1);
            p0 = quotient * q1 - p1;
            q0 += (p1 - p0) * quotient;
            if (q0 <= 0) { failed = true; break; }
            save(q1);
            const low = q0 & 15;
            if (low === 0 || low === 1 || low === 4 || low === 9) {
              const r = Math.floor(Math.sqrt(q0));
              if (r * r === q0 && !saved.includes(r)) { square = r; break; }
            }
            quotient = Math.floor((s + p0) / q0);
            p1 = quotient * q0 - p0;
            q1 += (p0 - p1) * quotient;
            if (q1 <= 0) { failed = true; break; }
            save(q0);
            if (failed) break;
          }
          if (square > 1) {
            q0 = square;
            p1 = p0 + square * Math.floor((s - p0) / square);
            q1 = (scaled - p1 * p1) / square;
            for (let i = 0; i < 20000 && q1 > 0; i++) {
              const quotient = Math.floor((s + p1) / q1);
              p0 = quotient * q1 - p1;
              const nextQ = q0 + (p1 - p0) * quotient;
              if (p0 === p1) {
                const g = gcdSmall(q1, n);
                if (g > 1 && g < n && n % g === 0) return g;
                break;
              }
              q0 = q1; q1 = nextQ; p1 = p0;
            }
            failed = true;
          }
          if (square === 1) failed = true;
          if (failed) states.splice(j, 1);
          else Object.assign(state, { q0, p1, q1 });
        }
      }
      return 0;
    }

    function splitLargePrimes(remainder, bound, maxDouble, minDouble) {
      if (remainder <= bound) return [1, Number(remainder)];
      if (remainder < minDouble || remainder > maxDouble) return null;
      // Prime residuals above the single bound cannot contribute. Fermat
      // rejection is only a speed filter; every accepted split is exact.
      if (powMod(2n, remainder - 1n, remainder) === 1n) return null;
      const factor = squfofSmall(Number(remainder));
      if (factor <= 1 || remainder % BigInt(factor) !== 0n) return null;
      const other = remainder / BigInt(factor);
      if (BigInt(factor) > bound || other > bound) return null;
      return factor < other ? [factor, Number(other)] : [Number(other), factor];
    }

    function packExponents(exponents) {
      return Uint32Array.from(exponents, ([column, count]) => {
        if (count >= 256 || column >= 2 ** 24) throw new Error("Relation exponent encoding overflow.");
        return column * 256 + count;
      });
    }

    function unpackRelation(row) {
      if (row.exponents instanceof Uint32Array) {
        return { ...row, exponents: Array.from(row.exponents, value => [value >>> 8, value & 255]) };
      }
      return row;
    }

    function createLargePrimeGraph(n) {
      // The spanning forest retains partials for reuse by multiple independent
      // cycles. Union-by-size bounds re-rooting work; no path search is needed
      // until an incoming edge actually closes a cycle.
      const vertices = new Map([[1, 0]]), primes = [1], parents = [-1], parentEdges = [-1];
      const roots = [0], sizes = [1], edges = [], seen = new Set();
      let crossWorkerMatches = 0;
      function vertex(p) {
        let v = vertices.get(p);
        if (v !== undefined) return v;
        v = primes.length;
        vertices.set(p, v); primes.push(p); parents.push(-1); parentEdges.push(-1);
        roots.push(v); sizes.push(1);
        return v;
      }
      function root(v) {
        while (roots[v] !== v) { roots[v] = roots[roots[v]]; v = roots[v]; }
        return v;
      }
      function add(p, q, row, worker = 0) {
        const key = row.x < n - row.x ? row.x : n - row.x;
        if (seen.has(key)) return { duplicate: true };
        seen.add(key);
        if (p === 1 && q === 1) return { row: unpackRelation(row), full: true };
        let u = vertex(p), v = vertex(q), ru = root(u), rv = root(v);
        if (ru !== rv) {
          if (sizes[ru] > sizes[rv]) { [u, v] = [v, u]; [ru, rv] = [rv, ru]; }
          let prev = -1, prevEdge = -1, current = u;
          while (current !== -1) {
            const next = parents[current], edge = parentEdges[current];
            parents[current] = prev; parentEdges[current] = prevEdge;
            prev = current; prevEdge = edge; current = next;
          }
          parents[u] = v; parentEdges[u] = edges.length;
          edges.push({ row, worker });
          roots[ru] = rv; sizes[rv] += sizes[ru];
          return { partial: true };
        }
        const ancestors = new Set();
        for (let a = u; a !== -1; a = parents[a]) ancestors.add(a);
        let common = v;
        while (!ancestors.has(common)) common = parents[common];
        const cycle = [{ row, worker }], large = [];
        for (let a = u; a !== common; a = parents[a]) {
          cycle.push(edges[parentEdges[a]]); large.push(primes[a]);
        }
        for (let a = v; a !== common; a = parents[a]) {
          cycle.push(edges[parentEdges[a]]); large.push(primes[a]);
        }
        large.push(primes[common]);
        let x = 1n, square = 1n;
        const counts = new Map();
        for (const edge of cycle) {
          x = x * edge.row.x % n;
          square = square * edge.row.square % n;
          const packed = edge.row.exponents;
          if (packed instanceof Uint32Array) {
            for (const value of packed) {
              const c = value >>> 8;
              counts.set(c, (counts.get(c) || 0) + (value & 255));
            }
          } else {
            for (const [c, count] of packed) counts.set(c, (counts.get(c) || 0) + count);
          }
        }
        for (const p of large) square = square * BigInt(p) % n;
        if (cycle.some(edge => edge.worker !== worker)) crossWorkerMatches++;
        return { row: { x, square, exponents: Array.from(counts) }, full: false };
      }
      return { add, get retained() { return edges.length; },
        get crossWorkerMatches() { return crossWorkerMatches; } };
    }

    function createSparseSolver(columns, rows, extract) {
      const pivots = new Array(columns), references = [];
      const threshold = 192;
      function xorSparse(a, b) {
        const out = [];
        let i = 0, j = 0;
        while (i < a.length && j < b.length) {
          if (a[i] < b[j]) out.push(a[i++]);
          else if (b[j] < a[i]) out.push(b[j++]);
          else { i++; j++; }
        }
        while (i < a.length) out.push(a[i++]);
        while (j < b.length) out.push(b[j++]);
        return out;
      }
      return function add(index) {
        let sparse = rows[index].exponents.filter(pair => pair[1] & 1).map(pair => pair[0]).sort((a, b) => a - b);
        let dense = null, high = -1;
        const used = [];
        function densify() {
          high = sparse[sparse.length - 1] >>> 5;
          dense = new Uint32Array(high + 1);
          for (const c of sparse) dense[c >>> 5] ^= 1 << (c & 31);
          sparse = null;
        }
        for (;;) {
          if (dense) while (high >= 0 && !dense[high]) high--;
          const column = dense ? high < 0 ? -1 : high * 32 + 31 - Math.clz32(dense[high]) :
            sparse.length ? sparse[sparse.length - 1] : -1;
          if (column < 0) break;
          const pivot = pivots[column];
          if (!pivot) {
            pivots[column] = { index, sparse: sparse ? Uint32Array.from(sparse) : null,
              dense: dense ? dense.slice(0, high + 1) : null };
            references[index] = Uint32Array.from(used);
            return null;
          }
          used.push(pivot.index);
          if (sparse && pivot.sparse) {
            sparse = xorSparse(sparse, pivot.sparse);
            if (sparse.length > threshold) densify();
          } else {
            if (!dense) densify();
            if (pivot.sparse) {
              for (const c of pivot.sparse) dense[c >>> 5] ^= 1 << (c & 31);
            } else {
              for (let w = 0; w <= high; w++) dense[w] ^= pivot.dense[w];
            }
          }
        }
        references[index] = Uint32Array.from(used);
        // Expand the elimination DAG only when there is a dependency. This
        // avoids a factor-base-sized combination bit vector at every pivot.
        const flags = new Uint8Array(index + 1);
        flags[index] = 1;
        for (let i = index; i >= 0; i--) {
          if (flags[i]) for (const j of references[i]) flags[j] ^= 1;
        }
        const combination = new Uint32Array(Math.ceil((index + 1) / 32));
        for (let i = 0; i <= index; i++) if (flags[i]) combination[i >>> 5] |= 1 << (i & 31);
        return extract(combination);
      };
    }

    function createSiqs(n, report, stats, threads = 1, overrides = {}, matrixWorkers = []) {
      const settings = { ...siqsParameters(n), ...overrides };
      const multiplier = settings.largeSieve ? chooseMultiplier(n) : 1;
      const sieveN = n * BigInt(multiplier);
      const M = settings.halfInterval, length = 2 * M;
      const digits = n.toString().length;
      beginSieve(n, report, stats);
      // This upper bound contains more than twice the requested prime count;
      // extend it if the quadratic-residue filter happens to leave too few.
      let bound = Math.ceil(3 * settings.size * Math.log(3 * settings.size));
      const fb = [];
      let lastPrime = 1;
      while (fb.length < settings.size) {
        if (bound >= 2 ** 26) throw new Error("The SIQS factor-base limit was exceeded.");
        for (const p of primesThrough(bound)) {
          if (p <= lastPrime) continue;
          const big = BigInt(p);
          if (n % big === 0n) return { divisor: big };
          const root = sqrtMod(sieveN, big);
          if (root !== null) fb.push({ p, big, root: Number(root), weight: Math.log2(p) });
          lastPrime = p;
          if (fb.length === settings.size) break;
        }
        bound *= 2;
      }
      const columns = fb.length + 1;
      let estimatedRelations = columns + (settings.largeSieve ? 96 : Math.max(16, Math.ceil(columns * 0.05)));
      const largeBound = BigInt(fb[fb.length - 1].p * settings.largeMultiplier);
      const doubleBound = settings.doubleLargePrimes ?
        largeBound * BigInt(Math.floor(Number(largeBound) ** 0.8)) : 0n;
      const candidateBound = doubleBound > 0n ?
        (doubleBound < (1n << 52n) ? doubleBound : (1n << 52n) - 1n) : largeBound;
      const minDouble = fb[fb.length - 1].big ** 2n;
      const targetA = isqrt(2n * sieveN) / BigInt(M);
      const numAFactors = Math.max(2, Math.min(12, Math.round(log2(targetA) / 11)));
      const desiredPrime = Math.pow(2, log2(targetA) / numAFactors);
      const pool = [];
      for (let j = 1; j < fb.length; j++) {
        if (fb[j].root && fb[j].p >= desiredPrime * 0.5 && fb[j].p <= desiredPrime * 1.5) pool.push(j);
      }
      if (pool.length < numAFactors + 1) throw new Error("Insufficient SIQS polynomial factors.");
      // Deterministic diversity of A values, independent of the MR random bases.
      let seed = Number(n & 0xffffffffn) ^ Number((n >> 32n) & 0xffffffffn) ^ 0x9e3779b9;
      if (!seed) seed = 1;
      function randomIndex(size) {
        seed ^= seed << 13; seed ^= seed >>> 17; seed ^= seed << 5;
        return Math.floor((seed >>> 0) / 4294967296 * size);
      }
      const usedA = new Set();

      // Modular square roots for small prime powers are shared by every A/B.
      const entries = [];
      for (let j = 0; j < fb.length; j++) {
        const { p, big, root, weight } = fb[j];
        let modulus = p;
        let roots = p === 2 || root === 0 ? [root] : [root, p - root];
        entries.push({ modulus, base: j, roots, weight, power: false });
        if (sieveN % big === 0n) continue; // Squarefree multiplier: no roots modulo p².
        while (modulus * p <= length) {
          const nextModulus = modulus * p, lifted = [];
          for (const r of roots) {
            if (p === 2) {
              for (const candidate of [r, r + modulus]) {
                if ((BigInt(candidate) ** 2n - sieveN) % BigInt(nextModulus) === 0n) lifted.push(candidate);
              }
            } else {
              const step = mod((sieveN - BigInt(r) ** 2n) / BigInt(modulus) * inverse(2n * BigInt(r), big), big);
              lifted.push(r + Number(step) * modulus);
            }
          }
          roots = lifted;
          modulus = nextModulus;
          if (!roots.length) break;
          entries.push({ modulus, base: j, roots, weight, power: true });
        }
      }

      const partials = new Map();
      const graph = settings.doubleLargePrimes ? createLargePrimeGraph(n) : null;
      const maxPartials = 50000;
      const maxRows = columns + 256;
      const parityWords = Math.ceil(columns / 32), combinationWords = Math.ceil(maxRows / 32);
      let rows = [], pivots = new Array(columns);
      let smoothRelations = 0, combinedRelations = 0, partialsFound = 0;
      const useBlockLanczos = settings.largeSieve;
      const search = { digits, factorBaseSize: fb.length, halfInterval: M,
        largePrimeBound: String(largeBound), polynomialGroups: 0, polynomials: 0,
        candidates: 0, linearAlgebraMs: 0, threads,
        workerPolynomials: Array(threads).fill(0), crossWorkerMatches: 0,
        multiplier, maxLargePrimes: graph ? 2 : 1, doubleLargePrimeBound: String(candidateBound),
        singleLargePrimeRelations: 0, doubleLargePrimeRelations: 0,
        linearAlgebraMethod: useBlockLanczos ? "block-lanczos" : "packed-gaussian",
        matrices: [] };
      stats.searches.push(search);

      function progress(stage = "Sieving polynomials", extra = {}) {
        const counts = { fullRelations: rows.length, smoothRelations, combinedRelations,
          partialRelations: graph ? graph.retained : partials.size, partialsFound, estimatedRelations,
          maxLargePrimes: graph ? 2 : 1 };
        Object.assign(stats, counts);
        Object.assign(search, counts);
        report({ stage, digits, threads, polynomials: search.polynomials,
          ...counts, factorBaseSize: fb.length, multiplier,
          doubleLargePrimeRelations: search.doubleLargePrimeRelations, linearAlgebraMs: search.linearAlgebraMs,
          ...extra });
      }
      progress();

      function extract(combination) {
        let X = 1n, Y = 1n;
        const powers = new Uint32Array(columns);
        for (let word = 0; word < combination.length; word++) {
          let bits = combination[word];
          while (bits) {
            const bit = 31 - Math.clz32(bits);
            const row = rows[word * 32 + bit];
            X = X * row.x % n;
            Y = Y * row.square % n;
            for (const [column, count] of row.exponents) powers[column] += count;
            bits = (bits ^ (1 << bit)) >>> 0;
          }
        }
        for (let column = 0; column < columns; column++) {
          if (powers[column] & 1) throw new Error("Invalid SIQS dependency.");
          if (column && powers[column]) {
            Y = Y * powMod(fb[column - 1].big, BigInt(powers[column] / 2), n) % n;
          }
        }
        for (const difference of [X - Y, X + Y]) {
          const divisor = gcd(difference, n);
          if (divisor > 1n && divisor < n) return divisor;
        }
        return null;
      }

      function addRelation(row) {
        const index = rows.length;
        rows.push(row);
        stats.relations++;
        // Large searches collect first, then filter and solve in a batch.
        // Helpers have only one message in flight while the coordinator is
        // solving; Clear can terminate the entire worker pool at any stage.
        if (useBlockLanczos) return null;
        const parity = new Uint32Array(parityWords);
        for (const [column, count] of row.exponents) {
          if (count & 1) parity[column >>> 5] ^= 1 << (column & 31);
        }
        const combination = new Uint32Array(combinationWords);
        combination[index >>> 5] = 1 << (index & 31);
        // No BigInt-to-string pivot searches or growing BigInt combinations.
        // Packed XORs and clz32 visit only the active portion of each vector.
        let highWord = parityWords - 1;
        while (highWord >= 0) {
          if (!parity[highWord]) { highWord--; continue; }
          const column = highWord * 32 + 31 - Math.clz32(parity[highWord]);
          const pivot = pivots[column];
          if (!pivot) {
            pivots[column] = { parity, combination, words: (index >>> 5) + 1 };
            return null;
          }
          for (let w = 0; w <= highWord; w++) parity[w] ^= pivot.parity[w];
          for (let w = 0; w < pivot.words; w++) combination[w] ^= pivot.combination[w];
        }
        return extract(combination);
      }

      function* matrixSteps() {
        const started = performance.now();
        progress("Filtering matrix", { matrixRows: columns, matrixColumns: rows.length });
        const original = rows.map(row => Uint32Array.from(
          row.exponents.filter(pair => pair[1] & 1).map(pair => pair[0]).sort((a, b) => a - b)));
        const filtered = matrixAlgebra.filterMatrix(columns, original);
        const record = { ...filtered.stats, filteringMs: performance.now() - started,
          solvingMs: 0, extractionMs: 0, method: "block-lanczos", dependenciesTried: 0 };
        search.matrices.push(record);
        function tryDependency(vector) {
          if (!matrixAlgebra.verifyDependency(columns, original, vector)) {
            throw new Error("Matrix filtering produced an invalid dependency.");
          }
          record.dependenciesTried++;
          const began = performance.now(), divisor = extract(vector);
          record.extractionMs += performance.now() - began;
          return divisor;
        }
        try {
          // Zero columns produced by filtering already give dependencies.
          for (const indices of filtered.dependencies) {
            const vector = new Uint32Array(Math.ceil(rows.length / 32));
            for (const index of indices) vector[index >>> 5] ^= 1 << (index & 31);
            const divisor = tryDependency(vector);
            if (divisor !== null) return divisor;
          }
          if (filtered.columns.length > filtered.rowCount) {
            const began = performance.now();
            const shape = { matrixRows: filtered.rowCount, matrixColumns: filtered.columns.length };
            const lanczosShape = { ...shape, matrixThreads: matrixWorkers.length + 1 };
            progress("Block Lanczos", lanczosShape);
            const solved = yield { rowCount: filtered.rowCount, columns: filtered.columns, options: {
              onProgress: state => progress("Block Lanczos", { ...lanczosShape,
                matrixIteration: state.iterations, matrixAttempt: state.attempt,
                matrixDimensionsSolved: state.dimensionsSolved, matrixDimensions: state.dimensions })
            } };
            record.solvingMs = performance.now() - began;
            Object.assign(record, solved.stats);
            progress("Extracting factors", shape);
            for (const dependency of solved.dependencies) {
              const originalVector = matrixAlgebra.expandDependency(dependency, filtered.groups, rows.length);
              const divisor = tryDependency(originalVector);
              if (divisor !== null) return divisor;
            }
          }
          // Filtering or an unlucky Krylov block can leave too little useful
          // rank. Keep every collected relation and continue sieving; never
          // recycle this large matrix or accept an unchecked dependency.
          estimatedRelations = rows.length + Math.max(64, Math.ceil(columns * 0.02));
          return null;
        } finally {
          search.linearAlgebraMs += performance.now() - started;
        }
      }
      function solveMatrix() {
        const steps = matrixSteps(), step = steps.next();
        if (step.done) return step.value; // Filtering alone found a factor or needs more relations.
        const { rowCount, columns: filteredColumns, options } = step.value;
        if (matrixWorkers.length) {
          return matrixAlgebra.blockLanczosParallel(rowCount, filteredColumns, { ...options, workers: matrixWorkers })
            .then(result => steps.next(result).value, error => steps.throw(error));
        }
        try {
          return steps.next(matrixAlgebra.blockLanczos(rowCount, filteredColumns, options)).value;
        } catch (error) { return steps.throw(error); }
      }

      function nextFamily() {
        let A, selected;
        for (;;) {
          selected = [];
          A = 1n;
          while (selected.length < numAFactors - 1) {
            const j = pool[randomIndex(pool.length)];
            if (!selected.includes(j)) { selected.push(j); A *= fb[j].big; }
          }
          const idealLast = Number(targetA / A); // A bounded factor-base prime estimate.
          let lo = 1, hi = fb.length - 1;
          while (lo < hi) {
            const mid = (lo + hi) >>> 1;
            if (fb[mid].p < idealLast) lo = mid + 1; else hi = mid;
          }
          let best = -1, distance = Infinity;
          for (let j = Math.max(1, lo - numAFactors); j <= Math.min(fb.length - 1, lo + numAFactors); j++) {
            if (fb[j].root && !selected.includes(j) && Math.abs(fb[j].p - idealLast) < distance) {
              best = j; distance = Math.abs(fb[j].p - idealLast);
            }
          }
          if (best < 0) continue;
          A *= fb[best].big;
          if (A * 2n < targetA || A > targetA * 2n || usedA.has(A)) continue;
          selected.push(best);
          if (usedA.size >= 10000) usedA.clear();
          usedA.add(A);
          break;
        }
        return selected;
      }

      function accept(batch, worker = 0) {
        search.polynomials += batch.polynomials;
        search.polynomialGroups += batch.polynomialGroups;
        search.candidates += batch.candidates;
        search.workerPolynomials[worker] += batch.polynomials;
        stats.polynomials += batch.polynomials;
        stats.polynomialGroups += batch.polynomialGroups;
        stats.reusedPolynomials += batch.polynomials - batch.polynomialGroups;
        for (const relation of batch.relations) {
          let { row, remainder } = relation;
          if (graph) {
            const [p, q] = relation.largePrimes;
            const combined = graph.add(p, q, row, worker);
            if (combined.duplicate) continue;
            if (p !== 1 || q !== 1) {
              partialsFound++;
              if (p === 1) search.singleLargePrimeRelations++;
              else search.doubleLargePrimeRelations++;
            }
            if (!combined.row) continue;
            row = combined.row;
            if (combined.full) smoothRelations++; else combinedRelations++;
            search.crossWorkerMatches = graph.crossWorkerMatches;
          } else if (remainder !== 1n) {
            row = unpackRelation(row);
            partialsFound++;
            const previous = partials.get(remainder);
            if (!previous) {
              if (partials.size >= maxPartials) partials.delete(partials.keys().next().value);
              partials.set(remainder, { row, worker });
              continue;
            }
            partials.delete(remainder);
            combinedRelations++;
            if (previous.worker !== worker) search.crossWorkerMatches++;
            row = { x: previous.row.x * row.x % n, square: remainder,
              exponents: mergeExponents(previous.row.exponents, row.exponents) };
          } else { row = unpackRelation(row); smoothRelations++; }
          const start = performance.now();
          const divisor = addRelation(row);
          search.linearAlgebraMs += performance.now() - start;
          if (divisor !== null) { progress(); return divisor; }
          if (!useBlockLanczos && rows.length === maxRows) {
            rows = [];
            pivots = new Array(columns);
            smoothRelations = combinedRelations = 0;
          }
        }
        if (useBlockLanczos && rows.length >= estimatedRelations) {
          const divisor = solveMatrix();
          if (divisor instanceof Promise) {
            return divisor.then(value => {
              progress(value !== null ? "Factors recovered" : "Sieving polynomials");
              return value;
            });
          }
          if (divisor !== null) { progress("Factors recovered"); return divisor; }
        }
        progress();
        return null;
      }

      return { context: { n, sieveN, M, fb, entries, numAFactors, largeBound, candidateBound, minDouble,
        largeSieve: settings.largeSieve, doubleLargePrimes: settings.doubleLargePrimes },
        nextFamily, accept };
    }

    // Each worker owns its sieve buffers and reuses them for every assigned
    // family. The coordinator assigns distinct A values; raw full and partial
    // relations all return to the same collector and matrix.
    function createSieve(context) {
      const { n, sieveN, M, fb, entries, numAFactors, largeBound, candidateBound, minDouble,
        largeSieve, doubleLargePrimes } = context;
      const length = 2 * M, largeLog = log2(candidateBound);
      const scores = new Uint16Array(length);
      const smallPattern = new Uint16Array(32 * 9 * 5);
      const candidateMap = largeSieve ? new Uint32Array(length) : null;
      const basePrimes = Uint32Array.from(fb, entry => entry.p);

      // Keep the hot loops outside the polynomial generator and use flat
      // typed arrays. This avoids object/iterator overhead at every sieve hit.
      function shiftRoots(roots, moduli, shifts, direction) {
        for (let i = 0; i < roots.length; i++) {
          const p = moduli[i];
          let value = roots[i] - direction * shifts[i];
          if (value < 0) value += p;
          else if (value >= p) value -= p;
          roots[i] = value;
        }
      }
      function markSieve(roots, moduli, weights, indices) {
        for (let j = 0; j < indices.length; j++) {
          const entry = indices[j], p = moduli[entry], weight = weights[entry];
          let index = roots[entry];
          const step4 = 4 * p, stop4 = length - 3 * p;
          for (; index < stop4; index += step4) {
            scores[index] += weight;
            scores[index + p] += weight;
            scores[index + 2 * p] += weight;
            scores[index + 3 * p] += weight;
          }
          for (; index < length; index += p) scores[index] += weight;
        }
      }
      function candidatePositions(a, b, c) {
        const positions = [];
        const coarseValue = Math.abs(c) / 64;
        const coarseScore = (Math.log2(coarseValue) - largeLog - 1.5) * 8;
        for (let i = 0; i < length; i++) {
          const offset = i - M;
          const approximate = Math.abs((a * offset + b) * offset + c);
          if (approximate >= coarseValue && scores[i] < coarseScore) continue;
          if (scores[i] / 8 + largeLog + 1.5 >= Math.log2(approximate)) positions.push(i);
        }
        return positions;
      }
      function divideCandidate(candidate, base, remainders, exponentLists) {
        const big = fb[base].big;
        let value = remainders[candidate], count = 0;
        while (value % big === 0n) { value /= big; count++; }
        if (count) {
          remainders[candidate] = value;
          exponentLists[candidate].push([base + 1, count]);
        }
      }
      function resieve(positions, remainders, exponentLists, roots, starts, rootCounts) {
        const directLimit = length / Math.max(1, positions.length);
        for (let base = 0; base < basePrimes.length; base++) {
          const p = basePrimes[base], start = starts[base];
          const r0 = roots[start], r1 = rootCounts[base] > 1 ? roots[start + 1] : -1;
          if (p < directLimit) {
            for (let i = 0; i < positions.length; i++) {
              const residue = positions[i] % p;
              if (residue === r0 || residue === r1) divideCandidate(i, base, remainders, exponentLists);
            }
          } else {
            for (let r = 0; r < rootCounts[base]; r++) {
              for (let index = roots[start + r]; index < length; index += p) {
                const candidate = candidateMap[index];
                if (candidate) divideCandidate(candidate - 1, base, remainders, exponentLists);
              }
            }
          }
        }
      }
      function* family(selected) {
        const A = selected.reduce((a, j) => a * fb[j].big, 1n);
        const selectedSet = new Set(selected);
        const aFactor = new Uint8Array(fb.length);
        for (const j of selected) aFactor[j] = 1;
        const pieces = selected.map(j => {
          const { big } = fb[j], quotient = A / big;
          let t = mod(BigInt(fb[j].root) * inverse(quotient, big), big);
          if (t > big / 2n) t = big - t;
          return quotient * t;
        });
        let B = pieces.reduce((a, b) => a + b, 0n);
        const baseRoots = new Array(fb.length);
        const channels = entries.map(entry => {
          const p = entry.modulus;
          const special = selectedSet.has(entry.base);
          const invA = special ? 0 : inverseSmall(Number(A % BigInt(p)), p);
          const b = Number(mod(B, BigInt(p)));
          const roots = special ? (entry.power ? [] : [0]) :
            entry.roots.map(r => modSmall((r - b) * invA + M, p));
          if (!entry.power) baseRoots[entry.base] = roots;
          return { ...entry, roots, special, integerWeight: Math.round(entry.weight * 8),
            shifts: special ? [] : pieces.slice(0, -1).map(piece => Number(2n * piece % BigInt(p)) * invA % p) };
        });
        const specialChannels = channels.filter(channel => channel.special && !channel.power);
        const baseChannels = channels.filter(channel => !channel.power);
        const inPattern = channel => {
          const prime = fb[channel.base].p;
          return (prime === 2 && channel.modulus <= 32) ||
            (prime === 3 && channel.modulus <= 9) || (prime === 5 && channel.modulus === 5);
        };
        const patternChannels = channels.filter(inPattern);
        const sieveChannels = channels.filter(channel => !inPattern(channel));
        let flatRoots, flatModuli, flatWeights, flatShifts, flatIndices, baseStarts, rootCounts;
        if (largeSieve) {
          const count = channels.reduce((sum, channel) => sum + channel.roots.length, 0);
          flatRoots = new Uint32Array(count);
          flatModuli = new Uint32Array(count);
          flatWeights = new Uint16Array(count);
          flatShifts = pieces.slice(0, -1).map(() => new Uint32Array(count));
          baseStarts = new Uint32Array(fb.length);
          rootCounts = new Uint8Array(fb.length);
          const indices = [];
          let at = 0;
          for (const channel of channels) {
            channel.start = at;
            if (!channel.power) {
              baseStarts[channel.base] = at;
              rootCounts[channel.base] = channel.roots.length;
            }
            for (const root of channel.roots) {
              flatRoots[at] = root;
              flatModuli[at] = channel.modulus;
              flatWeights[at] = channel.integerWeight;
              for (let bit = 0; bit < flatShifts.length; bit++) flatShifts[bit][at] = channel.shifts[bit] || 0;
              if (!inPattern(channel)) indices.push(at);
              at++;
            }
          }
          flatIndices = Uint32Array.from(indices);
        }
        const batchSize = 2 ** (numAFactors - 1);
        for (let polynomial = 0; polynomial < batchSize; polynomial++) {
          if (polynomial) {
            const bit = 31 - Math.clz32(polynomial & -polynomial);
            const direction = (polynomial & (1 << (bit + 1))) ? 1 : -1;
            B += BigInt(direction) * 2n * pieces[bit];
            if (largeSieve) shiftRoots(flatRoots, flatModuli, flatShifts[bit], direction);
            else for (const channel of channels) {
              if (channel.special) continue;
              const delta = direction * channel.shifts[bit], p = channel.modulus;
              for (let r = 0; r < channel.roots.length; r++) {
                let next = channel.roots[r] - delta;
                if (next < 0) next += p;
                else if (next >= p) next -= p;
                channel.roots[r] = next;
              }
            }
          }
          const numerator = B * B - sieveN;
          if (numerator % A !== 0n) throw new Error("Invalid SIQS polynomial.");
          const C = numerator / A;
          for (const channel of specialChannels) {
            const p = channel.modulus;
            const big = fb[channel.base].big;
            channel.roots[0] = modSmall(-Number(mod(C, big)) * inverseSmall(Number(mod(2n * B, big)), p) + M, p);
            if (largeSieve) flatRoots[channel.start] = channel.roots[0];
          }
          // The densest sieve updates repeat every 1440 positions. Build that
          // short pattern once and expand it with native typed-array copies.
          smallPattern.fill(0);
          for (const channel of patternChannels) {
            for (let r = 0; r < channel.roots.length; r++) {
              const root = largeSieve ? flatRoots[channel.start + r] : channel.roots[r];
              for (let i = root; i < smallPattern.length; i += channel.modulus) smallPattern[i] += channel.integerWeight;
            }
          }
          scores.set(smallPattern);
          for (let filled = smallPattern.length; filled < length;) {
            const count = Math.min(filled, length - filled);
            scores.copyWithin(filled, 0, count);
            filled += count;
          }
          if (largeSieve) markSieve(flatRoots, flatModuli, flatWeights, flatIndices);
          else for (const channel of sieveChannels) {
            const p = channel.modulus, weight = channel.integerWeight;
            const step4 = 4 * p, stop4 = length - 3 * p;
            for (const root of channel.roots) {
              let index = root;
              for (; index < stop4; index += step4) {
                scores[index] += weight;
                scores[index + p] += weight;
                scores[index + 2 * p] += weight;
                scores[index + 3 * p] += weight;
              }
              for (; index < length; index += p) scores[index] += weight;
            }
          }
          const relations = [];
          const positions = [], remainders = [], exponentLists = [];
          if (largeSieve) candidateMap.fill(0);
          let candidates = 0;
          const aApprox = Number(A), bApprox = Number(2n * B), cApprox = Number(C);
          if (largeSieve) {
            const offsets = candidatePositions(aApprox, bApprox, cApprox);
            candidates = offsets.length;
            for (const index of offsets) {
              const x = BigInt(index - M);
              let remainder = (A * x + 2n * B) * x + C;
              if (remainder === 0n) continue;
              const exponents = [];
              if (remainder < 0n) { remainder = -remainder; exponents.push([0, 1]); }
              positions.push(index);
              remainders.push(remainder);
              exponentLists.push(exponents);
              candidateMap[index] = positions.length;
            }
          } else {
          const coarseValue = Math.abs(cApprox) / 64;
          const coarseScore = (Math.log2(coarseValue) - largeLog - 1.5) * 8;
          for (let index = 0; index < length; index++) {
            const offset = index - M;
            // Floating point only selects candidates; exact division below
            // determines all accepted relations and factors.
            const approximate = Math.abs((aApprox * offset + bApprox) * offset + cApprox);
            // A conservative cheap check avoids a logarithm at almost every
            // sieve position, while retaining candidates close to a Q root.
            if (approximate >= coarseValue && scores[index] < coarseScore) continue;
            if (scores[index] / 8 + largeLog + 1.5 < Math.log2(approximate)) continue;
            candidates++;
            const x = BigInt(offset);
            let remainder = (A * x + 2n * B) * x + C;
            if (remainder === 0n) continue;
            const exponents = [];
            if (remainder < 0n) { remainder = -remainder; exponents.push([0, 1]); }
            for (let j = 0; j < fb.length; j++) {
              const { p, big } = fb[j];
              let count = aFactor[j]; // X² = A*Q(x), include A.
              const residue = index % p, roots = baseRoots[j];
              if (remainder !== 1n && (residue === roots[0] || residue === roots[1])) {
                while (remainder % big === 0n) { remainder /= big; count++; }
              }
              if (count) exponents.push([j + 1, count]);
            }
            if (remainder > largeBound) continue;
            relations.push({ remainder,
              row: { x: mod(A * x + B, n), square: 1n, exponents } });
          }
          }
          if (largeSieve) {
            // Re-sieve only the marked candidates. This avoids scanning all
            // 55,000 factor-base primes separately for every candidate. For
            // small primes, testing the few candidate offsets is cheaper than
            // revisiting their many thousands of sieve hits.
            resieve(positions, remainders, exponentLists, flatRoots, baseStarts, rootCounts);
            for (let i = 0; i < positions.length; i++) {
              const remainder = remainders[i];
              if (remainder > candidateBound) continue;
              const largePrimes = doubleLargePrimes ?
                splitLargePrimes(remainder, largeBound, candidateBound, minDouble) : null;
              if (doubleLargePrimes && largePrimes === null) continue;
              const exponents = exponentLists[i];
              for (const j of selected) {
                const pair = exponents.find(pair => pair[0] === j + 1);
                if (pair) pair[1]++; else exponents.push([j + 1, 1]);
              }
              relations.push({ remainder, largePrimes,
                row: { x: mod(A * BigInt(positions[i] - M) + B, n), square: 1n,
                  exponents: packExponents(exponents) } });
            }
          }
          yield { relations, polynomials: 1, polynomialGroups: polynomial === 0 ? 1 : 0,
            candidates, done: polynomial + 1 === batchSize };
        }
      }
      return { family };
    }

    // Limit both latency and message frequency. A helper sends only one batch
    // at a time, then waits for the coordinator's acknowledgement (backpressure).
    function sieveBatch(iterator) {
      const batch = { relations: [], polynomials: 0, polynomialGroups: 0, candidates: 0, done: false };
      const deadline = performance.now() + 20;
      do {
        const step = iterator.next();
        if (step.done) { batch.done = true; break; }
        batch.relations.push(...step.value.relations);
        batch.polynomials += step.value.polynomials;
        batch.polynomialGroups += step.value.polynomialGroups;
        batch.candidates += step.value.candidates;
        if (step.value.done) { batch.done = true; break; }
      } while (batch.polynomials < 16 && performance.now() < deadline);
      return batch;
    }

    function siqs(n, report, stats) {
      const state = createSiqs(n, report, stats);
      if (state.divisor) return state.divisor;
      const sieve = createSieve(state.context);
      for (;;) {
        for (const batch of sieve.family(state.nextFamily())) {
          const divisor = state.accept(batch);
          if (divisor !== null) return divisor;
        }
      }
    }

    function parallelSiqs(n, report, stats, ports, searchId, overrides = {}) {
      let nextRequestId = 0;
      const requests = new Map();
      const matrixWorkers = ports.map((port, index) => ({
        request(message, transfers) {
          return new Promise((resolve, reject) => {
            const requestId = ++nextRequestId;
            requests.set(requestId, { resolve, reject, index, matrixId: message.matrixId });
            try { port.postMessage({ ...message, requestId }, transfers); }
            catch (error) { requests.delete(requestId); reject(error); }
          });
        }
      }));
      const state = createSiqs(n, report, stats, ports.length + 1, overrides, matrixWorkers);
      if (state.divisor) return Promise.resolve(state.divisor);
      const sieve = createSieve(state.context);
      return new Promise((resolve, reject) => {
        let finished = false, localJob = null, timer = null, processing = false;
        const pending = [];
        // Alternate message and timer tasks: a self-posting MessagePort alone
        // can starve other ports in some runtimes. Starting the timer from a
        // message task also avoids nested-timer clamping in the browser.
        const scheduler = new MessageChannel();
        function finish(error, divisor) {
          if (finished) return;
          finished = true;
          clearTimeout(timer);
          scheduler.port1.close();
          scheduler.port2.close();
          for (const port of ports) port.onmessage = null;
          for (const request of requests.values()) request.reject(error || new Error("Matrix search ended."));
          requests.clear();
          pending.length = 0;
          if (error) reject(error); else resolve(divisor);
        }
        // Matrix products await helper replies. Serialize relation collection
        // while awaiting them, leaving at most one unacknowledged sieve batch
        // per helper. Matrix replies bypass this queue to avoid a deadlock.
        async function consume() {
          if (processing || finished) return;
          processing = true;
          try {
            while (pending.length && !finished) {
              const { batch, worker } = pending.shift();
              const divisor = await state.accept(batch, worker);
              if (finished) return;
              if (divisor !== null) { finish(null, divisor); return; }
              if (worker === 0) {
                if (batch.done) localJob = null;
                scheduler.port2.postMessage(null);
              } else {
                ports[worker - 1].postMessage(batch.done ?
                  { type: "family", searchId, selected: state.nextFamily() } :
                  { type: "continue", searchId });
              }
            }
          } catch (error) { finish(error); }
          finally { processing = false; }
        }
        function enqueue(batch, worker) {
          pending.push({ batch, worker });
          consume();
        }
        function pump() {
          if (finished) return;
          try {
            if (!localJob) localJob = sieve.family(state.nextFamily());
            const batch = sieveBatch(localJob);
            enqueue(batch, 0);
          } catch (error) { finish(error); }
        }
        scheduler.port1.onmessage = () => { timer = setTimeout(pump, 0); };
        try {
          ports.forEach((port, index) => {
            port.onmessage = ({ data }) => {
              if (data.type === "matrix-result" || data.type === "matrix-error") {
                const request = requests.get(data.requestId);
                if (!request || request.index !== index || request.matrixId !== data.matrixId) return;
                requests.delete(data.requestId);
                if (data.type === "matrix-error") request.reject(new Error(data.message));
                else request.resolve(data.result);
                return;
              }
              // A cofactor may finish with helper work still in flight.
              if (finished || data.searchId !== searchId) return;
              try {
                if (data.type === "error") throw new Error(data.message);
                enqueue(data.batch, index + 1);
              } catch (error) { finish(error); }
            };
            port.postMessage({ type: "init", searchId, context: state.context });
            port.postMessage({ type: "family", searchId, selected: state.nextFamily() });
          });
          scheduler.port2.postMessage(null);
        } catch (error) { finish(error); }
      });
    }

    function mpqs(n, report, stats) {
      return n.toString().length <= 25 ? tinyMpqs(n, report, stats) : siqs(n, report, stats);
    }

    function* factorSteps(value, onProgress, threads = 1, options = {}) {
      const validation = validateInput(value);
      if (!validation.valid) throw new Error(validation.message);
      const original = BigInt(validation.text);
      const stats = { polynomials: 0, relations: 0, ...emptyRelationCounts(), polynomialGroups: 0,
        reusedPolynomials: 0, searches: [],
        ecm: { enabled: options.ecmPretesting === true, curvesStarted: 0, curvesCompleted: 0,
          elapsedMs: 0, factorsFound: 0, searches: [] },
        squfof: { iterations: 0, forwardIterations: 0, reverseIterations: 0, multipliers: 0,
          failedMultipliers: 0, elapsedMs: 0, factorsFound: 0, searches: [] },
        rho: { budgetMs: 0, elapsedMs: 0, evaluations: 0, attempts: 0, factorsFound: 0, searches: [] } };
      if (original < 4n) {
        return { kind: "not-composite", message: original < 2n ?
          "0 and 1 are neither prime nor composite. Enter a composite number." :
          "This number is prime. Enter a composite number.", factors: [], stats };
      }
      const factors = [];
      onProgress({ stage: "Trial division through 1000", digits: original.toString().length });
      let remaining = original;
      for (const p of trialPrimes) {
        const prime = BigInt(p);
        while (remaining % prime === 0n) {
          factors.push(prime);
          remaining /= prime;
        }
        if (remaining === 1n) break;
      }
      if (factors.length > 0 && remaining > 1n) {
        // Send a snapshot before any expensive work on the remaining cofactor.
        // Decimal strings can be posted directly from the worker to the page.
        onProgress({ stage: "Trial division complete", factors: factors.map(String),
          remaining: remaining.toString(), digits: remaining.toString().length });
      }

      // All small factors are removed once, before primality testing or MPQS.
      // Every recursive divisor of the remainder is also coprime to them.
      function* split(n, ecmCredit = []) {
        if (n === 1n) return;
        onProgress({ stage: "Checking remaining factors", digits: n.toString().length });
        if (isProbablePrime(n)) { factors.push(n); return; }
        // Perfect powers are common in logs and do not need a sieve. Testing
        // prime exponents suffices; recursive splitting handles other powers.
        for (const k of trialPrimes) {
          if (k > n.toString(2).length) break;
          const r = integerRoot(n, k);
          if (r <= 1000n) break; // All factors <=1000 have already been removed.
          if (r ** BigInt(k) === n) {
            const start = factors.length;
            yield* split(r, ecmCredit);
            const rootFactors = factors.slice(start);
            for (let repeat = 1; repeat < k; repeat++) factors.push(...rootFactors);
            return;
          }
        }
        let divisor = null, factorSource = "Brent-Rho";
        if (n < SQUFOF_LIMIT) {
          factorSource = "SQUFOF";
          const digits = n.toString().length, record = { digits };
          stats.squfof.searches.push(record);
          Object.assign(stats, emptyRelationCounts());
          const reportSqufof = (state, extra = {}) => onProgress({ stage: "SQUFOF", digits,
            squfofIterations: state.iterations, squfofMultipliers: state.multipliers,
            squfofThreads: state.threads || threads, ...extra });
          reportSqufof({ iterations: 0, multipliers: 0 }, { ...emptyRelationCounts(), resetRelations: true });
          const result = yield { kind: "squfof", n, reportSqufof };
          divisor = result.divisor;
          if (divisor === null || divisor <= 1n || divisor >= n || n % divisor !== 0n) {
            throw new Error("SQUFOF did not return a valid divisor.");
          }
          Object.assign(record, result, { divisor: String(divisor) });
          for (const key of ["iterations", "forwardIterations", "reverseIterations", "multipliers",
            "failedMultipliers", "elapsedMs"]) stats.squfof[key] += result[key];
          stats.squfof.factorsFound++;
          reportSqufof(result);
        }
        const plan = brentRhoBudget(n, threads);
        if (divisor === null && plan.budgetMs > 0) {
          if (!stats.rho.searches.length) stats.rho.budgetMs = plan.budgetMs;
          const budgetMs = Math.min(plan.budgetMs, Math.max(0, stats.rho.budgetMs - stats.rho.elapsedMs));
          const maxEvaluations = RHO_MAX_EVALUATIONS - stats.rho.evaluations;
          if (budgetMs > 0 && maxEvaluations > 0) {
            const record = { digits: plan.digits, estimatedTotalMs: plan.estimatedTotalMs, budgetMs };
            stats.rho.searches.push(record);
            Object.assign(stats, emptyRelationCounts());
            const reportRho = (state, extra = {}) => onProgress({ stage: "Brent-Rho", digits: plan.digits,
              rhoAttempt: state.attempts, rhoEvaluations: state.evaluations,
              rhoElapsedMs: state.elapsedMs, rhoBudgetMs: budgetMs, rhoThreads: state.threads || threads, ...extra });
            reportRho({ attempts: 0, evaluations: 0, elapsedMs: 0 },
              { ...emptyRelationCounts(), resetRelations: true });
            const result = yield { kind: "rho", n, budgetMs, maxEvaluations, reportRho };
            Object.assign(record, result, { divisor: result.divisor === null ? null : String(result.divisor) });
            stats.rho.elapsedMs += result.elapsedMs;
            stats.rho.evaluations += result.evaluations;
            stats.rho.attempts += result.attempts;
            divisor = result.divisor;
            if (divisor !== null) stats.rho.factorsFound++;
            reportRho(result);
          }
        }
        if (divisor === null && stats.ecm.enabled) {
          const ecmPlan = ecm.normalPlan(n, ecmCredit);
          if (ecmPlan.totalCurves > 0) {
            factorSource = "ECM";
            Object.assign(stats, emptyRelationCounts());
            const reportECM = (state, extra = {}) => onProgress({
              stage: state.stage || "ECM pretesting", digits: ecmPlan.digits,
              ecmThreads: threads, ecmCurves: state.curvesCompleted,
              ecmTotalCurves: ecmPlan.totalCurves, ecmDepth: state.depth,
              ecmTarget: ecmPlan.targetDepth, ecmB1: state.B1, ecmB2: state.B2, ...extra
            });
            reportECM({ curvesCompleted: 0, depth: ecm.depth(ecmCredit) },
              { ...emptyRelationCounts(), resetRelations: true });
            const result = yield { kind: "ecm", n, plan: ecmPlan, reportECM };
            divisor = result.divisor;
            if (divisor !== null && (divisor <= 1n || divisor >= n || n % divisor !== 0n)) {
              throw new Error("ECM did not return a valid divisor.");
            }
            ecmCredit = result.counts;
            stats.ecm.searches.push({ ...result, digits: ecmPlan.digits,
              divisor: divisor === null ? null : String(divisor) });
            for (const key of ["curvesStarted", "curvesCompleted", "elapsedMs"]) stats.ecm[key] += result[key];
            if (divisor !== null) stats.ecm.factorsFound++;
            reportECM(result);
          }
        }
        if (divisor !== null) {
          // Resolve the smaller branch first. Only verified prime/probable-prime
          // factors go into the partial output; a raw divisor may be composite.
          const other = n / divisor, before = factors.length;
          yield* split(divisor < other ? divisor : other, ecmCredit);
          if (factors.length > before) {
            const known = factors.slice().sort((a, b) => a < b ? -1 : a > b ? 1 : 0);
            const unresolved = original / known.reduce((product, p) => product * p, 1n);
            if (unresolved > 1n) onProgress({ stage: `${factorSource} factor found`, factorSource,
              factors: known.map(String), remaining: String(unresolved), digits: unresolved.toString().length });
          }
          yield* split(divisor < other ? other : divisor, ecmCredit);
          return;
        }
        divisor = yield { n, stats };
        if (divisor <= 1n || divisor >= n || n % divisor !== 0n) {
          throw new Error("MPQS produced an invalid divisor.");
        }
        yield* split(divisor, ecmCredit);
        yield* split(n / divisor, ecmCredit);
      }
      yield* split(remaining);
      factors.sort((a, b) => a < b ? -1 : a > b ? 1 : 0);
      if (factors.reduce((product, p) => product * p, 1n) !== original) {
        throw new Error("The factors did not reconstruct the input.");
      }
      if (factors.length === 1) {
        return { kind: original <= 1000n ? "not-composite" : "probable-prime", message:
          original <= 1000n ? "This number is prime. Enter a composite number." :
            "This number passes the primality test. Enter a composite number.", factors: [], stats };
      }
      return { kind: "factored", factors: factors.map(String), stats };
    }


    function factor(value, onProgress = () => {}, options = {}) {
      const steps = factorSteps(value, onProgress, 1, options);
      let step = steps.next();
      while (!step.done) {
        const job = step.value;
        step = steps.next(job.kind === "rho" ?
          brentRho(job.n, job.budgetMs, { maxEvaluations: job.maxEvaluations, onProgress: job.reportRho }) :
          job.kind === "squfof" ? squfof(job.n, job.reportSqufof) :
          job.kind === "ecm" ? ecm.pretest(job.n, job.plan, job.reportECM) :
          mpqs(job.n, onProgress, job.stats));
      }
      return step.value;
    }

    async function factorParallel(value, onProgress, ports = [], options = {}) {
      if (!ports.length) return factor(value, onProgress, options);
      const steps = factorSteps(value, onProgress, ports.length + 1, options);
      let step = steps.next(), searchId = 0;
      while (!step.done) {
        const job = step.value;
        if (job.kind === "rho") {
          step = steps.next(await brentRhoParallel(job.n, job.budgetMs, ports, {
            maxEvaluations: job.maxEvaluations, onProgress: job.reportRho
          }));
        } else if (job.kind === "squfof") {
          step = steps.next(await squfofParallel(job.n, job.reportSqufof, ports));
        } else if (job.kind === "ecm") {
          step = steps.next(await ecm.pretestParallel(job.n, job.plan, ports, job.reportECM));
        } else {
          const divisor = job.n.toString().length <= 25 ? tinyMpqs(job.n, onProgress, job.stats) :
            await parallelSiqs(job.n, onProgress, job.stats, ports, ++searchId);
          step = steps.next(divisor);
        }
      }
      return step.value;
    }

    return { MAX_DIGITS, MR_ROUNDS, RHO_MIN_DIGITS, SQUFOF_LIMIT, ecm, createSqufofJob, squfof, squfofParallel,
      brentRhoBudget, brentAttempt, brentRho, brentRhoParallel, createRhoJob,
      validateInput, factor, factorParallel, createSieve, sieveBatch, isProbablePrime,
      squfofSmall, splitLargePrimes, createLargePrimeGraph, createSparseSolver, packExponents, chooseMultiplier, createSiqs, parallelSiqs,
      gcd, powMod, inverse, isqrt, integerRoot, sqrtMod, randomBelow };
  }

  function workerSource() {
    return `"use strict";
      const createMatrixAlgebra = ${createMatrixAlgebra.toString()};
      const createECM = ${createECM.toString()};
      const engine = (${createEngine.toString()})();
      self.onmessage = async function (event) {
        if (event.data.type === "connect") {
          const port = event.data.port;
          let sieve, job, searchId;
          const prepassJobs = new Map();
          const matrixWorker = createMatrixAlgebra().createMatrixWorker();
          port.onmessage = ({ data }) => {
            if (data.type === "rho" || data.type === "squfof" || data.type === "ecm") {
              try {
                let result = null;
                if (data.operation === "init") {
                  const job = data.type === "rho" ? engine.createRhoJob(data.n, data.budgetMs, data.options) :
                    data.type === "ecm" ? engine.ecm.createJob(data.n, data.options) :
                    engine.createSqufofJob(data.n, data.options);
                  prepassJobs.set(data.type, { job, id: data.jobId });
                } else {
                  const entry = prepassJobs.get(data.type);
                  if (!entry || data.jobId !== entry.id) throw new Error("Stale factoring worker request.");
                  if (data.operation === "step") result = entry.job.advance();
                  else if (data.operation === "release") prepassJobs.delete(data.type);
                  else throw new Error("Unknown factoring worker operation.");
                }
                port.postMessage({ type: data.type + "-result", jobId: data.jobId, requestId: data.requestId, result });
              } catch (error) {
                port.postMessage({ type: data.type + "-error", jobId: data.jobId, requestId: data.requestId, message: error.message });
              }
              return;
            }
            if (data.type === "matrix") {
              try {
                const result = matrixWorker(data);
                port.postMessage({ type: "matrix-result", matrixId: data.matrixId, requestId: data.requestId, result },
                  result ? [result.buffer] : []);
              } catch (error) {
                port.postMessage({ type: "matrix-error", matrixId: data.matrixId, requestId: data.requestId,
                  message: error.message });
              }
              return;
            }
            try {
              if (data.type === "init") {
                searchId = data.searchId;
                sieve = engine.createSieve(data.context);
                job = null;
                return;
              }
              if (data.searchId !== searchId) return;
              if (data.type === "family") job = sieve.family(data.selected);
              const batch = engine.sieveBatch(job);
              port.postMessage({ type: "batch", searchId, batch });
            } catch (error) {
              port.postMessage({ type: "error", searchId, message: error.message });
            }
          };
          return;
        }
        let lastProgress = 0;
        try {
          const result = await engine.factorParallel(event.data.value, function (progress) {
            const now = Date.now();
            if (now - lastProgress >= 100 || progress.stage !== "Sieving polynomials") {
              self.postMessage({ type: "progress", progress });
              lastProgress = now;
            }
          }, event.data.ports || [], { ecmPretesting: event.data.ecmPretesting === true });
          self.postMessage({ type: "result", result });
        } catch (error) {
          self.postMessage({ type: "error", message: error.message });
        }
      };`;
  }

  const MAX_THREADS = 32;

  function validateThreads(value) {
    const text = String(value).trim();
    return /^[0-9]+$/.test(text) && Number(text) >= 1 && Number(text) <= MAX_THREADS;
  }

  // The page owns every worker, including the coordinator, so cancellation
  // never leaves detached helper workers running. Ports keep relation traffic
  // off the UI thread. Injectable factories also exercise this exact protocol
  // in Node tests/benchmarks without a browser or external dependencies.
  function startFactorization(value, threads, onMessage, options = {}) {
    if (!validateThreads(threads)) throw new Error(`Threads must be a whole number from 1 to ${MAX_THREADS}.`);
    threads = Number(threads);
    const workers = [], channels = [];
    let workerURL = null, stopped = false;
    function terminate() {
      if (stopped) return;
      stopped = true;
      for (const worker of workers) worker.terminate();
      for (const channel of channels) {
        channel.port1.close();
        channel.port2.close();
      }
      if (workerURL) URL.revokeObjectURL(workerURL);
    }
    try {
      const source = workerSource();
      let makeWorker = options.createWorker;
      if (!makeWorker) {
        workerURL = URL.createObjectURL(new Blob([source], { type: "text/javascript" }));
        makeWorker = () => new Worker(workerURL);
      }
      const makeChannel = options.createChannel || (() => new MessageChannel());
      for (let index = 0; index < threads; index++) {
        const worker = makeWorker(source);
        workers.push(worker);
        worker.onerror = event => {
          if (stopped) return;
          if (event.preventDefault) event.preventDefault();
          terminate();
          onMessage({ type: "error", message: event.message || "A factoring worker stopped unexpectedly." });
        };
      }
      workers[0].onmessage = ({ data }) => {
        if (stopped) return;
        if (data.type === "result" || data.type === "error") terminate();
        onMessage(data);
      };
      const ports = [];
      for (let index = 1; index < threads; index++) {
        const channel = makeChannel();
        channels.push(channel);
        workers[index].postMessage({ type: "connect", port: channel.port2 }, [channel.port2]);
        ports.push(channel.port1);
      }
      workers[0].postMessage({ value, ports, ecmPretesting: options.ecmPretesting === true }, ports);
      return { terminate };
    } catch (error) {
      terminate();
      throw error;
    }
  }

  const api = { createEngine, workerSource, MAX_THREADS, validateThreads, startFactorization };
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.MPQS = api;
})(globalThis);
