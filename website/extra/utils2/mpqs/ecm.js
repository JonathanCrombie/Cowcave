/* Plain JavaScript two-stage ECM. Scheduling and Suyama parameterization
 * follow Ben Buhrow's public-domain yafu/factor/autofactor.c and
 * yafu/factor/avx-ecm/avxecm.c. No native arithmetic or external libraries. */
(function (root) {
  "use strict";
  function createECM() {
    "use strict";
    const bounds = [2000, 11000, 50000, 250000, 1000000];
    // YAFU's avx_ecm_data: Suyama curves, power-1 stage 2, B2 = 100*B1.
    // Rows target t15,t20,...,t40. Five bounds cover the <=110-digit UI.
    const expected = [
      [32, 8, 4, 3, 2],
      [1203, 98, 26, 11, 7],
      [88933, 2174, 281, 74, 32],
      [1.1e7, 73653, 4455, 671, 206],
      [2e9, 3495192, 95330, 7965, 1676],
      [9e99, 2.2e8, 2610023, 117616, 16589]
    ];
    const mod = (a, n) => (a % n + n) % n;
    function gcd(a, b) {
      while (b) [a, b] = [b, a % b];
      return a < 0n ? -a : a;
    }
    function inverse(a, n) {
      let b = n, x = 1n, y = 0n;
      a = mod(a, n);
      while (b) {
        const q = a / b;
        [a, b] = [b, a - q * b];
        [x, y] = [y, x - q * y];
      }
      return a === 1n ? mod(x, n) : null;
    }
    function depth(counts) {
      for (let row = 0; row < expected.length; row++) {
        const work = expected[row].reduce((sum, value, i) => sum + (counts[i] || 0) / value, 0);
        if (work < 1) return row === 0 ? 0 : 10 + 5 * row + 5 * work;
      }
      return 40;
    }
    function normalPlan(n, credit = []) {
      const digits = n.toString().length, targetDepth = 4 * digits / 13;
      if (digits > 110 || n < 2n) throw new Error("Invalid ECM input size.");
      const initialCounts = bounds.map((_, i) => credit[i] || 0);
      if (initialCounts.some(v => !Number.isSafeInteger(v) || v < 0)) throw new Error("Invalid ECM work credit.");
      const counts = initialCounts.slice(), levels = [];
      if (digits > 45) {
        for (let i = 0; i < bounds.length && depth(counts) < targetDepth; i++) {
          const cap = expected[i][i];
          let lo = 0, hi = Math.max(0, cap - counts[i]);
          while (lo < hi) {
            const mid = Math.floor((lo + hi) / 2), trial = counts.slice();
            trial[i] += mid;
            if (depth(trial) >= targetDepth) hi = mid;
            else lo = mid + 1;
          }
          if (lo) levels.push({ index: i, B1: bounds[i], B2: 100 * bounds[i], curves: lo });
          counts[i] += lo;
        }
      }
      return { digits, targetDepth, initialCounts, levels,
        totalCurves: levels.reduce((sum, level) => sum + level.curves, 0) };
    }

    // Montgomery x:z arithmetic. All coordinates and products remain BigInts.
    function double(P, a24, n) {
      const plus = (P.x + P.z) % n, minus = mod(P.x - P.z, n);
      const aa = plus * plus % n, bb = minus * minus % n, c = mod(aa - bb, n);
      return { x: aa * bb % n, z: c * ((bb + a24 * c) % n) % n };
    }
    function add(P, Q, difference, n) {
      const a = (P.x + P.z) * mod(Q.x - Q.z, n) % n;
      const b = mod(P.x - P.z, n) * ((Q.x + Q.z) % n) % n;
      const plus = (a + b) % n, minus = mod(a - b, n);
      return { x: difference.z * (plus * plus % n) % n,
        z: difference.x * (minus * minus % n) % n };
    }
    function multiply(P, scalar, a24, n) {
      if (scalar === 0) return { x: 1n, z: 0n };
      let lo = P, hi = double(P, a24, n);
      const bits = scalar.toString(2);
      for (let i = 1; i < bits.length; i++) {
        if (bits[i] === "0") { hi = add(lo, hi, P, n); lo = double(lo, a24, n); }
        else { lo = add(lo, hi, P, n); hi = double(hi, a24, n); }
      }
      return lo;
    }
    function makeCurve(n, sigma) {
      const u = mod(sigma * sigma - 5n, n), v = 4n * sigma % n;
      const x = u * u % n * u % n, z = v * v % n * v % n;
      const denominator = 16n * x % n * v % n, g = gcd(denominator, n);
      if (g !== 1n) return { divisor: g < n ? g : null, singular: true };
      const d = mod(v - u, n);
      const a24 = d * d % n * d % n * ((3n * u + v) % n) % n *
        inverse(denominator, n) % n;
      const discriminant = gcd(a24 * mod(a24 - 1n, n) % n, n);
      if (discriminant !== 1n) return { divisor: discriminant < n ? discriminant : null, singular: true };
      return { P: { x, z }, a24, divisor: null, singular: false };
    }
    function batchGcd(terms, n) {
      let product = 1n;
      for (const term of terms) product = product * term % n;
      const g = gcd(product, n);
      if (g > 1n && g < n) return { divisor: g, saturated: false };
      if (g === n) {
        // Different factors can divide different terms in the same batch.
        for (const term of terms) {
          const individual = gcd(term, n);
          if (individual > 1n && individual < n) return { divisor: individual, saturated: false };
        }
      }
      return { divisor: null, saturated: g === n };
    }
    function smallPrimes(limit) {
      const composite = new Uint8Array(limit + 1), primes = [];
      for (let p = 2; p <= limit; p++) {
        if (composite[p]) continue;
        primes.push(p);
        for (let j = p * p; j <= limit; j += p) composite[j] = 1;
      }
      return primes;
    }
    let cached = null;
    function* prepare(B1, B2) {
      if (cached?.B1 === B1 && cached.B2 === B2) return cached;
      // Pair q=kD +/- r. Both primes need just one projective x difference.
      const D = B1 >= 210 ? 210 : 30;
      const residues = [], slots = new Int8Array(D / 2 + 1).fill(-1);
      for (let r = 1; r < D / 2; r += 2) {
        if (gcd(BigInt(r), BigInt(D)) === 1n) { slots[r] = residues.length; residues.push(r); }
      }
      const first = Math.max(1, Math.floor((B1 + 1) / D + 0.5));
      const last = Math.floor(B2 / D + 0.5);
      const masks = new Uint32Array(Math.max(0, last - first + 1)), powers = [];
      const bases = smallPrimes(Math.floor(Math.sqrt(B2)));
      for (let low = 2; low <= B2; low += 32768) {
        const high = Math.min(B2, low + 32767), composite = new Uint8Array(high - low + 1);
        for (const p of bases) {
          if (p * p > high) break;
          for (let j = Math.max(p * p, Math.ceil(low / p) * p); j <= high; j += p) composite[j - low] = 1;
        }
        for (let p = low; p <= high; p++) {
          if (composite[p - low]) continue;
          if (p <= B1) {
            let power = p;
            while (power * p <= B1) power *= p;
            powers.push(power);
          } else {
            const k = Math.floor(p / D + 0.5), r = Math.abs(p - k * D);
            if (k >= first && slots[r] >= 0) masks[k - first] |= 1 << slots[r];
          }
        }
        yield "Preparing ECM";
      }
      cached = { B1, B2, D, first, masks, residues, powers };
      return cached;
    }
    function* curve(n, B1, B2, sigma) {
      const setup = makeCurve(n, sigma);
      if (setup.singular) return { divisor: setup.divisor, foundStage: 0 };
      const table = yield* prepare(B1, B2);
      const { a24 } = setup;
      let P = setup.P, terms = [];
      for (const power of table.powers) {
        P = multiply(P, power, a24, n);
        terms.push(P.z);
        if (terms.length >= 32 || P.z === 0n) {
          const result = batchGcd(terms, n); terms = [];
          if (result.divisor || result.saturated) return { divisor: result.divisor, foundStage: 1 };
        }
        yield "ECM stage 1";
      }
      if (terms.length) {
        const result = batchGcd(terms, n);
        if (result.divisor || result.saturated) return { divisor: result.divisor, foundStage: 1 };
      }
      terms = [];
      const babies = [];
      for (const r of table.residues) {
        babies.push(multiply(P, r, a24, n));
        yield "ECM stage 2";
      }
      const step = multiply(P, table.D, a24, n);
      let previous = multiply(step, table.first - 1, a24, n);
      let giant = multiply(step, table.first, a24, n);
      for (let index = 0; index < table.masks.length; index++) {
        let mask = table.masks[index];
        while (mask) {
          const bit = mask & -mask, baby = babies[31 - Math.clz32(bit)];
          terms.push(mod(giant.x * baby.z - baby.x * giant.z, n));
          mask ^= bit;
          if (terms.length >= 128) {
            const result = batchGcd(terms, n); terms = [];
            if (result.divisor || result.saturated) return { divisor: result.divisor, foundStage: 2 };
          }
        }
        const next = index + table.first === 1 ? double(giant, a24, n) : add(giant, step, previous, n);
        previous = giant; giant = next;
        yield "ECM stage 2";
      }
      return { divisor: terms.length ? batchGcd(terms, n).divisor : null, foundStage: 2 };
    }
    function createJob(n, options) {
      const { B1, B2, indices, firstSigma } = options;
      if (typeof n !== "bigint" || n < 7n || !(n & 1n) ||
          !Number.isSafeInteger(B1) || B1 < 30 || B1 > 1000000 ||
          !Number.isSafeInteger(B2) || B2 < B1 || B2 > 100000000 ||
          !Array.isArray(indices) || indices.some(i => !Number.isSafeInteger(i) || i < 0) ||
          new Set(indices).size !== indices.length || typeof firstSigma !== "bigint" || firstSigma < 6n) {
        throw new Error("Invalid ECM curve job.");
      }
      let current = null, cursor = 0;
      const state = { divisor: null, curvesStarted: 0, curvesCompleted: 0, operations: 0,
        elapsedMs: 0, stage: "Preparing ECM", foundStage: null, done: indices.length === 0 };
      return { advance() {
        if (state.done) return { ...state };
        const started = performance.now();
        for (let ticks = 0; ticks < 128; ticks++) {
          if (!current) {
            current = curve(n, B1, B2, firstSigma + BigInt(indices[cursor++]));
            state.curvesStarted++;
          }
          const step = current.next();
          state.operations++;
          if (step.done) {
            current = null;
            if (step.value.divisor !== null) {
              state.divisor = step.value.divisor;
              state.foundStage = step.value.foundStage;
              state.done = true;
            } else state.curvesCompleted++;
            if (cursor === indices.length) state.done = true;
          } else state.stage = step.value;
          if (state.done || performance.now() - started >= 4) break;
        }
        state.elapsedMs += Math.max(0, performance.now() - started);
        return { ...state };
      } };
    }
    function randomSigma() {
      const words = new Uint32Array(2);
      globalThis.crypto.getRandomValues(words);
      return 6n + (BigInt(words[0]) << 32n) + BigInt(words[1]);
    }
    function searchState(plan, threads) {
      return { divisor: null, threads, targetDepth: plan.targetDepth, totalCurves: plan.totalCurves,
        curvesStarted: 0, curvesCompleted: 0, operations: 0, elapsedMs: 0,
        counts: plan.initialCounts.slice(), depth: depth(plan.initialCounts), foundStage: null };
    }
    function combine(states, base, level) {
      const result = { ...base, counts: base.counts.slice(), B1: level.B1, B2: level.B2 };
      for (const key of ["curvesStarted", "curvesCompleted", "operations", "elapsedMs"]) {
        result[key] += states.reduce((sum, state) => sum + state[key], 0);
      }
      const found = states.find(state => state.divisor !== null);
      result.divisor = found?.divisor ?? null;
      result.foundStage = found?.foundStage ?? null;
      result.counts[level.index] += states.reduce((sum, state) => sum + state.curvesCompleted, 0);
      result.depth = depth(result.counts);
      result.stage = states.find(state => !state.done)?.stage || "ECM pretesting";
      return result;
    }
    function pretest(n, plan, onProgress = () => {}, firstSigma = randomSigma()) {
      let result = searchState(plan, 1), offset = 0, lastReport = -Infinity;
      for (const level of plan.levels) {
        const base = result, job = createJob(n, { ...level, firstSigma,
          indices: Array.from({ length: level.curves }, (_, i) => offset + i) });
        offset += level.curves;
        let state;
        do {
          state = job.advance();
          result = combine([state], base, level);
          if (state.done || performance.now() - lastReport >= 200) {
            onProgress(result); lastReport = performance.now();
          }
        } while (!state.done);
        if (result.divisor !== null) break;
      }
      return result;
    }
    let nextJobId = 0;
    async function pretestParallel(n, plan, ports, onProgress = () => {}, firstSigma = randomSigma()) {
      if (!ports.length) return pretest(n, plan, onProgress, firstSigma);
      if (ports.length > 31) throw new Error("Invalid ECM thread count.");
      const threads = ports.length + 1;
      let result = searchState(plan, threads), offset = 0;
      for (const level of plan.levels) {
        const base = result, jobId = ++nextJobId, previous = ports.map(port => port.onmessage);
        const assignments = Array.from({ length: threads }, (_, i) =>
          Array.from({ length: level.curves }, (_, id) => offset + id).filter((_, id) => id % threads === i));
        offset += level.curves;
        const states = assignments.map(indices => ({ divisor: null, curvesStarted: 0,
          curvesCompleted: 0, operations: 0, elapsedMs: 0, stage: "Preparing ECM",
          foundStage: null, done: !indices.length }));
        const pending = new Map();
        let requestId = 0;
        function request(index, operation, extra = {}) {
          return new Promise((resolve, reject) => {
            const id = ++requestId;
            pending.set(id, { index, resolve, reject });
            try { ports[index].postMessage({ type: "ecm", jobId, requestId: id, operation, ...extra }); }
            catch (error) { pending.delete(id); reject(error); }
          });
        }
        function accept(index, state) {
          const old = states[index];
          if (!state || ["curvesStarted", "curvesCompleted", "operations"].some(key =>
              !Number.isSafeInteger(state[key]) || state[key] < old[key]) ||
              state.curvesCompleted > state.curvesStarted || state.curvesStarted > assignments[index].length ||
              !Number.isFinite(state.elapsedMs) || state.elapsedMs < old.elapsedMs ||
              typeof state.done !== "boolean" ||
              (state.done && state.divisor === null && state.curvesCompleted !== assignments[index].length) ||
              (state.divisor !== null && !state.done) ||
              (state.divisor !== null && (typeof state.divisor !== "bigint" ||
                state.divisor <= 1n || state.divisor >= n || n % state.divisor !== 0n))) {
            throw new Error("Invalid ECM worker result.");
          }
          states[index] = state;
        }
        try {
          ports.forEach((port, index) => {
            port.onmessage = ({ data }) => {
              if (data.jobId !== jobId || (data.type !== "ecm-result" && data.type !== "ecm-error")) return;
              const call = pending.get(data.requestId);
              if (!call || call.index !== index) return;
              pending.delete(data.requestId);
              if (data.type === "ecm-error") call.reject(new Error(data.message));
              else call.resolve(data.result);
            };
          });
          const local = createJob(n, { ...level, firstSigma, indices: assignments[0] });
          await Promise.all(ports.map((_, i) => request(i, "init", { n,
            options: { ...level, firstSigma, indices: assignments[i + 1] } })));
          let lastReport = -Infinity;
          for (;;) {
            const replies = ports.map((_, i) => states[i + 1].done ? Promise.resolve() :
              request(i, "step").then(state => accept(i + 1, state)));
            replies.push(Promise.resolve().then(() => { if (!states[0].done) accept(0, local.advance()); }));
            // Drain this wave before changing methods or releasing helpers.
            const outcomes = await Promise.allSettled(replies);
            const failure = outcomes.find(outcome => outcome.status === "rejected");
            if (failure) throw failure.reason;
            result = combine(states, base, level);
            const done = result.divisor !== null || states.every(state => state.done);
            if (done || performance.now() - lastReport >= 200) {
              onProgress(result); lastReport = performance.now();
            }
            if (done) break;
          }
        } finally {
          await Promise.allSettled(ports.map((_, i) => request(i, "release")));
          ports.forEach((port, i) => { port.onmessage = previous[i]; });
        }
        if (result.divisor !== null) break;
      }
      return result;
    }
    return { normalPlan, depth, double, add, multiply, makeCurve, batchGcd,
      createJob, pretest, pretestParallel };
  }
  const api = { createECM };
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.MPQSECM = api;
})(globalThis);
