/* GF(2) matrix filtering and Montgomery block Lanczos, guided by the
 * public-domain common/lanczos code by Jason Papadopoulos in msieve.
 * A matrix column represents a relation; its rows are odd prime exponents.
 * 32-way blocks use JavaScript's exact bitwise operations, not floating-point
 * numerical linear algebra. Keep this factory self-contained for Blob workers.
 */
(function (root) {
  "use strict";
  function createMatrixAlgebra() {
    "use strict";
    const xorLists = (a, b) => {
      const out = [];
      let i = 0, j = 0;
      while (i < a.length && j < b.length) {
        if (a[i] < b[j]) out.push(a[i++]);
        else if (b[j] < a[i]) out.push(b[j++]);
        else { i++; j++; }
      }
      while (i < a.length) out.push(a[i++]);
      while (j < b.length) out.push(b[j++]);
      return Uint32Array.from(out);
    };
    function contains(sorted, value) {
      let lo = 0, hi = sorted.length;
      while (lo < hi) {
        const mid = (lo + hi) >>> 1;
        if (sorted[mid] < value) lo = mid + 1; else hi = mid;
      }
      return sorted[lo] === value;
    }
    function checkColumns(rowCount, columns) {
      if (!Number.isSafeInteger(rowCount) || rowCount < 0) throw new Error("Invalid matrix row count.");
      for (const column of columns) {
        let previous = -1;
        for (const r of column) {
          if (!Number.isInteger(r) || r <= previous || r >= rowCount) {
            throw new Error("Matrix columns must have sorted, distinct row indices.");
          }
          previous = r;
        }
      }
    }
    function filterMatrix(rowCount, input, options = {}) {
      checkColumns(rowCount, input);
      const excess = options.excess ?? 64, maxMergeWeight = options.maxMergeWeight ?? 1000;
      const columns = input.map(c => Uint32Array.from(c));
      const histories = input.map((_, i) => [i]);
      const counts = new Uint32Array(rowCount), ids = new Uint32Array(rowCount);
      const last = new Int32Array(rowCount), queue = [], dependencies = [];
      const stats = { originalRows: rowCount, originalColumns: input.length,
        singletonColumns: 0, mergedColumns: 0, heavyColumns: 0, passes: 0 };
      let live = columns.length, nonempty = 0;
      for (let c = 0; c < columns.length; c++) {
        for (const r of columns[c]) { if (!counts[r]++) nonempty++; ids[r] ^= c; }
      }
      function remove(c) {
        const column = columns[c];
        if (column === null) return;
        columns[c] = null; live--;
        for (const r of column) {
          ids[r] ^= c;
          if (--counts[r] === 0) nonempty--;
          else if (counts[r] === 1) queue.push(r);
        }
      }
      function singletons() {
        while (queue.length) {
          const r = queue.pop();
          if (counts[r] !== 1) continue;
          const c = ids[r];
          remove(c); histories[c] = null; stats.singletonColumns++;
        }
      }
      for (let r = 0; r < rowCount; r++) if (counts[r] === 1) queue.push(r);
      for (let c = 0; c < columns.length; c++) {
        if (columns[c].length === 0) {
          dependencies.push(Uint32Array.from(histories[c]));
          remove(c); histories[c] = null;
        }
      }
      for (;;) {
        stats.passes++;
        singletons();
        last.fill(-1);
        for (let c = 0; c < columns.length; c++) {
          if (columns[c]) for (const r of columns[c]) last[r] = c;
        }
        let merged = 0;
        for (let c = 0; c < columns.length; c++) {
          const col = columns[c];
          if (col === null) continue;
          for (const r of col) {
            if (counts[r] !== 2) continue;
            const other = last[r];
            // A previous merge may have invalidated the cached last column.
            if (other < 0 || other === c || columns[other] === null ||
                !contains(columns[other], r) ||
                col.length + columns[other].length >= maxMergeWeight) continue;
            const combined = xorLists(col, columns[other]);
            const history = histories[c].concat(histories[other]);
            remove(c); remove(other);
            histories[c] = null;
            if (!combined.length) {
              dependencies.push(Uint32Array.from(history));
              histories[other] = null;
            } else {
              columns[other] = combined; histories[other] = history; live++;
              for (const row of combined) {
                if (!counts[row]++) nonempty++;
                ids[row] ^= other; last[row] = other;
              }
            }
            merged++; stats.mergedColumns++;
            singletons();
            break;
          }
        }
        if (merged) continue;
        if (live <= nonempty + excess) break;
        // Retain a surplus of light columns; removing heavy columns may
        // expose new singletons and weight-two cliques, so filter again.
        const heavy = [];
        for (let c = 0; c < columns.length; c++) if (columns[c]) heavy.push(c);
        heavy.sort((a, b) => columns[b].length - columns[a].length || b - a);
        const toRemove = live - nonempty - excess;
        for (let j = 0; j < toRemove; j++) {
          remove(heavy[j]); histories[heavy[j]] = null; stats.heavyColumns++;
        }
      }
      const rowOrder = [];
      for (let r = 0; r < rowCount; r++) if (counts[r]) rowOrder.push(r);
      rowOrder.sort((a, b) => counts[b] - counts[a] || a - b);
      const rowMap = new Uint32Array(rowCount);
      rowOrder.forEach((r, i) => { rowMap[r] = i; });
      const reduced = [], groups = [];
      for (let c = 0; c < columns.length; c++) {
        if (!columns[c]) continue;
        reduced.push(Uint32Array.from(columns[c], r => rowMap[r]).sort());
        groups.push(Uint32Array.from(histories[c]));
      }
      Object.assign(stats, { rows: rowOrder.length, columns: reduced.length,
        nonzeros: reduced.reduce((sum, c) => sum + c.length, 0),
        emptyRows: rowCount - rowOrder.length, immediateDependencies: dependencies.length });
      return { rowCount: rowOrder.length, columns: reduced, groups, dependencies, stats };
    }

    // A pair of sparse layouts gives gather-only B*x and transpose(B)*x.
    // The densest rows can be omitted during Lanczos and restored for the
    // final kernel calculation, following msieve's post-Lanczos treatment.
    function packMatrix(rowCount, columns, postRows) {
      const n = columns.length, offsets = new Uint32Array(n + 1);
      const rowOffsets = new Uint32Array(rowCount + 1);
      for (let c = 0; c < n; c++) {
        offsets[c + 1] = offsets[c] + columns[c].length;
        for (const r of columns[c]) rowOffsets[r + 1]++;
      }
      for (let r = 0; r < rowCount; r++) rowOffsets[r + 1] += rowOffsets[r];
      const entries = new Uint32Array(offsets[n]), byRow = new Uint32Array(offsets[n]);
      const cursor = rowOffsets.slice(), starts = new Uint32Array(n);
      for (let c = 0; c < n; c++) {
        entries.set(columns[c], offsets[c]);
        let start = offsets[c];
        for (const r of columns[c]) {
          byRow[cursor[r]++] = c;
          if (r < postRows) start++;
        }
        starts[c] = start;
      }
      function multiply(x, out, full = false) {
        out.fill(0);
        for (let r = full ? 0 : postRows; r < rowCount; r++) {
          let value = 0;
          for (let j = rowOffsets[r]; j < rowOffsets[r + 1]; j++) value ^= x[byRow[j]];
          out[r] = value;
        }
      }
      function transpose(x, out) {
        for (let c = 0; c < n; c++) {
          let value = 0;
          for (let j = starts[c]; j < offsets[c + 1]; j++) value ^= x[entries[j]];
          out[c] = value;
        }
      }
      return { n, rowCount, activeRows: rowCount - postRows, multiply, transpose,
        postRows, offsets, starts, entries, rowOffsets, byRow };
    }

    // Give each worker contiguous row and column ranges, balanced by nonzeros.
    // Each worker stores only its two sparse slices, not a full matrix copy.
    function partitionMatrix(matrix, count) {
      const { n, rowCount, postRows, offsets, starts, entries, rowOffsets, byRow } = matrix;
      function cuts(prefix, first, last) {
        const result = [first], total = prefix[last] - prefix[first];
        for (let i = 1; i < count; i++) {
          const target = prefix[first] + total * i / count;
          let lo = result[i - 1], hi = last;
          while (lo < hi) {
            const mid = (lo + hi) >>> 1;
            if (prefix[mid] < target) lo = mid + 1; else hi = mid;
          }
          result.push(lo);
        }
        result.push(last);
        return result;
      }
      const activeOffsets = new Uint32Array(n + 1);
      for (let c = 0; c < n; c++) activeOffsets[c + 1] = activeOffsets[c] + offsets[c + 1] - starts[c];
      const rowCuts = cuts(rowOffsets, postRows, rowCount), columnCuts = cuts(activeOffsets, 0, n);
      return Array.from({ length: count }, (_, i) => {
        const rowStart = rowCuts[i], rowEnd = rowCuts[i + 1];
        const columnStart = columnCuts[i], columnEnd = columnCuts[i + 1];
        const rowIndex = rowOffsets.slice(rowStart, rowEnd + 1);
        const rowBase = rowIndex[0];
        for (let j = 0; j < rowIndex.length; j++) rowIndex[j] -= rowBase;
        const columnIndex = activeOffsets.slice(columnStart, columnEnd + 1);
        const columnBase = columnIndex[0];
        for (let j = 0; j < columnIndex.length; j++) columnIndex[j] -= columnBase;
        const columnEntries = new Uint32Array(columnIndex[columnIndex.length - 1]);
        for (let c = columnStart; c < columnEnd; c++) {
          columnEntries.set(entries.subarray(starts[c], offsets[c + 1]), columnIndex[c - columnStart]);
        }
        return { rowStart, columnStart, rowIndex, rowEntries: byRow.slice(rowBase, rowOffsets[rowEnd]),
          columnIndex, columnEntries };
      });
    }
    function partitionProduct(partition, operation, vector) {
      const index = operation === "multiply" ? partition.rowIndex : partition.columnIndex;
      const entries = operation === "multiply" ? partition.rowEntries : partition.columnEntries;
      const result = new Uint32Array(index.length - 1);
      for (let i = 0; i < result.length; i++) {
        let word = 0;
        for (let j = index[i]; j < index[i + 1]; j++) word ^= vector[entries[j]];
        result[i] = word;
      }
      return result;
    }
    // Shared by browser helpers and the real-worker regression tests. Matrix
    // messages are independent of the paused helper's sieve search/job state.
    function createMatrixWorker() {
      let partition = null, matrixId = null;
      return function handle(message) {
        if (message.operation === "init") {
          partition = message.partition; matrixId = message.matrixId;
          return null;
        }
        if (message.matrixId !== matrixId || !partition) throw new Error("Stale matrix worker request.");
        if (message.operation === "release") {
          partition = null; matrixId = null;
          return null;
        }
        if (message.operation !== "multiply" && message.operation !== "transpose") {
          throw new Error("Unknown matrix worker operation.");
        }
        return partitionProduct(partition, message.operation, message.vector);
      };
    }
    let nextMatrixId = 0;
    async function parallelOperator(matrix, workers) {
      const matrixId = ++nextMatrixId, partitions = partitionMatrix(matrix, workers.length + 1);
      const ranges = partitions.map(p => ({ rowStart: p.rowStart, columnStart: p.columnStart,
        rows: p.rowIndex.length - 1, columns: p.columnIndex.length - 1 }));
      const workerProducts = Array(partitions.length).fill(0);
      const request = (worker, operation, extra = {}, transfers = []) =>
        worker.request({ type: "matrix", matrixId, operation, ...extra }, transfers);
      async function close() {
        // Release sparse slices before helpers resume their existing sieve jobs.
        await Promise.allSettled(workers.map(worker => request(worker, "release")));
      }
      try {
        await Promise.all(workers.map((worker, i) => {
          const partition = partitions[i + 1];
          return request(worker, "init", { partition }, [partition.rowIndex.buffer, partition.rowEntries.buffer,
            partition.columnIndex.buffer, partition.columnEntries.buffer]);
        }));
      } catch (error) { await close(); throw error; }
      async function product(operation, vector, out) {
        const pending = workers.map(worker => request(worker, operation, { vector }));
        const local = partitionProduct(partitions[0], operation, vector);
        const results = [local, ...await Promise.all(pending)];
        out.fill(0);
        for (let i = 0; i < results.length; i++) {
          const expected = operation === "multiply" ? ranges[i].rows : ranges[i].columns;
          if (!(results[i] instanceof Uint32Array) || results[i].length !== expected) {
            throw new Error("Invalid matrix worker product.");
          }
          out.set(results[i], operation === "multiply" ? ranges[i].rowStart : ranges[i].columnStart);
          workerProducts[i]++;
        }
      }
      const scratch = new Uint32Array(matrix.rowCount);
      return {
        workerProducts, close,
        async apply(input, output) {
          await product("multiply", input, scratch);
          await product("transpose", scratch, output);
        }
      };
    }
    function multiplySmall(a, b) {
      const out = new Uint32Array(32);
      for (let r = 0; r < 32; r++) {
        let bits = a[r], value = 0;
        while (bits) {
          const bit = 31 - Math.clz32(bits);
          value ^= b[bit]; bits ^= 1 << bit;
        }
        out[r] = value;
      }
      return out;
    }
    function multiplyAccumulate(v, matrix, out) {
      const table = new Uint32Array(4 * 256);
      for (let byte = 0; byte < 4; byte++) {
        const start = byte * 256;
        for (let j = 1; j < 256; j++) {
          const bit = 31 - Math.clz32(j & -j);
          table[start + j] = table[start + (j & (j - 1))] ^ matrix[8 * byte + bit];
        }
      }
      for (let i = 0; i < v.length; i++) {
        const word = v[i];
        out[i] ^= table[word & 255] ^ table[256 + ((word >>> 8) & 255)] ^
          table[512 + ((word >>> 16) & 255)] ^ table[768 + (word >>> 24)];
      }
    }
    function innerProduct(x, y) {
      const table = new Uint32Array(4 * 256), out = new Uint32Array(32);
      for (let i = 0; i < x.length; i++) {
        const xi = x[i], yi = y[i];
        table[xi & 255] ^= yi; table[256 + ((xi >>> 8) & 255)] ^= yi;
        table[512 + ((xi >>> 16) & 255)] ^= yi; table[768 + (xi >>> 24)] ^= yi;
      }
      for (let byte = 0; byte < 4; byte++) {
        for (let bit = 0; bit < 8; bit++) {
          let value = 0;
          for (let j = 0; j < 256; j++) if (j & (1 << bit)) value ^= table[byte * 256 + j];
          out[byte * 8 + bit] = value;
        }
      }
      return out;
    }
    function nonsingularSubmatrix(t, previous) {
      const left = t.slice(), right = Uint32Array.from({ length: 32 }, (_, i) => 1 << i);
      const order = new Uint8Array(32), selected = [];
      let mask = 0, at = 0;
      for (let i = 0; i < previous.length; i++) {
        order[31 - i] = previous[i]; mask |= 1 << previous[i];
      }
      for (let i = 0; i < 32; i++) if (!(mask & (1 << i))) order[at++] = i;
      function swap(a, b) {
        let temp = left[a]; left[a] = left[b]; left[b] = temp;
        temp = right[a]; right[a] = right[b]; right[b] = temp;
      }
      for (let i = 0; i < 32; i++) {
        const row = order[i], bit = 1 << row;
        let j = i;
        while (j < 32 && !(left[order[j]] & bit)) j++;
        if (j < 32) {
          swap(row, order[j]);
          for (let r = 0; r < 32; r++) {
            if (r !== row && (left[r] & bit)) { left[r] ^= left[row]; right[r] ^= right[row]; }
          }
          selected.push(row);
        } else {
          j = i;
          while (j < 32 && !(right[order[j]] & bit)) j++;
          if (j === 32) return null;
          swap(row, order[j]);
          for (let r = 0; r < 32; r++) {
            if (r !== row && (right[r] & bit)) { left[r] ^= left[row]; right[r] ^= right[row]; }
          }
          left[row] = right[row] = 0;
        }
      }
      return { selected, inverse: right };
    }
    function unpackBlock(v) {
      const vectors = Array.from({ length: 32 }, () => new Uint32Array(Math.ceil(v.length / 32)));
      for (let i = 0; i < v.length; i++) {
        let bits = v[i];
        while (bits) {
          const bit = 31 - Math.clz32(bits);
          vectors[bit][i >>> 5] |= 1 << (i & 31);
          bits ^= 1 << bit;
        }
      }
      return vectors;
    }
    function highest(vector) {
      for (let w = vector.length - 1; w >= 0; w--) {
        if (vector[w]) return w * 32 + 31 - Math.clz32(vector[w]);
      }
      return -1;
    }
    function xorInto(a, b) { for (let i = 0; i < a.length; i++) a[i] ^= b[i]; }
    function combineNullspace(matrix, x, v) {
      const bx = new Uint32Array(matrix.rowCount), bv = bx.slice();
      matrix.multiply(x, bx, true); matrix.multiply(v, bv, true);
      const candidates = unpackBlock(x).concat(unpackBlock(v));
      const images = unpackBlock(bx).concat(unpackBlock(bv));
      const imagePivots = new Map(), dependencyPivots = new Map(), dependencies = [];
      for (let i = 0; i < candidates.length; i++) {
        const image = images[i], vector = candidates[i];
        let bit;
        while ((bit = highest(image)) >= 0) {
          const previous = imagePivots.get(bit);
          if (previous === undefined) { imagePivots.set(bit, i); break; }
          xorInto(image, images[previous]); xorInto(vector, candidates[previous]);
        }
        if (bit >= 0) continue;
        // The 64 candidates can include zero or dependent vectors. Return
        // only independent, nonzero dependencies of B, not just of B'B.
        const reduced = vector.slice();
        while ((bit = highest(reduced)) >= 0) {
          const previous = dependencyPivots.get(bit);
          if (!previous) {
            dependencyPivots.set(bit, reduced); dependencies.push(vector); break;
          }
          xorInto(reduced, previous);
        }
      }
      return dependencies;
    }
    function* lanczosCore(matrix, random, report) {
      const n = matrix.n, x = Uint32Array.from({ length: n }, random);
      let v = [new Uint32Array(n), new Uint32Array(n), new Uint32Array(n)];
      let next = new Uint32Array(n), inverses = Array.from({ length: 3 }, () => new Uint32Array(32));
      let previousT = new Uint32Array(32), previousU = new Uint32Array(32);
      let previousSelected = Array.from({ length: 32 }, (_, i) => i), previousMask = -1;
      yield { matrix, input: x, output: v[0] };
      const initial = v[0].slice();
      let solved = 0, iterations = 0, lastReport = -Infinity;
      for (; iterations <= matrix.activeRows + 32;) {
        iterations++;
        yield { matrix, input: v[0], output: next };
        const t = innerProduct(v[0], next), u = innerProduct(next, next);
        if (!t.some(value => value)) {
          return { dependencies: combineNullspace(matrix, x, v[0]), iterations, solved };
        }
        const sub = nonsingularSubmatrix(t, previousSelected);
        if (!sub || !sub.selected.length) break;
        let mask = 0;
        for (const bit of sub.selected) mask |= 1 << bit;
        if (solved < matrix.activeRows - 32 && (mask | previousMask) !== -1) break;
        inverses[0] = sub.inverse;
        for (let i = 0; i < n; i++) next[i] &= mask;
        const term = new Uint32Array(32);
        for (let i = 0; i < 32; i++) term[i] = (u[i] & mask) ^ t[i];
        const d = multiplySmall(inverses[0], term);
        for (let i = 0; i < 32; i++) d[i] ^= 1 << i;
        multiplyAccumulate(v[0], d, next);
        const e = multiplySmall(inverses[1], t);
        for (let i = 0; i < 32; i++) e[i] &= mask;
        multiplyAccumulate(v[1], e, next);
        if (previousMask !== -1) {
          let f = multiplySmall(previousT, inverses[1]);
          for (let i = 0; i < 32; i++) f[i] ^= 1 << i;
          f = multiplySmall(inverses[2], f);
          for (let i = 0; i < 32; i++) term[i] = ((previousU[i] & previousMask) ^ previousT[i]) & mask;
          multiplyAccumulate(v[2], multiplySmall(f, term), next);
        }
        // Direct inner products retain the initial block. This replaces
        // msieve's additional small-matrix recurrence with a simpler exact
        // calculation, while retaining Montgomery's three-block recurrence.
        const correction = multiplySmall(inverses[0], innerProduct(v[0], initial));
        multiplyAccumulate(v[0], correction, x);
        const spare = v[2]; v = [next, v[0], v[1]]; next = spare;
        inverses = [inverses[2], inverses[0], inverses[1]];
        previousT = t; previousU = u; previousSelected = sub.selected; previousMask = mask;
        solved += sub.selected.length;
        if (performance.now() - lastReport >= 200) {
          report({ iterations, dimensionsSolved: solved, dimensions: matrix.activeRows });
          lastReport = performance.now();
        }
        if (solved > matrix.activeRows + 32) break;
      }
      return { dependencies: [], iterations, solved, breakdown: true };
    }
    function verifyDependency(rowCount, columns, vector) {
      if (vector.length !== Math.ceil(columns.length / 32)) return false;
      const tail = columns.length & 31;
      if (tail && (vector[vector.length - 1] >>> tail)) return false;
      const parity = new Uint8Array(rowCount);
      let used = false;
      for (let c = 0; c < columns.length; c++) {
        if (!(vector[c >>> 5] & (1 << (c & 31)))) continue;
        used = true;
        for (const r of columns[c]) parity[r] ^= 1;
      }
      return used && parity.every(bit => bit === 0);
    }
    function* lanczosSteps(rowCount, columns, options = {}) {
      checkColumns(rowCount, columns);
      const postRows = Math.min(rowCount, options.postRows ?? (rowCount > 64 ? 16 : 0));
      const matrix = packMatrix(rowCount, columns, postRows);
      const stats = { method: "block-lanczos", blockSize: 32, postRows, attempts: 0,
        iterations: 0, dependencies: 0, breakdowns: 0, threads: (options.workers?.length || 0) + 1 };
      let seed = options.seed ?? globalThis.crypto.getRandomValues(new Uint32Array(1))[0];
      // Randomness only chooses a Krylov starting block. It is unrelated to
      // primality testing, which continues to use cryptographic uniform bases.
      function random() {
        seed = (seed + 0x6d2b79f5) | 0;
        let value = Math.imul(seed ^ (seed >>> 15), 1 | seed);
        value ^= value + Math.imul(value ^ (value >>> 7), 61 | value);
        return (value ^ (value >>> 14)) >>> 0;
      }
      for (let attempt = 0; attempt < (options.maxAttempts ?? 4); attempt++) {
        stats.attempts++;
        const result = yield* lanczosCore(matrix, random, state => options.onProgress?.({
          ...state, attempt: attempt + 1, threads: stats.threads }));
        stats.iterations += result.iterations;
        if (result.breakdown) stats.breakdowns++;
        if (!result.dependencies.length) continue;
        for (const vector of result.dependencies) {
          if (!verifyDependency(rowCount, columns, vector)) throw new Error("Invalid block Lanczos dependency.");
        }
        stats.dependencies = result.dependencies.length;
        return { dependencies: result.dependencies, stats };
      }
      return { dependencies: [], stats };
    }
    function blockLanczos(rowCount, columns, options = {}) {
      const steps = lanczosSteps(rowCount, columns, { ...options, workers: [] });
      const scratch = new Uint32Array(rowCount);
      let step = steps.next(), products = 0;
      while (!step.done) {
        const { matrix, input, output } = step.value;
        matrix.multiply(input, scratch); matrix.transpose(scratch, output);
        products += 2;
        step = steps.next();
      }
      step.value.stats.workerProducts = [products];
      return step.value;
    }
    async function blockLanczosParallel(rowCount, columns, options = {}) {
      if (!options.workers?.length) return blockLanczos(rowCount, columns, options);
      const steps = lanczosSteps(rowCount, columns, options);
      let operator;
      try {
        let step = steps.next();
        while (!step.done) {
          const { matrix, input, output } = step.value;
          if (!operator) operator = await parallelOperator(matrix, options.workers);
          await operator.apply(input, output);
          step = steps.next();
        }
        step.value.stats.workerProducts = operator?.workerProducts || Array(options.workers.length + 1).fill(0);
        return step.value;
      } finally {
        steps.return();
        if (operator) await operator.close();
      }
    }
    function expandDependency(vector, groups, originalColumns) {
      const result = new Uint32Array(Math.ceil(originalColumns / 32));
      for (let c = 0; c < groups.length; c++) {
        if (vector[c >>> 5] & (1 << (c & 31))) {
          for (const original of groups[c]) result[original >>> 5] ^= 1 << (original & 31);
        }
      }
      return result;
    }
    return { filterMatrix, blockLanczos, blockLanczosParallel, createMatrixWorker, expandDependency, verifyDependency,
      multiplySmall, multiplyAccumulate, innerProduct, nonsingularSubmatrix };
  }
  if (typeof module !== "undefined" && module.exports) module.exports = { createMatrixAlgebra };
  else root.MPQSMatrix = { createMatrixAlgebra };
})(typeof globalThis !== "undefined" ? globalThis : this);
