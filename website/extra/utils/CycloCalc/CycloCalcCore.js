(function (root, factory) {
  const api = factory();
  if (typeof module === "object" && module.exports) module.exports = api;
  root.CyclotomicCore = api;
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
  "use strict";

  const ZERO = 0n;
  const ONE = 1n;
  const TWO = 2n;
  const THREE = 3n;

  function absBigInt(value) {
    return value < ZERO ? -value : value;
  }

  function gcd(a, b) {
    a = absBigInt(a);
    b = absBigInt(b);
    while (b !== ZERO) {
      const remainder = a % b;
      a = b;
      b = remainder;
    }
    return a;
  }

  function pow(base, exponent) {
    if (exponent < ZERO) throw new RangeError("The exponent must be nonnegative.");
    let result = ONE;
    let factor = base;
    let remaining = exponent;
    while (remaining > ZERO) {
      if (remaining % TWO === ONE) result *= factor;
      remaining /= TWO;
      if (remaining > ZERO) factor *= factor;
    }
    return result;
  }

  function bitLength(value) {
    let v = absBigInt(value);
    let length = ZERO;
    while (v > ZERO) {
      v >>= ONE;
      length += ONE;
    }
    return length;
  }

  function factorInteger(value) {
    let remaining = absBigInt(value);
    const factors = [];
    if (remaining < TWO) return factors;
    let exponent = ZERO;
    while (remaining % TWO === ZERO) {
      remaining /= TWO;
      exponent += ONE;
    }
    if (exponent > ZERO) factors.push([TWO, exponent]);
    let candidate = THREE;
    while (candidate * candidate <= remaining) {
      exponent = ZERO;
      while (remaining % candidate === ZERO) {
        remaining /= candidate;
        exponent += ONE;
      }
      if (exponent > ZERO) factors.push([candidate, exponent]);
      candidate += TWO;
    }
    if (remaining > ONE) factors.push([remaining, ONE]);
    return factors;
  }

  function isPrime(value) {
    if (value < TWO) return false;
    const factors = factorInteger(value);
    return factors.length === 1 && factors[0][0] === value && factors[0][1] === ONE;
  }

  function modularPow(base, exponent, modulus) {
    let result = ONE;
    let currentBase = base % modulus;
    let currentExponent = exponent;
    while (currentExponent > ZERO) {
      if (currentExponent % TWO === ONE) result = (result * currentBase) % modulus;
      currentExponent /= TWO;
      currentBase = (currentBase * currentBase) % modulus;
    }
    return result;
  }

  const primalityCache = new Map();

  function classifyPrimality(value, maxIterations = 10n) {
    const n = absBigInt(value);
    if (n === ONE) return { kind: "unit", label: "1" };
    if (n < TWO) return { kind: "composite", label: "Composite" };
    const cacheKey = `${n}|${maxIterations}`;
    if (primalityCache.has(cacheKey)) return primalityCache.get(cacheKey);

    // Match the companion calculator's trial division, extended through 100.
    let divisor = TWO;
    while (divisor <= 100n && divisor * divisor <= n) {
      if (n % divisor === ZERO) {
        const result = { kind: "composite", label: "Composite" };
        primalityCache.set(cacheKey, result);
        return result;
      }
      divisor += ONE;
    }
    if (divisor * divisor > n) {
      const result = { kind: "prime", label: "Prime" };
      primalityCache.set(cacheKey, result);
      return result;
    }

    let d = n - ONE;
    let r = ZERO;
    while (d % TWO === ZERO) {
      d /= TWO;
      r += ONE;
    }

    for (let iteration = ONE; iteration <= maxIterations; iteration += ONE) {
      const witness = TWO + ((iteration - ONE) % (n - THREE));
      let x = modularPow(witness, d, n);
      if (x === ONE || x === n - ONE) continue;

      let mayBePrime = false;
      for (let round = ONE; round < r; round += ONE) {
        x = (x * x) % n;
        if (x === n - ONE) {
          mayBePrime = true;
          break;
        }
      }
      if (!mayBePrime) {
        const result = { kind: "composite", label: "Composite" };
        primalityCache.set(cacheKey, result);
        return result;
      }
    }

    const result = { kind: "probable-prime", label: "Prime / probable prime" };
    primalityCache.set(cacheKey, result);
    return result;
  }

  function divisors(value) {
    if (value < ONE) throw new RangeError("Divisors require a positive integer.");
    let result = [ONE];
    for (const [prime, exponent] of factorInteger(value)) {
      const previous = result.slice();
      let power = ONE;
      for (let e = ONE; e <= exponent; e += ONE) {
        power *= prime;
        for (const divisor of previous) result.push(divisor * power);
      }
    }
    result.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
    return result;
  }

  function totient(value) {
    let result = value;
    for (const [prime] of factorInteger(value)) result = (result / prime) * (prime - ONE);
    return result;
  }

  function primePowerBase(value) {
    const factors = factorInteger(value);
    return factors.length === 1 ? factors[0][0] : null;
  }

  const cyclotomicCache = new Map();

  function cyclotomic(index, a, b) {
    if (index < ONE) throw new RangeError("Cyclotomic indices start at 1.");
    const key = `${index}|${a}|${b}`;
    if (cyclotomicCache.has(key)) return cyclotomicCache.get(key);
    let value;
    if (index === ONE) {
      value = a - b;
    } else if (a === b) {
      const base = primePowerBase(index);
      const atOne = base === null ? ONE : base;
      value = pow(a, totient(index)) * atOne;
    } else {
      value = pow(a, index) - pow(b, index);
      for (const divisor of divisors(index)) {
        if (divisor !== index) value /= cyclotomic(divisor, a, b);
      }
    }
    cyclotomicCache.set(key, value);
    return value;
  }

  function standardIndices(exponent, sign) {
    if (sign === "-") return divisors(exponent);
    if (sign !== "+") throw new Error("The sign must be + or -.");
    return divisors(TWO * exponent).filter((index) => exponent % index !== ZERO);
  }

  function integerSquareRoot(value) {
    if (value < ZERO) return null;
    if (value < TWO) return value;
    let x = ONE << ((bitLength(value) + ONE) / TWO);
    let y = (x + value / x) / TWO;
    while (y < x) {
      x = y;
      y = (x + value / x) / TWO;
    }
    return x;
  }

  function exactSquareRoot(value) {
    const root = integerSquareRoot(value);
    return root !== null && root * root === value ? root : null;
  }

  function comparePower(base, exponent, target) {
    let result = ONE;
    let factor = base;
    let remaining = exponent;
    while (remaining > ZERO) {
      if (remaining % TWO === ONE) {
        result *= factor;
        if (result > target) return ONE;
      }
      remaining /= TWO;
      if (remaining > ZERO) {
        factor *= factor;
        if (factor > target) factor = target + ONE;
      }
    }
    return result < target ? -ONE : result > target ? ONE : ZERO;
  }

  function exactRoot(value, exponent) {
    if (value < ZERO || exponent < ONE) return null;
    if (value < TWO || exponent === ONE) return value;
    let low = ONE;
    let high = TWO;
    while (comparePower(high, exponent, value) < ZERO) high *= TWO;
    while (low <= high) {
      const middle = (low + high) / TWO;
      const comparison = comparePower(middle, exponent, value);
      if (comparison === ZERO) return middle;
      if (comparison < ZERO) low = middle + ONE;
      else high = middle - ONE;
    }
    return null;
  }

  function commonPerfectPower(a, b) {
    if (a < ONE || b < ONE || (a === ONE && b === ONE)) return null;
    const limit = a === ONE ? bitLength(b) : b === ONE ? bitLength(a) :
      (bitLength(a) < bitLength(b) ? bitLength(a) : bitLength(b));
    for (let exponent = limit; exponent >= TWO; exponent -= ONE) {
      const aRoot = a === ONE ? ONE : exactRoot(a, exponent);
      if (aRoot === null) continue;
      const bRoot = b === ONE ? ONE : exactRoot(b, exponent);
      if (bRoot !== null) return { exponent, aRoot, bRoot };
    }
    return null;
  }

  function compositionIndices(index, powerExponent) {
    return divisors(powerExponent).filter((e) => gcd(index, powerExponent / e) === ONE)
      .map((e) => index * e);
  }

  function phiThreeSpecials(a, b) {
    const results = [];
    const s = exactSquareRoot(b);
    if (s !== null && (a + b) % THREE === ZERO) {
      const r = exactSquareRoot((a + b) / THREE);
      if (r !== null) {
        const first = THREE * r * r - THREE * r * s + s * s;
        const second = THREE * r * r + THREE * r * s + s * s;
        results.push({ type: "balanced", r, s, factors: [first, second] });
      }
    }
    const k = exactSquareRoot(THREE * a * (a + b));
    if (k !== null) {
      const center = TWO * a + b;
      results.push({ type: "pell", k, factors: [center - k, center + k] });
    }
    return results;
  }

  function phiPrime(a, b, prime) {
    let total = ZERO;
    for (let i = ZERO; i < prime; i += ONE) total += pow(a, prime - ONE - i) * pow(b, i);
    return total;
  }

  function solvePhiPrime(target, b, prime) {
    let low = ONE;
    let high = TWO;
    while (phiPrime(high, b, prime) < target) high *= TWO;
    while (low <= high) {
      const middle = (low + high) / TWO;
      const value = phiPrime(middle, b, prime);
      if (value === target) return middle;
      if (value < target) low = middle + ONE;
      else high = middle - ONE;
    }
    return null;
  }

  function nestedRepresentations(a, b) {
    const candidates = [];
    const upper = bitLength(a) + ONE;
    for (let prime = THREE; prime <= upper; prime += TWO) {
      if (!isPrime(prime)) continue;
      const d = b === ONE ? ONE : exactRoot(b, prime - ONE);
      if (d === null) continue;
      const c = solvePhiPrime(a, d, prime);
      if (c !== null) candidates.push({ prime, c, d });
    }
    // The C=D=1 representation is mathematically valid but rarely informative;
    // cap its primality probe so a huge input cannot trigger an enormous trial division.
    if (b === ONE && a <= 100000n && isPrime(a) && a % TWO === ONE && !candidates.some((item) => item.prime === a && item.c === ONE)) {
      candidates.push({ prime: a, c: ONE, d: ONE });
    }
    return candidates;
  }

  function nestedFactors(index, a, b) {
    const relationships = [];
    let t = null;
    let companionIndex = null;
    let form = null;
    if (index % TWO === ONE && isPrime(index)) {
      t = index;
      companionIndex = TWO * t;
      form = "odd";
    } else if (index % TWO === ZERO && isPrime(index / TWO) && index / TWO % TWO === ONE) {
      t = index / TWO;
      companionIndex = t;
      form = "even";
    }
    if (t === null) return relationships;
    for (const representation of nestedRepresentations(a, b)) {
      const modulus = TWO * t;
      if (representation.prime % modulus !== modulus - ONE) continue;
      const factor = cyclotomic(companionIndex, representation.c, representation.d);
      const whole = cyclotomic(index, a, b);
      if (factor > ONE && whole % factor === ZERO) {
        relationships.push({
          type: "nested",
          form,
          t,
          r: representation.prime,
          c: representation.c,
          d: representation.d,
          factorIndex: companionIndex,
          factor,
          cofactor: whole / factor
        });
      }
    }
    return relationships;
  }

  function analyzePhiThree(a, b) {
    const analyses = [];
    const perfectPower = commonPerfectPower(a, b);
    if (perfectPower !== null) {
      const indices = compositionIndices(THREE, perfectPower.exponent);
      analyses.push({
        type: "perfect-power",
        index: THREE,
        ...perfectPower,
        indices,
        factors: indices.map((index) => cyclotomic(index, perfectPower.aRoot, perfectPower.bRoot))
      });
    }
    analyses.push(...phiThreeSpecials(a, b));
    return analyses;
  }

  function factorRelationships(index, a, b) {
    const relationships = [];
    const directPower = commonPerfectPower(a, b);
    if (directPower !== null) {
      const indices = compositionIndices(index, directPower.exponent);
      relationships.push({
        type: "perfect-power",
        index,
        ...directPower,
        indices,
        factors: indices.map((childIndex) => cyclotomic(childIndex, directPower.aRoot, directPower.bRoot))
      });
    }

    if (index === THREE) {
      relationships.push(...phiThreeSpecials(a, b));
      relationships.push({ type: "shift-3-to-6", a, b, shiftedA: a + b, shiftedB: b });
    } else if (index === 6n && a !== b) {
      const high = a > b ? a : b;
      const low = a > b ? b : a;
      const shiftedA = high - low;
      const shift = {
        type: "shift-6-to-3",
        a,
        b,
        usedSymmetry: a < b,
        shiftedA,
        shiftedB: low,
        analyses: analyzePhiThree(shiftedA, low)
      };
      relationships.push(shift);
    }
    relationships.push(...nestedFactors(index, a, b));
    return relationships;
  }

  function calculate(a, b, exponent, sign) {
    if (a < TWO || b < ONE || exponent < ONE) throw new RangeError("A, B, and N are outside the permitted range.");
    const common = gcd(a, b);
    const reducedA = a / common;
    const reducedB = b / common;
    const indices = standardIndices(exponent, sign);
    const factors = indices.map((index) => ({
      index,
      value: cyclotomic(index, a, b),
      reducedValue: cyclotomic(index, reducedA, reducedB),
      scale: pow(common, totient(index)),
      relationships: factorRelationships(index, reducedA, reducedB)
    }));
    const leftPower = pow(a, exponent);
    const rightPower = pow(b, exponent);
    const value = sign === "+" ? leftPower + rightPower : leftPower - rightPower;
    return { a, b, exponent, sign, common, reducedA, reducedB, indices, factors, leftPower, rightPower, value };
  }

  function product(values) {
    return values.reduce((accumulator, value) => accumulator * value, ONE);
  }

  return {
    gcd,
    pow,
    factorInteger,
    isPrime,
    modularPow,
    classifyPrimality,
    divisors,
    totient,
    cyclotomic,
    standardIndices,
    exactRoot,
    exactSquareRoot,
    commonPerfectPower,
    compositionIndices,
    phiThreeSpecials,
    nestedRepresentations,
    nestedFactors,
    factorRelationships,
    calculate,
    product
  };
});
