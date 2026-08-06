(function (root) {
  'use strict';

  const ZERO = 0n;
  const ONE = 1n;
  const TWO = 2n;
  const THREE = 3n;
  const FOUR = 4n;
  const SIX = 6n;
  const EIGHT = 8n;

  function absBigInt(value) {
    return value < ZERO ? -value : value;
  }

  function gcd(a, b) {
    let x = absBigInt(a);
    let y = absBigInt(b);
    while (y !== ZERO) {
      const remainder = x % y;
      x = y;
      y = remainder;
    }
    return x;
  }

  function pow(base, exponent) {
    if (exponent < ZERO) throw new RangeError('Negative BigInt exponent.');
    let result = ONE;
    let factor = base;
    let power = exponent;
    while (power > ZERO) {
      if (power % TWO === ONE) result *= factor;
      power /= TWO;
      if (power > ZERO) factor *= factor;
    }
    return result;
  }

  function factor(value) {
    if (value < ONE) throw new RangeError('Factorization requires a positive integer.');
    let remaining = value;
    const factors = [];

    function divideOut(prime) {
      let count = ZERO;
      while (remaining % prime === ZERO) {
        remaining /= prime;
        count += ONE;
      }
      if (count > ZERO) factors.push({ prime, count });
    }

    divideOut(TWO);
    divideOut(THREE);
    let candidate = 5n;
    let step = TWO;
    while (candidate * candidate <= remaining) {
      divideOut(candidate);
      candidate += step;
      step = SIX - step;
    }
    if (remaining > ONE) factors.push({ prime: remaining, count: ONE });
    return factors;
  }

  function mergeFactorizations(left, right) {
    const merged = new Map();
    for (const item of [...left, ...right]) {
      const key = item.prime.toString();
      const current = merged.get(key);
      if (current) current.count += item.count;
      else merged.set(key, { prime: item.prime, count: item.count });
    }
    return [...merged.values()].sort((a, b) => a.prime < b.prime ? -1 : a.prime > b.prime ? 1 : 0);
  }

  function squareFreeFromFactors(factors) {
    let squareRootPart = ONE;
    let squareFreePart = ONE;
    for (const item of factors) {
      squareRootPart *= pow(item.prime, item.count / TWO);
      if (item.count % TWO === ONE) squareFreePart *= item.prime;
    }
    return { squareRootPart, squareFreePart };
  }

  function mobius(value) {
    if (value < ONE) throw new RangeError('Mobius requires a positive integer.');
    const factors = factor(value);
    let parity = ZERO;
    for (const item of factors) {
      if (item.count > ONE) return ZERO;
      parity = ONE - parity;
    }
    return parity === ZERO ? ONE : -ONE;
  }

  function totient(value) {
    if (value < ONE) throw new RangeError('Totient requires a positive integer.');
    let result = value;
    for (const item of factor(value)) result = (result / item.prime) * (item.prime - ONE);
    return result;
  }

  function jacobi(aValue, nValue) {
    if (nValue <= ZERO || nValue % TWO === ZERO) throw new RangeError('Jacobi denominator must be positive and odd.');
    let a = aValue % nValue;
    if (a < ZERO) a += nValue;
    let n = nValue;
    let result = ONE;

    while (a !== ZERO) {
      while (a % TWO === ZERO) {
        a /= TWO;
        const residue = n % EIGHT;
        if (residue === THREE || residue === 5n) result = -result;
      }
      const swap = a;
      a = n;
      n = swap;
      if (a % FOUR === THREE && n % FOUR === THREE) result = -result;
      a %= n;
    }
    return n === ONE ? result : ZERO;
  }

  function divisors(value) {
    let values = [ONE];
    for (const item of factor(value)) {
      const expanded = [];
      let primePower = ONE;
      let count = ZERO;
      while (count <= item.count) {
        for (const existing of values) expanded.push(existing * primePower);
        primePower *= item.prime;
        count += ONE;
      }
      values = expanded;
    }
    values.sort((a, b) => a < b ? -1 : a > b ? 1 : 0);
    return values;
  }

  function mapGet(map, key) {
    return map.get(key.toString()) ?? ZERO;
  }

  const lucasCache = new Map();

  function computeLucasCD(squareFreePart) {
    const cacheKey = squareFreePart.toString();
    if (lucasCache.has(cacheKey)) return lucasCache.get(cacheKey);
    if (squareFreePart <= ONE) throw new RangeError('Lucas Algorithm L requires a square-free index greater than one.');

    const modFour = squareFreePart % FOUR;
    const nPrime = modFour === ONE ? squareFreePart : TWO * squareFreePart;
    const degree = totient(nPrime) / TWO;
    const q = new Map();

    let k = ONE;
    while (k <= degree) {
      let value;
      if (k % TWO === ONE) {
        value = jacobi(squareFreePart, k);
      } else {
        const common = gcd(k, nPrime);
        const cosineResidue = ((squareFreePart - ONE) * k) % EIGHT;
        let cosine = ZERO;
        if (cosineResidue === ZERO) cosine = ONE;
        else if (cosineResidue === FOUR) cosine = -ONE;
        value = cosine * totient(common) * mobius(nPrime / common);
      }
      q.set(k.toString(), value);
      k += ONE;
    }

    const gamma = new Map([['0', ONE]]);
    const delta = new Map([['0', ONE]]);
    const maxGamma = degree / TWO;
    const maxDelta = (degree - ONE) / TWO;
    k = ONE;

    while (k <= maxGamma) {
      let gammaValue = ZERO;
      let j = ZERO;
      while (j < k) {
        const oddIndex = TWO * k - TWO * j - ONE;
        const evenIndex = oddIndex + ONE;
        gammaValue += squareFreePart * mapGet(q, oddIndex) * mapGet(delta, j);
        gammaValue -= mapGet(q, evenIndex) * mapGet(gamma, j);
        j += ONE;
      }
      const gammaDivisor = TWO * k;
      if (gammaValue % gammaDivisor !== ZERO) throw new Error('Algorithm L produced a non-integral C coefficient.');
      gammaValue /= gammaDivisor;
      gamma.set(k.toString(), gammaValue);

      if (k <= maxDelta) {
        let deltaValue = ZERO;
        j = ZERO;
        while (j < k) {
          const oddIndex = TWO * k + ONE - TWO * j;
          const evenIndex = oddIndex - ONE;
          deltaValue += mapGet(q, oddIndex) * mapGet(gamma, j);
          deltaValue -= mapGet(q, evenIndex) * mapGet(delta, j);
          j += ONE;
        }
        deltaValue += gammaValue;
        const deltaDivisor = TWO * k + ONE;
        if (deltaValue % deltaDivisor !== ZERO) throw new Error('Algorithm L produced a non-integral D coefficient.');
        delta.set(k.toString(), deltaValue / deltaDivisor);
      }
      k += ONE;
    }

    k = maxGamma + ONE;
    while (k <= degree) {
      gamma.set(k.toString(), mapGet(gamma, degree - k));
      k += ONE;
    }
    k = maxDelta + ONE;
    while (k < degree) {
      delta.set(k.toString(), mapGet(delta, degree - ONE - k));
      k += ONE;
    }

    const result = { squareFreePart, nPrime, degree, gamma, delta };
    lucasCache.set(cacheKey, result);
    return result;
  }

  function evaluateHomogeneous(coefficients, degree, x, y) {
    let result = mapGet(coefficients, ZERO);
    let yPower = ONE;
    let index = ONE;
    while (index <= degree) {
      yPower *= y;
      result = result * x + mapGet(coefficients, index) * yPower;
      index += ONE;
    }
    return result;
  }

  function rawFactors(a, b, exponent, decomposition, lucas) {
    const r = decomposition.squareFreePart;
    if (exponent % r !== ZERO) throw new Error('The exponent is not divisible by the square-free index.');
    const h = exponent / r;
    if (h % TWO !== ONE) throw new Error('The exponent quotient must be odd.');

    const x = pow(a, h);
    const y = pow(b, h);
    const cValue = evaluateHomogeneous(lucas.gamma, lucas.degree, x, y);
    const dDegree = lucas.degree - ONE;
    const dValue = evaluateHomogeneous(lucas.delta, dDegree, x, y);
    const radical = pow(decomposition.squareRootPart, h) * pow(r, (h + ONE) / TWO);
    const radicalTerm = radical * dValue;
    return {
      exponent,
      h,
      x,
      y,
      cValue,
      dValue,
      radical,
      radicalTerm,
      L: cValue - radicalTerm,
      M: cValue + radicalTerm
    };
  }

  function oddPart(value) {
    let result = value;
    while (result % TWO === ZERO) result /= TWO;
    return result;
  }

  function primitiveFactors(a, b, exponent, decomposition, lucas) {
    const r = decomposition.squareFreePart;
    const relevantDivisors = divisors(oddPart(exponent));
    let lNumerator = ONE;
    let lDenominator = ONE;
    let mNumerator = ONE;
    let mDenominator = ONE;
    const rows = [];

    for (const divisor of relevantDivisors) {
      if (gcd(r, divisor) !== ONE) continue;
      const mu = mobius(divisor);
      if (mu === ZERO) continue;
      const chi = jacobi(r, divisor);
      const raw = rawFactors(a, b, exponent / divisor, decomposition, lucas);
      const selectedL = chi === ONE ? raw.L : raw.M;
      const selectedM = chi === ONE ? raw.M : raw.L;

      if (mu === ONE) {
        lNumerator *= selectedL;
        mNumerator *= selectedM;
      } else {
        lDenominator *= selectedL;
        mDenominator *= selectedM;
      }
      rows.push({ divisor, mu, chi, raw, selectedL, selectedM });
    }

    if (lNumerator % lDenominator !== ZERO || mNumerator % mDenominator !== ZERO) {
      throw new Error('Primitive extraction did not divide exactly.');
    }
    return {
      L: lNumerator / lDenominator,
      M: mNumerator / mDenominator,
      lNumerator,
      lDenominator,
      mNumerator,
      mDenominator,
      rows
    };
  }

  function cyclotomicHomogeneous(a, b, index) {
    let numerator = ONE;
    let denominator = ONE;
    for (const divisor of divisors(index)) {
      const mu = mobius(index / divisor);
      if (mu === ZERO) continue;
      const term = pow(a, divisor) - pow(b, divisor);
      if (mu === ONE) numerator *= term;
      else denominator *= term;
    }
    if (denominator === ZERO || numerator % denominator !== ZERO) throw new Error('Cyclotomic evaluation did not divide exactly.');
    return numerator / denominator;
  }

  function analyze(a, b, exponent, sign) {
    const factorsA = factor(a);
    const factorsB = factor(b);
    const productFactors = mergeFactorizations(factorsA, factorsB);
    const decomposition = squareFreeFromFactors(productFactors);
    const r = decomposition.squareFreePart;
    const target = pow(a, exponent) + (sign === '+' ? pow(b, exponent) : -pow(b, exponent));
    const expectedSign = r % FOUR === ONE ? '-' : '+';
    const divisible = r > ONE && exponent % r === ZERO;
    const h = divisible ? exponent / r : null;
    const eligibleExponent = divisible && h % TWO === ONE;
    const eligible = eligibleExponent && sign === expectedSign;

    const result = {
      a,
      b,
      exponent,
      sign,
      target,
      factorsA,
      factorsB,
      productFactors,
      decomposition,
      expectedSign,
      h,
      eligibleExponent,
      eligible
    };

    if (!eligible) return result;

    const lucas = computeLucasCD(r);
    const raw = rawFactors(a, b, exponent, decomposition, lucas);
    const primitive = primitiveFactors(a, b, exponent, decomposition, lucas);
    const cyclotomicIndex = sign === '+' ? TWO * exponent : exponent;
    const cyclotomicValue = cyclotomicHomogeneous(a, b, cyclotomicIndex);
    const rawProduct = raw.L * raw.M;
    const primitiveProduct = primitive.L * primitive.M;

    if (primitiveProduct !== cyclotomicValue) throw new Error('Primitive L and M do not match the expected cyclotomic value.');
    if (target % rawProduct !== ZERO || target % primitiveProduct !== ZERO) throw new Error('Calculated factors do not divide the requested number exactly.');

    Object.assign(result, {
      lucas,
      raw,
      primitive,
      cyclotomicIndex,
      cyclotomicValue,
      rawProduct,
      rawCofactor: target / rawProduct,
      primitiveProduct,
      primitiveCofactor: target / primitiveProduct
    });
    return result;
  }

  root.AurifeuilleanMath = Object.freeze({
    absBigInt,
    gcd,
    pow,
    factor,
    mergeFactorizations,
    squareFreeFromFactors,
    mobius,
    totient,
    jacobi,
    divisors,
    computeLucasCD,
    evaluateHomogeneous,
    rawFactors,
    primitiveFactors,
    cyclotomicHomogeneous,
    analyze
  });
})(globalThis);
