// Cymatix physics: Kirchhoff thin-plate eigenmodes and their steady-state response to point forces.
// Circular plates use the exact Bessel-function solutions; rectangular plates use a Rayleigh–Ritz
// solution on orthonormalised polynomials. Units are SI. Λ is the dimensionless frequency ω·L²·√(ρh/D).

const TAU = 2 * Math.PI;
export const NR = 257; // radial samples per circular mode shape, ρ = r/a from 0 to 1
export const DRIVER_LIMIT = 0.97; // how far the driver can sit from the centre, as a share of the half-width

export const clamp = (x, lo, hi) => Math.max(lo, Math.min(hi, x));

const byFrequency = (p, q) => p.L - q.L;
const flexuralRigidity = ({ E, h, nu }) => E * h ** 3 / (12 * (1 - nu * nu));
const trapezoidWeight = (k, N) => (k === 0 || k === N ? 0.5 : 1);
const lerpAt = (values, i, t) => values[i] + t * (values[i + 1] - values[i]);
const toUnitInterval = (u, halfWidth) => (u / halfWidth + 1) / 2; // plate coordinate → ξ ∈ [0, 1]
const range = n => [...Array(n).keys()];

// ---------- circular plates ----------

// The integral for I_n sums terms about e^(n²/2x) times larger than the answer, so past this point it has too few digits
// left, and its rounding noise would look like extra modes. The power series has only positive terms, so it stays exact.
const integralLosesDigits = (n, x) => n * n > 30 * Math.abs(x);

// I_n(x) = (x/2)^n/n! · Σ_k (x²/4)^k / (k!·(n+1)_k)
function besselISeries(n, x) {
  const half = x / 2;
  let term = 1;
  for (let k = 1; k <= n; k++) term *= half / k;
  const step = half * half;
  let sum = term;
  for (let k = 1; term > 1e-17 * Math.abs(sum); k++) {
    term *= step / (k * (n + k));
    sum += term;
  }
  return sum;
}

/** [J_n(x), J_n'(x), I_n(x), I_n'(x)] for integer n ≥ 0, from Bessel's integrals. */
export function bessel(n, x) {
  // The trapezoid rule on these smooth periodic integrands converges exponentially (≈1e-15 here).
  const N = Math.ceil(n + Math.abs(x)) + 40;
  const h = Math.PI / N;
  let J = 0, Jp = 0, I = 0, Ip = 0;
  for (let k = 0; k <= N; k++) {
    const t = k * h;
    const s = Math.sin(t);
    const c = Math.cos(t);
    const w = trapezoidWeight(k, N);
    const phase = n * t - x * s;
    const e = w * Math.exp(x * c) * Math.cos(n * t);
    J += w * Math.cos(phase);
    Jp += w * s * Math.sin(phase);
    I += e;
    Ip += e * c;
  }
  if (x !== 0 && integralLosesDigits(n, x)) {
    const In = besselISeries(n, x);
    return [J / N, Jp / N, In, besselISeries(n + 1, x) + n / x * In];
  }
  return [J / N, Jp / N, I / N, Ip / N];
}

// Edge-condition rows at ρ = 1 for w = J_n(λρ) + C·I_n(λρ); a mode exists where their 2×2 determinant is zero.
function edgeRows(bc, n, lam, nu) {
  const [J, Jp, I, Ip] = bessel(n, lam);
  if (bc === 'clamped') return [J, I, Jp, Ip]; // w = 0, ∂w/∂r = 0

  const n2 = n * n;
  const l2 = lam * lam;
  // second derivatives, from Bessel's equations
  const Jpp = -Jp / lam - (1 - n2 / l2) * J;
  const Ipp = -Ip / lam + (1 + n2 / l2) * I;
  // bending moment M_r
  const mJ = l2 * Jpp + nu * (lam * Jp - n2 * J);
  const mI = l2 * Ipp + nu * (lam * Ip - n2 * I);
  if (bc === 'simply') return [J, I, mJ, mI]; // w = 0, M_r = 0

  // Kirchhoff shear V_r
  const vJ = -l2 * lam * Jp - (1 - nu) * n2 * (lam * Jp - J);
  const vI = l2 * lam * Ip - (1 - nu) * n2 * (lam * Ip - I);
  return [mJ, mI, vJ, vI]; // free: M_r = 0, V_r = 0
}

function radialShape(n, lam, C) {
  const R = new Float32Array(NR);
  if (!lam) { // rigid translation / rocking
    for (let j = 0; j < NR; j++) R[j] = (j / (NR - 1)) ** n;
    return R;
  }
  // Bessel's integrals as in bessel(), with the per-angle terms computed once for all NR radii
  const N = Math.ceil(n + lam) + 40;
  const terms = Array.from({ length: N + 1 }, (_, k) => {
    const t = k * Math.PI / N;
    return { sin: Math.sin(t), cos: Math.cos(t), nt: n * t, weight: trapezoidWeight(k, N) };
  });
  for (let j = 0; j < NR; j++) {
    const x = lam * j / (NR - 1);
    const useSeries = integralLosesDigits(n, x);
    let J = 0, I = 0;
    for (const { sin, cos, nt, weight } of terms) {
      J += weight * Math.cos(nt - x * sin);
      if (!useSeries) I += weight * Math.exp(x * cos) * Math.cos(nt);
    }
    R[j] = useSeries ? J / N + C * besselISeries(n, x) : (J + C * I) / N;
  }
  return R;
}

const isNegligible = (value, Rmax) => Math.abs(value) < 1e-6 * Rmax;

function circleMode(n, lam, C) {
  const R = radialShape(n, lam, C);
  let norm = 0, Rmax = 0;
  for (let j = 0; j < NR; j++) {
    const simpson = j === 0 || j === NR - 1 ? 1 : j % 2 ? 4 : 2;
    norm += simpson * R[j] * R[j] * j / (NR - 1); // ∫R²ρdρ
    Rmax = Math.max(Rmax, Math.abs(R[j]));
  }
  norm /= 3 * (NR - 1);

  let nodalCircles = 0, sign = 0;
  for (let j = 1; j < NR - 1; j++) {
    if (isNegligible(R[j], Rmax)) continue;
    if (sign && Math.sign(R[j]) !== sign) nodalCircles++;
    sign = Math.sign(R[j]);
  }
  return {
    n,
    m: nodalCircles,
    L: lam * lam,
    C,
    R,
    Rmax,
    // inv folds in ∫cos²(nθ)dθ so that a cos/sin degenerate pair sums to R(ρ)R(ρ₀)cos n(θ−θ₀)·inv
    inv: 1 / (Math.PI * (n ? 1 : 2) * norm),
    label: lam ? `(${n},${nodalCircles})` : 'rigid',
  };
}

// 50 bisection steps on [a, b], across which f changes sign; fa is f(a).
function bisect(f, a, b, fa = f(a)) {
  for (let i = 0; i < 50; i++) {
    const c = (a + b) / 2;
    const fc = f(c);
    if (fa * fc <= 0) {
      b = c;
    } else {
      a = c;
      fa = fc;
    }
  }
  return (a + b) / 2;
}

/** Radii r/a of a circular mode's nodal circles: sign changes of the sampled shape, refined on J_n(λρ) + C·I_n(λρ). */
export function nodalRadii(mode) {
  if (!mode.L) return [];
  const { n, C, R, Rmax } = mode;
  const lam = Math.sqrt(mode.L);
  const shape = r => {
    const [J, , I] = bessel(n, lam * r);
    return J + C * I;
  };
  const radii = [];
  for (let j = 1; j < NR - 2; j++) {
    if (isNegligible(R[j], Rmax) || Math.sign(R[j]) === Math.sign(R[j + 1])) continue;
    radii.push(bisect(shape, j / (NR - 1), (j + 1) / (NR - 1)));
  }
  return radii;
}

/** Exact modes of a circular plate (radius 1) up to λ = lamMax; edge bc is 'free' | 'clamped' | 'simply'. */
export function circleModes(bc, nu, lamMax) {
  const determinant = (n, lam) => {
    const [a, b, c, d] = edgeRows(bc, n, lam, nu);
    return a * d - b * c;
  };
  const modes = bc === 'free' ? [circleMode(0, 0, 0), circleMode(1, 0, 0)] : [];
  // the first root grows with n, so stop at the first n with none
  for (let n = 0, found = true; found && n < lamMax; n++) {
    found = false;
    // roots of order n lie above λ = n; starting at 0.8n also keeps the scan out of the range where J_n(λ) is rounding noise
    let lo = Math.max(0.5, 0.8 * n);
    let fLo = determinant(n, lo);
    while (lo < lamMax) {
      const hi = lo + 0.1; // roots of one n are ~π apart, so a 0.1 scan cannot skip one
      const fHi = determinant(n, hi);
      if (fLo * fHi < 0) {
        const lam = bisect(l => determinant(n, l), lo, hi, fLo);
        const [a, b] = edgeRows(bc, n, lam, nu);
        modes.push(circleMode(n, lam, -a / b));
        found = true;
      }
      lo = hi;
      fLo = fHi;
    }
  }
  return modes.sort(byFrequency);
}

// ---------- rectangular plates ----------

export function gaussLegendre(n) { // nodes and weights on [0, 1]
  const nodes = new Float64Array(n);
  const weights = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let z = Math.cos(Math.PI * (i + 0.75) / (n + 0.5));
    let slope = 1;
    for (let iter = 0; iter < 100; iter++) {
      let p1 = 1, p2 = 0;
      for (let j = 1; j <= n; j++) {
        const p3 = p2;
        p2 = p1;
        p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
      }
      slope = n * (z * p1 - p2) / (z * z - 1);
      const dz = p1 / slope;
      z -= dz;
      if (Math.abs(dz) < 1e-15) break;
    }
    nodes[i] = (1 - z) / 2;
    weights[i] = 1 / ((1 - z * z) * slope * slope);
  }
  return [nodes, weights];
}

export const ORDER = { free: 0, simply: 1, clamped: 2 }; // vanishing order of w at an edge: none, w = 0, w = w' = 0

// Raw 1D Ritz basis on ξ ∈ [0,1]: shifted Legendre polynomials times ξ^p·(1−ξ)^q, which enforces the geometric
// condition at each end (p, q from ORDER). Fills value, d/dξ, d²/dξ².
export function rawBasis(p, q, nb, xi, f, f1, f2) {
  const t = 2 * xi - 1;
  const a = xi;
  const b = 1 - xi;
  const pow = (x, k) => (k > 0 ? x ** k : 1);

  // the weight ξ^p·(1−ξ)^q and its derivatives
  const w = pow(a, p) * pow(b, q);
  const w1 = (p ? p * pow(a, p - 1) * pow(b, q) : 0) - (q ? q * pow(a, p) * pow(b, q - 1) : 0);
  const w2 = (p > 1 ? p * (p - 1) * pow(a, p - 2) * pow(b, q) : 0)
    - (p && q ? 2 * p * q * pow(a, p - 1) * pow(b, q - 1) : 0)
    + (q > 1 ? q * (q - 1) * pow(a, p) * pow(b, q - 2) : 0);

  // P_k and its first two t-derivatives, with the previous order's values for the recurrence
  let P = 1, d = 0, s = 0, prevP = 0, prevD = 0, prevS = 0;
  for (let k = 0; k < nb; k++) {
    f[k] = w * P;
    f1[k] = w1 * P + 2 * w * d;
    f2[k] = w2 * P + 4 * w1 * d + 4 * w * s;
    const nextP = ((2 * k + 1) * t * P - k * prevP) / (k + 1);
    const nextD = ((2 * k + 1) * (P + t * d) - k * prevD) / (k + 1);
    const nextS = ((2 * k + 1) * (2 * d + t * s) - k * prevS) / (k + 1);
    prevP = P;
    prevD = d;
    prevS = s;
    P = nextP;
    d = nextD;
    s = nextS;
  }
}

// 1D basis orthonormalised on [0,1] (Gram–Schmidt, twice for stability) with its stiffness integrals
// E11 = ∫U_i'U_k', E22 = ∫U_i''U_k'', E20 = ∫U_i''U_k, and an evaluator at any ξ.
export function basis1D(p, q, nb) {
  const nq = nb + 12;
  const [qx, qw] = gaussLegendre(nq);
  const integrate = (a, b) => {
    let total = 0;
    for (let j = 0; j < nq; j++) total += qw[j] * a[j] * b[j];
    return total;
  };

  // the raw functions and their derivatives at the quadrature points
  const f = new Float64Array(nb);
  const f1 = new Float64Array(nb);
  const f2 = new Float64Array(nb);
  const newTable = () => Array.from({ length: nb }, () => new Float64Array(nq));
  const V = newTable(), V1 = newTable(), V2 = newTable();
  for (let j = 0; j < nq; j++) {
    rawBasis(p, q, nb, qx[j], f, f1, f2);
    for (let k = 0; k < nb; k++) {
      V[k][j] = f[k];
      V1[k][j] = f1[k];
      V2[k][j] = f2[k];
    }
  }

  // Gram–Schmidt; row k of C maps raw → orthonormal
  const C = [], U = [];
  for (let k = 0; k < nb; k++) {
    const u = Float64Array.from(V[k]);
    const ck = new Float64Array(nb);
    ck[k] = 1;
    for (let pass = 0; pass < 2; pass++) {
      for (let j = 0; j < k; j++) {
        const projection = integrate(u, U[j]);
        for (let i = 0; i < nq; i++) u[i] -= projection * U[j][i];
        for (let i = 0; i < nb; i++) ck[i] -= projection * C[j][i];
      }
    }
    const scale = 1 / Math.sqrt(integrate(u, u));
    for (let i = 0; i < nq; i++) u[i] *= scale;
    for (let i = 0; i < nb; i++) ck[i] *= scale;
    U.push(u);
    C.push(ck);
  }

  const orthonormalised = W => C.map(ck => {
    const out = new Float64Array(nq);
    for (let i = 0; i < nb; i++) {
      for (let j = 0; j < nq; j++) out[j] += ck[i] * W[i][j];
    }
    return out;
  });
  const integrals = (A, B) => A.map(a => B.map(b => integrate(a, b)));
  const U1 = orthonormalised(V1);
  const U2 = orthonormalised(V2);
  const at = x => {
    rawBasis(p, q, nb, x, f, f1, f2);
    return C.map(ck => ck.reduce((total, c, i) => total + c * f[i], 0));
  };
  return {
    at,
    E11: integrals(U1, U1),
    E22: integrals(U2, U2),
    E20: integrals(U2, U),
    sym: p === q, // equal ends: Legendre parity survives the symmetric weight, so functions alternate even/odd about ξ = ½
  };
}

// Householder reduction of the symmetric matrix in V to tridiagonal form (diagonal d, off-diagonal e). From JAMA's tred2.
function tridiagonalise(V, d, e, n) {
  for (let j = 0; j < n; j++) d[j] = V[(n - 1) * n + j];
  for (let i = n - 1; i > 0; i--) {
    let scale = 0, h = 0;
    for (let k = 0; k < i; k++) scale += Math.abs(d[k]);
    if (scale === 0) {
      e[i] = d[i - 1];
      for (let j = 0; j < i; j++) {
        d[j] = V[(i - 1) * n + j];
        V[i * n + j] = 0;
        V[j * n + i] = 0;
      }
    } else {
      for (let k = 0; k < i; k++) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      let f = d[i - 1];
      let g = Math.sqrt(h);
      if (f > 0) g = -g;
      e[i] = scale * g;
      h -= f * g;
      d[i - 1] = f - g;
      for (let j = 0; j < i; j++) e[j] = 0;
      for (let j = 0; j < i; j++) {
        f = d[j];
        V[j * n + i] = f;
        g = e[j] + V[j * n + j] * f;
        for (let k = j + 1; k <= i - 1; k++) {
          g += V[k * n + j] * d[k];
          e[k] += V[k * n + j] * f;
        }
        e[j] = g;
      }
      f = 0;
      for (let j = 0; j < i; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      const hh = f / (h + h);
      for (let j = 0; j < i; j++) e[j] -= hh * d[j];
      for (let j = 0; j < i; j++) {
        f = d[j];
        g = e[j];
        for (let k = j; k <= i - 1; k++) V[k * n + j] -= f * e[k] + g * d[k];
        d[j] = V[(i - 1) * n + j];
        V[i * n + j] = 0;
      }
    }
    d[i] = h;
  }

  // accumulate the transformations into V
  for (let i = 0; i < n - 1; i++) {
    V[(n - 1) * n + i] = V[i * n + i];
    V[i * n + i] = 1;
    const h = d[i + 1];
    if (h !== 0) {
      for (let k = 0; k <= i; k++) d[k] = V[k * n + i + 1] / h;
      for (let j = 0; j <= i; j++) {
        let g = 0;
        for (let k = 0; k <= i; k++) g += V[k * n + i + 1] * V[k * n + j];
        for (let k = 0; k <= i; k++) V[k * n + j] -= g * d[k];
      }
    }
    for (let k = 0; k <= i; k++) V[k * n + i + 1] = 0;
  }
  for (let j = 0; j < n; j++) {
    d[j] = V[(n - 1) * n + j];
    V[(n - 1) * n + j] = 0;
  }
  V[(n - 1) * n + n - 1] = 1;
  for (let i = 1; i < n; i++) e[i - 1] = e[i];
  e[n - 1] = 0;
}

// Implicit QL iteration on the tridiagonal form: d becomes the eigenvalues and V's columns the eigenvectors. From JAMA's tql2.
function implicitQL(V, d, e, n) {
  const eps = 2 ** -52;
  let f = 0, tst1 = 0;
  for (let l = 0; l < n; l++) {
    tst1 = Math.max(tst1, Math.abs(d[l]) + Math.abs(e[l]));
    let m = l;
    while (m < n - 1 && Math.abs(e[m]) > eps * tst1) m++;
    if (m > l) {
      do {
        let g = d[l];
        let p = (d[l + 1] - g) / (2 * e[l]);
        let r = Math.hypot(p, 1);
        if (p < 0) r = -r;
        d[l] = e[l] / (p + r);
        d[l + 1] = e[l] * (p + r);
        const dl1 = d[l + 1];
        let h = g - d[l];
        for (let i = l + 2; i < n; i++) d[i] -= h;
        f += h;
        p = d[m];
        let c = 1, c2 = 1, c3 = 1, s = 0, s2 = 0;
        const el1 = e[l + 1];
        for (let i = m - 1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = Math.hypot(p, e[i]);
          e[i + 1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i + 1] = h + s * (c * g + s * d[i]);
          for (let k = 0; k < n; k++) {
            h = V[k * n + i + 1];
            V[k * n + i + 1] = s * V[k * n + i] + c * h;
            V[k * n + i] = c * V[k * n + i] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;
      } while (Math.abs(e[l]) > eps * tst1);
    }
    d[l] += f;
    e[l] = 0;
  }
}

// Eigen-decomposition of a symmetric n×n matrix (JAMA, public domain). Eigenvectors are the columns of `vectors`.
export function symEig(A, n) {
  const V = Float64Array.from(A);
  const d = new Float64Array(n);
  const e = new Float64Array(n);
  tridiagonalise(V, d, e, n);
  implicitQL(V, d, e, n);
  return { values: d, vectors: V };
}

// Each symmetric axis splits its functions into even and odd ones, which the stiffness never couples.
export const parityClasses = (basis, n) => (basis.sym
  ? [range(n).filter(i => i % 2 === 0), range(n).filter(i => i % 2 === 1)]
  : [range(n)]);

// Basis size by the number of symmetric axes: each one halves the block size, so it affords a richer basis.
const DEFAULT_BASIS = [22, 28, 36];
const MAX_BASIS = [48, 68, 96]; // blocks of about 2300 unknowns: tens of seconds in worker threads

/**
 * Rayleigh–Ritz modes of a rectangular plate with sides 1 (ξ) × beta (η). `edges` is one condition for all four
 * edges or [ξ = 0, ξ = 1, η = 0, η = 1], each 'free' | 'simply' | 'clamped'.
 */
export function rectModes(edges, nu, beta, nb) {
  const problem = rectProblem(edges, nu, beta, nb);
  return assembleRectModes(problem, problem.blocks.map((_, block) => solveRectBlock(problem, block)));
}

// The basis and the independent parity blocks of a rectangle. Each block can be solved separately, even in another thread.
function rectProblem(edges, nu, beta, nb) {
  const edgeList = typeof edges === 'string' ? [edges, edges, edges, edges] : edges;
  const [x0, x1, y0, y1] = edgeList.map(edge => ORDER[edge]);
  const symAxes = (x0 === x1) + (y0 === y1);
  nb ??= DEFAULT_BASIS[symAxes];

  // The same number of functions per unit length on both sides: an elongated plate gets more along its long side
  // (b/a = 8 needs ~100 there for its first 15 modes), keeping nbx·nby ≈ nb².
  const r = Math.sqrt(beta);
  const nbx = Math.max(8, Math.round(nb / r));
  const nby = Math.max(8, Math.round(nb * r));
  const X = basis1D(x0, x1, nbx);
  const Y = x0 === y0 && x1 === y1 && nbx === nby ? X : basis1D(y0, y1, nby);
  const blocks = parityClasses(X, nbx).flatMap(I => parityClasses(Y, nby).map(J => ({ I, J })));
  return { nu, beta, nbx, nby, X, Y, blocks };
}

function solveRectBlock({ nu, beta, X, Y, blocks }, block) {
  const { I, J } = blocks[block];
  const b2 = beta * beta;
  const b4 = b2 * b2;
  const stiffness = (i, j, k, l) => (j === l ? X.E22[i][k] : 0) + (i === k ? Y.E22[j][l] / b4 : 0)
    + (nu * (X.E20[i][k] * Y.E20[l][j] + X.E20[k][i] * Y.E20[j][l]) + 2 * (1 - nu) * X.E11[i][k] * Y.E11[j][l]) / b2;
  const pairs = I.flatMap(i => J.map(j => [i, j]));
  const m = pairs.length;
  const K = new Float64Array(m * m);
  pairs.forEach(([i, j], a) => {
    pairs.forEach(([k, l], b) => { K[a * m + b] = stiffness(i, j, k, l); });
  });
  return symEig(K, m);
}

function assembleRectModes({ nbx, nby, X, Y, blocks }, solutions) {
  const modes = [];
  solutions.forEach(({ values, vectors }, block) => {
    const { I, J } = blocks[block];
    const pairs = I.flatMap(i => J.map(j => [i, j]));
    const m = pairs.length;
    for (let e = 0; e < m; e++) {
      modes.push({
        L: Math.sqrt(Math.max(values[e], 0)),
        solution: { pairs, vectors, m, e },
        inv: 1,
        Rmax: 1,
        // parity about the plate's centre lines: 0 even, 1 odd, −1 none (that axis has unequal ends)
        px: X.sym ? I[0] : -1,
        py: Y.sym ? J[0] : -1,
      });
    }
  });
  modes.sort(byFrequency);

  // shortcut: only the lowest ~nbx·nby/4 Ritz modes are converged; raise nb for denser high-frequency patterns.
  const kept = modes.slice(0, Math.round(nbx * nby / 4));
  let elastic = 0;
  for (const mode of kept) {
    const { pairs, vectors, m, e } = mode.solution;
    const vec = new Float64Array(nbx * nby);
    pairs.forEach(([i, j], a) => { vec[i * nby + j] = vectors[a * m + e]; });
    delete mode.solution;
    mode.vec = vec;
    mode.label = mode.L < 1e-3 ? 'rigid' : `#${++elastic}`;
  }
  return { modes: kept, nbx, nby, atX: X.at, atY: Y.at };
}

/** How many parity blocks a rectangle splits into, for solving them in parallel with rectBlock. */
export const rectBlockCount = (edges, nu, beta, nb) => rectProblem(edges, nu, beta, nb).blocks.length;

/** Solves one parity block: the numbers a worker thread sends back to rectModesFromBlocks. */
export function rectBlock(edges, nu, beta, nb, block) {
  return solveRectBlock(rectProblem(edges, nu, beta, nb), block);
}

/** The same result as rectModes, from blocks solved elsewhere, in block order. */
export function rectModesFromBlocks(edges, nu, beta, nb, solutions) {
  return assembleRectModes(rectProblem(edges, nu, beta, nb), solutions);
}

/**
 * Mindlin correction (shear deformation + rotary inertia) of a Kirchhoff angular frequency wK. The mode's bending
 * wavenumber k (k⁴ = wK²ρh/D) is kept and Mindlin's dispersion relation
 *   ρh·I·ω⁴ − (S·I·k² + ρh·D·k² + ρh·S)·ω² + S·D·k⁴ = 0,  S = κ²Gh, I = ρh³/12,
 * is solved for the flexural root. Exact for simply supported rectangles; a local approximation for other edges.
 */
export function mindlinOmega(wK, { E, rho, nu, h }, kappa2 = Math.PI ** 2 / 12) {
  if (!wK) return 0;
  const D = flexuralRigidity({ E, h, nu });
  const rhoH = rho * h;
  const k2 = wK * Math.sqrt(rhoH / D);
  const S = kappa2 * E / (2 * (1 + nu)) * h;
  const I = rho * h ** 3 / 12;
  const a = rhoH * I;
  const b = S * I * k2 + rhoH * D * k2 + rhoH * S;
  const c = S * D * k2 * k2;
  return Math.sqrt(2 * c / (b + Math.sqrt(b * b - 4 * a * c))); // lower root, written to avoid cancellation
}

// ---------- how many modes a plate needs ----------

const LAM_MAX_INLINE = 50; // above this λ a circular plate is worth solving in a worker thread
const circleReach = (wScale, fmax) => Math.max(10, 5 * Math.ceil(Math.sqrt(1.3 * TAU * fmax / wScale) / 5));

/** The λ a circular plate needs to reach 30 % past fmax, up to lamLimit. */
export function circleLamFor(cfg, lamLimit = 120) {
  const D = flexuralRigidity(cfg);
  const wScale = Math.sqrt(D / (cfg.rho * cfg.h)) / (cfg.d / 2) ** 2;
  return Math.min(lamLimit, circleReach(wScale, cfg.fmax));
}

/** The basis size a rectangle needs to reach 20 % past fmax, or null if the default basis already does. */
export function rectBasisFor(cfg) {
  const edges = cfg.edges ?? [cfg.bc, cfg.bc, cfg.bc, cfg.bc];
  const [x0, x1, y0, y1] = edges.map(edge => ORDER[edge]);
  const symAxes = (x0 === x1) + (y0 === y1);
  const beta = cfg.Ly / cfg.Lx;
  const D = flexuralRigidity(cfg);
  const targetL = 1.2 * TAU * cfg.fmax * cfg.Lx ** 2 * Math.sqrt(cfg.rho * cfg.h / D);
  // about β·Λ/(4π) modes lie below Λ (Weyl's law), and the solver keeps nb²/4 of them
  const wanted = Math.ceil(2 * Math.sqrt(1.1 * beta * targetL / (4 * Math.PI)));
  if (wanted <= DEFAULT_BASIS[symAxes]) return null;
  return Math.min(wanted, MAX_BASIS[symAxes]);
}

// ---------- the driven plate ----------

/** A plate with its modes, a G×G field grid over [-1,1]², a point driver, and the plate-acceleration field. */
export class Plate {
  // cfg: { shape: 'circle'|'rect', bc, E, rho, nu, h, d, Lx, Ly (all SI), zeta, fmax (Hz), G,
  //        mount: 'force'|'shaker'|'rigid', ma (shaker moving mass, kg), mindlin (default true) }
  configure(cfg) {
    const circle = cfg.shape === 'circle';
    const D = flexuralRigidity(cfg);
    const rhoH = cfg.rho * cfg.h;
    const L = circle ? cfg.d / 2 : cfg.Lx;
    this.circle = circle;
    this.wScale = Math.sqrt(D / rhoH) / (L * L); // ω = Λ·wScale
    this.zInf = 8 * Math.sqrt(D * rhoH); // an infinite plate's point impedance F/v, used for the modes above the cut-off
    this.mass = circle ? rhoH * L * L : rhoH * cfg.Lx * cfg.Ly; // modal mass scale for unit-normalised shapes
    this.zeta = cfg.zeta;
    this.mount = cfg.mount ?? 'force';
    this.ma = cfg.ma ?? 0;

    const key = circle ? this.#loadCircleModes(cfg) : this.#loadRectModes(cfg);

    // solved.thick: those modes already include shear and rotary inertia
    const withThickness = w => (cfg.mindlin === false || cfg.solved?.thick ? w : mindlinOmega(w, cfg));
    for (const mode of this.solvedModes) {
      mode.w = withThickness(mode.L * this.wScale);
      mode.f = mode.w / TAU;
    }
    const topL = circle ? this.lamMax ** 2 : this.solvedModes.at(-1).L;
    this.fTop = withThickness(topL * this.wScale) / TAU;

    const longSide = Math.max(cfg.Lx, cfg.Ly);
    this.hx = circle ? 1 : cfg.Lx / longSide;
    this.hy = circle ? 1 : cfg.Ly / longSide;
    // solvedModes are reused between configurations; a fluid or similar returns adjusted copies through cfg.adjustModes
    this.modes = cfg.adjustModes?.(this.solvedModes) ?? this.solvedModes;

    if (key !== this.key || cfg.G !== this.G) {
      this.key = key;
      this.#buildGrid(cfg.G);
    }
    this.setDrive(this.du ?? 0, this.dv ?? 0);
  }

  // cfg.solved carries modes solved elsewhere (a worker thread, a larger basis); without it they are solved here.
  #loadCircleModes({ bc, nu, fmax, solved }) {
    // modes reach 30 % past the top analysed frequency; λ is rounded up so small edits reuse them
    this.lamMax = solved?.lamMax ?? Math.min(LAM_MAX_INLINE, circleReach(this.wScale, fmax));
    const key = `c|${bc}|${nu}|${this.lamMax}|${solved?.id ?? ''}`;
    if (key !== this.key) this.solvedModes = solved?.modes ?? circleModes(bc, nu, this.lamMax);
    return key;
  }

  #loadRectModes({ edges, bc, nu, Lx, Ly, solved }) {
    this.edges = edges ?? [bc, bc, bc, bc];
    const key = `r|${this.edges}|${nu}|${Ly / Lx}|${solved?.id ?? ''}`;
    if (key !== this.key) {
      this.ritz = solved?.ritz ?? rectModes(this.edges, nu, Ly / Lx);
      this.solvedModes = this.ritz.modes;
    }
    return key;
  }

  contains(u, v) {
    return this.circle ? u * u + v * v < 1 : Math.abs(u) < this.hx && Math.abs(v) < this.hy;
  }

  #gridCoord(i) {
    return -1 + (i + 0.5) * 2 / this.G;
  }

  #buildGrid(G) {
    this.G = G;
    const n = G * G;
    this.gamma = new Float32Array(n);
    this.g2 = new Float64Array(n);
    const inside = new Uint8Array(n);
    for (let j = 0; j < G; j++) {
      for (let i = 0; i < G; i++) {
        inside[j * G + i] = this.contains(this.#gridCoord(i), this.#gridCoord(j)) ? 1 : 0;
      }
    }
    this.edge = this.#rimCopies(inside);
    this.cr = new Float64Array(this.modes.length);
    this.ci = new Float64Array(this.modes.length);
    if (this.circle) this.#buildCircleGrid(inside);
    else this.#buildRectGrid();
  }

  // Cells just outside the edge copy the nearest plate cell, so interpolation at the rim never sees a fake zero.
  // Returns flat (outside cell, plate cell) pairs.
  #rimCopies(inside) {
    const { G } = this;
    const pairs = [];
    for (let j = 0; j < G; j++) {
      for (let i = 0; i < G; i++) {
        if (inside[j * G + i]) continue;
        let nearest = -1, nearestDist = 9;
        for (let dj = -2; dj <= 2; dj++) {
          for (let di = -2; di <= 2; di++) {
            const ii = i + di, jj = j + dj, dist = di * di + dj * dj;
            if (ii >= 0 && jj >= 0 && ii < G && jj < G && inside[jj * G + ii] && dist < nearestDist) {
              nearestDist = dist;
              nearest = jj * G + ii;
            }
          }
        }
        if (nearest >= 0) pairs.push(j * G + i, nearest);
      }
    }
    return Int32Array.from(pairs);
  }

  #buildCircleGrid(inside) {
    const { G } = this;
    const cells = [], radialIndex = [];
    inside.forEach((isInside, q) => {
      if (!isInside) return;
      cells.push(q);
      radialIndex.push(Math.hypot(this.#gridCoord(q % G), this.#gridCoord((q / G) | 0)) * (NR - 1));
    });
    this.cells = Int32Array.from(cells);
    this.rIdx = Float32Array.from(radialIndex);
    this.cosD = new Float32Array(cells.length);
    this.nTop = Math.max(...this.modes.map(mode => mode.n));
    this.Gre = new Float64Array((this.nTop + 1) * NR);
    this.Gim = new Float64Array((this.nTop + 1) * NR);
    this.act = new Uint8Array(this.nTop + 1);
  }

  #buildRectGrid() {
    const { G } = this;
    const { nbx, nby, atX, atY } = this.ritz;
    // the basis functions along one axis at every grid line on the plate, and which lines those are
    const tabulate = (halfWidth, basisAt, count) => {
      const table = new Float32Array(count * G);
      const lines = [];
      for (let i = 0; i < G; i++) {
        const u = this.#gridCoord(i);
        if (Math.abs(u) >= halfWidth) continue;
        lines.push(i);
        basisAt(toUnitInterval(u, halfWidth)).forEach((value, a) => { table[a * G + i] = value; });
      }
      return [table, Int32Array.from(lines)];
    };
    [this.Xg, this.cols] = tabulate(this.hx, atX, nbx);
    [this.Yg, this.rows] = tabulate(this.hy, atY, nby);
    this.Bre = new Float64Array(nbx * nby);
    this.Bim = new Float64Array(nbx * nby);
    this.Tre = new Float64Array(nbx * G);
    this.Tim = new Float64Array(nbx * G);
  }

  /** Point driver at normalised plate coordinates (u, v), kept just inside the edge. */
  setDrive(u, v) {
    if (this.circle) {
      const r = Math.hypot(u, v);
      if (r > DRIVER_LIMIT) {
        u *= DRIVER_LIMIT / r;
        v *= DRIVER_LIMIT / r;
      }
    } else {
      u = clamp(u, -DRIVER_LIMIT * this.hx, DRIVER_LIMIT * this.hx);
      v = clamp(v, -DRIVER_LIMIT * this.hy, DRIVER_LIMIT * this.hy);
    }
    this.du = u;
    this.dv = v;
    if (this.circle) this.#driveCircle(u, v);
    else this.#driveRect(u, v);
    this.patch?.setDrive(u, v);
  }

  #driveCircle(u, v) {
    const { G } = this;
    const r = Math.hypot(u, v);
    const x = r * (NR - 1);
    const j = Math.min(NR - 2, x | 0);
    const fraction = x - j;
    const cosTheta = r ? u / r : 1;
    const sinTheta = r ? v / r : 0;
    for (const mode of this.modes) mode.d0 = lerpAt(mode.R, j, fraction);
    this.cells.forEach((q, k) => {
      const cx = this.#gridCoord(q % G);
      const cy = this.#gridCoord((q / G) | 0);
      const cellR = Math.hypot(cx, cy);
      this.cosD[k] = cellR ? (cx * cosTheta + cy * sinTheta) / cellR : 1;
    });
  }

  #driveRect(u, v) {
    const { nbx, nby, atX, atY } = this.ritz;
    const X = atX(toUnitInterval(u, this.hx));
    const Y = atY(toUnitInterval(v, this.hy));
    for (const mode of this.modes) {
      let shape = 0;
      for (let i = 0; i < nbx; i++) {
        for (let j = 0; j < nby; j++) shape += mode.vec[i * nby + j] * X[i] * Y[j];
      }
      mode.d0 = shape;
    }
  }

  /**
   * Plate acceleration amplitude in units of g on the grid (this.gamma) for drives [{f: Hz, F}] at the driver.
   * Point force (F in N): W(x) = F Σ φ_k(x)φ_k(x₀) / (m(ω_k² − ω² + 2iζω_kω)) = F·G(x, x₀).
   * Shaker (F = coil force in N) whose moving mass m_a is bolted to the plate: the force reaching the plate is
   * F / (1 − m_a·ω²·G(x₀, x₀)), so resonances move from the free plate's (m_a → 0) towards the pinned plate's.
   * Rigid shaker (F = bolt acceleration in g, the m_a → ∞ limit): W(x) = W₀·G(x, x₀)/G(x₀, x₀), which resonates
   * where G(x₀, x₀) = 0, i.e. at the modes of the plate pinned at the bolt. Different frequencies add in power.
   */
  field(peaks, g) {
    this.g2.fill(0);
    this.driveCount = 0;
    for (const { f, F } of peaks) {
      const w = TAU * f;
      this.#nextDrive(w);
      if (this.patch) {
        this.patch.coefficients(w, F, g, this.cr, this.ci);
      } else {
        const [gr, gi] = this.#solveCoefficients(w, F, g);
        const [qr, qi] = this.#mountFactor(w, F, g, gr, gi);
        if (qr !== 1 || qi !== 0) this.#scaleCoefficients(qr, qi);
      }
      const tolerance = this.#strongestCoefficient() * 1e-3; // skip modes below 0.1 % of the strongest
      if (this.circle) this.#addCircleField(tolerance);
      else this.#addRectField(tolerance);
    }

    let max = 0;
    for (let i = 0; i < this.g2.length; i++) {
      const value = Math.sqrt(this.g2[i]);
      this.gamma[i] = value;
      if (value > max) max = value;
    }
    for (let k = 0; k < this.edge.length; k += 2) this.gamma[this.edge[k]] = this.gamma[this.edge[k + 1]];
    return max;
  }

  // With keepDrives set, each drive's complex field is kept for deflections().
  #nextDrive(omega) {
    if (!this.keepDrives) return;
    const cells = this.G * this.G;
    const store = this.driveStore ??= [];
    const drive = store[this.driveCount] ??= { omega, re: new Float32Array(cells), im: new Float32Array(cells) };
    if (drive.re.length !== cells) {
      drive.re = new Float32Array(cells);
      drive.im = new Float32Array(cells);
    } else {
      drive.re.fill(0);
      drive.im.fill(0);
    }
    drive.omega = omega;
    this.currentDrive = drive;
    this.driveCount++;
  }

  /** Each drive's complex deflection in metres from the last field() call, for a force multiplied by `gain`. */
  deflections(gain, g) {
    const drives = (this.driveStore ?? []).slice(0, this.driveCount);
    for (const { omega, re, im } of drives) {
      const scale = gain * g / (omega * omega);
      for (let q = 0; q < re.length; q++) {
        re[q] *= scale;
        im[q] *= scale;
      }
    }
    return drives;
  }

  // Sets the modal coefficients cr + i·ci for a point force and returns the response [re, im] at the driver.
  #solveCoefficients(w, F, g) {
    const { cr, ci, modes, zeta, mass, detune } = this;
    const scale = w * w * F / (mass * g);
    const forcePerMass = F / mass;
    let gr = 0, gi = 0;
    modes.forEach((mode, k) => {
      const wk = mode.w;
      let a = wk * wk - w * w;
      const b = 2 * (mode.zeta ?? zeta) * wk * w;
      if (detune) a += detune(k, this.driveCount, a, b, forcePerMass * Math.abs(mode.d0) * mode.inv);
      const amplitude = scale * mode.d0 * mode.inv / (a * a + b * b);
      cr[k] = amplitude * a;
      ci[k] = -amplitude * b;
      gr += cr[k] * mode.d0;
      gi += ci[k] * mode.d0;
    });
    return [gr, gi];
  }

  // The complex factor [re, im] by which the mounting scales every coefficient, from the response at the driver.
  #mountFactor(w, F, g, gr, gi) {
    if (this.mount === 'force') return [1, 0];
    // Modes above the cut-off ω_K still add flexibility at the driver. Up there a plate's point receptance tends to the
    // infinite plate's, a constant density 2/(π·zInf) per unit ω, so their sum is (1/(π·zInf·ω))·ln((ω_K+ω)/(ω_K−ω)).
    const wK = TAU * this.fTop;
    if (w < wK) gr += w * F / (g * Math.PI * this.zInf) * Math.log((wK + w) / (wK - w));

    if (this.mount === 'rigid') { // bolt point accelerates at exactly F g
      const d = gr * gr + gi * gi;
      if (d > 0) return [F * gr / d, -F * gi / d];
    }
    if (this.mount === 'shaker') { // q = 1/(1 − m_a ω² G₀), with G₀ = (gr + i·gi)·g/(ω²F) in m/N
      const zr = 1 - this.ma * gr * g / F;
      const zi = -this.ma * gi * g / F;
      const d = zr * zr + zi * zi;
      return [zr / d, -zi / d];
    }
    return [1, 0];
  }

  #scaleCoefficients(qr, qi) {
    const { cr, ci } = this;
    for (let k = 0; k < cr.length; k++) {
      const a = cr[k], b = ci[k];
      cr[k] = a * qr - b * qi;
      ci[k] = a * qi + b * qr;
    }
  }

  #strongestCoefficient() {
    const { cr, ci, modes } = this;
    let strongest = 0;
    modes.forEach((mode, k) => { strongest = Math.max(strongest, Math.hypot(cr[k], ci[k]) * mode.Rmax); });
    return strongest;
  }

  #addCircleField(tolerance) {
    const { modes, Gre, Gim, act, cr, ci, cells, rIdx, cosD, g2 } = this;
    const drive = this.keepDrives ? this.currentDrive : null;
    Gre.fill(0);
    Gim.fill(0);
    act.fill(0);
    let nHighest = -1;
    modes.forEach((mode, k) => { // G_n(ρ) = Σ_m c_nm R_nm(ρ)
      if (Math.hypot(cr[k], ci[k]) * mode.Rmax < tolerance) return;
      const offset = mode.n * NR, R = mode.R, a = cr[k], b = ci[k];
      for (let j = 0; j < NR; j++) {
        Gre[offset + j] += a * R[j];
        Gim[offset + j] += b * R[j];
      }
      act[mode.n] = 1;
      nHighest = Math.max(nHighest, mode.n);
    });
    for (let q = 0; q < cells.length; q++) { // W = Σ_n G_n(ρ) cos n(θ−θ₀), cosines by Chebyshev recurrence
      const x = rIdx[q];
      const j = x | 0;
      const fraction = x - j;
      const cosDelta = cosD[q];
      let previous = cosDelta, cosine = 1, wr = 0, wi = 0;
      for (let n = 0; n <= nHighest; n++) {
        if (act[n]) {
          const offset = n * NR + j;
          wr += cosine * lerpAt(Gre, offset, fraction);
          wi += cosine * lerpAt(Gim, offset, fraction);
        }
        const next = 2 * cosDelta * cosine - previous;
        previous = cosine;
        cosine = next;
      }
      g2[cells[q]] += wr * wr + wi * wi;
      if (drive) {
        drive.re[cells[q]] = wr;
        drive.im[cells[q]] = wi;
      }
    }
  }

  #addRectField(tolerance) {
    const { modes, cr, ci, Bre, Bim, Tre, Tim, Xg, Yg, rows, cols, g2, G } = this;
    const drive = this.keepDrives ? this.currentDrive : null;
    const { nbx, nby } = this.ritz;
    const nn = nbx * nby;
    Bre.fill(0);
    Bim.fill(0);
    modes.forEach((mode, k) => { // B_ab = Σ_k c_k a^k_ab, so W = Σ_ab B_ab X_a(u) Y_b(v)
      if (Math.hypot(cr[k], ci[k]) < tolerance) return;
      for (let t = 0; t < nn; t++) {
        Bre[t] += cr[k] * mode.vec[t];
        Bim[t] += ci[k] * mode.vec[t];
      }
    });
    for (const j of rows) {
      for (let a = 0; a < nbx; a++) {
        let sr = 0, si = 0;
        for (let b = 0; b < nby; b++) {
          const y = Yg[b * G + j];
          sr += Bre[a * nby + b] * y;
          si += Bim[a * nby + b] * y;
        }
        Tre[a * G + j] = sr;
        Tim[a * G + j] = si;
      }
      for (const i of cols) {
        let wr = 0, wi = 0;
        for (let a = 0; a < nbx; a++) {
          const x = Xg[a * G + i];
          wr += x * Tre[a * G + j];
          wi += x * Tim[a * G + j];
        }
        g2[j * G + i] += wr * wr + wi * wi;
        if (drive) {
          drive.re[j * G + i] = wr;
          drive.im[j * G + i] = wi;
        }
      }
    }
  }
}

// ---------- spectrum ----------

/**
 * Strongest local maxima of analyser dB bins in [fmin, fmax] above floorDb, refined by parabolic interpolation.
 * Peaks more than 55 dB below the strongest are dropped: the analyser's Blackman sidelobes sit 58 dB down.
 */
export function findPeaks(db, binHz, fmin, fmax, k, floorDb) {
  const peaks = [];
  const lo = Math.max(1, Math.ceil(fmin / binHz));
  const hi = Math.min(db.length - 2, Math.floor(fmax / binHz));
  for (let i = lo; i <= hi; i++) {
    const a = db[i - 1], b = db[i], c = db[i + 1];
    if (!(b > floorDb && b > a && b >= c)) continue;
    const curvature = a - 2 * b + c;
    const offset = Number.isFinite(curvature) && curvature < 0 ? 0.5 * (a - c) / curvature : 0;
    peaks.push({ f: (i + offset) * binHz, db: b - 0.25 * (a - c) * offset });
  }
  peaks.sort((p, q) => q.db - p.db);
  return peaks.filter(p => p.db > peaks[0].db - 55).slice(0, k);
}
