// Cymatix physics: Kirchhoff thin-plate eigenmodes and their steady-state response to point forces.
// Circular plates use the exact Bessel-function solutions; rectangular plates use a Rayleigh–Ritz
// solution on orthonormalised polynomials. Units are SI. Λ is the dimensionless frequency ω·L²·√(ρh/D).

const TAU = 2 * Math.PI;
export const NR = 257; // radial samples per circular mode shape, ρ = r/a from 0 to 1

/** [J_n(x), J_n'(x), I_n(x), I_n'(x)] for integer n ≥ 0, from Bessel's integrals. */
export function bessel(n, x) {
  // The trapezoid rule on these smooth periodic integrands converges exponentially (≈1e-15 here).
  const N = Math.ceil(n + Math.abs(x)) + 40, h = Math.PI / N;
  let J = 0, Jp = 0, I = 0, Ip = 0;
  for (let k = 0; k <= N; k++) {
    const t = k * h, s = Math.sin(t), c = Math.cos(t), w = k === 0 || k === N ? 0.5 : 1;
    const ph = n * t - x * s, e = w * Math.exp(x * c) * Math.cos(n * t);
    J += w * Math.cos(ph); Jp += w * s * Math.sin(ph); I += e; Ip += e * c;
  }
  return [J / N, Jp / N, I / N, Ip / N];
}

// Edge-condition rows at ρ = 1 for w = J_n(λρ) + C·I_n(λρ); a mode exists where their 2×2 determinant is zero.
function edgeRows(bc, n, lam, nu) {
  const [J, Jp, I, Ip] = bessel(n, lam), n2 = n * n, l2 = lam * lam;
  if (bc === 'clamped') return [J, I, Jp, Ip]; // w = 0, ∂w/∂r = 0
  const Jpp = -Jp / lam - (1 - n2 / l2) * J, Ipp = -Ip / lam + (1 + n2 / l2) * I; // from Bessel's equations
  const mJ = l2 * Jpp + nu * (lam * Jp - n2 * J), mI = l2 * Ipp + nu * (lam * Ip - n2 * I); // bending moment M_r
  if (bc === 'simply') return [J, I, mJ, mI]; // w = 0, M_r = 0
  const vJ = -l2 * lam * Jp - (1 - nu) * n2 * (lam * Jp - J), vI = l2 * lam * Ip - (1 - nu) * n2 * (lam * Ip - I); // Kirchhoff shear V_r
  return [mJ, mI, vJ, vI]; // free: M_r = 0, V_r = 0
}

function radialShape(n, lam, C) {
  const R = new Float32Array(NR);
  if (!lam) { for (let j = 0; j < NR; j++) R[j] = (j / (NR - 1)) ** n; return R; } // rigid translation / rocking
  const N = Math.ceil(n + lam) + 40, s = [], c = [], nt = [], wt = [];
  for (let k = 0; k <= N; k++) { const t = k * Math.PI / N; s.push(Math.sin(t)); c.push(Math.cos(t)); nt.push(n * t); wt.push(k === 0 || k === N ? 0.5 : 1); }
  for (let j = 0; j < NR; j++) {
    const x = lam * j / (NR - 1);
    let J = 0, I = 0;
    for (let k = 0; k <= N; k++) { J += wt[k] * Math.cos(nt[k] - x * s[k]); I += wt[k] * Math.exp(x * c[k]) * Math.cos(nt[k]); }
    R[j] = (J + C * I) / N;
  }
  return R;
}

function circleMode(n, lam, C) {
  const R = radialShape(n, lam, C);
  let norm = 0, Rmax = 0, m = 0, sign = 0;
  for (let j = 0; j < NR; j++) {
    norm += (j === 0 || j === NR - 1 ? 1 : j % 2 ? 4 : 2) * R[j] * R[j] * j / (NR - 1); // Simpson for ∫R²ρdρ
    Rmax = Math.max(Rmax, Math.abs(R[j]));
  }
  norm /= 3 * (NR - 1);
  for (let j = 1; j < NR - 1; j++) {
    if (Math.abs(R[j]) < 1e-6 * Rmax) continue;
    if (sign && Math.sign(R[j]) !== sign) m++; // nodal circles
    sign = Math.sign(R[j]);
  }
  // inv folds in ∫cos²(nθ)dθ so that a cos/sin degenerate pair sums to R(ρ)R(ρ₀)cos n(θ−θ₀)·inv
  return { n, m, L: lam * lam, C, R, Rmax, inv: 1 / (Math.PI * (n ? 1 : 2) * norm), label: lam ? `(${n},${m})` : 'rigid' };
}

/** Radii r/a of a circular mode's nodal circles: sign changes of the sampled shape, refined on J_n(λρ) + C·I_n(λρ). */
export function nodalRadii(m) {
  if (!m.L) return [];
  const lam = Math.sqrt(m.L), R = r => { const [J, , I] = bessel(m.n, lam * r); return J + m.C * I; }, out = [];
  for (let j = 1; j < NR - 2; j++) {
    if (Math.abs(m.R[j]) < 1e-6 * m.Rmax || Math.sign(m.R[j]) === Math.sign(m.R[j + 1])) continue;
    let a = j / (NR - 1), b = (j + 1) / (NR - 1), fa = R(a);
    for (let i = 0; i < 50; i++) { const c = (a + b) / 2, fc = R(c); if (fa * fc <= 0) b = c; else { a = c; fa = fc; } }
    out.push((a + b) / 2);
  }
  return out;
}

/** Exact modes of a circular plate (radius 1) up to λ = lamMax; edge bc is 'free' | 'clamped' | 'simply'. */
export function circleModes(bc, nu, lamMax) {
  const det = (n, l) => { const r = edgeRows(bc, n, l, nu); return r[0] * r[3] - r[1] * r[2]; };
  const modes = bc === 'free' ? [circleMode(0, 0, 0), circleMode(1, 0, 0)] : [];
  for (let n = 0, found = true; found && n < lamMax; n++) { // the first root grows with n, so stop at the first n with none
    found = false;
    // roots of order n lie above λ = n; starting at 0.8n also keeps the scan out of the range where J_n(λ) is rounding noise
    for (let lo = Math.max(0.5, 0.8 * n), flo = det(n, lo); lo < lamMax;) {
      const hi = lo + 0.1, fhi = det(n, hi); // roots of one n are ~π apart, so a 0.1 scan cannot skip one
      if (flo * fhi < 0) {
        let a = lo, b = hi, fa = flo;
        for (let i = 0; i < 50; i++) { const c = (a + b) / 2, fc = det(n, c); if (fa * fc <= 0) b = c; else { a = c; fa = fc; } }
        const lam = (a + b) / 2, r = edgeRows(bc, n, lam, nu);
        modes.push(circleMode(n, lam, -r[0] / r[1]));
        found = true;
      }
      lo = hi; flo = fhi;
    }
  }
  return modes.sort((p, q) => p.L - q.L);
}

function gaussLegendre(n) { // nodes and weights on [0, 1]
  const x = new Float64Array(n), w = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let z = Math.cos(Math.PI * (i + 0.75) / (n + 0.5)), pp = 1;
    for (let it = 0; it < 100; it++) {
      let p1 = 1, p2 = 0;
      for (let j = 1; j <= n; j++) { const p3 = p2; p2 = p1; p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j; }
      pp = n * (z * p1 - p2) / (z * z - 1);
      const dz = p1 / pp; z -= dz;
      if (Math.abs(dz) < 1e-15) break;
    }
    x[i] = (1 - z) / 2; w[i] = 1 / ((1 - z * z) * pp * pp);
  }
  return [x, w];
}

const ORDER = { free: 0, simply: 1, clamped: 2 }; // vanishing order of w at an edge: none, w = 0, w = w' = 0

// Raw 1D Ritz basis on ξ ∈ [0,1]: shifted Legendre polynomials times ξ^p·(1−ξ)^q, which enforces the geometric
// condition at each end (p, q from ORDER). Fills value, d/dξ, d²/dξ².
function rawBasis(p, q, nb, xi, f, f1, f2) {
  const t = 2 * xi - 1, a = xi, b = 1 - xi, pw = (x, k) => (k > 0 ? x ** k : 1);
  const w = pw(a, p) * pw(b, q);
  const w1 = (p ? p * pw(a, p - 1) * pw(b, q) : 0) - (q ? q * pw(a, p) * pw(b, q - 1) : 0);
  const w2 = (p > 1 ? p * (p - 1) * pw(a, p - 2) * pw(b, q) : 0) - (p && q ? 2 * p * q * pw(a, p - 1) * pw(b, q - 1) : 0)
    + (q > 1 ? q * (q - 1) * pw(a, p) * pw(b, q - 2) : 0);
  let p0 = 0, d0 = 0, s0 = 0, P = 1, d = 0, s = 0; // P_{k-1}, P_k and their t-derivatives
  for (let k = 0; k < nb; k++) {
    f[k] = w * P; f1[k] = w1 * P + 2 * w * d; f2[k] = w2 * P + 4 * w1 * d + 4 * w * s;
    const pn = ((2 * k + 1) * t * P - k * p0) / (k + 1);
    const dn = ((2 * k + 1) * (P + t * d) - k * d0) / (k + 1);
    const sn = ((2 * k + 1) * (2 * d + t * s) - k * s0) / (k + 1);
    p0 = P; d0 = d; s0 = s; P = pn; d = dn; s = sn;
  }
}

// 1D basis orthonormalised on [0,1] (Gram–Schmidt, twice for stability) with its stiffness integrals
// E11 = ∫U_i'U_k', E22 = ∫U_i''U_k'', E20 = ∫U_i''U_k, and an evaluator at any ξ.
function basis1D(p, q, nb) {
  const nq = nb + 12, [qx, qw] = gaussLegendre(nq);
  const f = new Float64Array(nb), f1 = new Float64Array(nb), f2 = new Float64Array(nb);
  const V = [], V1 = [], V2 = [];
  for (let k = 0; k < nb; k++) { V.push(new Float64Array(nq)); V1.push(new Float64Array(nq)); V2.push(new Float64Array(nq)); }
  for (let j = 0; j < nq; j++) { rawBasis(p, q, nb, qx[j], f, f1, f2); for (let k = 0; k < nb; k++) { V[k][j] = f[k]; V1[k][j] = f1[k]; V2[k][j] = f2[k]; } }
  const C = [], U = []; // row k of C maps raw → orthonormal
  for (let k = 0; k < nb; k++) {
    const u = Float64Array.from(V[k]), ck = new Float64Array(nb);
    ck[k] = 1;
    for (let pass = 0; pass < 2; pass++) for (let j = 0; j < k; j++) {
      let pr = 0;
      for (let i = 0; i < nq; i++) pr += qw[i] * u[i] * U[j][i];
      for (let i = 0; i < nq; i++) u[i] -= pr * U[j][i];
      for (let i = 0; i < nb; i++) ck[i] -= pr * C[j][i];
    }
    let nrm = 0;
    for (let i = 0; i < nq; i++) nrm += qw[i] * u[i] * u[i];
    nrm = 1 / Math.sqrt(nrm);
    for (let i = 0; i < nq; i++) u[i] *= nrm;
    for (let i = 0; i < nb; i++) ck[i] *= nrm;
    U.push(u); C.push(ck);
  }
  const mix = W => C.map(ck => { const o = new Float64Array(nq); for (let i = 0; i < nb; i++) for (let j = 0; j < nq; j++) o[j] += ck[i] * W[i][j]; return o; });
  const inner = (A, B) => A.map(a => B.map(b => { let t = 0; for (let j = 0; j < nq; j++) t += qw[j] * a[j] * b[j]; return t; }));
  const U1 = mix(V1), U2 = mix(V2);
  const at = xi => { rawBasis(p, q, nb, xi, f, f1, f2); return C.map(ck => ck.reduce((t, c, i) => t + c * f[i], 0)); };
  // equal ends: Legendre parity survives the symmetric weight, so functions alternate even/odd about ξ = ½
  return { at, E11: inner(U1, U1), E22: inner(U2, U2), E20: inner(U2, U), sym: p === q };
}

// Eigen-decomposition of a symmetric n×n matrix: Householder tridiagonalisation, then implicit QL
// (tred2/tql2 from JAMA, public domain). Eigenvectors are the columns of `vectors`.
function symEig(A, n) {
  const V = Float64Array.from(A), d = new Float64Array(n), e = new Float64Array(n);
  for (let j = 0; j < n; j++) d[j] = V[(n - 1) * n + j];
  for (let i = n - 1; i > 0; i--) {
    let scale = 0, h = 0;
    for (let k = 0; k < i; k++) scale += Math.abs(d[k]);
    if (scale === 0) {
      e[i] = d[i - 1];
      for (let j = 0; j < i; j++) { d[j] = V[(i - 1) * n + j]; V[i * n + j] = 0; V[j * n + i] = 0; }
    } else {
      for (let k = 0; k < i; k++) { d[k] /= scale; h += d[k] * d[k]; }
      let f = d[i - 1], g = Math.sqrt(h);
      if (f > 0) g = -g;
      e[i] = scale * g; h -= f * g; d[i - 1] = f - g;
      for (let j = 0; j < i; j++) e[j] = 0;
      for (let j = 0; j < i; j++) {
        f = d[j]; V[j * n + i] = f; g = e[j] + V[j * n + j] * f;
        for (let k = j + 1; k <= i - 1; k++) { g += V[k * n + j] * d[k]; e[k] += V[k * n + j] * f; }
        e[j] = g;
      }
      f = 0;
      for (let j = 0; j < i; j++) { e[j] /= h; f += e[j] * d[j]; }
      const hh = f / (h + h);
      for (let j = 0; j < i; j++) e[j] -= hh * d[j];
      for (let j = 0; j < i; j++) {
        f = d[j]; g = e[j];
        for (let k = j; k <= i - 1; k++) V[k * n + j] -= f * e[k] + g * d[k];
        d[j] = V[(i - 1) * n + j]; V[i * n + j] = 0;
      }
    }
    d[i] = h;
  }
  for (let i = 0; i < n - 1; i++) {
    V[(n - 1) * n + i] = V[i * n + i]; V[i * n + i] = 1;
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
  for (let j = 0; j < n; j++) { d[j] = V[(n - 1) * n + j]; V[(n - 1) * n + j] = 0; }
  V[(n - 1) * n + n - 1] = 1;
  for (let i = 1; i < n; i++) e[i - 1] = e[i];
  e[n - 1] = 0;
  let f = 0, tst1 = 0;
  const eps = 2 ** -52;
  for (let l = 0; l < n; l++) {
    tst1 = Math.max(tst1, Math.abs(d[l]) + Math.abs(e[l]));
    let m = l;
    while (m < n - 1 && Math.abs(e[m]) > eps * tst1) m++;
    if (m > l) {
      do {
        let g = d[l], p = (d[l + 1] - g) / (2 * e[l]), r = Math.hypot(p, 1);
        if (p < 0) r = -r;
        d[l] = e[l] / (p + r); d[l + 1] = e[l] * (p + r);
        const dl1 = d[l + 1];
        let h = g - d[l];
        for (let i = l + 2; i < n; i++) d[i] -= h;
        f += h;
        p = d[m];
        let c = 1, c2 = 1, c3 = 1, s = 0, s2 = 0;
        const el1 = e[l + 1];
        for (let i = m - 1; i >= l; i--) {
          c3 = c2; c2 = c; s2 = s;
          g = c * e[i]; h = c * p; r = Math.hypot(p, e[i]);
          e[i + 1] = s * r; s = e[i] / r; c = p / r; p = c * d[i] - s * g;
          d[i + 1] = h + s * (c * g + s * d[i]);
          for (let k = 0; k < n; k++) { h = V[k * n + i + 1]; V[k * n + i + 1] = s * V[k * n + i] + c * h; V[k * n + i] = c * V[k * n + i] - s * h; }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1; e[l] = s * p; d[l] = c * p;
      } while (Math.abs(e[l]) > eps * tst1);
    }
    d[l] += f; e[l] = 0;
  }
  return { values: d, vectors: V };
}

/**
 * Rayleigh–Ritz modes of a rectangular plate with sides 1 (ξ) × beta (η). `edges` is one condition for all four
 * edges or [ξ = 0, ξ = 1, η = 0, η = 1], each 'free' | 'simply' | 'clamped'.
 */
export function rectModes(edges, nu, beta, nb) {
  const [x0, x1, y0, y1] = (typeof edges === 'string' ? [edges, edges, edges, edges] : edges).map(e => ORDER[e]);
  const symAxes = (x0 === x1) + (y0 === y1);
  nb ??= [22, 28, 36][symAxes]; // each symmetric axis halves the block size, so it affords a richer basis
  // The same number of functions per unit length on both sides: an elongated plate gets more along its long side
  // (b/a = 8 needs ~100 there for its first 15 modes), keeping nbx·nby ≈ nb².
  const r = Math.sqrt(beta), nbx = Math.max(8, Math.round(nb / r)), nby = Math.max(8, Math.round(nb * r));
  const X = basis1D(x0, x1, nbx), Y = x0 === y0 && x1 === y1 && nbx === nby ? X : basis1D(y0, y1, nby);
  const b2 = beta * beta, b4 = b2 * b2, modes = [], range = n => [...Array(n).keys()];
  // each symmetric axis splits its functions into even and odd ones, which the stiffness never couples
  const classes = (B, n) => B.sym ? [range(n).filter(i => i % 2 === 0), range(n).filter(i => i % 2 === 1)] : [range(n)];
  for (const I of classes(X, nbx)) for (const J of classes(Y, nby)) {
    const idx = I.flatMap(i => J.map(j => [i, j])), m = idx.length, K = new Float64Array(m * m);
    for (let a = 0; a < m; a++) for (let b = 0; b < m; b++) {
      const [i, j] = idx[a], [k, l] = idx[b];
      K[a * m + b] = (j === l ? X.E22[i][k] : 0) + (i === k ? Y.E22[j][l] / b4 : 0)
        + (nu * (X.E20[i][k] * Y.E20[l][j] + X.E20[k][i] * Y.E20[j][l]) + 2 * (1 - nu) * X.E11[i][k] * Y.E11[j][l]) / b2;
    }
    const { values, vectors } = symEig(K, m);
    for (let e = 0; e < m; e++) {
      const vec = new Float64Array(nbx * nby);
      idx.forEach(([i, j], a) => { vec[i * nby + j] = vectors[a * m + e]; });
      // parity about the plate's centre lines: 0 even, 1 odd, −1 none (that axis has unequal ends)
      modes.push({ L: Math.sqrt(Math.max(values[e], 0)), vec, inv: 1, Rmax: 1, px: X.sym ? I[0] : -1, py: Y.sym ? J[0] : -1 });
    }
  }
  modes.sort((p, q) => p.L - q.L);
  // shortcut: only the lowest ~nbx·nby/4 Ritz modes are converged; raise nb for denser high-frequency patterns.
  const kept = modes.slice(0, Math.round(nbx * nby / 4));
  let elastic = 0;
  for (const m of kept) m.label = m.L < 1e-3 ? 'rigid' : `#${++elastic}`;
  return { modes: kept, nbx, nby, atX: X.at, atY: Y.at };
}

/**
 * Mindlin correction (shear deformation + rotary inertia) of a Kirchhoff angular frequency wK. The mode's bending
 * wavenumber k (k⁴ = wK²ρh/D) is kept and Mindlin's dispersion relation
 *   ρh·I·ω⁴ − (S·I·k² + ρh·D·k² + ρh·S)·ω² + S·D·k⁴ = 0,  S = κ²Gh, I = ρh³/12,
 * is solved for the flexural root. Exact for simply supported rectangles; a local approximation for other edges.
 */
export function mindlinOmega(wK, { E, rho, nu, h }, kappa2 = Math.PI ** 2 / 12) {
  if (!wK) return 0;
  const D = E * h ** 3 / (12 * (1 - nu * nu)), rhoH = rho * h, k2 = wK * Math.sqrt(rhoH / D);
  const S = kappa2 * E / (2 * (1 + nu)) * h, I = rho * h ** 3 / 12;
  const a = rhoH * I, b = S * I * k2 + rhoH * D * k2 + rhoH * S, c = S * D * k2 * k2;
  return Math.sqrt(2 * c / (b + Math.sqrt(b * b - 4 * a * c))); // lower root, written to avoid cancellation
}

/** A plate with its modes, a G×G field grid over [-1,1]², a point driver, and the plate-acceleration field. */
export class Plate {
  // cfg: { shape: 'circle'|'rect', bc, E, rho, nu, h, d, Lx, Ly (all SI), zeta, fmax (Hz), G,
  //        mount: 'force'|'shaker'|'rigid', ma (shaker moving mass, kg), mindlin (default true) }
  configure(cfg) {
    const circle = cfg.shape === 'circle';
    const D = cfg.E * cfg.h ** 3 / (12 * (1 - cfg.nu * cfg.nu)), rhoH = cfg.rho * cfg.h, L = circle ? cfg.d / 2 : cfg.Lx;
    this.circle = circle;
    this.wScale = Math.sqrt(D / rhoH) / (L * L); // ω = Λ·wScale
    this.zInf = 8 * Math.sqrt(D * rhoH); // an infinite plate's point impedance F/v, used for the modes above the cut-off
    this.mass = circle ? rhoH * L * L : rhoH * cfg.Lx * cfg.Ly; // modal mass scale for unit-normalised shapes
    this.zeta = cfg.zeta;
    this.mount = cfg.mount ?? 'force'; this.ma = cfg.ma ?? 0;
    let key;
    if (circle) {
      // modes reach 30 % past the top analysed frequency; λ is rounded up so small edits reuse them
      this.lamMax = Math.min(50, Math.max(10, 5 * Math.ceil(Math.sqrt(1.3 * TAU * cfg.fmax / this.wScale) / 5)));
      key = `c|${cfg.bc}|${cfg.nu}|${this.lamMax}`;
      if (key !== this.key) this.modes = circleModes(cfg.bc, cfg.nu, this.lamMax);
    } else {
      const edges = cfg.edges ?? [cfg.bc, cfg.bc, cfg.bc, cfg.bc];
      key = `r|${edges}|${cfg.nu}|${cfg.Ly / cfg.Lx}`;
      this.edges = edges;
      if (key !== this.key) { this.ritz = rectModes(edges, cfg.nu, cfg.Ly / cfg.Lx); this.modes = this.ritz.modes; }
    }
    const thick = w => cfg.mindlin === false ? w : mindlinOmega(w, cfg);
    for (const m of this.modes) { m.w = thick(m.L * this.wScale); m.f = m.w / TAU; }
    this.fTop = thick((circle ? this.lamMax ** 2 : this.modes.at(-1).L) * this.wScale) / TAU;
    const s = Math.max(cfg.Lx, cfg.Ly);
    this.hx = circle ? 1 : cfg.Lx / s; this.hy = circle ? 1 : cfg.Ly / s;
    if (key !== this.key || cfg.G !== this.G) { this.key = key; this.buildGrid(cfg.G); }
    this.setDrive(this.du ?? 0, this.dv ?? 0);
  }

  contains(u, v) { return this.circle ? u * u + v * v < 1 : Math.abs(u) < this.hx && Math.abs(v) < this.hy; }

  buildGrid(G) {
    this.G = G;
    const n = G * G, at = i => -1 + (i + 0.5) * 2 / G, inside = new Uint8Array(n);
    this.at = at;
    this.gamma = new Float32Array(n); this.g2 = new Float64Array(n);
    for (let j = 0; j < G; j++) for (let i = 0; i < G; i++) inside[j * G + i] = this.contains(at(i), at(j)) ? 1 : 0;
    // Cells just outside the edge copy the nearest plate cell, so interpolation at the rim never sees a fake zero.
    const edge = [];
    for (let j = 0; j < G; j++) for (let i = 0; i < G; i++) {
      if (inside[j * G + i]) continue;
      let best = -1, bd = 9;
      for (let dj = -2; dj <= 2; dj++) for (let di = -2; di <= 2; di++) {
        const ii = i + di, jj = j + dj, d = di * di + dj * dj;
        if (ii >= 0 && jj >= 0 && ii < G && jj < G && inside[jj * G + ii] && d < bd) { bd = d; best = jj * G + ii; }
      }
      if (best >= 0) edge.push(j * G + i, best);
    }
    this.edge = Int32Array.from(edge);
    this.cr = new Float64Array(this.modes.length); this.ci = new Float64Array(this.modes.length);
    if (this.circle) {
      const cells = [], rIdx = [];
      for (let q = 0; q < n; q++) if (inside[q]) { cells.push(q); rIdx.push(Math.hypot(at(q % G), at((q / G) | 0)) * (NR - 1)); }
      this.cells = Int32Array.from(cells); this.rIdx = Float32Array.from(rIdx); this.cosD = new Float32Array(cells.length);
      this.nTop = Math.max(...this.modes.map(m => m.n));
      this.Gre = new Float64Array((this.nTop + 1) * NR); this.Gim = new Float64Array((this.nTop + 1) * NR); this.act = new Uint8Array(this.nTop + 1);
    } else {
      const { nbx, nby, atX, atY } = this.ritz, cols = [], rows = [];
      this.Xg = new Float32Array(nbx * G); this.Yg = new Float32Array(nby * G);
      for (let i = 0; i < G; i++) {
        const u = at(i);
        if (Math.abs(u) < this.hx) { cols.push(i); atX((u / this.hx + 1) / 2).forEach((b, a) => { this.Xg[a * G + i] = b; }); }
        if (Math.abs(u) < this.hy) { rows.push(i); atY((u / this.hy + 1) / 2).forEach((b, a) => { this.Yg[a * G + i] = b; }); }
      }
      this.cols = Int32Array.from(cols); this.rows = Int32Array.from(rows);
      this.Bre = new Float64Array(nbx * nby); this.Bim = new Float64Array(nbx * nby);
      this.Tre = new Float64Array(nbx * G); this.Tim = new Float64Array(nbx * G);
    }
  }

  /** Point driver at normalised plate coordinates (u, v), kept just inside the edge. */
  setDrive(u, v) {
    if (this.circle) { const r = Math.hypot(u, v); if (r > 0.97) { u *= 0.97 / r; v *= 0.97 / r; } }
    else { u = Math.max(-0.97 * this.hx, Math.min(0.97 * this.hx, u)); v = Math.max(-0.97 * this.hy, Math.min(0.97 * this.hy, v)); }
    this.du = u; this.dv = v;
    if (this.circle) {
      const r = Math.hypot(u, v), x = r * (NR - 1), j = Math.min(NR - 2, x | 0), fr = x - j, cu = r ? u / r : 1, sv = r ? v / r : 0, G = this.G;
      for (const m of this.modes) m.d0 = m.R[j] + fr * (m.R[j + 1] - m.R[j]);
      this.cells.forEach((q, k) => { const cx = this.at(q % G), cy = this.at((q / G) | 0), rr = Math.hypot(cx, cy); this.cosD[k] = rr ? (cx * cu + cy * sv) / rr : 1; });
    } else {
      const { nbx, nby, atX, atY } = this.ritz, X = atX((u / this.hx + 1) / 2), Y = atY((v / this.hy + 1) / 2);
      for (const m of this.modes) { let s = 0; for (let i = 0; i < nbx; i++) for (let j = 0; j < nby; j++) s += m.vec[i * nby + j] * X[i] * Y[j]; m.d0 = s; }
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
    const { cr, ci, modes } = this;
    for (const { f, F } of peaks) {
      const w = TAU * f, s = w * w * F / (this.mass * g);
      let cmax = 0, gr = 0, gi = 0;
      modes.forEach((m, k) => {
        const wk = m.w, a = wk * wk - w * w, b = 2 * this.zeta * wk * w, amp = s * m.d0 * m.inv / (a * a + b * b);
        cr[k] = amp * a; ci[k] = -amp * b;
        gr += cr[k] * m.d0; gi += ci[k] * m.d0; // response at the driver
      });
      // Modes above the cut-off ω_K still add flexibility at the driver. Up there a plate's point receptance tends to the
      // infinite plate's, a constant density 2/(π·zInf) per unit ω, so their sum is (1/(π·zInf·ω))·ln((ω_K+ω)/(ω_K−ω)).
      const wK = TAU * this.fTop;
      if (this.mount !== 'force' && w < wK) gr += w * F / (g * Math.PI * this.zInf) * Math.log((wK + w) / (wK - w));
      let qr = 1, qi = 0, d = gr * gr + gi * gi;
      if (this.mount === 'rigid' && d > 0) { qr = F * gr / d; qi = -F * gi / d; } // bolt point accelerates at exactly F g
      if (this.mount === 'shaker') { // q = 1/(1 − m_a ω² G₀), with G₀ = (gr + i·gi)·g/(ω²F) in m/N
        const zr = 1 - this.ma * gr * g / F, zi = -this.ma * gi * g / F;
        d = zr * zr + zi * zi; qr = zr / d; qi = -zi / d;
      }
      if (qr !== 1 || qi !== 0) for (let k = 0; k < modes.length; k++) { const a = cr[k], b = ci[k]; cr[k] = a * qr - b * qi; ci[k] = a * qi + b * qr; }
      modes.forEach((m, k) => { cmax = Math.max(cmax, Math.hypot(cr[k], ci[k]) * m.Rmax); });
      if (this.circle) this.circleField(cmax * 1e-3); else this.rectField(cmax * 1e-3); // skip modes below 0.1 % of the strongest
    }
    let max = 0;
    for (let i = 0; i < this.g2.length; i++) { const v = Math.sqrt(this.g2[i]); this.gamma[i] = v; if (v > max) max = v; }
    for (let k = 0; k < this.edge.length; k += 2) this.gamma[this.edge[k]] = this.gamma[this.edge[k + 1]];
    return max;
  }

  circleField(tol) {
    const { modes, Gre, Gim, act, cr, ci, cells, rIdx, cosD, g2 } = this;
    Gre.fill(0); Gim.fill(0); act.fill(0);
    let nHi = -1;
    modes.forEach((m, k) => { // G_n(ρ) = Σ_m c_nm R_nm(ρ)
      if (Math.hypot(cr[k], ci[k]) * m.Rmax < tol) return;
      const o = m.n * NR, R = m.R, a = cr[k], b = ci[k];
      for (let j = 0; j < NR; j++) { Gre[o + j] += a * R[j]; Gim[o + j] += b * R[j]; }
      act[m.n] = 1; nHi = Math.max(nHi, m.n);
    });
    for (let q = 0; q < cells.length; q++) { // W = Σ_n G_n(ρ) cos n(θ−θ₀), cosines by Chebyshev recurrence
      const x = rIdx[q], j = x | 0, fr = x - j, cd = cosD[q];
      let tp = cd, t = 1, wr = 0, wi = 0;
      for (let n = 0; n <= nHi; n++) {
        if (act[n]) { const o = n * NR + j; wr += t * (Gre[o] + fr * (Gre[o + 1] - Gre[o])); wi += t * (Gim[o] + fr * (Gim[o + 1] - Gim[o])); }
        const tn = 2 * cd * t - tp; tp = t; t = tn;
      }
      g2[cells[q]] += wr * wr + wi * wi;
    }
  }

  rectField(tol) {
    const { modes, cr, ci, Bre, Bim, Tre, Tim, Xg, Yg, rows, cols, g2, G } = this, { nbx, nby } = this.ritz, nn = nbx * nby;
    Bre.fill(0); Bim.fill(0);
    modes.forEach((m, k) => { // B_ab = Σ_k c_k a^k_ab, so W = Σ_ab B_ab X_a(u) Y_b(v)
      if (Math.hypot(cr[k], ci[k]) < tol) return;
      for (let t = 0; t < nn; t++) { Bre[t] += cr[k] * m.vec[t]; Bim[t] += ci[k] * m.vec[t]; }
    });
    for (let r = 0; r < rows.length; r++) {
      const j = rows[r];
      for (let a = 0; a < nbx; a++) {
        let sr = 0, si = 0;
        for (let b = 0; b < nby; b++) { const y = Yg[b * G + j]; sr += Bre[a * nby + b] * y; si += Bim[a * nby + b] * y; }
        Tre[a * G + j] = sr; Tim[a * G + j] = si;
      }
      for (let c = 0; c < cols.length; c++) {
        const i = cols[c];
        let wr = 0, wi = 0;
        for (let a = 0; a < nbx; a++) { const x = Xg[a * G + i]; wr += x * Tre[a * G + j]; wi += x * Tim[a * G + j]; }
        g2[j * G + i] += wr * wr + wi * wi;
      }
    }
  }
}

/**
 * Strongest local maxima of analyser dB bins in [fmin, fmax] above floorDb, refined by parabolic interpolation.
 * Peaks more than 55 dB below the strongest are dropped: the analyser's Blackman sidelobes sit 58 dB down.
 */
export function findPeaks(db, binHz, fmin, fmax, k, floorDb) {
  const out = [], lo = Math.max(1, Math.ceil(fmin / binHz)), hi = Math.min(db.length - 2, Math.floor(fmax / binHz));
  for (let i = lo; i <= hi; i++) {
    const a = db[i - 1], b = db[i], c = db[i + 1];
    if (!(b > floorDb && b > a && b >= c)) continue;
    const den = a - 2 * b + c, d = Number.isFinite(den) && den < 0 ? 0.5 * (a - c) / den : 0;
    out.push({ f: (i + d) * binHz, db: b - 0.25 * (a - c) * d });
  }
  out.sort((p, q) => q.db - p.db);
  return out.filter(p => p.db > out[0].db - 55).slice(0, k);
}
