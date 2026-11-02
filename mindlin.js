// True Mindlin (first-order shear deformation) modes: w and the two rotations are solved together, instead of a thin-plate
// shape with a corrected frequency (physics.js, mindlinOmega). Mode objects match physics.js. Λ = L stays ωa²√(ρh/D) with the
// Kirchhoff D, and a = 1 throughout, so hOverA is the only thickness input.
//
// Energy per D/2, with s = S·a²/D = 6κ²(1−ν)/τ², τ = h/a, α = τ²/12:
//   U = ∫[ψx,x² + ψy,y² + 2νψx,xψy,y + ½(1−ν)(ψx,y + ψy,x)²] + s∫[(w,x + ψx)² + (w,y + ψy)²]
//   T = Λ²∫[w² + α(ψx² + ψy²)]
// Simply supported means the "hard" support w = 0, tangential rotation = 0, normal moment = 0: Mindlin's closed form exists
// only for it, and the soft support (rotation free) adds a boundary layer that smooth bases and the Kirchhoff limit don't need.

import { NR, bessel, gaussLegendre, rawBasis, symEig } from './physics.js';

const KAPPA2 = Math.PI ** 2 / 12; // makes the thickness-shear cut-off exact (Mindlin 1951)
const byFrequency = (p, q) => p.L - q.L;
const range = n => [...Array(n).keys()];
const shearRatio = (nu, hOverA, kappa2) => 6 * kappa2 * (1 - nu) / hOverA ** 2; // s above

// ---------- rectangles: three-field Rayleigh–Ritz ----------

// Vanishing order at an edge. w is only first order at a clamped edge: the slope there is carried by shear, unlike Kirchhoff.
const W_ORDER = { free: 0, simply: 1, clamped: 1 }; // also the rotation along the edge (tangential)
const NORMAL_ORDER = { free: 0, simply: 0, clamped: 1 }; // the rotation that bends across the edge

// Orthonormal polynomials on [0,1] with ends of order p, q (rawBasis), tabulated on the quadrature nodes with their slopes.
function fieldBasis(p, q, nb, [nodes, weights]) {
  const f = new Float64Array(nb), f1 = new Float64Array(nb), f2 = new Float64Array(nb);
  const nq = nodes.length;
  const V = range(nb).map(() => new Float64Array(nq));
  const V1 = range(nb).map(() => new Float64Array(nq));
  for (let j = 0; j < nq; j++) {
    rawBasis(p, q, nb, nodes[j], f, f1, f2);
    for (let k = 0; k < nb; k++) {
      V[k][j] = f[k];
      V1[k][j] = f1[k];
    }
  }
  const dot = (a, b) => a.reduce((total, x, j) => total + weights[j] * x * b[j], 0);
  const U = [], U1 = [], C = [];
  for (let k = 0; k < nb; k++) { // Gram–Schmidt, twice for stability, carrying slopes and the raw → orthonormal map along
    const u = V[k], d = V1[k], c = new Float64Array(nb);
    c[k] = 1;
    for (let pass = 0; pass < 2; pass++) {
      for (let j = 0; j < k; j++) {
        const projection = dot(u, U[j]);
        for (let i = 0; i < nq; i++) {
          u[i] -= projection * U[j][i];
          d[i] -= projection * U1[j][i];
        }
        for (let i = 0; i < nb; i++) c[i] -= projection * C[j][i];
      }
    }
    const scale = 1 / Math.sqrt(dot(u, u));
    for (const array of [u, d, c]) array.forEach((x, i) => { array[i] = x * scale; });
    U.push(u);
    U1.push(d);
    C.push(c);
  }
  const at = x => {
    rawBasis(p, q, nb, x, f, f1, f2);
    return C.map(ck => ck.reduce((total, c, i) => total + c * f[i], 0));
  };
  return { U, U1, at, dot };
}

// The bases of one axis with the integrals the energy needs. `w` also carries the tangential rotation; `n` carries the
// normal rotation and reaches the same polynomial degree as w, so ψ = −∇w is representable and the thin limit cannot lock.
function axisBases(lo, hi, nb) {
  const nn = nb + W_ORDER[lo] + W_ORDER[hi] - NORMAL_ORDER[lo] - NORMAL_ORDER[hi];
  const quadrature = gaussLegendre(nn + 12);
  const w = fieldBasis(W_ORDER[lo], W_ORDER[hi], nb, quadrature);
  const n = fieldBasis(NORMAL_ORDER[lo], NORMAL_ORDER[hi], nn, quadrature);
  const integrals = (A, B) => A.map(a => B.map(b => w.dot(a, b)));
  return {
    at: w.at,
    nb,
    nn,
    sym: lo === hi, // equal ends: both bases alternate even/odd about the centre
    ww: integrals(w.U1, w.U1), // ∫w_i′w_k′
    nn2: integrals(n.U1, n.U1), // ∫n_i′n_k′
    wn: integrals(w.U1, n.U), // ∫w_i′n_k
    nw: integrals(n.U1, w.U), // ∫n_i′w_k
  };
}

const parityFilter = (axis, parity) => i => !axis.sym || i % 2 === parity;

// Rigid motions a plate with these edges keeps: w = A + Bx + Cy needs w = 0 along each supported edge, and no slope at a clamped one.
const rigidCount = edges => {
  const held = edges.filter(edge => edge !== 'free');
  if (held.length === 0) return 3;
  return held.length === 1 && held[0] === 'simply' ? 1 : 0;
};

/**
 * Three-field Ritz modes of a rectangle with sides 1 (ξ) × beta (η), thickness hOverA of the ξ side. `edges` as in rectModes.
 * Returns { modes, nbx, nby, atX, atY } like rectModes: `vec` holds the w coefficients, and each mode has
 * ∫∫(w² + α(ψx² + ψy²))dξdη = 1 on the unit domain, so the modal sum in Plate stays correct. mode.L already includes
 * thickness: don't apply mindlinOmega to it. Below h/a ≈ 0.002 rounding reaches 0.05 %, and free edges' boundary layer
 * (width ≈ h/π) is under-resolved by up to 0.1 % at h/a = 0.005: use rectModes under h/a = 0.01. Takes 0.6–0.9 s at default nb.
 */
/** The w-basis evaluators of mindlinRectModes, for modes solved in a worker thread, where functions cannot cross. */
export function mindlinRectEvaluators(edges, beta, nb) {
  const [x0, x1, y0, y1] = typeof edges === 'string' ? [edges, edges, edges, edges] : edges;
  nb ??= MINDLIN_BASIS[(x0 === x1) + (y0 === y1)];
  const r = Math.sqrt(beta);
  return {
    atX: axisBases(x0, x1, Math.max(8, Math.round(nb / r))).at,
    atY: axisBases(y0, y1, Math.max(8, Math.round(nb * r))).at,
  };
}

export const MINDLIN_BASIS = [14, 18, 24]; // by symmetric axes

export function mindlinRectModes(edges, nu, beta, hOverA, { kappa2 = KAPPA2, nb } = {}) {
  const edgeList = typeof edges === 'string' ? [edges, edges, edges, edges] : edges;
  const [x0, x1, y0, y1] = edgeList;
  nb ??= MINDLIN_BASIS[(x0 === x1) + (y0 === y1)]; // 0.6–0.9 s; each symmetric axis halves the blocks of ≈ 3·nb² unknowns
  const r = Math.sqrt(beta);
  const X = axisBases(x0, x1, Math.max(8, Math.round(nb / r)));
  const Y = axisBases(y0, y1, Math.max(8, Math.round(nb * r)));

  const s = shearRatio(nu, hOverA, kappa2);
  const twist = (1 - nu) / 2;
  const b2 = beta * beta;
  const rot = Math.sqrt(12) / hOverA; // ψ = rot·ψ̃ turns the mass matrix into the identity, so symEig applies directly
  const rot2 = rot * rot;
  // stiffness blocks by field pair, for w = Xw_i·Yw_j, ψx = Xn_i·Yw_j, ψy = Xw_i·Yn_j
  const delta = (i, k) => (i === k ? 1 : 0);
  const block = [
    [(i, j, k, l) => s * (delta(j, l) * X.ww[i][k] + delta(i, k) * Y.ww[j][l] / b2),
      (i, j, k, l) => s * rot * delta(j, l) * X.wn[i][k],
      (i, j, k, l) => s * rot / beta * delta(i, k) * Y.wn[j][l]],
    [null,
      (i, j, k, l) => rot2 * (delta(j, l) * X.nn2[i][k] + delta(i, k) * (twist * Y.ww[j][l] / b2 + s * delta(j, l))),
      (i, j, k, l) => rot2 / beta * (nu * X.nw[i][k] * Y.nw[l][j] + twist * X.wn[k][i] * Y.wn[j][l])],
    [null, null,
      (i, j, k, l) => rot2 * (delta(i, k) * Y.nn2[j][l] / b2 + delta(j, l) * (twist * X.ww[i][k] + s * delta(i, k)))],
  ];

  const modes = [];
  const unknowns = (field, I, J) => I.flatMap(i => J.map(j => [field, i, j]));
  for (const px of X.sym ? [0, 1] : [-1]) {
    for (const py of Y.sym ? [0, 1] : [-1]) {
      // w even across an axis makes the rotation along that axis odd, and the other rotation even
      const Iw = range(X.nb).filter(parityFilter(X, px));
      const Jw = range(Y.nb).filter(parityFilter(Y, py));
      const In = range(X.nn).filter(parityFilter(X, 1 - px));
      const Jn = range(Y.nn).filter(parityFilter(Y, 1 - py));
      const dofs = [...unknowns(0, Iw, Jw), ...unknowns(1, In, Jw), ...unknowns(2, Iw, Jn)];
      const m = dofs.length;
      const K = new Float64Array(m * m);
      for (let a = 0; a < m; a++) {
        const [fa, i, j] = dofs[a];
        for (let b = a; b < m; b++) {
          const [fb, k, l] = dofs[b];
          K[a * m + b] = K[b * m + a] = block[fa][fb](i, j, k, l);
        }
      }
      const { values, vectors } = symEig(K, m);
      const nw = Iw.length * Jw.length;
      for (let e = 0; e < m; e++) {
        const vec = new Float64Array(X.nb * Y.nb);
        for (let a = 0; a < nw; a++) vec[dofs[a][1] * Y.nb + dofs[a][2]] = vectors[a * m + e];
        modes.push({ L: Math.sqrt(Math.max(values[e], 0)), vec, inv: 1, Rmax: 1, px, py });
      }
    }
  }
  modes.sort(byFrequency);

  // shortcut: keeps the lowest ~nbx·nby/4 like rectModes; raise nb for denser high-frequency patterns.
  const kept = modes.slice(0, Math.round(X.nb * Y.nb / 4));
  const rigid = rigidCount(edgeList); // counted, not thresholded: at small τ rounding lifts a rigid mode's Λ² above a long plate's first
  kept.forEach((mode, k) => {
    if (k < rigid) mode.L = 0;
    mode.label = k < rigid ? 'rigid' : `#${k - rigid + 1}`;
  });
  return { modes: kept, nbx: X.nb, nby: Y.nb, atX: X.at, atY: Y.at };
}

// ---------- circles: exact three-wavenumber solution ----------

// Miller's backward recurrence: [f_n(x), f_n′(x)] with relative accuracy at any x, where f = J (sign −1) or f = e^(−x)·I (sign +1).
// The trapezoid in physics.js's bessel() is only absolute, and its I_n overflows at the |k₃| ≈ π/τ of a thin plate.
function miller(n, x, sign) {
  if (!x) return [n ? 0 : 1, n === 1 ? 0.5 : 0];
  const big = Math.max(n, x);
  const top = sign < 0 ? Math.ceil(big + 15 + 12 * Math.cbrt(big)) : Math.ceil(n + 30 + 9 * Math.sqrt(x));
  let f = 1e-300, above = 0, sum = 0, fn = 0, below = 0, over = 0;
  for (let k = top; k >= 0; k--) {
    if (k === n + 1) over = f;
    if (k === n) fn = f;
    if (k === n - 1) below = f;
    if (sign > 0 || k % 2 === 0) sum += (k ? 2 : 1) * f; // J₀ + 2ΣJ₂ₖ = 1, I₀ + 2ΣIₖ = e^x
    const next = 2 * k / x * f + sign * above;
    above = f;
    f = next;
    if (Math.abs(f) > 1e250) {
      f *= 1e-250;
      above *= 1e-250;
      sum *= 1e-250;
      fn *= 1e-250;
      below *= 1e-250;
      over *= 1e-250;
    }
  }
  if (n === 0) below = sign * over;
  return [fn / sum, (below + sign * over) / (2 * sum)];
}

// [Z, dZ/dρ] of the radial solution for signed k²: J_n(kρ), or I_n(|k|ρ)/I_n(|k|) so a huge |k| cannot overflow.
function radial(n, k2, rho) {
  const k = Math.sqrt(Math.abs(k2));
  const x = k * rho;
  if (k2 > 0) {
    const [J, Jp] = miller(n, x, -1);
    return [J, k * Jp];
  }
  const [I, Ip] = miller(n, x, 1);
  const edge = miller(n, k, 1)[0];
  const scale = Math.exp(x - k) / edge;
  return [I * scale, k * Ip * scale];
}

// Signed squared wavenumbers: dilatational k₁² > 0, k₂² and rotational k₃² < 0 below the thickness-shear cut-off Λ_c² = 12s/τ².
function wavenumbers(L, { alpha, s, nu }) {
  const P = alpha * L * L, Q = L * L / s;
  const root = P + Q + Math.sqrt((P - Q) ** 2 + 4 * L * L);
  return [root / 2, 2 * (P * Q - L * L) / root, 2 * (P - s) / (1 - nu)]; // k₂² from k₁²k₂² = PQ − Λ², free of cancellation
}

// ψ = ∇Φ + curl H with Φ = Σ B_i Z_i(ρ)cos nθ, H = B₃Z₃(ρ)sin nθ and w = Σ σ_i B_i Z_i(ρ)cos nθ, σ_i = (P − s − k_i²)/s.
// Values at ρ = 1 per unit B: w, ψr, ψθ, M_r, M_rθ and s·(ψr + w,r) (shear scaled by s so a thin plate's rows stay comparable).
const [W, PR, PT, MR, MT, QR] = range(6);
function edgeColumns(n, L, wave, { alpha, s, nu }) {
  const P = alpha * L * L;
  const potential = k2 => {
    const [Z, dZ] = radial(n, k2, 1);
    return [(P - s - k2) / s * Z, dZ, -n * Z, -(1 - nu) * dZ + ((1 - nu) * n * n - k2) * Z, (1 - nu) * n * (Z - dZ), (P - k2) * dZ];
  };
  const [Z3, dZ3] = radial(n, wave[2], 1);
  return [potential(wave[0]), potential(wave[1]),
    [0, n * Z3, -dZ3, (1 - nu) * n * (dZ3 - Z3), (1 - nu) / 2 * (2 * dZ3 + (wave[2] - 2 * n * n) * Z3), s * n * Z3]];
}

// n = 0 splits off pure torsion (H alone, w = 0), which the plate's sand cannot see: only Φ's two columns remain.
function edgeMatrix(bc, n, L, ctx) {
  const rows = { clamped: [W, PR, PT], simply: [W, MR, PT], free: [MR, MT, QR] }[bc].filter(row => n || (row !== PT && row !== MT));
  const columns = edgeColumns(n, L, wavenumbers(L, ctx), ctx);
  return rows.map(row => columns.slice(0, rows.length).map(column => column[row]));
}

const det = ([a, b, c]) => (c
  ? a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])
  : a[0] * b[1] - a[1] * b[0]);

// The null vector at a root: from the first row when n = 0, else the best-conditioned pair of unit rows.
function nullVector(matrix) {
  const rows = matrix.map(row => {
    const size = Math.hypot(...row) || 1;
    return row.map(x => x / size);
  });
  if (rows.length === 2) return [-rows[0][1], rows[0][0], 0];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  return [cross(rows[0], rows[1]), cross(rows[1], rows[2]), cross(rows[0], rows[2])]
    .reduce((best, v) => (Math.hypot(...v) > Math.hypot(...best) ? v : best));
}

/** [w, ψr, ψθ] of a circular mode at ρ (w·cos nθ, ψr·cos nθ, ψθ·sin nθ). */
export function mindlinCircleFields(mode, rho) {
  const { n, B, wave, sigma } = mode;
  if (!mode.L) return [rho ** n, n ? -1 : 0, n ? 1 : 0]; // rigid: ψ = −∇w
  const [Z1, dZ1] = radial(n, wave[0], rho);
  const [Z2, dZ2] = radial(n, wave[1], rho);
  const [Z3, dZ3] = radial(n, wave[2], rho);
  const over = rho ? n / rho : 0;
  return [
    B[0] * sigma[0] * Z1 + B[1] * sigma[1] * Z2,
    B[0] * dZ1 + B[1] * dZ2 + B[2] * over * Z3,
    -over * (B[0] * Z1 + B[1] * Z2) - B[2] * dZ3,
  ];
}

const gaussCache = new Map();
const gauss = n => {
  if (!gaussCache.has(n)) gaussCache.set(n, gaussLegendre(n));
  return gaussCache.get(n);
};

// ∫(w² + α(ψr² + ψθ²))ρdρ, with a panel over the edge layer of width 1/|k₃|.
function modalNorm(mode, alpha) {
  const k1 = Math.sqrt(mode.wave[0]), k3 = Math.sqrt(Math.abs(mode.wave[2]));
  const split = Math.max(0.5, 1 - 12 / k3);
  let total = 0;
  for (const [lo, hi] of [[0, split], [split, 1]]) {
    const [nodes, weights] = gauss(Math.ceil(24 + (hi - lo) * (mode.n + k1 + (mode.wave[2] > 0 ? k3 : 0))));
    nodes.forEach((t, q) => {
      const rho = lo + (hi - lo) * t;
      const [w, pr, pt] = mindlinCircleFields(mode, rho);
      total += (hi - lo) * weights[q] * (w * w + alpha * (pr * pr + pt * pt)) * rho;
    });
  }
  return total;
}

const isNegligible = (value, Rmax) => Math.abs(value) < 1e-6 * Rmax;

function circleMode(n, L, wave, B, ctx) {
  const P = ctx.alpha * L * L;
  const sigma = wave.slice(0, 2).map(k2 => (P - ctx.s - k2) / ctx.s);
  const toUnit = 1 / (B[0] * sigma[0] || 1); // scale so the J_n(k₁ρ) term has coefficient 1, like J_n + C·I_n in physics.js
  const mode = { n, L, wave, sigma, B: B.map(x => x * toUnit) };

  const R = new Float32Array(NR);
  let Rmax = 0, nodalCircles = 0, sign = 0;
  for (let j = 0; j < NR; j++) {
    R[j] = mindlinCircleFields(mode, j / (NR - 1))[0];
    Rmax = Math.max(Rmax, Math.abs(R[j]));
  }
  for (let j = 1; j < NR - 1; j++) {
    if (isNegligible(R[j], Rmax)) continue;
    if (sign && Math.sign(R[j]) !== sign) nodalCircles++;
    sign = Math.sign(R[j]);
  }
  return {
    ...mode,
    m: nodalCircles,
    R,
    Rmax,
    inv: 1 / (Math.PI * (n ? 1 : 2) * modalNorm(mode, ctx.alpha)),
    label: `(${n},${nodalCircles})`,
  };
}

function rigidCircleMode(n, alpha) {
  const R = Float32Array.from({ length: NR }, (_, j) => (j / (NR - 1)) ** n);
  // n = 1 rocks, so rotary inertia adds to the plate's moment of inertia: ∫ρ³dρ + α·∫2ρdρ
  return { n, m: 0, L: 0, R, Rmax: 1, inv: n ? 1 / (Math.PI * (0.25 + alpha)) : 1 / Math.PI, label: 'rigid' };
}

// Bisection on [a, b], across which f changes sign.
function bisect(f, a, b) {
  let fa = f(a);
  for (let i = 0; i < 50; i++) {
    const c = (a + b) / 2, fc = f(c);
    if (fa * fc <= 0) {
      b = c;
    } else {
      a = c;
      fa = fc;
    }
  }
  return (a + b) / 2;
}

// Λ at which a wavenumber on each branch reaches k, from the dispersion αΛ⁴/s − (k²(α + 1/s) + 1)Λ² + k⁴ = 0 and k₃ = 2(αΛ² − s)/(1 − ν).
function branchFrequencies(k2, { alpha, s, nu }) {
  const a = alpha / s, b = k2 * (alpha + 1 / s) + 1;
  const root = b + Math.sqrt(b * b - 4 * a * k2 * k2);
  return [Math.sqrt(2 * k2 * k2 / root), Math.sqrt(root / (2 * a)), Math.sqrt((s + (1 - nu) * k2 / 2) / alpha)];
}

// Scan points: below the cut-off the evanescent k₂, k₃ don't oscillate, so k₁ moves by 0.2 between neighbours, against ~π between
// roots. Above it three families interleave and can nearly coincide, so every wavenumber moves by 0.05.
function scanGrid(n, lamMax, ctx) {
  const Lmax = lamMax * lamMax;
  const cutoff = Math.sqrt(ctx.s / ctx.alpha);
  const grid = [Lmax];
  for (let k = Math.max(0.5, 0.8 * n), more = true; more; k += 0.2) { // roots of order n lie above k = n, as in physics.js
    const [flexural] = branchFrequencies(k * k, ctx);
    grid.push(flexural);
    more = flexural < Math.min(Lmax, cutoff);
  }
  if (cutoff < Lmax) {
    grid.push(cutoff * (1 - 1e-9), cutoff * (1 + 1e-9));
    for (let k = 0.05, more = true; more; k += 0.05) {
      const branches = branchFrequencies(k * k, ctx);
      grid.push(...branches);
      more = Math.min(...branches) < Lmax;
    }
  }
  return { grid: grid.filter(L => L <= Lmax).sort((p, q) => p - q), cutoff };
}

/**
 * Exact Mindlin modes of a circular plate (radius 1) up to Λ = lamMax²; bc is 'free' | 'clamped' | 'simply' (hard: ψθ = 0).
 * Each mode has n, m (nodal circles of w), L, R (w at ρ = j/(NR−1)), Rmax, label, and inv including rotary inertia, so
 * Plate's modal sum applies unchanged; plus B, wave, sigma for mindlinCircleFields and mindlinNodalRadii.
 */
export function mindlinCircleModes(bc, nu, hOverA, lamMax, { kappa2 = KAPPA2 } = {}) {
  const ctx = { nu, alpha: hOverA * hOverA / 12, s: shearRatio(nu, hOverA, kappa2) };
  const modes = bc === 'free' ? [rigidCircleMode(0, ctx.alpha), rigidCircleMode(1, ctx.alpha)] : [];
  for (let n = 0, found = true; found; n++) {
    found = false;
    const determinant = L => det(edgeMatrix(bc, n, L, ctx));
    const { grid, cutoff } = scanGrid(n, lamMax, ctx);
    let lo = grid[0], fLo = determinant(lo);
    for (const hi of grid.slice(1)) {
      const fHi = determinant(hi);
      // k₂ and k₃ change regime at the cut-off: the determinant jumps there and isn't a root
      if (fLo * fHi < 0 && (lo - cutoff) * (hi - cutoff) > 0) {
        const L = bisect(determinant, lo, hi);
        modes.push(circleMode(n, L, wavenumbers(L, ctx), nullVector(edgeMatrix(bc, n, L, ctx)), ctx));
        found = true;
      }
      lo = hi;
      fLo = fHi;
    }
  }
  return modes.sort(byFrequency);
}

/** Radii r/a of a Mindlin circular mode's nodal circles, inner to outer; use instead of physics.js's nodalRadii. */
export function mindlinNodalRadii(mode) {
  if (!mode.L) return [];
  const { R, Rmax } = mode;
  const radii = [];
  for (let j = 1; j < NR - 2; j++) {
    if (isNegligible(R[j], Rmax) || Math.sign(R[j]) === Math.sign(R[j + 1])) continue;
    radii.push(bisect(rho => radialShapeAt(mode, rho), j / (NR - 1), (j + 1) / (NR - 1)));
  }
  return radii;
}

/**
 * A circular mode's radial shape at ρ, exact rather than from the sampled R, for every kind of circular mode: thin
 * (J_n + C·I_n), Mindlin, and a fluid-loaded mode, which is a mix of the vacuum modes of the same n.
 */
export function radialShapeAt(mode, rho) {
  if (mode.mix) return mode.mix.reduce((sum, { part, scale }) => sum + scale * radialShapeAt(part, rho), 0);
  if (mode.B) return mindlinCircleFields(mode, rho)[0];
  if (!mode.L) return rho ** mode.n;
  const [J, , I] = bessel(mode.n, Math.sqrt(mode.L) * rho);
  return J + mode.C * I;
}
