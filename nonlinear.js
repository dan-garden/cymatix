// Cymatix large-amplitude plate vibration: von Kármán membrane stiffening reduced to one cubic per mode,
//   q̈ + 2ζω q̇ + ω²q + γq³ = f cos ωt,   γ = κ·ω²/h²  (ω the Kirchhoff frequency, w = q·φ, κ dimensionless),
// solved by harmonic balance: ω² → ω² + ¾γ|q|². κ depends only on the mode shape, ν and the in-plane edge condition.
//
// γ comes from the membrane energy left after the in-plane displacements relax, found from the Airy function Φ,
// ∇⁴Φ = −(E/2)[w,w], in the complementary-energy form (Reissner). For w = qφ and Φ = Eq²Ψ, Ψ maximises
//   R(Ψ) − ½B(Ψ,Ψ),   R(ψ) = ∫ ½ψ,yy φ,x² + ½ψ,xx φ,y² − ψ,xy φ,x φ,y dA,
//   B(Ψ,ψ) = ∫ ∇²Ψ∇²ψ − (1+ν)(Ψ,xxψ,yy + Ψ,yyψ,xx − 2Ψ,xyψ,xy) dA,
// and ¼γq⁴·(modal mass) = (Ehq⁴/2L²)·R(Ψ). Displacement conditions are natural in this form: an unconstrained Ψ is a
// fully fixed edge (u = v = 0), Ψ = Ψ,n = 0 a stress-free edge. B is the plate stiffness with ν → −ν (Thomas & Bilbao 2008).

import { gaussLegendre, NR, rawBasis } from './physics.js';

/**
 * In-plane condition assumed for each transverse edge condition. A free edge is stress-free by definition. A knife-edge or
 * pinned support only reacts transverse force, so it is taken as stress-free too. A clamp grips the plate, so it is taken as
 * fully fixed; a real clamp lies below this limit, so pass 'free' for the other bound.
 */
export const DEFAULT_IN_PLANE = { free: 'free', simply: 'free', clamped: 'fixed' };

const inPlaneOf = (bc, override) => (bc === 'free' ? 'free' : override ?? DEFAULT_IN_PLANE[bc]);

// Cholesky of the symmetric band matrix in a (row i, offsets 0..w from the diagonal, a[i*(w+1)+d]), in place; L·Lᵀ.
function bandCholesky(a, n, w) {
  const stride = w + 1;
  for (let i = 0; i < n; i++) {
    for (let d = 0; d <= w && i + d < n; d++) {
      let sum = a[i * stride + d];
      for (let k = 1; k <= w - d && i - k >= 0; k++) sum -= a[(i - k) * stride + k] * a[(i - k) * stride + k + d];
      a[i * stride + d] = d === 0 ? Math.sqrt(sum) : sum / a[i * stride];
    }
  }
}

function bandSolve(a, n, w, b) {
  const stride = w + 1;
  for (let i = 0; i < n; i++) {
    let sum = b[i];
    for (let k = 1; k <= w && i - k >= 0; k++) sum -= a[(i - k) * stride + k] * b[i - k];
    b[i] = sum / a[i * stride];
  }
  for (let i = n - 1; i >= 0; i--) {
    let sum = b[i];
    for (let k = 1; k <= w && i + k < n; k++) sum -= a[i * stride + k] * b[i + k];
    b[i] = sum / a[i * stride];
  }
}

// ---------- circular plates ----------

// Sixth-order first-derivative weights ×60 on samples 0..6, evaluated at samples 0, 1 and 2.
const EDGE_WEIGHTS = [
  [-147, 360, -450, 400, -225, 72, -10],
  [-10, -77, 150, -100, 50, -15, 2],
  [2, -24, -35, 80, -30, 8, -1],
];

// R' of a uniformly sampled smooth shape: sixth-order differences, central inside and one-sided near the ends.
function radialSlope(R, slope) {
  const s = R.length - 1;
  for (let j = 3; j <= s - 3; j++) {
    slope[j] = (45 * (R[j + 1] - R[j - 1]) - 9 * (R[j + 2] - R[j - 2]) + (R[j + 3] - R[j - 3])) * s / 60;
  }
  for (let j = 0; j < 3; j++) {
    let front = 0, back = 0;
    for (let i = 0; i < 7; i++) {
      front += EDGE_WEIGHTS[j][i] * R[i];
      back += EDGE_WEIGHTS[j][i] * R[s - i];
    }
    slope[j] = front * s / 60;
    slope[s - j] = -back * s / 60;
  }
}

/** Hermite cubic element (nodes at the ends, length l): value, slope and curvature of its four shape functions at x ∈ [0,1]. */
function hermite(x, l, v, d1, d2) {
  const x2 = x * x, x3 = x2 * x;
  v[0] = 1 - 3 * x2 + 2 * x3;
  v[1] = l * (x - 2 * x2 + x3);
  v[2] = 3 * x2 - 2 * x3;
  v[3] = l * (x3 - x2);
  d1[0] = (-6 * x + 6 * x2) / l;
  d1[1] = 1 - 4 * x + 3 * x2;
  d1[2] = (6 * x - 6 * x2) / l;
  d1[3] = 3 * x2 - 2 * x;
  d2[0] = (-6 + 12 * x) / (l * l);
  d2[1] = (-4 + 6 * x) / l;
  d2[2] = (6 - 12 * x) / (l * l);
  d2[3] = (6 * x - 2) / l;
}

const SAMPLES_PER_ELEMENT = 2; // Hermite elements span this many intervals of the radial samples; kept even for Simpson's rule

/**
 * The factorised radial stiffness of Ψ = F(ρ)cos pθ, on Hermite elements with F(0) = F'(0) = 0 (regular at the centre; p is
 * 0 or even). Stress-free edge: F(1) = F'(1) = 0. Fixed edge: nothing, the natural condition is u = v = 0.
 */
function radialFactor(p, nu, inPlane, elements) {
  // nodal (F, F') without the centre node; at a free edge F' = 0 and, unless p = 0 where F(0) = 0 already fixes the gauge, F = 0
  const dofs = 2 * elements - (inPlane === 'free' ? (p ? 2 : 1) : 0);
  const w = 3;
  const K = new Float64Array(dofs * (w + 1));
  const [gx, gw] = gaussLegendre(5);
  const v = new Float64Array(4), d1 = new Float64Array(4), d2 = new Float64Array(4);
  const l = 1 / elements;
  const angular = p ? Math.PI : 2 * Math.PI;
  const radial = new Float64Array(4), hoop = new Float64Array(4), twist = new Float64Array(4);
  for (let e = 0; e < elements; e++) {
    for (let g = 0; g < 5; g++) {
      const rho = (e + gx[g]) * l;
      hermite(gx[g], l, v, d1, d2);
      for (let a = 0; a < 4; a++) {
        radial[a] = d1[a] / rho - p * p * v[a] / (rho * rho); // σ_rr
        twist[a] = p * (d1[a] / rho - v[a] / (rho * rho)); // σ_rθ
        hoop[a] = d2[a]; // σ_θθ
      }
      const weight = gw[g] * l * rho * angular;
      for (let a = 0; a < 4; a++) {
        const i = 2 * e + a - 2;
        if (i < 0 || i >= dofs) continue;
        for (let c = a; c < 4; c++) {
          const j = 2 * e + c - 2;
          if (j < 0 || j >= dofs) continue;
          // ∇²F∇²G − (1+ν)(F''σ_G + G''σ_F − 2 p²(F/ρ)'(G/ρ)')
          const entry = (hoop[a] + radial[a]) * (hoop[c] + radial[c])
            - (1 + nu) * (hoop[a] * radial[c] + hoop[c] * radial[a] - 2 * twist[a] * twist[c]);
          K[i * (w + 1) + j - i] += weight * entry;
        }
      }
    }
  }
  bandCholesky(K, dofs, w);
  return { K, dofs, w };
}

// R(F) for the harmonic p of w = R(ρ)cos nθ, assembled into load from samples of R and R' (Simpson within each element),
// using the θ-integrals of cos²nθ, sin²nθ and sin nθ cos nθ with cos pθ and sin pθ.
function radialLoad(p, n, R, slope, elements, dofs, load) {
  load.fill(0);
  const v = new Float64Array(4), d1 = new Float64Array(4), d2 = new Float64Array(4);
  const l = 1 / elements;
  const s = R.length - 1;
  const span = s / elements;
  const simpson = (i, m) => (i === 0 || i === m ? 1 : i % 2 ? 4 : 2) / (3 * m);
  for (let e = 0; e < elements; e++) {
    for (let i = 0; i <= span; i++) {
      const j = e * span + i;
      if (j === 0) continue; // every integrand vanishes at the centre for the regular F chosen
      const rho = j / s;
      const dR = slope[j], Rj = R[j];
      hermite(i / span, l, v, d1, d2);
      const weight = simpson(i, span) * l;
      for (let a = 0; a < 4; a++) {
        const dof = 2 * e + a - 2;
        if (dof < 0 || dof >= dofs) continue;
        let term;
        if (n === 0) {
          term = Math.PI * d1[a] * dR * dR; // only p = 0 exists: ∫F'R'²
        } else if (p === 0) {
          term = (Math.PI / 2) * (d1[a] * dR * dR + d2[a] * n * n * Rj * Rj / rho);
        } else {
          const radial = d1[a] - p * p * v[a] / rho;
          const twist = d1[a] / rho - v[a] / (rho * rho);
          term = Math.PI * (0.25 * radial * dR * dR - 0.25 * d2[a] * n * n * Rj * Rj / rho - 0.5 * p * n * twist * dR * Rj);
        }
        load[dof] += weight * term;
      }
    }
  }
}

/**
 * κ of each circular mode from circleModes, for w = q·R(ρ)cos nθ (the plate's modal coordinate). The plate's peak deflection is
 * A = |q|·Rmax, and in terms of it the backbone is ω_nl² = ω²(1 + ¾κ(A/h)²/Rmax²).
 * Cost: see the report; no Worker needed.
 */
export function circleKappas(modes, bc, nu, inPlane) {
  const type = inPlaneOf(bc, inPlane);
  const elements = (NR - 1) / SAMPLES_PER_ELEMENT;
  const factors = new Map();
  const factor = p => {
    if (!factors.has(p)) factors.set(p, radialFactor(p, nu, type, elements));
    return factors.get(p);
  };
  const slope = new Float64Array(NR);
  const load = new Float64Array(2 * elements);
  const copy = new Float64Array(2 * elements);
  const kappa = new Float64Array(modes.length);
  modes.forEach((mode, k) => {
    if (!mode.L) return; // rigid modes do not stretch
    radialSlope(mode.R, slope);
    let R = 0;
    for (const p of mode.n ? [0, 2 * mode.n] : [0]) {
      const { K, dofs, w } = factor(p);
      radialLoad(p, mode.n, mode.R, slope, elements, dofs, load);
      copy.set(load);
      bandSolve(K, dofs, w, copy);
      for (let i = 0; i < dofs; i++) R += load[i] * copy[i];
    }
    kappa[k] = 24 * (1 - nu * nu) * R * mode.inv / (mode.L * mode.L);
  });
  return kappa;
}

// ---------- rectangular plates ----------

// Fills P[k] = P_k(t) for k < count.
function legendre(t, P, count) {
  P[0] = 1;
  if (count > 1) P[1] = t;
  for (let k = 1; k + 1 < count; k++) P[k + 1] = ((2 * k + 1) * t * P[k] - k * P[k - 1]) / (k + 1);
}

// The 1D Airy functions on ξ ∈ [0,1] with orthonormal second derivatives, so that the stiffness stays well conditioned:
// X'' = √(2k+1)·P_k(t), t = 2ξ − 1. For k ≥ 2 both integrals vanish at both ends (X = X' = 0 there, the stress-free edge).
// The low-order functions complete the set as the end conditions allow:
//   'free': none; 'fixed': whatever X(end), X'(end) a fixed end leaves free; 'straight'/'sliding' (both ends): only the constant,
// and t² joins 'sliding' through airyPairs. `even` keeps even functions only, all that a mode with definite parity can excite.
function airyAxis(end0, end1, count, nodes, even) {
  const low = [];
  const add = (v, d1, d2) => low.push({ v, d1, d2 });
  const fixed0 = end0 === 'fixed', fixed1 = end1 === 'fixed';
  if (fixed0 && fixed1) {
    add(() => 1, () => 0, () => 0);
    add(t => t, () => 2, () => 0);
    add(t => t * t, t => 4 * t, () => 8);
    add(t => t * t * t, t => 6 * t * t, t => 24 * t);
  } else if (fixed0 !== fixed1) { // cubic and quadratic with X = X' = 0 at the free end
    const s = fixed0 ? 1 : -1; // u = 0 at the free end
    add(t => (1 - s * t) ** 2, t => -4 * s * (1 - s * t), () => 8);
    add(t => (1 - s * t) ** 3, t => -6 * s * (1 - s * t) ** 2, t => 24 * (1 - s * t));
  } else if (end0 === 'straight' || end0 === 'sliding') {
    add(() => 1, () => 0, () => 0);
    add(t => t * t, t => 4 * t, () => 8); // pairs only with the constant, see airyPairs
  }
  const orders = [];
  for (let k = 2; orders.length < count; k++) if (!even || k % 2 === 0) orders.push(k);
  const kept = even ? low.filter((_, i) => low[i].v(1) === low[i].v(-1)) : low;
  const size = kept.length + orders.length;
  const nq = nodes.length;
  const V = new Float64Array(size * nq), V1 = new Float64Array(size * nq), V2 = new Float64Array(size * nq);
  const P = new Float64Array(orders.at(-1) + 3);
  for (let q = 0; q < nq; q++) {
    const t = 2 * nodes[q] - 1;
    legendre(t, P, P.length);
    kept.forEach((fn, a) => {
      const scale = 1 / 8; // keeps the low functions on the scale of the others
      V[a * nq + q] = fn.v(t) * scale;
      V1[a * nq + q] = fn.d1(t) * scale;
      V2[a * nq + q] = fn.d2(t) * scale;
    });
    orders.forEach((k, b) => {
      const a = kept.length + b, c = Math.sqrt(2 * k + 1);
      const i1 = (P[k + 1] - P[k - 1]) / (2 * k + 1);
      const i2 = (P[k + 2] - P[k]) / ((2 * k + 1) * (2 * k + 3)) - (P[k] - P[k - 2]) / ((2 * k + 1) * (2 * k - 1));
      V[a * nq + q] = c * i2 / 4;
      V1[a * nq + q] = c * i1 / 2;
      V2[a * nq + q] = c * P[k];
    });
  }
  return { size, low: kept.length, V, V1, V2 };
}

// ∫A_i B_k dξ for all i, k.
function crossIntegrals(A, B, size, nq, weights) {
  const out = new Float64Array(size * size);
  for (let i = 0; i < size; i++) {
    for (let k = 0; k < size; k++) {
      let sum = 0;
      for (let q = 0; q < nq; q++) sum += weights[q] * A[i * nq + q] * B[k * nq + q];
      out[i * size + k] = sum;
    }
  }
  return out;
}

// The tensor-product pairs (i, j) of the Airy basis: the excluded ones carry no stress, or break the shear-free edge.
function airyPairs(X, Y, type) {
  const pairs = [];
  const uniform = type === 'straight' || type === 'sliding';
  for (let i = 0; i < X.size; i++) {
    for (let j = 0; j < Y.size; j++) {
      if (uniform && ((i === 1 && j !== 0) || (j === 1 && i !== 0))) continue;
      if (type === 'straight' && (i === 1 || j === 1)) continue;
      if (uniform && i === 0 && j === 0) continue;
      pairs.push(i, j);
    }
  }
  return Int32Array.from(pairs);
}

// Dense Cholesky in place, lower triangle; a pivot lost to rounding marks a stress-free combination, whose row is dropped.
function denseCholesky(A, n) {
  const dropped = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    const diag = A[i * n + i];
    for (let j = i; j < n; j++) {
      let sum = A[j * n + i];
      for (let k = 0; k < i; k++) sum -= A[j * n + k] * A[i * n + k];
      A[j * n + i] = sum;
    }
    if (A[i * n + i] <= 1e-11 * Math.max(diag, 1e-300)) {
      dropped[i] = 1;
      for (let j = i; j < n; j++) A[j * n + i] = 0;
      A[i * n + i] = 1;
    } else {
      const pivot = Math.sqrt(A[i * n + i]);
      for (let j = i; j < n; j++) A[j * n + i] /= pivot;
    }
  }
  return dropped;
}

// bᵀK⁻¹b from K's Cholesky factor, by forward substitution alone: it is |L⁻¹b|².
function energyOf(L, n, dropped, b) {
  let total = 0;
  for (let i = 0; i < n; i++) {
    if (dropped[i]) {
      b[i] = 0;
      continue;
    }
    let sum = b[i];
    for (let k = 0; k < i; k++) sum -= L[i * n + k] * b[k];
    b[i] = sum / L[i * n + i];
    total += b[i] * b[i];
  }
  return total;
}

// Values and ξ-slopes of a polynomial basis, given only as an evaluator `at`, on the nodes: expand in Legendre, differentiate.
function basisTables(at, nb, nodes) {
  const degree = nb + 24; // exact for the Ritz polynomials (P_k times a weight of degree ≤ 4), and enough for a smooth test shape
  const [gx, gw] = gaussLegendre(degree + 2);
  const f = new Float64Array(degree), f1 = new Float64Array(degree), f2 = new Float64Array(degree);
  const coefficients = new Float64Array(nb * degree);
  gx.forEach((x, q) => {
    rawBasis(0, 0, degree, x, f, f1, f2);
    const u = at(x);
    for (let i = 0; i < nb; i++) {
      for (let m = 0; m < degree; m++) coefficients[i * degree + m] += gw[q] * (2 * m + 1) * u[i] * f[m];
    }
  });
  const nq = nodes.length;
  const V = new Float64Array(nb * nq), D = new Float64Array(nb * nq);
  nodes.forEach((x, q) => {
    rawBasis(0, 0, degree, x, f, f1, f2);
    for (let i = 0; i < nb; i++) {
      for (let m = 0; m < degree; m++) {
        V[i * nq + q] += coefficients[i * degree + m] * f[m];
        D[i * nq + q] += coefficients[i * degree + m] * f1[m];
      }
    }
  });
  return { V, D };
}

const uniformType = ends => (ends.every(e => e === ends[0]) && (ends[0] === 'straight' || ends[0] === 'sliding') ? ends[0] : null);

/**
 * κ of each rectangular mode from rectModes (w = q·φ, ∫∫φ²dξdη = 1), plus the shape's peak |φ|: the plate's peak deflection is
 * A = |q|·peak. inPlane is one type or [ξ=0, ξ=1, η=0, η=1], each 'free' | 'fixed'; for all four edges alike also 'straight'
 * (edges stay straight, zero net force, no shear) or 'sliding' (no normal displacement, no shear: Chu & Herrmann's hinged plate).
 * nAiry is the number of Airy functions per axis beyond the low-order ones.
 * Cost: as rectModes itself, so run it where rectModes runs, or use a Worker. κ holds for every thickness and material.
 */
export function rectKappas(ritz, edges, nu, beta, inPlane, nAiry = 44) {
  const { modes, nbx, nby, atX, atY } = ritz;
  const edgeList = typeof edges === 'string' ? [edges, edges, edges, edges] : edges;
  const inPlaneList = edgeList.map((bc, e) => inPlaneOf(bc, Array.isArray(inPlane) ? inPlane[e] : inPlane));
  const type = uniformType(inPlaneList);
  const evenX = edgeList[0] === edgeList[1] && inPlaneList[0] === inPlaneList[1];
  const evenY = edgeList[2] === edgeList[3] && inPlaneList[2] === inPlaneList[3];
  const countX = evenX ? nAiry >> 1 : nAiry, countY = evenY ? nAiry >> 1 : nAiry;

  const nq = (Math.max(nbx, nby) + (nAiry >> 1) + 10) | 1; // odd, so that the centre, where low modes peak, is a node
  const [nodes, weights] = gaussLegendre(nq);
  const X = airyAxis(inPlaneList[0], inPlaneList[1], countX, nodes, evenX);
  const Y = airyAxis(inPlaneList[2], inPlaneList[3], countY, nodes, evenY);

  // B(Ψ,ψ) in the tensor basis: the plate stiffness of rectModes with ν → −ν
  const b2 = beta * beta, b4 = b2 * b2;
  const integrate = (T, A, B) => crossIntegrals(A, B, T.size, nq, weights);
  const [ax, cx, dx, fx] = [integrate(X, X.V, X.V), integrate(X, X.V1, X.V1), integrate(X, X.V2, X.V2), integrate(X, X.V2, X.V)];
  const [ay, cy, dy, fy] = [integrate(Y, Y.V, Y.V), integrate(Y, Y.V1, Y.V1), integrate(Y, Y.V2, Y.V2), integrate(Y, Y.V2, Y.V)];
  const pairs = airyPairs(X, Y, type);
  const n = pairs.length / 2;
  const K = new Float64Array(n * n);
  const sx = X.size, sy = Y.size;
  for (let a = 0; a < n; a++) {
    const i = pairs[2 * a], j = pairs[2 * a + 1];
    for (let b = 0; b <= a; b++) {
      const k = pairs[2 * b], l = pairs[2 * b + 1];
      K[a * n + b] = dx[i * sx + k] * ay[j * sy + l] + ax[i * sx + k] * dy[j * sy + l] / b4
        + (-nu * (fx[i * sx + k] * fy[l * sy + j] + fx[k * sx + i] * fy[j * sy + l]) + 2 * (1 + nu) * cx[i * sx + k] * cy[j * sy + l]) / b2;
    }
  }
  const dropped = denseCholesky(K, n);

  const px = basisTables(atX, nbx, nodes), py = basisTables(atY, nby, nodes);
  const kappa = new Float64Array(modes.length), peak = new Float64Array(modes.length);

  // scratch for the per-mode contractions
  const tA = new Float64Array(nbx * nq), tB = new Float64Array(nbx * nq);
  const sxField = new Float64Array(nq * nq), syField = new Float64Array(nq * nq), value = new Float64Array(nq * nq);
  const fieldXX = new Float64Array(nq * nq), fieldYY = new Float64Array(nq * nq), fieldXY = new Float64Array(nq * nq);
  const load = new Float64Array(sx * sy), rhs = new Float64Array(n);
  const partial = new Float64Array(sx * nq);

  // load[i,j] += ∫∫ A_i(ξ) F(ξ,η) B_j(η), through the intermediate partial[i,η]
  const addLoad = (A, F, B, scale) => {
    partial.fill(0);
    for (let i = 0; i < sx; i++) {
      for (let q = 0; q < nq; q++) {
        const w = weights[q] * A[i * nq + q];
        for (let r = 0; r < nq; r++) partial[i * nq + r] += w * F[q * nq + r];
      }
    }
    for (let i = 0; i < sx; i++) {
      for (let j = 0; j < sy; j++) {
        let sum = 0;
        for (let r = 0; r < nq; r++) sum += weights[r] * partial[i * nq + r] * B[j * nq + r];
        load[i * sy + j] += scale * sum;
      }
    }
  };

  modes.forEach((mode, m) => {
    if (mode.L < 1e-3) return;
    const v = mode.vec;
    // along η first: tA[i,η] = Σ_j v_ij Y_j(η), tB[i,η] = Σ_j v_ij Y_j'(η)
    tA.fill(0);
    tB.fill(0);
    for (let i = 0; i < nbx; i++) {
      for (let j = 0; j < nby; j++) {
        const c = v[i * nby + j];
        if (c === 0) continue;
        for (let r = 0; r < nq; r++) {
          tA[i * nq + r] += c * py.V[j * nq + r];
          tB[i * nq + r] += c * py.D[j * nq + r];
        }
      }
    }
    sxField.fill(0);
    syField.fill(0);
    value.fill(0);
    for (let i = 0; i < nbx; i++) {
      for (let q = 0; q < nq; q++) {
        const d = px.D[i * nq + q], u = px.V[i * nq + q];
        for (let r = 0; r < nq; r++) {
          sxField[q * nq + r] += d * tA[i * nq + r];
          syField[q * nq + r] += u * tB[i * nq + r];
          value[q * nq + r] += u * tA[i * nq + r];
        }
      }
    }
    let top = 0;
    for (let t = 0; t < nq * nq; t++) {
      top = Math.max(top, Math.abs(value[t]));
      fieldXX[t] = sxField[t] * sxField[t];
      fieldYY[t] = syField[t] * syField[t];
      fieldXY[t] = sxField[t] * syField[t];
    }
    peak[m] = top;
    // R(ψ) = ∫∫ ½ψ,ηη φ,ξ² + ½ψ,ξξ φ,η² − ψ,ξη φ,ξ φ,η
    load.fill(0);
    addLoad(X.V, fieldXX, Y.V2, 0.5);
    addLoad(X.V2, fieldYY, Y.V, 0.5);
    addLoad(X.V1, fieldXY, Y.V1, -1);
    for (let a = 0; a < n; a++) rhs[a] = load[pairs[2 * a] * sy + pairs[2 * a + 1]] / b2;
    // stationary: B(Ψ,·) = R(·)/β², so R(Ψ) = β²·rhsᵀK⁻¹rhs
    const R = b2 * energyOf(K, n, dropped, rhs);
    kappa[m] = 24 * (1 - nu * nu) * R / (b2 * mode.L * mode.L);
  });
  return { kappa, peak };
}

// ---------- the plate bolted at the driver ----------

/**
 * Modes of the plate pinned at the driver, as combinations of its free modes. A rigid bolt resonates where G(x₀,x₀) = 0, between
 * consecutive free modes that move the bolt point, and these resonant shapes mix several free modes (largest coefficients about
 * 0.7, 0.5, 0.3). Their stretching energy is dominated by cross terms the per-mode κ cannot see: Σκ_k c_k⁴ falls short of the true
 * value by a factor of 2 to 100. So a rigid mount needs κ of these shapes, from circleKappas (driver at the centre only, as
 * off-centre shapes are not one angular harmonic) or rectKappas.
 * L: the free modes' Λ. d0: each mode's shape at the driver for unit modal norm. Returns [{ L, c }], c the free-mode coefficients,
 * Σc² = 1. Modes that do not move the driver stay as they are and are not returned.
 */
export function pinnedModes(L, d0, count = Infinity) {
  const active = L.map((_, k) => k).filter(k => Math.abs(d0[k]) > 1e-9).sort((p, q) => L[p] - L[q]);
  const receptance = x => active.reduce((sum, k) => sum + d0[k] * d0[k] / (L[k] * L[k] - x), 0);
  const modes = [];
  for (let i = 0; i + 1 < active.length && modes.length < count; i++) {
    let lo = L[active[i]] ** 2, hi = L[active[i + 1]] ** 2;
    if (hi - lo < 1e-9 * hi) continue; // a degenerate pair leaves no root between its members
    lo += 1e-10 * (hi - lo);
    hi -= 1e-10 * (hi - lo);
    for (let iter = 0; iter < 80; iter++) {
      const mid = (lo + hi) / 2;
      if (receptance(mid) > 0) hi = mid;
      else lo = mid;
    }
    const x = (lo + hi) / 2;
    const c = new Float64Array(L.length);
    let norm = 0;
    for (const k of active) {
      c[k] = d0[k] / (L[k] * L[k] - x);
      norm += c[k] * c[k];
    }
    norm = Math.sqrt(norm);
    for (const k of active) c[k] /= norm;
    modes.push({ L: Math.sqrt(x), c });
  }
  return modes;
}

// ---------- steady state, per drive frequency ----------

/**
 * Harmonic balance of q̈ + 2ζω_k q̇ + ω_k²q + γq³ = f cos ωt: z = ¾γ|q|² is a root of
 *   z[(Ω + z)² + b²] = S,   Ω = ω_k² − ω² (detuning), b = 2ζω_kω (damping), S = ¾γf² (drive),
 * and q = f/(Ω + z + ib). The cubic is solved in closed form, so there is nothing to diverge. Where it has three positive roots
 * the middle one is unstable and separates the basins: the state stays on the branch `previous` (last frame's z) is on. When its
 * branch ends at a fold the root simply moves to the other branch within that call, which is the jump of a Duffing sweep; a real
 * plate takes about 1/(ζω_k) seconds to get there. Allocation-free, about 100 ns per call.
 */
export function duffingShift(detuning, damping, drive, previous = 0) {
  if (!(drive > 0)) return 0;
  const O = detuning, b2 = damping * damping;
  // z = t − 2Ω/3 gives t³ + Pt + Q = 0
  const P = b2 - O * O / 3;
  const Q = -(2 / 27) * O * O * O - (2 / 3) * O * b2 - drive;
  const disc = Q * Q / 4 + P * P * P / 27;
  let z;
  if (disc >= 0) {
    const u = Math.cbrt(-Q / 2 - Math.sign(Q || 1) * Math.sqrt(disc));
    z = u - P / (3 * u) - 2 * O / 3;
  } else {
    const radius = 2 * Math.sqrt(-P / 3);
    const angle = Math.acos(Math.max(-1, Math.min(1, 3 * Q / (P * radius)))) / 3;
    const upper = radius * Math.cos(angle) - 2 * O / 3;
    const middle = radius * Math.cos(angle + 4 * Math.PI / 3) - 2 * O / 3;
    z = previous > middle ? upper : radius * Math.cos(angle + 2 * Math.PI / 3) - 2 * O / 3;
  }
  // one Newton step restores the digits lost to z ≪ |Ω|
  const a = O + z;
  const slope = a * a + b2 + 2 * z * a;
  if (slope > 0) z -= (z * (a * a + b2) - drive) / slope;
  return Math.max(z, 0);
}

// ---------- when to warn ----------

/**
 * s = shift/half-width = (⅜κ(q/h)²)/ζ = z/(2ζω_k²), the frequency shift in half-power half-widths, sets the regime:
 *  s < 0.1  'negligible': the peak moves by under a tenth of its half-width, and the amplitude changes by under 0.5 % at the peak.
 *  s < 1.5  'corrected':  harmonic balance within 1 % at ζ ≤ 0.005 and 6 % at ζ ≤ 0.05 of a multi-mode integration. Folds, and so
 *                          jumps in a sweep, begin at s = 2/√3 ≈ 1.15.
 *  s < 3    'caution':    4 runs in 16 lost the branch at s = 3 in the integration.
 *  s ≥ 3    'broken':     other modes take over: 14 of 16 runs at s = 6 ended below half the one-mode amplitude.
 * Measured on a simply supported plate with 19 modes (one mode repeated with 46), ζ = 0.002–0.05. s = 3 is also where a perfect
 * circular plate's unforced companion mode goes unstable. In peak deflection A: A/h = √(8ζs/(3κ_A)), κ_A = κ/peak².
 */
export function regime(s) {
  return s < 0.1 ? 'negligible' : s < 1.5 ? 'corrected' : s < 3 ? 'caution' : 'broken';
}

/** The largest A/h in a regime: the s of regime() for a mode with κ_A = κ/peak² and damping ratio ζ. */
export const amplitudeLimit = (s, kappaA, zeta) => Math.sqrt(8 * zeta * s / (3 * kappaA));
