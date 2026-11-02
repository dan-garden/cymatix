// Fluid loading of a thin plate surrounded by fluid on both sides, without a baffle.
// The flow is antisymmetric about the plate, so the potential jump μ across it obeys the hypersingular equation
//   S μ = w on the plate,  with Fourier symbol S(k) = ½√(k² − k_a²),
// and the fluid adds the kinetic energy ½ρ_f ω² ∫ w S⁻¹w dS. Then Γ_k = ∫ w S⁻¹w dS / (L ∫ w² dS) and
//   β_k = (ρ_f L / ρh) Γ_k,  ω_fluid = ω_vacuum / √(1 + β_k).
// Γ is Kwak's NAVMI factor doubled, because φ = 0 outside the plate in his free-surface problem is this same flow.
// The solution needs μ ~ √(distance) at the edge, which the bases below build in; that is the flow around the edge.
// For large k, Γ → 2/√(k²−k_a²) (infinite plate, here Γ̃ with L = 1), the limit in the compressible sections.
//
// To use it on a configured Plate: added = circleAddedMass(modes) or rectAddedMass(ritz, Ly/Lx), cached by shape; then
// plate.modes = loadModes(vacuumModes, added, fluid.rho * L / (rho * h)) followed by setDrive, L being the radius or Lx.

import { basis1D, bessel, gaussLegendre, NR, ORDER, symEig } from './physics.js';
import { radialShapeAt } from './mindlin.js';

const TAU = 2 * Math.PI;

// ρ kg/m³, c m/s, μ Pa·s: Escaler & De La Torre 2018, Table 4
export const VACUUM = { rho: 0, c: Infinity, mu: 0 };
export const AIR = { rho: 1.2041, c: 343.24, mu: 1.821e-5 };
export const WATER = { rho: 998, c: 1482, mu: 1.002e-3 };

// ---------- circular plates ----------

// ψ_s(ρ) = s!/(√2 Γ(s+3/2)) · ρⁿ √(1−ρ²) P_s^(n,½)(1−2ρ²) has Hankel transform J_ν(k)/k^{3/2} with ν = n + 2s + 3/2
// (Weber–Schafheitlin, DLMF 10.22.56), so ⟨ψ_s, Sψ_t⟩ = ½∫J_ν_s J_ν_t dk/k = δ_st/(4ν_s) (DLMF 10.22.57): S is diagonal on them.
// Fills table[s][q] with ψ_s at ρ² = 1 − t_q², without the factor √(1−ρ²) = t_q.
function jacobiTable(n, S, t) {
  const table = Array.from({ length: S }, () => new Float64Array(t.length));
  const beta = 0.5;
  t.forEach((tq, q) => {
    const y = 2 * tq * tq - 1; // 1 − 2ρ²
    const rhoN = (1 - tq * tq) ** (n / 2);
    let prev = 0, P = 1, front = Math.SQRT2 / Math.sqrt(Math.PI); // 1/(√2 Γ(3/2))
    for (let s = 0; s < S; s++) {
      table[s][q] = front * rhoN * P;
      const a = 2 * s + n + beta; // three-term recurrence for P_s^(n,β)
      const next = s === 0
        ? ((n + beta + 2) * y + n - beta) / 2 * P
        : ((a + 1) * ((a + 2) * a * y + n * n - beta * beta) * P - 2 * (s + n) * (s + beta) * (a + 2) * prev)
          / (2 * (s + 1) * (s + 1 + n + beta) * a);
      prev = P;
      P = next;
      front *= (s + 1) / (s + 1.5);
    }
  });
  return table;
}

// √(1−ρ²) = t turns the edge singularity into a polynomial: every integrand below is analytic in t.
const circleQuadrature = (lam, n) => gaussLegendre(Math.ceil(lam + n) + 48);

const circleOrder = (n, top) => Math.max(10, Math.ceil((top - n) / 2) + 10); // beyond ν ≈ λ the projections fall exponentially

// Projections c_s = ∫ψ_s R ρdρ of one mode, and ∫R²ρdρ.
function circleProjection(mode, S, table, tq, wq) {
  const c = new Float64Array(S);
  let norm = 0;
  for (let q = 0; q < tq.length; q++) {
    const t = tq[q];
    const R = radialShapeAt(mode, Math.sqrt(1 - t * t));
    norm += wq[q] * R * R * t; // ρdρ = −t dt
    const weight = wq[q] * R * t * t;
    for (let s = 0; s < S; s++) c[s] += weight * table[s][q];
  }
  return [c, norm];
}

/**
 * Added-mass matrices of circular modes, L = radius. One block per n, because the fluid couples only equal n. With shapes
 * normalised to ∫φ²dS = 1, C_ij = ∫φ_i S⁻¹φ_j dS / L, so Γ_k = C_kk. About 30 ms for 300 modes.
 */
export function circleAddedMass(modes) {
  const groups = new Map();
  modes.forEach((mode, k) => {
    if (!groups.has(mode.n)) groups.set(mode.n, []);
    groups.get(mode.n).push(k);
  });
  const blocks = [...groups].map(([n, members]) => {
    const top = Math.max(...members.map(k => Math.sqrt(modes[k].L)));
    const S = circleOrder(n, top);
    const [tq, wq] = circleQuadrature(top, n);
    const table = jacobiTable(n, S, tq);
    const projections = members.map(k => circleProjection(modes[k], S, table, tq, wq));
    const N = members.length;
    const C = new Float64Array(N * N);
    projections.forEach(([ci, ni], i) => projections.slice(0, i + 1).forEach(([cj, nj], j) => {
      let sum = 0;
      for (let s = 0; s < S; s++) sum += 4 * (n + 2 * s + 1.5) * ci[s] * cj[s];
      C[i * N + j] = C[j * N + i] = sum / Math.sqrt(ni * nj);
    }));
    return { members, C, norms: projections.map(([, norm]) => norm) };
  });
  return { blocks, circle: true };
}

// ---------- compressible loading of circular plates ----------

/** Spherical Bessel j_0(x)…j_top(x) into out: downward recurrence, normalised on the larger of j_0, j_1. */
function sphericalRow(x, top, out) {
  const start = 2 * Math.ceil((Math.max(top, x) + 20 + Math.sqrt(40 * Math.max(top, x))) / 2) + 1;
  let above = 0, here = 1e-300;
  for (let m = start; m > 0; m--) {
    const below = (2 * m + 1) / x * here - above;
    above = here;
    here = below;
    if (m - 1 <= top) out[m - 1] = here;
    if (Math.abs(here) > 1e250) {
      above *= 1e-250;
      here *= 1e-250;
      for (let i = m - 1; i <= top; i++) out[i] *= 1e-250;
    }
  }
  const j0 = Math.sin(x) / x, j1 = (j0 - Math.cos(x)) / x;
  const scale = Math.abs(j0) > Math.abs(j1) ? j0 / out[0] : j1 / out[1];
  for (let i = 0; i <= top; i++) out[i] *= scale;
}

// Nodes of ∫dk on either side of k = κ, each with the substitution that removes the square root there: k = κ sinθ below
// and s = √(k²−κ²) above. Returns [k, weight, symbol part] triples, the part being γ = √(κ²−k²) or s.
function ringNodes(kappa, top) {
  const [gx, gw] = gaussLegendre(8);
  const below = [], above = [];
  const nb = Math.max(1, Math.ceil(kappa / 1.5));
  for (let p = 0; p < nb && kappa > 0; p++) {
    for (let q = 0; q < 8; q++) {
      const theta = (p + gx[q]) / nb * Math.PI / 2;
      const gamma = kappa * Math.cos(theta);
      below.push([kappa * Math.sin(theta), gw[q] / nb * Math.PI / 2 * gamma, gamma]);
    }
  }
  const sTop = Math.sqrt(top * top - kappa * kappa);
  const na = Math.ceil(sTop / 1.5);
  for (let p = 0; p < na; p++) {
    for (let q = 0; q < 8; q++) {
      const s = (p + gx[q]) / na * sTop;
      const k = Math.hypot(kappa, s);
      above.push([k, gw[q] / na * sTop * s / k, s]);
    }
  }
  return [below, above];
}

// Solves the complex symmetric system (Are + iAim)x = b by elimination; the operator is never near singular, because a
// screen has no interior resonances. Returns [xr, xi].
function solveComplex(Are, Aim, bre, bim, n) {
  const ar = Float64Array.from(Are), ai = Float64Array.from(Aim), xr = Float64Array.from(bre), xi = Float64Array.from(bim);
  for (let c = 0; c < n; c++) {
    const pr = ar[c * n + c], pi = ai[c * n + c], d = pr * pr + pi * pi;
    for (let i = c + 1; i < n; i++) {
      const fr = (ar[i * n + c] * pr + ai[i * n + c] * pi) / d, fi = (ai[i * n + c] * pr - ar[i * n + c] * pi) / d;
      for (let j = c; j < n; j++) {
        const rr = ar[c * n + j], ri = ai[c * n + j];
        ar[i * n + j] -= fr * rr - fi * ri;
        ai[i * n + j] -= fr * ri + fi * rr;
      }
      const rr = xr[c], ri = xi[c];
      xr[i] -= fr * rr - fi * ri;
      xi[i] -= fr * ri + fi * rr;
    }
  }
  for (let i = n - 1; i >= 0; i--) {
    let sr = xr[i], si = xi[i];
    for (let j = i + 1; j < n; j++) {
      sr -= ar[i * n + j] * xr[j] - ai[i * n + j] * xi[j];
      si -= ar[i * n + j] * xi[j] + ai[i * n + j] * xr[j];
    }
    const pr = ar[i * n + i], pi = ai[i * n + i], d = pr * pr + pi * pi;
    xr[i] = (sr * pr + si * pi) / d;
    xi[i] = (si * pr - sr * pi) / d;
  }
  return [xr, xi];
}

/**
 * Complex loading Γ̃ = [Γ_R, Γ_I] of one circular mode, meaning Γ_R − iΓ_I under e^{iωt}, at κ = k_a·a = ω·a/c. Solves
 * the full compressible problem on the ψ_s basis, whose matrix is
 *   S_st(κ) = δ_st/(4ν_s) + (1/π)∫[√(k²−κ²) − k] j_ms(k) j_mt(k) dk/k,  m_s = n + 2s + 1,  √(k²−κ²) = iγ below κ.
 * About 1 ms per mode. Limits: Γ_I → 16κ³/(27π²) for a rigid disc (compact dipole), Γ̃ → 2/√(k²−κ²) for large k.
 */
export function circleRadiation(mode, kappa) {
  const n = mode.n, lam = Math.sqrt(mode.L);
  const S = circleOrder(n, Math.max(lam, kappa));
  const [tq, wq] = circleQuadrature(lam, n);
  const [c, norm] = circleProjection(mode, S, jacobiTable(n, S, tq), tq, wq);
  const top = n + 2 * S + 1;
  const row = new Float64Array(top + 1);
  const re = new Float64Array(S * S), im = new Float64Array(S * S);
  const [below, above] = ringNodes(kappa, Math.max(lam, kappa) + top + 80);
  const add = (k, weight, symbolRe, symbolIm) => {
    sphericalRow(k, top, row);
    for (let s = 0; s < S; s++) {
      const a = row[n + 2 * s + 1] * weight / (Math.PI * k);
      for (let t = s; t < S; t++) {
        const value = a * row[n + 2 * t + 1];
        re[s * S + t] += value * symbolRe;
        im[s * S + t] += value * symbolIm;
      }
    }
  };
  for (const [k, weight, gamma] of below) add(k, weight, -k, gamma);
  for (const [k, weight, root] of above) add(k, weight, root - k, 0);
  for (let s = 0; s < S; s++) {
    re[s * S + s] += 1 / (4 * (n + 2 * s + 1.5));
    for (let t = 0; t < s; t++) {
      re[s * S + t] = re[t * S + s];
      im[s * S + t] = im[t * S + s];
    }
  }
  const [xr, xi] = solveComplex(re, im, c, new Float64Array(S), S);
  let gr = 0, gi = 0;
  for (let s = 0; s < S; s++) {
    gr += c[s] * xr[s];
    gi += c[s] * xi[s];
  }
  return [gr / norm, -gi / norm];
}

// ---------- rectangular plates ----------

/** J_0(u)…J_top(u) into out, by Miller's downward recurrence, normalised with J_0 + 2ΣJ_2k = 1. */
function besselRow(u, top, out) {
  const start = 2 * Math.ceil((Math.max(top, u) + 20 + Math.sqrt(40 * Math.max(top, u))) / 2);
  let above = 0, here = 1e-300, even = 0;
  out.fill(0);
  for (let k = start; k > 0; k--) {
    const below = 2 * k / u * here - above;
    above = here;
    here = below;
    if (k - 1 <= top) out[k - 1] = here;
    if (k > 1 && (k - 1) % 2 === 0) even += here;
    if (Math.abs(here) > 1e250) {
      above *= 1e-250;
      here *= 1e-250;
      even *= 1e-250;
      for (let i = k - 1; i <= top; i++) out[i] *= 1e-250;
    }
  }
  const scale = 1 / (here + 2 * even);
  for (let i = 0; i <= top; i++) out[i] *= scale;
}

// Composite Gauss–Legendre on [0, top]: graded towards the origin, where the symbol has a conical point, then panels
// of width 2, which give nine points per period π of J_iJ_k.
function spectralNodes(top) {
  const [gx, gw] = gaussLegendre(6);
  const edges = [0, 1 / 64, 1 / 32, 1 / 16, 1 / 8, 1 / 4, 1 / 2, 1];
  while (edges.at(-1) < top) edges.push(edges.at(-1) + 2);
  const count = (edges.length - 1) * 6;
  const u = new Float64Array(count), w = new Float64Array(count);
  for (let p = 0; p + 1 < edges.length; p++) {
    const width = edges[p + 1] - edges[p];
    for (let q = 0; q < 6; q++) {
      u[p * 6 + q] = edges[p] + width * gx[q];
      w[p * 6 + q] = width * gw[q];
    }
  }
  return [u, w];
}

// With p_i(t) = √(1−t²)U_i(t), p̂_i(u) = π(i+1)iⁱJ_{i+1}(u)/u. Rows of g_i(u) = (−1)^⌊i/2⌋(i+1)J_{i+1}(u)/u on the nodes,
// so that within one parity p̂_i p̂_k* = π² g_i g_k.
function spectralRows(P, u) {
  const row = new Float64Array(P + 1);
  const g = Array.from({ length: P }, () => new Float64Array(u.length));
  u.forEach((value, m) => {
    besselRow(value, P, row);
    for (let i = 0; i < P; i++) g[i][m] = (i >> 1) % 2 ? -(i + 1) * row[i + 1] / value : (i + 1) * row[i + 1] / value;
  });
  return g;
}

function choleskySolver(A, n) { // A is overwritten by its lower factor
  for (let j = 0; j < n; j++) {
    for (let i = j; i < n; i++) {
      let sum = A[i * n + j];
      for (let k = 0; k < j; k++) sum -= A[i * n + k] * A[j * n + k];
      A[i * n + j] = i === j ? Math.sqrt(sum) : sum / A[j * n + j];
    }
  }
  return b => { // b ← A⁻¹b
    for (let i = 0; i < n; i++) {
      let sum = b[i];
      for (let k = 0; k < i; k++) sum -= A[i * n + k] * b[k];
      b[i] = sum / A[i * n + i];
    }
    for (let i = n - 1; i >= 0; i--) {
      let sum = b[i];
      for (let k = i + 1; k < n; k++) sum -= A[k * n + i] * b[k];
      b[i] = sum / A[i * n + i];
    }
    return b;
  };
}

/**
 * The hypersingular operator on a 1 × beta rectangle in the Galerkin basis p_i(t)p_j(s): one factorised block per parity
 * of (i, j), which the geometry decouples whatever the edge conditions. With a = 1/hx, b = 1/hy, its symbol splits as
 *   √(a²u² + b²v²) = au + bv + E(u,v).
 * The au and bv terms are exact in one axis (∫J_iJ_k du/u = δ_ik/2i, DLMF 10.22.57). E decays fast enough to truncate at
 * U, and its tails (E → −bv as u → ∞, and likewise in v) turn the one-axis integrals left over into the truncated ones.
 */
function rectOperator(beta, Px, Py) {
  const hx = 0.5, hy = beta / 2, a = 1 / hx, b = 1 / hy;
  const [u, wu] = spectralNodes(3 * Px + 40);
  const [v, wv] = spectralNodes(3 * Py + 40);
  const gx = spectralRows(Px, u), gy = spectralRows(Py, v);
  const Nu = u.length, Nv = v.length, front = Math.PI ** 2 * hx * hy / 2;

  const pairsOf = P => { // every i ≤ k of equal parity, with its index in a P×P lookup
    const index = new Int32Array(P * P), list = [];
    for (let i = 0; i < P; i++) {
      for (let k = i; k < P; k += 2) {
        index[i * P + k] = index[k * P + i] = list.length;
        list.push([i, k]);
      }
    }
    return [list, index];
  };
  const tabulate = (pairs, g, w) => pairs.map(([i, k]) => {
    const G = Float64Array.from(w, (weight, m) => weight * g[i][m] * g[k][m]);
    return { G, total: G.reduce((sum, value) => sum + value, 0), D: i === k ? (i + 1) / 2 : 0 };
  });
  const [pairsX, indexX] = pairsOf(Px), [pairsY, indexY] = pairsOf(Py);
  const X = tabulate(pairsX, gx, wu), Y = tabulate(pairsY, gy, wv);

  const E = new Float64Array(Nu * Nv);
  for (let m = 0; m < Nu; m++) {
    for (let n = 0; n < Nv; n++) E[m * Nv + n] = Math.hypot(a * u[m], b * v[n]) - a * u[m] - b * v[n];
  }

  // A_ij,kl = (π²hxhy/2)·[a·D_ik·I_jl + b·I_ik·D_jl + ∫∫G_ik G_jl E], for every pair of pairs
  const scratch = new Float64Array(Nu);
  const Q = new Float64Array(X.length * Y.length);
  Y.forEach((y, yIndex) => {
    for (let m = 0; m < Nu; m++) {
      let sum = 0;
      const offset = m * Nv;
      for (let n = 0; n < Nv; n++) sum += y.G[n] * E[offset + n];
      scratch[m] = sum;
    }
    X.forEach((x, xIndex) => {
      let sum = 0;
      for (let m = 0; m < Nu; m++) sum += x.G[m] * scratch[m];
      Q[xIndex * Y.length + yIndex] = front * (a * x.D * y.total + b * x.total * y.D + sum);
    });
  });

  const blocks = [];
  for (let px = 0; px < 2; px++) {
    for (let py = 0; py < 2; py++) {
      const pairs = [];
      for (let i = px; i < Px; i += 2) for (let j = py; j < Py; j += 2) pairs.push([i, j]);
      const n = pairs.length;
      const A = new Float64Array(n * n);
      pairs.forEach(([i, j], r) => pairs.forEach(([k, l], c) => {
        A[r * n + c] = Q[indexX[i * Px + k] * Y.length + indexY[j * Py + l]];
      }));
      blocks.push({ px, py, pairs, solve: choleskySolver(A, n) });
    }
  }
  return { blocks, Px, Py, hx, hy, grids: new Map() };
}

// ∫p_i(t) f(ξ(t)) dt for the Ritz functions f, exact by Gauss–Chebyshev of the second kind.
function ritzOverlap(at, count, P) {
  const N = count + P + 8;
  const overlap = Array.from({ length: count }, () => new Float64Array(P));
  for (let q = 1; q <= N; q++) {
    const theta = q * Math.PI / (N + 1);
    const t = Math.cos(theta), s = Math.sin(theta);
    const weight = Math.PI / (N + 1) * s * s;
    const values = at((t + 1) / 2);
    let previous = 0, U = 1;
    for (let i = 0; i < P; i++) {
      for (let a = 0; a < count; a++) overlap[a][i] += weight * values[a] * U;
      const next = 2 * t * U - previous;
      previous = U;
      U = next;
    }
  }
  return overlap;
}

/**
 * Added-mass matrices of rectangular Ritz modes, L = Lx, from what rectModes returns and beta = Ly/Lx. The fluid couples
 * modes of equal parity about the centre lines. Px, Py are fluid basis sizes per axis; the Ritz sizes are the default, which
 * converge Γ to 0.1 % or better over the kept modes. About 0.5 s for a 1 × 1 plate (324 modes): use a Worker and cache by aspect.
 */
export function rectAddedMass(ritz, beta, { Px = ritz.nbx, Py = ritz.nby } = {}) {
  const { modes, nbx, nby, atX, atY } = ritz;
  const operator = rectOperator(beta, Px, Py);
  const TX = ritzOverlap(atX, nbx, Px), TY = ritzOverlap(atY, nby, Py);
  const area = operator.hx * operator.hy;
  const half = new Float64Array(Px * nby);

  // ⟨p_i p_j, φ⟩ of each mode
  const moments = modes.map(({ vec }) => {
    half.fill(0);
    for (let a = 0; a < nbx; a++) {
      for (let j = 0; j < nby; j++) {
        const value = vec[a * nby + j];
        if (value) for (let i = 0; i < Px; i++) half[i * nby + j] += value * TX[a][i];
      }
    }
    const b = new Float64Array(Px * Py);
    for (let i = 0; i < Px; i++) {
      for (let j = 0; j < nby; j++) {
        const value = half[i * nby + j];
        if (value) for (let l = 0; l < Py; l++) b[i * Py + l] += value * area * TY[j][l];
      }
    }
    return b;
  });

  // A mode lives in the parity blocks it touches; with unequal edges (px or py = −1) it spans two or four.
  const groups = new Map();
  modes.forEach((mode, k) => {
    const key = `${mode.px}|${mode.py}`;
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push(k);
  });
  const blocks = [...groups.values()].map(members => {
    const N = members.length;
    const C = new Float64Array(N * N);
    for (const { pairs, solve } of operator.blocks) {
      const rhs = members.map(k => Float64Array.from(pairs, ([i, j]) => moments[k][i * Py + j]));
      if (!rhs.some(b => b.some(value => Math.abs(value) > 1e-13))) continue;
      const solved = rhs.map(b => solve(b.slice()));
      for (let r = 0; r < N; r++) {
        for (let c = 0; c <= r; c++) {
          let sum = 0;
          for (let q = 0; q < pairs.length; q++) sum += rhs[r][q] * solved[c][q];
          C[r * N + c] += sum / beta; // ∫φ²dS = beta with Lx = 1 and ∫∫φ²dξdη = 1
        }
      }
    }
    for (let r = 0; r < N; r++) for (let c = 0; c < r; c++) C[c * N + r] = C[r * N + c];
    return { members, C };
  });
  return { blocks, circle: false, operator, moments, beta, waves: modes.map(mode => Math.sqrt(mode.L)) };
}

/** Γ_k = C_kk of every mode. */
export function gammas(added, count) {
  const out = new Float64Array(count);
  for (const { members, C } of added.blocks) members.forEach((k, i) => { out[k] = C[i * members.length + i]; });
  return out;
}

/**
 * Modes loaded by a fluid: new mode objects, keeping the vacuum ones untouched, with the same fields as Plate expects.
 * ratio = ρ_f L/(ρh), L as in the added-mass call; the vacuum frequencies are mode.w. Each block solves
 *   diag(ω²)x = ω_f²(I + ratio·C)x,
 * which changes the shapes, and so nodal lines, too. coupled: false keeps the shapes: ω_f = ω/√(1 + β_k), mass × (1 + β_k).
 * ratio = 0 returns the same modes. A loaded circular shape is no longer J_n + C·I_n, so nodalRadii does not apply to it.
 */
export function loadModes(modes, added, ratio, { coupled = true } = {}) {
  if (!ratio) return modes;
  const loaded = modes.slice();
  for (const block of added.blocks) {
    const { members, C } = block;
    const N = members.length;
    if (!coupled) {
      members.forEach((k, i) => {
        const inertia = 1 + ratio * C[i * N + i];
        const w = modes[k].w / Math.sqrt(inertia);
        loaded[k] = { ...modes[k], w, f: w / TAU, inv: modes[k].inv / inertia };
      });
      continue;
    }
    const inverseFactor = invertCholesky(C, ratio, N);
    // L⁻¹·diag(ω²)·L⁻ᵀ, whose eigenvectors y give x = L⁻ᵀy
    const A = new Float64Array(N * N);
    for (let i = 0; i < N; i++) {
      for (let j = 0; j <= i; j++) {
        let sum = 0;
        for (let k = 0; k <= j; k++) sum += inverseFactor[i * N + k] * modes[members[k]].w ** 2 * inverseFactor[j * N + k];
        A[i * N + j] = A[j * N + i] = sum;
      }
    }
    const { values, vectors } = symEig(A, N);
    const order = Array.from(values.keys()).sort((p, q) => values[p] - values[q]);
    order.forEach((e, index) => {
      const x = new Float64Array(N);
      for (let i = 0; i < N; i++) {
        for (let k = i; k < N; k++) x[i] += inverseFactor[k * N + i] * vectors[k * N + e];
      }
      const k = members[index];
      const w = Math.sqrt(Math.max(values[e], 0));
      loaded[k] = { ...modes[k], w, f: w / TAU, ...(added.circle ? circleShape(modes, block, x) : rectShape(modes, members, x)) };
    });
  }
  return loaded;
}

function invertCholesky(C, ratio, N) { // inverse of the lower factor of I + ratio·C
  const L = new Float64Array(N * N), inverse = new Float64Array(N * N);
  for (let j = 0; j < N; j++) {
    for (let i = j; i < N; i++) {
      let sum = ratio * C[i * N + j] + (i === j ? 1 : 0);
      for (let k = 0; k < j; k++) sum -= L[i * N + k] * L[j * N + k];
      L[i * N + j] = i === j ? Math.sqrt(sum) : sum / L[j * N + j];
    }
  }
  for (let c = 0; c < N; c++) {
    for (let i = c; i < N; i++) {
      let sum = i === c ? 1 : 0;
      for (let k = c; k < i; k++) sum -= L[i * N + k] * inverse[k * N + c];
      inverse[i * N + c] = sum / L[i * N + i];
    }
  }
  return inverse;
}

function circleShape(modes, { members, norms }, x) { // xᵀ(I + ratio·C)x = 1, so ∫R²ρdρ counts as 1 in inv
  const R = new Float32Array(NR);
  const sum = new Float64Array(NR);
  members.forEach((k, i) => {
    const scale = x[i] / Math.sqrt(norms[i]), source = modes[k].R;
    for (let j = 0; j < NR; j++) sum[j] += scale * source[j];
  });
  let Rmax = 0;
  for (let j = 0; j < NR; j++) {
    R[j] = sum[j];
    Rmax = Math.max(Rmax, Math.abs(sum[j]));
  }
  const mix = members.map((k, i) => ({ part: modes[k], scale: x[i] / Math.sqrt(norms[i]) }));
  return { R, Rmax, inv: 1 / (Math.PI * (modes[members[0]].n ? 1 : 2)), mix };
}

function rectShape(modes, members, x) {
  const vec = new Float64Array(modes[0].vec.length);
  members.forEach((k, i) => {
    const source = modes[k].vec;
    for (let t = 0; t < vec.length; t++) vec[t] += x[i] * source[t];
  });
  return { vec };
}

// ---------- compressible loading of rectangular plates ----------

function radiationGrid(operator, reach) { // spectral rows on a tensor grid that reaches wavenumber `reach`
  const { Px, Py, hx, hy, grids } = operator;
  const key = Math.ceil(reach / 10);
  if (!grids.has(key)) {
    const [u, wu] = spectralNodes(Math.max(30, 1.5 * key * 10 * hx));
    const [v, wv] = spectralNodes(Math.max(30, 1.5 * key * 10 * hy));
    const K = new Float64Array(u.length * v.length);
    u.forEach((x, m) => v.forEach((y, n) => { K[m * v.length + n] = Math.hypot(x / hx, y / hy); }));
    grids.set(key, { wu, wv, gx: spectralRows(Px, u), gy: spectralRows(Py, v), K, Nu: u.length, Nv: v.length });
  }
  return grids.get(key);
}

// (A(κ) − A(0))·x on one parity block, in place into yr, yi: from coefficients to the wavenumber grid, times the symbol
// change ½(√(k²−κ²) − k), and back. Sum factorisation keeps this near 1e6 operations.
function applySymbolChange(grid, kappa, front, pairs, x, y, workspace) {
  const { wu, wv, gx, gy, K, Nu, Nv } = grid;
  const k2 = kappa * kappa;
  const rows = [...new Set(pairs.map(p => p[1]))];
  const { T, F } = workspace;
  for (const part of [0, 1]) {
    const source = x[part];
    T[part].fill(0);
    pairs.forEach(([i, j], q) => {
      const value = source[q], gi = gx[i], offset = rows.indexOf(j) * Nu;
      for (let m = 0; m < Nu; m++) T[part][offset + m] += value * gi[m];
    });
    F[part].fill(0);
    rows.forEach((j, r) => {
      const gj = gy[j], offset = r * Nu;
      for (let m = 0; m < Nu; m++) {
        const value = T[part][offset + m], base = m * Nv;
        for (let n = 0; n < Nv; n++) F[part][base + n] += value * gj[n];
      }
    });
  }
  const [Fr, Fi] = F;
  for (let t = 0; t < K.length; t++) {
    const k = K[t], d = k * k - k2;
    const re = d > 0 ? Math.sqrt(d) - k : -k, im = d > 0 ? 0 : Math.sqrt(-d);
    const fr = Fr[t], fi = Fi[t];
    Fr[t] = re * fr - im * fi;
    Fi[t] = re * fi + im * fr;
  }
  for (const part of [0, 1]) {
    T[part].fill(0);
    rows.forEach((j, r) => {
      const gj = gy[j], offset = r * Nu;
      for (let m = 0; m < Nu; m++) {
        const base = m * Nv;
        let sum = 0;
        for (let n = 0; n < Nv; n++) sum += wv[n] * gj[n] * F[part][base + n];
        T[part][offset + m] = sum;
      }
    });
    y[part].fill(0);
    pairs.forEach(([i, j], q) => {
      const gi = gx[i], offset = rows.indexOf(j) * Nu;
      let sum = 0;
      for (let m = 0; m < Nu; m++) sum += wu[m] * gi[m] * T[part][offset + m];
      y[part][q] = front * sum;
    });
  }
}

/**
 * Γ̃ = [Γ_R, Γ_I] of rectangular mode k at κ = ω·Lx/c, by Galerkin on the Krylov space of A(0)⁻¹(A(κ) − A(0)) from the
 * incompressible solution, which needs no complex matrix. Against an independent vortex-lattice solve this agrees within about
 * 2 % on low modes. On circles, where the full solve exists, six vectors keep Γ_I within 40 % near coincidence on high modes.
 * About 30 ms per mode: Worker only. Well above coincidence 2/√(κ²−k²) is within 5 % and costs nothing.
 */
export function rectRadiation(added, k, kappa, { vectors = 6 } = {}) {
  const { operator, moments } = added;
  const { blocks, Py, hx, hy } = operator;
  const grid = radiationGrid(operator, Math.max(kappa, added.waves[k]));
  const front = Math.PI ** 2 * hx * hy / 2; // as in rectOperator, which also carries the symbol's ½
  let gammaRe = 0, gammaIm = 0;
  for (const { pairs, solve } of blocks) {
    const b = Float64Array.from(pairs, ([i, j]) => moments[k][i * Py + j]);
    if (!b.some(value => Math.abs(value) > 1e-13)) continue;
    const n = pairs.length;
    const rows = new Set(pairs.map(p => p[1])).size;
    const workspace = {
      T: [new Float64Array(rows * grid.Nu), new Float64Array(rows * grid.Nu)],
      F: [new Float64Array(grid.Nu * grid.Nv), new Float64Array(grid.Nu * grid.Nv)],
    };
    // v_1 = A₀⁻¹b, v_{j+1} ∝ A₀⁻¹ΔA v_j, so A₀v_j is known without a product by A₀
    const v = [], A0v = [], dv = [];
    let current = [solve(b.slice()), new Float64Array(n)];
    let image = [b.slice(), new Float64Array(n)];
    for (let j = 0; j < vectors; j++) {
      const applied = [new Float64Array(n), new Float64Array(n)];
      applySymbolChange(grid, kappa, front, pairs, current, applied, workspace);
      v.push(current);
      A0v.push(image);
      dv.push(applied);
      const nr = solve(applied[0].slice()), ni = solve(applied[1].slice());
      let scale = 0;
      for (let q = 0; q < n; q++) scale += nr[q] * nr[q] + ni[q] * ni[q];
      scale = 1 / Math.sqrt(scale || 1);
      current = [nr.map(value => value * scale), ni.map(value => value * scale)];
      image = [applied[0].map(value => value * scale), applied[1].map(value => value * scale)];
    }
    const m = vectors;
    const Tr = new Float64Array(m * m), Ti = new Float64Array(m * m), br = new Float64Array(m), bi = new Float64Array(m);
    for (let r = 0; r < m; r++) {
      for (let q = 0; q < n; q++) {
        br[r] += v[r][0][q] * b[q];
        bi[r] += v[r][1][q] * b[q];
      }
      for (let c = 0; c < m; c++) {
        for (let q = 0; q < n; q++) {
          const ar = A0v[c][0][q] + dv[c][0][q], ai = A0v[c][1][q] + dv[c][1][q];
          Tr[r * m + c] += v[r][0][q] * ar - v[r][1][q] * ai;
          Ti[r * m + c] += v[r][0][q] * ai + v[r][1][q] * ar;
        }
      }
    }
    const [xr, xi] = solveComplex(Tr, Ti, br, bi, m);
    for (let r = 0; r < m; r++) {
      gammaRe += br[r] * xr[r] - bi[r] * xi[r];
      gammaIm += br[r] * xi[r] + bi[r] * xr[r];
    }
  }
  return [gammaRe / added.beta, -gammaIm / added.beta];
}

// ---------- modal damping ----------

/**
 * Thermoelastic relaxation of aluminium from Lambourg, Chaigne & Matignon 2001, Table II: bending stiffness × (1 + sR/(s + c/h²)).
 * Zener's Eα²T(1+ν)/((1−ν)ρc_p) and π²k/(ρc_p) with handbook 6061-T6 values give R 6 % higher and c 23 % lower.
 */
export const ALUMINIUM_THERMAL = { R: 8.45e-3, c: 8.0e-4 };

/** Zener's two constants from handbook values, SI: expansion alpha, conductivity k, heat capacity cp, temperature T. */
export const thermalConstants = ({ E, nu, rho, alpha, k, cp, T = 293 }) => ({
  R: E * alpha * alpha * T * (1 + nu) / ((1 - nu) * rho * cp),
  c: Math.PI ** 2 * k / (rho * cp),
});

/**
 * Share of a circular mode's bending energy in (∇²w)², the only part that makes heat: ∫(∇²w)² / ∫[(∇²w)² − 2(1−ν)(w_xxw_yy − w_xy²)].
 * ∇²J_n(λρ) = −λ²J_n and ∇²I_n(λρ) = λ²I_n, and the denominator is λ⁴∫w² for a mode.
 */
export function circleHeatShare(mode) {
  if (!mode.L) return 0;
  const lam = Math.sqrt(mode.L);
  const [x, w] = gaussLegendre(Math.ceil(2 * lam) + mode.n + 40);
  let top = 0, bottom = 0;
  x.forEach((rho, q) => {
    const [J, , I] = bessel(mode.n, lam * rho);
    top += w[q] * rho * (mode.C * I - J) ** 2;
    bottom += w[q] * rho * (mode.C * I + J) ** 2;
  });
  return top / bottom;
}

/** The same share for every Ritz mode, from the stiffness with ν = 1, which is exactly the (∇²w)² form. */
export function rectHeatShares(ritz, edges, beta) {
  const { modes, nbx, nby } = ritz;
  const [x0, x1, y0, y1] = (typeof edges === 'string' ? [edges, edges, edges, edges] : edges).map(edge => ORDER[edge]);
  const X = basis1D(x0, x1, nbx), Y = basis1D(y0, y1, nby);
  const b2 = beta * beta;
  const scratch = new Float64Array(nbx * nby), out = new Float64Array(nbx * nby);
  // Σ V_ij A_ik B_lj V_kl, with A or B the identity if null
  const form = (vec, A, B) => {
    for (let i = 0; i < nbx; i++) {
      for (let l = 0; l < nby; l++) {
        let sum = 0;
        if (A) for (let k = 0; k < nbx; k++) sum += A[i][k] * vec[k * nby + l];
        else sum = vec[i * nby + l];
        scratch[i * nby + l] = sum;
      }
    }
    for (let i = 0; i < nbx; i++) {
      for (let j = 0; j < nby; j++) {
        let sum = 0;
        if (B) for (let l = 0; l < nby; l++) sum += scratch[i * nby + l] * B[l][j];
        else sum = scratch[i * nby + j];
        out[i * nby + j] = sum;
      }
    }
    return vec.reduce((total, value, t) => total + value * out[t], 0);
  };
  return Float64Array.from(modes, ({ L, vec }) => (L < 1e-3 ? 0
    : (form(vec, X.E22, null) + form(vec, null, Y.E22) / (b2 * b2) + 2 * form(vec, X.E20, Y.E20) / b2) / (L * L)));
}

/** Thermoelastic damping ratio of a mode. Zener's single relaxation; the exact plate series differs by under 2 %. */
export function thermoelasticZeta(w, h, { R, c }, heatShare) {
  const s = c / (h * h);
  return 0.5 * R * heatShare * w * s / (w * w + s * s);
}

/** Viscous damping ratio from the Stokes layers on both faces, where the flow slips past the plate at about the plate's own speed. */
export const viscousZeta = (w, rhoH, { rho, mu }) => (w ? Math.sqrt(mu * rho * w / 2) / (rhoH * w) : 0);

/** Radiation damping ratio from Γ̃ = [Γ_R, Γ_I] and ratio = ρ_f L/(ρh). */
export const radiationZeta = ([gammaR, gammaI], ratio) => ratio * gammaI / (2 * (1 + ratio * gammaR));

/**
 * Physical damping ratio of every mode: thermoelastic + viscous + radiation. Shapes come from the vacuum `modes` (for
 * air they barely change) and frequencies from `loaded`. options: { h, rhoH, L, fluid, thermal, ritz, beta, edges }, the last three for
 * rectangles. Circles: about 3 ms per mode. Rectangles: about 30 ms per mode, so 10 s for 324 modes: run in a Worker.
 */
export function modalDamping(modes, loaded, added, { h, rhoH, L, fluid, thermal, ritz, beta, edges }) {
  const ratio = fluid.rho * L / rhoH;
  const shares = added.circle ? modes.map(circleHeatShare) : rectHeatShares(ritz, edges, beta);
  return Float64Array.from(modes, (mode, k) => {
    const w = loaded[k].w;
    if (w < 1e-6) return 0;
    const kappa = w * L / fluid.c;
    const gamma = !ratio ? [0, 0] : added.circle ? circleRadiation(mode, kappa) : rectRadiation(added, k, kappa);
    return thermoelasticZeta(w, h, thermal, shares[k]) + viscousZeta(w, rhoH, fluid) + radiationZeta(gamma, ratio);
  });
}

// ---------- sand ----------

/**
 * Fractional frequency shift δω/ω = −½Σm_g φ²/m_k of mode k under grains at plate coordinates (flat [u, v, …], as Plate
 * uses them, in [−1, 1]), each of mass grainMass (kg), by first-order perturbation. `inertia` is the share of its mass a grain
 * adds: 1 for grains at rest on the plate; a bouncing layer tends to about 1/Γ (see grainInertia), Γ being the local acceleration in g.
 */
export function sandShift(plate, k, grains, grainMass, inertia = 1) {
  const mode = plate.modes[k];
  let sum = 0;
  if (plate.circle) {
    const theta0 = Math.atan2(plate.dv, plate.du), { R, n } = mode;
    for (let g = 0; g < grains.length; g += 2) {
      const u = grains[g], v = grains[g + 1], rho = Math.min(Math.hypot(u, v), 1) * (NR - 1);
      const j = Math.min(NR - 2, rho | 0);
      sum += (R[j] + (rho - j) * (R[j + 1] - R[j])) ** 2 * Math.cos(n * (Math.atan2(v, u) - theta0)) ** 2;
    }
  } else {
    const { nbx, nby, atX, atY } = plate.ritz;
    for (let g = 0; g < grains.length; g += 2) {
      const X = atX((grains[g] / plate.hx + 1) / 2), Y = atY((grains[g + 1] / plate.hy + 1) / 2);
      let shape = 0;
      for (let i = 0; i < nbx; i++) for (let j = 0; j < nby; j++) shape += mode.vec[i * nby + j] * X[i] * Y[j];
      sum += shape * shape;
    }
  }
  return -0.5 * grainMass * inertia * sum * mode.inv / plate.mass;
}

/** Inelastic bouncing ball limit (Masmoudi et al. 2016, via Ferreyra et al. 2021): apparent mass ∝ 1/Γ above lift-off at Γ = 1. */
export const grainInertia = gamma => (gamma > 1 ? 1 / gamma : 1);
