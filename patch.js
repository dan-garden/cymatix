// Cymatix finite rigid patch: a bolt, washer or shaker head holds a small disk of the plate, not a point.
// The plate couples to a rigid disk through its receptance G(x, x′, ω), written on two rings of radius ≈ r₀ with smooth
// ring forces f = Σ u_i S_i (Galerkin collocation). With T̃ the rigid-body motions seen by each S_i:
//   G̃ u = T̃ q,  q = (w₀, rim displacement of each tilt),  (T̃ᵀ G̃⁻¹ T̃ + B(ω)) q = (F, 0, 0).
// G is the modal sum plus the closed-form tail of the modes above the cut-off (tailSeries). The plate's modal
// coefficients follow from the ring forces. Units are SI; time dependence e^{+iωt}, as in physics.js.

import { radialShapeAt } from './mindlin.js';

const EULER = 0.5772156649015329;
const TAU = 2 * Math.PI;

// A Kirchhoff plate outside a rigid disk feels only w and ∂w/∂r at the rim, so two rings are enough (nothing inside the
// disk is needed). Between rings δ apart the plate can still bend, a slope compliance δ/(3D) per unit rim moment, the
// same as clamping δ/3 further in: placing the rings at r₀(1 + δ/3) and r₀(1 − 2δ/3) cancels it. The error is then O(δ²).
const RING_GAP = 0.1;
const SAMPLES = 16; // points per ring: resolves mode shapes up to n = 13 against harmonics h ≤ 2 without aliasing
const TAIL_SAMPLES = 256; // the tail is sharp over the ring gap, so its ring integrals need arcs well below δ
const TAIL_TERMS = 41; // reaches κ_K·r = 16, where rounding is still about 1e-10 and the tail is down to 2 % of its value at r = 0
const FAR_TAIL_X = 16 * 16 / 4;

/**
 * The tail τ(x, Ω) of the modes above ω_K between two points r apart, where T(r, ω) = τ / (4π√(Dρh)·ω_K),
 * x = κ_K²r²/4, κ_K⁴ = ω_K²ρh/D and Ω = ω/ω_K < 1. A plate's high modes act like an infinite plate's, so
 *   T(r, ω) = 1/(4π√(Dρh)) ∫_{ω_K}^∞ J₀(κ(u)·r)/(u² − ω²) du = Re G∞(r, ω) − (the same integral over [0, ω_K]),
 * with Re G∞ from Y₀ + K₀. Expanding all of it in x, the logarithms of Ω cancel exactly:
 *   τ = Σ_j c_j·(α_j + β_j·ln x),  c_j = x^j/(j!)²,
 *   even j: α_j = −P_j, β_j = 0;  odd j: α_j = Q_j + Ω^{j−1}(2γ − 2H_j), β_j = Ω^{j−1};
 *   P₀ = −atanh(Ω)/Ω, Q₁ = ½ln(1 − Ω²), both obeying P_{j+2} = 1/(j+1) + Ω²·P_j.
 * At x = 0 this is atanh(Ω)/Ω, the point tail ln((ω_K+ω)/(ω_K−ω))/(8πω√(Dρh)) of Plate's #mountFactor.
 * The leading term x·ln x is the biharmonic r²ln r/(8πD).
 */
function tailCoefficients(omega, alpha, beta) {
  const o2 = omega * omega;
  let P = omega > 1e-8 ? -Math.atanh(omega) / omega : -1;
  let Q = 0.5 * Math.log1p(-o2);
  let power = 1, harmonic = 0;
  alpha[0] = -P;
  for (let j = 1; j < TAIL_TERMS; j += 2) {
    harmonic += 1 / j;
    alpha[j] = Q + power * (2 * EULER - 2 * harmonic);
    beta[j] = power;
    Q = 1 / (j + 1) + o2 * Q;
    power *= o2;
    P = 1 / j + o2 * P;
    alpha[j + 1] = -P;
    harmonic += 1 / (j + 1);
  }
}

// c_j(x) added with weight to A, and c_j(x)·ln x to B: the series' dependence on distance, apart from its dependence on Ω.
function addTailMoments(x, weight, A, B) {
  if (!x) {
    A[0] += weight;
    return;
  }
  const logX = Math.log(x);
  let c = weight;
  for (let j = 0; j < TAIL_TERMS; j++) {
    A[j] += c;
    B[j] += c * logX;
    c *= x / ((j + 1) * (j + 1));
  }
}

const dotTail = (A, B, alpha, beta) => {
  let sum = 0;
  for (let j = 0; j < TAIL_TERMS; j++) sum += A[j] * alpha[j] + B[j] * beta[j];
  return sum;
};

/** The tail series at one separation; see tailCoefficients. */
export function tailSeries(x, omega) {
  const alpha = new Float64Array(TAIL_TERMS), beta = new Float64Array(TAIL_TERMS);
  const A = new Float64Array(TAIL_TERMS), B = new Float64Array(TAIL_TERMS);
  tailCoefficients(omega, alpha, beta);
  addTailMoments(x, 1, A, B);
  return dotTail(A, B, alpha, beta);
}

// Complex LU with partial pivoting, in place on split n×n arrays; `pivots` records the row swaps.
function luFactor(re, im, n, pivots) {
  for (let c = 0; c < n; c++) {
    let best = -1, pivot = c;
    for (let r = c; r < n; r++) {
      const size = re[r * n + c] ** 2 + im[r * n + c] ** 2;
      if (size > best) {
        best = size;
        pivot = r;
      }
    }
    pivots[c] = pivot;
    if (pivot !== c) {
      for (let k = 0; k < n; k++) {
        const a = re[c * n + k], b = im[c * n + k];
        re[c * n + k] = re[pivot * n + k];
        im[c * n + k] = im[pivot * n + k];
        re[pivot * n + k] = a;
        im[pivot * n + k] = b;
      }
    }
    const dr = re[c * n + c] / best, di = -im[c * n + c] / best; // 1/pivot
    for (let r = c + 1; r < n; r++) {
      const a = re[r * n + c], b = im[r * n + c];
      const mr = a * dr - b * di, mi = a * di + b * dr;
      re[r * n + c] = mr;
      im[r * n + c] = mi;
      for (let k = c + 1; k < n; k++) {
        const ur = re[c * n + k], ui = im[c * n + k];
        re[r * n + k] -= mr * ur - mi * ui;
        im[r * n + k] -= mr * ui + mi * ur;
      }
    }
  }
}

// Solves LU·x = b in place for the right-hand side stored at x[offset…], y[offset…].
function luSolve(re, im, n, pivots, x, y, offset) {
  for (let c = 0; c < n; c++) {
    const p = pivots[c];
    if (p !== c) {
      const a = x[offset + c], b = y[offset + c];
      x[offset + c] = x[offset + p];
      y[offset + c] = y[offset + p];
      x[offset + p] = a;
      y[offset + p] = b;
    }
    for (let r = c + 1; r < n; r++) {
      const mr = re[r * n + c], mi = im[r * n + c];
      x[offset + r] -= mr * x[offset + c] - mi * y[offset + c];
      y[offset + r] -= mr * y[offset + c] + mi * x[offset + c];
    }
  }
  for (let c = n - 1; c >= 0; c--) {
    const size = re[c * n + c] ** 2 + im[c * n + c] ** 2;
    const dr = re[c * n + c] / size, di = -im[c * n + c] / size;
    const a = x[offset + c], b = y[offset + c];
    const vr = a * dr - b * di, vi = a * di + b * dr;
    x[offset + c] = vr;
    y[offset + c] = vi;
    for (let r = 0; r < c; r++) {
      const ur = re[r * n + c], ui = im[r * n + c];
      x[offset + r] -= ur * vr - ui * vi;
      y[offset + r] -= ur * vi + ui * vr;
    }
  }
}

// Ring-force harmonics the drive can excite, [h, phase]: phase 0 is cos hα and 1 is sin hα about the patch axis.
function excitedHarmonics(plate, u, v) {
  const centred = !u && !v;
  if (plate.circle) return centred ? [[0, 0]] : [[0, 0], [1, 0], [2, 0]]; // the plate is mirror-symmetric about its radius through x₀
  const symmetric = plate.modes.every(mode => mode.px >= 0 && mode.py >= 0);
  if (centred && symmetric) return [[0, 0], [2, 0]]; // even about both centre lines
  return [[0, 0], [1, 0], [1, 1], [2, 0], [2, 1]];
}

/**
 * A rigid disk on a configured Plate. `cfg` is the Plate's configuration. Options: r0 (m; 0 is the point model),
 * mount 'rigid' | 'shaker', ma (kg), inertia (kg·m² about an in-plane axis through the centre), lockTilt (default true
 * for 'rigid', false for 'shaker'), ks (N/m), cs (N·s/m) for a suspended shaker. For tests: tail (default true),
 * harmonics, ringGap, samples.
 *
 * U = rings × harmonics unknowns: 2 on a centred circle, 4 on a centred symmetric rectangle, 6 on any other circle and 10
 * on any other rectangle (half as many again with a third ring). Measured in Node on this machine: coefficients() is
 * O(K·U²/2) + O(U³) with no allocation, about 70 µs at K = 1300, U = 6, so under 2 ms per frame for 24 drives, and 2 µs
 * on a centred circle. setDrive is 0.12 ms per mode on an off-centre circle (K·2·SAMPLES Bessel integrals, 160 ms at
 * K = 1300), 6 ms on a centred one and 40 ms on a rectangle: fine for a click, no Worker needed, not for every frame.
 */
export class RigidPatch {
  constructor(plate, cfg, options) {
    this.plate = plate;
    this.options = options;
    const rhoH = cfg.rho * cfg.h;
    const D = cfg.E * cfg.h ** 3 / (12 * (1 - cfg.nu * cfg.nu));
    this.tailScale = 1 / (4 * Math.PI * Math.sqrt(D * rhoH));
    this.kappa2PerOmega = Math.sqrt(rhoH / D); // κ² = ω·√(ρh/D)
    this.unit = plate.circle ? cfg.d / 2 : Math.max(cfg.Lx, cfg.Ly) / 2; // metres per normalised plate unit
    this.setDrive(plate.du ?? 0, plate.dv ?? 0);
  }

  /** Radius that fits inside the plate at (u, v), so that no ring leaves it. */
  #fittingRadius(u, v) {
    const { plate, unit } = this;
    const { r0 = 0, ringGap = RING_GAP } = this.options;
    const room = plate.circle ? 1 - Math.hypot(u, v) : Math.min(plate.hx - Math.abs(u), plate.hy - Math.abs(v));
    // shortcut: a patch overhanging an edge shrinks to fit; a real washer past the edge is not modelled
    return Math.min(r0, room * unit / (1 + ringGap / 3));
  }

  /** Rebuilds everything that depends on the driver position (u, v), normalised plate coordinates as in Plate.setDrive. */
  setDrive(u, v) {
    const { plate, options } = this;
    this.r0 = this.#fittingRadius(u, v);
    const point = !this.r0;
    const axis = this.axis = plate.circle && (u || v) ? Math.atan2(v, u) : 0;
    const harmonics = point ? [[0, 0]] : options.harmonics ?? excitedHarmonics(plate, u, v);
    // Held at the rim only, the plate inside the disk keeps its own clamped-disk modes, λ² = 10.2 at the inner ring's radius.
    // A real washer has none, so if the modes reach that high a third ring moves them about four times higher.
    const { ringGap = RING_GAP } = options;
    const pair = [1 + ringGap / 3, 1 - 2 * ringGap / 3];
    const innerLambda2 = this.kappa2PerOmega * TAU * plate.fTop * (this.r0 * pair[1]) ** 2;
    const rings = point ? [0] : innerLambda2 > 8 ? [...pair, 0.5] : pair;
    const samples = point ? 1 : options.samples ?? SAMPLES;

    // the sample points, and the force shapes S_i on them: S_i = c_h·cos(hα)/M so that T̃ below is simple
    this.points = [];
    for (const radius of rings) {
      for (let l = 0; l < samples; l++) {
        const alpha = TAU * l / samples;
        const x = (this.r0 / this.unit) * radius * Math.cos(alpha), y = (this.r0 / this.unit) * radius * Math.sin(alpha);
        const pu = u + x * Math.cos(axis) - y * Math.sin(axis), pv = v + x * Math.sin(axis) + y * Math.cos(axis);
        this.points.push({ alpha, radius, u: pu, v: pv });
      }
    }
    this.basis = rings.flatMap((radius, ring) => harmonics.map(([h, phase]) => ({ ring, radius, h, phase })));
    const U = this.basis.length;
    this.U = U;
    this.samples = samples;
    this.rings = rings;
    this.weights = this.basis.map(({ ring, h, phase }) => this.points.map((p, l) => (
      (l / samples | 0) === ring ? (h ? 2 : 1) * (phase ? Math.sin(h * p.alpha) : Math.cos(h * p.alpha)) / samples : 0)));

    // rigid-body motions the disk may take: translation, plus each tilt that a ring force can feel, unless tilts are locked
    const locked = point || (options.lockTilt ?? options.mount === 'rigid');
    this.motions = [0, ...(locked ? [] : harmonics.filter(([h]) => h === 1).map(([, phase]) => 1 + phase))];
    this.dofs = this.motions.length;

    this.#buildShapes(u, v, axis);
    this.#buildTail();
    this.#allocate();
  }

  // ψ_ki = ∫φ_k S_i, for the modes that any S_i can excite, and their products for the packed receptance.
  #buildShapes(u, v, axis) {
    const { plate, U } = this;
    const all = plate.modes.map(plate.circle ? this.#circleProjector(u, v, axis) : this.#rectProjector());
    this.active = [];
    all.forEach((psi, k) => {
      const biggest = psi.reduce((best, x) => Math.max(best, Math.abs(x)), 0);
      if (biggest > 1e-13 * plate.modes[k].Rmax) this.active.push(k);
    });
    const A = this.active.length, pairs = U * (U + 1) / 2;
    this.psi = new Float64Array(A * U);
    this.products = new Float64Array(A * pairs);
    this.active.forEach((k, row) => {
      const psi = all[k], { inv } = plate.modes[k];
      this.psi.set(psi, row * U);
      let t = row * pairs;
      for (let i = 0; i < U; i++) {
        for (let j = i; j < U; j++, t++) this.products[t] = inv * psi[i] * psi[j];
      }
    });
    this.active = Int32Array.from(this.active);
  }

  #circleProjector(u, v, axis) {
    const { points, weights, U } = this;
    const radii = points.map(p => Math.hypot(p.u, p.v));
    const angles = points.map(p => Math.atan2(p.v, p.u) - axis);
    const centred = !u && !v;
    return mode => {
      const psi = new Float64Array(U);
      if (centred && mode.n % this.samples && !this.options.harmonics) return psi; // only h = 0: the integral vanishes
      points.forEach((_, l) => {
        // exact shape, not the 257-point table: nearby points differ by less than its interpolation error
        const value = radialShapeAt(mode, radii[l]) * Math.cos(mode.n * angles[l]);
        for (let i = 0; i < U; i++) psi[i] += weights[i][l] * value;
      });
      return psi;
    };
  }

  #rectProjector() {
    const { plate, points, weights, U } = this;
    const { nbx, nby, atX, atY } = plate.ritz;
    const moments = weights.map(() => new Float64Array(nbx * nby));
    points.forEach((p, l) => {
      const X = atX((p.u / plate.hx + 1) / 2), Y = atY((p.v / plate.hy + 1) / 2);
      for (let i = 0; i < U; i++) {
        if (!weights[i][l]) continue;
        for (let a = 0; a < nbx; a++) {
          for (let b = 0; b < nby; b++) moments[i][a * nby + b] += weights[i][l] * X[a] * Y[b];
        }
      }
    });
    return mode => {
      const psi = new Float64Array(U);
      for (let i = 0; i < U; i++) {
        for (let t = 0; t < nbx * nby; t++) psi[i] += moments[i][t] * mode.vec[t];
      }
      return psi;
    };
  }

  // The tail is rotation-invariant, so S_i and S_j couple only when h and phase agree:
  // ∫∫S_i T S_j = (c_h/M) Σ_Δ T(d_Δ) cos(2πhΔ/M). The moments of the series over Δ are frequency-independent.
  #buildTail() {
    const { basis, rings, r0, plate } = this;
    const perR2 = this.kappa2PerOmega * TAU * plate.fTop / 4;
    const groups = new Map();
    this.tailGroup = [];
    basis.forEach((bi, i) => {
      for (let j = i; j < basis.length; j++) {
        const bj = basis[j];
        if (bi.h !== bj.h || bi.phase !== bj.phase) {
          this.tailGroup.push(-1);
          continue;
        }
        const key = `${bi.ring}|${bj.ring}|${bi.h}`;
        if (!groups.has(key)) {
          const A = new Float64Array(TAIL_TERMS), B = new Float64Array(TAIL_TERMS);
          const a = rings[bi.ring] * r0, b = rings[bj.ring] * r0;
          const count = r0 ? TAIL_SAMPLES : 1;
          for (let delta = 0; delta < count; delta++) {
            const angle = TAU * delta / count;
            const x = perR2 * (a * a + b * b - 2 * a * b * Math.cos(angle));
            if (x < FAR_TAIL_X) addTailMoments(x, (bi.h ? 2 : 1) * Math.cos(bi.h * angle) / count, A, B);
          }
          groups.set(key, { index: groups.size, A, B });
        }
        this.tailGroup.push(groups.get(key).index);
      }
    });
    this.tailMoments = [...groups.values()];
    this.tailValues = new Float64Array(groups.size);
    this.tailAlpha = new Float64Array(TAIL_TERMS);
    this.tailBeta = new Float64Array(TAIL_TERMS);
  }

  #allocate() {
    const { U, dofs, basis, motions } = this;
    const pairs = U * (U + 1) / 2, A = this.active.length;
    this.packedRe = new Float64Array(pairs);
    this.packedIm = new Float64Array(pairs);
    this.Gr = new Float64Array(U * U);
    this.Gi = new Float64Array(U * U);
    this.pivots = new Int32Array(U);
    this.Xr = new Float64Array(dofs * U); // X = G̃⁻¹T̃
    this.Xi = new Float64Array(dofs * U);
    this.T = new Float64Array(dofs * U); // T̃: ring force h = 0 feels the translation, h = 1 the tilts
    basis.forEach(({ radius, h, phase }, i) => {
      motions.forEach((motion, d) => {
        if (motion === 0 && h === 0) this.T[d * U + i] = 1;
        if (motion > 0 && h === 1 && phase === motion - 1) this.T[d * U + i] = radius;
      });
    });
    this.ur = new Float64Array(U); // ring-force amplitudes (N)
    this.ui = new Float64Array(U);
    this.modalRe = new Float64Array(A); // 1/(mass·(ω_k² − ω² + 2iζω_kω))
    this.modalIm = new Float64Array(A);
    this.Ar = new Float64Array(dofs * dofs); // A = T̃ᵀG̃⁻¹T̃, the plate's dynamic stiffness at the disk
    this.Ai = new Float64Array(dofs * dofs);
    this.Sr = new Float64Array(dofs * dofs);
    this.Si = new Float64Array(dofs * dofs);
    this.qr = new Float64Array(dofs);
    this.qi = new Float64Array(dofs);
    this.dofPivots = new Int32Array(dofs);
  }

  /** Fills Gr, Gi with G̃(ω): modal sum over the active modes plus the tail. */
  #assemble(w) {
    const { plate, U, Gr, Gi, packedRe, packedIm, products, active, modalRe, modalIm } = this;
    const { modes, zeta, mass } = plate;
    const pairs = U * (U + 1) / 2;
    packedRe.fill(0);
    packedIm.fill(0);
    for (let row = 0, base = 0; row < active.length; row++, base += pairs) {
      const mode = modes[active[row]];
      const wk = mode.w;
      const a = wk * wk - w * w;
      const b = 2 * (mode.zeta ?? zeta) * wk * w;
      const scale = 1 / (mass * (a * a + b * b));
      const gr = modalRe[row] = a * scale;
      const gi = modalIm[row] = -b * scale;
      for (let t = 0; t < pairs; t++) {
        const product = products[base + t];
        packedRe[t] += gr * product;
        packedIm[t] += gi * product;
      }
    }
    this.#addTail(w);
    for (let i = 0, t = 0; i < U; i++) {
      for (let j = i; j < U; j++, t++) {
        Gr[i * U + j] = Gr[j * U + i] = packedRe[t];
        Gi[i * U + j] = Gi[j * U + i] = packedIm[t];
      }
    }
  }

  #addTail(w) {
    const { plate, options, packedRe, tailMoments, tailValues, tailGroup, tailAlpha, tailBeta } = this;
    if (options.tail === false) return;
    const wK = TAU * plate.fTop;
    if (w >= wK) return; // as for the point mount: the expansion only holds below the cut-off
    tailCoefficients(w / wK, tailAlpha, tailBeta);
    const scale = this.tailScale / wK;
    for (let d = 0; d < tailMoments.length; d++) {
      tailValues[d] = scale * dotTail(tailMoments[d].A, tailMoments[d].B, tailAlpha, tailBeta);
    }
    for (let t = 0; t < tailGroup.length; t++) {
      if (tailGroup[t] >= 0) packedRe[t] += tailValues[tailGroup[t]];
    }
  }

  /** Factors G̃(ω), then forms X = G̃⁻¹T̃ and A = T̃ᵀX. */
  #factor(w) {
    const { U, Gr, Gi, pivots, Xr, Xi, T, Ar, Ai, dofs } = this;
    this.#assemble(w);
    luFactor(Gr, Gi, U, pivots);
    Xr.set(T);
    Xi.fill(0);
    for (let d = 0; d < dofs; d++) luSolve(Gr, Gi, U, pivots, Xr, Xi, d * U);
    for (let a = 0; a < dofs; a++) {
      for (let b = 0; b < dofs; b++) {
        let sr = 0, si = 0;
        for (let i = 0; i < U; i++) {
          sr += T[a * U + i] * Xr[b * U + i];
          si += T[a * U + i] * Xi[b * U + i];
        }
        Ar[a * dofs + b] = sr;
        Ai[a * dofs + b] = si;
      }
    }
  }

  /**
   * Rigid-body motion q. The disk's own dynamic stiffness B is added to A; tilts always obey their row of (A + B)q = 0.
   * 'shaker': the translation row is driven by F. 'rigid': the translation is prescribed, q₀ = F·g/ω².
   */
  #solveMotion(w, F, g) {
    const { options, Ar, Ai, Sr, Si, qr, qi, dofs, dofPivots, r0 } = this;
    const { mount = 'shaker', ma = 0, inertia = 0, ks = 0, cs = 0 } = options;
    Sr.set(Ar);
    Si.set(Ai);
    for (let d = 1; d < dofs; d++) Sr[d * dofs + d] -= w * w * inertia / (r0 * r0); // a rim displacement q is an angle q/r₀
    qr.fill(0);
    qi.fill(0);
    if (mount === 'rigid') {
      Sr.fill(0, 0, dofs);
      Si.fill(0, 0, dofs);
      Sr[0] = 1;
      qr[0] = F * g / (w * w);
    } else {
      Sr[0] += ks - w * w * ma;
      Si[0] += w * cs;
      qr[0] = F;
    }
    luFactor(Sr, Si, dofs, dofPivots);
    luSolve(Sr, Si, dofs, dofPivots, qr, qi, 0);
  }

  // Ring forces u = X·q (N).
  #ringForces() {
    const { U, Xr, Xi, qr, qi, ur, ui, dofs } = this;
    ur.fill(0);
    ui.fill(0);
    for (let d = 0; d < dofs; d++) {
      for (let i = 0; i < U; i++) {
        ur[i] += Xr[d * U + i] * qr[d] - Xi[d * U + i] * qi[d];
        ui[i] += Xr[d * U + i] * qi[d] + Xi[d * U + i] * qr[d];
      }
    }
  }

  /**
   * Fills Plate-style modal coefficients cr + i·ci for one drive at w (rad/s): c_k = (ω²/g)·inv_k·Σ_i u_i ψ_ki /
   * (mass·(ω_k² − ω² + 2iζω_kω)), the acceleration in g per unit mode shape. F is as in Plate.field: newtons for a
   * shaker, the bolt's acceleration in g for 'rigid'. Allocation-free.
   */
  coefficients(w, F, g, cr, ci) {
    const { plate, U, psi, active, modalRe, modalIm, ur, ui } = this;
    this.#factor(w);
    this.#solveMotion(w, F, g);
    this.#ringForces();
    cr.fill(0);
    ci.fill(0);
    const scale = w * w / g;
    for (let row = 0, base = 0; row < active.length; row++, base += U) {
      let sr = 0, si = 0;
      for (let i = 0; i < U; i++) {
        sr += ur[i] * psi[base + i];
        si += ui[i] * psi[base + i];
      }
      const k = active[row];
      const factor = scale * plate.modes[k].inv;
      cr[k] = factor * (modalRe[row] * sr - modalIm[row] * si);
      ci[k] = factor * (modalRe[row] * si + modalIm[row] * sr);
    }
  }

  /** The plate's dynamic stiffness at the disk, A = T̃ᵀG̃⁻¹T̃ (N/m), row-major, dofs × dofs. A pole is a mode of the held plate. */
  stiffness(w) {
    this.#factor(w);
    return [this.Ar, this.Ai];
  }

  /** G̃(ω) between the ring forces (m/N), row-major, U × U: its zeros in ω are the modes of the plate clamped along the disk. */
  receptance(w) {
    this.#assemble(w);
    return [this.Gr, this.Gi];
  }

  /**
   * Response at normalised plate coordinates (u, v) for coefficients from the last coefficients() call, [re, im] in g,
   * including the tail, which matters within about 2/κ_K of the disk and which the plate's own field() leaves out.
   */
  responseAt(u, v, w, g, cr, ci) {
    const { plate, points, weights, ur, ui, U, axis } = this;
    let wr = 0, wi = 0;
    plate.modes.forEach((mode, k) => {
      if (!cr[k] && !ci[k]) return;
      const value = this.#shapeAt(mode, u, v, axis);
      wr += cr[k] * value;
      wi += ci[k] * value;
    });
    const wK = TAU * plate.fTop;
    if (w < wK && this.options.tail !== false) {
      const scale = (w * w / g) * this.tailScale / wK;
      points.forEach((p, l) => {
        const r2 = this.unit ** 2 * ((u - p.u) ** 2 + (v - p.v) ** 2);
        const x = this.kappa2PerOmega * wK * r2 / 4;
        if (x >= FAR_TAIL_X) return;
        // fades out from κ_K·r = 10 to 16, so that the many large, opposed ring forces never cut off at different radii
        const fade = Math.min(1, (16 - 2 * Math.sqrt(x)) / 6);
        const tail = fade * scale * tailSeries(x, w / wK); // shortcut: allocates; only analysis calls this, field() does not
        for (let i = 0; i < U; i++) {
          wr += ur[i] * weights[i][l] * tail;
          wi += ui[i] * weights[i][l] * tail;
        }
      });
    }
    return [wr, wi];
  }

  #shapeAt(mode, u, v, axis) {
    const { plate } = this;
    if (plate.circle) {
      const angle = Math.atan2(v, u) - axis;
      return radialShapeAt(mode, Math.hypot(u, v)) * Math.cos(mode.n * angle);
    }
    const { nbx, nby, atX, atY } = plate.ritz;
    const X = atX((u / plate.hx + 1) / 2), Y = atY((v / plate.hy + 1) / 2);
    let value = 0;
    for (let a = 0; a < nbx; a++) {
      for (let b = 0; b < nby; b++) value += mode.vec[a * nby + b] * X[a] * Y[b];
    }
    return value;
  }
}
