// Grain transport in real units on a vibrating plate: drift and diffusivity fields from the plate deflection, and a
// portable per-grain update. 1 simulated second = 1 real second. Laws and fits come from event-driven simulations of a
// rigid grain and a boundary-layer solve; see the header of each constant for its origin and how far it was checked.

const G0 = 9.81;

/** Air at room conditions: viscosity, density, sound speed and the mean free path at 1 atm. */
export const AIR = { rho: 1.2, mu: 1.81e-5, c: 343, meanFreePath: 68e-9 };

/**
 * Grain presets. D diameter (m), rho density (kg/m³), e restitution, mu friction, scatter multiplies the diffusivity,
 * arm (m) adds to h/2 for flat-bottomed grains, gammaStick is the Γ below which the grain stays put (> 1 with adhesion),
 * repose angle of repose (degrees), packing solid fraction of a heap.
 * Sources: quartz 2650 kg/m³ and repose 32.5 ± 2.4° (arXiv 1702.05980); NaCl 2170 kg/m³ (CRC Handbook via Wikipedia);
 * lycopodium 35 µm, settling 4.2 cm/s (Loubet et al. 2007, abstract), density 0.9–1.3 g/cm³ (Gregory, secondary).
 * e and mu for sand: DEM values for sand grains, e 0.64–0.94, mu 0.26–0.35 (Heim & Henke 2026). No measurement of sand on
 * plate metal at the mm/s impact speeds found here was located, so e and mu are assumed for that contact; salt and
 * lycopodium copy them. Solder values are Zhou et al. 2016's nominal e, mu, D, H; gammaStick 5 is fitted to the dead zone.
 */
export const SPECIES = {
  sand: { name: 'Sand, 0.3 mm', D: 300e-6, rho: 2650, e: 0.7, mu: 0.35, scatter: 1, arm: 0, gammaStick: 1, repose: 32.5, packing: 0.6 },
  fineSand: { name: 'Fine sand, 0.125 mm', D: 125e-6, rho: 2650, e: 0.7, mu: 0.35, scatter: 1, arm: 0, gammaStick: 1, repose: 32.5, packing: 0.6 },
  salt: { name: 'Table salt, 0.3 mm', D: 300e-6, rho: 2170, e: 0.7, mu: 0.35, scatter: 1, arm: 0, gammaStick: 1, repose: 32.5, packing: 0.6 },
  lycopodium: { name: 'Lycopodium powder, 35 µm', D: 35e-6, rho: 1100, e: 0.5, mu: 0.35, scatter: 1, arm: 0, gammaStick: 1, repose: 32.5, packing: 0.4 },
  flattenedSolder: { name: 'Flattened solder ball, 0.6 mm', D: 600e-6, rho: 7400, e: 0.5, mu: 0.5, scatter: 0.04, arm: 1.28e-4, gammaStick: 5, repose: 32.5, packing: 0.6 },
};

/** Mass of one grain, and how many grains make a given mass of solid material. */
export const grainMass = ({ D, rho }) => (Math.PI / 6) * rho * D ** 3;
export const grainCount = (massKg, grain) => massKg / grainMass(grain);

/**
 * Width of a heap on a nodal line carrying massPerLength (kg/m) at its angle of repose: triangular section, slope tan θ.
 * Derived, not checked against any measured line width.
 */
export const heapWidth = (massPerLength, { rho, packing, repose }) =>
  2 * Math.sqrt(massPerLength / (rho * packing * Math.tan(repose * Math.PI / 180)));

/** Particle relaxation time τ = ρD²C/(18μ), with the Cunningham slip factor C. */
export function relaxationTime({ D, rho }, fluid = AIR) {
  const kn = 2 * fluid.meanFreePath / D;
  const slip = 1 + kn * (1.257 + 0.4 * Math.exp(-1.1 / kn));
  return rho * D * D * slip / (18 * fluid.mu);
}

/** Stokes terminal speed, valid while the Reynolds number on D is below about 1. */
export const settlingSpeed = (grain, fluid = AIR, g = G0) => relaxationTime(grain, fluid) * g;

/** Bending wavenumber k⁴ = ω²ρh/D of a Kirchhoff plate. */
export function bendingWavenumber(omega, { E, rho, nu, h }) {
  const flex = E * h ** 3 / (12 * (1 - nu * nu));
  return Math.sqrt(omega) * (rho * h / flex) ** 0.25;
}

// In-phase impact factor: with sticking friction a grain leaves each bounce with the top-fibre velocity −(h/2)∂ẇ/∂x at the
// moment of impact, and impacts fall while the plate rises. Fitted to rigid polygon runs at Γ = 2 and 5 (0.82, 0.77, 0.68,
// 0.56, 0.54 for e = 0.1…0.9). The quadrature part of ∇W gave no drift within error (|Ψ| < 0.1, Γ ≤ 12).
const impactFactor = e => 0.87 - 0.35 * e;
// Below μ ≈ 0.3 impacts slip and the fibre grips less: 0.46 at μ = 0.05, 0.68 at μ ≥ 0.3 (same runs).
const frictionFactor = mu => Math.min(1, 0.65 + 1.2 * mu);
// Activity rises with Γ from just below 1 to 1 at about 1.65 (Φ/Φ∞ in the same runs, thin and thick plates).
const ramp = gamma => Math.max(0, Math.min(1, (gamma - 0.95) / 0.7));
/** The share of full motion at a drive level: nothing below Γ ≈ 0.95, full above about 1.65. */
export const activityAt = ramp;
// How much the fibre's own acceleration lowers the threshold; a sliding estimate is hk/(2μ), the polygon runs show 0.35 of it.
const FIBRE_THRESHOLD = 0.35;
// Below this power-sum Γ a mix of tones does nothing even when its rare peak passes 1 (8 tones at Γ_peak 1.2, Γ_power 0.42).
const POWER_GATE = 0.75;
// Spread per axis: D = C(Θ)·f(e)·V³/g, Θ = 2V²/(gR) tumbling per flight. C from polygon runs, f(e) from the e sweep.
// Real scatter: ±×3 between runs at a point, ×10 across μ. Not validated against any measurement for sand.
const tumbleFactor = theta => 1 + 8 * theta * theta / (theta * theta + 0.0625);
const hopFactor = e => Math.max(0.35, ((1 + e) / (3 * (1 - e))) ** 3);
const SCATTER_SPEED_SHARE = 0.2; // derived estimate, not measured
// Air: flying grains relax toward the air. x = T_flight/τ from powder runs; brake = mean of e^(−t/τ) over a flight.
const flightsPerTau = (gamma, S) => 1.4 * (gamma / S) ** 0.65;
const brakeOf = x => (x > 1e-4 ? (1 - Math.exp(-x)) / x : 1 - x / 2);
// Rectified first-order air drift (Faraday 1831 ¶25), powder runs: u₁/U ≈ 0.15 (1+S)^−0.6 within about ×3 when lofted.
const RECTIFIED = 0.15;
const RECTIFIED_POWER = -0.6;
const LOFT_SCALE = 1.5; // first-order and streaming drag need the hop to leave the Stokes layer: H = x²/(x²+1.5²), x = hop/δ
const HOP_SCALE = 0.3; // mean hop height ≈ 0.3 min(Vτ, V²/2g), powder runs
const STREAM_COUPLING = 0.5; // grains feel half the slip at their mean hop height, powder runs

let work = null;
const scratchFor = cells => {
  if (work && work.cells === cells) return work;
  const grid = () => new Float32Array(cells);
  work = {
    cells, A: grid(), gradSq: grid(), V2: grid(), G2: grid(), peak: grid(), fibre2: grid(), nx: grid(), ny: grid(),
    ax: grid(), ay: grid(), sx: grid(), sy: grid(),
  };
  return work;
};

/** Central differences of src into the x and y gradient of cell value q, one-sided on the grid boundary. */
const gradX = (src, q, i, G, dx) => (i === 0 ? src[q + 1] - src[q] : i === G - 1 ? src[q] - src[q - 1] : (src[q + 1] - src[q - 1]) / 2) / dx;
const gradY = (src, q, j, G, dx) => (j === 0 ? src[q + G] - src[q] : j === G - 1 ? src[q] - src[q - G] : (src[q + G] - src[q - G]) / 2) / dx;

/**
 * Streaming slip of one drive added into (sx, sy), m/s, from u_s = −(1/4ω)Re[(U*·∇)U + (2−3i)U*(∇·U)], U = iω∇W/κ, ∇²W = −k²W:
 * u_s = (ω/4)(k²/κ²)[∇|W|² − 3 Im(W*∇W)] − (ω/8κ²)∇|∇W|².
 * Coefficients from a scipy solve of the second-order layer in the wall frame and a non-expanded time-domain solve. The
 * Taylor wall term, up to 84× larger, cancels against the Reynolds stress: the wall only translates the layer.
 * Reduces to Rayleigh's (3/8)ωkW² on a standing wave and (3/4)U₀²/c on a travelling one; points to antinodes as Faraday saw.
 * shortcut: ∇²W = −k²W ignores the edge near-field of the modes; upgrade by forming ∇·U from the discrete Laplacian if edges matter.
 */
export function addStreamingSlip({ G, dx, omega, k, kappa, re, im }, sx, sy, gradSq) {
  const c1 = (omega / 4) * (k * k) / (kappa * kappa);
  const c2 = omega / (8 * kappa * kappa);
  for (let j = 0, q = 0; j < G; j++) {
    for (let i = 0; i < G; i++, q++) {
      const rx = gradX(re, q, i, G, dx), ry = gradY(re, q, j, G, dx);
      const ix = gradX(im, q, i, G, dx), iy = gradY(im, q, j, G, dx);
      gradSq[q] = rx * rx + ry * ry + ix * ix + iy * iy;
    }
  }
  for (let j = 0, q = 0; j < G; j++) {
    for (let i = 0; i < G; i++, q++) {
      const rx = gradX(re, q, i, G, dx), ry = gradY(re, q, j, G, dx);
      const ix = gradX(im, q, i, G, dx), iy = gradY(im, q, j, G, dx);
      const a = re[q], b = im[q];
      // ∇|W|² = 2 Re(W*∇W), and the phase term is Im(W*∇W) = a∇b − b∇a
      sx[q] += c1 * (2 * (a * rx + b * ix) - 3 * (a * ix - b * rx)) - c2 * gradX(gradSq, q, i, G, dx);
      sy[q] += c1 * (2 * (a * ry + b * iy) - 3 * (a * iy - b * ry)) - c2 * gradY(gradSq, q, j, G, dx);
    }
  }
}

/**
 * Per-cell transport fields, once per frame. Pure CPU, no allocation after the first call at a grid size.
 *   input: {
 *     G            grid size, cells span [−1, 1]² in plate units,
 *     drives       [{ omega (rad/s), re, im: Float32Array(G²) of complex deflection W in metres, k?: bending wavenumber 1/m }],
 *     plate        { h (m), unit (m per plate unit, half the longest side or the radius), E, rho, nu for k if drives lack it },
 *     grain        a SPECIES entry, fluid (default AIR, pass null for vacuum), g (default 9.81),
 *     edge         optional Plate.edge pairs [outsideCell, insideCell, ...]: amplitude copied outward so edges have no false gradient,
 *   }
 *   out: { ux, uy: Float32Array(G²) drift in plate units per second; D: per-axis diffusivity in plate units²/s; active: Uint8Array }
 * Drift for one tone: u = −Φ(Γ)·(h/2 + arm)·ω·∇|W|·B + air; several tones share by velocity amplitude.
 * Run time for G = 160 and 8 drives: about 2 ms for sand, 7 ms for powder (its streaming needs second differences). No Web Worker needed.
 */
export function transportFields({ G, drives, plate, grain, fluid = AIR, g = G0, edge = null }, out) {
  const cells = G * G;
  const w = scratchFor(cells);
  const { A, gradSq, V2, G2, peak, fibre2, nx, ny, ax, ay, sx, sy } = w;
  V2.fill(0); G2.fill(0); peak.fill(0); fibre2.fill(0); nx.fill(0); ny.fill(0);
  const unit = plate.unit, dx = 2 * unit / G;
  const tau = fluid ? relaxationTime(grain, fluid) : Infinity;
  const dragged = tau < 0.05; // τ above this (quartz over about 75 µm): air currents stay small against friction, so skip their cost
  if (dragged) { ax.fill(0); ay.fill(0); sx.fill(0); sy.fill(0); }
  const h2 = plate.h / 2;
  const fibreGain = FIBRE_THRESHOLD * h2 / grain.mu;

  for (const drive of drives) {
    const { omega, re, im } = drive;
    const om2g = omega * omega / g;
    for (let q = 0; q < cells; q++) A[q] = Math.sqrt(re[q] * re[q] + im[q] * im[q]);
    if (edge) for (let k = 0; k < edge.length; k += 2) A[edge[k]] = A[edge[k + 1]];
    let k = 0, kappa = 0, rect = 0;
    if (dragged) {
      k = drive.k ?? bendingWavenumber(omega, plate);
      // shortcut: κ² = k² − k₀² fails near coincidence, clamped at 0.2k; upgrade with a radiation model for thick plates at high f
      kappa = Math.sqrt(Math.max(k * k - (omega / fluid.c) ** 2, 0.04 * k * k));
      rect = RECTIFIED * (omega / kappa) * (1 + omega * tau) ** RECTIFIED_POWER;
    }
    for (let j = 0, q = 0; j < G; j++) {
      for (let i = 0; i < G; i++, q++) {
        const a = A[q];
        const gx = gradX(A, q, i, G, dx), gy = gradY(A, q, j, G, dx);
        const v = omega * a;
        const gamma = om2g * a;
        V2[q] += v * v;
        G2[q] += gamma * gamma;
        peak[q] += gamma;
        fibre2[q] += om2g * om2g * (gx * gx + gy * gy);
        nx[q] += omega * v * gx;
        ny[q] += omega * v * gy;
        if (dragged) { ax[q] += rect * gx; ay[q] += rect * gy; }
      }
    }
    if (dragged) addStreamingSlip({ G, dx, omega, k, kappa, re, im }, sx, sy, gradSq);
  }

  const { e, mu, D: grainD, gammaStick = 1, scatter = 1, arm = 0 } = grain;
  const phi = impactFactor(e) * frictionFactor(mu);
  const hop = hopFactor(e);
  const radius = grainD / 2;
  const nu = fluid ? fluid.mu / fluid.rho : 0;
  const { ux, uy, D, level } = out;
  for (let q = 0; q < cells; q++) {
    const v2 = V2[q];
    if (v2 === 0) { ux[q] = uy[q] = D[q] = level[q] = 0; continue; }
    const v = Math.sqrt(v2);
    const power = Math.sqrt(G2[q]);
    const fibre = fibreGain * Math.sqrt(fibre2[q]);
    const eff = (power + fibre) / gammaStick;
    const top = (peak[q] + fibre) / gammaStick;
    // A rare peak helps a mix of tones near threshold; the taper keeps a single tone (peak = power) unchanged.
    const taper = Math.max(0, Math.min(1, (1.4 - eff) / 0.55));
    level[q] = eff <= POWER_GATE ? 0 : Math.min(top, eff + Math.min(0.7 * (peak[q] - power), 0.25 * taper));
    let brake = 1, dragU = 0, dragV = 0;
    if (fluid) {
      const omEff = power * g / v;
      brake = brakeOf(flightsPerTau(eff, omEff * tau));
    }
    if (dragged) {
      const omEff = power * g / v;
      const hopHeight = HOP_SCALE * Math.min(v * tau, v2 / (2 * g));
      const loft = hopHeight * hopHeight / (hopHeight * hopHeight + (LOFT_SCALE * LOFT_SCALE) * (2 * nu / omEff));
      const stream = STREAM_COUPLING * (1 - brake) * loft;
      dragU = loft * ax[q] + stream * sx[q];
      dragV = loft * ay[q] + stream * sy[q];
    }
    const drift = -(phi * h2 + arm) * brake / v;
    ux[q] = (drift * nx[q] + dragU) / unit;
    uy[q] = (drift * ny[q] + dragV) / unit;
    const theta = 2 * v2 / (g * radius);
    const hopping = scatter * tumbleFactor(theta) * hop * brake * brake * (v * v2 / g);
    // Air stops a grain's sideways motion within τ, whatever its flight time, so a random walk cannot beat σ²τ, with a
    // sideways speed σ about a fifth of the bounce speed. It only limits grains light enough for air to matter.
    D[q] = Math.min(hopping, SCATTER_SPEED_SHARE ** 2 * v2 * tau) / (unit * unit);
  }
  return out;
}

// Grains cannot overlap. A heap can hold this many layers before its sides slump: the height at which a line one grid cell
// wide stands at about its angle of repose, for the 1–2 mm cells used here. A limit, not a measured value.
const PILE_LAYERS = 3;
const PACKING_FRACTION = 0.8;
const MOST_CROWDED = 4; // above this many times over, a cell spreads no faster: it keeps the random step within a cell

/** How many grains of diameter D fit in one grid cell before they must spread sideways. */
export function cellCapacity({ D }, G, unit) {
  const cellSide = 2 * unit / G;
  return PILE_LAYERS * PACKING_FRACTION * cellSide * cellSide / (Math.PI * D * D / 4);
}

/**
 * How fast a crowd spreads, in plate units²/s per unit of overfilling: about a grid cell's area each second. Grains in an
 * overfull cell jostle at this rate until it is no longer overfull, so only the capacity sets where they end up.
 */
export const crowdDiffusivity = G => (2 / G) ** 2;

/** Counts grains per cell (CPU path) and sets how overfull each cell is. */
export function updateCrowding(x, y, count, { counts, crowd }, G, capacity) {
  counts.fill(0);
  const cellsPerUnit = G / 2;
  for (let p = 0; p < count; p++) {
    const fx = (x[p] + 1) * cellsPerUnit, fy = (y[p] + 1) * cellsPerUnit;
    const i = Math.max(0, Math.min(G - 1, fx | 0)), j = Math.max(0, Math.min(G - 1, fy | 0));
    counts[j * G + i]++;
  }
  for (let q = 0; q < counts.length; q++) crowd[q] = Math.min(MOST_CROWDED, Math.max(0, counts[q] / capacity - 1));
}

/** Allocates the output of transportFields for a G×G grid. */
export const newFields = G => ({
  ux: new Float32Array(G * G), uy: new Float32Array(G * G), D: new Float32Array(G * G), level: new Float32Array(G * G),
  crowd: new Float32Array(G * G), counts: new Float32Array(G * G),
});

const uniform = rng => { // xorshift32
  let s = rng[0];
  s ^= s << 13; s ^= s >>> 17; s ^= s << 5;
  rng[0] = s;
  return (s >>> 0) / 4294967296;
};

const sample = (values, G, ax, ay, o) => // bilinear with weights from stepGrains
  (values[o] + ax * (values[o + 1] - values[o])) * (1 - ay) + (values[o + G] + ax * (values[o + G + 1] - values[o + G])) * ay;

/**
 * Advances count grains by dt seconds, in place, in plate units. x, y are Float32Array in [−1, 1]; rng is a Uint32Array(1)
 * xorshift state (non-zero). Allocation-free and written to translate line by line into a WGSL compute shader:
 * sample the fields, advect (capped at half a cell per step, so a grain cannot jump across a node), add noise, stay inside.
 */
export function stepGrains(x, y, count, { ux, uy, D, level, crowd }, G, hx, hy, circle, dt, rng) {
  const cellsPerUnit = G / 2;
  const end = G - 1.001;
  const maxStep = 0.5 / cellsPerUnit;
  const limitX = hx * 0.995, limitY = hy * 0.995;
  const crowdSpread = crowdDiffusivity(G);
  for (let p = 0; p < count; p++) {
    const px = x[p], py = y[p];
    const fx = Math.max(0, Math.min(end, (px + 1) * cellsPerUnit - 0.5));
    const fy = Math.max(0, Math.min(end, (py + 1) * cellsPerUnit - 0.5));
    const i = fx | 0, j = fy | 0, ax = fx - i, ay = fy - j, o = j * G + i;
    // the gate is decided here for each grain, so a quiet zone has a smooth edge rather than the grid's staircase
    const activity = ramp(sample(level, G, ax, ay, o));
    const crowded = sample(crowd, G, ax, ay, o);
    if (activity === 0 && crowded === 0) continue;
    const dif = activity * sample(D, G, ax, ay, o) + crowded * crowdSpread;
    const u = activity * sample(ux, G, ax, ay, o), v = activity * sample(uy, G, ax, ay, o);
    const kick = Math.sqrt(2 * Math.max(dif, 0) * dt) * 2; // sum of three uniforms has variance 1/4
    const n1 = uniform(rng) + uniform(rng) + uniform(rng) - 1.5;
    const n2 = uniform(rng) + uniform(rng) + uniform(rng) - 1.5;
    let nx = px + Math.max(-maxStep, Math.min(maxStep, u * dt)) + kick * n1;
    let ny = py + Math.max(-maxStep, Math.min(maxStep, v * dt)) + kick * n2;
    if (circle) {
      const r = Math.hypot(nx, ny);
      if (r > 0.995) { const shrink = Math.max(0, 1.99 - r) / r; nx *= shrink; ny *= shrink; }
    } else {
      if (nx > limitX) nx = 2 * limitX - nx; else if (nx < -limitX) nx = -2 * limitX - nx;
      if (ny > limitY) ny = 2 * limitY - ny; else if (ny < -limitY) ny = -2 * limitY - ny;
      nx = Math.max(-limitX, Math.min(limitX, nx));
      ny = Math.max(-limitY, Math.min(limitY, ny));
    }
    x[p] = nx;
    y[p] = ny;
  }
}
