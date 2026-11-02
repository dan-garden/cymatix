import test from 'node:test';
import assert from 'node:assert/strict';
import { Plate } from './physics.js';
import {
  AIR, SPECIES, addStreamingSlip, bendingWavenumber, grainCount, heapWidth, newFields, settlingSpeed, stepGrains, transportFields, activityAt } from './grains.js';

const TAU = 2 * Math.PI;
const g = 9.81;
const G = 64;
const UNIT = 0.02; // 1 plate unit = 2 cm, so one cell is about 0.6 mm
const near = (got, want, factor, msg) => assert.ok(got / want <= factor && got / want >= 1 / factor, `${msg}: got ${got}, want ${want}`);

/** A 1-D standing wave across the grid: W = A0 cos(kx), one drive. */
function standing(A0, f, k, im = null) {
  const re = new Float32Array(G * G), imag = new Float32Array(G * G);
  for (let j = 0, q = 0; j < G; j++) {
    for (let i = 0; i < G; i++, q++) {
      const x = ((i + 0.5) * 2 / G - 1) * UNIT;
      re[q] = A0 * Math.cos(k * x);
      if (im) imag[q] = im(x);
    }
  }
  return { omega: TAU * f, re, im: imag, k };
}
/** The cell-centre position nearest x (metres), and x-drift in m/s and D in m²/s there. */
const cellX = xMetres => {
  const i = Math.round((xMetres / UNIT + 1) * G / 2 - 0.5);
  return { i, x: ((i + 0.5) * 2 / G - 1) * UNIT };
};
function at(fields, xMetres) {
  const q = (G >> 1) * G + cellX(xMetres).i;
  const activity = activityAt(fields.level[q]);
  return { u: activity * fields.ux[q] * UNIT, D: activity * fields.D[q] * UNIT * UNIT, active: activity > 0 ? 1 : 0 };
}
const steel = { E: 200e9, rho: 7800, nu: 0.3 };
const grainOf = (e, mu, extra = {}) => ({ D: 300e-6, rho: 2650, e, mu, scatter: 1, arm: 0, gammaStick: 1, ...extra });
const run = (drives, h, grain, fluid = AIR) => transportFields({ G, drives, plate: { h, unit: UNIT, ...steel }, grain, fluid }, newFields(G));

// Reference: rigid-polygon event simulation in the scratch folder (sim2d.mjs; D = 300 µm, N = 7, roughness 0.2, 64–96 grains, dt = T/128).
// Φ = −u / ((h/2) ω ∇|W|). This is an independent code, not the fits in grains.js.
test('sand drift follows the polygon simulation in Γ, restitution, friction, thickness and frequency', () => {
  const k = 500, { x } = cellX(Math.PI / (4 * k));
  const uRef = (h, omega, G_) => (h / 2) * omega * G_;
  const cases = [ // [Γ local at x, f, h, e, μ, Φ simulated, its standard error], from the sweep2 and tables runs
    [1.2, 1000, 1e-3, 0.5, 0.5, 0.40, 0.023], [1.5, 1000, 1e-3, 0.5, 0.5, 0.62, 0.019], [2, 1000, 1e-3, 0.5, 0.5, 0.65, 0.029],
    [3, 1000, 1e-3, 0.5, 0.5, 0.74, 0.034], [5, 1000, 1e-3, 0.5, 0.5, 0.72, 0.041], [8, 1000, 1e-3, 0.5, 0.5, 0.77, 0.088],
    [5, 4000, 1e-3, 0.5, 0.5, 0.66, 0.034], [5, 400, 1e-3, 0.5, 0.5, 0.90, 0.116],
    [5, 1000, 1e-3, 0.1, 0.5, 0.855, 0.020], [5, 1000, 1e-3, 0.3, 0.5, 0.772, 0.021], [5, 1000, 1e-3, 0.7, 0.5, 0.585, 0.104],
    [5, 1000, 1e-3, 0.5, 0.05, 0.46, 0.032], [5, 1000, 1e-3, 0.5, 0.8, 0.73, 0.021],
    [5, 1000, 0.25e-3, 0.5, 0.5, 0.79, 0.144], [5, 1000, 4e-3, 0.5, 0.5, 0.70, 0.012],
  ];
  for (const [gamma, f, h, e, mu, phiSim, sigma] of cases) {
    const om = TAU * f;
    const A0 = gamma * g / (om * om) / Math.cos(k * x);
    const { u } = at(run([standing(A0, f, k)], h, grainOf(e, mu), null), x);
    const phi = Math.abs(u) / uRef(h, om, A0 * k * Math.sin(k * x));
    assert.ok(Math.abs(phi - phiSim) <= 0.2 * phiSim + 2 * sigma, `Φ at Γ=${gamma} f=${f} h=${h} e=${e} μ=${mu}: got ${phi}, simulated ${phiSim} ± ${sigma}`);
    assert.ok(u > 0, 'drift points down the amplitude gradient');
  }
});

test('sand drift is independent of grain size, as simulated', () => {
  const k = 500, { x } = cellX(Math.PI / (4 * k)), A0 = 5 * g / (TAU * 1000) ** 2 / Math.cos(k * x);
  const a = at(run([standing(A0, 1000, k)], 1e-3, grainOf(0.5, 0.5, { D: 100e-6 }), null), x).u;
  const b = at(run([standing(A0, 1000, k)], 1e-3, grainOf(0.5, 0.5, { D: 1000e-6 }), null), x).u;
  near(a, b, 1.15, 'u at 100 µm and 1 mm');
});

test('the rigid-body result of Zhou et al. 2016 sets the flattened-grain arm', () => {
  // Their Fig. 6a simulator (rigid 600×300 µm block, A = 2.5 µm sin(2πx/6 mm), f = 1 kHz, h = 0), reproduced by sim2d.mjs: arrival times agree
  // within 5–15 % from 0.4–1.4 mm. Read from their figure, the curve starting at 1.2 mm reaches 0.6 mm at about 0.5 s; the closed form
  // arm·ω·∇|W| with no fitted constant beyond the arm measured separately at Γ = 7 should give the same time.
  const k = TAU / 6e-3, om = TAU * 1000, A0 = 2.5e-6;
  const u = x => SPECIES.flattenedSolder.arm * om * A0 * k * Math.abs(Math.cos(k * x));
  let t = 0, x = 1.2e-3;
  while (x > 0.6e-3) { x -= u(x) * 1e-3; t += 1e-3; }
  near(t, 0.5, 1.3, '1.2 mm to 0.6 mm');
});

// Zhou et al. 2016 Fig. 3b, digitised from the published JPEG: note voltage for a 75 % displacement quantile of 175 µm in a 0.5 s
// triangular note, silicon 50×50×0.525 mm, 33 nm/V piezo at the centre. Nothing here is fitted to these numbers.
const NOTES = Array.from({ length: 59 }, (_, n) => 1046.5 * 2 ** (n / 12));
const MEASURED_INV_VOLT = [
  .027, .019, .014, .012, .011, .011, .011, .011, .012, .013, .025, .015, .013, .014, .016, .018, .021, .026, .046, .028, .035, .065,
  .143, .128, .051, .049, .094, .073, .035, .038, .061, .042, .050, .075, .205, .230, .084, .056, .062, .043, .050, .085, .654, .082,
  .048, .043, .045, .050, .062, .122, .205, .442, .130, .107, .574, .136, .219, .091, .073,
];

test('absolute time scale: voltage for a fixed displacement across 59 Chladni notes on silicon (Zhou et al. 2016, Fig. 3b)', () => {
  const gridN = 48;
  const plate = new Plate();
  plate.configure({ shape: 'rect', bc: 'free', edges: ['free', 'free', 'free', 'free'], E: 130e9, rho: 2330, nu: 0.28, h: 525e-6, d: 0.05, Lx: 0.05, Ly: 0.05, zeta: 0.01, fmax: 31000, G: gridN, mount: 'rigid' });
  plate.setDrive(0, 0);
  const grain = SPECIES.flattenedSolder;
  const out = newFields(gridN), re = new Float32Array(gridN * gridN), im = new Float32Array(gridN * gridN);
  const q75 = (f, volts) => {
    const om = TAU * f;
    plate.field([{ f, F: om * om * 33e-9 * volts / g }], g);
    const distances = new Float32Array(gridN * gridN);
    for (const e of [1 / 12, 3 / 12, 5 / 12, 7 / 12, 9 / 12, 11 / 12]) { // midpoints of the triangular envelope, 6 equal bands
      for (let q = 0; q < re.length; q++) re[q] = plate.gamma[q] * e * g / (om * om);
      transportFields({ G: gridN, drives: [{ omega: om, re, im }], plate: { h: 525e-6, unit: 0.025, E: 130e9, rho: 2330, nu: 0.28 }, grain, fluid: null, edge: plate.edge }, out);
      for (let q = 0; q < re.length; q++) distances[q] += activityAt(out.level[q]) * Math.hypot(out.ux[q], out.uy[q]) * 0.025 * (0.5 / 6);
    }
    const inner = [];
    for (let j = 2; j < gridN - 2; j++) for (let i = 2; i < gridN - 2; i++) inner.push(distances[j * gridN + i]);
    inner.sort((a, b) => a - b);
    return inner[Math.floor(inner.length * 0.75)];
  };
  const logRatios = [];
  NOTES.forEach((f, n) => {
    if (MEASURED_INV_VOLT[n] < 0.0115) return; // the amplifier saturated at 90 V on these notes
    let lo = 0.05, hi = 2000;
    for (let it = 0; it < 18; it++) { const u = Math.sqrt(lo * hi); (q75(f, u) > 175e-6) ? hi = u : lo = u; }
    logRatios.push(Math.log10(1 / Math.sqrt(lo * hi) / MEASURED_INV_VOLT[n]));
  });
  const mean = logRatios.reduce((a, b) => a + b, 0) / logRatios.length;
  const spread = Math.sqrt(logRatios.reduce((a, b) => a + (b - mean) ** 2, 0) / logRatios.length);
  assert.ok(Math.abs(mean) < Math.log10(2), `mean ratio ×${(10 ** mean).toFixed(2)}`);
  assert.ok(spread < Math.log10(3.5), `scatter ×${(10 ** spread).toFixed(2)}`);
  // A tilt-only law (arm ≈ A rather than h/2) would miss this by about ×150, measured at 3990 Hz in the same script.
});

test('Zhou et al. 2016 Fig. 3a: displacement grows linearly with note length, so drift dominates scatter', () => {
  // Measured (Fig. 3a left, 3990 Hz, 13.4 V): 0.010, 0.020, 0.044, 0.105 mm at 62.5, 125, 250, 500 ms: exponent 1.13.
  const fit = Math.log(0.105 / 0.010) / Math.log(8);
  // Median plate point there has Γ ≈ 20 (physics.js, 33 nm/V), so V = Γg/ω.
  const om = TAU * 3990, gamma = 19.5, A = gamma * g / (om * om);
  const k = bendingWavenumber(om, { E: 130e9, rho: 2330, nu: 0.28, h: 525e-6 });
  const re = Float32Array.from({ length: G * G }, (_, q) => (A / Math.SQRT1_2) * Math.cos(k * ((q % G + 0.5) * 2 / G - 1) * UNIT));
  const fields = transportFields({ G, drives: [{ omega: om, re, im: new Float32Array(G * G) }], plate: { h: 525e-6, unit: UNIT, E: 130e9, rho: 2330, nu: 0.28 }, grain: SPECIES.flattenedSolder, fluid: null }, newFields(G));
  const { u, D } = at(fields, Math.PI / (4 * k));
  // Péclet over the shortest note: with drift-law exponent 1 the measured 1.13 needs noise well below drift at 62.5 ms
  const t = 0.0625 / 4; // ∫e³ dt for a triangular envelope
  assert.ok(Math.sqrt(2 * D * t) * 1.18 < 15e-6, `scatter over the shortest note is ${Math.sqrt(2 * D * t) * 1.18} m, measured median 10 µm`);
  assert.ok(Math.abs(u) > 0 && fit > 0.8 && fit < 1.3, 'measured exponent near the drift law');
});

test('streaming slip reproduces the scipy boundary-layer solve and the classic limits', () => {
  // streaming_bl.py: wall-frame and fixed-frame (Taylor term −61× Rayleigh cancels) both give ratio 1.0000; time-domain 1.011.
  const f = 1000, om = TAU * f, k = 157, W = 5e-6;
  const dx = 2 * UNIT / G; // k dx ≈ 0.1, so second differences carry about 1 %
  const sx = new Float32Array(G * G), sy = new Float32Array(G * G), tmp = new Float32Array(G * G);
  addStreamingSlip({ G, dx, omega: om, k, kappa: k, re: standing(W, f, k).re, im: new Float32Array(G * G) }, sx, sy, tmp);
  const { i, x } = cellX(Math.PI / (4 * k));
  near(-sx[(G >> 1) * G + i], (3 / 8) * om * k * W * W * Math.sin(2 * k * x), 1.03, 'standing wave: (3/8)ωkW² sin 2kx toward the antinode');
  // Travelling wave W = A e^{−ikx}: classic Rayleigh value (3/4)U₀²/c along the wave, U₀ = ωA
  const re = standing(W, f, k).re;
  const im = standing(0, f, k, xM => -W * Math.sin(k * xM)).im;
  sx.fill(0); sy.fill(0);
  addStreamingSlip({ G, dx, omega: om, k, kappa: k, re, im }, sx, sy, tmp);
  near(sx[(G >> 1) * G + (G >> 1)], (3 / 4) * (om * W) ** 2 * k / om, 1.03, 'travelling wave along +x');
});

// Point-grain powder model, powder.mjs in the scratch folder (first-order Stokes layer, streaming profile, Schiller–Naumann, e = 0.3,
// μ = 0.5, h = 1 mm). Each row: D µm, f, Γ, simulated u at kx = π/4 (negative is toward the antinode). The closed form is a rough fit to
// these runs, so factor tolerances are wide, and sign is the main check.
test('fine powder reverses direction with size and amplitude, as in the powder simulation', () => {
  const rows = [
    [35, 300, 100, -1.1e-2], [35, 300, 300, -2.6e-2], [20, 300, 100, -2.8e-2], [10, 300, 100, -3.9e-2],
    [35, 1000, 300, -2.0e-3], [20, 1000, 300, -1.1e-2],
    [35, 1000, 10, 5.4e-4], [35, 3000, 30, 9.4e-4], [120, 1000, 100, 5.3e-3], [120, 300, 30, 1.6e-3],
  ];
  let scoreLogSum = 0;
  for (const [Dum, f, gam, uSim] of rows) {
    const om = TAU * f, k = 157 * Math.sqrt(f / 1000), { x } = cellX(Math.PI / (4 * k));
    const W = gam * g / (om * om);
    const u = at(run([standing(W, f, k)], 1e-3, { D: Dum * 1e-6, rho: 1100, e: 0.3, mu: 0.5, scatter: 1, arm: 0, gammaStick: 1 }), x).u;
    assert.equal(Math.sign(u), Math.sign(uSim), `sign at D=${Dum} f=${f} Γ=${gam}: got ${u}, sim ${uSim}`);
    scoreLogSum += Math.abs(Math.log10(Math.abs(u / uSim)));
  }
  assert.ok(scoreLogSum / rows.length < Math.log10(4), `mean error ×${(10 ** (scoreLogSum / rows.length)).toFixed(2)}`);
  // Sand on the same plate is friction-driven and goes to the node: simulated +2.9e-3 at Γ = 40, 1 kHz.
  const k = 157, { x } = cellX(Math.PI / (4 * k)), W = 40 * g / (TAU * 1000) ** 2;
  near(at(run([standing(W, 1000, k)], 1e-3, grainOf(0.3, 0.5)), x).u, 2.95e-3, 1.5, 'sand with air');
});

test('Faraday 1831 ¶36: the antinode drift of lycopodium weakens as air thins, reversing at low pressure', () => {
  // Measured: no collection at 6.8–13.5 kPa, collection at 34 kPa. Closed form: Stokes layer grows, the hop stays inside it.
  const f = 300, k = 63, { x } = cellX(Math.PI / (4 * k)), W = 100e-6, grain = SPECIES.lycopodium;
  const atPressure = kPa => at(run([standing(W, f, k)], 1e-3, grain, { ...AIR, rho: AIR.rho * kPa / 101, meanFreePath: AIR.meanFreePath * 101 / kPa }), x).u;
  assert.ok(atPressure(101) < 0, 'toward the antinode at 1 atm');
  assert.ok(atPressure(34) < 0, 'toward the antinode at 34 kPa');
  assert.ok(atPressure(6.8) > atPressure(101), 'weaker or reversed at 6.8 kPa');
});

test('several tones: velocity-amplitude share and thresholds follow the two-tone and eight-tone simulations', () => {
  const k = 500, { x } = cellX(Math.PI / (4 * k));
  const tone = (f, gamma) => standing(gamma * g / (TAU * f) ** 2 / Math.cos(k * x), f, k);
  const u = tones => Math.abs(at(run(tones, 1e-3, grainOf(0.5, 0.5), null), x).u);
  const single = u([tone(1000, 5)]);
  // simulated ratios to one tone Γ = 5 at 1 kHz (−1.74e-3): tones_rolling.mjs
  near(u([tone(1000, 5), tone(1370, 5)]) / single, 2.27 / 1.74, 1.25, 'equal mix');
  near(u([tone(1000, 5), tone(1370, 1)]) / single, 1.68 / 1.74, 1.25, 'weak second tone');
  near(u([tone(1000, 3), tone(2700, 3)]) / single, 1.13 / 1.74, 1.25, 'high second tone');
  // Peak sum 0.8: simulated no motion within error. Peak sum 1.2, power 0.85: slow motion. Eight tones, power 0.42: no motion.
  assert.equal(u([tone(1000, 0.4), tone(1370, 0.4)]), 0, 'below the peak sum');
  assert.ok(u([tone(1000, 0.6), tone(1370, 0.6)]) > 0, 'above the peak sum though the power sum is below 1');
  assert.equal(u([1000, 1137, 1291, 1423, 1619, 1787, 1973, 2111].map(f => tone(f, 0.15))), 0, 'many weak tones: rare peaks do nothing');
});

test('grain update: drift, spread and walls match the fields it was given', () => {
  const gridN = 8, fields = newFields(gridN);
  fields.ux.fill(0.02); fields.D.fill(0.001);
  fields.level.fill(2);
  const n = 40000, dt = 0.01, steps = 20;
  const x = new Float32Array(n), y = new Float32Array(n), rng = new Uint32Array([12345]);
  for (let s = 0; s < steps; s++) stepGrains(x, y, n, fields, gridN, 0.995, 0.995, false, dt, rng);
  const t = steps * dt;
  const mean = x.reduce((a, b) => a + b, 0) / n;
  const variance = x.reduce((a, b) => a + (b - mean) ** 2, 0) / n;
  near(mean, 0.02 * t, 1.1, 'mean displacement = u t');
  near(variance, 2 * 0.001 * t, 1.1, 'variance = 2 D t');
  fields.ux.fill(10); fields.D.fill(0.5);
  fields.level.fill(2);
  for (let s = 0; s < 200; s++) stepGrains(x, y, n, fields, gridN, 0.5, 0.25, false, dt, rng);
  assert.ok(x.every(v => Math.abs(v) <= 0.5) && y.every(v => Math.abs(v) <= 0.25), 'grains stay on a rectangle');
  for (let s = 0; s < 200; s++) stepGrains(x, y, n, fields, gridN, 1, 1, true, dt, rng);
  assert.ok(x.every((v, i) => Math.hypot(v, y[i]) <= 1), 'grains stay on a disc');
});

test('grains collect on the nodal line of a mode from the fields alone', () => {
  // W = W0 sin(kx), node at the centre, antinodes at x = ±1 plate unit. Expected scale: u = Φ (h/2) ω k W0 ≈ 5 × 10⁻⁴ m/s, so a minute crosses the plate.
  const k = Math.PI / (2 * UNIT), W0 = 6 * g / (TAU * 1000) ** 2;
  const drive = standing(W0, 1000, k);
  for (let j = 0, q = 0; j < G; j++) for (let i = 0; i < G; i++, q++) drive.re[q] = W0 * Math.sin(k * ((i + 0.5) * 2 / G - 1) * UNIT);
  const fields = run([drive], 2e-3, grainOf(0.5, 0.5, { scatter: 0.05 }), null);
  const n = 4000, x = Float32Array.from({ length: n }, (_, i) => -0.95 + 1.9 * i / n), y = new Float32Array(n), rng = new Uint32Array([7]);
  const before = x.reduce((s, v) => s + Math.abs(v), 0) / n;
  for (let s = 0; s < 3600; s++) stepGrains(x, y, n, fields, G, 1, 1, false, 1 / 60, rng);
  const after = x.reduce((s, v) => s + Math.abs(v), 0) / n;
  assert.ok(after < 0.5 * before, `mean |x| fell from ${before} to ${after}`);
  // At rest, within the band where Γ < 1 about the node: half-width g/(ω² k W0) in plate units, ±1/6.
  assert.ok(after < 2 / 6, `mean |x| ${after} near the Γ = 1 band edge`);
});

test('presets and counting rules against sourced numbers', () => {
  // Loubet et al. 2007 abstract: lycopodium mean settling 4.2 ± 0.7 cm/s for 35 µm. Density from Gregory's range, not fitted.
  near(settlingSpeed(SPECIES.lycopodium), 0.042, 1.17, 'lycopodium settling');
  // 10 g of 0.3 mm quartz spheres: 10e-3 / (π/6 × 2650 × (3e-4)³), worked by hand: 2.67e5.
  near(grainCount(10e-3, SPECIES.sand), 2.67e5, 1.01, 'grain count');
  // Triangle area = ρ-corrected mass per length; width² tan θ / 4 = M/(ρφ), rearranged by hand for M = 1 g/m on quartz, φ = 0.6, 32.5°.
  near(heapWidth(1e-3, SPECIES.sand), 2 * Math.sqrt(1e-3 / (2650 * 0.6 * Math.tan(32.5 * Math.PI / 180))), 1.0001, 'heap width');
});

test('run time: one field update for G = 160 and 8 drives, and a million grains', () => {
  const gridN = 160, re = new Float32Array(gridN * gridN), im = new Float32Array(gridN * gridN);
  for (let q = 0; q < re.length; q++) { re[q] = 2e-6 * Math.sin(q * 0.01); im[q] = 1e-6 * Math.cos(q * 0.007); }
  const drives = Array.from({ length: 8 }, (_, i) => ({ omega: TAU * (500 + 200 * i), re, im, k: 300 }));
  const plate = { h: 1e-3, unit: 0.1, ...steel };
  const out = newFields(gridN);
  const timeOf = grain => {
    for (let i = 0; i < 20; i++) transportFields({ G: gridN, drives, plate, grain }, out);
    const t0 = performance.now();
    for (let i = 0; i < 20; i++) transportFields({ G: gridN, drives, plate, grain }, out);
    return (performance.now() - t0) / 20;
  };
  const sand = timeOf(SPECIES.sand), powder = timeOf(SPECIES.lycopodium);
  const n = 1_000_000, x = new Float32Array(n), y = new Float32Array(n), rng = new Uint32Array([3]);
  const t0 = performance.now();
  stepGrains(x, y, n, out, gridN, 1, 1, false, 1 / 60, rng);
  const grains = performance.now() - t0;
  console.log(`transportFields G=160, 8 drives: sand ${sand.toFixed(2)} ms, powder ${powder.toFixed(2)} ms; stepGrains 1M: ${grains.toFixed(1)} ms`);
  assert.ok(sand < 15, `sand ${sand} ms`);
});
