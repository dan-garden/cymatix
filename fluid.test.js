import test from 'node:test';
import assert from 'node:assert/strict';
import { circleModes, rectModes, Plate, nodalRadii } from './physics.js';
import { nearest, median, sweepResonances } from './testkit.js';
import {
  AIR, WATER, ALUMINIUM_THERMAL, circleAddedMass, rectAddedMass, gammas, loadModes, circleRadiation, rectRadiation,
  circleHeatShare, rectHeatShares, thermalConstants, thermoelasticZeta, viscousZeta, radiationZeta, modalDamping,
  sandShift, grainInertia,
} from './fluid.js';

const close = (got, want, rel, msg) =>
  assert.ok(Math.abs(got - want) <= rel * Math.abs(want), `${msg}: got ${got}, want ${want}`);
const rms = errors => Math.sqrt(errors.reduce((sum, e) => sum + e * e, 0) / errors.length);
const circleGammas = modes => gammas(circleAddedMass(modes), modes.length);
const rectGammas = (edges, nu, beta) => {
  const ritz = rectModes(edges, nu, beta);
  return [ritz, gammas(rectAddedMass(ritz, beta), ritz.modes.length)];
};
const axisymmetric = (modes, s) => modes.findIndex(mode => mode.n === 0 && mode.m === s && mode.L > 0);

test('rigid disc: closed-form added mass and added moment of inertia', () => {
  // μ = (4/π)√(a²−r²) gives (8/3)ρa³ for translation; ring solver below (independent) reaches 0.8488 too.
  const modes = circleModes('free', 0.3, 12);
  const G = circleGammas(modes);
  close(G[modes.findIndex(m => m.L === 0 && m.n === 0)], 8 / (3 * Math.PI), 1e-12, 'translation');
  close(G[modes.findIndex(m => m.L === 0 && m.n === 1)], (16 / 45) / (Math.PI / 4), 1e-12, 'rocking');
});

test('free disc: Kwak (1991) NAVMI factors as printed in Escaler & De La Torre 2018, Table 3 (ν = 1/3, one side)', () => {
  const modes = circleModes('free', 1 / 3, 25);
  const G = circleGammas(modes);
  [0.218247, 0.137003, 0.0962644, 0.07401, 0.06008, 0.05056].forEach((one, i) => {
    close(G[axisymmetric(modes, i + 1)] / 2, one, 1.5e-4, `s${i + 1}`);
  });
});

test('disc modes on all three edges: independent real-space vortex-ring solver in scipy (scratch rings.py), ν = 0.3', () => {
  const ref = {
    clamped: [[0, 0.93348], [1, 0.40649], [2, 0.25290], [3, 0.18242]],
    simply: [[0, 0.98915], [1, 0.39168], [2, 0.24343], [3, 0.17636]],
    free: [[1, 0.43616], [2, 0.27419], [3, 0.19259], [4, 0.14805]],
  };
  for (const [bc, rows] of Object.entries(ref)) {
    const modes = circleModes(bc, 0.3, 20);
    const G = circleGammas(modes);
    rows.forEach(([s, want]) => close(G[axisymmetric(modes, s)], want, 1e-4, `${bc} s=${s}`));
  }
});

test('rectangles: independent vortex-lattice solver in numpy (scratch vlm.py)', () => {
  const [free, Gf] = rectGammas('free', 0.3, 1);
  close(Gf[free.modes.findIndex(m => m.L < 1e-3 && m.px === 0 && m.py === 0)], 0.45473, 2e-4, 'rigid square');
  close(Gf[free.modes.findIndex(m => m.L < 1e-3 && m.px === 1 && m.py === 0)], 0.24076, 3e-4, 'rocking square');
  const [long, Gl] = rectGammas('free', 0.3, 2);
  close(Gl[long.modes.findIndex(m => m.L < 1e-3 && m.px === 0 && m.py === 0)], 0.59441, 3e-4, 'rigid 1 × 2');
  const [ss, Gs] = rectGammas('simply', 0.3, 1); // sin·sin shapes, exact on the lattice side
  const bySum = sum => ss.modes.findIndex(m => Math.abs(m.L - sum * Math.PI ** 2) < 1e-2);
  [[2, 0.53836], [5, 0.31499], [8, 0.24240], [18, 0.15744]].forEach(([sum, want]) => close(Gs[bySum(sum)], want, 3e-4, `simply m² + n² = ${sum}`));
});

test('cantilever: unequal edges against vortex lattice and the strip limit', () => {
  // first bending shape, Γ on the chord b. vlm wall.py with the exact beam shape gives 0.3570 at a/b = 1, 0.6431 at a/b = 5.
  const first = ab => {
    const [ritz, G] = rectGammas(['clamped', 'free', 'free', 'free'], 0.3, 1 / ab);
    return G[ritz.modes.findIndex(m => m.L > 1e-3 && m.py === 0)] * ab;
  };
  close(first(1), 0.3570, 0.01, 'a/b = 1');
  close(first(5), 0.6431, 0.002, 'a/b = 5');
  assert.ok(first(8) < Math.PI / 4, 'stays below the strip limit π/4');
});

test('infinite-plate limit Γ → 2/k: confirmed, approached as 1/k and from different sides for different edges', () => {
  const free = circleModes('free', 0.33, 50), clamped = circleModes('clamped', 0.33, 30);
  const ratio = (modes, G, s) => {
    const k = axisymmetric(modes, s);
    return G[k] * Math.sqrt(modes[k].L) / 2;
  };
  const Gf = circleGammas(free), Gc = circleGammas(clamped);
  close(ratio(free, Gf, 15), 1, 0.025, 'free, k ≈ 47'); // 0.981
  close(ratio(clamped, Gc, 8), 1, 0.08, 'clamped, k ≈ 28'); // 1.07
  assert.ok(ratio(free, Gf, 1) < 0.7, 'the free fundamental is 34 % below it');
  assert.ok(ratio(clamped, Gc, 0) > 1.45, 'the clamped fundamental is 49 % above it');
});

test('compressible loading of discs: compact dipole, infinite-plate limits either side of coincidence', () => {
  const modes = circleModes('free', 0.33, 30);
  const rigid = modes.find(m => m.L === 0 && m.n === 0);
  // far field of a force dipole F = ω·(8/3)ρa³·v gives Γ_I = 16κ³/(27π²)
  for (const kappa of [0.02, 0.05, 0.1]) close(circleRadiation(rigid, kappa)[1], 16 * kappa ** 3 / (27 * Math.PI ** 2), 0.005, `κ = ${kappa}`);
  const high = modes.find(m => m.n === 5 && m.m === 5);
  const k = Math.sqrt(high.L);
  close(circleRadiation(high, 3 * k)[1], 2 / Math.sqrt(9 * k * k - k * k), 0.02, 'resistive above coincidence');
  close(circleRadiation(high, 1.5 * k)[1], 2 / Math.sqrt(2.25 * k * k - k * k), 0.03, 'resistive at 1.5 k');
  close(circleRadiation(high, 0.3 * k)[0], 2 / Math.sqrt(k * k - 0.09 * k * k), 0.05, 'reactive below coincidence');
  assert.ok(Math.abs(circleRadiation(high, 3 * k)[0]) < 0.05 * circleRadiation(high, 3 * k)[1], 'loading has turned resistive');
});

test('compressible loading of rectangles: independent Helmholtz vortex lattice (scratch vlmh.py), unit square', () => {
  const ss = rectModes('simply', 0.3, 1), added = rectAddedMass(ss, 1);
  const index = sum => ss.modes.findIndex(m => Math.abs(m.L - sum * Math.PI ** 2) < 1e-2);
  // [m² + n², κ, Γ_R, Γ_I]
  for (const [sum, kappa, re, im] of [[2, 3, 0.77800, 0.32354], [8, 6, 0.38400, 0.045809], [8, 10, 0.10259, 0.35781], [18, 10, 0.27486, 0.030416], [18, 14, 0.12683, 0.30045]]) {
    const [gr, gi] = rectRadiation(added, index(sum), kappa);
    close(gr, re, 0.015, `Γ_R, ${sum} at κ ${kappa}`);
    close(gi, im, 0.03, `Γ_I, ${sum} at κ ${kappa}`);
  }
});

test('Escaler free plate in water: coupling explains the 2–11 % that NAVMI missed (their Table 5, FEM and theory)', () => {
  const nu = 0.35, h = 8e-4, a = 0.12, rhoH = 2700 * h;
  const D = 69e9 * h ** 3 / (12 * (1 - nu * nu));
  const modes = circleModes('free', nu, 40);
  modes.forEach(mode => { mode.w = mode.L * Math.sqrt(D / rhoH) / (a * a); });
  const added = circleAddedMass(modes);
  const ratio = WATER.rho * a / rhoH;
  const f = loaded => [1, 2, 3, 4, 5, 6].map(s => loaded[axisymmetric(modes, s)].f);
  const fem = [27.2, 139.2, 365.7, 727.5, 1240.8, 1918.7];
  const coupled = f(loadModes(modes, added, ratio)), plain = f(loadModes(modes, added, ratio, { coupled: false }));
  const errors = coupled.map((value, i) => value / fem[i] - 1);
  errors.forEach((e, i) => assert.ok(Math.abs(e) < 0.03, `coupled s${i + 1}: ${(100 * e).toFixed(1)} %`)); // +2.4, +1.3 … +1 %
  plain.forEach((value, i) => assert.ok(value < fem[i], 'uncoupled sits below the FEM, as their "th" row does'));
  close(plain[0], 24.7, 0.02, 'their th row, s1'); // nu differs slightly: they use 1/3 for Γ
});

test('Escaler centre-driven plate, measured in their Table 7: coupling halves the water error, air trims the air error', () => {
  const cfg = { shape: 'circle', bc: 'free', mount: 'shaker', ma: 0.08, E: 69e9, rho: 2700, nu: 0.35, h: 8e-4, d: 0.24, Lx: 0.24, Ly: 0.24, zeta: 0.005, fmax: 5000, G: 25 };
  const plate = new Plate();
  plate.configure(cfg);
  const vacuum = plate.modes, added = circleAddedMass(vacuum), rhoH = 2700 * 8e-4;
  const resonances = (fluid, coupled) => {
    plate.modes = fluid ? loadModes(vacuum, added, fluid.rho * 0.12 / rhoH, { coupled }) : vacuum;
    plate.setDrive(0, 0);
    return sweepResonances(plate, 15, 5000);
  };
  const water = [119, 275, 558, 997, 1594];
  const errors = (found, measured) => measured.map(f => nearest(found, f) / f - 1);
  const coupled = errors(resonances(WATER, true), water), plain = errors(resonances(WATER, false), water);
  assert.ok(rms(coupled) < 0.036, `coupled water rms ${rms(coupled)}`); // 3.4 %: model low on every mode, see the report
  assert.ok(rms(plain) > 0.065, `uncoupled water rms ${rms(plain)}`); // 7.0 %
  const air = [70, 305, 844, 1657, 2747, 4129];
  const before = errors(resonances(null), air), after = errors(resonances(AIR, true), air);
  assert.ok(rms(before) < 0.0095 && rms(after) < 0.0066, `air rms ${rms(before)} → ${rms(after)}`); // 0.92 % → 0.64 %
});

test('Shu et al. with air: no improvement, and the moving mass refit stays at 15–20 g', () => {
  const measured = [263, 362, 585, 697, 775, 990, 1247, 1373, 1545, 1678, 1938, 2255, 2485, 2808, 2917, 3101, 3202, 3409, 3646, 3879, 4174, 4369, 4481, 4733, 4817];
  const cfg = { shape: 'rect', bc: 'free', mount: 'shaker', E: 70e9, rho: 2700, nu: 0.33, h: 1e-3, d: 0.32, Lx: 0.32, Ly: 0.32, zeta: 0.005, fmax: 5000, G: 25, ma: 0.02 };
  const plate = new Plate();
  plate.configure(cfg);
  const added = rectAddedMass(plate.ritz, 1), vacuum = plate.modes, ratio = AIR.rho * 0.32 / 2.7;
  const score = (fluid, ma) => {
    plate.configure({ ...cfg, ma });
    plate.modes = fluid ? loadModes(vacuum, added, ratio) : vacuum;
    plate.setDrive(0, 0);
    const found = sweepResonances(plate, 150, 5200).filter(f => f > 200 && f < 5000);
    return rms(measured.map(f => nearest(found, f) / f - 1));
  };
  close(score(false, 0.02), 0.0174, 0.02, 'vacuum, 20 g');
  close(score(true, 0.015), 0.0186, 0.02, 'air, best moving mass 15 g');
});

test('heat share of bending energy: 1 where the Gauss term vanishes, near 0 for twisting', () => {
  // simply supported sin·sin and clamped shapes have ∫(w_xxw_yy − w_xy²) = 0, so all bending energy is in (∇²w)²
  const ss = rectModes('simply', 0.3, 1);
  rectHeatShares(ss, 'simply', 1).slice(0, 20).forEach((share, i) => close(share, 1, 2e-3, `simply ${i}`));
  const clamped = circleModes('clamped', 0.3, 15);
  clamped.slice(0, 8).forEach(mode => close(circleHeatShare(mode), 1, 1e-6, mode.label));
  // the first elastic mode of a free square is close to w = xy, with ∇²w = 0
  const free = rectModes('free', 0.3, 1);
  const twist = rectHeatShares(free, 'free', 1)[free.modes.findIndex(m => m.L > 1e-3)];
  assert.ok(twist < 0.1, `twist ${twist}`);
});

test('thermoelastic constants: Zener from handbook 6061-T6 against the values Lambourg et al. 2001 used', () => {
  // handbook values from en.wikipedia.org/wiki/6061_aluminium_alloy: k 152 W/mK, c_p 897 J/kgK, α 23.2e-6 /K
  const { R, c } = thermalConstants({ E: 69e9, nu: 0.33, rho: 2700, alpha: 23.2e-6, k: 152, cp: 897 });
  close(R, ALUMINIUM_THERMAL.R, 0.07, 'R'); // +5.6 %
  close(c, ALUMINIUM_THERMAL.c, 0.25, 'c'); // −23 %: their alloy and its conductivity are not stated
});

// Lambourg, Chaigne & Matignon 2001, Fig. 10a, read off a 300 dpi crop of the circles (kHz, 1/s), aluminium plate a3.
const DECAY = [[.0756, .030], [.1105, .0667], [.1455, .2274], [.196, .084], [.2065, .0886], [.343, .2174], [.368, .145], [.378, .230], [.3956, .104], [.4538, .208], [.602, .239],
  [.6196, .252], [.662, .264], [.7375, .244], [.743, .322], [.80, .284], [.881, .403], [.925, .451], [.995, .326], [1.101, .474], [1.156, .376]];

test('measured decay of a free aluminium plate: thermoelastic + radiation + viscous, no parameter fitted here', () => {
  const lx = 0.4195, ly = 0.4, h = 0.004, rhoH = 2660 * h, D1 = 6160e6, D4 = 8600e6; // Table I, plate a3
  const nu = 1 - D4 / (2 * D1), beta = ly / lx, D = D1 * h ** 3;
  const fluid = { rho: 1.2, c: 344, mu: AIR.mu }; // Table II
  const ritz = rectModes('free', nu, beta), added = rectAddedMass(ritz, beta);
  const ratio = fluid.rho * lx / rhoH;
  const G = gammas(added, ritz.modes.length);
  const shares = rectHeatShares(ritz, 'free', beta);
  const modes = ritz.modes.map((mode, k) => ({ k, L: mode.L, w: mode.L * Math.sqrt(D / rhoH) / (lx * lx) / Math.sqrt(1 + ratio * G[k]) })).filter(m => m.L > 1e-3).slice(0, 26);
  modes.forEach(m => {
    m.th = thermoelasticZeta(m.w, h, ALUMINIUM_THERMAL, shares[m.k]) * m.w;
    m.rest = (viscousZeta(m.w, rhoH, fluid) + radiationZeta(rectRadiation(added, m.k, m.w * lx / fluid.c, { vectors: 4 }), ratio)) * m.w;
  });
  const logs = pick => DECAY.map(([f, d]) => {
    const m = modes.reduce((p, q) => (Math.abs(q.w / 6283.2 - f) < Math.abs(p.w / 6283.2 - f) ? q : p));
    return Math.log(pick(m) / d);
  });
  const model = logs(m => m.th + m.rest);
  assert.ok(rms(model) < 0.2, `rms log ${rms(model)}`); // 0.16, a factor 1.17; median 5 %
  assert.ok(median(model.map(Math.abs)) < 0.08, 'median');
  // Chance: the same values shuffled over the modes reach rms log 1.0 (5th percentile 0.77); the best single ζ reaches 0.36.
  const single = logs(m => 7.15e-5 * m.w);
  assert.ok(rms(model) < 0.6 * rms(single), 'beats the best single ζ');
  assert.ok(rms(logs(m => 0.005 * m.w)) > 4, 'ζ = 0.005 is about 70 times too much for this free plate');
});

test('modalDamping on a disc adds up its parts and stays off in vacuum', () => {
  const modes = circleModes('free', 0.33, 12);
  const h = 1e-3, rhoH = 2.7, D = 70e9 * h ** 3 / (12 * (1 - 0.33 ** 2));
  modes.forEach(mode => { mode.w = mode.L * Math.sqrt(D / rhoH) / 0.01; });
  const options = { h, rhoH, L: 0.1, thermal: ALUMINIUM_THERMAL };
  const added = circleAddedMass(modes);
  const zeta = modalDamping(modes, modes, added, { ...options, fluid: AIR });
  const k = axisymmetric(modes, 1), w = modes[k].w;
  const parts = thermoelasticZeta(w, h, ALUMINIUM_THERMAL, circleHeatShare(modes[k])) + viscousZeta(w, rhoH, AIR)
    + radiationZeta(circleRadiation(modes[k], w * 0.1 / AIR.c), AIR.rho * 0.1 / rhoH);
  close(zeta[k], parts, 1e-12, 'sum');
  assert.equal(modalDamping(modes, modes, added, { ...options, fluid: { rho: 0, c: Infinity, mu: 0 } })[k], thermoelasticZeta(w, h, ALUMINIUM_THERMAL, circleHeatShare(modes[k])));
});

test('sand: first-order shift against closed forms', () => {
  const plate = new Plate();
  plate.configure({ shape: 'circle', bc: 'free', E: 69e9, rho: 2700, nu: 0.33, h: 1e-3, d: 0.3, Lx: 0.3, Ly: 0.3, zeta: 0.005, fmax: 3000, G: 25 });
  plate.setDrive(0, 0);
  const k = plate.modes.findIndex(mode => mode.n === 0 && mode.m === 2);
  const plateMass = 2700 * 1e-3 * Math.PI * 0.15 ** 2, sandMass = 0.01;
  const grid = [];
  for (let i = -300; i <= 300; i++) {
    for (let j = -300; j <= 300; j++) {
      const u = (i + 0.37) / 300, v = (j + 0.11) / 300;
      if (u * u + v * v < 1) grid.push(u, v);
    }
  }
  // a uniform layer adds mass everywhere: δω/ω = −½ m_sand/m_plate
  close(sandShift(plate, k, grid, sandMass / (grid.length / 2)), -0.5 * sandMass / plateMass, 0.01, 'uniform layer');
  // grains on a nodal circle add nothing
  const radii = nodalRadii(plate.modes[k]);
  const ring = [];
  for (let q = 0; q < 200; q++) ring.push(radii[0] * Math.cos(q / 32), radii[0] * Math.sin(q / 32));
  assert.ok(Math.abs(sandShift(plate, k, ring, sandMass / 200)) < 1e-6, 'nodal circle');
  assert.equal(grainInertia(0.5), 1);
  assert.equal(grainInertia(4), 0.25);
});

test('fluid off: loadModes hands back the same modes, so nothing in Plate changes', () => {
  const modes = circleModes('free', 0.33, 10);
  assert.equal(loadModes(modes, circleAddedMass(modes), 0), modes);
});
