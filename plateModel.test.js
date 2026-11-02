// The full plate model working together: each module has its own tests against sources, and these check the wiring.
import test from 'node:test';
import assert from 'node:assert/strict';
import { Plate } from './physics.js';
import { PlateModel, stiffeningRegime } from './plateModel.js';
import { AIR, WATER } from './fluid.js';
import { sweepResonances } from './calibrate.js';
import { solveNow, solvePlan } from './solvePlan.js';
import { AIR as GRAIN_AIR, SPECIES, newFields, stepGrains, transportFields } from './grains.js';
import { mindlinOmega } from './physics.js';
import { nearest, pct } from './testkit.js';

const G0 = 9.81;
const ALUMINIUM = { E: 69e9, rho: 2700, nu: 0.33 };
// Escaler & De La Torre 2018: 240 × 0.8 mm, driven through a bar with 0.08 kg moving mass
const ESCALER = { shape: 'circle', bc: 'free', mount: 'shaker', ma: 0.08, E: 69e9, rho: 2700, nu: 0.35, h: 0.0008, d: 0.24, Lx: 0.24, Ly: 0.24, zeta: 0.005, fmax: 5000, G: 25 };
const ESCALER_MEASURED = [70, 305, 844, 1657, 2747, 4129];
// Shu et al. 2022: 320 × 1 mm on a shaker, material values from the paper
const SHU = { shape: 'rect', bc: 'free', mount: 'shaker', ma: 0.015, E: 70e9, rho: 2700, nu: 0.33, h: 0.001, d: 0.32, Lx: 0.32, Ly: 0.32, zeta: 0.005, fmax: 5000, G: 25 };
const SHU_MEASURED = [263, 362, 585, 697, 775, 990, 1247, 1373, 1545, 1678, 1938, 2255, 2485, 2808, 2917, 3101, 3202, 3409, 3646, 3879, 4174, 4369, 4481, 4733, 4817];

const rms = values => Math.sqrt(values.reduce((sum, x) => sum + x * x, 0) / values.length);
const configured = cfg => {
  const plate = new PlateModel();
  plate.configure(cfg);
  plate.setDrive(0, 0);
  return plate;
};

test('with every option off, the model is exactly the core plate', () => {
  const cfg = { shape: 'rect', bc: 'free', ...ALUMINIUM, h: 0.001, d: 0.3, Lx: 0.3, Ly: 0.2, zeta: 0.01, fmax: 4000, G: 32, mount: 'shaker', ma: 0.02 };
  const core = new Plate(), model = new PlateModel();
  core.configure(cfg);
  model.configure(cfg);
  for (const plate of [core, model]) plate.setDrive(0.2, -0.1);
  const drives = [{ f: 700, F: 1 }, { f: 2100, F: 0.3 }];
  assert.equal(model.field(drives, G0), core.field(drives, G0));
  assert.deepEqual(model.gamma, core.gamma);
});

test('air lowers each mode by a little, and water by a lot, as the fluid module predicts', () => {
  const cfg = { ...ESCALER, mount: 'force', fmax: 3000 };
  const vacuum = configured(cfg).modes.map(mode => mode.f);
  const inAir = configured({ ...cfg, fluid: AIR }).modes.map(mode => mode.f);
  const inWater = configured({ ...cfg, fluid: WATER }).modes.map(mode => mode.f);
  for (let k = 2; k < 20; k++) {
    assert.ok(inAir[k] < vacuum[k] && inAir[k] > 0.97 * vacuum[k], `air, mode ${k}: ${inAir[k]} vs ${vacuum[k]}`);
    assert.ok(inWater[k] < 0.75 * vacuum[k], `water, mode ${k}: ${inWater[k]} vs ${vacuum[k]}`);
  }
});

test('Escaler & De La Torre stays within 2 % with air included and nothing fitted', () => {
  const resonances = sweepResonances(configured({ ...ESCALER, fluid: AIR }), 40, 5000);
  ESCALER_MEASURED.forEach((f, i) => assert.ok(pct(resonances[i], f) < 2, `s=${i + 1}: ${resonances[i]} vs ${f}`));
});

test('Shu et al. stays below 2 % RMS with air and the moving mass refitted to 15 g', () => {
  const resonances = sweepResonances(configured({ ...SHU, fluid: AIR }), 150, 5200);
  const errors = SHU_MEASURED.map(f => nearest(resonances, f) / f - 1);
  assert.ok(rms(errors) < 0.02, `rms ${rms(errors)}`);
});

test('a finite clamp and air work together, and the clamp raises the bolted resonances', async () => {
  const point = sweepResonances(configured({ ...ESCALER, fluid: AIR }), 40, 5000);
  const clamped = sweepResonances(configured({ ...ESCALER, fluid: AIR, clampRadius: 0.004 }), 40, 5000);
  assert.equal(clamped.filter(f => f < 4500).length, 6);
  // fixed in patch.js by exact annular plates; here it must stay the same size when air is added
  ESCALER_MEASURED.forEach((_, i) => {
    const rise = clamped[i] / point[i] - 1;
    assert.ok(rise > 0.01 && rise < 0.08, `s=${i + 1}: ${(100 * rise).toFixed(1)} %`);
  });
});

test('a plate solved as true Mindlin replaces the per-mode correction', async () => {
  const cfg = { shape: 'rect', bc: 'clamped', ...ALUMINIUM, nu: 0.3, h: 0.02, d: 0.2, Lx: 0.2, Ly: 0.2, zeta: 0.005, fmax: 8000, G: 16 };
  const plan = solvePlan(cfg);
  assert.equal(plan.thick, true);
  const plate = new PlateModel();
  plate.configure({ ...cfg, solved: await solveNow(cfg) });
  // Λ₁ = 32.49 for clamped h/a = 0.1 in true Mindlin, 3D 32.74 (mindlin.test.js); the old correction gives about 31.4
  const D = cfg.E * cfg.h ** 3 / (12 * (1 - cfg.nu ** 2));
  const scale = Math.sqrt(D / (cfg.rho * cfg.h)) / cfg.Lx ** 2;
  assert.ok(pct(plate.modes[0].w / scale, 32.49) < 0.2, `Λ₁ ${plate.modes[0].w / scale}`);
  const corrected = mindlinOmega(35.985 * scale, cfg) / scale;
  assert.ok(pct(corrected, 32.74) > pct(plate.modes[0].w / scale, 32.74), `old ${corrected} vs new ${plate.modes[0].w / scale}, 3D 32.74`);
});

test('plan: thin plates stay on the fast path, and thick or large ones go to a worker', () => {
  const thin = { shape: 'rect', bc: 'free', ...ALUMINIUM, h: 0.001, d: 0.3, Lx: 0.3, Ly: 0.3, zeta: 0.01, fmax: 8000, G: 32 };
  assert.equal(solvePlan(thin), null);
  assert.equal(solvePlan({ ...thin, h: 0.005 }).thick, true);
  assert.equal(solvePlan({ ...thin, Lx: 1, Ly: 1 }).task, 'rectBlock');
  assert.equal(solvePlan({ ...thin, shape: 'circle', d: 1, h: 0.0005 }).task, 'circleModes');
  assert.equal(solvePlan({ ...thin, h: 0.005, mindlin: false }), null);
});

test('stiffening lifts the resonance of a force-driven plate only when it is on, and only at large amplitude', () => {
  const cfg = { shape: 'circle', bc: 'free', mount: 'force', ...ALUMINIUM, h: 0.0003, d: 0.24, Lx: 0.24, Ly: 0.24, zeta: 0.005, fmax: 600, G: 16 };
  const mode = configured(cfg).modes.find(m => m.n === 2 && m.m === 0);
  const peakFrequency = (forceGain, stiffening) => {
    const plate = new PlateModel();
    plate.configure({ ...cfg, stiffening });
    plate.setDrive(0.5, 0);
    if (stiffening) plate.prepareStiffening(cfg);
    plate.forceGain = forceGain;
    let best = 0, bestF = 0;
    for (let f = mode.f * 0.97; f < mode.f * 1.2; f *= 1.0005) {
      const value = plate.field([{ f, F: forceGain }], G0);
      if (value > best) { best = value; bestF = f; }
    }
    return { bestF, deflection: best * G0 / (2 * Math.PI * bestF) ** 2 };
  };
  const base = peakFrequency(1e-3, false).bestF;
  assert.ok(pct(peakFrequency(1e-3, true).bestF, base) < 0.1, 'small forces do not stiffen');
  const strong = peakFrequency(0.2, true);
  assert.ok(strong.bestF > base * 1.005, `${strong.bestF} vs ${base}, deflection ${strong.deflection}`);
  // the state does not leak: off means linear at any force
  assert.ok(pct(peakFrequency(0.2, false).bestF, base) < 0.1);
});

test('the stiffening regime follows the thresholds in nonlinear.js', () => {
  const omega = 2 * Math.PI * 100;
  assert.equal(stiffeningRegime(5 * G0, 2 * Math.PI * 1000, 0.001, 0.005), 'negligible');
  assert.equal(stiffeningRegime(100 * G0, omega, 0.001, 0.005), 'broken');
});

test('deflections are the same field as gamma, in metres', () => {
  const plate = new PlateModel();
  plate.keepDrives = true;
  plate.configure({ shape: 'rect', bc: 'free', mount: 'force', ...ALUMINIUM, h: 0.001, d: 0.3, Lx: 0.3, Ly: 0.3, zeta: 0.01, fmax: 3000, G: 32 });
  plate.setDrive(0, 0);
  const f = 987;
  plate.field([{ f, F: 1 }], G0);
  const gamma = Float32Array.from(plate.gamma);
  const omega = 2 * Math.PI * f;
  const [drive] = plate.deflections(3, G0);
  for (let q = 0; q < gamma.length; q += 37) {
    if (plate.edge.includes(q)) continue;
    const metres = Math.hypot(drive.re[q], drive.im[q]);
    assert.ok(Math.abs(metres * omega * omega / G0 - 3 * gamma[q]) <= 1e-3 * (3 * gamma[q]) + 1e-9, `cell ${q}`);
  }
});

test('the whole chain collects sand onto the nodal lines, in real time units', () => {
  const plate = new PlateModel();
  plate.keepDrives = true;
  const cfg = { shape: 'rect', bc: 'free', mount: 'force', ...ALUMINIUM, h: 0.001, d: 0.3, Lx: 0.3, Ly: 0.3, zeta: 0.01, fmax: 3000, G: 64, fluid: AIR };
  plate.configure(cfg);
  plate.setDrive(0, 0);
  const f = plate.modes.find(m => m.label === '#20').f;
  const unit = plate.field([{ f, F: 1 }], G0);
  const gain = 12 / unit;
  const fields = newFields(plate.G);
  const quietShare = (x, y) => {
    let quiet = 0;
    for (let p = 0; p < x.length; p++) {
      const i = Math.min(63, Math.floor((x[p] + 1) * 32)), j = Math.min(63, Math.floor((y[p] + 1) * 32));
      if (plate.gamma[j * 64 + i] * gain < 2) quiet++;
    }
    return quiet / x.length;
  };
  const count = 4000, x = new Float32Array(count), y = new Float32Array(count);
  let seed = 11;
  const random = () => (seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296;
  for (let p = 0; p < count; p++) {
    x[p] = random() * 1.9 - 0.95;
    y[p] = random() * 1.9 - 0.95;
  }
  const before = quietShare(x, y);
  const grain = SPECIES.sand;
  const rng = Uint32Array.of(99);
  for (let t = 0; t < 400; t++) { // 200 s of plate time
    plate.field([{ f, F: 1 }], G0);
    transportFields({ G: 64, drives: plate.deflections(gain, G0), plate: { h: cfg.h, unit: 0.15, E: cfg.E, rho: cfg.rho, nu: cfg.nu }, grain, fluid: GRAIN_AIR, g: G0, edge: plate.edge }, fields);
    stepGrains(x, y, count, fields, 64, plate.hx, plate.hy, false, 0.5, rng);
  }
  const after = quietShare(x, y);
  assert.ok(after > 0.8 && after > before + 0.3, `after 200 s, ${after} of grains are where Γ < 2 (started at ${before})`);
});
