import test from 'node:test';
import assert from 'node:assert/strict';
import { Plate } from './physics.js';
import { autoExaggeration, deflectionAt, heapHeightPerGrain, peakDeflection } from './plateSurface.js';

const G0 = 9.81;
const TAU = 2 * Math.PI;
const clampedDisc = { shape: 'circle', bc: 'clamped', E: 69e9, rho: 2700, nu: 0.3, h: 0.001, d: 0.24, Lx: 0.24, Ly: 0.24, fmax: 40000, G: 41 };
const CENTRE = (41 * 41 - 1) / 2;

function drivenAt(f, zeta) {
  const plate = new Plate();
  plate.keepDrives = true;
  plate.configure({ ...clampedDisc, zeta });
  plate.setDrive(0, 0);
  plate.field([{ f, F: 1 }], G0);
  return { plate, drives: plate.deflections(1, G0) };
}

test('slowly pushed, the plate follows the force: w = Fa²/16πD (Timoshenko) in phase, and nothing a quarter cycle later', () => {
  const f = 0.5;
  const { drives } = drivenAt(f, 0);
  const D = clampedDisc.E * clampedDisc.h ** 3 / (12 * (1 - clampedDisc.nu ** 2));
  const exact = 0.12 ** 2 / (16 * Math.PI * D);
  const out = new Float32Array(41 * 41);
  const atStart = deflectionAt(drives, 0, out)[CENTRE];
  assert.ok(Math.abs(atStart / exact - 1) < 5e-3, `${atStart} against ${exact}`);
  assert.ok(Math.abs(deflectionAt(drives, 1 / (4 * f), out)[CENTRE]) < 1e-4 * exact);
});

test('at resonance the plate lags the force by a quarter cycle, and the phase turns with damping', () => {
  const plate = new Plate();
  plate.configure({ ...clampedDisc, zeta: 0.01 });
  const f = plate.modes[0].f;
  const { drives } = drivenAt(f, 0.01);
  const out = new Float32Array(41 * 41);
  const inPhase = deflectionAt(drives, 0, out)[CENTRE];
  const lagging = deflectionAt(drives, 1 / (4 * f), out)[CENTRE];
  assert.ok(lagging > 0, 'w = A·sin ωt for a force F·cos ωt');
  assert.ok(Math.abs(inPhase) < 0.05 * lagging, `in phase ${inPhase}, lagging ${lagging}`);
});

test('shape repeats every period, and tones add linearly', () => {
  const { drives } = drivenAt(700, 0.01);
  const [one] = drives;
  const other = { omega: TAU * 1100, re: one.re.map(x => 0.5 * x), im: one.im.map(x => -0.3 * x) };
  const size = one.re.length;
  const period = TAU / one.omega;
  const a = deflectionAt([one], 0.0003, new Float32Array(size));
  const b = deflectionAt([one], 0.0003 + period, new Float32Array(size));
  a.forEach((value, q) => assert.ok(Math.abs(value - b[q]) <= 1e-5 * Math.abs(value) + 1e-15));
  const sum = deflectionAt([one, other], 0.0002, new Float32Array(size));
  const first = deflectionAt([one], 0.0002, new Float32Array(size));
  const second = deflectionAt([other], 0.0002, new Float32Array(size));
  sum.forEach((value, q) => assert.ok(Math.abs(value - first[q] - second[q]) <= 1e-5 * Math.abs(value) + 1e-15));
});

test('peakDeflection bounds the shape at every time and is reached by a single tone', () => {
  const { drives } = drivenAt(700, 0.01);
  const peak = peakDeflection(drives);
  const out = new Float32Array(drives[0].re.length);
  let largest = 0;
  for (let t = 0; t < 64; t++) {
    for (const value of deflectionAt(drives, t / (64 * 700), out)) {
      largest = Math.max(largest, Math.abs(value));
    }
  }
  assert.ok(largest <= peak * (1 + 1e-6));
  assert.ok(largest > 0.99 * peak, `reached ${largest} of ${peak}`);
  assert.equal(peakDeflection([]), 0);
});

test('rim cells copy the plate, so the surface does not drop to zero at the edge', () => {
  const { plate, drives } = drivenAt(700, 0.01);
  const out = deflectionAt(drives, 0.0001, new Float32Array(drives[0].re.length), plate.edge);
  for (let k = 0; k < plate.edge.length; k += 2) assert.equal(out[plate.edge[k]], out[plate.edge[k + 1]]);
});

test('auto exaggeration fills about the share asked, in readable steps, and never shrinks a real shape', () => {
  for (const peak of [3e-6, 7.7e-6, 1e-5, 2.6e-4]) {
    const scale = autoExaggeration(peak, 0.12);
    assert.ok(scale * peak <= 0.15 * 0.12 * (1 + 1e-9), `${peak}: ×${scale}`);
    assert.ok(scale * peak > 0.15 * 0.12 / 2.6, `${peak}: ×${scale}`);
    const leading = scale / 10 ** Math.floor(Math.log10(scale));
    assert.ok([1, 2, 5].some(step => Math.abs(leading - step) < 1e-9), `×${scale}`);
  }
  assert.equal(autoExaggeration(0, 0.12), 1);
  assert.equal(autoExaggeration(0.5, 0.12), 1);
});

test('heap height conserves the grains\' volume', () => {
  const sand = { D: 300e-6, packing: 0.6 };
  const G = 160, unit = 0.12, grains = 100000;
  const perGrain = heapHeightPerGrain(sand, G, unit);
  const cell = (2 * unit / G) ** 2;
  const solid = grains * (Math.PI / 6) * sand.D ** 3;
  assert.ok(Math.abs(grains * perGrain * cell * sand.packing / solid - 1) < 1e-12);
  // 60 grains in a 1.5 mm cell stand about two grains high, a sane loose heap
  const height = 60 * perGrain;
  assert.ok(height > 1.5 * sand.D && height < 2.5 * sand.D, `${height} m`);
});
