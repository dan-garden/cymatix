import test from 'node:test';
import assert from 'node:assert/strict';
import { Plate } from './physics.js';
import { SPECIES, newFields, stepGrains, transportFields } from './grains.js';
import { autoExaggeration, heapHeightPerGrain, peakDeflection, plateHeight, slowMotion, strongestTones } from './plateSurface.js';

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
  return { plate, drive: plate.deflections(1, G0)[0] };
}

const heightsAt = (drive, phase) => plateHeight([{ drive, phase }], new Float32Array(drive.re.length));

test('slowly pushed, the plate follows the force: w = Fa²/16πD (Timoshenko) in phase, and nothing a quarter cycle later', () => {
  const { drive } = drivenAt(0.5, 0);
  const D = clampedDisc.E * clampedDisc.h ** 3 / (12 * (1 - clampedDisc.nu ** 2));
  const exact = 0.12 ** 2 / (16 * Math.PI * D);
  const atStart = heightsAt(drive, 0)[CENTRE];
  assert.ok(Math.abs(atStart / exact - 1) < 5e-3, `${atStart} against ${exact}`);
  assert.ok(Math.abs(heightsAt(drive, TAU / 4)[CENTRE]) < 1e-4 * exact);
});

test('at resonance the plate lags the force by a quarter cycle, and the phase turns with damping', () => {
  const plate = new Plate();
  plate.configure({ ...clampedDisc, zeta: 0.01 });
  const { drive } = drivenAt(plate.modes[0].f, 0.01);
  const inPhase = heightsAt(drive, 0)[CENTRE];
  const lagging = heightsAt(drive, TAU / 4)[CENTRE];
  assert.ok(lagging > 0, 'w = A·sin ωt for a force F·cos ωt');
  assert.ok(Math.abs(inPhase) < 0.05 * lagging, `in phase ${inPhase}, lagging ${lagging}`);
});

test('each tone agrees with the exact Navier series for a simply supported plate, in size and phase, on and off resonance', () => {
  const [a, b, unit] = [0.24, 0.16, 0.12];
  const [E, rho, nu, h, zeta, G] = [69e9, 2700, 0.3, 0.001, 0.01, 41];
  const [u0, v0] = [0.3, -0.2];
  const plate = new Plate();
  plate.keepDrives = true;
  plate.configure({ shape: 'rect', bc: 'simply', E, rho, nu, h, d: a, Lx: a, Ly: b, fmax: 40000, G, zeta });
  plate.setDrive(u0, v0);
  const tones = [220, 660, 1760];
  plate.field(tones.map(f => ({ f, F: 1 })), G0);
  const drives = plate.deflections(1, G0);
  const stiffness = Math.sqrt((E * h ** 3 / (12 * (1 - nu ** 2))) / (rho * h));
  const toX = u => (u + 1) * unit, toY = v => (v + plate.hy) * unit;

  function navier(omega, x, y) { // w = Σ (4F/ρhab) φ(x₀,y₀) φ(x,y) / (ω_mn² − ω² + 2iζω_mn·ω)
    let re = 0, im = 0;
    for (let m = 1; m <= 120; m++) {
      for (let n = 1; n <= 120; n++) {
        const [p, q] = [m * Math.PI / a, n * Math.PI / b];
        const wmn = stiffness * (p * p + q * q);
        const shape = Math.sin(p * toX(u0)) * Math.sin(q * toY(v0)) * Math.sin(p * x) * Math.sin(q * y);
        const [c, d] = [wmn * wmn - omega * omega, 2 * zeta * wmn * omega];
        const scale = 4 * shape / (rho * h * a * b * (c * c + d * d));
        re += scale * c;
        im -= scale * d;
      }
    }
    return [re, im];
  }

  drives.forEach(({ omega, re, im }, tone) => {
    const [largest] = strongestTones([{ omega, re, im }], G0, 1);
    for (const [i, j] of [[10, 20], [30, 14], [20, 26], [33, 22]]) {
      const [u, v] = [(i + 0.5) * 2 / G - 1, (j + 0.5) * 2 / G - 1];
      assert.ok(Math.abs(v) < plate.hy, 'the test point is on the plate');
      const [wr, wi] = navier(omega, toX(u), toY(v));
      const q = j * G + i;
      // relative to the tone's own peak, so a point near a still line is not judged on a tiny number
      const error = Math.hypot(re[q] - wr, im[q] - wi) / largest.peak;
      assert.ok(error < 0.01, `${tones[tone]} Hz at (${u}, ${v}): ${(100 * error).toFixed(2)} % of the peak`);
    }
  });
});

test('the height repeats every cycle and never passes the amplitude', () => {
  const { drive } = drivenAt(700, 0.01);
  const a = heightsAt(drive, 0.9);
  const b = heightsAt(drive, 0.9 + TAU);
  a.forEach((value, q) => assert.ok(Math.abs(value - b[q]) <= 1e-5 * Math.abs(value) + 1e-15));
  a.forEach((value, q) => {
    const amplitude = Math.hypot(drive.re[q], drive.im[q]);
    assert.ok(Math.abs(value) <= amplitude * (1 + 1e-5) + 1e-15);
  });
});

test('several tones add, each at its own phase, and never pass the peak deflection', () => {
  const { drive } = drivenAt(700, 0.01);
  const other = drivenAt(1500, 0.01).drive;
  const tones = [{ drive, phase: 0.4 }, { drive: other, phase: 2.1 }];
  const both = plateHeight(tones, new Float32Array(drive.re.length));
  const [a, b] = [heightsAt(drive, 0.4), heightsAt(other, 2.1)];
  both.forEach((value, q) => assert.ok(Math.abs(value - (a[q] + b[q])) <= 1e-5 * (Math.abs(a[q]) + Math.abs(b[q])) + 1e-15));
  const limit = peakDeflection([drive, other]);
  both.forEach(value => assert.ok(Math.abs(value) <= limit * (1 + 1e-5)));
  assert.ok(plateHeight([], new Float32Array(4)).every(value => value === 0));
});

test('tones are chosen by how hard they shake, not how tall they stand, and listed by frequency', () => {
  const tone = (f, height) => ({ omega: TAU * f, re: Float32Array.of(0, height), im: Float32Array.of(0, 0) });
  // a saw wave's real shape: the tallest tone is not the one shaking hardest
  const chosen = strongestTones([tone(220, 22.7e-6), tone(1760, 0.2e-6), tone(660, 4.7e-6), tone(3000, 1e-9)], G0, 2);
  assert.deepEqual(chosen.map(t => t.frequency.toFixed(0)), ['220', '660']);
  assert.ok(chosen[1].shaking > chosen[0].shaking);
  assert.ok(Math.abs(chosen[1].shaking - (TAU * 660) ** 2 * 4.7e-6 / G0) < 1e-3);
  assert.ok(Math.abs(chosen[0].peak - 22.7e-6) < 1e-12);
  assert.deepEqual(strongestTones([], G0, 4), []);
});

test('the peak phase brings the tallest point to its full height', () => {
  const { drive } = drivenAt(700, 0.01);
  const [{ peak, peakPhase }] = strongestTones([drive], G0, 1);
  const heights = heightsAt(drive, peakPhase);
  let tallest = 0;
  for (const height of heights) tallest = Math.max(tallest, height);
  assert.ok(Math.abs(tallest / peak - 1) < 1e-5, `${tallest} of ${peak}`);
});

test('peak deflection is a tone\'s real peak, and bounds tones that happen to line up', () => {
  const { drive } = drivenAt(700, 0.01);
  assert.ok(Math.abs(peakDeflection([drive]) / strongestTones([drive], G0, 1)[0].peak - 1) < 1e-6);
  assert.equal(peakDeflection([]), 0);
  const size = drive.re.length;
  const other = { omega: TAU * 1100, re: drive.re.map(x => 0.5 * x), im: drive.im.map(x => 0.5 * x) };
  const both = peakDeflection([drive, other]);
  assert.ok(Math.abs(both / (1.5 * peakDeflection([drive])) - 1) < 1e-5, `${both} from ${size} cells`);
});

test('sand settles where the plate stands almost still, so a field line barely moves there, far above chance', () => {
  const G = 81, unit = 0.12, cell = 2 / G;
  const plate = new Plate();
  plate.keepDrives = true;
  plate.configure({ shape: 'rect', bc: 'free', E: 68.9e9, rho: 2700, nu: 0.33, h: 0.001, d: 0.24, Lx: 0.24, Ly: 0.24, fmax: 8000, G, zeta: 0.01 });
  plate.setDrive(0, 0);
  const unitPeak = plate.field([{ f: 981.93, F: 1 }], G0);
  const [drive] = plate.deflections(12 / unitPeak, G0);
  const [{ peak, peakPhase }] = strongestTones([drive], G0, 1);

  // still: nodes of the standing shape, where the tone moves under 10 % of its peak
  const c = Math.cos(peakPhase), s = Math.sin(peakPhase);
  const standing = drive.re.map((re, q) => re * c - drive.im[q] * s);
  const isStill = new Uint8Array(G * G);
  for (let j = 1; j < G - 1; j++) {
    for (let i = 1; i < G - 1; i++) {
      const q = j * G + i;
      const crosses = standing[q] * standing[q + 1] < 0 || standing[q] * standing[q + G] < 0;
      if (crosses && Math.hypot(drive.re[q], drive.im[q]) < 0.1 * peak) isStill[q] = 1;
    }
  }
  const near = new Uint8Array(G * G); // within one cell of a still line
  isStill.forEach((still, q) => {
    if (!still) return;
    const [i, j] = [q % G, Math.floor(q / G)];
    for (let b = Math.max(0, j - 1); b <= Math.min(G - 1, j + 1); b++) {
      for (let a = Math.max(0, i - 1); a <= Math.min(G - 1, i + 1); a++) near[b * G + a] = 1;
    }
  });
  const chance = near.reduce((sum, value) => sum + value, 0) / (G * G);

  const fields = transportFields({ G, drives: [drive], plate: { h: 0.001, unit, E: 68.9e9, rho: 2700, nu: 0.33 }, grain: SPECIES.sand, fluid: null, edge: plate.edge }, newFields(G));
  const count = 3000;
  const rng = Uint32Array.of(2463534242);
  const random = () => { // xorshift32
    rng[0] ^= rng[0] << 13; rng[0] ^= rng[0] >>> 17; rng[0] ^= rng[0] << 5;
    return (rng[0] >>> 0) / 4294967296;
  };
  const x = Float32Array.from({ length: count }, () => 2 * random() - 1), y = Float32Array.from({ length: count }, () => 2 * random() - 1);
  for (let step = 0; step < 1600; step++) stepGrains(x, y, count, fields, G, 1, 1, false, 0.25, rng);
  let onLine = 0;
  for (let p = 0; p < count; p++) {
    const q = Math.floor((y[p] + 1) / cell) * G + Math.floor((x[p] + 1) / cell);
    onLine += near[q];
  }
  const share = onLine / count;
  assert.ok(share > 0.85 && share > 2.5 * chance, `${(100 * share).toFixed(1)} % of sand on lines covering ${(100 * chance).toFixed(1)} % of the plate`);
});

test('slow motion keeps true speed ratios and never outruns the screen', () => {
  assert.equal(slowMotion([1433]), 1433);
  assert.equal(slowMotion([220, 660]), 220, '220 Hz at 1 cycle/s and 660 Hz at 3');
  assert.equal(slowMotion([200, 8000]), 1000, 'fastest at 8 cycles/s, lowest slower than 1');
});

test('rim cells copy the plate, so the height does not drop to zero at the edge', () => {
  const { plate, drive } = drivenAt(700, 0.01);
  const out = plateHeight([{ drive, phase: 0.2 }], new Float32Array(drive.re.length), plate.edge);
  for (let k = 0; k < plate.edge.length; k += 2) assert.equal(out[plate.edge[k]], out[plate.edge[k + 1]]);
});

test('auto exaggeration fills about the share asked, in readable steps, and never shrinks a real shape', () => {
  for (const peak of [3e-6, 7.7e-6, 1e-5, 2.6e-4]) {
    const scale = autoExaggeration(peak, 0.12);
    assert.ok(scale * peak <= 0.3 * 0.12 * (1 + 1e-9), `${peak}: ×${scale}`);
    assert.ok(scale * peak > 0.3 * 0.12 / 2.6, `${peak}: ×${scale}`);
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
