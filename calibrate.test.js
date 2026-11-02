import test from 'node:test';
import assert from 'node:assert/strict';
import { fitPlate } from './calibrate.js';

// Shu et al., Entropy 24:215 (2022), Fig. 2: 25 measured resonances of a 320 × 320 × 1 mm aluminium square on a shaker.
const SHU_MEASURED = [263, 362, 585, 697, 775, 990, 1247, 1373, 1545, 1678, 1938, 2255, 2485, 2808, 2917, 3101, 3202, 3409, 3646, 3879, 4174, 4369, 4481, 4733, 4817];
const SHU_PLATE = {
  shape: 'rect', bc: 'free', mount: 'shaker', E: 70e9, rho: 2700, nu: 0.33, h: 0.001,
  d: 0.32, Lx: 0.32, Ly: 0.32, zeta: 0.005, fmax: 5000, G: 25,
};

test('a fit to Shu et al. recovers the plate and a moving mass well above chance', () => {
  const fit = fitPlate(SHU_MEASURED, SHU_PLATE, { massRatios: [0.03, 0.08, 0.18] });
  assert.equal(fit.massRatio, 0.08, 'about 22 g on a 276 g plate; 16–25 g fits best in physics.test.js');
  assert.ok(Math.abs(fit.scale - 1) < 0.03, `the stated plate should need almost no scaling: ${fit.scale}`);
  assert.ok(fit.rms < 0.02, `rms ${fit.rms}`);
  assert.ok(fit.chance > 2 * fit.rms, `rms ${fit.rms} against chance ${fit.chance}`);
});

test('a plate with the wrong stiffness is recovered through E', () => {
  const wrong = { ...SHU_PLATE, E: 70e9 / 1.21 }; // 10 % low in frequency
  const fit = fitPlate(SHU_MEASURED, wrong, { massRatios: [0.08] });
  assert.ok(Math.abs(fit.scale / 1.1 - 1) < 0.03, `scale ${fit.scale}`);
  assert.ok(Math.abs(fit.E / 70e9 - 1) < 0.06, `E ${fit.E}`);
});

test('resonances of a different plate do not fit better than chance', () => {
  // the same frequencies scrambled by a fixed irregular factor are not a square plate's sequence
  const scrambled = SHU_MEASURED.map((f, i) => f * (1 + 0.09 * Math.sin(2.4 * i)));
  const honest = fitPlate(SHU_MEASURED, SHU_PLATE, { massRatios: [0.08] });
  const fit = fitPlate(scrambled, SHU_PLATE, { massRatios: [0.08] });
  assert.ok(fit.rms > 1.3 * honest.rms, `scrambled ${fit.rms} against real ${honest.rms}`);
});

test('a fit with nothing to match returns null', () => {
  assert.equal(fitPlate([1e7, 2e7], SHU_PLATE, { massRatios: [0.08] }), null);
});
