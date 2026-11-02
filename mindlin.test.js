// True Mindlin modes against closed forms, independent exact and 3D solutions (scipy, see each test) and published tables.
import test from 'node:test';
import assert from 'node:assert/strict';
import { Plate, circleModes, gaussLegendre, mindlinOmega, nodalRadii, rectModes, bessel, circleLamFor } from './physics.js';
import { mindlinCircleFields, mindlinCircleModes, mindlinNodalRadii, mindlinRectModes } from './mindlin.js';
import { findMode, median, pct, sweepResonances } from './testkit.js';

const NU = 0.3;
const elastic = modes => modes.filter(mode => mode.L > 0).map(mode => mode.L);

// the existing per-mode correction, as Λ
const corrected = (L, hOverA, nu = NU) => {
  const mat = { E: 1e9, rho: 1000, nu, h: hOverA };
  const c = Math.sqrt(mat.E * hOverA ** 3 / (12 * (1 - nu * nu)) / (mat.rho * hOverA));
  return mindlinOmega(L * c, mat) / c;
};

test('simply supported rectangle: the closed form, itself checked against Hashemi & Arsanjani 2005 Table 9', () => {
  for (const [beta, hOverA] of [[1, 0.1], [1, 0.2], [2, 0.1], [0.5, 0.2]]) {
    const thin = rectModes('simply', NU, beta).modes.map(mode => corrected(mode.L, hOverA));
    const modes = mindlinRectModes('simply', NU, beta, hOverA, { nb: 18 }).modes;
    // up to Λ = 150: from the cut-off (Λ = 161 at h/a = 0.2) the thickness-shear family, absent from the closed form's flexural root, interleaves
    thin.filter(L => L < 150).slice(0, 20).forEach((L, k) => assert.ok(pct(modes[k].L, L) < 1e-3, `β ${beta} h/a ${hOverA} #${k + 1}: ${modes[k].L} vs ${L}`));
  }
  // Table 9, κ² = π²/12, h/a = 0.1, as read for physics.test.js
  const square = mindlinRectModes('simply', NU, 1, 0.1, { nb: 18 }).modes;
  [19.0584, 45.4478, 45.4478, 69.7166, 84.9263, 84.9263, 106.5154, 106.5154, 133.3705]
    .forEach((want, k) => assert.ok(pct(square[k].L, want) < 3e-4, `Table 9 #${k + 1}`));
});

test('Lévy rectangles: exact state-space solution, scipy, y edges simply supported', () => {
  // Euler–Lagrange ODEs in x with numerically computed exponents, then the 6 × 6 edge determinant (levy.py); ν = 0.3, κ² = π²/12
  const EDGE = { S: 'simply', C: 'clamped', F: 'free' };
  const exact = {
    'CC|0.1|1': [26.6449, 49.0625, 59.1183, 78.6831, 86.7203, 101.1524, 111.8482, 118.6575],
    'FF|0.1|1': [9.4388, 15.3838, 33.8408, 36.3339, 42.7604, 62.084, 66.1367, 76.5405],
    'CF|0.1|1': [12.2454, 30.3862, 38.6074, 55.7431, 62.6487, 78.4292, 85.7649, 93.5436],
    'SC|0.1|1': [22.3759, 47.0628, 52.0899, 74.0041, 85.759, 93.0643, 109.072, 112.5271],
    'SF|0.1|1': [11.3704, 26.142, 38.2605, 53.1988, 55.576, 78.27, 80.9756, 92.0781],
    'CC|0.2|1': [22.3085, 39.7558, 44.4669, 58.3574, 65.5675, 70.2846, 80.006, 81.8396],
    'FF|0.2|1': [8.9774, 14.0803, 29.0958, 31.2149, 35.8905, 49.4472, 52.3701, 59.4499],
    'CF|0.2|1': [11.3514, 25.706, 32.8323, 44.6161, 48.3715, 60.6236, 63.8489, 69.8578],
    'SF|0.2|1': [10.6901, 23.1226, 32.6558, 43.4778, 45.21, 60.5664, 62.1585, 69.3679],
    'CF|0.1|2': [5.5568, 12.2454, 23.3255, 23.4442, 30.3862, 38.6074, 41.324, 55.7431],
    'FF|0.2|2': [2.331, 6.2277, 8.9774, 14.0803, 18.9645, 22.3183, 24.0684, 29.0958],
  };
  for (const [key, want] of Object.entries(exact)) {
    const [ends, hOverA, beta] = key.split('|');
    const edges = [EDGE[ends[0]], EDGE[ends[1]], 'simply', 'simply'];
    const got = elastic(mindlinRectModes(edges, NU, +beta, +hOverA, { nb: 18 }).modes);
    want.forEach((value, k) => assert.ok(pct(got[k], value) < 5e-3, `${key} #${k + 1}: ${got[k]} vs ${value}`));
  }
});

test('thick squares against 3D elasticity (Ritz, scipy): old per-mode correction against true Mindlin', () => {
  // ritz3d.py, Legendre products, flexural class, ν = 0.3, first three modes of each parity class; bases grown until ≤ 0.03 % moved
  const threeD = {
    'free|0.1': { '00': [18.9544, 23.3449, 55.8192], '01': [31.9542, 55.4893, 87.9509], '11': [12.723, 60.7587, 67.8736] },
    'free|0.2': { '00': [17.4332, 21.2519, 45.3072], '01': [27.6468, 45.5363, 66.7733], '11': [11.7098, 48.4536, 54.1176] },
    'clamped|0.1': { '00': [32.7386, 103.4988, 104.4793], '01': [62.5535, 125.2745, 152.8051], '11': [87.772, 170.3394, 171.886] },
    'clamped|0.2': { '00': [26.8839, 72.25, 73.2401], '01': [47.0715, 85.8075, 101.0824], '11': [63.3124, 111.4627, 112.3392] },
  };
  const bounds = { // [new, old] worst error in % over these nine modes
    'free|0.1': [0.5, 3.4], 'free|0.2': [1, 5.3], 'clamped|0.1': [1.6, 5.3], 'clamped|0.2': [2.8, 10.6],
  };
  for (const [key, classes] of Object.entries(threeD)) {
    const [edge, tau] = key.split('|');
    const fresh = mindlinRectModes(edge, NU, 1, +tau).modes;
    const old = rectModes(edge, NU, 1).modes;
    const inClass = (modes, cls) => modes.filter(mode => mode.L > 1e-3 && `${mode.px}${mode.py}` === cls).map(mode => mode.L);
    let worstNew = 0, worstOld = 0;
    for (const [cls, values] of Object.entries(classes)) {
      const n = inClass(fresh, cls), o = inClass(old, cls).map(L => corrected(L, +tau));
      values.forEach((value, k) => {
        worstNew = Math.max(worstNew, pct(n[k], value));
        worstOld = Math.max(worstOld, pct(o[k], value));
      });
    }
    assert.ok(worstNew < bounds[key][0], `${key}: Mindlin within ${worstNew} %`);
    assert.ok(worstOld > 2 * worstNew && worstOld < bounds[key][1], `${key}: old ${worstOld} %`);
  }
});

test('thin limit: rectangles recover rectModes without shear locking', () => {
  for (const [edge, hOverA, tol] of [['clamped', 0.005, 0.06], ['free', 0.005, 0.2], ['simply', 0.005, 1e-3], ['clamped', 0.002, 0.01], ['free', 0.002, 0.04]]) {
    const thin = elastic(rectModes(edge, NU, 1).modes.map(mode => ({ L: mode.L < 1e-3 ? 0 : mode.L })));
    const got = elastic(mindlinRectModes(edge, NU, 1, hOverA, { nb: 18 }).modes);
    // locking would drive Ritz far above Kirchhoff; only the physical drop (and the free edge's O(h) layer) may remain
    for (let k = 0; k < 20; k++) assert.ok(pct(got[k], corrected(thin[k], hOverA)) < tol, `${edge} ${hOverA} #${k + 1}: ${got[k]} vs ${thin[k]}`);
  }
  // free-square fundamental at h/a = 0.01 sits on 3D (13.4172), not on Kirchhoff + dispersion (13.4648): the edge layer is real, O(h)
  const fundamental = elastic(mindlinRectModes('free', NU, 1, 0.01).modes)[0];
  assert.ok(pct(fundamental, 13.4172) < 0.03, `${fundamental}`);
  assert.ok(pct(corrected(elastic(rectModes('free', NU, 1).modes.map(mode => ({ L: mode.L < 1e-3 ? 0 : mode.L })))[0], 0.01), 13.4172) > 0.3);
});

test('rectangle convergence with nb: exact Lévy values reached from above', () => {
  const edges = ['free', 'free', 'simply', 'simply'];
  const exact = [9.4388, 15.3838, 33.8408, 36.3339, 42.7604, 62.084, 66.1367, 76.5405];
  const worst = nb => {
    const got = elastic(mindlinRectModes(edges, NU, 1, 0.1, { nb }).modes);
    return Math.max(...exact.map((value, k) => pct(got[k], value)));
  };
  const [coarse, medium, fine] = [8, 12, 18].map(worst);
  assert.ok(coarse < 0.25 && medium < 5e-3 && fine < 5e-4, `${coarse}, ${medium}, ${fine} %`);
  assert.ok(medium < coarse / 30 && fine < medium / 3);
});

test('modal sum: uniform-load static deflection matches closed forms, only with rotary inertia in the normalisation', () => {
  // Square, simply supported: Navier series with shear, w̄ = Σ 16/(π²ij)·(±1)/k⁴·(1 + k²/s) per unit qa³/D
  const tau = 0.2, s = 6 * (Math.PI ** 2 / 12) * (1 - NU) / tau ** 2;
  let navier = 0;
  for (let i = 1; i < 400; i += 2) {
    for (let j = 1; j < 400; j += 2) {
      const k2 = Math.PI ** 2 * (i * i + j * j);
      navier += 16 / (Math.PI ** 2 * i * j) * (-1) ** ((i + j) / 2 - 1) / k2 ** 2 * (1 + k2 / s);
    }
  }
  const { modes, nbx, nby, atX, atY } = mindlinRectModes('simply', NU, 1, tau);
  const [nodes, weights] = gaussLegendre(nbx + 12);
  const meanX = new Float64Array(nbx), meanY = new Float64Array(nby);
  nodes.forEach((x, q) => {
    atX(x).forEach((value, i) => { meanX[i] += weights[q] * value; });
    atY(x).forEach((value, j) => { meanY[j] += weights[q] * value; });
  });
  const X = atX(0.5), Y = atY(0.5);
  let sum = 0;
  for (const mode of modes) {
    let centre = 0, mean = 0;
    for (let i = 0; i < nbx; i++) {
      for (let j = 0; j < nby; j++) {
        centre += mode.vec[i * nby + j] * X[i] * Y[j];
        mean += mode.vec[i * nby + j] * meanX[i] * meanY[j];
      }
    }
    sum += centre * mean / mode.L ** 2;
  }
  assert.ok(pct(sum, navier) < 0.3, `square: ${sum} vs ${navier}`);

  // Clamped circle: bending as Kirchhoff, plus the shear part from S(w′ + ψ) = −qr/2: w̄(0) = 1/64 + 1/(4s)
  const circles = mindlinCircleModes('clamped', NU, tau, 22).filter(mode => mode.n === 0);
  let circleSum = 0;
  for (const mode of circles) {
    const [nodesR, weightsR] = gaussLegendre(40 + Math.ceil(Math.sqrt(mode.wave[0])));
    let mean = 0;
    nodesR.forEach((rho, q) => { mean += weightsR[q] * mindlinCircleFields(mode, rho)[0] * rho; });
    circleSum += mode.R[0] * mode.inv * 2 * Math.PI * mean / mode.L ** 2;
  }
  assert.ok(pct(circleSum, 1 / 64 + 1 / (4 * s)) < 0.2, `circle: ${circleSum}`);
});

test('circles: exact solution matches a polar Ritz with no Bessel functions, below and above the thickness-shear cut-off', () => {
  // circle_ritz.py: Ritz on Mindlin's energy with polynomial fields in scipy, ν = 0.3, κ² = π²/12; cut-off at h/a = 0.2 is Λ = 160.97
  const ritz = {
    'free|0.1|0': [8.868, 36.041, 76.676, 126.274], 'free|0.1|1': [19.711, 54.257, 99.935, 152.75], 'free|0.1|5': [31.27, 82.722, 137.951, 197.055],
    'free|0.2|0': [8.505, 31.111, 59.645, 90.059, 120.57, 149.63, 171.181, 183.356, 199.042], 'free|0.2|1': [17.978, 44.434, 74.331, 105.033, 135.12, 160.67, 163.18, 167.273, 175.087],
    'clamped|0.1|0': [9.941, 36.479, 75.664, 123.32], 'clamped|0.1|1': [20.177, 53.847, 97.877, 148.678], 'clamped|0.1|5': [76.588, 131.757, 189.686, 249.584],
    'clamped|0.2|0': [9.24, 30.211, 56.682, 85.571], 'clamped|0.2|1': [17.758, 42.374, 70.459, 100.11, 130.43, 160.71, 165.367, 170.243],
    'simply|0.1|0': [4.894, 28.24, 65.942, 113.574], 'simply|0.1|1': [13.562, 44.731, 88.027, 139.297], 'simply|0.1|5': [67.572, 122.605, 181.41, 242.556],
    'simply|0.2|0': [4.777, 24.994, 52.514, 82.766], 'simply|0.2|1': [12.708, 37.591, 66.984, 97.902, 129.2, 160.41, 161.48, 169.924, 176.02, 183.177],
  };
  const solved = new Map();
  for (const [key, want] of Object.entries(ritz)) {
    const [bc, tau, n] = key.split('|');
    if (!solved.has(bc + tau)) solved.set(bc + tau, mindlinCircleModes(bc, NU, +tau, 16.5));
    const got = elastic(solved.get(bc + tau).filter(mode => mode.n === +n));
    want.forEach((value, k) => assert.ok(Math.abs(got[k] - value) < 0.012, `${key} #${k + 1}: ${got[k]} vs ${value}`));
  }
});

test('circles: published Mindlin tables (Senjanović et al., Arch. Mech. 66, 2014, Tables 1–4; n = 0)', () => {
  // Their columns for n ≥ 1 use a "modified" theory: at h/R = 0.2 on free edges they sit up to 5 % above ours, where the
  // polar Ritz and 3D elasticity side with ours, so only n = 0 (where the rotational field decouples) is a fair test.
  const table = { // Tables 1–3, κ² = 5/6, ν = 0.3 inferred from their thin values: 5.358 and 10.216 are ν = 0.3's
    simply: { 0.1: [4.894, 28.253, 66.011, 113.761, 167.902], 0.2: [4.778, 25.031, 52.649, 83.050, 114.340] },
    clamped: { 0.1: [9.944, 36.511, 75.780, 123.581, 176.880], 0.2: [9.250, 30.283, 56.870, 85.906, 116.059] },
    free: { 0.1: [8.869, 36.059, 76.758, 126.482, 181.862], 0.2: [8.508, 31.156, 59.794, 90.360, 121.047] },
  };
  for (const [bc, rows] of Object.entries(table)) {
    for (const [tau, want] of Object.entries(rows)) {
      const got = elastic(mindlinCircleModes(bc, NU, +tau, 14, { kappa2: 5 / 6 }).filter(mode => mode.n === 0));
      want.forEach((value, k) => assert.ok(Math.abs(got[k] - value) < 0.002, `${bc} ${tau} #${k + 1}: ${got[k]} vs ${value}`));
    }
  }
  // Table 4, κ² = π²/12, h/R = 0.25
  for (const [bc, want] of [['simply', [4.696, 23.254, 46.775, 71.603]], ['clamped', [8.807, 27.253, 49.420, 73.054]], ['free', [8.267, 28.605, 52.584, 76.936]]]) {
    const got = elastic(mindlinCircleModes(bc, NU, 0.25, 9).filter(mode => mode.n === 0));
    want.forEach((value, k) => assert.ok(Math.abs(got[k] - value) < 0.002, `${bc} h/R 0.25 #${k + 1}: ${got[k]}`));
  }
});

test('thin limit: circles recover circleModes in frequency, shape, nodal radii and modal mass', () => {
  for (const bc of ['free', 'clamped', 'simply']) {
    const thin = circleModes(bc, NU, 12);
    for (const hOverA of [0.001, 0.005]) {
      const modes = mindlinCircleModes(bc, NU, hOverA, 12);
      for (const mode of thin.filter(m => m.L > 0)) {
        const got = findMode(modes, mode.n, mode.m);
        // the gap is the physical drop, O(h²) from dispersion plus O(h) from a free edge's layer
        assert.ok(pct(got.L, mode.L) < 40 * hOverA, `${bc} (${mode.n},${mode.m}) ${got.L} vs ${mode.L}`);
        assert.ok(pct(got.inv, mode.inv) < 1, `${bc} inv ${got.inv} vs ${mode.inv}`);
        mindlinNodalRadii(got).forEach((r, k) => assert.ok(Math.abs(r - nodalRadii(mode)[k]) < 3e-4, `${bc} radius`));
        for (let j = 0; j < 257; j += 16) assert.ok(Math.abs(got.R[j] - mode.R[j]) < 0.01 * mode.Rmax, `${bc} shape at ${j}`);
      }
    }
  }
  // Rocking keeps a free disc's true moment of inertia: m(a²/4 + h²/12)
  const rocking = mindlinCircleModes('free', NU, 0.2, 3).find(mode => mode.n === 1 && !mode.L);
  assert.ok(pct(1 / rocking.inv, Math.PI * (0.25 + 0.04 / 12)) < 1e-12);
});

test('circular modes are orthogonal under the full Mindlin mass, and not under w alone', () => {
  const tau = 0.2, alpha = tau * tau / 12;
  const modes = mindlinCircleModes('free', NU, tau, 12).filter(mode => mode.n === 3);
  const [nodes, weights] = gaussLegendre(120);
  const inner = (a, b, rotary) => nodes.reduce((total, rho, q) => {
    const [wa, ra, ta] = mindlinCircleFields(a, rho), [wb, rb, tb] = mindlinCircleFields(b, rho);
    return total + weights[q] * rho * (wa * wb + rotary * (ra * rb + ta * tb));
  }, 0);
  const norm = mode => inner(mode, mode, alpha);
  const overlap = (rotary) => inner(modes[0], modes[1], rotary) / Math.sqrt(norm(modes[0]) * norm(modes[1]));
  assert.ok(Math.abs(overlap(alpha)) < 1e-6, `with rotary inertia ${overlap(alpha)}`);
  assert.ok(Math.abs(overlap(0)) > 0.01, `w alone ${overlap(0)}`);
  assert.ok(pct(modes[0].inv, 1 / (Math.PI * norm(modes[0]))) < 1e-4);
});

test('circles against 3D elasticity (Ritz, scipy): Mindlin and the free edge\'s O(h) layer are real', () => {
  // disc3d.py, ν = 0.3, κ² = π²/12; two lowest modes at n = 0 and n = 5; bases grown until ≤ 0.03 % moved
  const threeD = {
    'free|0.1|0': [8.872, 36.132], 'free|0.1|5': [31.327, 83.184], 'free|0.2|0': [8.519, 31.344], 'free|0.2|5': [27.391, 64.095],
    'clamped|0.1|0': [9.973, 36.683], 'clamped|0.1|5': [77.245, 133.333], 'clamped|0.2|0': [9.312, 30.620], 'clamped|0.2|5': [58.389, 92.716],
  };
  for (const [key, want] of Object.entries(threeD)) {
    const [bc, tau, n] = key.split('|');
    const got = elastic(mindlinCircleModes(bc, NU, +tau, 12).filter(mode => mode.n === +n));
    want.forEach((value, k) => assert.ok(pct(got[k], value) < 13 * +tau, `${key} #${k + 1}: ${got[k]} vs ${value}`));
  }
  // (9,0), free, h/a = 0.005: 3D 101.9516; thin theory with per-mode dispersion misses the edge layer by ten times the Mindlin gap
  const high = findMode(mindlinCircleModes('free', NU, 0.005, 11), 9, 0).L;
  const old = corrected(findMode(circleModes('free', NU, 11), 9, 0).L, 0.005);
  assert.ok(pct(high, 101.9516) < 0.01 && pct(old, 101.9516) > 0.1, `${high}, old ${old}`);
});

test('Waller\'s outer circles: thickness cannot explain them, at her thickness or any plausible one', () => {
  // Leissa Table 2.8 footnote: t = 0.0788 in, plates 12.0 and 10.1 in across, rough readings on one of 18.36 in (t assumed the same).
  // Table 2.9, n = 0: s = 4 outer .918, s = 5 outer .956, s = 6 outer .958. ν = 0.33.
  const measured = [[4, 0.918], [5, 0.956], [6, 0.958]];
  const thin = circleModes('free', 0.33, 20);
  for (const hOverA of [0.0788 / 9.18, 0.0788 / 6, 0.0788 / 5.05]) {
    const modes = mindlinCircleModes('free', 0.33, hOverA, 20);
    for (const [s] of measured) {
      const shift = mindlinNodalRadii(findMode(modes, 0, s)).at(-1) - nodalRadii(findMode(thin, 0, s)).at(-1);
      assert.ok(Math.abs(shift) < 2e-4, `h/a ${hOverA}, s = ${s}: Mindlin moves the circle by ${shift}`);
    }
  }
  // 3D elasticity (disc3d.py, scipy) at h/a = 0.01313: mid-plane and top surface sit on thin theory too, four digits
  const threeD = { 4: 0.9194, 5: 0.9352, 6: 0.9458 };
  for (const [s, r] of Object.entries(threeD)) assert.ok(Math.abs(nodalRadii(findMode(thin, 0, +s)).at(-1) - r) < 2e-4);

  // 3D at h/a = 0.2 moves the (0,5), (0,6) mid-plane circles by +0.012, +0.028: the size needed, but Λ falls by half and inner circles move too
  // What does fit them: Colwell & Hardy's approximation, eq. 2.19, ρ = j/λ with j a zero of J_n, which Leissa's own Table 2.7 prints for s = 5
  const zeros = { 4: 11.7915344, 5: 14.9309177, 6: 18.0710640 }; // the 4th–6th zeros of J₀, checked just below
  Object.values(zeros).forEach(z => assert.ok(Math.abs(bessel(0, z)[0]) < 1e-6));
  const outerByApproximation = s => zeros[s] / Math.sqrt(findMode(thin, 0, s).L);
  assert.ok(Math.abs(outerByApproximation(5) - 0.956) < 0.004 && Math.abs(outerByApproximation(6) - 0.958) < 0.003);
  assert.ok(Math.abs(nodalRadii(findMode(thin, 0, 5)).at(-1) - 0.956) > 0.02 && Math.abs(nodalRadii(findMode(thin, 0, 6)).at(-1) - 0.958) > 0.012);
  assert.ok(Math.abs(outerByApproximation(4) - 0.918) > 0.02 && Math.abs(nodalRadii(findMode(thin, 0, 4)).at(-1) - 0.918) < 0.002); // s = 4 follows theory
});

test('Plate runs on these modes unchanged: centre-driven circle still matches Escaler & De La Torre 2018', () => {
  const cfg = { shape: 'circle', bc: 'free', mount: 'shaker', ma: 0.08, E: 69e9, rho: 2700, nu: 0.35, h: 0.0008, d: 0.24, Lx: 0.24, Ly: 0.24, zeta: 0.005, fmax: 5000, G: 25, mindlin: false };
  const plate = new Plate();
  const lamMax = circleLamFor(cfg);
  plate.configure({ ...cfg, solved: { id: 'mindlin', thick: true, lamMax, modes: mindlinCircleModes('free', 0.35, 0.0008 / 0.12, lamMax) } });
  plate.setDrive(0, 0);
  const measured = [70, 305, 844, 1657, 2747, 4129];
  const resonances = sweepResonances(plate, 40, 5000);
  measured.forEach((f, i) => assert.ok(pct(resonances[i], f) < 2, `${resonances[i]} Hz vs ${f}`));
  assert.ok(median(measured.map((f, i) => pct(resonances[i], f))) < 0.6);

  // a rectangle needs nothing else either: the top of the kept range is a real frequency and the field peaks at the first resonance
  const square = new Plate();
  const rectCfg = { shape: 'rect', bc: 'clamped', E: 70e9, rho: 2700, nu: NU, h: 0.02, d: 0.2, Lx: 0.2, Ly: 0.2, zeta: 0.005, fmax: 5000, G: 25, mindlin: false };
  square.configure({ ...rectCfg, solved: { id: 'mindlin', thick: true, ritz: mindlinRectModes('clamped', NU, 1, 0.1) } });
  square.setDrive(0, 0);
  const first = square.modes[0].f;
  // Λ₁ = 32.49 for clamped h/a = 0.1 in true Mindlin, below 3D's 32.74
  assert.ok(pct(square.modes[0].L, 32.74) < 1);
  assert.ok(square.field([{ f: first, F: 1 }], 9.81) > 5 * square.field([{ f: first * 1.2, F: 1 }], 9.81));
});
