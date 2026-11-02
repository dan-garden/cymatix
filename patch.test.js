// Finite rigid patch against exact theory and measurements. Every number marked "scipy" comes from an independent solver
// (J, Y, I, K characteristic determinants, or a Navier series with no tail), not from this code.
import test from 'node:test';
import assert from 'node:assert/strict';
import { Plate } from './physics.js';
import { RigidPatch, tailSeries } from './patch.js';

const TAU = 2 * Math.PI, g = 9.81;
const percent = (got, want) => 100 * (got / want - 1);
const rms = values => Math.sqrt(values.reduce((sum, x) => sum + x * x, 0) / values.length);
const nearest = (values, target) => values.reduce((best, x) => (Math.abs(x - target) < Math.abs(best - target) ? x : best));
const A = 0.12; // plate radius for the theory tests

function circle(bc, nu, lamMax, extra = {}) {
  const base = { E: 69e9, rho: 2700, h: 0.001 };
  const wScale = Math.sqrt(base.E * base.h ** 2 / (12 * (1 - nu * nu) * base.rho)) / (A * A);
  const cfg = { ...base, nu, shape: 'circle', bc, d: 2 * A, Lx: 2 * A, Ly: 2 * A, zeta: 0, G: 9, mindlin: false, fmax: 0.999 * lamMax ** 2 * wScale / (1.3 * TAU), ...extra };
  const plate = new Plate();
  plate.configure(cfg);
  plate.setDrive(0, 0);
  return { plate, cfg, wScale };
}

function determinant(matrix, n) {
  const m = Float64Array.from(matrix);
  let product = 1;
  for (let c = 0; c < n; c++) {
    let pivot = c;
    for (let r = c; r < n; r++) if (Math.abs(m[r * n + c]) > Math.abs(m[pivot * n + c])) pivot = r;
    if (pivot !== c) {
      product = -product;
      for (let k = 0; k < n; k++) [m[c * n + k], m[pivot * n + k]] = [m[pivot * n + k], m[c * n + k]];
    }
    product *= m[c * n + c];
    for (let r = c + 1; r < n; r++) {
      const factor = m[r * n + c] / m[c * n + c];
      for (let k = c; k < n; k++) m[r * n + k] -= factor * m[c * n + k];
    }
  }
  return product;
}

const bisect = (fn, lo, hi, flo) => {
  for (let i = 0; i < 60; i++) {
    const mid = (lo + hi) / 2, fm = fn(mid);
    if (flo * fm <= 0) hi = mid;
    else {
      lo = mid;
      flo = fm;
    }
  }
  return (lo + hi) / 2;
};

// Zeros of fn between consecutive poles, with samples crowded towards the poles because a small patch moves a zero very little.
function zerosBetween(fn, poles, top) {
  const edges = [1e-3, ...poles.filter(L => L > 0 && L < top), top];
  const found = [];
  for (let e = 0; e + 1 < edges.length; e++) {
    let a = 0, fa = NaN;
    for (let i = 1; i < 60; i++) {
      const t = i / 60;
      const b = edges[e] + (edges[e + 1] - edges[e]) * (0.5 - 0.5 * Math.cos(Math.PI * t * t * (3 - 2 * t)));
      const fb = fn(b);
      if (fa * fb < 0) found.push(bisect(fn, a, b, fa));
      a = b;
      fa = fb;
    }
  }
  return found;
}

/** λ² of the plate clamped along the disk, harmonic n: where the receptance block between the ring forces is singular. */
function clampedRoots(plate, cfg, wScale, r0, n, top) {
  const patch = new RigidPatch(plate, cfg, { r0, mount: 'rigid', harmonics: [[n, 0]] });
  const poles = plate.modes.filter(m => m.n === n).map(m => m.L).sort((a, b) => a - b);
  return zerosBetween(L => determinant(patch.receptance(L * wScale)[0], patch.U), poles, top);
}

test('tail series equals direct quadrature of the infinite-plate integral (scipy)', () => {
  for (const [x, omega, want] of [[0.01, 0.1, 0.94877109], [0.01, 0.7, 1.18111672], [0.25, 0.3, 0.44786338], [1, 0.5, -0.07757378], [1, 0.9, -0.02047515],
    [4, 0.2, -0.08334904], [4, 0.77, -0.16981837], [16, 0.5, -0.03483492], [49, 0.4, -0.01213143]]) {
    assert.ok(Math.abs(tailSeries(x, omega) - want) < 2e-8, `τ(${x}, ${omega}) = ${tailSeries(x, omega)}, scipy ${want}`);
  }
  // r → 0 is exactly the point tail of Plate's #mountFactor: ln((1+Ω)/(1−Ω))/(2Ω)
  for (const omega of [0.1, 0.5, 0.9]) assert.ok(Math.abs(tailSeries(0, omega) - Math.log((1 + omega) / (1 - omega)) / (2 * omega)) < 1e-14);
});

test('modal sum plus tail reproduces the static Green function of a clamped disk, which the modal sum alone cannot', () => {
  // closed form, a = 1: (a²/16πD)[d²·ln(d²/q) + (1−r²)(1−s²)], d² = r²+s²−2rs·cos φ, q = 1+r²s²−2rs·cos φ; at r = s = 0 it is Fa²/16πD (Timoshenko)
  const { plate, cfg } = circle('clamped', 0.3, 30);
  const D = cfg.E * cfg.h ** 3 / (12 * (1 - cfg.nu ** 2));
  const exact = (r, s, phi) => {
    const d2 = r * r + s * s - 2 * r * s * Math.cos(phi), q = 1 + r * r * s * s - 2 * r * s * Math.cos(phi);
    return A * A / (16 * Math.PI * D) * ((d2 > 0 ? d2 * Math.log(d2 / q) : 0) + (1 - r * r) * (1 - s * s));
  };
  const beta = 0.02;
  const check = (tail, maxPercent, minPercent) => {
    const patch = new RigidPatch(plate, cfg, { r0: beta * A, mount: 'rigid', tail, harmonics: [[0, 0], [1, 0], [2, 0]] });
    const [G] = patch.receptance(TAU * 0.5);
    const errors = [];
    patch.basis.forEach((bi, i) => patch.basis.forEach((bj, j) => {
      if (bi.h !== bj.h) return;
      let want = 0;
      for (let l = 0; l < 4096; l++) {
        const phi = TAU * l / 4096;
        want += (bi.h ? 2 : 1) * Math.cos(bi.h * phi) * exact(bi.radius * beta, bj.radius * beta, phi) / 4096;
      }
      errors.push(Math.abs(percent(G[i * patch.U + j], want)));
    }));
    const worst = Math.max(...errors);
    assert.ok(worst < maxPercent && worst > minPercent, `tail ${tail}: worst ${worst.toFixed(2)} %`);
  };
  check(true, 1.5, 0); // 1.3 % at λ_K = 30
  check(false, 100, 50); // 94 % without it: the ring-to-ring differences live almost entirely in the tail
});

test('r₀ → 0 recovers the point model', () => {
  const material = { E: 69e9, rho: 2700, nu: 0.33, h: 0.001, d: 0.24, Lx: 0.24, Ly: 0.24, zeta: 0.005, fmax: 5000, G: 25 };
  const worstDifference = (cfg, u, v, options) => {
    const plate = new Plate();
    plate.configure(cfg);
    plate.setDrive(u, v);
    const patch = new RigidPatch(plate, cfg, options);
    const cr = new Float64Array(plate.modes.length), ci = new Float64Array(plate.modes.length);
    let worst = 0, biggest = 0;
    for (const f of [97, 500, 1234.5, 3300]) {
      plate.field([{ f, F: 1 }], g);
      patch.coefficients(TAU * f, 1, g, cr, ci);
      for (let k = 0; k < cr.length; k++) {
        worst = Math.max(worst, Math.hypot(cr[k] - plate.cr[k], ci[k] - plate.ci[k]));
        biggest = Math.max(biggest, Math.hypot(plate.cr[k], plate.ci[k]));
      }
    }
    return worst / biggest;
  };
  for (const [mount, ma] of [['rigid', 0], ['shaker', 0.02], ['force', 0]]) {
    const options = { r0: 0, mount: mount === 'rigid' ? 'rigid' : 'shaker', ma };
    const rect = { ...material, shape: 'rect', bc: 'free', mount, ma };
    const round = { ...material, shape: 'circle', bc: 'free', mount, ma };
    assert.ok(worstDifference(rect, 0.3, 0.1, options) < 1e-12, `${mount} rectangle at r₀ = 0`);
    assert.ok(worstDifference(round, 0, 0, options) < 1e-6, `${mount} circle at r₀ = 0, centred`); // R[] is stored as Float32
    assert.ok(worstDifference(round, 0.3, 0.1, options) < 5e-4, `${mount} circle at r₀ = 0`); // Plate interpolates d0 linearly in R[]; the patch uses the exact shape
  }
  // small but finite: centred any mount; off centre with the tilt free, since a locked tilt is a different plate (next test)
  const small = { r0: 1e-5, mount: 'shaker', ma: 0.02 };
  assert.ok(worstDifference({ ...material, shape: 'rect', bc: 'free', mount: 'shaker', ma: 0.02 }, 0.3, 0.1, small) < 1e-4);
  const bolt = r0 => worstDifference({ ...material, shape: 'circle', bc: 'free', mount: 'rigid', ma: 0 }, 0, 0, { r0, mount: 'rigid' });
  assert.ok(bolt(1e-5) < 1e-3 && bolt(1e-4) < 0.05 && bolt(1e-4) > 10 * bolt(1e-5), 'centred bolt converges on the point as r₀ falls');
  assert.ok(worstDifference({ ...material, shape: 'rect', bc: 'free', mount: 'rigid', ma: 0 }, 0.3, 0.1, { r0: 1e-5, mount: 'rigid', lockTilt: false }) < 1e-3);
});

test('exact annular plates, clamped inside (scipy): all three outer edges, r₀/a from 0.02 to 0.3', () => {
  const exact = {
    free: { 0.02: [[3.7941, 21.4003, 62.3236], [2.347, 23.3756, 65.9402], [5.3688, 35.3283, 84.6142]], 0.05: [[3.9234, 22.6422, 66.1765], [2.8309, 24.9721, 70.0606], [5.424, 35.6846, 85.8865]], 0.1: [[4.2374, 25.2618, 73.9006], [3.4781, 27.6728, 77.4868], [5.6227, 36.9409, 90.1835]], 0.3: [[6.6604, 42.6142, 123.4661], [6.5523, 44.6311, 125.8234], [7.9565, 50.9473, 133.0848]] },
    clamped: { 0.02: [[23.2503, 63.6568, 124.8619], [24.3723, 66.985, 130.1255], [34.9476, 84.8327, 154.4173]], 0.05: [[24.5465, 67.5581, 132.7759], [26.066, 71.1402, 137.9423], [35.3163, 86.1158, 157.4343]], 0.1: [[27.2805, 75.3663, 148.2137], [28.9158, 78.6352, 152.5418], [36.6173, 90.4483, 167.1351]], 0.3: [[45.3462, 125.3621, 246.1573], [46.6435, 127.381, 248.503], [51.1388, 133.6743, 255.669]] },
    simply: { 0.02: [[15.1368, 50.798, 106.8707], [16.2386, 53.7694, 111.7626], [25.6576, 70.3069, 134.7909]], 0.05: [[15.9777, 53.8899, 113.6277], [17.432, 57.1772, 118.5317], [25.8899, 71.2857, 137.2776]], 0.1: [[17.7894, 60.1441, 126.8805], [19.3942, 63.2432, 131.0855], [26.717, 74.631, 145.3746]], 0.3: [[29.9777, 100.4228, 211.1294], [31.4025, 102.4828, 213.5194], [36.2434, 108.9216, 220.826]] },
  };
  for (const [bc, table] of Object.entries(exact)) {
    const { plate, cfg, wScale } = circle(bc, 0.3, 32);
    for (const [beta, rows] of Object.entries(table)) {
      rows.forEach((want, n) => {
        const got = clampedRoots(plate, cfg, wScale, beta * A, n, want.at(-1) * 1.05);
        want.forEach((L, s) => assert.ok(Math.abs(percent(got[s], L)) < (beta < 0.2 ? 0.35 : 1), `${bc} b/a ${beta} n ${n} s ${s}: ${got[s]} vs ${L}`));
      });
    }
  }
});

test('published annular tables, Leissa 1969 Tables 2.30, 2.24 and 2.31 (ν = 0.3)', () => {
  // Read from the book's text. Rows [n, s, λ² at b/a = 0.1, 0.3]. The patch must land within the table's rounding (up to 0.5 %) plus its own error.
  const free = [[0, 0, 4.23, 6.66], [2, 0, 5.62, 7.96], [3, 0, 12.4, 13.27], [0, 1, 25.3, 42.6], [2, 1, 37.0, 50.9], [3, 1, 53.2, 62.1], [1, 1, null, 44.6], [1, 0, null, null]];
  const simply = [[0, 0, 17.8, 29.9], [2, 0, 26.8, 36.2], [3, 0, 40.0, 45.4], [0, 1, 60.1, 100], [2, 1, 74.7, 109], [3, 1, 95.3, 120], [1, 0, null, 31.4], [1, 1, null, 102]];
  const check = (bc, rows) => {
    const { plate, cfg, wScale } = circle(bc, 0.3, 32);
    for (const [n, s, at1, at3] of rows) {
      [[0.1, at1], [0.3, at3]].forEach(([beta, want]) => {
        if (want === null) return;
        const got = clampedRoots(plate, cfg, wScale, beta * A, n, want * 1.1);
        assert.ok(Math.abs(percent(got[s], want)) < 1.3, `${bc} b/a ${beta} (${n},${s}): ${got[s]} vs ${want}`);
      });
    }
  };
  check('free', free);
  check('simply', simply);
  // Table 2.31 (Southwell), b/a ≤ 0.4: [n, b/a, λ²]
  const { plate, cfg, wScale } = circle('free', 0.3, 32);
  for (const [n, beta, want] of [[0, 0.276, 6.25], [1, 0.397, 9.0], [2, 0.186, 6.25], [2, 0.349, 9.0]]) {
    const got = clampedRoots(plate, cfg, wScale, beta * A, n, want * 1.2);
    assert.ok(Math.abs(percent(nearest(got, want), want)) < 1, `Southwell n ${n} b/a ${beta}: ${nearest(got, want)} vs ${want}`);
  }
  // Table cells that disagree with exact theory, where the patch follows scipy (J, Y, I, K) and Raju's Table 2.29 (3.482 at ν = 1/3),
  // not the printed value: free n = 1 b/a 0.1 prints 3.14, b/a 0.3 prints 6.33, (1,1) at 0.1 prints 27.3; simply (1,0) at 0.1 prints 19.0.
  assert.ok(Math.abs(percent(clampedRoots(plate, cfg, wScale, 0.1 * A, 1, 5)[0], 3.4781)) < 0.5);
  assert.ok(Math.abs(percent(clampedRoots(plate, cfg, wScale, 0.3 * A, 1, 9)[0], 6.5523)) < 0.6);
});

test('a clamp of any size still holds slope: n = 1 does not return to the free plate as r₀ → 0', () => {
  // scipy: the lowest n = 1 mode with the hub clamped is λ² = 2.347 at b/a = 0.02 and 3.478 at 0.1. A point pin leaves that mode at 0.
  const { plate, cfg, wScale } = circle('free', 0.3, 32);
  assert.ok(Math.abs(percent(clampedRoots(plate, cfg, wScale, 0.02 * A, 1, 3)[0], 2.347)) < 0.5);
});

test('a disk free to tilt, with and without inertia (scipy, n = 1 rigid hub)', () => {
  const { plate, cfg, wScale } = circle('free', 0.3, 30);
  const rhoH = cfg.rho * cfg.h;
  for (const [beta, j, want] of [[0.05, 0, [20.4758, 59.8283, 119.0477]], [0.05, 0.01, [13.0036, 27.2707, 70.783]], [0.1, 0, [20.4939, 60.0578, 120.1765]], [0.1, 0.01, [14.5685, 32.2404, 79.3431]]]) {
    const r0 = beta * A, inertia = j * rhoH * Math.PI * A ** 4; // j = I/(ρhπa⁴), the plate's own hub added in the exact solution
    const patch = new RigidPatch(plate, cfg, { r0, mount: 'shaker', inertia, harmonics: [[1, 0]] });
    const tiltRow = L => {
      const w = L * wScale;
      const [stiffness] = patch.stiffness(w);
      return stiffness[patch.dofs ** 2 - 1] - w * w * inertia / (r0 * r0);
    };
    const roots = [];
    let f0 = tiltRow(1);
    for (let L = 1.05; L < want.at(-1) * 1.05; L += 0.05) {
      const f = tiltRow(L);
      if (f0 * f < 0) {
        const root = bisect(tiltRow, L - 0.05, L, f0);
        if (Math.abs(tiltRow(root + 1e-7)) + Math.abs(tiltRow(root - 1e-7)) < Math.abs(f0) + Math.abs(f)) roots.push(root); // not a pole
      }
      f0 = f;
    }
    want.forEach((L, s) => assert.ok(Math.abs(percent(roots[s], L)) < 0.25, `b/a ${beta} j ${j} root ${s}: ${roots[s]} vs ${L}`));
  }
});

test('off-centre and rectangle code paths: consistent with the centred paths and an independent Navier series', () => {
  // A circle driven a hair off centre goes through the three-harmonic path and must still find the exact n = 0 mode (free, b/a 0.1)
  const { plate, cfg, wScale } = circle('free', 0.3, 32);
  plate.setDrive(1e-5, 0);
  const patch = new RigidPatch(plate, cfg, { r0: 0.1 * A, mount: 'rigid' });
  assert.ok(patch.basis.some(({ h }) => h === 1) && patch.dofs === 1, 'three-harmonic path, tilt locked');
  const poles = plate.modes.filter(m => m.n === 0).map(m => m.L).sort((a, b) => a - b);
  const bolted = L => {
    const [stiffness] = patch.stiffness(L * wScale);
    return 1 / stiffness[0];
  };
  const found = zerosBetween(bolted, poles, 30); // also holds the anti-resonances, where 1/A₀₀ changes sign through infinity
  for (const L of [4.2374, 25.2618]) assert.ok(Math.abs(percent(nearest(found, L), L)) < 0.5, `n = 0 off centre: ${nearest(found, L)} vs ${L}`);

  // Simply supported unit square clamped over a disk at its centre. scipy: Navier modes, m, n odd ≤ 801, no tail at all.
  const square = { E: 69e9, rho: 2700, nu: 0.3, h: 0.001, shape: 'rect', bc: 'simply', d: 1, Lx: 1, Ly: 1, zeta: 0, G: 9, mindlin: false, fmax: 100 };
  const s = new Plate();
  s.configure(square);
  s.setDrive(0, 0);
  const scale = s.wScale;
  const evenPoles = [...new Set(s.modes.filter(m => m.px === 0 && m.py === 0 && m.L > 0).map(m => m.L.toFixed(6)))].map(Number).sort((a, b) => a - b);
  for (const [r0, want] of [[0.02, [55.163, 99.334, 151.383, 213.842]], [0.05, [61.387, 102.664, 160.021, 227.868]]]) {
    const held = new RigidPatch(s, square, { r0, mount: 'rigid' });
    assert.ok(held.basis.every(({ h }) => h === 0 || h === 2), 'even harmonics only');
    const roots = zerosBetween(L => determinant(held.receptance(L * scale)[0], held.U), evenPoles, 240);
    want.forEach(L => assert.ok(Math.abs(percent(nearest(roots, L), L)) < 0.4, `square r₀ ${r0}: ${nearest(roots, L)} vs ${L}`));
  }
});

// Escaler & De La Torre 2018 (doi 10.1016/j.jfluidstructs.2018.07.017): values as in physics.test.js. Paper not opened here.
const ESCALER = { shape: 'circle', bc: 'free', E: 69e9, rho: 2700, nu: 0.35, h: 0.0008, d: 0.24, Lx: 0.24, Ly: 0.24, fmax: 5000, G: 25 };
const ESCALER_MA = 0.08;
const ESCALER_F = [70, 305, 844, 1657, 2747, 4129];
const SUSPENDED_EXACT = [35.72, 85.01, 323.67, 894.45]; // scipy: hub r₀ = 4 mm, m_a = 0.08 kg, k_s = 12.3 N/mm
const ESCALER_RADII = [[0.456], [0.813, 0.132], [0.873, 0.526, 0.063], [0.907, 0.652, 0.375, 0.05], [0.926, 0.727, 0.513, 0.295, 0.044], [0.939, 0.778, 0.602, 0.422, 0.245, 0.041]];

function escalerResonances(r0, mindlin) {
  const cfg = { ...ESCALER, zeta: 0, mindlin };
  const plate = new Plate();
  plate.configure(cfg);
  plate.setDrive(0, 0);
  const patch = new RigidPatch(plate, cfg, { r0, mount: 'shaker', ma: ESCALER_MA });
  const Z = f => {
    const w = TAU * f;
    return patch.stiffness(w)[0][0] - w * w * ESCALER_MA;
  };
  const found = [];
  let f0 = 20, z0 = Z(f0);
  for (let f = f0 * 1.0005; f < 5000; f *= 1.0005) { // finer than the gap between a mode of the held plate and its neighbour
    const z = Z(f);
    if (z0 * z < 0) {
      const root = bisect(Z, f / 1.0005, f, z0);
      if (Math.abs(Z(root * (1 + 1e-9))) + Math.abs(Z(root * (1 - 1e-9))) < Math.abs(z0) + Math.abs(z)) found.push(root); // a pole would blow up
    }
    z0 = z;
  }
  return found;
}

function nodalCircles(r0, zeta, f, mindlin) {
  const cfg = { ...ESCALER, zeta, mindlin };
  const plate = new Plate();
  plate.configure(cfg);
  plate.setDrive(0, 0);
  const patch = new RigidPatch(plate, cfg, { r0, mount: 'shaker', ma: ESCALER_MA });
  const cr = new Float64Array(plate.modes.length), ci = new Float64Array(plate.modes.length);
  const w = TAU * f;
  patch.coefficients(w, 1, g, cr, ci);
  const start = Math.max(0.003, 1.05 * r0 / 0.12), N = 1000;
  const radius = j => start + (0.998 - start) * j / N;
  const y = Array.from({ length: N + 1 }, (_, j) => Math.hypot(...patch.responseAt(radius(j), 0, w, g, cr, ci)));
  const minima = [];
  for (let j = N - 1; j > 0; j--) {
    if (y[j] < y[j - 1] && y[j] <= y[j + 1]) minima.push(radius(j + 0.5 * (y[j - 1] - y[j + 1]) / (y[j - 1] - 2 * y[j] + y[j + 1])));
  }
  return minima;
}

test('Escaler plate: the patch reproduces the exact hub with its moving mass (scipy), including the nodal circles', () => {
  const exact = { // r₀: frequencies in Hz and nodal radii, thin plate
    0.004: [[72.67, 322.97, 894.39, 1757.48, 2913.86, 4362.42], [[0.4427], [0.8034, 0.1285], [0.8761, 0.5276, 0.066], [0.9108, 0.659, 0.3819, 0.0503], [0.9303, 0.7342, 0.5174, 0.3034, 0.044], [0.9428, 0.7823, 0.605, 0.4287, 0.2537, 0.0408]]],
    0.008: [[75.36, 345.82, 960.47, 1888.56, 3131.43, 4687.78], [[0.4571], [0.8095, 0.1626], [0.8803, 0.5436, 0.102], [0.9139, 0.6709, 0.4033, 0.0859], [0.9327, 0.7435, 0.5343, 0.3277, 0.0791], [0.9448, 0.7899, 0.6189, 0.4488, 0.2798, 0.0755]]],
  };
  for (const [key, [frequencies, radii]] of Object.entries(exact)) {
    const r0 = Number(key);
    const found = escalerResonances(r0, false);
    assert.equal(found.length, 6, `r₀ ${r0}: ${found}`);
    frequencies.forEach((f, i) => assert.ok(Math.abs(percent(found[i], f)) < 0.3, `r₀ ${r0} s ${i + 1}: ${found[i]} vs ${f}`));
    // nodal circles on the mindlin-free model, lightly damped, compared where the tail has settled
    radii.forEach((row, i) => {
      const circles = nodalCircles(r0, 5e-4, found[i], false);
      row.forEach(r => assert.ok(Math.abs(nearest(circles, r) - r) < 0.003, `r₀ ${r0} s ${i + 1}: nodal ${r}, found ${circles}`));
    });
  }
});

test('Escaler plate against measurements: innermost circles now match, absolute frequencies get worse', () => {
  const inner = rows => rows.filter(row => row.length > 1).map(row => row.at(-1)); // the innermost circle of s = 2…6
  const errors = r0 => {
    const f = escalerResonances(r0, true);
    return {
      frequency: f.map((x, i) => percent(x, ESCALER_F[i])),
      circles: inner(ESCALER_RADII).map((r, i) => Math.abs(nearest(nodalCircles(r0, 0.005, f[i + 1], true), r) - r)),
    };
  };
  const point = errors(0), patch = errors(0.004);
  // measured s = 2…6 innermost: 0.132 0.063 0.050 0.044 0.041. Point: 0.101 0.032 0.015 0.009 0.006. r₀ = 4 mm: 0.129 0.066 0.051 0.045 0.041.
  point.circles.forEach((e, i) => assert.ok(e > 0.025, `point s ${i + 2}: ${e}`));
  patch.circles.forEach((e, i) => assert.ok(e < 0.005, `patch s ${i + 2}: ${e}`));
  // Frequencies at nominal thickness and modulus: point +1.8 … −0.3 %; r₀ = 4 mm +3.8 … +6.0 %, an almost uniform offset for s ≥ 2.
  assert.ok(Math.max(...point.frequency) < 1.9 && Math.min(...point.frequency) > -0.5, `point ${point.frequency}`);
  assert.ok(Math.min(...patch.frequency) > 3.5 && Math.max(...patch.frequency) < 6.3, `patch ${patch.frequency}`);
  const spread = patch.frequency.slice(1);
  assert.ok(Math.max(...spread) - Math.min(...spread) < 0.5, `spread above s = 1: ${spread}`);
});

test('Shu et al. 2022 (doi 10.3390/e24020215): no washer radius helps, and none replaces the fitted moving mass', () => {
  const measured = [263, 362, 585, 697, 775, 990, 1247, 1373, 1545, 1678, 1938, 2255, 2485, 2808, 2917, 3101, 3202, 3409, 3646, 3879, 4174, 4369, 4481, 4733, 4817];
  const cfg = { shape: 'rect', bc: 'free', E: 70e9, rho: 2700, nu: 0.33, h: 0.001, d: 0.32, Lx: 0.32, Ly: 0.32, zeta: 0.005, fmax: 5000, G: 25 };
  const plate = new Plate();
  plate.configure(cfg);
  plate.setDrive(0, 0);
  const cr = new Float64Array(plate.modes.length), ci = new Float64Array(plate.modes.length);
  const errors = options => {
    const patch = new RigidPatch(plate, cfg, options);
    const response = [];
    for (let f = 150; f < 5200; f *= 1.001) {
      patch.coefficients(TAU * f, 1, g, cr, ci);
      let energy = 0;
      for (let k = 0; k < cr.length; k++) energy += cr[k] * cr[k] + ci[k] * ci[k];
      response.push([f, Math.sqrt(energy)]);
    }
    const peaks = [];
    for (let i = 10; i + 10 < response.length; i++) {
      const y = response[i][1];
      if (y > response[i - 1][1] && y >= response[i + 1][1] && y > 1.2 * response[i - 10][1] && y > 1.2 * response[i + 10][1]) peaks.push(response[i][0]);
    }
    return measured.map(f => percent(nearest(peaks, f), f));
  };
  const fitted = errors({ r0: 0, mount: 'shaker', ma: 0.02 });
  [-2.4, 0, 0.5].forEach((want, i) => assert.ok(Math.abs(fitted[i] - want) < 0.2, `point f${i + 1}: ${fitted[i]}`)); // as physics.test.js
  const washer = errors({ r0: 0.004, mount: 'shaker', ma: 0.02 });
  assert.ok(rms(washer) > rms(fitted) + 0.2, `RMS ${rms(washer)} vs ${rms(fitted)}`); // 1.85 vs 1.45 %
  // Without a moving mass, even a 10 mm clamp leaves f1 far too low: the mass, not the clamp, sets it.
  assert.ok(errors({ r0: 0.01, mount: 'rigid' })[0] < -14);
});

test('suspended shaker: the first resonance on its suspension matches the exact hub with mass and spring (scipy)', () => {
  // k_s = 12.3 N/mm: LDS V406/8 datasheet (0.2 kg moving element), applied here to the Escaler plate with its own 0.08 kg as a typical value.
  const ks = 12300, cfg = { ...ESCALER, zeta: 1e-4, mindlin: false };
  const plate = new Plate();
  plate.configure(cfg);
  plate.setDrive(0, 0);
  const patch = new RigidPatch(plate, cfg, { r0: 0.004, mount: 'shaker', ma: ESCALER_MA, ks });
  const cr = new Float64Array(plate.modes.length), ci = new Float64Array(plate.modes.length);
  const motion = f => {
    patch.coefficients(TAU * f, 1, g, cr, ci);
    return Math.hypot(patch.qr[0], patch.qi[0]);
  };
  for (const want of SUSPENDED_EXACT) {
    let best = want, top = 0;
    for (let f = want * 0.97; f < want * 1.03; f *= 1.0002) {
      const y = motion(f);
      if (y > top) {
        top = y;
        best = f;
      }
    }
    assert.ok(Math.abs(percent(best, want)) < 0.3, `${best} vs ${want}`);
  }
});
