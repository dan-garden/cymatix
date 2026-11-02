// Verification against A. W. Leissa, "Vibration of Plates", NASA SP-160 (1969): every table in chapters 2 and 4 that
// covers what the app models (isotropic, uniform circular and rectangular plates with free, simply supported or clamped
// edges, point supports and central masses). Values were read from the scanned pages. Where a printed value is
// inconsistent with the rest of the book, the test says so and checks it against the corroborating source.
import test from 'node:test';
import assert from 'node:assert/strict';
import { circleModes, nodalRadii, Plate } from './physics.js';
import { WALLER_RADII, elasticRectModes, findMode, forEachCell, median, nearest, pct, within } from './testkit.js';

const L = (ms, n, s) => findMode(ms, n, s).L;
const radii = (ms, n, s) => nodalRadii(findMode(ms, n, s)).reverse(); // outside in, as printed
// Errors of a printed table of radii, indexed [s − 1][n], against ours.
const radiusErrors = (ms, table) => {
  const errors = [];
  table.forEach((cols, s) => cols.forEach((printed, n) => {
    const r = radii(ms, n, s + 1);
    printed.forEach((x, k) => errors.push(Math.abs(r[k] - x)));
  }));
  return errors;
};

test('Leissa ch. 2: clamped circular plate, Tables 2.1 and 2.2', () => {
  const ms = circleModes('clamped', 0.3, 40);
  const pairs = [];
  const t21 = [[10.2158, 21.26, 34.88, 51.04], [39.771, 60.82, 84.58, 111.01], [89.104, 120.08, 153.81, 190.30], [158.183, 199.06, 242.71, 289.17], [247.005, 297.77, 351.38, 407.72],
    [355.568, 416.20, 479.65, 545.97], [483.872, 554.37, 627.75, 703.95], [631.914, 712.30, 795.52, 881.67], [799.702, 889.95, 983.07, 1079.0], [987.216, 1087.4, 1190.4, 1296.2]];
  forEachCell(t21, (v, s, n) => pairs.push([L(ms, n, s), v, `(${n},${s})`]));
  const t22 = [[69.6659, 90.7390, 114.2126, 140.0561, 168.2445, 198.7561, 231.5732, 266.6790, 304.0601, 343.7038, 385.5996], [140.1079, 171.8029, 206.0706, 242.8782, 282.1977, 324.0036, 368.2734],
    [229.5186, 271.4283, 316.0015, 363.2097], [338.4113, 390.3896]];
  forEachCell(t22, (v, s, i) => pairs.push([L(ms, i + 4, s), v, `(${i + 4},${s})`]));
  assert.equal(pairs.length, 64);
  within(pairs, 0.03, 'λ²'); // median 0.001 %
  const d = radiusErrors(ms, [[[.379], [.4899], [.559], [.606], [.641], [.669]], [[.583, .255], [.640, .350], [.679, .414], [.708, .462], [.730, .501], [.749, .532]],
    [[.688, .439, .191], [.721, .497, .272], [.746, .540, .330], [.765, .574, .375], [.781, .601, .412], [.787, .618, .439]],
    [[.749, .550, .351, .153], [.767, .589, .407, .222], [.789, .620, .449, .274], [.803, .645, .488, .316]],
    [[.791, .625, .459, .293, .127], [.807, .653, .499, .344, .188]], [[.822, .678, .535, .393, .251, .109], [.833, .699, .566, .432, .298, .163]],
    [[.844, .720, .593, .469, .344, .220, .096], [.853, .735, .617, .499, .381, .263, .144]]]);
  assert.equal(d.length, 88);
  assert.ok(median(d) < 0.001 && Math.max(...d) < 0.008, `radii: median ${median(d)}, max ${Math.max(...d)}`); // max is (5,3), whose printed .787 breaks its row's trend
});

test('Leissa ch. 2: simply supported circular plate, Tables 2.3 and 2.4', () => {
  // Table 2.3 is labelled ν = 0.3 but its values are the ν = 1/3 roots: 0.01 % median at 1/3, up to 0.84 % at 0.3
  const third = circleModes('simply', 1 / 3, 20);
  const ms = circleModes('simply', 0.3, 20);
  const pairs = [];
  forEachCell([[4.977, 13.94, 25.65], [29.76, 48.51, 70.14], [74.20, 102.80, 134.33], [138.34, 176.84, 218.24]], (v, s, n) => {
    pairs.push([L(third, n, s), v, `(${n},${s})`]);
  });
  within(pairs, 0.15, 'λ² at ν = 1/3');
  assert.ok(pct(L(ms, 0, 0), 4.977) > 0.8, 'the ν = 0.3 fundamental would be 4.935');
  const d = radiusErrors(ms, [[[.441], [.550], [.613]], [[.644, .279], [.692, .378], [.726, .443]], [[.736, .469, .204], [.765, .528, .288], [.787, .570, .348]]]);
  assert.ok(Math.max(...d) < 0.003, `radii max ${Math.max(...d)}`);
});

test('Leissa ch. 2: completely free circular plate, Tables 2.5–2.7', () => {
  const m33 = circleModes('free', 0.33, 42);
  const m25 = circleModes('free', 0.25, 34);
  const t25 = [[null, null, 5.253, 12.23, 21.6, 33.1, 46.2], [9.084, 20.52, 35.25, 52.91, 73.1, 95.8, 121.0], [38.55, 59.86, 83.9, 111.3, 142.8, 175.0, 210.3],
    [87.80, 119.0, 154.0, 192.1, 232.3, 274.6, 319.7], [157.0, 198.2, 242.7, 290.7, 340.4, 392.4, 447.3], [245.9, 296.9, 350.8, 408.4, 467.9, 529.5, 593.9],
    [354.6, 415.3, 479.2, 546.2, 615.0, 686.4, 760.1], [483.1, 651.8, 627.0, 703.3, 781.8, 864.4, 952.3], [631.0, 711.3, 794.7, 880.3, 968.5, 1061, 1158.7],
    [798.6, 888.6, 981.6, 1076, 1175, 1277, 1384], [986.0, 1086, 1188, 1292, 1401, 1513, 1631]];
  const exact = [], rough = []; // footnote a: "values true within 2 percent" for n = 4–6, s = 0–1
  forEachCell(t25, (v, s, n) => {
    if (s === 7 && n === 1) return;
    (n >= 4 && s <= 1 ? rough : exact).push([L(m33, n, s), v, `(${n},${s})`]);
  });
  within(rough, 2, 'footnote-a cell');
  // this older table is good to a few tenths of a percent (median 0.13 %, worst 0.74 %): at (3,3) it even rises from
  // Table 2.6's ν = 0.25 value, while the exact roots fall with ν
  within(exact, 0.75, 'λ² ν = 0.33');
  assert.ok(median(exact.map(([g, w]) => pct(g, w))) < 0.15);
  // (1,7) is printed 651.8, out of sequence; ours sits with Table 2.6's 553.00 for the same mode (ν = 0.25)
  assert.ok(pct(L(m33, 1, 7), 553.00) < 0.1, `(1,7): ${L(m33, 1, 7)}`);
  const t26 = [[null, null, 5.513, 12.75], [8.892, 20.41, 35.28, 53.16], [38.34, 59.74, 84.38, 112.36], [87.65, 118.88, 153.29, 191.02], [156.73, 196.67, 241.99, 289.51],
    [245.52, 296.46, 350.48, 408.16], [354.08, 414.86, 478.73, 545.83], [482.37, 553.00, 626.75, 703.63], [630.41, 710.92, 794.51, 881.20], [798.23, 888.58, 982.01, 1078.5]];
  const p26 = [];
  forEachCell(t26, (v, s, n) => p26.push([L(m25, n, s), v, `(${n},${s})`]));
  within(p26, 0.6, 'λ² ν = 0.25');
  assert.ok(median(p26.map(([g, w]) => pct(g, w))) < 0.02);
  // where Table 2.6 is off (its n = 3 column, 0.1–0.3 %), ours equals the exact frequency equation (eq. 2.14, scipy)
  assert.ok(pct(L(m25, 3, 2), 112.0142) < 0.005 && pct(L(m25, 3, 3), 190.7208) < 0.005);
  // Table 2.7 radii (ν = 0.33). Where ours and the printed value differ by more than 0.005, Waller's measured radius
  // (Table 2.9) must be closer to ours: the printed outer circles of s = 2–3, n ≥ 1 run high.
  const t27 = [[[.680], [.781], [.822], [.847], [.863], [.881]], [[.841, .391], [.871, .4972], [.8897, .562], [.925, .605], [.926, .635], [.993, .663]],
    [[.893, .591, .257], [.932, .643, .351], [.936, .678, .414], [.939, .704, .460], [.943, .726, .498], [.947, .745, .529]],
    [[.941, .691, .441, .192], [.946, .723, .498, .272], [.950, .746, .540, .330], [.951, .763, .572, .374], [.955, .779, .600, .411], [.958, .793, .623, .443]],
    [[.952, .752, null, .352, .154], [.956, .773, .590, .407, .222], [.959, .790, .620, .449, .274], [.960, .803, .644, .483, .316], [.963, .814, .644, .512, .351], [.966, .825, .682, .536, .381]]];
  let agree = 0, sided = 0, outer5 = 0;
  t27.forEach((cols, s) => cols.forEach((printed, n) => {
    const r = radii(m33, n, s + 1);
    printed.forEach((x, k) => {
      if (x == null) return;
      if (Math.abs(r[k] - x) <= 0.005) {
        agree++;
        return;
      }
      const w = WALLER_RADII[n]?.[s]?.slice().reverse()[k]; // Waller lists inside out
      // s = 5 has no measurement for n ≥ 1; its printed outer circles sit 0.012–0.017 beyond exact theory, like the rest
      // of that column (and like Waller's s = 5–6, n = 0 outermost circles, while all his inner ones match within 0.004)
      if (w === undefined && s === 4 && k === 0) {
        assert.ok(x - r[k] > 0 && x - r[k] < 0.02, `(${n},5) outer: ${r[k]} vs ${x}`);
        outer5++;
        return;
      }
      if (w === undefined) { // printed .644 repeats n = 3
        assert.ok(n === 4 && s === 4 && k === 2, `unexplained (${n},${s + 1}) circle ${k + 1}: ${r[k]} vs ${x}`);
        return;
      }
      assert.ok(Math.abs(r[k] - w) < Math.abs(x - w), `(${n},${s + 1}) circle ${k + 1}: ours ${r[k]}, printed ${x}, measured ${w}`);
      sided++;
    });
  }));
  // of 88 printed radii: 67 agree within 0.005, 15 differ with Waller's measurement siding with ours, 6 are s = 5 outer
  assert.ok(agree >= 67 && sided + outer5 <= 21, `${agree} agree, ${sided} settled by Waller, ${outer5} outer s = 5`);
  // the exact radii, independently from eq. 2.14 in scipy, for the (0,5) mode
  radii(m33, 0, 5).reverse().forEach((r, k) => assert.ok(Math.abs(r - [0.153, 0.352, 0.552, 0.753, 0.935][k]) < 0.001));
});

// resonance nearest each target from a fine sweep of the plate's own response (thin theory, so λ² maps exactly to Hz)
const resonances = (cfg, targets) => {
  const p = new Plate();
  p.configure({ E: 70e9, rho: 2700, h: 0.001, d: 0.24, Lx: 0.24, Ly: 0.24, zeta: 1e-4, fmax: 9e4, G: 9, mindlin: false, ...cfg });
  p.setDrive(0, 0);
  return targets.map(L2 => {
    const f0 = L2 * p.wScale / (2 * Math.PI);
    let best = 0, bestF = f0;
    for (let f = f0 * 0.97; f < f0 * 1.03; f *= 1.00005) {
      const y = p.field([{ f, F: 1 }], 9.81);
      if (y > best) {
        best = y;
        bestF = f;
      }
    }
    return bestF * 2 * Math.PI / p.wScale;
  });
};

test('Leissa ch. 2: point supports and a central mass, through the app\'s driver models', () => {
  const compare = (cfg, targets, maxPct, name) => within(resonances(cfg, targets).map((g, i) => [g, targets[i]]), maxPct, name);
  const heldAtCentre = bc => ({ shape: 'circle', bc, nu: 0.3, mount: 'rigid' });
  // free plate held at its centre (rigid bolt): exact roots of the point-support equation (solved in scipy, ν = 0.3)
  compare(heldAtCentre('free'), [3.752, 20.908, 60.645, 119.814], 0.05, 'centre-fixed free plate');
  // Table 2.14 (exact column, ν = 0.3) and the text's centre-supported simply supported and clamped plates (3 figures)
  compare(heldAtCentre('free'), [3.73, 20.9, 60.5, 119.7], 0.7, 'Table 2.14');
  compare(heldAtCentre('simply'), [14.8, 49.4], 0.5, 'simply supported + centre support');
  compare(heldAtCentre('clamped'), [22.7, 61.9], 0.5, 'clamped + centre support');
  // Table 2.15: clamped plate with a mass μ·πρa² at the centre = the shaker model with that moving mass and no force.
  // Exact roots (point-mass Green's function J0, I0, Y0 + (2/π)K0, solved in scipy) are matched to 0.01 %; Leissa's
  // "precise values" are up to 0.34 % low on the second root.
  for (const [mu, exact, printed] of [[0.05, [9.0100, 32.9454, 72.0069, 129.5534], [9.0120, 32.833, 72.012, 129.39]],
    [0.10, [8.1093, 29.7733, 67.7337, 125.6151], [8.1111, 29.681, 67.733, 125.69]]]) {
    const got = resonances({ shape: 'circle', bc: 'clamped', nu: 0.3, mount: 'shaker', ma: mu * Math.PI * 2700 * 0.001 * 0.12 ** 2 }, exact);
    within(got.map((g, i) => [g, exact[i]]), 0.03, `central mass μ = ${mu}, exact`);
    within(got.map((g, i) => [g, printed[i]]), 0.35, `central mass μ = ${mu}, Table 2.15`);
  }
});

// Leissa's rectangle: 0 ≤ x ≤ a, 0 ≤ y ≤ b, x = 0 and x = a simply supported; λ = ωa²√(ρh/D), λ* = ωb²√(ρh/D) = λ(b/a)².
const S = 'simply', C = 'clamped', F = 'free';
const EDGES = { SCSC: [S, S, C, C], SCSS: [S, S, S, C], SCSF: [S, S, C, F], SSSF: [S, S, S, F], SFSF: [S, S, F, F] };
const rect = (k, nu, beta) => elasticRectModes(EDGES[k], nu, beta);

test('Leissa ch. 4: SS–C–SS–C rectangles, Tables 4.3–4.8', () => {
  const sq = rect('SCSC', 0.3, 1).map(m => m.L);
  within([28.946, 54.743, 69.320, 94.584, 102.213, 129.086].map((v, i) => [sq[i], v]), 0.03, 'T4.3');
  within([140.189, 154.765, 199.797, 208.373, 234.578, 279.627, 307.300, 333.926, 379.274, 425.885, 452.877, 498.501].map(v => [nearest(sq, v), v]), 0.03, 'T4.4');
  // T4.5/T4.8: rows m (half-waves along x), columns n (n-th root across y), printed to 0.1: matched within the mode's
  // symmetry class (x-parity of m, y-parity of n)
  const grid = (beta, rows, scale, name) => {
    const ms = rect('SCSC', 0.3, beta);
    const out = [];
    forEachCell(rows, (v, m, n) => {
      const sameClass = ms.filter(q => q.px === m % 2 && q.py === n % 2).map(q => q.L * scale);
      out.push([nearest(sameClass, v), v, `${name} (${m + 1},${n + 1})`]);
    });
    return out;
  };
  const g5 = grid(1, [[28.9, 69.2, 129.1, 208.6, 307.4, 426.1], [54.8, 94.6, 154.8, 234.5, 333.9, 452.9], [102.2, 140.2, 199.9, 279.5, 379.1, 498.4],
    [170.3, 206.6, 265.2, 344.6, 443.8, 563.5], [258.5, 293.8, 351.1, 429.8, 529.0, 647.9], [366.8, 400.9, 457.4, 535.1, 633.7, 752.2]], 1, 'T4.5');
  const g8 = [[0.5, [[54.8, 94.6, 154.8, 234.5, 333.9, 452.9], [170.3, 206.6, 265.2, 344.6, 443.8, 563.5], [366.8, 400.9, 457.4, 535.1, 633.7, 752.2], [642.8, 675.9, 730.5, 806.9, 904.2, 1021], [997.7, 1030, 1084, 1159, 1257, 1375], [1432, 1464, 1517, 1592, 1686, 1802]]],
    [1.5, [[25.0, 64.9, 124.5, 203.7, 302.4, 420.9], [35.1, 75.6, 135.7, 215.1, 314.1, 432.8], [54.8, 94.6, 154.8, 234.5, 333.9, 452.9], [84.1, 122.3, 182.6, 262.5, null, 481.1], [122.6, 160.0, 219.3, 298.9, 398.5, 518.0], [170.3, 206.6, 265.2, 344.6, 443.8, 563.5]]],
    [2.0, [[23.8, 63.4, 123.0, 202.1, 300.7, 419.0], [28.9, 69.2, 129.1, 208.6, 307.4, 426.1], [39.0, 79.5, 139.7, 219.3, 318.2, 437.1], [54.8, 94.6, 154.8, 234.5, 333.9, 452.9], [75.9, 114.7, 174.6, 254.7, 354.1, 473.3], [102.2, 140.2, 199.9, 279.5, 379.1, 498.4]]]]
    .flatMap(([ab, rows]) => grid(1 / ab, rows, 1 / ab ** 2, `T4.8 a/b ${ab}`));
  assert.equal(g5.length + g8.length, 143);
  within([...g5, ...g8], 0.35, ''); // printed to 0.1: median 0.03 %
  within([1, 1.5, 2, 2.5, 3].map((ba, i) => [rect('SCSC', 0.3, ba)[0].L, [28.946, 17.369, 13.688, 12.129, 11.359][i], `b/a ${ba}`]), 0.06, 'T4.6 λ');
  // T4.6 λ* at a/b = 1.5 is printed 24.047, off a smooth curve through its neighbours; ours is 25.044 (a 4/5 misprint)
  within([1, 2, 2.5, 3].map((ab, i) => [rect('SCSC', 0.3, 1 / ab)[0].L / ab ** 2, [28.946, 23.814, 23.271, 22.985][i], `a/b ${ab}`]), 0.06, 'T4.6 λ*');
  assert.ok(pct(rect('SCSC', 0.3, 1 / 1.5)[0].L / 2.25, 25.047) < 0.05);
  // T4.7, the lowest mode odd across y: a/b = 1 is printed 68.181, but Table 4.3 gives this same mode (ω12) as 69.320
  const odd = ab => rect('SCSC', 0.3, 1 / ab).filter(m => m.px === 0 && m.py === 1)[0].L / ab ** 2;
  within([1.5, 2, 2.5, 3].map((ab, i) => [odd(ab), [65.118, 63.641, 62.967, 62.602][i], `a/b ${ab}`]), 0.5, 'T4.7');
  assert.ok(pct(odd(1), 69.320) < 0.02);
});

test('Leissa ch. 4: SS–C–SS–SS, SS–C–SS–F, SS–SS–SS–F and SS–F–SS–F rectangles, Tables 4.10–4.20', () => {
  const sq = k => rect(k, 0.3, 1).map(m => m.L);
  within([23.646, 51.674, 58.641, 86.126, 100.259, 113.217].map((v, i) => [sq('SCSS')[i], v]), 0.02, 'T4.10');
  within([133.784, 140.840, 188.102].map(v => [nearest(sq('SCSS'), v), v]), 0.02, 'T4.11');
  within([1, 1.5, 2, 2.5, 3].map((ba, i) => [rect('SCSS', 0.3, ba)[0].L, [23.646, 15.573, 12.918, 11.754, 11.142][i]]), 0.05, 'T4.13 λ');
  within([1, 1.5, 2, 2.5, 3].map((ab, i) => [rect('SCSS', 0.3, 1 / ab)[0].L / ab ** 2, [23.646, 18.899, 17.330, 16.629, 16.254][i]]), 0.05, 'T4.13 λ*');
  within([12.69, 33.06, 41.70, 63.01, 72.40, 90.61].map((v, i) => [sq('SCSF')[i], v]), 0.05, 'T4.14');
  within([11.68, 27.76, 41.20, 59.07, 61.86, 90.29].map((v, i) => [sq('SSSF')[i], v]), 0.05, 'T4.16');
  within([9.631, 16.13, 36.72, 38.94, 46.74, 70.75].map((v, i) => [sq('SFSF')[i], v]), 0.05, 'T4.18');
  // fundamentals against aspect ratio, ν = 0.25 (older calculations, ~0.25 % high on average)
  within([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2].map((ab, i) => [rect('SCSF', 0.25, 1 / ab)[0].L, [12.859, 13.520, 14.310, 15.198, 16.086, 17.172, 18.258, 19.343, 20.527, 21.910, 23.192, 26.153][i], `a/b ${ab}`]), 0.4, 'T4.15');
  // T4.17 against the exact Lévy solution (two opposite edges simply supported, solved in scipy): ours matches it to
  // 0.01 %. The printed table is an older calculation, ~0.2 % high, with three wrong cells: 11.349 at a/b = 0.6 and
  // 11.547 at 0.8 (exact 10.5800, 11.1367), and 35.133 at 5 (exact 34.6977).
  const levy = [[0.5, 10.3556], [0.6, 10.5800], [0.8, 11.1367], [1.0, 11.8196], [1.2, 12.6084], [1.4, 13.4846], [2.0, 16.4907], [3.0, 22.2171], [4.0, 28.3603], [5.0, 34.6977]];
  within(levy.map(([ab, v]) => [rect('SSSF', 0.25, 1 / ab)[0].L, v, `a/b ${ab}`]), 0.01, 'T4.17 exact');
  within([0.5, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0].map((ab, i) => [rect('SSSF', 0.25, 1 / ab)[0].L, [10.362, 11.843, 12.632, 13.520, 14.409, 15.396, 16.481, 19.244, 22.205, 28.324][i], `a/b ${ab}`]), 0.5, 'T4.17');
  // T4.19 (ν = 0.16) and T4.20 (ν = 0.3), Jankovic: every mode with n ≥ 2. Leissa notes their (m,1) values are wrong
  // (k² > α² assumed); those are left out, and T4.18 above covers the square's correct (1,1) and (2,1).
  const tables = [[0.16, { 0.5: '12.01(12) 18.21(13) 28.37(14) 41.69(22) 43.03(15) 48.66(23) 60.01(24) 62.43(16) 75.67(25) 86.81(17) 91.02(32) 95.77(26) 98.17(33) 109.99(34) 116.15(18) 120.44(27) 126.34(35) 147.21(36) 149.85(28) 150.47(19)',
    0.6: '12.87(12) 21.39(13) 35.69(14) 42.68(22) 52.49(23) 56.62(15) 68.42(24) 84.60(16) 90.49(25) 92.02(32) 102.22(33) 118.97(34) 119.01(26) 119.76(17) 142.12(35) 154.24(27)',
    0.667: '13.51(12) 23.21(13) 41.19(14) 43.42(22) 55.36(23) 67.00(15) 74.70(24) 92.78(32) 101.57(16) 101.61(25) 105.27(33) 125.74(34) 136.58(26) 145.23(17) 153.99(35)',
    0.75: '14.33(12) 26.90(13) 44.44(22) 48.77(14) 59.26(23) 81.41(15) 83.26(24) 93.83(32) 109.49(33) 116.91(25) 125.24(16) 135.03(34)',
    0.8: '14.85(12) 28.92(13) 45.10(22) 53.70(14) 61.76(23) 88.78(24) 90.86(15) 94.53(32) 112.22(33) 126.82(25) 140.77(16) 141.03(34)',
    1: '17.06(12) 37.95(13) 48.09(22) 72.88(23) 76.36(14) 97.69(32) 113.60(24) 124.56(33) 134.58(15)',
    1.333: '21.10(12) 53.96(22) 56.46(13) 94.99(23) 104.28(32) 124.77(14) 149.60(33)', 1.5: '23.21(12) 57.27(22) 67.42(13) 107.63(23) 108.11(32) 154.00(14)',
    2: '29.77(12) 68.24(22) 107.35(13) 121.40(32) 151.81(23)' }],
  [0.3, { 0.5: '11.66(12) 17.66(13) 27.73(14) 41.18(22) 42.35(15) 47.73(23) 59.06(24) 61.85(16) 74.51(25) 86.27(17) 90.26(32) 94.46(26) 97.33(33) 108.91(34) 115.67(18) 119.10(27) 125.01(35) 145.63(36) 148.49(28) 150.05(19)',
    0.6: '12.43(12) 20.69(13) 34.93(14) 42.04(22) 51.54(23) 55.90(15) 67.14(24) 83.97(16) 88.99(25) 91.12(32) 101.11(33) 117.39(26) 117.50(34) 119.19(17) 140.29(35) 152.63(27)',
    0.667: '12.98(12) 22.95(13) 40.33(14) 42.67(22) 54.23(23) 66.21(15) 73.19(24) 91.76(32) 99.88(25) 100.91(16) 103.98(33) 123.99(34) 134.78(26) 144.44(17) 151.81(35)',
    0.75: '13.71(12) 25.99(13) 43.56(22) 47.85(14) 57.89(23) 80.59(15) 81.47(24) 92.67(32) 107.92(33) 114.93(25) 124.55(16) 132.89(34) 158.91(26)',
    0.8: '14.17(12) 27.94(13) 44.15(22) 52.73(14) 60.24(23) 86.82(24) 90.03(15) 93.28(32) 110.47(33) 124.70(25) 138.65(16) 140.07(34)',
    1: '16.13(12) 36.71(13) 46.73(22) 70.73(23) 75.72(14) 96.29(32) 111.02(24) 122.04(33) 133.69(15)',
    1.333: '19.73(12) 51.93(22) 54.62(13) 91.81(23) 101.80(32) 123.56(14) 145.57(33)', 1.5: '21.62(12) 54.84(22) 65.79(13) 103.97(23) 105.16(32) 152.77(14) 159.16(33)',
    2: '27.52(12) 64.54(22) 105.48(13) 116.83(32) 146.90(23)' }]];
  const errs = [];
  for (const [nu, byAspect] of tables) {
    for (const [ab, row] of Object.entries(byAspect)) {
      const ms = rect('SFSF', nu, 1 / +ab);
      for (const tok of row.split(' ')) {
        const [, v, m] = tok.match(/([\d.]+)\((\d)(\d)\)/);
        const w = +v;
        const sameClass = ms.filter(q => q.px === (+m - 1) % 2).map(q => q.L);
        errs.push([pct(nearest(sameClass, w), w), `ν ${nu} a/b ${ab} ${tok}`]);
      }
    }
  }
  assert.ok(errs.length > 200, `${errs.length} values`);
  assert.ok(median(errs.map(e => e[0])) < 0.02, 'median');
  // only four cells exceed 0.3 %: three near 0.6 %, and ν = 0.16 a/b = 0.667 (1,3) at 2.4 %, which breaks the ν-trend
  // between the two tables (ours 23.76; the same cell is 22.95 at ν = 0.3, and the neighbours differ by 0.7–0.9)
  const big = errs.filter(e => e[0] > 0.3);
  assert.ok(big.length <= 4 && big.every(e => e[0] < 0.7 || e[1] === 'ν 0.16 a/b 0.667 23.21(13)'), big.map(e => `${e[1]} ${e[0].toFixed(2)} %`).join('; '));
});
