import test from 'node:test';
import assert from 'node:assert/strict';
import { bessel, circleModes, rectModes, Plate, findPeaks, mindlinOmega, nodalRadii, NR } from './physics.js';

const close = (got, want, rel, msg) =>
  assert.ok(Math.abs(got - want) <= rel * Math.max(Math.abs(want), 1e-12), `${msg}: got ${got}, want ${want}`);

// Power series as an independent reference: J_n(x) = Σ (−1)^k (x/2)^(2k+n) / (k!(n+k)!), I_n without the sign.
const series = (n, x, sign) => {
  let term = (x / 2) ** n, sum = 0;
  for (let i = 1; i <= n; i++) term /= i;
  for (let k = 0; k < 80; k++) { sum += term; term *= sign * (x / 2) ** 2 / ((k + 1) * (n + k + 1)); }
  return sum;
};

test('Bessel functions match their power series and derivative identities', () => {
  close(bessel(0, 1)[0], 0.7651976865579666, 1e-13, 'J0(1)');
  for (const [n, x] of [[0, 2.4], [1, 2.5], [5, 10], [12, 7], [3, 0.3]]) {
    const [J, , I] = bessel(n, x);
    close(J, series(n, x, -1), 1e-9, `J${n}(${x})`);
    close(I, series(n, x, 1), 1e-12, `I${n}(${x})`);
  }
  for (const [n, x] of [[3, 7.3], [12, 40]]) {
    const [, Jp, , Ip] = bessel(n, x), a = bessel(n - 1, x), b = bessel(n + 1, x);
    close(Jp, (a[0] - b[0]) / 2, 1e-10, "J'");
    close(Ip, (a[2] + b[2]) / 2, 1e-10, "I'");
  }
});

test('circular plate eigenvalues λ² match the closed-form frequency equations', () => {
  const pick = (modes, n, m) => modes.find(q => q.n === n && q.m === m && q.L > 0).L;
  const check = (modes, rows, tol, name) => rows.forEach(([n, m, L]) => close(pick(modes, n, m), L, tol, `${name} (${n},${m})`));
  check(circleModes('clamped', 0.3, 12), [[0, 0, 10.2158], [1, 0, 21.2604], [2, 0, 34.8770], [0, 1, 39.7711], [3, 0, 51.0300], [1, 1, 60.8287]], 1e-4, 'clamped');
  // simply supported: roots of Leissa's J_{n+1}/J_n + I_{n+1}/I_n = 2λ/(1−ν), solved independently with scipy
  check(circleModes('simply', 0.3, 12), [[0, 0, 4.935149], [1, 0, 13.898165], [2, 0, 25.613297], [0, 1, 29.720005], [3, 0, 39.957314], [1, 1, 48.478928]], 1e-6, 'simply supported');
  const free = circleModes('free', 0.33, 12);
  // free: roots of Leissa's eq. 2.14 solved with scipy (his 1969 table rounds these with ~0.3 % error)
  check(free, [[2, 0, 5.262037], [0, 1, 9.068899], [3, 0, 12.243894], [1, 1, 20.513087], [4, 0, 21.527217], [2, 1, 35.242521], [0, 2, 38.506969]], 1e-6, 'free');
  assert.equal(free.filter(m => m.L === 0).length, 2, 'free plate has two rigid-body modes');
});

test('square plate Ritz eigenvalues match exact and published values', () => {
  const L = (bc, nu) => rectModes(bc, nu, 1).modes.map(m => m.L);
  const exact = [];
  for (let m = 1; m < 40; m++) for (let n = 1; n < 40; n++) exact.push(Math.PI ** 2 * (m * m + n * n));
  exact.sort((a, b) => a - b);
  const ss = L('simply', 0.3);
  ss.forEach((v, k) => close(v, exact[k], 2e-3, `simply supported #${k + 1}`)); // every kept mode is converged
  [35.985, 73.394, 73.394, 108.22, 131.58, 132.20].forEach((v, k) => close(L('clamped', 0.3)[k], v, 1e-3, `clamped #${k + 1}`));
  const free = L('free', 0.3);
  free.slice(0, 3).forEach(v => assert.ok(v < 1e-3, 'three rigid-body modes'));
  [13.468, 19.596, 24.270, 34.801, 34.801, 61.093, 61.093, 63.686].forEach((v, k) => close(free[k + 3], v, 2e-3, `free #${k + 1}`));
});

test('quasi-static centre deflection under a centre point load matches plate theory', () => {
  const mat = { E: 69e9, rho: 2700, nu: 0.3, h: 0.001, zeta: 0, G: 41, fmax: 40000 }, g = 9.81, f = 0.5, w = 2 * Math.PI * f;
  const D = mat.E * mat.h ** 3 / (12 * (1 - mat.nu ** 2)), centre = (41 * 41 - 1) / 2;
  const deflection = shape => {
    const p = new Plate();
    p.configure({ ...mat, ...shape });
    p.setDrive(0, 0);
    p.field([{ f, F: 1 }], g);
    return p.gamma[centre] * g / (w * w);
  };
  // clamped circle: w = F a² / (16πD); simply supported square: w = 0.01160 F a² / D (Timoshenko, ν = 0.3)
  close(deflection({ shape: 'circle', bc: 'clamped', d: 0.24, Lx: 0.24, Ly: 0.24 }), 0.12 ** 2 / (16 * Math.PI * D), 5e-3, 'clamped circle');
  close(deflection({ shape: 'rect', bc: 'simply', d: 0.24, Lx: 0.24, Ly: 0.24 }), 0.0116 * 0.24 ** 2 / D, 1.5e-2, 'simply supported square');
});

test('real Chladni experiment: sand-pattern resonances of a bolted aluminium square (Shu et al., Entropy 24:215, 2022)', () => {
  // Fig. 2: 320 × 320 × 1 mm aluminium, free edges, screwed to a shaker at its centre, 0.3 mm silica sand.
  // Material values are the paper's (E 70 GPa, ρ 2700, ν 0.33).
  const measured = [263, 362, 585, 697, 775, 990, 1247, 1373, 1545, 1678, 1938, 2255, 2485, 2808, 2917, 3101, 3202, 3409, 3646, 3879, 4174, 4369, 4481, 4733, 4817];
  const resonances = (mount, ma, G) => {
    const p = new Plate(), peak = [], fs = [];
    p.configure({ shape: 'rect', bc: 'free', mount, ma, E: 70e9, rho: 2700, nu: 0.33, h: 0.001, d: 0.32, Lx: 0.32, Ly: 0.32, zeta: 0.005, fmax: 5000, G });
    p.setDrive(0, 0);
    for (let f = 150; f < 5200; f *= 1.001) fs.push([f, p.field([{ f, F: 1 }], 9.81)]);
    for (let i = 10; i + 10 < fs.length; i++) { // a resonance: local maximum standing 20 % above the response 1 % either side
      const y = fs[i][1];
      if (y > fs[i - 1][1] && y >= fs[i + 1][1] && y > 1.2 * fs[i - 10][1] && y > 1.2 * fs[i + 10][1]) peak.push(fs[i][0]);
    }
    // the experiment found 25 resonances between 200 and 5000 Hz; the model finds 24 (one close pair near 2.2 kHz merges)
    const count = peak.filter(f => f > 200 && f < 5000).length;
    assert.ok(count >= 24 && count <= 25, `${mount}: resonance count ${count}`);
    return measured.map(f => peak.reduce((b, x) => Math.abs(x - f) < Math.abs(b - f) ? x : b) / f - 1);
  };
  // No free parameter: rigid bolt. Median over f4–f25 is 0.4 %; the same model scaled by ±3–10 % scores 1.7–3.2 %.
  const rigid = resonances('rigid', 0, 41);
  rigid.forEach((e, i) => assert.ok(Math.abs(e) < (i < 3 ? 0.2 : 0.045), `rigid f${i + 1}: ${(100 * e).toFixed(1)} %`));
  assert.ok(rigid.slice(3).map(Math.abs).sort((a, b) => a - b)[11] < 0.01, 'rigid median error above 1 %');
  // One parameter: the shaker's moving mass (armature + screw), not given in the paper; 20 g fits best (16–25 g within 0.1 % RMS).
  // It moves the three lowest modes from −20/−5/−7 % (rigid) to −2.4/0.0/+0.5 %.
  const shaker = resonances('shaker', 0.02, 25);
  shaker.forEach((e, i) => assert.ok(Math.abs(e) < (i < 3 ? 0.025 : 0.045), `shaker f${i + 1}: ${(100 * e).toFixed(1)} %`));
  assert.ok(Math.sqrt(shaker.reduce((s, e) => s + e * e, 0) / 25) < 0.018, 'shaker RMS error above 1.8 %'); // 1.7 %; rigid: 4.6 %
});

test('Mindlin thick-plate correction matches exact Mindlin solutions (Hashemi & Arsanjani 2005, Table 9)', () => {
  // simply supported square, ν = 0.3, β = ωa²√(ρh/D); the thin-plate wavenumber is exact here, so the correction is too
  const thin = rectModes('simply', 0.3, 1).modes.slice(0, 9).map(m => m.L);
  for (const [d, kappa2, ref] of [
    [0.1, Math.PI ** 2 / 12, [19.0584, 45.4478, 45.4478, 69.7166, 84.9263, 84.9263, 106.5154, 106.5154, 133.3705]],
    [0.2, 5 / 6, [17.4485, 38.1521, 38.1521, 55.1500, 65.1452, 65.1452, 78.6969, 78.6969, 94.7658]],
  ]) {
    const mat = { E: 1e9, rho: 1000, nu: 0.3, h: d }, D = mat.E * d ** 3 / (12 * (1 - 0.09)), c = Math.sqrt(D / (mat.rho * d)); // side a = 1
    thin.forEach((L, i) => close(mindlinOmega(L * c, mat, kappa2) / c, ref[i], 2e-5, `h/a ${d} mode ${i + 1}`));
  }
});

test('nodal circles of free circular plates: Kirchhoff, Strehlke and Waller', () => {
  const pick = (ms, n, s) => ms.find(m => m.n === n && m.m === s && m.L > 0);
  // Rayleigh, Theory of Sound §220: radii calculated by Kirchhoff (μ = 1/4) and observed by Strehlke (1855), n = 0
  const quarter = circleModes('free', 0.25, 20);
  [[[0.68062], [0.67815]], [[0.39151, 0.84200], [0.39133, 0.84149]], [[0.25679, 0.59147, 0.89381], [0.25631, 0.59107, 0.89360]]].forEach(([calc, seen], i) => {
    const r = nodalRadii(pick(quarter, 0, i + 1)); // inner to outer
    r.forEach((x, k) => { assert.ok(Math.abs(x - calc[k]) < 2e-4, `Kirchhoff s=${i + 1}: ${x} vs ${calc[k]}`); assert.ok(Math.abs(x - seen[k]) < 3e-3, `Strehlke s=${i + 1}`); });
  });
  // Waller (1938), via Leissa NASA SP-160 Table 2.9: 60 measured radii on free brass plates, n = 0–5, 1–4 circles
  const waller = [[[.680], [.391, .843], [.257, .591, .895], [.190, .441, .692, .918]], [[.781], [.497, .867], [.349, .643, .902], [.269, .495, .726, .928]],
    [[.823], [.562, .887], [.415, .681, .913], [.328, .540, .748, .934]], [[.843], [.604, .898], [.461, .706, .919], [.374, .571, .764, .938]],
    [[.859], [.635, .906], [.505, .728, .925], [.411, .596, .779, .941]], [[.871], [.662, .915], [.531, .745, .933], [.443, .623, .794, .944]]];
  const brass = circleModes('free', 0.33, 30), d = [];
  waller.forEach((rows, n) => rows.forEach((seen, s) => { const r = nodalRadii(pick(brass, n, s + 1)); seen.forEach((x, k) => d.push(Math.abs(r[k] - x))); }));
  d.sort((a, b) => a - b);
  assert.equal(d.length, 60);
  assert.ok(d[30] < 0.002 && d[59] < 0.007, `median ${d[30]}, max ${d[59]} (r/a)`); // 0.0014 and 0.0057
});

test('free circular plate frequency ratios: Waller (1938), 76 measured modes', () => {
  // Leissa NASA SP-160 Table 2.8: ω/ω(n=2, s=0) for free brass plates; rows are s = 0…6, from n = 2 (s = 0) or n = 0
  const rows = [[1, 2.29, 4.10, 6.19, 8.80, 11.7, 15.1, 18.8, 23.0, 27.6, 33.0, 38.0, 44.5, 51.0, 58.0, 65.5, 74.3, 83.5, 90.9, 100, 109, 118, 127],
    [1.70, 3.99, 6.79, 10.3, 13.8, 18.1, 23.1, 28.6, 34.7, 40.8, 47.5, 54.7, 63.0, 72.5], [7.51, 11.7, 16.1, 21.2, 27.1, 33.4, 40.5, 47.7, 56.0, 65, 75, 86, 100],
    [16.4, 22.7, 29.4, 36.3, 43.0, 50.5, 59, 69, 80, 90], [29.1, 37.0, 43.2, 53.1, 63.0, 74, 86, 99], [47.0, 56.0, 66, 77, 92], [72.7, 83, 100]];
  const ms = circleModes('free', 0.33, 40), L0 = ms.find(m => m.n === 2 && m.m === 0 && m.L > 0).L, err = [];
  rows.forEach((row, s) => row.forEach((ratio, i) => { const m = ms.find(q => q.n === (s ? i : i + 2) && q.m === s && q.L > 0); err.push(Math.abs(m.L / L0 / ratio - 1)); }));
  err.sort((a, b) => a - b);
  assert.equal(err.length, 76);
  assert.ok(err[38] < 0.02, `median ${err[38]}`); // 1.6 % (brass ν not stated; 0.30–0.35 give 1.6–1.9 %)
});

test('centre-driven circular Chladni plate (Escaler & De La Torre, J. Fluids Struct. 2018)', () => {
  // 6061-T6, a = 0.12 m, h = 0.8 mm, E 69 GPa, ν 0.35; centre driven through a bar, moving mass ≈ 0.08 kg (authors' value)
  const measured = [70, 305, 844, 1657, 2747, 4129], radii = [[0.456], [0.813, 0.132], [0.873, 0.526, 0.063], [0.907, 0.652, 0.375, 0.050], [0.926, 0.727, 0.513, 0.295, 0.044], [0.939, 0.778, 0.602, 0.422, 0.245, 0.041]];
  const p = new Plate(), fs = [], peak = [];
  p.configure({ shape: 'circle', bc: 'free', mount: 'shaker', ma: 0.08, E: 69e9, rho: 2700, nu: 0.35, h: 0.0008, d: 0.24, Lx: 0.24, Ly: 0.24, zeta: 0.005, fmax: 5000, G: 25 });
  p.setDrive(0, 0);
  for (let f = 40; f < 5000; f *= 1.001) fs.push([f, p.field([{ f, F: 1 }], 9.81)]);
  for (let i = 10; i + 10 < fs.length; i++) { const y = fs[i][1]; if (y > fs[i - 1][1] && y >= fs[i + 1][1] && y > 1.2 * fs[i - 10][1] && y > 1.2 * fs[i + 10][1]) peak.push(fs[i][0]); }
  assert.equal(peak.filter(f => f < 4500).length, 6, 'resonance count');
  measured.forEach((f, i) => {
    const sim = peak[i];
    assert.ok(Math.abs(sim / f - 1) < 0.02, `s=${i + 1}: ${sim.toFixed(0)} Hz vs ${f}`); // +1.8, +1.1, +0.5, +0.3, +0.1, −0.4 %
    p.field([{ f: sim, F: 1 }], 9.81); // nodal circles: minima of |W(ρ)| along a radius (the response is axisymmetric)
    const W = new Float64Array(2 * NR);
    p.modes.forEach((m, k) => { if (m.n === 0) for (let j = 0; j < NR; j++) { W[2 * j] += p.cr[k] * m.R[j]; W[2 * j + 1] += p.ci[k] * m.R[j]; } });
    const A = j => Math.hypot(W[2 * j], W[2 * j + 1]), found = [];
    for (let j = NR - 2; j > 0; j--) if (A(j) < A(j - 1) && A(j) <= A(j + 1)) { const a = A(j - 1), b = A(j), c = A(j + 1); found.push((j + 0.5 * (a - c) / (a - 2 * b + c)) / (NR - 1)); }
    // circles beyond r/a = 0.2 agree within 0.022 (mean 0.008); the ones inside lie on the 8 mm drive bar, which a point
    // mount does not resolve, so they are not checked
    radii[i].forEach((r, k) => { if (r > 0.2) assert.ok(Math.abs(found[k] - r) < 0.025, `s=${i + 1} circle ${k + 1}: ${found[k]} vs ${r}`); });
  });
});

test('cantilever plate eigenvalues match Leissa NASA SP-160 Tables 4.51–4.52 (ν = 0.3)', () => {
  // λ = ωa²√(ρh/D), a = length from the clamped edge; Leissa's 1969 values are upper bounds from a coarser method
  const ref = { 1: [[3.472, 21.29, 27.2, 54.3, 61.3], [8.55, 31.1, 64.2, 71.1]], 0.8: [[3.479, 18.49, 22.2, 44.5, 61.1], [7.27, 28.0, 43.0, 66.6]],
    0.7: [[3.484, 15.48, 22.0, 40.0, 60.3], [6.63, 26.5, 34.1, 59.8]] };
  for (const [ab, [sym, anti]] of Object.entries(ref)) {
    const modes = rectModes(['clamped', 'free', 'free', 'free'], 0.3, 1 / ab).modes;
    [[sym, 0], [anti, 1]].forEach(([want, py]) => {
      const got = modes.filter(m => m.py === py && m.L > 1e-3).map(m => m.L);
      want.forEach((v, i) => assert.ok(got[i] <= v * 1.002 && got[i] >= v * 0.99, `a/b ${ab} ${py ? 'anti' : 'sym'}metric ${i + 1}: ${got[i]} vs ${v}`));
    });
  }
});

test('cantilever plates: measured frequencies (Leissa NASA SP-160 Tables 4.55–4.58)', () => {
  const C = ['clamped', 'free', 'free', 'free'], err = [];
  const model = (ab, nu, hOverA) => rectModes(C, nu, 1 / ab).modes.filter(m => m.L > 1e-3).map(m => {
    if (!hOverA) return m.L;
    const mat = { E: 1, rho: 1, nu, h: hOverA }, c = Math.sqrt(hOverA * hOverA / (12 * (1 - nu * nu))); // a = 1
    return mindlinOmega(m.L * c, mat) / c;
  });
  const inOrder = (meas, mod) => meas.forEach((v, i) => err.push(v / mod[i] - 1)); // the first modes, in order
  const nearest = (meas, mod) => meas.forEach(v => err.push(v / mod.reduce((b, x) => Math.abs(x - v) < Math.abs(b - v) ? x : b) - 1));
  // 4.55: aluminium (ν taken as 0.33), four aspect ratios
  for (const [ab, meas] of [[0.5, [3.34, 5.38, 10.31]], [1, [3.37, 8.26, 20.55, 27.15, 29.75]], [2, [3.36, 14.43, 20.86]], [5, [3.32, 20.84, 32.40]]]) inOrder(meas, model(ab, 0.33));
  inOrder([3.34, 8.23, 20.56, 25.97], model(1, 0.33)); // 4.57: 6061-T6, 7.5 × 7.5 × 0.125 in (Craig, Plass & Caughfield)
  inOrder([3.35, 8.53, 20.90, 26.72, 30.61], model(1, 0.29)); // 4.58: steel, 10 × 10 × 0.0627 in (Gustafson, Stokey & Zorowski)
  // 4.56: steel strips, a = 5 in, h = 0.1 in (Plunkett & Wilson), labelled modes up to λ ≈ 290; thick-plate corrected
  for (const [ab, meas] of [[2, [3.50, 21.7, 60.5, 118.7, 196.0, 292.0, 14.5, 48.1, 92.3, 154.0, 228.0, 92.8, 125.1, 176.0, 244.0, 246.0, 274.0]],
    [2.5, [3.50, 21.6, 60.4, 117.5, 17.3, 54.8, 101.5, 139.1]], [3.33, [3.50, 21.5, 59.8, 116.5, 195.0, 22.5, 69.6, 125.0, 187.0, 246.0]],
    [5, [3.45, 21.1, 59.3, 115.2, 190.0, 281.0, 32.0, 98.0, 169.0, 248.0]]]) nearest(meas, model(ab, 0.3, 0.02));
  assert.equal(err.length, 68);
  // real clamps are a little soft, so experiments sit up to a few % below ideal clamping (Leissa notes the same)
  err.forEach((e, i) => assert.ok(e > -0.06 && e < 0.03, `mode ${i}: ${(100 * e).toFixed(1)} %`));
  const mean = err.reduce((a, b) => a + b, 0) / err.length, median = err.map(Math.abs).sort((a, b) => a - b)[34];
  assert.ok(mean > -0.03 && mean < 0 && median < 0.025, `mean ${(100 * mean).toFixed(2)} %, median |e| ${(100 * median).toFixed(2)} %`);
});

test('findPeaks on a real analyser spectrum: frequency within 0.05 Hz, level within 0.05 dB', () => {
  // Web Audio AnalyserNode: Blackman window (α = 0.16), X[k] = (1/N) Σ x[n]w[n]e^{−2πikn/N}, dB = 20·log10|X|
  const N = 8192, sr = 44100, bin = sr / N;
  for (const f of [441.7, 1002.49, 6000.9]) { // a full-scale sine at awkward between-bin frequencies
    const db = new Float32Array(N / 2).fill(-200), k0 = Math.round(f / bin);
    for (let k = k0 - 12; k <= k0 + 12; k++) { // main lobe plus the first sidelobes
      let re = 0, im = 0;
      for (let n = 0; n < N; n++) {
        const x = Math.sin(2 * Math.PI * f * n / sr) * (0.42 - 0.5 * Math.cos(2 * Math.PI * n / N) + 0.08 * Math.cos(4 * Math.PI * n / N));
        re += x * Math.cos(2 * Math.PI * k * n / N); im -= x * Math.sin(2 * Math.PI * k * n / N);
      }
      db[k] = 20 * Math.log10(Math.hypot(re, im) / N);
    }
    const p = findPeaks(db, bin, 20, 20000, 4, -100);
    assert.equal(p.length, 1, 'one peak, no sidelobes');
    assert.ok(Math.abs(p[0].f - f) < 0.05, `frequency ${p[0].f} vs ${f}`);
    assert.ok(Math.abs(p[0].db - 20 * Math.log10(0.21)) < 0.05, `level ${p[0].db} dB`); // app.js maps 0.21 → 0 dBFS
  }
});
