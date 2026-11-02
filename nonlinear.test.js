// Checks of the large-amplitude module against sources and independent code. "scipy" below means scripts that share no code with
// nonlinear.js: the Airy problem solved in strong form by exact integrals (circles), or a displacement Ritz for (u, v), which bounds
// the membrane energy from above while this module bounds it from below (rectangles).
import test from 'node:test';
import assert from 'node:assert/strict';
import { circleModes, rectModes, Plate, NR } from './physics.js';
import { amplitudeLimit, circleKappas, duffingShift, pinnedModes, rectKappas, regime } from './nonlinear.js';
import { findMode, sweepResonances } from './testkit.js';

const close = (got, want, rel, msg) => assert.ok(Math.abs(got / want - 1) <= rel, `${msg}: got ${got}, want ${want}`);

// a one-term trial shape in the form of a circleModes mode
const radialTrial = (shape, L) => {
  const R = Float32Array.from({ length: NR }, (_, j) => shape(j / (NR - 1)));
  let norm = 0;
  for (let j = 0; j < NR; j++) norm += (j === 0 || j === NR - 1 ? 1 : j % 2 ? 4 : 2) * R[j] * R[j] * j / (NR - 1);
  return { n: 0, m: 0, L, R, Rmax: Math.max(...R.map(Math.abs)), inv: 1 / (2 * Math.PI * norm / (3 * (NR - 1))) };
};

test('circles reproduce Yamaki (1961) as tabulated in Leissa SP-160 Table 12.4, with his one-term shapes', () => {
  const nu = 0.3;
  const c2 = (1 + nu) / (5 + nu), c1 = -(6 + 2 * nu) / (5 + nu);
  const shapes = { simply: r => 1 + c1 * r * r + c2 * r ** 4, clamped: r => (1 - r * r) ** 2 };
  const rows = [ // edge, in-plane, ω²ρa⁴/Eh², β
    ['simply', 'free', 2.242, 0.591],
    ['simply', 'fixed', 2.242, 4.148],
    ['clamped', 'free', 9.768, 1.429],
    ['clamped', 'fixed', 9.768, 4.602],
  ];
  for (const [bc, inPlane, alpha, beta] of rows) {
    const trial = radialTrial(shapes[bc], Math.sqrt(alpha * 12 * (1 - nu * nu)));
    close(circleKappas([trial], bc, nu, inPlane)[0] / trial.Rmax ** 2, beta / alpha, 5e-3, `${bc} ${inPlane}`);
  }
});

test('free-edge and clamped circles match the Γ tables of the VK-Gong manual (Touzé group), ν = 0.38', () => {
  const nu = 0.38;
  const gamma = (modes, kappa, n, s) => { // Γ of unit-norm modes, q̈ + ω²q + 12(1−ν²)Γq³ = 0 with w/h = qΦ
    const k = modes.indexOf(findMode(modes, n, s));
    return kappa[k] * modes[k].L ** 2 * modes[k].inv / (12 * (1 - nu * nu));
  };
  const free = circleModes('free', nu, 12);
  const kf = circleKappas(free, 'free', nu);
  close(gamma(free, kf, 2, 0), 1.8983, 1e-4, 'free (2,0)');
  close(gamma(free, kf, 0, 1), 8.5747, 1e-4, 'free (0,1)');
  const clamped = circleModes('clamped', nu, 12);
  const kc = circleKappas(clamped, 'clamped', nu, 'fixed');
  close(gamma(clamped, kc, 0, 0), 9.0374, 1e-4, 'clamped (0,0)');
  // their 64.413 for (1,0) sits 0.7 % above both this and the scipy solution of their own eq. 3.9, 63.943
  close(gamma(clamped, kc, 1, 0), 63.943, 1e-4, 'clamped (1,0)');
});

test('circles match an independent strong-form solution of the Airy equation (scipy)', () => {
  const rows = [ // ν, edge, in-plane, n, s, κ
    [0.33, 'free', 'free', [[2, 0, 0.17180], [0, 1, 0.72111], [3, 0, 0.18451], [4, 0, 0.17873], [5, 0, 0.16980], [6, 0, 0.16080]]],
    [0.38, 'free', 'free', [[1, 1, 0.16465], [0, 2, 0.36100]]],
    [0.38, 'clamped', 'fixed', [[2, 0, 0.22176], [0, 1, 0.62972]]],
    [0.38, 'clamped', 'free', [[0, 0, 0.18119], [1, 0, 0.08255]]],
    [0.3, 'simply', 'fixed', [[0, 0, 1.67790], [1, 0, 0.54775], [2, 0, 0.33854]]],
    [0.3, 'simply', 'free', [[0, 0, 0.25345], [1, 0, 0.07939]]],
    [0.3, 'clamped', 'fixed', [[3, 0, 0.19080]]],
  ];
  for (const [nu, bc, inPlane, cases] of rows) {
    const modes = circleModes(bc, nu, 14);
    const kappa = circleKappas(modes, bc, nu, inPlane);
    for (const [n, s, want] of cases) close(kappa[modes.indexOf(findMode(modes, n, s))], want, 5e-4, `${bc} ${inPlane} (${n},${s})`);
  }
});

test('square plates match Yamaki (1961) via Leissa SP-160 Sec. 12.2 and the closed form', () => {
  const nu = 0.3;
  const ritz = rectModes('simply', nu, 1);
  const kappaA = inPlane => {
    const { kappa } = rectKappas(ritz, 'simply', nu, 1, inPlane);
    return kappa[0] / 4; // the fundamental is 2·sinπξ·sinπη, peak 2
  };
  const prefactor = 3 * (1 - nu * nu);
  close(kappaA('free'), prefactor * 0.06492, 1e-3, 'stress-free, case (a)');
  close(kappaA('straight'), prefactor / 8, 1e-5, 'straight edges, case (c)');
  // Chu & Herrmann's hinged plate, case (b). The scan of Leissa reads (3−ν)/8(1+ν), which cannot be: the clamped case has 1/(1−ν),
  // the closed form gives 1/8 + 1/(4(1−ν)), and the scipy displacement Ritz gives 5.2650/4.
  close(kappaA('sliding'), prefactor * (3 - nu) / (8 * (1 - nu)), 1e-5, 'hinged, case (b)');
  close(kappaA('sliding'), 5.2650 / 4, 1e-4, 'hinged, scipy');

  // Yamaki's clamped trial shape sin²πξ·sin²πη: one constant, fixed by his normalisation, links all three published cases.
  const scale = Math.sqrt(64 / 9), at = x => [Math.sin(Math.PI * x) ** 2];
  const trial = { modes: [{ L: 37.22060906, vec: Float64Array.of(scale) }], nbx: 1, nby: 1, atX: at, atY: at };
  const clamped = inPlane => {
    const { kappa, peak } = rectKappas(trial, 'clamped', nu, 1, inPlane, 60);
    return kappa[0] / peak[0] ** 2;
  };
  const constant = clamped('straight') / 0.16656;
  close(clamped('sliding') / constant, 0.16656 + 0.14063 / (1 - nu), 1e-4, 'clamped case (b)');
  close(clamped('free') / constant, 0.14903, 6e-3, 'clamped case (a), a lower bound on his series value');
});

test('every simply supported mode matches the closed form, up to mode 100 and beyond', () => {
  const nu = 0.3, beta = 1.3;
  const ritz = rectModes('simply', nu, beta);
  const exact = (L, inPlane) => { // from the sine modes: ∇⁴Ψ has only cos 2mπξ and cos 2nπη, plus the mean strain if 'sliding'
    let best = [1, 1];
    for (let m = 1; m < 30; m++) {
      for (let n = 1; n < 30; n++) {
        if (Math.abs(Math.PI ** 2 * (m * m + n * n / beta ** 2) - L) < Math.abs(Math.PI ** 2 * (best[0] ** 2 + best[1] ** 2 / beta ** 2) - L)) best = [m, n];
      }
    }
    const a = (best[0] * Math.PI) ** 2, b = (best[1] * Math.PI / beta) ** 2;
    const move = 3 * (1 - nu * nu) * (a * a + b * b) / (4 * (a + b) ** 2);
    const mean = inPlane === 'sliding' ? 3 * (a * a + b * b + 2 * nu * a * b) / (2 * (a + b) ** 2) : 0;
    return 4 * (move + mean);
  };
  for (const inPlane of ['straight', 'sliding']) {
    const { kappa } = rectKappas(ritz, 'simply', nu, beta, inPlane);
    ritz.modes.slice(0, 100).forEach((mode, k) => close(kappa[k], exact(mode.L, inPlane), 1e-3, `${inPlane} #${k + 1}`));
    // the Airy basis is a lower bound and it does run out at the very top, far above any frequency that reaches large amplitude
    ritz.modes.forEach((mode, k) => assert.ok(kappa[k] <= exact(mode.L, inPlane) * (1 + 1e-6), `${inPlane} lower bound #${k + 1}`));
  }
});

test('rectangles are bracketed by an independent displacement Ritz upper bound (scipy)', () => {
  const nu = 0.3;
  const rows = [ // edge, β, in-plane, [mode index among elastic modes, upper bound on κ]
    ['simply', 1.3, 'free', [[0, 0.59833], [1, 0.83944], [2, 0.39241]]],
    ['simply', 1, 'fixed', [[0, 5.29179]]],
    ['clamped', 1, 'free', [[0, 0.95333], [3, 1.27868], [4, 1.41692]]],
    ['clamped', 1, 'fixed', [[0, 2.92590], [3, 4.07842], [4, 5.10395]]],
    ['clamped', 1, 'sliding', [[0, 2.87264]]],
    ['clamped', 1.5, 'fixed', [[0, 2.95545], [1, 3.86997], [4, 4.31007]]],
    ['free', 1, 'free', [[0, 1.12510], [1, 2.79180], [2, 1.54700]]],
    ['free', 1.5, 'free', [[0, 0.92482], [1, 0.25744], [2, 1.48453], [3, 0.06068], [4, 2.23470]]],
  ];
  for (const [bc, beta, inPlane, cases] of rows) {
    const ritz = rectModes(bc, nu, beta);
    const { kappa } = rectKappas(ritz, bc, nu, beta, inPlane);
    const elastic = ritz.modes.map((m, k) => k).filter(k => ritz.modes[k].L > 1e-3);
    for (const [index, upper] of cases) {
      const got = kappa[elastic[index]];
      assert.ok(got <= upper * (1 + 1e-4), `${bc} ${beta} ${inPlane} #${index}: ${got} above the upper bound ${upper}`);
      close(got, upper, 2e-3, `${bc} ${beta} ${inPlane} #${index}`);
    }
  }
});

// Forced response by RK4, the reference for harmonic balance: the amplitude at the drive frequency after settling.
function forcedResponse(state, { w, zeta, gamma, f }, cycles) {
  const steps = 96, dt = 2 * Math.PI / w / steps;
  const accel = (q, v, t) => f * Math.cos(w * t + state.phase) - 2 * zeta * v - q - gamma * q * q * q;
  let { q, v } = state, t = 0, a = 0, b = 0, n = 0;
  for (let c = 0; c < cycles; c++) {
    for (let i = 0; i < steps; i++) {
      const k1 = accel(q, v, t);
      const k2 = accel(q + 0.5 * dt * v, v + 0.5 * dt * k1, t + 0.5 * dt);
      const k3 = accel(q + 0.5 * dt * v + 0.25 * dt * dt * k1, v + 0.5 * dt * k2, t + 0.5 * dt);
      const k4 = accel(q + dt * v + 0.5 * dt * dt * k2, v + dt * k3, t + dt);
      q += dt * v + dt * dt * (k1 + k2 + k3) / 6;
      v += dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
      t += dt;
      if (c >= cycles - 8) {
        a += q * Math.cos(w * t + state.phase);
        b += q * Math.sin(w * t + state.phase);
        n++;
      }
    }
  }
  state.q = q;
  state.v = v;
  state.phase = (state.phase + w * t) % (2 * Math.PI); // the drive stays continuous across a frequency step, as in a real sweep
  return 2 * Math.hypot(a, b) / n;
}

test('harmonic balance follows a time integration through a sweep, in both directions, including the hysteresis', () => {
  const zeta = 0.01, gamma = 1, f = 0.06;
  const sweep = frequencies => {
    const state = { q: 0, v: 0, phase: 0 };
    let z = 0;
    return frequencies.map(w => {
      z = duffingShift(1 - w * w, 2 * zeta * w, 0.75 * gamma * f * f, z);
      return [w, Math.sqrt(z / (0.75 * gamma)), forcedResponse(state, { w, zeta, gamma, f }, 300)];
    });
  };
  const up = Array.from({ length: 27 }, (_, i) => 0.95 + 0.02 * i);
  const upward = sweep(up), downward = sweep(up.slice().reverse());
  for (const [w, balance, integrated] of [...upward, ...downward]) close(balance, integrated, 0.015, `ω = ${w}`);
  const at = (rows, w) => rows.find(row => Math.abs(row[0] - w) < 1e-9);
  const middle = up[16];
  assert.ok(at(upward, middle)[1] > 5 * at(downward, middle)[1], 'the two sweeps sit on different branches');
  assert.ok(at(upward, middle)[2] > 5 * at(downward, middle)[2], 'so does the integration');
});

test('the fold appears at the shift predicted by amplitudeLimit, and regime() names it', () => {
  const zeta = 0.005, kappaA = 0.5;
  const hysteresis = force => { // largest gap between the two branches over a fine sweep, by the module's own branch tracking
    let up = 0, down = 0, gap = 0;
    const b = w => 2 * zeta * w;
    const list = Array.from({ length: 4000 }, (_, i) => 0.99 + i * 2e-5);
    const uppers = list.map(w => (up = duffingShift(1 - w * w, b(w), 0.75 * kappaA * force * force, up)));
    list.slice().reverse().forEach((w, i) => {
      down = duffingShift(1 - w * w, b(w), 0.75 * kappaA * force * force, down);
      gap = Math.max(gap, Math.abs(down - uppers[uppers.length - 1 - i]));
    });
    return gap;
  };
  // the cusp of the cubic: Ω = −√3b, z = 2b/√3, so f = x·|Ω + z + ib| = x·2b/√3 with b = 2ζ
  const forceAtCusp = amplitudeLimit(2 / Math.sqrt(3), kappaA, zeta) * 4 * zeta / Math.sqrt(3);
  assert.ok(hysteresis(forceAtCusp * 0.97) < 1e-4, 'no jump below the cusp');
  assert.ok(hysteresis(forceAtCusp * 1.1) > 1e-3, 'a jump above it');
  assert.equal(regime(1), 'corrected');
  assert.equal(regime(0.05), 'negligible');
  assert.equal(regime(2), 'caution');
  assert.equal(regime(4), 'broken');
});

test('the unforced companion of a degenerate pair goes unstable where the Mathieu analysis says', () => {
  // q̈₁ + 2ζq̇₁ + q₁ + Γq₁(q₁²+q₂²) = f cos ωt, the same for q₂ with no force (rotation invariance gives the one Γ; VK-Gong's table
  // confirms it with Γ¹₁₂₂₁ = Γ¹₁₁₁₁). On the driven branch S = 1 − ω² + ¾Γa², and q₂ grows where S² + b² < ½SΓa².
  const zeta = 0.005, Gamma = 1;
  const growth = (a, S) => {
    const w = Math.sqrt(1 + 0.75 * Gamma * a * a - S), b = 2 * zeta * w;
    const f = a * Math.hypot(S, b);
    const dt = 2 * Math.PI / w / 64;
    const rate = (y, t) => {
      const e = y[0] * y[0] + y[2] * y[2];
      return [y[1], f * Math.cos(w * t) - 2 * zeta * y[1] - y[0] - Gamma * y[0] * e, y[3], -2 * zeta * y[3] - y[2] - Gamma * y[2] * e];
    };
    // start on the driven solution, q₁ = a cos(ωt − φ), with a tiny companion
    const phase = Math.atan2(b, S);
    let y = [a * Math.cos(phase), a * w * Math.sin(phase), 1e-6, 0], t = 0, top = 0;
    const periods = 1500;
    for (let i = 0; i < periods * 64; i++) {
      const k1 = rate(y, t);
      const k2 = rate(y.map((v, j) => v + 0.5 * dt * k1[j]), t + 0.5 * dt);
      const k3 = rate(y.map((v, j) => v + 0.5 * dt * k2[j]), t + 0.5 * dt);
      const k4 = rate(y.map((v, j) => v + dt * k3[j]), t + dt);
      y = y.map((v, j) => v + dt * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6);
      t += dt;
      if (i > (periods - 20) * 64) top = Math.max(top, Math.abs(y[2]));
    }
    return top;
  };
  const b0 = 2 * zeta;
  // inside the predicted tongue at S = b: threshold Γa² = 4b
  const edge = Math.sqrt(4 * b0 / Gamma);
  assert.ok(growth(edge * 1.25, b0 * 1.05) > 1e-3, 'above the threshold the companion grows');
  assert.ok(growth(edge * 0.8, b0 * 1.05) < 1e-5, 'below the threshold it decays');
});

test('a bolted plate has its own modes, and the free modes’ κ alone misses their stiffening', () => {
  const nu = 0.33;
  const all = circleModes('free', nu, 40).filter(mode => mode.n === 0);
  const pinned = pinnedModes(all.map(mode => mode.L), all.map(mode => mode.R[0] * Math.sqrt(mode.inv)), 3);
  const kappas = circleKappas(all, 'free', nu);

  // the same plate through Plate's rigid mount: resonances where the bolt's receptance vanishes
  const plate = new Plate();
  plate.configure({ E: 69e9, rho: 2700, nu, h: 0.001, shape: 'circle', bc: 'free', d: 0.24, Lx: 0.24, Ly: 0.24, zeta: 0.002, G: 21, fmax: 40000, mount: 'rigid' });
  plate.setDrive(0, 0);
  const found = sweepResonances(plate, 20, 2500);

  pinned.forEach(({ L, c }, index) => {
    const f = L * plate.wScale / (2 * Math.PI);
    if (f < 2400) close(found.reduce((p, x) => (Math.abs(x - f) < Math.abs(p - f) ? x : p)), f, 0.01, `pinned mode ${index}`);
    const R = new Float32Array(NR);
    all.forEach((mode, k) => mode.R.forEach((value, j) => { R[j] += c[k] * Math.sqrt(mode.inv) * value; }));
    const shape = radialTrial(r => R[Math.round(r * (NR - 1))], L);
    const own = circleKappas([shape], 'free', nu)[0] * L * L / shape.inv; // γ per ω² scale, unit norm
    const fromFree = all.reduce((sum, mode, k) => sum + (mode.L ? kappas[k] * mode.L ** 2 / mode.inv * c[k] ** 4 : 0), 0);
    assert.ok(fromFree < 0.5 * own, `pinned ${index}: the sum over free modes gives ${fromFree}, the shape itself ${own}`);
  });
});
