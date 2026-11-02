import { Plate, findPeaks, nodalRadii } from './physics.js';

const $ = id => document.getElementById(id), num = id => +$(id).value;
const MATERIALS = { aluminum: [69, 2700, 0.33], steel: [200, 7850, 0.29], brass: [100, 8500, 0.34], glass: [70, 2500, 0.22], acrylic: [3.2, 1180, 0.37] };
const NOTES = ['C', 'C♯', 'D', 'D♯', 'E', 'F', 'F♯', 'G', 'G♯', 'A', 'A♯', 'B'];
const fmtHz = f => f >= 1000 ? `${(f / 1000).toFixed(2)} kHz` : `${f.toFixed(1)} Hz`;
const cents = c => `${c >= 0 ? '+' : '−'}${Math.abs(c).toFixed(0)}¢`;
const note = f => { const m = 69 + 12 * Math.log2(f / 440), r = Math.round(m); return `${NOTES[(r % 12 + 12) % 12]}${Math.floor(r / 12) - 1} ${cents((m - r) * 100)}`; };
const status = (msg, warn = false) => { $('status').textContent = msg; $('status').classList.toggle('warn', warn); };
const sliderHz = () => 20 * 1000 ** (num('tone') / 1000); // log slider, 20 Hz – 20 kHz
const form = $('controls');
form.reset(); // every page load starts from the defaults written in the HTML, whatever the browser restored
form.addEventListener('submit', e => e.preventDefault());

const outputs = () => {
  for (const el of form.querySelectorAll('input[type=range]')) {
    const out = el.nextElementSibling;
    if (out?.tagName === 'OUTPUT') out.value = el.id === 'grains' ? `${el.value / 1000}k` : el.value;
  }
};
form.addEventListener('input', e => { if (e.target.type === 'range') outputs(); });

// ---------- plate ----------
const plate = new Plate();
let rebuildTimer = 0;

function plateConfig() {
  const c = { shape: $('shape').value, bc: $('bc').value, edges: ['eL', 'eR', 'eT', 'eB'].map(id => $(id).value), E: num('E') * 1e9, rho: num('rho'), nu: num('nu'), h: num('h') / 1000,
    d: num('d') / 100, Lx: num('lx') / 100, Ly: num('ly') / 100, zeta: num('zeta'), fmax: num('fmax'), G: num('res'), mount: $('mount').value, ma: num('ma') / 1000, mindlin: $('mindlin').checked };
  return c.E > 0 && c.rho > 0 && c.nu >= 0 && c.nu < 0.5 && c.fmax > 0 && c.ma >= 0 ? c : null;
}

function rebuildNow() {
  clearTimeout(rebuildTimer);
  const c = plateConfig();
  if (!c) return status("Check the inputs: modulus and density must be positive, Poisson's ratio in [0, 0.5), moving mass ≥ 0.", true);
  const t = performance.now(), geom = `${plate.circle}|${plate.hx}|${plate.hy}`;
  document.body.dataset.shape = c.shape;
  document.body.dataset.mount = c.mount;
  $('forceUnit').textContent = c.mount === 'rigid' ? 'g' : 'N';
  plate.configure(c);
  setDrive();
  if (geom !== `${plate.circle}|${plate.hx}|${plate.hy}`) scatter();
  buildBase();
  listModes(performance.now() - t);
}

function rebuild() {
  clearTimeout(rebuildTimer);
  status('Solving plate modes…');
  rebuildTimer = setTimeout(rebuildNow, 80);
}

function listModes(ms) {
  const shown = plate.modes.filter(m => m.f > 0).slice(0, 150);
  $('modes').replaceChildren(...shown.map(m => {
    const b = document.createElement('button');
    b.type = 'button'; b.textContent = `${m.label} ${fmtHz(m.f)}`; b.onclick = () => playTone(m.f);
    if (m.R) b.onmouseenter = () => { b.title ||= `${m.f.toFixed(2)} Hz · nodal circles at r/a = ${nodalRadii(m).map(r => r.toFixed(4)).join(', ') || 'none'}`; }; // computed on hover
    return b;
  }));
  const c = plateConfig(), short = num('fmax') > plate.fTop * 1.01;
  // Mindlin's correction is exact for simply supported edges; for free/clamped ones it is within ~1–5 % at h/L = 0.1
  const thick = c.bc !== 'simply' && c.h / (c.shape === 'circle' ? c.d : Math.min(c.Lx, c.Ly)) > 0.05;
  status(`${plate.modes.length} modes, valid to ${fmtHz(plate.fTop)}${short ? ' — below "Highest"; a thicker or smaller plate covers more' : ''}`
    + `${thick ? ' · thick plate: free/clamped frequencies carry a few % model error' : ''} · solved in ${ms.toFixed(0)} ms`, short || thick);
}

function setDrive() { const r = num('dr'), a = num('da') * Math.PI / 180; plate.setDrive(r * Math.cos(a), r * Math.sin(a)); }

// ---------- audio ----------
let actx, analyser, out, buses, osc, mic, spec, source = null, toneHz = num('toneHz');

function audio() {
  if (actx) return;
  actx = new AudioContext();
  analyser = actx.createAnalyser();
  out = actx.createGain(); out.gain.value = $('hear').checked ? 1 : 0; out.connect(actx.destination);
  buses = {};
  for (const k of ['file', 'mic', 'tone']) {
    const g = buses[k] = actx.createGain();
    g.gain.value = 0; g.connect(analyser);
    if (k !== 'mic') g.connect(out); // never monitor the microphone, to avoid feedback
  }
  actx.createMediaElementSource($('player')).connect(buses.file);
  osc = actx.createOscillator(); osc.type = $('wave').value; osc.frequency.value = toneHz; osc.connect(buses.tone); osc.start();
  setAnalysis();
}

function setAnalysis() {
  if (!analyser) return;
  analyser.fftSize = num('fft'); analyser.smoothingTimeConstant = num('smooth');
  spec = new Float32Array(analyser.frequencyBinCount);
}

// source: 'file' | 'mic' | 'tone' | null; only that source's controls are shown
function use(src) {
  if (src || actx) { audio(); actx.resume(); }
  source = src;
  document.body.dataset.source = src ?? 'none';
  if (buses) for (const k in buses) buses[k].gain.value = k === src ? (k === 'tone' ? num('level') : 1) : 0;
  if (src !== 'file') $('player').pause();
  if (src !== 'mic' && mic) { mic.stream.getTracks().forEach(t => t.stop()); mic.node.disconnect(); mic = null; }
  for (const [id, k] of [['srcFile', 'file'], ['srcMic', 'mic'], ['srcTone', 'tone']]) $(id).classList.toggle('on', src === k);
}

function setTone(f) { // one frequency for the number box, the log slider and the oscillator
  toneHz = Math.min(20000, Math.max(1, f));
  $('toneHz').value = +toneHz.toFixed(2);
  $('tone').value = 1000 * Math.log(Math.max(20, toneHz) / 20) / Math.log(1000);
  if (osc) osc.frequency.setTargetAtTime(toneHz, actx.currentTime, 0.01);
}

function playTone(f) { setTone(f); use('tone'); }

$('file').addEventListener('change', () => {
  const f = $('file').files[0], p = $('player');
  if (!f) return;
  if (p.src.startsWith('blob:')) URL.revokeObjectURL(p.src);
  p.src = URL.createObjectURL(f);
  use('file'); p.play();
});
$('player').addEventListener('play', () => source !== 'file' && use('file'));
$('srcFile').addEventListener('click', () => {
  if (source === 'file') return use(null);
  use('file');
  if ($('player').src) $('player').play(); else $('file').click();
});
$('srcMic').addEventListener('click', async () => {
  if (source === 'mic') return use(null);
  use('mic');
  try {
    const stream = await navigator.mediaDevices.getUserMedia({ audio: { echoCancellation: false, noiseSuppression: false, autoGainControl: false } });
    if (source !== 'mic') return stream.getTracks().forEach(t => t.stop()); // switched away while asking
    const node = actx.createMediaStreamSource(stream);
    node.connect(buses.mic); mic = { stream, node };
  } catch (e) { use(null); status(`Microphone unavailable: ${e.message}`, true); }
});
$('micStop').addEventListener('click', () => use(null));
$('srcTone').addEventListener('click', () => use(source === 'tone' ? null : 'tone'));
$('tone').addEventListener('input', () => setTone(sliderHz()));
$('toneHz').addEventListener('input', () => { const f = num('toneHz'); if (f >= 1 && f <= 20000) setTone(f); });
$('toneHz').addEventListener('change', () => setTone(num('toneHz') || toneHz)); // tidy an out-of-range entry on blur
$('wave').addEventListener('input', () => { if (osc) osc.type = $('wave').value; });
$('freqSet').addEventListener('input', () => { if ($('freqSet').value) playTone(num('freqSet')); });
$('level').addEventListener('input', () => { if (source === 'tone') buses.tone.gain.value = num('level'); });
$('hear').addEventListener('input', () => { if (out) out.gain.value = $('hear').checked ? 1 : 0; });
for (const id of ['fft', 'smooth']) $(id).addEventListener('input', setAnalysis);

// ---------- sand ----------
let N = 0, px = new Float32Array(0), py = new Float32Array(0);

function scatter() {
  N = num('grains'); px = new Float32Array(N); py = new Float32Array(N);
  for (let p = 0; p < N; p++) {
    let u, v;
    do { u = Math.random() * 2 - 1; v = Math.random() * 2 - 1; } while (!plate.contains(u, v));
    px[p] = u; py[p] = v;
  }
}

// Grains hop only where the plate's acceleration exceeds gravity (Γ > 1). Each hop drifts down the gradient of Γ with a
// random scatter, so sand drains from antinodes and collects where the plate is quiet: the nodal lines.
// shortcut: grain transport is this phenomenological model, not ballistic per-cycle flight; the plate field itself is exact.
function stepSand(dt) {
  const G = plate.G, gam = plate.gamma, k = G / 2, top = G - 1.001, h = 1.5 / k, mob = num('mob') * dt, jit = num('jit');
  const circle = plate.circle, hx = plate.hx * 0.995, hy = plate.hy * 0.995;
  const sample = (x, y) => {
    let fx = (x + 1) * k - 0.5, fy = (y + 1) * k - 0.5;
    fx = fx < 0 ? 0 : fx > top ? top : fx; fy = fy < 0 ? 0 : fy > top ? top : fy;
    const i = fx | 0, j = fy | 0, ax = fx - i, ay = fy - j, o = j * G + i;
    return (gam[o] + ax * (gam[o + 1] - gam[o])) * (1 - ay) + (gam[o + G] + ax * (gam[o + G + 1] - gam[o + G])) * ay;
  };
  for (let p = 0; p < N; p++) {
    let x = px[p], y = py[p];
    const gm = sample(x, y);
    if (gm <= 1) continue;
    const e = Math.min(gm - 1, 4), gx = sample(x + h, y) - sample(x - h, y), gy = sample(x, y + h) - sample(x, y - h);
    const drift = mob * e / (Math.hypot(gx, gy) + 1e-9), kick = 2 * jit * Math.sqrt(e * dt);
    x += -gx * drift + kick * (Math.random() + Math.random() + Math.random() - 1.5);
    y += -gy * drift + kick * (Math.random() + Math.random() + Math.random() - 1.5);
    if (circle) { const r = Math.hypot(x, y); if (r > 0.995) { const s = Math.max(0, 1.99 - r) / r; x *= s; y *= s; } }
    else {
      if (x > hx) x = 2 * hx - x; else if (x < -hx) x = -2 * hx - x;
      if (y > hy) y = 2 * hy - y; else if (y < -hy) y = -2 * hy - y;
      x = Math.max(-hx, Math.min(hx, x)); y = Math.max(-hy, Math.min(hy, y));
    }
    px[p] = x; py[p] = y;
  }
}

// ---------- rendering ----------
const cv = $('plate'), c2 = cv.getContext('2d'), sc = $('spectrum'), s2 = sc.getContext('2d');
const rgb = hex => [1, 3, 5].map(i => parseInt(hex.slice(i, i + 2), 16));
const pack = ([r, g, b]) => (255 << 24 | b << 16 | g << 8 | r) >>> 0; // ImageData is RGBA bytes, little-endian
const mix = (a, b, t) => a.map((v, i) => Math.round(v + (b[i] - v) * t));
const MAGMA = [[0, 0, 4], [59, 15, 112], [140, 41, 129], [222, 73, 104], [254, 159, 109], [252, 253, 191]];
const mapLut = Uint32Array.from({ length: 512 }, (_, i) => { // dimmed half for Γ < 1, where sand can rest
  const t = (i & 255) / 255 * (MAGMA.length - 1), a = Math.floor(t), c = mix(MAGMA[a], MAGMA[Math.min(a + 1, MAGMA.length - 1)], t - a);
  return pack(i < 256 ? mix(c, [0, 0, 0], 0.55) : c);
});
let S = 0, img, buf, dens, base, onPlate, sandLut;

function resize() {
  const box = $('stage').getBoundingClientRect(), css = Math.max(120, Math.floor(Math.min(box.width, box.height) - 20));
  S = Math.min(1400, Math.round(css * devicePixelRatio));
  cv.width = cv.height = S; cv.style.width = cv.style.height = `${css}px`;
  img = c2.createImageData(S, S); buf = new Uint32Array(img.data.buffer); dens = new Uint8Array(S * S);
  buildBase();
}

function buildBase() {
  if (!S || !plate.G) return;
  const bg = rgb($('bg').value), pl = rgb($('plateCol').value), sand = rgb($('sandCol').value);
  base = new Uint32Array(S * S); onPlate = new Uint8Array(S * S);
  for (let y = 0; y < S; y++) for (let x = 0; x < S; x++) {
    const u = (x + 0.5) / S * 2 - 1, v = (y + 0.5) / S * 2 - 1, k = y * S + x;
    if (!plate.contains(u, v)) { base[k] = pack(bg); continue; }
    const edge = plate.circle ? 1 - Math.hypot(u, v) : Math.min(plate.hx - Math.abs(u), plate.hy - Math.abs(v));
    base[k] = pack(edge < 0.008 ? mix(pl, [255, 255, 255], 0.25) : mix(pl, [0, 0, 0], 0.15 * (u * u + v * v))); // rim, dish shading
    onPlate[k] = 1;
  }
  sandLut = Uint32Array.from({ length: 32 }, (_, d) => pack(mix(mix(pl, sand, 1 - 0.55 ** d), [255, 255, 255], Math.max(0, d - 6) / 80)));
}

function render() {
  dens.fill(0);
  const half = S / 2, sz = Math.max(1, Math.round(num('grain') * S / 600)), lim = S - sz; // grain size scales with the plate, not the screen
  for (let p = 0; p < N; p++) {
    const x = Math.min(lim, ((px[p] + 1) * half) | 0), y = Math.min(lim, ((py[p] + 1) * half) | 0);
    for (let dy = 0; dy < sz; dy++) for (let dx = 0; dx < sz; dx++) { const k = (y + dy) * S + x + dx; if (dens[k] < 255) dens[k]++; }
  }
  const ov = $('overlay').checked, gam = plate.gamma, top = Math.max(2, maxG), G = plate.G, cs = G / S, end = G - 1.001;
  for (let y = 0, k = 0; y < S; y++) {
    const fy = Math.min(end, Math.max(0, (y + 0.5) * cs - 0.5)), j = fy | 0, ay = fy - j;
    for (let x = 0; x < S; x++, k++) {
      const d = dens[k];
      if (d) buf[k] = sandLut[d > 31 ? 31 : d];
      else if (ov && onPlate[k]) { // bilinear Γ
        const fx = Math.min(end, Math.max(0, (x + 0.5) * cs - 0.5)), i = fx | 0, ax = fx - i, o = j * G + i;
        const g = (gam[o] + ax * (gam[o + 1] - gam[o])) * (1 - ay) + (gam[o + G] + ax * (gam[o + G + 1] - gam[o + G])) * ay;
        buf[k] = mapLut[(g < 1 ? 0 : 256) + Math.min(255, (g / top * 255) | 0)];
      } else buf[k] = base[k];
    }
  }
  c2.putImageData(img, 0, 0);
  if (!plate.circle) { // supports: clamped edges solid, simply supported dashed, free edges bare
    const r = devicePixelRatio, x0 = (1 - plate.hx) * half, x1 = (1 + plate.hx) * half, y0 = (1 - plate.hy) * half, y1 = (1 + plate.hy) * half;
    [[x0, y0, x0, y1], [x1, y0, x1, y1], [x0, y0, x1, y0], [x0, y1, x1, y1]].forEach(([a, b, c, d], i) => {
      const e = plate.edges[i];
      if (e === 'free') return;
      c2.strokeStyle = '#9aa3ad'; c2.lineWidth = (e === 'clamped' ? 5 : 2) * r; c2.setLineDash(e === 'simply' ? [6 * r, 5 * r] : []);
      c2.beginPath(); c2.moveTo(a, b); c2.lineTo(c, d); c2.stroke();
    });
    c2.setLineDash([]);
  }
  if ($('showDrive').checked) {
    const r = devicePixelRatio;
    c2.strokeStyle = '#ff6b5a'; c2.lineWidth = 2 * r;
    c2.beginPath(); c2.arc((plate.du + 1) * half, (plate.dv + 1) * half, 6 * r, 0, 2 * Math.PI); c2.stroke();
  }
}

function drawSpectrum() {
  const r = devicePixelRatio, W = Math.round(sc.clientWidth * r), H = Math.round(sc.clientHeight * r);
  if (sc.width !== W || sc.height !== H) { sc.width = W; sc.height = H; }
  s2.clearRect(0, 0, W, H);
  const f0 = Math.max(1, num('fmin')), f1 = Math.max(f0 * 1.01, num('fmax')), lr = Math.log(f1 / f0), lo = -120, hi = -10;
  const X = f => Math.log(f / f0) / lr * W, Y = db => H - 18 * r - (Math.max(lo, db) - lo) / (hi - lo) * (H - 30 * r);
  s2.lineWidth = r; s2.font = `${11 * r}px system-ui, sans-serif`;
  for (const m of plate.modes) { // the plate's eigenfrequencies, coloured by nodal diameters
    if (m.f < f0 || m.f > f1) continue;
    s2.strokeStyle = `hsl(${(m.n ?? 1) * 47 % 360} 70% 60% / 0.35)`;
    s2.beginPath(); s2.moveTo(X(m.f), 0); s2.lineTo(X(m.f), H - 18 * r); s2.stroke();
  }
  s2.fillStyle = '#8b919b';
  for (const f of [20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000]) if (f >= f0 && f <= f1) s2.fillText(f >= 1000 ? `${f / 1000}k` : `${f}`, X(f) + 2 * r, H - 5 * r);
  s2.strokeStyle = '#8b919b80'; s2.setLineDash([4 * r, 4 * r]);
  s2.beginPath(); s2.moveTo(0, Y(num('floor'))); s2.lineTo(W, Y(num('floor'))); s2.stroke(); s2.setLineDash([]);
  if (!source || !spec) return;
  const bin = actx.sampleRate / analyser.fftSize;
  s2.strokeStyle = $('sandCol').value; s2.beginPath();
  for (let x = 0; x < W; x++) {
    const i0 = Math.floor(f0 * Math.exp(lr * x / W) / bin), i1 = Math.max(i0, Math.floor(f0 * Math.exp(lr * (x + 1) / W) / bin));
    let m = -Infinity;
    for (let i = i0; i <= i1 && i < spec.length; i++) if (spec[i] > m) m = spec[i];
    if (x) s2.lineTo(x, Y(m)); else s2.moveTo(x, Y(m));
  }
  s2.stroke();
  peaks.forEach((p, i) => {
    s2.fillStyle = '#ff6b5a'; s2.beginPath(); s2.arc(X(p.f), Y(p.db), 3 * r, 0, 2 * Math.PI); s2.fill();
    if (i < 5) { s2.fillStyle = '#d9dce1'; s2.fillText(fmtHz(p.f), X(p.f) + 5 * r, Y(p.db) - 4 * r); }
  });
}

function readout() {
  const lines = [];
  if (!source) lines.push('Pick a sound: audio file, microphone or test tone');
  else if (!peaks.length) lines.push('No peaks above the floor — the sand is at rest');
  else {
    const p = peaks[0], m = plate.modes.reduce((b, q) => q.f > 0 && (!b || Math.abs(Math.log(q.f / p.f)) < Math.abs(Math.log(b.f / p.f))) ? q : b, null);
    lines.push(`strongest     ${fmtHz(p.f)}  ${note(p.f)}  ${p.db.toFixed(1)} dB`);
    if (m) lines.push(`nearest mode  ${m.label}  ${fmtHz(m.f)}  ${cents(1200 * Math.log2(p.f / m.f))}`);
  }
  lines.push(`peak plate acceleration ${maxG.toFixed(2)} g`);
  if (peaks.length) { // amplitude = a/ω² at the strongest peak; linear plate theory needs it well below the thickness
    const amp = maxG * num('g') / (2 * Math.PI * peaks[0].f) ** 2, ratio = amp / (num('h') / 1000);
    lines.push(`peak deflection ≈ ${(amp * 1e6).toPrecision(3)} µm = ${(100 * ratio).toFixed(1)} % of thickness${ratio > 0.2 ? '  ⚠ beyond linear theory' : ''}`);
  }
  lines.push(`drive ${gain.toPrecision(3)} ${plate.mount === 'rigid' ? 'g' : 'N'} at 0 dBFS`, `field ${fieldMs.toFixed(1)} ms · ${N / 1000}k grains`);
  $('readout').textContent = lines.join('\n');
}

// ---------- main loop ----------
let peaks = [], gain = num('force'), maxG = 0, fieldMs = 0, paused = false, last = performance.now(), frameNo = 0;

function frame(now) {
  const dt = Math.min(0.05, (now - last) / 1000);
  last = now;
  peaks = [];
  if (source && spec) {
    analyser.getFloatFrequencyData(spec);
    peaks = findPeaks(spec, actx.sampleRate / analyser.fftSize, num('fmin'), num('fmax'), num('peaks'), num('floor'));
  }
  // A full-scale sine reads −13.5 dB: the analyser's Blackman window has coherent gain 0.42, halved for one side.
  const t0 = performance.now(), raw = plate.field(peaks.map(p => ({ f: p.f, F: 10 ** (p.db / 20) / 0.21 })), num('g'));
  fieldMs += (performance.now() - t0 - fieldMs) * 0.05;
  if (!$('auto').checked) gain = num('force');
  else if (raw > 0) { const want = Math.min(1e4, num('target') / raw); gain = want < gain ? want : gain + (want - gain) * Math.min(1, 2 * dt); } // fast attack, slow release
  for (let i = 0; i < plate.gamma.length; i++) plate.gamma[i] *= gain; // Γ is linear in the force
  maxG = raw * gain;
  if (!paused) stepSand(dt);
  render();
  drawSpectrum();
  if (++frameNo % 6 === 0) readout();
  requestAnimationFrame(frame);
}

// ---------- wiring ----------
for (const id of ['shape', 'bc', 'E', 'rho', 'nu', 'd', 'lx', 'ly', 'h', 'zeta', 'res', 'fmax', 'mount', 'ma', 'mindlin', 'eL', 'eR', 'eT', 'eB']) $(id).addEventListener('input', rebuild);
for (const el of document.querySelectorAll('.edge')) { el.innerHTML = $('bc').innerHTML; el.value = $('bc').value; }
$('bc').addEventListener('input', () => { for (const el of document.querySelectorAll('.edge')) el.value = $('bc').value; }); // "Edge" sets all four
for (const id of ['E', 'rho', 'nu']) $(id).addEventListener('input', () => { $('material').value = 'custom'; });
$('material').addEventListener('input', () => {
  const m = MATERIALS[$('material').value];
  if (m) [$('E').value, $('rho').value, $('nu').value] = m;
  rebuild();
});
for (const id of ['dr', 'da']) $(id).addEventListener('input', setDrive);
for (const id of ['sandCol', 'plateCol', 'bg']) $(id).addEventListener('input', buildBase);
$('grains').addEventListener('change', scatter);
$('resprinkle').addEventListener('click', scatter);
$('pause').addEventListener('click', () => { paused = !paused; $('pause').textContent = paused ? 'Resume' : 'Pause'; $('pause').classList.toggle('on', paused); });
cv.addEventListener('click', e => {
  const b = cv.getBoundingClientRect(), u = (e.clientX - b.left) / b.width * 2 - 1, v = (e.clientY - b.top) / b.height * 2 - 1;
  $('dr').value = Math.min(0.97, Math.hypot(u, v)); $('da').value = Math.atan2(v, u) * 180 / Math.PI;
  for (const id of ['dr', 'da']) $(id).dispatchEvent(new Event('input', { bubbles: true }));
});
new ResizeObserver(resize).observe($('stage'));

// ---------- reset and presets ----------
// Each preset starts from the defaults, applies its settings in order ('bc' before single edges), clears the sand,
// and plays a tone at a named eigenmode or a fixed frequency.
const PRESETS = [
  { name: 'Eight-point star', note: 'free circle, (8,1) mode', set: { dr: 0.6, da: -90, target: 12 }, mode: '(8,1)' },
  { name: 'Nested rings', note: 'circle bolted to a shaker', set: { mount: 'shaker', ma: 80, E: 69, nu: 0.35, h: 0.8, target: 15 }, f: 1700 },
  { name: 'Chladni square', note: 'free square, centre drive', set: { shape: 'rect', target: 12 }, mode: '#38' },
  { name: 'Shu et al. 1247 Hz', note: 'bolted square, as photographed', set: { shape: 'rect', mount: 'shaker', ma: 20, E: 70, lx: 32, ly: 32, fmax: 5000, target: 15 }, f: 1247 },
  // fitted to a widely shared 12-pattern photo (345–5907 Hz): only h/L² and the mass ratio are determined, ±1 %
  { name: '12-pattern photo plate', note: '24 cm × 0.75 mm Al, 7 g shaker', set: { shape: 'rect', lx: 24, ly: 24, h: 0.75, mount: 'shaker', ma: 7, target: 12 }, f: 2041 },
  { name: 'Cantilever cross', note: 'top edge clamped', set: { shape: 'rect', eT: 'clamped', lx: 19, ly: 19, h: 3.15, dr: 0.75, da: 40, fmax: 3000, target: 12 }, mode: '#5' },
  { name: 'Clamped drum', note: 'clamped circle, (2,1) mode', set: { bc: 'clamped', dr: 0.3, da: -90, target: 12 }, mode: '(2,1)' },
  { name: 'Square-wave harmonics', note: 'free circle, 220 Hz', set: { dr: 0.45, da: -90, target: 10 }, f: 220, wave: 'square' },
  { name: 'Grid of cells', note: 'simply supported 2:1 rectangle', set: { shape: 'rect', bc: 'simply', lx: 40, ly: 20, dr: 0.3, da: 30, target: 12 }, mode: '#9' },
];

function syncForm() { // derived state after the whole form changed (page load, reset, preset)
  outputs();
  rebuildNow();
  setTone(num('toneHz'));
  if (osc) osc.type = $('wave').value;
  if (out) out.gain.value = $('hear').checked ? 1 : 0;
  setAnalysis();
  scatter();
  gain = num('force');
}

function applyPreset(p) {
  form.reset();
  for (const [id, v] of Object.entries(p.set)) {
    if ($(id).type === 'checkbox') $(id).checked = v; else $(id).value = v;
    if (id === 'bc') for (const el of document.querySelectorAll('.edge')) el.value = v;
  }
  if (p.wave) $('wave').value = p.wave;
  syncForm();
  const f = p.f ?? plate.modes.find(m => m.label === p.mode)?.f;
  if (f) playTone(f); else status(`Preset mode ${p.mode} is not in this plate's mode list`, true);
}

$('presets').replaceChildren(...PRESETS.map(p => {
  const b = document.createElement('button'), small = document.createElement('small');
  b.type = 'button'; small.textContent = p.note; b.append(p.name, small); b.onclick = () => applyPreset(p);
  return b;
}));
$('resetAll').addEventListener('click', () => { form.reset(); syncForm(); use(null); });

syncForm();
resize();
requestAnimationFrame(frame);
