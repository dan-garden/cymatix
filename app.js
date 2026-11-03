import { DRIVER_LIMIT, clamp, findPeaks } from './physics.js';
import { PlateModel, stiffeningRegime } from './plateModel.js';
import { AIR, VACUUM, WATER } from './fluid.js';
import { mindlinNodalRadii } from './mindlin.js';
import { AIR as GRAIN_AIR, SPECIES, cellCapacity, grainMass, newFields, stepGrains, transportFields, updateCrowding } from './grains.js';
import { createSandGL } from './sandGL.js';
import { createOrbitCamera, orbitControls, pickPlane } from './orbitCamera.js';
import { createPlateView3D } from './plateView3D.js';
import { traceFieldLines } from './fieldLines.js';
import { autoExaggeration, heapHeightPerGrain, peakDeflection, plateHeight, slowMotion, strongestTones } from './plateSurface.js';
import { runPlan, solvePlan } from './solvePlan.js';
import { MATERIALS } from './materials.js';
import { cancelGroup, runInGroup, runTask } from './solverPool.js';

// ---------- helpers ----------

const $ = id => document.getElementById(id);
const num = id => +$(id).value;
const listen = (ids, event, handler) => {
  for (const id of [ids].flat()) $(id).addEventListener(event, handler);
};
const makeButton = (onClick, ...content) => {
  const button = document.createElement('button');
  button.type = 'button';
  button.append(...content);
  button.onclick = onClick;
  return button;
};
const strokeLine = (ctx, x0, y0, x1, y1) => {
  ctx.beginPath();
  ctx.moveTo(x0, y0);
  ctx.lineTo(x1, y1);
  ctx.stroke();
};

const NOTES = ['C', 'C♯', 'D', 'D♯', 'E', 'F', 'F♯', 'G', 'G♯', 'A', 'A♯', 'B'];
const fmtHz = f => (f >= 1000 ? `${(f / 1000).toFixed(2)} kHz` : `${f.toFixed(1)} Hz`);
const fmtCents = c => `${c >= 0 ? '+' : '−'}${Math.abs(c).toFixed(0)}¢`;
const fmtNote = f => {
  const midi = 69 + 12 * Math.log2(f / 440);
  const nearest = Math.round(midi);
  return `${NOTES[(nearest % 12 + 12) % 12]}${Math.floor(nearest / 12) - 1} ${fmtCents((midi - nearest) * 100)}`;
};

const showStatus = (message, warn = false) => {
  $('status').textContent = message;
  $('status').classList.toggle('warn', warn);
};

// ---------- form ----------

const form = $('controls');
form.reset(); // every page load starts from the defaults written in the HTML, whatever the browser restored
form.addEventListener('submit', e => e.preventDefault());

const updateRangeOutputs = () => {
  for (const el of form.querySelectorAll('input[type=range]')) {
    const output = el.nextElementSibling;
    if (output?.tagName === 'OUTPUT') output.value = el.id === 'grains' ? `${el.value / 1000}k` : el.value;
  }
};
const syncLineControls = () => { $('heightScale').disabled = !$('fieldLines').checked; }; // Waves only moves the lines
form.addEventListener('input', e => {
  if (e.target.type === 'range') updateRangeOutputs();
});

// ---------- plate ----------

const FLUIDS = { air: AIR, vacuum: VACUUM, water: WATER };
const plate = new PlateModel();
plate.keepDrives = true;
let rebuildTimer = 0;
let activeConfig = null; // the last valid plate settings

function plateConfig() {
  const config = {
    shape: $('shape').value,
    bc: $('bc').value,
    edges: ['eL', 'eR', 'eT', 'eB'].map(id => $(id).value),
    E: num('E') * 1e9,
    rho: num('rho'),
    nu: num('nu'),
    h: num('h') / 1000,
    d: num('d') / 100,
    Lx: num('lx') / 100,
    Ly: num('ly') / 100,
    zeta: num('zeta'),
    fmax: num('fmax'),
    G: num('res'),
    mount: $('mount').value,
    ma: num('ma') / 1000,
    mindlin: $('mindlin').checked,
    fluid: FLUIDS[$('fluid').value],
    clampRadius: num('clamp') / 1000,
    stiffening: $('stiffening').checked,
    modalZeta: modalDampingFor(),
  };
  const { E, rho, nu, fmax, ma } = config;
  return E > 0 && rho > 0 && nu >= 0 && nu < 0.5 && fmax > 0 && ma >= 0 ? config : null;
}

const plateGeometry = () => `${plate.circle}|${plate.hx}|${plate.hy}`;

function rebuildNow() {
  clearTimeout(rebuildTimer);
  const config = plateConfig();
  if (!config) {
    showStatus("Check the inputs: modulus and density must be positive, Poisson's ratio in [0, 0.5), moving mass ≥ 0.", true);
    return;
  }
  activeConfig = config;
  const startedAt = performance.now();
  const previousGeometry = plateGeometry();
  document.body.dataset.shape = config.shape;
  document.body.dataset.mount = config.mount;
  $('forceUnit').textContent = config.mount === 'rigid' ? 'g' : 'N';
  requestPhysicalDamping(config);
  const plan = solvePlan(config);
  const solved = plan && solvedPlates.get(plan.id);
  if (!plan || solved) {
    cancelGroup('plate');
    activeSolveId = null;
  }
  plate.configure({ ...config, solved });
  if (plan && !solved) startSolve(plan);
  applyDriver();
  if (previousGeometry !== plateGeometry()) scatterSand();
  buildBase();
  listModes(performance.now() - startedAt);
}

// ---------- physical damping, worked out in a worker thread ----------

let dampingCache = { key: null, zetas: null };

const dampingKey = ({ shape, bc, edges, E, rho, nu, h, d, Lx, Ly, fmax, fluid }) =>
  JSON.stringify([shape, bc, edges, E, rho, nu, h, d, Lx, Ly, fmax, fluid.rho]);

function modalDampingFor() {
  return $('dampingModel').value === 'physical' ? dampingCache.zetas : null;
}

function requestPhysicalDamping(config) {
  if ($('dampingModel').value !== 'physical') return;
  const key = dampingKey(config);
  if (dampingCache.key === key) return;
  dampingCache = { key, zetas: null };
  const { modalZeta, ...request } = config;
  runTask('physicalDamping', request).then(zetas => {
    if (dampingCache.key !== key) return; // the plate changed while this was solving
    dampingCache.zetas = zetas;
    rebuildNow();
  });
}

// ---------- better and larger solves in worker threads ----------

const solvedPlates = new Map(); // the last few, so moving a slider back is instant
const SOLVED_KEPT = 3;
let activeSolveId = null;

const runPlateSolve = (task, ...args) => runInGroup('plate', task, ...args);

async function startSolve(plan) {
  if (activeSolveId === plan.id) return;
  cancelGroup('plate'); // a solve for a plate that has since changed only holds the threads up
  activeSolveId = plan.id;
  const solved = await runPlan(plan, runPlateSolve);
  solvedPlates.set(plan.id, solved);
  if (solvedPlates.size > SOLVED_KEPT) solvedPlates.delete(solvedPlates.keys().next().value);
  if (activeSolveId !== plan.id) return; // the plate changed again while this one solved
  activeSolveId = null;
  rebuildNow();
}

function scheduleRebuild() {
  clearTimeout(rebuildTimer);
  showStatus('Solving plate modes…');
  rebuildTimer = setTimeout(rebuildNow, 80);
}

const nodalTooltip = mode => {
  const radii = mindlinNodalRadii(mode).map(r => r.toFixed(4)).join(', ') || 'none';
  return `${mode.f.toFixed(2)} Hz · nodal circles at r/a = ${radii}`;
};

function modeButton(mode) {
  const button = makeButton(() => playTone(mode.f), `${mode.label} ${fmtHz(mode.f)}`);
  if (mode.R) {
    button.onmouseenter = () => { button.title ||= nodalTooltip(mode); }; // computed on hover
  }
  return button;
}

function listModes(solveMs) {
  const shown = plate.modes.filter(mode => mode.f > 0).slice(0, 150);
  $('modes').replaceChildren(...shown.map(modeButton));

  const config = plateConfig();
  const tooShort = num('fmax') > plate.fTop * 1.01;
  // Mindlin's correction is exact for simply supported edges; for free/clamped ones it is within ~1–5 % at h/L = 0.1
  const span = config.shape === 'circle' ? config.d : Math.min(config.Lx, config.Ly);
  const tooThick = config.bc !== 'simply' && config.h / span > 0.05;
  const shortNote = !tooShort ? ''
    : activeSolveId ? ' — solving more modes to reach "Highest"…'
    : ' — below "Highest"; a thicker or smaller plate covers more';
  const thickNote = tooThick && !plate.solvedModes[0]?.B ? ' · thick plate: solving true Mindlin modes…' : '';
  const waterNote = $('fluid').value === 'water' ? ' · sand transport is modelled in air only' : '';
  showStatus(
    `${plate.modes.length} modes, valid to ${fmtHz(plate.fTop)}${shortNote}${thickNote}${waterNote} · solved in ${solveMs.toFixed(0)} ms`,
    tooShort,
  );
}

function showMaterialSource() {
  const material = MATERIALS[$('material').value];
  const hint = $('materialSource');
  hint.replaceChildren();
  if (!material) return;
  const link = document.createElement('a');
  link.href = material.url;
  link.target = '_blank';
  link.rel = 'noopener';
  link.textContent = material.source;
  hint.append('Values from ', link, '. A real sheet can differ by a few percent.');
}

function applyDriver() {
  const r = num('dr');
  const angle = num('da') * Math.PI / 180;
  plate.setDrive(r * Math.cos(angle), r * Math.sin(angle));
}

// ---------- calibration ----------

let lastFit = null;

const measuredFrequencies = () => $('measured').value
  .split(/[\s,;]+/)
  .map(Number)
  .filter(f => f > 0);

async function fitToMeasurements() {
  const measured = measuredFrequencies();
  const config = plateConfig();
  lastFit = null;
  $('applyFit').hidden = true;
  if (measured.length < 2 || !config) {
    $('fitResult').textContent = 'Enter at least two frequencies, and check the plate inputs.';
    return;
  }
  $('fit').disabled = true;
  $('fitResult').textContent = 'Fitting…';
  const fit = await runTask('fitPlate', measured, config);
  $('fit').disabled = false;
  if (!fit) {
    $('fitResult').textContent = 'No fit: the model has no resonances in that range.';
    return;
  }
  lastFit = fit;
  $('applyFit').hidden = false;
  const massNote = config.mount === 'shaker' ? `, moving mass ${(fit.ma * 1000).toFixed(1)} g` : '';
  const verdict = fit.chance > 2 * fit.rms ? 'well above chance' : 'not clearly better than chance, so treat this with caution';
  $('fitResult').textContent = `Frequency scale ×${fit.scale.toFixed(3)}${massNote}. Error ${(100 * fit.rms).toFixed(1)} % RMS `
    + `against ${(100 * fit.chance).toFixed(1)} % for a wrong model: ${verdict}.`;
}

function applyFit() {
  if (!lastFit) return;
  $('E').value = +(lastFit.E / 1e9).toPrecision(4);
  $('ma').value = +(lastFit.ma * 1000).toPrecision(4);
  $('material').value = 'custom';
  showMaterialSource();
  rebuildNow();
}

// ---------- audio ----------

const SOURCE_BUTTONS = { file: 'srcFile', mic: 'srcMic', tone: 'srcTone' };
let actx, analyser, outputGain, buses, oscillator, mic, spectrumDb;
let source = null; // 'file' | 'mic' | 'tone' | null
let toneHz = num('toneHz');

const applyHearOutput = () => {
  if (outputGain) outputGain.gain.value = $('hear').checked ? 1 : 0;
};
const applyWaveform = () => {
  if (oscillator) oscillator.type = $('wave').value;
};

function setUpAudio() {
  if (actx) return;
  actx = new AudioContext();
  analyser = actx.createAnalyser();
  outputGain = actx.createGain();
  applyHearOutput();
  outputGain.connect(actx.destination);
  buses = {};
  for (const name of Object.keys(SOURCE_BUTTONS)) {
    const bus = actx.createGain();
    buses[name] = bus;
    bus.gain.value = 0;
    bus.connect(analyser);
    if (name !== 'mic') bus.connect(outputGain); // never monitor the microphone, to avoid feedback
  }
  actx.createMediaElementSource($('player')).connect(buses.file);
  oscillator = actx.createOscillator();
  oscillator.type = $('wave').value;
  oscillator.frequency.value = toneHz;
  oscillator.connect(buses.tone);
  oscillator.start();
  applyAnalysisSettings();
}

function applyAnalysisSettings() {
  if (!analyser) return;
  analyser.fftSize = num('fft');
  analyser.smoothingTimeConstant = num('smooth');
  spectrumDb = new Float32Array(analyser.frequencyBinCount);
}

const stopTracks = stream => stream.getTracks().forEach(track => track.stop());

// Only the chosen source's bus is open and only its controls are shown.
function useSource(next) {
  if (next || actx) {
    setUpAudio();
    actx.resume();
  }
  source = next;
  document.body.dataset.source = next ?? 'none';
  if (buses) {
    for (const [name, bus] of Object.entries(buses)) {
      const open = name === 'tone' ? num('level') : 1;
      bus.gain.value = name === next ? open : 0;
    }
  }
  if (next !== 'file') $('player').pause();
  if (next !== 'mic' && mic) {
    stopTracks(mic.stream);
    mic.node.disconnect();
    mic = null;
  }
  for (const [name, id] of Object.entries(SOURCE_BUTTONS)) $(id).classList.toggle('on', next === name);
}

const toggleSource = name => useSource(source === name ? null : name);

function setTone(f) { // one frequency for the number box, the log slider and the oscillator
  toneHz = clamp(f, 1, 20000);
  $('toneHz').value = +toneHz.toFixed(2);
  $('tone').value = 1000 * Math.log(Math.max(20, toneHz) / 20) / Math.log(1000);
  if (oscillator) oscillator.frequency.setTargetAtTime(toneHz, actx.currentTime, 0.01);
}

function playTone(f) {
  setTone(f);
  useSource('tone');
}

async function startMic() {
  useSource('mic');
  try {
    const stream = await navigator.mediaDevices.getUserMedia({
      audio: { echoCancellation: false, noiseSuppression: false, autoGainControl: false },
    });
    if (source !== 'mic') { // switched away while asking
      stopTracks(stream);
      return;
    }
    const node = actx.createMediaStreamSource(stream);
    node.connect(buses.mic);
    mic = { stream, node };
  } catch (error) {
    useSource(null);
    showStatus(`Microphone unavailable: ${error.message}`, true);
  }
}

listen('file', 'change', () => {
  const file = $('file').files[0];
  const player = $('player');
  if (!file) return;
  if (player.src.startsWith('blob:')) URL.revokeObjectURL(player.src);
  player.src = URL.createObjectURL(file);
  useSource('file');
  player.play();
});
listen('player', 'play', () => {
  if (source !== 'file') useSource('file');
});
listen('srcFile', 'click', () => {
  if (source === 'file') {
    useSource(null);
    return;
  }
  useSource('file');
  if ($('player').src) $('player').play();
  else $('file').click();
});
listen('srcMic', 'click', () => {
  if (source === 'mic') useSource(null);
  else startMic();
});
listen('micStop', 'click', () => useSource(null));
listen('srcTone', 'click', () => toggleSource('tone'));
listen('tone', 'input', () => setTone(20 * 1000 ** (num('tone') / 1000))); // log slider, 20 Hz – 20 kHz
listen('toneHz', 'input', () => {
  const f = num('toneHz');
  if (f >= 1 && f <= 20000) setTone(f);
});
listen('toneHz', 'change', () => setTone(num('toneHz') || toneHz)); // tidy an out-of-range entry on blur
listen('wave', 'input', applyWaveform);
listen('freqSet', 'input', () => {
  if ($('freqSet').value) playTone(num('freqSet'));
});
listen('level', 'input', () => {
  if (source === 'tone') buses.tone.gain.value = num('level');
});
listen('hear', 'input', applyHearOutput);
listen(['fft', 'smooth'], 'input', applyAnalysisSettings);

// ---------- sand ----------

let grainCount = 0;
let grainX = new Float32Array(0);
let grainY = new Float32Array(0);

function scatterSand() {
  grainCount = num('grains');
  grainX = new Float32Array(grainCount);
  grainY = new Float32Array(grainCount);
  for (let p = 0; p < grainCount; p++) {
    let u, v;
    do {
      u = Math.random() * 2 - 1;
      v = Math.random() * 2 - 1;
    } while (!plate.contains(u, v));
    grainX[p] = u;
    grainY[p] = v;
  }
  sandGL?.setGrains(grainX, grainY);
}

// Bilinear interpolation on a G×G grid at continuous grid coordinates, which must already be clamped to the grid.
const bilinear = (values, G, fx, fy) => {
  const i = fx | 0;
  const j = fy | 0;
  const ax = fx - i;
  const ay = fy - j;
  const o = j * G + i;
  return (values[o] + ax * (values[o + 1] - values[o])) * (1 - ay)
    + (values[o + G] + ax * (values[o + G + 1] - values[o + G])) * ay;
};
const gridEnd = G => G - 1.001; // largest coordinate that still has a cell to its right and below

// Grains move by the drift and diffusivity calculated from the plate's motion in grains.js, in real units. Each step here
// covers no more simulated time than this, so that a fast time-lapse does not jump across the grid.
const MAX_SIM_STEP = 0.05;
const MAX_STEPS_PER_FRAME = { gpu: 16, cpu: 4 };
const sandRandom = Uint32Array.of(12345);
let transport = null;
let gpuFields = null;
let currentDrives = []; // this frame's per-drive deflection in metres

const grainSpecies = () => SPECIES[$('grainType').value];
const metresPerUnit = () => (plate.circle ? activeConfig.d / 2 : Math.max(activeConfig.Lx, activeConfig.Ly) / 2);

function updateTransport() {
  const { G } = plate;
  if (transport?.ux.length !== G * G) {
    transport = newFields(G);
    gpuFields = new Float32Array(G * G * 4);
  }
  const { h, E, rho, nu } = activeConfig;
  currentDrives = plate.deflections(gain, num('g'));
  transportFields({
    G,
    drives: currentDrives,
    plate: { h, unit: metresPerUnit(), E, rho, nu },
    grain: grainSpecies(),
    fluid: $('fluid').value === 'vacuum' ? null : GRAIN_AIR,
    g: num('g'),
    edge: plate.edge,
  }, transport);
}

function uploadFields() {
  const { ux, uy, D, level } = transport;
  const { gamma } = plate;
  for (let q = 0, o = 0; q < gamma.length; q++, o += 4) {
    gpuFields[o] = ux[q];
    gpuFields[o + 1] = uy[q];
    gpuFields[o + 2] = D[q];
    gpuFields[o + 3] = level[q];
  }
  sandGL.setFields(plate.G, gpuFields, gamma);
}

function moveSand(dt) {
  const simTime = dt * num('timeLapse');
  const limit = sandGL ? MAX_STEPS_PER_FRAME.gpu : MAX_STEPS_PER_FRAME.cpu;
  const steps = clamp(Math.ceil(simTime / MAX_SIM_STEP), 1, limit);
  const { circle, hx, hy, G } = plate;
  const capacity = cellCapacity(grainSpecies(), G, metresPerUnit());
  if (sandGL) sandGL.updateCrowding(grainCount, capacity);
  else updateCrowding(grainX, grainY, grainCount, transport, G, capacity);
  for (let i = 0; i < steps; i++) {
    if (sandGL) sandGL.step(simTime / steps, { circle, hx, hy });
    else stepGrains(grainX, grainY, grainCount, transport, G, hx, hy, circle, simTime / steps, sandRandom);
  }
}

// ---------- plate rendering ----------

const plateCanvas = $('plate');
const sandGL = new URLSearchParams(location.search).has('cpu') ? null : createSandGL(plateCanvas);
const plateCtx = sandGL ? null : plateCanvas.getContext('2d');
const marksCanvas = $('plateMarks');
const marksCtx = marksCanvas.getContext('2d');
const BLACK = [0, 0, 0];
const WHITE = [255, 255, 255];
const hexToRgb = hex => [1, 3, 5].map(i => parseInt(hex.slice(i, i + 2), 16));
const packRgba = ([r, g, b]) => (255 << 24 | b << 16 | g << 8 | r) >>> 0; // ImageData is RGBA bytes, little-endian
const mixRgb = (a, b, t) => a.map((value, i) => Math.round(value + (b[i] - value) * t));
const MAGMA = [[0, 0, 4], [59, 15, 112], [140, 41, 129], [222, 73, 104], [254, 159, 109], [252, 253, 191]];
const GAMMA_COLORS = Uint32Array.from({ length: 512 }, (_, i) => { // dimmed half for Γ < 1, where sand can rest
  const t = (i & 255) / 255 * (MAGMA.length - 1);
  const a = Math.floor(t);
  const color = mixRgb(MAGMA[a], MAGMA[Math.min(a + 1, MAGMA.length - 1)], t - a);
  return packRgba(i < 256 ? mixRgb(color, BLACK, 0.55) : color);
});
let size = 0; // canvas side in device pixels
let image, pixels, grainDensity, basePixels, onPlate, sandColors;
let marksKey = null; // what the marks canvas shows now

function resize() {
  const box = $('stage').getBoundingClientRect();
  const cssSize = Math.max(120, Math.floor(Math.min(box.width, box.height) - 20));
  size = Math.min(1400, Math.round(cssSize * devicePixelRatio));
  marksKey = null; // resizing a canvas clears it
  for (const canvas of [plateCanvas, marksCanvas]) {
    canvas.width = size;
    canvas.height = size;
    canvas.style.width = `${cssSize}px`;
    canvas.style.height = `${cssSize}px`;
  }
  if (sandGL) {
    sandGL.resize(size);
  } else {
    image = plateCtx.createImageData(size, size);
    pixels = new Uint32Array(image.data.buffer);
    grainDensity = new Uint8Array(size * size);
  }
  buildBase();
}

function buildBase() {
  if (!size || !plate.G) return;
  const background = hexToRgb($('bg').value);
  const plateColor = hexToRgb($('plateCol').value);
  const sandColor = hexToRgb($('sandCol').value);
  basePixels = new Uint32Array(size * size);
  onPlate = new Uint8Array(size * size);
  for (let y = 0; y < size; y++) {
    for (let x = 0; x < size; x++) {
      const u = (x + 0.5) / size * 2 - 1;
      const v = (y + 0.5) / size * 2 - 1;
      const k = y * size + x;
      if (!plate.contains(u, v)) {
        basePixels[k] = packRgba(background);
        continue;
      }
      const toEdge = plate.circle ? 1 - Math.hypot(u, v) : Math.min(plate.hx - Math.abs(u), plate.hy - Math.abs(v));
      const shaded = toEdge < 0.008
        ? mixRgb(plateColor, WHITE, 0.25) // rim
        : mixRgb(plateColor, BLACK, 0.15 * (u * u + v * v)); // dish shading
      basePixels[k] = packRgba(shaded);
      onPlate[k] = 1;
    }
  }
  sandColors = Uint32Array.from({ length: 32 }, (_, density) => {
    const covered = mixRgb(plateColor, sandColor, 1 - 0.55 ** density);
    return packRgba(mixRgb(covered, WHITE, Math.max(0, density - 6) / 80));
  });
  if (sandGL) {
    // the shader reads alpha as "this pixel is on the plate"
    sandGL.setBase(basePixels.map((pixel, k) => (onPlate[k] ? pixel : pixel & 0x00ffffff)));
    sandGL.setSandColors(sandColors);
    sandGL.setMapColors(GAMMA_COLORS);
  }
}

const grainPixels = () => Math.max(1, Math.round(num('grain') * size / 600)); // scales with the plate, not the screen

function countGrains(half) {
  grainDensity.fill(0);
  const grainSize = grainPixels();
  const limit = size - grainSize;
  for (let p = 0; p < grainCount; p++) {
    const x = Math.min(limit, ((grainX[p] + 1) * half) | 0);
    const y = Math.min(limit, ((grainY[p] + 1) * half) | 0);
    for (let dy = 0; dy < grainSize; dy++) {
      for (let dx = 0; dx < grainSize; dx++) {
        const k = (y + dy) * size + x + dx;
        if (grainDensity[k] < 255) grainDensity[k]++;
      }
    }
  }
}

function paintPixels() {
  const showMap = $('overlay').checked;
  const { gamma, G } = plate;
  const mapTop = Math.max(2, peakGamma);
  const cellsPerPixel = G / size;
  const end = gridEnd(G);
  for (let y = 0, k = 0; y < size; y++) {
    const fy = clamp((y + 0.5) * cellsPerPixel - 0.5, 0, end);
    for (let x = 0; x < size; x++, k++) {
      const density = grainDensity[k];
      if (density) {
        pixels[k] = sandColors[Math.min(31, density)];
      } else if (showMap && onPlate[k]) {
        const fx = clamp((x + 0.5) * cellsPerPixel - 0.5, 0, end);
        const g = bilinear(gamma, G, fx, fy);
        pixels[k] = GAMMA_COLORS[(g < 1 ? 0 : 256) + Math.min(255, (g / mapTop * 255) | 0)];
      } else {
        pixels[k] = basePixels[k];
      }
    }
  }
  plateCtx.putImageData(image, 0, 0);
}

// Supports: clamped edges solid, simply supported dashed, free edges bare.
function drawSupports(half) {
  const dpr = devicePixelRatio;
  const x0 = (1 - plate.hx) * half;
  const x1 = (1 + plate.hx) * half;
  const y0 = (1 - plate.hy) * half;
  const y1 = (1 + plate.hy) * half;
  const edgeLines = [[x0, y0, x0, y1], [x1, y0, x1, y1], [x0, y0, x1, y0], [x0, y1, x1, y1]];
  edgeLines.forEach((line, i) => {
    const edge = plate.edges[i];
    if (edge === 'free') return;
    marksCtx.strokeStyle = '#9aa3ad';
    marksCtx.lineWidth = (edge === 'clamped' ? 5 : 2) * dpr;
    marksCtx.setLineDash(edge === 'simply' ? [6 * dpr, 5 * dpr] : []);
    strokeLine(marksCtx, ...line);
  });
  marksCtx.setLineDash([]);
}

function drawDriver(half) {
  const dpr = devicePixelRatio;
  marksCtx.strokeStyle = '#ff6b5a';
  marksCtx.lineWidth = 2 * dpr;
  marksCtx.beginPath();
  marksCtx.arc((plate.du + 1) * half, (plate.dv + 1) * half, 6 * dpr, 0, 2 * Math.PI);
  marksCtx.stroke();
}

function updateMarks(half) {
  const key = [size, plate.circle, plate.hx, plate.hy, plate.edges, plate.du, plate.dv, $('showDrive').checked].join();
  if (key === marksKey) return;
  marksKey = key;
  marksCtx.clearRect(0, 0, size, size);
  if (!plate.circle) drawSupports(half);
  if ($('showDrive').checked) drawDriver(half);
  plateView?.setMarks(marksCanvas);
}

function renderPlate(dt) {
  const half = size / 2;
  updateMarks(half);
  if (sandGL) {
    const showMap = $('overlay').checked;
    sandGL.draw(grainCount, { grainSize: grainPixels(), showMap, mapTop: Math.max(2, peakGamma), offscreen: viewIn3D });
    if (viewIn3D) drawPlateIn3D(dt);
  } else {
    countGrains(half);
    paintPixels();
  }
}

// ---------- 3D view ----------

const MAX_LISTED = 4;
const NO_LINES = new Float32Array(0);
const camera = createOrbitCamera();
const plateView = sandGL && createPlateView3D(sandGL.gl);
const stillMotion = matchMedia('(prefers-reduced-motion: reduce)');
const wavePhases = []; // one per tone, lowest first, advancing at true speed ratios from an arbitrary start: the spectrum gives no phase between tones
let plateHeights = new Float32Array(0);
let soloTone = null; // the tone picked in the legend for the lines to ride alone
// shortcut: a tone is known by its rounded frequency, so one that wanders in live audio comes back; fine for a viewing aid
const toneKey = tone => Math.round(tone.frequency);
let viewIn3D = false;

function setView(in3D) {
  viewIn3D = in3D;
  document.body.dataset.view = in3D ? '3d' : '2d';
  $('view2D').classList.toggle('on', !in3D);
  $('view3D').classList.toggle('on', in3D);
  plateCanvas.tabIndex = in3D ? 0 : -1;
  plateCanvas.setAttribute('aria-label', in3D
    ? 'Simulated vibrating plate covered in sand, in 3D. Drag or use the arrow keys to rotate, scroll to zoom, click to move the driver.'
    : 'Simulated vibrating plate covered in sand. Click to move the driver.');
}

const fmtLength = metres => {
  if (metres >= 1) return `${metres.toPrecision(3)} m`;
  if (metres >= 1e-3) return `${(metres * 1e3).toPrecision(3)} mm`;
  return metres >= 1e-6 ? `${(metres * 1e6).toPrecision(3)} µm` : `${(metres * 1e9).toPrecision(3)} nm`;
};

function updateLegend(magnified, slowdown, driftTop) {
  const showLines = $('fieldLines').checked;
  const tones = strongestTones(currentDrives, num('g'), MAX_LISTED).filter(tone => tone.peak > 0);
  const solo = showLines ? tones.find(tone => toneKey(tone) === soloTone) : undefined;
  const rows = tones.map(tone => {
    const on = !solo || tone === solo;
    const pick = showLines ? `<button type="button" aria-pressed="${tone === solo}" aria-label="Lines ride only the ${fmtHz(tone.frequency)} tone">${on ? '■' : '□'}</button>` : '';
    return `<tr${showLines ? ` data-tone="${toneKey(tone)}"` : ''}${on ? '' : ' class="off"'}><td>${pick}</td>`
      + `<td>${fmtHz(tone.frequency)}</td><td>${fmtLength(tone.peak)}</td><td>${tone.shaking.toPrecision(3)} g</td></tr>`;
  }).reverse();
  const header = `<tr${showLines ? ' title="Click a tone for the lines to ride it alone, and again for all tones"' : ''}><th></th><th>Tone</th><th title="The largest real up-and-down motion of this tone alone">Height</th>`
    + '<th title="Peak acceleration of this tone alone. Above about 1 g, grains hop">Shaking</th></tr>';
  const notes = [];
  if (tones.length) {
    const shown = currentDrives.length > tones.length ? [`${tones.length} of ${currentDrives.length} tones`] : [];
    if (showLines) shown.push(`lines ${magnified ? `×${magnified.toLocaleString('en')}` : 'flat'}`, `slowed ×${Math.round(slowdown).toLocaleString('en')}`);
    if (shown.length) notes.push(shown.join(' · '));
    if (showLines && magnified && !solo && currentDrives.length > 1) notes.push('together: one possible timing');
  }
  if (driftTop > 0) notes.push(`field lines: grain drift, up to ${fmtLength(driftTop)}/s`);
  else if (tones.length && showLines) notes.push('no field lines: shaking too weak to move grains');
  const scale = num('grain');
  notes.push(`grains ${fmtLength(grainSpecies().D)}${scale === 1 ? ' true size' : `, drawn ×${scale}`}, at least 1 pixel`);
  const legend = $('waveLegend');
  const html = (tones.length ? header + rows.join('') : '') + notes.map(note => `<tr><td colspan="4">${note}</td></tr>`).join('');
  if (legend.innerHTML !== html) legend.innerHTML = html;
  legend.hidden = !html;
}

function drawPlateIn3D(dt) {
  const { G, circle, hx, hy } = plate;
  const unit = metresPerUnit();
  const species = grainSpecies();
  const tones = strongestTones(currentDrives, num('g'), currentDrives.length).filter(tone => tone.peak > 0);
  const solo = tones.find(tone => toneKey(tone) === soloTone);
  const ridden = solo ? [solo] : tones;
  const slowdown = slowMotion(ridden.map(tone => tone.frequency));
  tones.forEach((tone, index) => {
    wavePhases[index] = stillMotion.matches ? tone.peakPhase : (wavePhases[index] ?? 0) + tone.drive.omega * dt / slowdown;
    tone.phase = wavePhases[index];
  });
  if (plateHeights.length !== G * G) plateHeights = new Float32Array(G * G);
  const chosen = $('heightScale').value;
  const magnified = chosen === 'off' ? 0 : chosen === 'auto' ? autoExaggeration(solo ? solo.peak : peakDeflection(currentDrives), unit) : +chosen;
  const traced = $('fieldLines').checked && transport ? traceFieldLines(transport, plate) : { vertices: NO_LINES, topSpeed: 0 };

  plateView.draw({
    camera,
    grid: G,
    image: sandGL.image,
    grains: {
      positions: sandGL.positions,
      side: sandGL.side,
      count: grainCount,
      counts: sandGL.counts,
      heapScale: heapHeightPerGrain(species, G, unit) / unit,
      diameter: num('grain') * species.D / unit,
      colour: hexToRgb($('sandCol').value).map(value => value / 255),
    },
    lines: { vertices: traced.vertices, heights: traced.vertices.length ? plateHeight(ridden, plateHeights, plate.edge) : plateHeights, scale: magnified / unit },
    circle,
    hx,
    hy,
    background: hexToRgb($('bg').value).map(value => value / 255),
    size,
  });
  if (frameCount % 6 === 0) updateLegend(magnified, slowdown, traced.topSpeed * unit);
}

// ---------- spectrum rendering ----------

const spectrumCanvas = $('spectrum');
const spectrumCtx = spectrumCanvas.getContext('2d');
const DB_MIN = -120;
const DB_MAX = -10;
const AXIS_LABELS_HZ = [20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000];

function drawSpectrum() {
  const ctx = spectrumCtx;
  const dpr = devicePixelRatio;
  const W = Math.round(spectrumCanvas.clientWidth * dpr);
  const H = Math.round(spectrumCanvas.clientHeight * dpr);
  if (spectrumCanvas.width !== W || spectrumCanvas.height !== H) {
    spectrumCanvas.width = W;
    spectrumCanvas.height = H;
  }
  ctx.clearRect(0, 0, W, H);

  const f0 = Math.max(1, num('fmin'));
  const f1 = Math.max(f0 * 1.01, num('fmax'));
  const logSpan = Math.log(f1 / f0);
  const toX = f => Math.log(f / f0) / logSpan * W;
  const toY = db => H - 18 * dpr - (Math.max(DB_MIN, db) - DB_MIN) / (DB_MAX - DB_MIN) * (H - 30 * dpr);
  ctx.lineWidth = dpr;
  ctx.font = `${11 * dpr}px system-ui, sans-serif`;

  // the plate's eigenfrequencies, coloured by nodal diameters
  for (const mode of plate.modes) {
    if (mode.f < f0 || mode.f > f1) continue;
    ctx.strokeStyle = `hsl(${(mode.n ?? 1) * 47 % 360} 70% 60% / 0.35)`;
    strokeLine(ctx, toX(mode.f), 0, toX(mode.f), H - 18 * dpr);
  }

  ctx.fillStyle = '#8b919b';
  for (const f of AXIS_LABELS_HZ) {
    if (f >= f0 && f <= f1) ctx.fillText(f >= 1000 ? `${f / 1000}k` : `${f}`, toX(f) + 2 * dpr, H - 5 * dpr);
  }

  // the detection floor
  ctx.strokeStyle = '#8b919b80';
  ctx.setLineDash([4 * dpr, 4 * dpr]);
  strokeLine(ctx, 0, toY(num('floor')), W, toY(num('floor')));
  ctx.setLineDash([]);

  if (!source || !spectrumDb) return;

  // the live spectrum: each pixel column shows the loudest bin it covers
  const binHz = actx.sampleRate / analyser.fftSize;
  ctx.strokeStyle = $('sandCol').value;
  ctx.beginPath();
  for (let x = 0; x < W; x++) {
    const firstBin = Math.floor(f0 * Math.exp(logSpan * x / W) / binHz);
    const lastBin = Math.max(firstBin, Math.floor(f0 * Math.exp(logSpan * (x + 1) / W) / binHz));
    let loudest = -Infinity;
    for (let i = firstBin; i <= lastBin && i < spectrumDb.length; i++) {
      if (spectrumDb[i] > loudest) loudest = spectrumDb[i];
    }
    if (x) ctx.lineTo(x, toY(loudest));
    else ctx.moveTo(x, toY(loudest));
  }
  ctx.stroke();

  peaks.forEach((peak, i) => {
    ctx.fillStyle = '#ff6b5a';
    ctx.beginPath();
    ctx.arc(toX(peak.f), toY(peak.db), 3 * dpr, 0, 2 * Math.PI);
    ctx.fill();
    if (i < 5) {
      ctx.fillStyle = '#d9dce1';
      ctx.fillText(fmtHz(peak.f), toX(peak.f) + 5 * dpr, toY(peak.db) - 4 * dpr);
    }
  });
}

// ---------- readout ----------

const logDistance = (a, b) => Math.abs(Math.log(a / b));
const nearestMode = f => plate.modes.reduce(
  (best, mode) => (mode.f > 0 && (!best || logDistance(mode.f, f) < logDistance(best.f, f)) ? mode : best),
  null,
);

// Stiffening depends on absolute amplitude, so its coefficients are only worked out once a drive is strong enough.
let kappasRequestedFor = null;

function prepareStiffeningIfNeeded() {
  if (!plate.stiffeningOn || plate.hasKappas || !peaks.length) return;
  const omega = 2 * Math.PI * peaks[0].f;
  if (stiffeningRegime(peakGamma * num('g'), omega, plate.h, num('zeta')) === 'negligible') return;
  const config = activeConfig;
  if (plate.circle) {
    plate.prepareStiffening(config);
    return;
  }
  if (plate.ritz.modes !== plate.solvedModes || kappasRequestedFor === plate.solvedModes) return;
  kappasRequestedFor = plate.solvedModes;
  const asked = plate.solvedModes;
  // default basis only, so a worker thread's modes line up with these
  if (solvePlan(config)) return;
  runTask('rectKappasDefault', config.edges, config.nu, config.Ly / config.Lx).then(kappas => {
    if (plate.solvedModes === asked) plate.setKappas(kappas);
  });
}

const STIFFENING_NOTES = {
  corrected: '  stiffening included',
  caution: '  ⚠ stiffening beyond its tested range',
  broken: '  ⚠ one-mode stiffening has broken down',
};

function stiffeningNote(a, omega) {
  const state = stiffeningRegime(a, omega, plate.h, num('zeta'));
  if (state === 'negligible') return '';
  if (state === 'corrected' && !plate.detune) return '  ⚠ stiffening not modelled on this mount';
  return STIFFENING_NOTES[state];
}

// Heaping and Faraday waves in a deep vibrated layer (Melo et al. 1995): Γ above about 2.2, more than two grains deep,
// up to about 110 Hz. This simulation has no grain contact, so it says so.
function layerNote(frequency) {
  const species = grainSpecies();
  const area = plate.circle ? Math.PI * (num('d') / 200) ** 2 : (num('lx') / 100) * (num('ly') / 100);
  const layersDeep = grainCount * Math.PI * (species.D / 2) ** 2 / area;
  return peakGamma > 2.2 && frequency < 110 && layersDeep > 2 ? '⚠ deep layer at low frequency: granular waves are not modelled' : '';
}

function updateReadout() {
  const lines = [];
  const strongest = peaks[0];
  if (!source) {
    lines.push('Pick a sound: audio file, microphone or test tone');
  } else if (!strongest) {
    lines.push('No peaks above the floor — the sand is at rest');
  } else {
    lines.push(`strongest     ${fmtHz(strongest.f)}  ${fmtNote(strongest.f)}  ${strongest.db.toFixed(1)} dB`);
    const mode = nearestMode(strongest.f);
    if (mode) lines.push(`nearest mode  ${mode.label}  ${fmtHz(mode.f)}  ${fmtCents(1200 * Math.log2(strongest.f / mode.f))}`);
  }
  lines.push(`peak plate acceleration ${peakGamma.toFixed(2)} g`);
  if (strongest) { // linear plate theory needs the deflection well below the thickness
    const amplitude = peakDeflection(currentDrives);
    const ratio = amplitude / (num('h') / 1000);
    const omega = 2 * Math.PI * strongest.f;
    lines.push(`peak deflection ${currentDrives.length > 1 ? 'up to' : '≈'} ${(amplitude * 1e6).toPrecision(3)} µm = ${(100 * ratio).toFixed(1)} % of thickness${stiffeningNote(peakGamma * num('g'), omega)}`);
    const note = layerNote(strongest.f);
    if (note) lines.push(note);
  }
  lines.push(
    `drive ${gain.toPrecision(3)} ${plate.mount === 'rigid' ? 'g' : 'N'} at 0 dBFS`,
    `field ${fieldMs.toFixed(1)} ms · ${grainCount / 1000}k grains ≈ ${(grainCount * grainMass(grainSpecies()) * 1000).toPrecision(2)} g · ${sandGL ? 'GPU' : 'CPU'} · time ×${num('timeLapse')}`,
  );
  $('readout').textContent = lines.join('\n');
}

// ---------- main loop ----------

// A full-scale sine reads −13.5 dB: the analyser's Blackman window has coherent gain 0.42, halved for one side.
const FULL_SCALE_AMPLITUDE = 0.21;
let peaks = [];
let gain = num('force');
let peakGamma = 0;
let fieldMs = 0;
let paused = false;
let lastFrameAt = performance.now();
let frameCount = 0;

function readPeaks() {
  if (!source || !spectrumDb) return [];
  analyser.getFloatFrequencyData(spectrumDb);
  return findPeaks(spectrumDb, actx.sampleRate / analyser.fftSize, num('fmin'), num('fmax'), num('peaks'), num('floor'));
}

function updateGain(unitPeak, dt) {
  if (!$('auto').checked) {
    gain = num('force');
  } else if (unitPeak > 0) { // fast attack, slow release
    const wanted = Math.min(1e4, num('target') / unitPeak);
    gain = wanted < gain ? wanted : gain + (wanted - gain) * Math.min(1, 2 * dt);
  }
}

function frame(now) {
  const dt = Math.min(0.05, (now - lastFrameAt) / 1000);
  lastFrameAt = now;
  peaks = readPeaks();

  const fieldStartedAt = performance.now();
  const drives = peaks.map(peak => ({ f: peak.f, F: 10 ** (peak.db / 20) / FULL_SCALE_AMPLITUDE }));
  const unitPeak = plate.field(drives, num('g'));
  fieldMs += (performance.now() - fieldStartedAt - fieldMs) * 0.05;

  updateGain(unitPeak, dt);
  plate.forceGain = gain;
  for (let i = 0; i < plate.gamma.length; i++) plate.gamma[i] *= gain; // Γ is linear in the force
  peakGamma = unitPeak * gain;

  updateTransport();
  if (sandGL) uploadFields();
  if (!paused) moveSand(dt);
  prepareStiffeningIfNeeded();
  renderPlate(dt);
  drawSpectrum();
  if (++frameCount % 6 === 0) updateReadout();
  requestAnimationFrame(frame);
}

// ---------- wiring ----------

const edgeSelects = document.querySelectorAll('.edge');
const setAllEdges = value => {
  for (const select of edgeSelects) select.value = value;
};

listen(
  ['shape', 'bc', 'E', 'rho', 'nu', 'd', 'lx', 'ly', 'h', 'zeta', 'res', 'fmax', 'mount', 'ma', 'mindlin', 'eL', 'eR', 'eT', 'eB', 'fluid', 'clamp', 'stiffening', 'dampingModel'],
  'input',
  scheduleRebuild,
);
for (const select of edgeSelects) select.innerHTML = $('bc').innerHTML;
setAllEdges($('bc').value);
listen('bc', 'input', () => setAllEdges($('bc').value)); // "Edge" sets all four
listen(['E', 'rho', 'nu'], 'input', () => {
  $('material').value = 'custom';
  showMaterialSource();
});
listen('material', 'input', () => {
  const material = MATERIALS[$('material').value];
  if (material) {
    $('E').value = material.E;
    $('rho').value = material.rho;
    $('nu').value = material.nu;
  }
  showMaterialSource();
  scheduleRebuild();
});
listen(['dr', 'da'], 'input', applyDriver);
listen(['sandCol', 'plateCol', 'bg'], 'input', buildBase);
listen('grains', 'change', scatterSand);
listen('resprinkle', 'click', scatterSand);
$('grainType').replaceChildren(...Object.entries(SPECIES).map(([key, { name }]) => new Option(name, key)));
form.reset(); // the grain list only now exists, so its default can be restored
listen('fit', 'click', fitToMeasurements);
listen('applyFit', 'click', applyFit);
listen('pause', 'click', () => {
  paused = !paused;
  $('pause').textContent = paused ? 'Resume' : 'Pause';
  $('pause').classList.toggle('on', paused);
});
function moveDriverTo(u, v) {
  $('dr').value = Math.min(DRIVER_LIMIT, Math.hypot(u, v));
  $('da').value = Math.atan2(v, u) * 180 / Math.PI;
  for (const id of ['dr', 'da']) $(id).dispatchEvent(new Event('input', { bubbles: true }));
}
plateCanvas.addEventListener('click', e => {
  if (viewIn3D) return; // orbitControls handles clicks there, to tell them from drags
  const box = plateCanvas.getBoundingClientRect();
  moveDriverTo((e.clientX - box.left) / box.width * 2 - 1, (e.clientY - box.top) / box.height * 2 - 1);
});
if (plateView) {
  $('viewBar').hidden = false;
  orbitControls(plateCanvas, camera, {
    isActive: () => viewIn3D,
    onClick: (x, y) => {
      const hit = pickPlane(camera, x, y);
      if (hit && plate.contains(...hit)) moveDriverTo(...hit);
    },
  });
  listen('view2D', 'click', () => setView(false));
  listen('view3D', 'click', () => setView(true));
  listen('waveLegend', 'click', event => {
    const key = event.target.closest('[data-tone]')?.dataset.tone;
    if (key) soloTone = +key === soloTone ? null : +key;
  });
}
new ResizeObserver(resize).observe($('stage'));

// ---------- reset and presets ----------

// Each preset starts from the defaults, applies its settings in order ('bc' before single edges), clears the sand,
// and plays a tone at a named eigenmode or a fixed frequency.
const PRESETS = [
  {
    name: 'Eight-point star',
    note: 'free circle, (8,1) mode',
    set: { dr: 0.6, da: -90, target: 12 },
    mode: '(8,1)',
  },
  {
    name: 'Nested rings',
    note: 'circle bolted to a shaker',
    set: { mount: 'shaker', ma: 80, E: 69, nu: 0.35, h: 0.8, target: 15 },
    f: 1700,
  },
  {
    name: 'Chladni square',
    note: 'free square, centre drive',
    set: { shape: 'rect', target: 12 },
    mode: '#38',
  },
  {
    name: 'Shu et al. 1247 Hz',
    note: 'bolted square, as photographed',
    // moving mass fitted to their 25 resonances with air included (1.9 % RMS); without air the same fit gives 20 g
    set: { shape: 'rect', mount: 'shaker', ma: 15, E: 70, lx: 32, ly: 32, fmax: 5000, target: 15 },
    f: 1247,
  },
  // fitted to a widely shared 12-pattern photo (345–5907 Hz): only h/L² and the mass ratio are determined, ±1 %
  {
    name: '12-pattern photo plate',
    note: '24 cm × 0.75 mm Al, 7 g shaker',
    set: { shape: 'rect', lx: 24, ly: 24, h: 0.75, mount: 'shaker', ma: 7, target: 12 },
    f: 2041,
  },
  {
    name: 'Cantilever cross',
    note: 'top edge clamped',
    set: { shape: 'rect', eT: 'clamped', lx: 19, ly: 19, h: 3.15, dr: 0.75, da: 40, fmax: 3000, target: 12 },
    mode: '#5',
  },
  {
    name: 'Clamped drum',
    note: 'clamped circle, (2,1) mode',
    set: { bc: 'clamped', dr: 0.3, da: -90, target: 12 },
    mode: '(2,1)',
  },
  {
    name: 'Square-wave harmonics',
    note: 'free circle, 220 Hz',
    set: { dr: 0.45, da: -90, target: 10 },
    f: 220,
    wave: 'square',
  },
  // Faraday 1831: under strong drive, fine powder is carried by air to the antinodes, the opposite of sand
  {
    name: 'Faraday powder',
    note: 'lycopodium gathers at antinodes',
    set: { grainType: 'lycopodium', h: 2, target: 100 },
    mode: '(0,1)',
  },
  {
    name: 'Grid of cells',
    note: 'simply supported 2:1 rectangle',
    set: { shape: 'rect', bc: 'simply', lx: 40, ly: 20, dr: 0.3, da: 30, target: 12 },
    mode: '#9',
  },
];

function syncForm() { // derived state after the whole form changed (page load, reset, preset)
  updateRangeOutputs();
  syncLineControls();
  showMaterialSource();
  rebuildNow();
  setTone(num('toneHz'));
  applyWaveform();
  applyHearOutput();
  applyAnalysisSettings();
  scatterSand();
  gain = num('force');
}

function applyPreset(preset) {
  form.reset();
  for (const [id, value] of Object.entries(preset.set)) {
    if ($(id).type === 'checkbox') $(id).checked = value;
    else $(id).value = value;
    if (id === 'bc') setAllEdges(value);
  }
  if (preset.wave) $('wave').value = preset.wave;
  syncForm();
  const f = preset.f ?? plate.modes.find(mode => mode.label === preset.mode)?.f;
  if (f) playTone(f);
  else showStatus(`Preset mode ${preset.mode} is not in this plate's mode list`, true);
}

function presetButton(preset) {
  const note = document.createElement('small');
  note.textContent = preset.note;
  return makeButton(() => applyPreset(preset), preset.name, note);
}

$('presets').replaceChildren(...PRESETS.map(presetButton));
listen('fieldLines', 'input', syncLineControls);
listen('resetAll', 'click', () => {
  form.reset();
  // the view bar sits outside the form, so presets leave it alone and only this resets it
  $('fieldLines').checked = $('fieldLines').defaultChecked;
  $('heightScale').value = 'auto';
  soloTone = null;
  syncForm();
  useSource(null);
});

if (sandGL) $('grains').max = 3000000;
syncForm();
resize();
requestAnimationFrame(frame);
