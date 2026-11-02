// Fit a modelled plate to measured resonances. Frequencies scale with h·√(E/ρ)/L², and with the moving mass held as a
// fraction of the plate's mass the whole response only stretches along the frequency axis. So one sweep per mass ratio
// is enough, and the frequency scale can then be searched without solving anything again.
import { PlateModel } from './plateModel.js';

const SWEEP_STEP = 1.001;
const SCALE_STEP = 1.002;
const MASS_RATIOS = [0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.5];

/** Sweeps a configured plate; a resonance is a local maximum standing 20 % above the response 1 % either side. */
export function sweepResonances(plate, fStart, fEnd) {
  const response = [];
  for (let f = fStart; f < fEnd; f *= SWEEP_STEP) response.push([f, plate.field([{ f, F: 1 }], 9.81)]);
  const resonances = [];
  for (let i = 10; i + 10 < response.length; i++) {
    const y = response[i][1];
    const isPeak = y > response[i - 1][1] && y >= response[i + 1][1];
    const standsOut = y > 1.2 * response[i - 10][1] && y > 1.2 * response[i + 10][1];
    if (isPeak && standsOut) resonances.push(response[i][0]);
  }
  return resonances;
}

// RMS relative distance from each measured frequency to the nearest scaled resonance. Both lists are sorted.
function rmsError(measured, resonances, scale) {
  let sum = 0, j = 0;
  for (const f of measured) {
    while (j + 1 < resonances.length && Math.abs(resonances[j + 1] * scale - f) <= Math.abs(resonances[j] * scale - f)) j++;
    sum += (resonances[j] * scale / f - 1) ** 2;
  }
  return Math.sqrt(sum / measured.length);
}

const scalesBetween = (lo, hi) => {
  const scales = [];
  for (let s = lo; s <= hi; s *= SCALE_STEP) scales.push(s);
  return scales;
};

// The error a wrong model reaches by luck: the median over scale factors 0.7–1.3 about the best one.
function chanceError(measured, resonances, bestScale) {
  const errors = scalesBetween(0.7 * bestScale, 1.3 * bestScale).map(s => rmsError(measured, resonances, s));
  return errors.sort((a, b) => a - b)[errors.length >> 1];
}

function bestScale(measured, resonances, minScale) {
  let best = { scale: 1, rms: Infinity };
  for (const scale of scalesBetween(Math.max(0.4, minScale), 2.5)) {
    const rms = rmsError(measured, resonances, scale);
    if (rms < best.rms) best = { scale, rms };
  }
  return best;
}

/**
 * Fits the frequency scale and, for a shaker mount, the moving mass to measured resonances (Hz).
 * config is a Plate configuration in SI units. Returns the best fit: { scale, ma, massRatio, rms, chance, E, resonances },
 * where E is the Young's modulus giving that scale with h, ρ and the sizes kept, and chance is the error a wrong model
 * reaches by luck. A fit means little unless rms is well below chance.
 */
export function fitPlate(measured, config, { massRatios = MASS_RATIOS } = {}) {
  const sorted = [...measured].sort((a, b) => a - b);
  const area = config.shape === 'circle' ? Math.PI * (config.d / 2) ** 2 : config.Lx * config.Ly;
  const plateMass = config.rho * config.h * area;
  const ratios = config.mount === 'shaker' ? massRatios : [0];
  let best = { rms: Infinity };
  for (const massRatio of ratios) {
    const plate = new PlateModel();
    const fTopWanted = 2.5 * sorted.at(-1);
    plate.configure({ ...config, G: 25, fmax: fTopWanted, zeta: Math.min(config.zeta, 0.005), ma: massRatio * plateMass });
    plate.setDrive(0, 0);
    const fEnd = Math.min(fTopWanted, plate.fTop);
    const resonances = sweepResonances(plate, sorted[0] / 2.5, fEnd);
    if (resonances.length < 2) continue;
    // scales below this would need modes the solver does not have
    const fit = bestScale(sorted, resonances, 1.05 * sorted.at(-1) / fEnd);
    if (fit.rms < best.rms) best = { ...fit, massRatio, resonances };
  }
  if (!best.resonances) return null;
  const { scale, massRatio, rms, resonances } = best;
  return {
    scale,
    massRatio,
    rms,
    ma: massRatio * plateMass,
    E: config.E * scale * scale,
    chance: chanceError(sorted, resonances, scale),
    resonances: resonances.map(f => f * scale),
  };
}
