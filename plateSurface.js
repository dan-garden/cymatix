// The plate's true motion, tone by tone, for drawing in 3D. Pure physics in metres, with no drawing here.

const TAU = 2 * Math.PI;

/**
 * The plate's height in metres at each cell, as the sum of Re[W·e^{iφ}] over tones given as { drive: { re, im }, phase }.
 * W comes from Plate.deflections. Cells just outside the edge copy the nearest plate cell (edge, from Plate), so the height
 * does not fall to zero at the rim.
 */
export function plateHeight(tones, out, edge = []) {
  out.fill(0);
  for (const { drive: { re, im }, phase } of tones) {
    const c = Math.cos(phase), s = Math.sin(phase);
    for (let q = 0; q < out.length; q++) out[q] += re[q] * c - im[q] * s;
  }
  for (let k = 0; k < edge.length; k += 2) out[edge[k]] = out[edge[k + 1]];
  return out;
}

/**
 * The tones that shake the plate hardest, lowest frequency first so their order holds as strengths change. For each: peak
 * height |W| in metres, shaking in g (peak acceleration ω²|W|/g, which decides whether grains hop), and the phase at which
 * its tallest point reaches full height.
 */
export function strongestTones(drives, g, count) {
  return drives
    .map(drive => {
      const { omega, re, im } = drive;
      let largest = 0, cell = 0;
      for (let q = 0; q < re.length; q++) {
        const size = re[q] * re[q] + im[q] * im[q];
        if (size > largest) {
          largest = size;
          cell = q;
        }
      }
      const peak = Math.sqrt(largest);
      return { drive, frequency: omega / TAU, peak, shaking: omega * omega * peak / g, peakPhase: -Math.atan2(im[cell], re[cell]) };
    })
    .sort((a, b) => b.shaking - a.shaking)
    .slice(0, count)
    .sort((a, b) => a.frequency - b.frequency);
}

/** The largest |w| the plate can reach at any cell and time: the tones' amplitudes, all lined up. Exact for a single tone. */
export function peakDeflection(drives) {
  if (!drives.length) return 0;
  let peak = 0;
  for (let q = 0; q < drives[0].re.length; q++) {
    let sum = 0;
    for (const { re, im } of drives) sum += Math.hypot(re[q], im[q]);
    if (sum > peak) peak = sum;
  }
  return peak;
}

/**
 * How many times to slow time: the lowest tone shows one cycle per second, unless that would push the highest past `fastest`
 * cycles per second, where a screen could no longer show it truly. All tones share it, so their true speed ratios stay.
 */
export function slowMotion(frequencies, fastest = 8) {
  return Math.max(Math.min(...frequencies), Math.max(...frequencies) / fastest);
}

/** How far to magnify heights so the peak fills about `share` of the half-width, as 1, 2 or 5 × 10ⁿ so the label stays readable. */
export function autoExaggeration(peak, unit, share = 0.3) {
  if (!(peak > 0)) return 1;
  const wanted = share * unit / peak;
  if (wanted <= 1) return 1; // already that large: never shrink a real shape
  const power = 10 ** Math.floor(Math.log10(wanted));
  return [5, 2, 1].map(step => step * power).find(value => value <= wanted);
}

/**
 * Height in metres that one grain adds to its grid cell, if the cell's grains lie as a loose heap of solid fraction
 * `packing`: volume conservation, n·(π/6)D³ = height·packing·cell². Grain positions have no height of their own.
 */
export function heapHeightPerGrain({ D, packing }, G, unit) {
  const cellSide = 2 * unit / G;
  return (Math.PI / 6) * D ** 3 / (packing * cellSide * cellSide);
}
