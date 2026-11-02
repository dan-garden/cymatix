// The plate's true shape at an instant, for drawing it in 3D. Pure physics in metres, with no drawing here.

/**
 * w(x, t) = Σ_j Re[W_j·e^{iω_j·t}] on the grid, from each drive's complex deflection in metres (Plate.deflections).
 * Cells just outside the edge copy the nearest plate cell (edge, from Plate), so the surface does not fall to zero at the rim.
 * shortcut: the spectrum gives no phase between tones, so tones all start in phase; a single tone is exact.
 */
export function deflectionAt(drives, time, out, edge = []) {
  out.fill(0);
  for (const { omega, re, im } of drives) {
    const c = Math.cos(omega * time), s = Math.sin(omega * time);
    for (let q = 0; q < out.length; q++) out[q] += re[q] * c - im[q] * s;
  }
  for (let k = 0; k < edge.length; k += 2) out[edge[k]] = out[edge[k + 1]];
  return out;
}

/** The largest |w| the plate can reach at any cell and time: the tones' amplitudes, all lined up. */
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

/** How far to magnify heights so the peak fills about `share` of the half-width, as 1, 2 or 5 × 10ⁿ so the label stays readable. */
export function autoExaggeration(peak, unit, share = 0.15) {
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
