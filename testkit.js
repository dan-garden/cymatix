// Helpers shared by the test files.
import assert from 'node:assert/strict';
import { rectModes } from './physics.js';

export const pct = (got, want) => 100 * Math.abs(got / want - 1);

export const median = values => values.slice().sort((x, y) => x - y)[values.length >> 1];

export const nearest = (values, target) =>
  values.reduce((best, x) => (Math.abs(x - target) < Math.abs(best - target) ? x : best));

export const within = (pairs, maxPct, name) => {
  for (const [got, want, tag] of pairs) {
    assert.ok(pct(got, want) <= maxPct, `${name} ${tag ?? ''}: ours ${got} vs ${want}`);
  }
};

/** Calls visit(value, row, column) for every non-null cell of a table. */
export const forEachCell = (table, visit) => {
  table.forEach((cells, row) => cells.forEach((value, column) => {
    if (value != null) visit(value, row, column);
  }));
};

/** The elastic mode of a circular plate with n nodal diameters and s nodal circles. */
export const findMode = (modes, n, s) => modes.find(mode => mode.n === n && mode.m === s && mode.L > 0);

const rectCache = new Map();

/** Elastic modes of a rectangle, cached because the tables ask for the same plate many times. */
export const elasticRectModes = (edges, nu, beta, nb) => {
  const key = `${edges}|${nu}|${beta}|${nb}`;
  if (!rectCache.has(key)) rectCache.set(key, rectModes(edges, nu, beta, nb).modes.filter(mode => mode.L > 1e-3));
  return rectCache.get(key);
};

export { sweepResonances } from './calibrate.js';

/** Waller (1938), via Leissa NASA SP-160 Table 2.9: measured radii on free brass plates, [n][s − 1], inside out. */
export const WALLER_RADII = [
  [[.680], [.391, .843], [.257, .591, .895], [.190, .441, .692, .918]],
  [[.781], [.497, .867], [.349, .643, .902], [.269, .495, .726, .928]],
  [[.823], [.562, .887], [.415, .681, .913], [.328, .540, .748, .934]],
  [[.843], [.604, .898], [.461, .706, .919], [.374, .571, .764, .938]],
  [[.859], [.635, .906], [.505, .728, .925], [.411, .596, .779, .941]],
  [[.871], [.662, .915], [.531, .745, .933], [.443, .623, .794, .944]],
];
