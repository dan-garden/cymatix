import test from 'node:test';
import assert from 'node:assert/strict';
import { Plate } from './physics.js';
import { AIR, SPECIES, newFields, stepGrains, transportFields } from './grains.js';
import { traceFieldLines } from './fieldLines.js';

const G0 = 9.81;
const TAU = 2 * Math.PI;
const G = 81;
const UNIT = 0.12;
const SQUARE = { G, circle: false, hx: 1, hy: 1 };
const lightPlate = { h: 0.001, unit: UNIT, E: 69e9, rho: 2700, nu: 0.3 };

// The exact (2,2) mode of a simply supported square: still along its edges and along x = 0 and y = 0.
const shape = (x, y) => Math.sin(Math.PI * (x + 1)) * Math.sin(Math.PI * (y + 1));
const slopeOf = (x, y) => [
  Math.PI * Math.cos(Math.PI * (x + 1)) * Math.sin(Math.PI * (y + 1)),
  Math.PI * Math.sin(Math.PI * (x + 1)) * Math.cos(Math.PI * (y + 1)),
];

function exactMode(omega, gamma) {
  const re = new Float32Array(G * G), im = new Float32Array(G * G);
  const height = gamma * G0 / (omega * omega);
  for (let j = 0, q = 0; j < G; j++) {
    for (let i = 0; i < G; i++, q++) re[q] = height * shape((i + 0.5) * 2 / G - 1, (j + 0.5) * 2 / G - 1);
  }
  return { omega, re, im };
}

function fieldsFor(drive, grain = SPECIES.sand, fluid = null, plate = lightPlate) {
  return transportFields({ G, drives: [drive], plate, grain, fluid }, newFields(G));
}

/** Splits drawing segments back into lines, each a list of [x, y] points. */
function linesOf(vertices) {
  const lines = [];
  for (let o = 0; o < vertices.length; o += 8) {
    if (vertices[o + 2] === 0) lines.push([[vertices[o], vertices[o + 1]]]);
    lines.at(-1).push([vertices[o + 4], vertices[o + 5]]);
  }
  return lines;
}

const OMEGA = 2097; // the (2,2) mode of a 24 cm, 1 mm aluminium square

test('sand lines run down the slope of the amplitude, as the drift law u ∝ −∇|W| says', () => {
  const lines = linesOf(traceFieldLines(fieldsFor(exactMode(OMEGA, 5)), SQUARE).vertices);
  assert.ok(lines.length > 100, `${lines.length} lines`);
  for (const line of lines) {
    for (let p = 0; p + 1 < line.length; p++) {
      const [[x0, y0], [x1, y1]] = [line[p], line[p + 1]];
      if (Math.abs(shape(x0, y0)) > 0.97) continue; // on top of an antinode the slope is flat, within a cell of it the grid decides
      const size = Math.sign(shape(x0, y0));
      const [sx, sy] = slopeOf(x0, y0).map(value => -size * value); // −∇|W|
      const along = ((x1 - x0) * sx + (y1 - y0) * sy) / (Math.hypot(x1 - x0, y1 - y0) * Math.hypot(sx, sy));
      assert.ok(along > 0.95, `at (${x0}, ${y0}): cos angle ${along}`);
    }
    const heights = line.map(([x, y]) => Math.abs(shape(x, y)));
    assert.ok(heights.at(-1) < heights[0], `amplitude rose from ${heights[0]} to ${heights.at(-1)}`);
  }
});

test('sand lines stop where grains stop: where the shaking falls to about 1 g beside the still lines', () => {
  const gamma = 5;
  const lines = linesOf(traceFieldLines(fieldsFor(exactMode(OMEGA, gamma)), SQUARE).vertices);
  const oneStep = gamma * Math.PI * (2 / G); // the most the shaking changes in one step of a line
  for (const line of lines) {
    const [x, y] = line.at(-1);
    assert.ok(gamma * Math.abs(shape(x, y)) < 0.95 + oneStep, `ends at ${gamma * Math.abs(shape(x, y))} g`);
  }
});

test('a grain moved by the real sand step follows its line to the end', () => {
  const fields = fieldsFor(exactMode(OMEGA, 5));
  fields.D.fill(0);
  const line = linesOf(traceFieldLines(fields, SQUARE).vertices).sort((p, q) => q.length - p.length)[0];
  assert.ok(line.length > 10, `${line.length} points`);
  // from partway along, since right beside an antinode any tiny difference decides which way a grain leaves
  const [startX, startY] = line[Math.floor(line.length / 3)];
  const x = Float32Array.of(startX), y = Float32Array.of(startY);
  const rng = Uint32Array.of(12345);
  const cell = 2 / G;
  for (let step = 0; step < 4000; step++) {
    stepGrains(x, y, 1, fields, G, 1, 1, false, 0.25, rng);
    const closest = Math.min(...line.map(([px, py]) => Math.hypot(px - x[0], py - y[0])));
    assert.ok(closest < cell, `step ${step}: ${closest} from the line`);
  }
  const [endX, endY] = line.at(-1);
  assert.ok(Math.hypot(x[0] - endX, y[0] - endY) < 2 * cell, 'and stops at its end');
});

test('Faraday powder reverses as the drive grows: gently it follows sand, strongly its lines climb to the antinodes', () => {
  const thick = { ...lightPlate, h: 0.002 };
  const climbers = gamma => {
    const lines = linesOf(traceFieldLines(fieldsFor(exactMode(2 * OMEGA, gamma), SPECIES.lycopodium, AIR, thick), SQUARE).vertices);
    const climbing = lines.filter(line => {
      const heights = line.map(([x, y]) => Math.abs(shape(x, y)));
      return heights.at(-1) > heights[0];
    });
    for (const line of climbing) assert.ok(Math.abs(shape(...line.at(-1))) > 0.9, 'a climbing line ends at an antinode');
    return climbing.length / lines.length;
  };
  assert.equal(climbers(5), 0);
  const [strong, stronger] = [climbers(100), climbers(200)];
  assert.ok(strong > 0.4 && stronger > 0.7 && stronger > strong, `${strong} of lines climb at 100 g, ${stronger} at 200 g`);
});

test('on a real clamped drum, sand lines end beside still grains or at the clamped rim', () => {
  const plate = new Plate();
  plate.keepDrives = true;
  plate.configure({ shape: 'circle', bc: 'clamped', E: 68.9e9, rho: 2700, nu: 0.33, h: 0.001, d: 0.24, Lx: 0.24, Ly: 0.24, fmax: 8000, G, zeta: 0.01 });
  plate.setDrive(0, -0.5);
  const unitPeak = plate.field([{ f: 1433, F: 1 }], G0);
  const [drive] = plate.deflections(10 / unitPeak, G0);
  const fields = transportFields({ G, drives: [drive], plate: { ...lightPlate, nu: 0.33, E: 68.9e9 }, grain: SPECIES.sand, fluid: null, edge: plate.edge }, newFields(G));
  const lines = linesOf(traceFieldLines(fields, { G, circle: true, hx: 1, hy: 1 }).vertices);
  assert.ok(lines.length > 50, `${lines.length} lines`);
  const shakingNear = (x, y) => { // the least shaking within two cells
    const [ci, cj] = [Math.floor((x + 1) * G / 2), Math.floor((y + 1) * G / 2)];
    let least = Infinity;
    for (let j = cj - 2; j <= cj + 2; j++) {
      for (let i = ci - 2; i <= ci + 2; i++) {
        const q = j * G + i;
        const [u, v] = [(i + 0.5) * 2 / G - 1, (j + 0.5) * 2 / G - 1];
        if (u * u + v * v < 1) least = Math.min(least, Math.hypot(drive.re[q], drive.im[q]) * drive.omega ** 2 / G0);
      }
    }
    return least;
  };
  // 10 g falls about 2 g per cell here, so a still line between cells can show as 1.5 g
  const ended = lines.filter(line => {
    const [x, y] = line.at(-1);
    return x * x + y * y > 0.9 || shakingNear(x, y) < 1.5;
  });
  // the few others stop at saddles of the amplitude, where the drift itself vanishes
  assert.ok(ended.length > 0.95 * lines.length, `${ended.length} of ${lines.length}`);
  for (const [x, y] of lines.flat()) assert.ok(x * x + y * y < 1.0001);
});

test('lines stay on a narrow plate, and there are none when nothing moves', () => {
  const narrow = { G, circle: false, hx: 1, hy: 0.5 };
  for (const [x, y] of linesOf(traceFieldLines(fieldsFor(exactMode(OMEGA, 5)), narrow).vertices).flat()) {
    assert.ok(Math.abs(x) < 1 && Math.abs(y) < 0.5);
  }
  const still = traceFieldLines(fieldsFor(exactMode(OMEGA, 0.5)), SQUARE);
  assert.equal(still.topSpeed, 0, 'below 1 g grains do not hop');
  assert.equal(still.vertices.length, 0);
});

test('each vertex carries how far along its line it is and how fast grains go there', () => {
  const { vertices, topSpeed } = traceFieldLines(fieldsFor(exactMode(OMEGA, 5)), SQUARE);
  assert.ok(topSpeed > 0);
  for (let o = 0; o < vertices.length; o += 4) {
    assert.ok(vertices[o + 2] >= 0 && vertices[o + 2] <= 1 + 1e-6, `along ${vertices[o + 2]}`);
    assert.ok(vertices[o + 3] > 0 && vertices[o + 3] <= 1.5, `speed share ${vertices[o + 3]}`);
  }
});
