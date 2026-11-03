// Field lines of the grains' drift: the path a grain follows on average, from the fields that really move them. Pure maths,
// with no drawing here. A line runs from where the grains start to where they stop. Its spacing means nothing, unlike the
// field lines of a magnet, because the lines are simply started from points spread evenly over the plate.
import { activityAt, driftAt } from './grains.js';

const STOP_SHARE = 0.005; // a line ends where drift falls below this share of the fastest on the plate
const SPREAD_X = 0.7548776662, SPREAD_Y = 0.5698402910; // an even, repeatable scatter (the R2 sequence), so lines hold still
const FLOATS_PER_VERTEX = 4;

let buffer = new Float32Array(0);
let trail = new Float32Array(0);
const velocity = [0, 0];

function fastestDrift({ ux, uy, level }) {
  let fastest = 0;
  for (let q = 0; q < ux.length; q++) fastest = Math.max(fastest, activityAt(level[q]) * Math.hypot(ux[q], uy[q]));
  return fastest;
}

/**
 * Traces up to `lines` lines, each from a point on the plate both ways, one grid cell per step by the midpoint rule: forward
 * until it leaves the plate, reaches still grains, or crosses a line where the drift reverses, and back to where its grains
 * come from. Returns segments as pairs of [x, y, along, speed] vertices, where along runs 0–1 from start to end, speed
 * is a share of the fastest, and topSpeed is that fastest drift in plate units/s.
 */
export function traceFieldLines(fields, { G, circle, hx, hy }, lines = 400, reach = 1) {
  const topSpeed = fastestDrift(fields);
  if (topSpeed === 0) return { vertices: buffer.subarray(0, 0), topSpeed };
  const cell = 2 / G;
  const maxSteps = Math.ceil(reach / cell);
  const needed = lines * 2 * maxSteps * 2 * FLOATS_PER_VERTEX;
  if (buffer.length < needed) buffer = new Float32Array(needed);
  if (trail.length < 3 * (2 * maxSteps + 1)) trail = new Float32Array(3 * (2 * maxSteps + 1));
  const minSpeed = STOP_SHARE * topSpeed;
  const inside = (x, y) => (circle ? x * x + y * y < 1 : Math.abs(x) < hx && Math.abs(y) < hy);
  const speedAt = (x, y) => {
    driftAt(fields, G, x, y, velocity);
    return Math.hypot(velocity[0], velocity[1]);
  };
  const at = { x: 0, y: 0, dx: 0, dy: 0, speed: 0 };
  const settle = (x, y) => {
    at.x = x;
    at.y = y;
    at.speed = speedAt(x, y);
    at.dx = velocity[0] / at.speed;
    at.dy = velocity[1] / at.speed;
  };
  // One step forward (sign 1) or back (sign −1), or false if the line ends here.
  const step = sign => {
    const middle = speedAt(at.x + sign * at.dx * cell / 2, at.y + sign * at.dy * cell / 2);
    if (middle < minSpeed) return false;
    const stepX = velocity[0] / middle, stepY = velocity[1] / middle;
    if (stepX * at.dx + stepY * at.dy <= 0) return false;
    const nextX = at.x + sign * stepX * cell, nextY = at.y + sign * stepY * cell;
    if (!inside(nextX, nextY) || speedAt(nextX, nextY) < minSpeed) return false;
    settle(nextX, nextY);
    return true;
  };
  const record = point => {
    trail[3 * point] = at.x;
    trail[3 * point + 1] = at.y;
    trail[3 * point + 2] = at.speed / topSpeed;
  };
  const emit = (point, along, o) => {
    buffer[o] = trail[3 * point];
    buffer[o + 1] = trail[3 * point + 1];
    buffer[o + 2] = along;
    buffer[o + 3] = trail[3 * point + 2];
  };

  let used = 0;
  for (let line = 0; line < lines; line++) {
    const startX = (((0.5 + line * SPREAD_X) % 1) * 2 - 1) * hx;
    const startY = (((0.5 + line * SPREAD_Y) % 1) * 2 - 1) * hy;
    if (!inside(startX, startY) || speedAt(startX, startY) < minSpeed) continue;
    let first = maxSteps, last = maxSteps;
    settle(startX, startY);
    record(maxSteps);
    while (maxSteps - first < maxSteps && step(-1)) record(--first);
    settle(startX, startY);
    while (last - maxSteps < maxSteps && step(1)) record(++last);
    for (let point = first; point < last; point++) {
      emit(point, (point - first) / (last - first), used);
      emit(point + 1, (point + 1 - first) / (last - first), used + FLOATS_PER_VERTEX);
      used += 2 * FLOATS_PER_VERTEX;
    }
  }
  return { vertices: buffer.subarray(0, used), topSpeed };
}
