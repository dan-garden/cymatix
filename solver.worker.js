// Runs the slow solves off the main thread, so the page stays responsive while they finish.
import { fitPlate } from './calibrate.js';
import { TASKS } from './solvePlan.js';

const ALL_TASKS = { ...TASKS, fitPlate };

// Typed arrays in the result are handed over rather than copied.
const buffersIn = result => Object.values(result ?? {})
  .filter(value => ArrayBuffer.isView(value))
  .map(value => value.buffer);

self.onmessage = async ({ data: { id, task, args } }) => {
  const result = await ALL_TASKS[task](...args);
  self.postMessage({ id, result }, buffersIn(result));
};
