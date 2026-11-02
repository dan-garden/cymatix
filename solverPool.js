// A small pool of solver threads. Rectangles split into independent parity blocks, so they solve in parallel.
const THREADS = Math.max(1, Math.min(4, (navigator.hardwareConcurrency ?? 2) - 1));

const idle = [];
const busy = new Map(); // worker → the group its task belongs to
let queue = [];
let nextId = 1;

function startWorker() {
  return new Worker(new URL('./solver.worker.js', import.meta.url), { type: 'module' });
}

function dispatch() {
  while (queue.length && (idle.length || busy.size < THREADS)) {
    const worker = idle.pop() ?? startWorker();
    const { task, args, group, resolve, reject } = queue.shift();
    const id = nextId++;
    busy.set(worker, group);
    worker.onmessage = ({ data }) => {
      if (data.id !== id) return;
      busy.delete(worker);
      idle.push(worker);
      resolve(data.result);
      dispatch();
    };
    worker.onerror = error => {
      busy.delete(worker); // a crashed worker is not reused
      reject(error);
      dispatch();
    };
    worker.postMessage({ id, task, args });
  }
}

/** Runs a task in a worker thread. `group` names tasks that cancelGroup may abandon together. */
export const runInGroup = (group, task, ...args) => new Promise((resolve, reject) => {
  queue.push({ task, args, group, resolve, reject });
  dispatch();
});

export const runTask = (task, ...args) => runInGroup(null, task, ...args);

/** Abandons a group's queued and running tasks: their promises never settle, and their threads are freed at once. */
export function cancelGroup(group) {
  queue = queue.filter(entry => entry.group !== group);
  for (const [worker, workerGroup] of busy) {
    if (workerGroup !== group) continue;
    worker.terminate();
    busy.delete(worker);
  }
  dispatch();
}
