// Decides which solve a plate deserves. The main thread answers at once with thin-plate modes, and anything better or
// larger is planned here and run in worker threads (or directly, in tests), then replaces them.
import { circleLamFor, circleModes, rectBasisFor, rectBlock, rectBlockCount, rectModes, rectModesFromBlocks } from './physics.js';
import { MINDLIN_BASIS, mindlinCircleModes, mindlinRectEvaluators, mindlinRectModes } from './mindlin.js';
import { rectKappas } from './nonlinear.js';
import { PlateModel } from './plateModel.js';

// Below these thickness ratios, the per-mode correction is as close as the plate's own inputs (0.3–0.5 % at 0.01).
const THICK_CIRCLE = 0.005;
const THICK_RECT = 0.01;
const LAM_INLINE = 50;
const MINDLIN_MAX_BASIS = [22, 28, 34]; // by symmetric axes: about 8 s at the top

function circlePlan(cfg) {
  const { bc, nu, h, d, mindlin } = cfg;
  const lamMax = circleLamFor(cfg);
  const hOverA = h / (d / 2);
  if (mindlin !== false && hOverA >= THICK_CIRCLE) {
    return { id: `m${hOverA}|${lamMax}`, thick: true, lamMax, task: 'mindlinCircleModes', args: [bc, nu, hOverA, lamMax] };
  }
  if (lamMax <= LAM_INLINE) return null;
  return { id: `t${lamMax}`, lamMax, task: 'circleModes', args: [bc, nu, lamMax] };
}

function rectPlan(cfg) {
  const { bc, nu, h, Lx, Ly, mindlin } = cfg;
  const edges = cfg.edges ?? [bc, bc, bc, bc];
  const beta = Ly / Lx;
  const wanted = rectBasisFor(cfg);
  if (mindlin !== false && h / Math.min(Lx, Ly) >= THICK_RECT) {
    const symAxes = (edges[0] === edges[1]) + (edges[2] === edges[3]);
    const nb = Math.min(MINDLIN_MAX_BASIS[symAxes], Math.max(wanted ?? 0, MINDLIN_BASIS[symAxes]));
    return { id: `m${h / Lx}|${nb}`, thick: true, task: 'mindlinRectData', args: [edges, nu, beta, h / Lx, nb], edges, beta, nb };
  }
  if (!wanted) return null;
  return { id: `t${wanted}`, task: 'rectBlock', args: [edges, nu, beta, wanted], edges, nu, beta, nb: wanted };
}

/** The solve a plate should have, or null if the main thread's own answer is already right. */
export const solvePlan = cfg => (cfg.shape === 'circle' ? circlePlan(cfg) : rectPlan(cfg));

/** Solves that a worker thread can run: arguments and results must cross between threads. */
export const TASKS = {
  circleModes,
  mindlinCircleModes,
  rectBlock,
  rectKappasDefault(edges, nu, beta) {
    return rectKappas(rectModes(edges, nu, beta), edges, nu, beta).kappa;
  },
  async physicalDamping(cfg) {
    const plate = new PlateModel();
    plate.configure({ ...cfg, solved: await solveNow(cfg) });
    return plate.physicalDamping(cfg);
  },
  mindlinRectData(edges, nu, beta, hOverA, nb) {
    const { modes, nbx, nby } = mindlinRectModes(edges, nu, beta, hOverA, { nb });
    return { modes, nbx, nby };
  },
};

/** Runs a plan through `run(task, ...args)`, which may be a worker pool, and returns cfg.solved for Plate.configure. */
export async function runPlan(plan, run) {
  const { id, thick, lamMax, task, args } = plan;
  if (task === 'rectBlock') {
    const { edges, nu, beta, nb } = plan;
    const blocks = Array.from({ length: rectBlockCount(edges, nu, beta, nb) }, (_, block) => run('rectBlock', ...args, block));
    return { id, ritz: rectModesFromBlocks(edges, nu, beta, nb, await Promise.all(blocks)) };
  }
  const result = await run(task, ...args);
  if (task !== 'mindlinRectData') return { id, thick, lamMax, modes: result };
  return { id, thick, ritz: { ...result, ...mindlinRectEvaluators(plan.edges, plan.beta, plan.nb) } };
}

/** The same solve on this thread, for tests and scripts. */
export const solveNow = cfg => {
  const plan = solvePlan(cfg);
  return plan && runPlan(plan, (task, ...args) => TASKS[task](...args));
};
