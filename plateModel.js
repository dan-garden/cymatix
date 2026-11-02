// The full plate: the core modal model plus the effects that live in their own modules. Each one is optional in the
// configuration, and with all of them absent this is exactly Plate.
import { Plate } from './physics.js';
import { ALUMINIUM_THERMAL, circleAddedMass, loadModes, modalDamping, rectAddedMass } from './fluid.js';
import { RigidPatch } from './patch.js';
import { circleKappas, duffingShift, regime } from './nonlinear.js';

// Exact coupled fluid loading is solved for this many modes. Higher ones use the infinite-plate limit, Γ → 2/k, which
// fluid.test.js confirms, because the exact solve grows too slow for the main thread and they barely couple.
const EXACT_FLUID_MODES = 400;
const FLUID_BASIS_LIMIT = 36;
const MAX_DRIVES = 24;

const loadedByLocalLimit = (mode, ratio) => {
  const inertia = 1 + ratio * 2 / Math.sqrt(mode.L);
  const w = mode.w / Math.sqrt(inertia);
  return { ...mode, w, f: w / (2 * Math.PI), inv: mode.inv / inertia };
};

export class PlateModel extends Plate {
  #added = null;
  #addedFor = null; // the solved modes #added and #kappas belong to
  #kappas = null;
  #kappasFor = null;
  #shifts = null;

  // Extra cfg, all optional: fluid { rho, c, mu } around the plate; clampRadius (m) of a rigid bolt or shaker clamp;
  // stiffening, for von Kármán large-amplitude stiffening on a force mount; modalZeta, a per-mode damping ratio.
  configure(cfg) {
    this.patch = null;
    this.detune = null;
    this.h = cfg.h;
    super.configure({ ...cfg, adjustModes: modes => this.#adjust(modes, cfg) });
    const { mount, clampRadius, ma } = cfg;
    if (clampRadius > 0 && (mount === 'rigid' || mount === 'shaker')) {
      this.patch = new RigidPatch(this, cfg, { r0: clampRadius, mount, ma });
    }
    this.stiffeningOn = !!cfg.stiffening && this.mount === 'force';
    if (this.#kappasFor === this.solvedModes) this.#armStiffening();
  }

  #adjust(modes, cfg) {
    const inFluid = this.#inFluid(modes, cfg);
    const { modalZeta } = cfg;
    if (!modalZeta) return inFluid;
    return inFluid.map((mode, k) => ({ ...mode, zeta: modalZeta[k] }));
  }

  #inFluid(modes, cfg) {
    if (!cfg.fluid?.rho) return modes;
    const exact = modes.slice(0, EXACT_FLUID_MODES);
    if (this.#addedFor !== modes) {
      this.#added = this.circle ? circleAddedMass(exact) : this.#rectAddedMass(exact, cfg);
      this.#addedFor = modes;
    }
    const length = this.circle ? cfg.d / 2 : cfg.Lx;
    const ratio = cfg.fluid.rho * length / (cfg.rho * cfg.h);
    const loaded = loadModes(exact, this.#added, ratio);
    if (modes.length <= exact.length) return loaded;
    return loaded.concat(modes.slice(exact.length).map(mode => loadedByLocalLimit(mode, ratio)));
  }

  #rectAddedMass(exact, cfg) {
    const { ritz } = this;
    return rectAddedMass({ ...ritz, modes: exact }, cfg.Ly / cfg.Lx, {
      Px: Math.min(FLUID_BASIS_LIMIT, ritz.nbx),
      Py: Math.min(FLUID_BASIS_LIMIT, ritz.nby),
    });
  }

  /**
   * Per-mode damping from thermoelastic, viscous and radiation losses, for aluminium (the one material with sourced
   * thermal constants). A real rig adds losses at the bolt, the shaker and the sand, so this is a lower bound. Slow for
   * rectangles: run it in a worker thread.
   */
  physicalDamping(cfg) {
    const { fluid, h, rho, Ly, Lx, d, edges } = cfg;
    return modalDamping(this.solvedModes, this.modes, this.#added, {
      h,
      rhoH: rho * h,
      L: this.circle ? d / 2 : Lx,
      fluid,
      thermal: ALUMINIUM_THERMAL,
      ritz: this.ritz,
      beta: Ly / Lx,
      edges,
    });
  }

  /** Circular plates compute their stiffening coefficients at once; rectangles bring theirs from a worker thread. */
  prepareStiffening(cfg) {
    if (this.circle) this.setKappas(circleKappas(this.solvedModes, cfg.bc, cfg.nu));
  }

  setKappas(kappas) {
    this.#kappas = kappas;
    this.#kappasFor = this.solvedModes;
    this.#armStiffening();
  }

  get hasKappas() {
    return this.#kappasFor === this.solvedModes;
  }

  #armStiffening() {
    this.detune = null;
    if (!this.stiffeningOn) return;
    const { modes, wScale, h } = this;
    const gammas = Float64Array.from(modes, (mode, k) => this.#kappas[k] * (mode.L * wScale) ** 2 / (h * h));
    if (this.#shifts?.length !== modes.length * MAX_DRIVES) this.#shifts = new Float64Array(modes.length * MAX_DRIVES);
    const shifts = this.#shifts;
    this.forceGain = this.forceGain ?? 1;
    this.detune = (k, drive, detuning, damping, forcePerMass) => {
      const f = forcePerMass * this.forceGain;
      const slot = k * MAX_DRIVES + Math.min(drive, MAX_DRIVES) - 1;
      const shift = duffingShift(detuning, damping, 0.75 * gammas[k] * f * f, shifts[slot]);
      shifts[slot] = shift;
      return shift;
    };
  }
}

/**
 * How far large-amplitude stiffening has gone, from the plate acceleration a (m/s²) at ω, as the regime in nonlinear.js.
 * κ_A = 1 is a cautious value: over all modes of a free circle it runs 0.35–1.09, and on a clamped square about 1.
 */
export function stiffeningRegime(a, omega, h, zeta, kappaA = 1) {
  const x = a / (omega * omega * h);
  return regime(0.375 * kappaA * x * x / zeta);
}
