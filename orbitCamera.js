// An orbiting camera for any 3D scene, looking at a point on the plane y = 0 with y up. Pure maths, plus optional pointer and key controls.

const DEGREES = Math.PI / 180;
const MIN_ELEVATION = 5 * DEGREES;
const MAX_ELEVATION = 89 * DEGREES; // straight down has no defined "up"
export const HOME_DISTANCE = 4.2;
const MIN_DISTANCE = 0.2; // near enough to see single grains
const MAX_DISTANCE = 20;
const KEY_STEP = 3 * DEGREES;
const CLICK_SLOP = 4; // CSS pixels a pointer may move and still count as a click

const cross = ([a, b, c], [d, e, f]) => [b * f - c * e, c * d - a * f, a * e - b * d];
const dot = (p, q) => p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
const normalize = p => {
  const length = Math.hypot(...p);
  return p.map(value => value / length);
};

export const createOrbitCamera = ({ azimuth = 0, elevation = 40 * DEGREES, distance = HOME_DISTANCE, fov = 35 * DEGREES } = {}) => ({
  azimuth,
  elevation,
  distance,
  fov,
  target: [0, 0], // the [x, z] looked at
});

const toEye = ({ azimuth, elevation, distance }) => [
  distance * Math.cos(elevation) * Math.sin(azimuth),
  distance * Math.sin(elevation),
  distance * Math.cos(elevation) * Math.cos(azimuth),
];

export const eyeOf = camera => {
  const [x, y, z] = toEye(camera);
  return [x + camera.target[0], y, z + camera.target[1]];
};

/** The camera's right, up and forward directions. */
export function axesOf(camera) {
  const forward = normalize(toEye(camera).map(value => -value));
  const right = normalize(cross(forward, [0, 1, 0]));
  return { right, up: cross(right, forward), forward };
}

/** The view-projection matrix as 16 numbers, column by column, ready for a shader. */
export function viewProjection(camera, aspect = 1, near = 0.02, far = 50) {
  const { right, up, forward } = axesOf(camera);
  const eye = eyeOf(camera);
  const f = 1 / Math.tan(camera.fov / 2);
  const depthA = (far + near) / (near - far);
  const depthB = 2 * far * near / (near - far);
  const row = axis => [axis[0], axis[1], axis[2], -dot(axis, eye)];
  const [r, u, b] = [row(right), row(up), row(forward.map(value => -value))];
  return Float32Array.from({ length: 16 }, (_, index) => {
    const column = index >> 2;
    switch (index & 3) {
      case 0: return f / aspect * r[column];
      case 1: return f * u[column];
      case 2: return depthA * b[column] + (column === 3 ? depthB : 0);
      default: return -b[column];
    }
  });
}

/** Where a world point lands on screen, as x and y in [-1, 1] with y up. */
export function projectPoint(camera, [x, y, z], aspect = 1) {
  const m = viewProjection(camera, aspect);
  const w = m[3] * x + m[7] * y + m[11] * z + m[15];
  return [(m[0] * x + m[4] * y + m[8] * z + m[12]) / w, (m[1] * x + m[5] * y + m[9] * z + m[13]) / w];
}

/** The point [x, z] on the plane y = 0 under a screen position, or null if that position looks above the horizon. */
export function pickPlane(camera, screenX, screenY, aspect = 1) {
  const { right, up, forward } = axesOf(camera);
  const spread = Math.tan(camera.fov / 2);
  const ray = forward.map((value, i) => value + right[i] * screenX * spread * aspect + up[i] * screenY * spread);
  if (ray[1] >= 0) return null;
  const eye = eyeOf(camera);
  const reach = -eye[1] / ray[1];
  return [eye[0] + ray[0] * reach, eye[2] + ray[2] * reach];
}

export function orbitBy(camera, azimuthChange, elevationChange) {
  camera.azimuth += azimuthChange;
  camera.elevation = Math.min(MAX_ELEVATION, Math.max(MIN_ELEVATION, camera.elevation + elevationChange));
}

export function zoomBy(camera, factor) {
  camera.distance = Math.min(MAX_DISTANCE, Math.max(MIN_DISTANCE, camera.distance * factor));
}

/**
 * Zoom keeping the plane point under a screen position where it is, so any part of the plate can be looked at closely.
 * The target may only leave the origin as far as the zoom has come in, so zooming back out always returns home.
 */
export function zoomAt(camera, factor, screenX, screenY) {
  const before = pickPlane(camera, screenX, screenY);
  zoomBy(camera, factor);
  const after = pickPlane(camera, screenX, screenY);
  const room = Math.max(0, 1 - camera.distance / HOME_DISTANCE);
  camera.target = camera.target.map((value, axis) => {
    const moved = before && after ? value + before[axis] - after[axis] : value;
    return Math.min(room, Math.max(-room, moved));
  });
}

/**
 * Drag to orbit, scroll to zoom towards the pointer, arrow keys to orbit, and a click without a drag reports its screen position.
 * Controls do nothing while isActive() is false, so a canvas can share them with another view.
 */
export function orbitControls(canvas, camera, { isActive, onClick }) {
  let drag = null;
  const screenPosition = e => {
    const box = canvas.getBoundingClientRect();
    return [(e.clientX - box.left) / box.width * 2 - 1, 1 - (e.clientY - box.top) / box.height * 2];
  };
  canvas.addEventListener('pointerdown', e => {
    if (!isActive()) return;
    canvas.setPointerCapture(e.pointerId);
    drag = { startX: e.clientX, startY: e.clientY, x: e.clientX, y: e.clientY };
  });
  canvas.addEventListener('pointermove', e => {
    if (!drag) return;
    orbitBy(camera, -(e.clientX - drag.x) * 0.3 * DEGREES, (e.clientY - drag.y) * 0.3 * DEGREES);
    drag.x = e.clientX;
    drag.y = e.clientY;
  });
  canvas.addEventListener('pointerup', e => {
    if (!drag) return;
    const moved = Math.hypot(e.clientX - drag.startX, e.clientY - drag.startY);
    drag = null;
    if (moved <= CLICK_SLOP) onClick?.(...screenPosition(e));
  });
  canvas.addEventListener('pointercancel', () => { drag = null; });
  canvas.addEventListener('wheel', e => {
    if (!isActive()) return;
    e.preventDefault();
    zoomAt(camera, Math.exp(e.deltaY * 0.001), ...screenPosition(e));
  }, { passive: false });
  const KEYS = {
    ArrowLeft: [KEY_STEP, 0],
    ArrowRight: [-KEY_STEP, 0],
    ArrowUp: [0, KEY_STEP],
    ArrowDown: [0, -KEY_STEP],
  };
  canvas.addEventListener('keydown', e => {
    if (!isActive() || !KEYS[e.key]) return;
    e.preventDefault();
    orbitBy(camera, ...KEYS[e.key]);
  });
}
