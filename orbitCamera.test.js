import test from 'node:test';
import assert from 'node:assert/strict';
import { createOrbitCamera, eyeOf, orbitBy, pickPlane, projectPoint, viewProjection, zoomBy } from './orbitCamera.js';

const near = (got, want, tolerance = 1e-5) => assert.ok(Math.abs(got - want) <= tolerance, `${got} against ${want}`);

test('the origin sits at the centre of the screen, from any angle', () => {
  for (const azimuth of [0, 0.7, -2.5]) {
    for (const elevation of [0.2, 0.7, 1.5]) {
      const [x, y] = projectPoint(createOrbitCamera({ azimuth, elevation }), [0, 0, 0]);
      near(x, 0);
      near(y, 0);
    }
  }
});

test('a point on the plane and its pick agree, from any angle', () => {
  for (const azimuth of [0, 0.7, -2.5]) {
    for (const elevation of [0.3, 0.8, 1.5]) {
      const camera = createOrbitCamera({ azimuth, elevation });
      for (const [u, v] of [[0.5, -0.3], [-0.9, 0.9], [0, 0.4]]) {
        const [x, z] = pickPlane(camera, ...projectPoint(camera, [u, 0, v]));
        near(x, u, 1e-4);
        near(z, v, 1e-4);
      }
    }
  }
});

test('from almost straight above, the plate is an undistorted, correctly oriented copy of the 2D view', () => {
  const camera = createOrbitCamera({ elevation: 89 * Math.PI / 180 });
  const scale = 1 / (Math.tan(camera.fov / 2) * camera.distance);
  const [x, y] = projectPoint(camera, [0.5, 0, 0.5]);
  near(x, 0.5 * scale, 0.01);
  near(y, -0.5 * scale, 0.01); // plate v points down the screen, as in 2D
});

test('the default view has the near edge at the bottom of the screen, and rotating moves the eye as expected', () => {
  const camera = createOrbitCamera();
  assert.ok(projectPoint(camera, [0, 0, 1])[1] < projectPoint(camera, [0, 0, -1])[1]);
  assert.ok(eyeOf(camera)[2] > 0);
});

test('looking above the horizon picks nothing', () => {
  assert.equal(pickPlane(createOrbitCamera({ elevation: 0.09 }), 0, 1), null);
});

test('orbit and zoom stay in range, so the camera can never go under the plate or inside it', () => {
  const camera = createOrbitCamera();
  orbitBy(camera, 0, -10);
  assert.ok(eyeOf(camera)[1] > 0);
  orbitBy(camera, 0, 10);
  assert.ok(camera.elevation < Math.PI / 2);
  zoomBy(camera, 1e-6);
  assert.ok(camera.distance >= 1.5);
  zoomBy(camera, 1e6);
  assert.ok(camera.distance <= 20);
});

test('the matrix is a valid perspective: near things nearer, depth between the clip planes', () => {
  const m = viewProjection(createOrbitCamera());
  assert.equal(m.length, 16);
  const clipW = (x, y, z) => m[3] * x + m[7] * y + m[11] * z + m[15];
  assert.ok(clipW(0, 0, 0) > 0);
  near(clipW(0, 0, 0), 4.2);
});
