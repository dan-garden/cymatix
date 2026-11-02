// Draws the plate as a surface in 3D: the same image as the 2D view, carried on the true shape from plateSurface.js and
// seen through a camera from orbitCamera.js. Only drawing happens here. World axes: x = plate u, z = plate v, y = height,
// all in plate units, so the plate's outline and size are to scale.
import { GRID_VALUE, compile, makeTexture } from './sandGL.js';
import { viewProjection } from './orbitCamera.js';

const SURFACE_VS = `#version 300 es
precision highp float; precision highp int; precision highp sampler2D;
uniform sampler2D uDeflection; // R32F, metres
uniform sampler2D uCounts;     // R32F, grains per cell
uniform mat4 uViewProjection;
uniform float uDeflectionScale; // plate units per metre, magnified
uniform float uHeapScale;       // plate units per grain, true to scale
uniform int uG;
out vec2 vPlate;
out vec3 vNormal;
${GRID_VALUE}
float heightAt(vec2 p) {
  return uDeflectionScale * gridValue(uDeflection, p) + uHeapScale * gridValue(uCounts, p);
}
void main() {
  int nodes = uG + 1;
  vec2 p = vec2(gl_VertexID % nodes, gl_VertexID / nodes) / float(uG) * 2.0 - 1.0;
  float step = 2.0 / float(uG);
  vec2 dx = vec2(step, 0.0), dy = vec2(0.0, step);
  vNormal = normalize(vec3(heightAt(p - dx) - heightAt(p + dx), 2.0 * step, heightAt(p - dy) - heightAt(p + dy)));
  vPlate = p;
  gl_Position = uViewProjection * vec4(p.x, heightAt(p), p.y, 1.0);
}`;

const SURFACE_FS = `#version 300 es
precision highp float; precision highp sampler2D;
uniform sampler2D uImage; // the 2D view, bottom row at plate v = +1
uniform sampler2D uMarks; // driver and supports, top row at plate v = −1
uniform bool uCircle;
uniform vec2 uHalf;
uniform vec3 uBackground;
in vec2 vPlate;
in vec3 vNormal;
out vec4 outColor;
const vec3 LIGHT = normalize(vec3(-0.45, 0.8, 0.4));
void main() {
  float toEdge = uCircle ? 1.0 - length(vPlate) : min(uHalf.x - abs(vPlate.x), uHalf.y - abs(vPlate.y));
  float coverage = clamp(toEdge / fwidth(toEdge) + 0.5, 0.0, 1.0);
  if (coverage <= 0.0) discard;
  vec2 uv = (vPlate + 1.0) / 2.0;
  vec3 colour = texture(uImage, vec2(uv.x, 1.0 - uv.y)).rgb;
  vec4 marks = texture(uMarks, uv);
  colour = mix(colour, marks.rgb, marks.a);
  // relief only: a flat plate keeps exactly the colours of the 2D view
  float relief = 1.0 + 1.1 * (dot(normalize(vNormal), LIGHT) - LIGHT.y);
  outColor = vec4(mix(uBackground, colour * max(relief, 0.15), coverage), 1.0);
}`;

/** Returns a view that draws into the canvas of the given WebGL2 context, on request. */
export function createPlateView3D(gl) {
  const surface = compile(gl, SURFACE_VS, SURFACE_FS);
  const vertexArray = gl.createVertexArray();
  const indexBuffer = gl.createBuffer();
  let G = 0, indexCount = 0, deflectionTexture = null, marksTexture = null;

  function buildGrid(grid) {
    G = grid;
    const nodes = G + 1;
    const indices = new Uint32Array(6 * G * G);
    for (let j = 0, o = 0; j < G; j++) {
      for (let i = 0; i < G; i++) {
        const a = j * nodes + i, b = a + 1, c = a + nodes, d = c + 1;
        indices.set([a, c, b, b, c, d], o);
        o += 6;
      }
    }
    indexCount = indices.length;
    gl.bindVertexArray(vertexArray);
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indexBuffer);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, indices, gl.STATIC_DRAW);
    if (deflectionTexture) gl.deleteTexture(deflectionTexture);
    deflectionTexture = makeTexture(gl, gl.R32F, gl.RED, gl.FLOAT, G, G);
  }

  const bind = (unit, texture, name) => {
    gl.activeTexture(gl.TEXTURE0 + unit);
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.uniform1i(surface.uniforms.get(name), unit);
  };

  return {
    /** Take the supports and driver marker from a transparent 2D canvas. */
    setMarks(canvas) {
      if (marksTexture) gl.deleteTexture(marksTexture);
      marksTexture = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, marksTexture);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, canvas);
    },

    /**
     * Draw one frame. deflection: G×G Float32Array in metres. Scales turn metres and grains per cell into plate units.
     * image and counts: textures from the sand engine. background: [r, g, b] in 0–1.
     */
    draw({ camera, grid, deflection, deflectionScale, heapScale, image, counts, circle, hx, hy, background, size }) {
      const sandVertexArray = gl.getParameter(gl.VERTEX_ARRAY_BINDING);
      if (grid !== G) buildGrid(grid);
      gl.bindVertexArray(vertexArray);
      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, deflectionTexture);
      gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, G, G, gl.RED, gl.FLOAT, deflection);

      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      gl.viewport(0, 0, size, size);
      gl.disable(gl.BLEND);
      gl.enable(gl.DEPTH_TEST);
      gl.clearColor(...background, 1);
      gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

      gl.useProgram(surface.program);
      const u = surface.uniforms;
      bind(0, deflectionTexture, 'uDeflection');
      bind(1, counts, 'uCounts');
      bind(2, image, 'uImage');
      bind(3, marksTexture, 'uMarks');
      gl.uniformMatrix4fv(u.get('uViewProjection'), false, viewProjection(camera));
      gl.uniform1f(u.get('uDeflectionScale'), deflectionScale);
      gl.uniform1f(u.get('uHeapScale'), heapScale);
      gl.uniform1i(u.get('uG'), G);
      gl.uniform1i(u.get('uCircle'), circle ? 1 : 0);
      gl.uniform2f(u.get('uHalf'), hx, hy);
      gl.uniform3f(u.get('uBackground'), ...background);
      gl.drawElements(gl.TRIANGLES, indexCount, gl.UNSIGNED_INT, 0);

      gl.disable(gl.DEPTH_TEST);
      gl.bindVertexArray(sandVertexArray);
    },
  };
}
