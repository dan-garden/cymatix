// Draws the plate in 3D: a flat plate carrying the 2D view's colours, every grain as a true-size ball where the simulation has
// it, and the grains' field lines rising and falling as the plate moves, magnified. Only drawing happens here.
// World axes: x = plate u, z = plate v, y = height, all in plate units, so the plate and grains are to scale.
import { GRID_VALUE, compile, makeTexture } from './sandGL.js';
import { HOME_DISTANCE, axesOf, viewProjection } from './orbitCamera.js';

const SURFACE_VS = `#version 300 es
precision highp float;
uniform mat4 uViewProjection;
out vec2 vPlate;
void main() {
  vPlate = vec2(gl_VertexID & 1, gl_VertexID >> 1) * 2.0 - 1.0;
  gl_Position = uViewProjection * vec4(vPlate.x, 0.0, vPlate.y, 1.0);
}`;

const SURFACE_FS = `#version 300 es
precision highp float; precision highp sampler2D;
uniform sampler2D uImage; // the 2D view without grains, bottom row at plate v = +1
uniform sampler2D uMarks; // driver and supports, top row at plate v = −1
uniform vec3 uBackground;
uniform bool uCircle;
uniform vec2 uHalf;
in vec2 vPlate;
out vec4 outColor;
void main() {
  float toEdge = uCircle ? 1.0 - length(vPlate) : min(uHalf.x - abs(vPlate.x), uHalf.y - abs(vPlate.y));
  float coverage = clamp(toEdge / fwidth(toEdge) + 0.5, 0.0, 1.0);
  if (coverage <= 0.0) discard;
  vec2 uv = (vPlate + 1.0) / 2.0;
  vec3 colour = texture(uImage, vec2(uv.x, 1.0 - uv.y)).rgb;
  vec4 marks = texture(uMarks, uv);
  outColor = vec4(mix(uBackground, mix(colour, marks.rgb, marks.a), coverage), 1.0);
}`;

const GRAIN_VS = `#version 300 es
precision highp float; precision highp int; precision highp sampler2D;
uniform sampler2D uPositions; // RGBA32F, plate coordinates
uniform sampler2D uCounts;    // R32F, grains per cell
uniform mat4 uViewProjection;
uniform int uSide;
uniform int uG;
uniform float uHeapScale;   // plate units of pile per grain in a cell, true to scale
uniform float uDiameter;    // plate units
uniform float uPixelsAtOne; // pixels per plate unit, one plate unit from the eye
${GRID_VALUE}
float shareOf(uint grain) { // a fixed number 0–1 for each grain (the PCG hash, Jarzynski & Olano 2020)
  uint state = grain * 747796405u + 2891336453u;
  uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
  return float((word >> 22u) ^ word) / 4294967296.0;
}
void main() {
  vec2 p = texelFetch(uPositions, ivec2(gl_VertexID % uSide, gl_VertexID / uSide), 0).xy;
  // shortcut: a grain sits at a fixed share of its pile's height, from the volume the cell's grains fill; true stacking needs grain contact
  float pile = uHeapScale * gridValue(uCounts, p);
  float height = 0.5 * uDiameter + shareOf(uint(gl_VertexID)) * max(pile - uDiameter, 0.0);
  gl_Position = uViewProjection * vec4(p.x, height, p.y, 1.0);
  gl_PointSize = max(uDiameter * uPixelsAtOne / gl_Position.w, 1.0); // a grain under a pixel still shows as a pixel, as in 2D
}`;

const GRAIN_FS = `#version 300 es
precision highp float;
uniform vec3 uSand;
uniform mat3 uToWorld; // from the sprite's right, up and out of the screen to world directions
out vec4 outColor;
const vec3 LIGHT = normalize(vec3(-0.45, 0.8, 0.4));
void main() {
  vec2 c = gl_PointCoord * 2.0 - 1.0;
  float out2 = dot(c, c);
  if (out2 > 1.0) discard;
  vec3 normal = uToWorld * vec3(c.x, -c.y, sqrt(1.0 - out2));
  outColor = vec4(uSand * (0.4 + 0.8 * max(dot(normal, LIGHT), 0.0)), 1.0);
}`;

const LINE_VS = `#version 300 es
precision highp float; precision highp int; precision highp sampler2D;
layout(location = 0) in vec4 aFrom; // one segment of a line, drawn as a ribbon: x, z, how far along its line 0–1, speed as a share of the fastest
layout(location = 1) in vec4 aTo;
uniform sampler2D uHeights; // R32F, the plate's height now in metres
uniform mat4 uViewProjection;
uniform float uScale; // plate units per metre, magnified
uniform float uWidthAtOne; // ribbon width in pixels, one plate unit from the eye, so nearer lines are wider
uniform float uSize; // canvas side in pixels
uniform int uG;
out float vStrength;
out float vHeight;
out float vAcross;
const float CLOSEST = 0.03; // a segment nearer the eye than this is dropped, since it has no sound position on screen
${GRID_VALUE}
vec4 clipOf(vec4 point) {
  return uViewProjection * vec4(point.x, uScale * gridValue(uHeights, point.xy), point.y, 1.0);
}
void main() {
  float along = float(gl_VertexID & 1);
  vAcross = float(gl_VertexID >> 1) * 2.0 - 1.0;
  vec4 point = mix(aFrom, aTo, along);
  vHeight = uScale * gridValue(uHeights, point.xy);
  // each line brightens the way the grains travel, faster drift is brighter, and it is fainter where it dips below the plate
  vStrength = (0.3 + 0.7 * point.z) * (0.45 + 0.55 * sqrt(point.w)) * (vHeight < 0.0 ? 0.6 : 1.0);
  vec4 from = clipOf(aFrom), to = clipOf(aTo);
  gl_Position = mix(from, to, along);
  if (min(from.w, to.w) < CLOSEST) {
    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
    return;
  }
  vec2 onScreen = to.xy / to.w - from.xy / from.w;
  vec2 across = vec2(-onScreen.y, onScreen.x) / max(length(onScreen), 1e-6);
  gl_Position.xy += across * vAcross * clamp(uWidthAtOne / gl_Position.w, 1.0, 6.0) / uSize * gl_Position.w;
}`;

const LINE_FS = `#version 300 es
precision highp float;
in float vStrength;
in float vHeight;
in float vAcross;
out vec4 outColor;
const vec3 LINE_COLOUR = vec3(0.45, 0.95, 1.0);
const float TALLEST = 0.3; // plate units, about the most Auto draws
void main() {
  // crests lighter, so the shape reads even from straight above
  vec3 colour = mix(LINE_COLOUR, vec3(1.0), 0.7 * clamp(vHeight / TALLEST, 0.0, 1.0));
  outColor = vec4(colour, vStrength * (1.0 - smoothstep(0.4, 1.0, abs(vAcross))));
}`;

const FLOATS_PER_LINE_VERTEX = 4;
const BYTES_PER_LINE_VERTEX = 4 * FLOATS_PER_LINE_VERTEX;
const LINE_WIDTH = 2; // CSS pixels, at the home distance

/** Returns a view that draws into the canvas of the given WebGL2 context, on request. */
export function createPlateView3D(gl) {
  const surface = compile(gl, SURFACE_VS, SURFACE_FS);
  const grain = compile(gl, GRAIN_VS, GRAIN_FS);
  const line = compile(gl, LINE_VS, LINE_FS);
  const lineVertexArray = gl.createVertexArray();
  const lineBuffer = gl.createBuffer();
  let G = 0, heightsTexture = null, marksTexture = null;

  gl.bindVertexArray(lineVertexArray);
  gl.bindBuffer(gl.ARRAY_BUFFER, lineBuffer);
  [0, 1].forEach(end => {
    gl.enableVertexAttribArray(end);
    gl.vertexAttribPointer(end, FLOATS_PER_LINE_VERTEX, gl.FLOAT, false, 2 * BYTES_PER_LINE_VERTEX, end * BYTES_PER_LINE_VERTEX);
    gl.vertexAttribDivisor(end, 1);
  });
  gl.bindBuffer(gl.ARRAY_BUFFER, null);
  gl.bindVertexArray(null);

  const bind = (program, unit, texture, name) => {
    gl.activeTexture(gl.TEXTURE0 + unit);
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.uniform1i(program.uniforms.get(name), unit);
  };

  function drawSurface(matrix, { image, circle, hx, hy, background }) {
    const u = surface.uniforms;
    gl.useProgram(surface.program);
    bind(surface, 0, image, 'uImage');
    bind(surface, 1, marksTexture, 'uMarks');
    gl.uniformMatrix4fv(u.get('uViewProjection'), false, matrix);
    gl.uniform1i(u.get('uCircle'), circle ? 1 : 0);
    gl.uniform2f(u.get('uHalf'), hx, hy);
    gl.uniform3f(u.get('uBackground'), ...background);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
  }

  function drawGrains(matrix, camera, size, { positions, side, count, counts, heapScale, diameter, colour }) {
    const u = grain.uniforms;
    const { right, up, forward } = axesOf(camera);
    gl.useProgram(grain.program);
    bind(grain, 0, positions, 'uPositions');
    bind(grain, 1, counts, 'uCounts');
    gl.uniformMatrix4fv(u.get('uViewProjection'), false, matrix);
    gl.uniformMatrix3fv(u.get('uToWorld'), false, [...right, ...up, ...forward.map(value => -value)]);
    gl.uniform1i(u.get('uSide'), side);
    gl.uniform1i(u.get('uG'), G);
    gl.uniform1f(u.get('uHeapScale'), heapScale);
    gl.uniform1f(u.get('uDiameter'), diameter);
    gl.uniform1f(u.get('uPixelsAtOne'), size / 2 / Math.tan(camera.fov / 2));
    gl.uniform3f(u.get('uSand'), ...colour);
    gl.drawArrays(gl.POINTS, 0, count);
  }

  function drawLines(matrix, size, { vertices, heights, scale }) {
    const u = line.uniforms;
    gl.bindVertexArray(lineVertexArray);
    gl.useProgram(line.program);
    bind(line, 0, heightsTexture, 'uHeights');
    gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, G, G, gl.RED, gl.FLOAT, heights);
    gl.uniformMatrix4fv(u.get('uViewProjection'), false, matrix);
    gl.uniform1f(u.get('uScale'), scale);
    gl.uniform1f(u.get('uWidthAtOne'), LINE_WIDTH * HOME_DISTANCE * size / gl.canvas.clientWidth);
    gl.uniform1f(u.get('uSize'), size);
    gl.uniform1i(u.get('uG'), G);
    gl.bindBuffer(gl.ARRAY_BUFFER, lineBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, vertices, gl.DYNAMIC_DRAW);
    gl.bindBuffer(gl.ARRAY_BUFFER, null);
    gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, vertices.length / (2 * FLOATS_PER_LINE_VERTEX));
  }

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
     * Draw one frame. image: the plate's colours from the sand engine. background: [r, g, b] in 0–1.
     * grains: { positions, side, count: from the sand engine, counts: its grains per cell, heapScale: plate units per grain, diameter: plate units, colour: [r, g, b] in 0–1 }.
     * lines: { vertices: from traceFieldLines, heights: from plateHeight, scale: plate units per metre }.
     */
    draw({ camera, grid, image, grains, lines, circle, hx, hy, background, size }) {
      const sandVertexArray = gl.getParameter(gl.VERTEX_ARRAY_BINDING);
      if (grid !== G) {
        G = grid;
        if (heightsTexture) gl.deleteTexture(heightsTexture);
        heightsTexture = makeTexture(gl, gl.R32F, gl.RED, gl.FLOAT, G, G);
      }
      const matrix = viewProjection(camera);

      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      gl.viewport(0, 0, size, size);
      gl.disable(gl.BLEND);
      gl.clearColor(...background, 1);
      gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

      drawSurface(matrix, { image, circle, hx, hy, background });

      // the plate writes no depth, so grains only hide each other: a grain's height is too small to trust against the plate's
      gl.enable(gl.DEPTH_TEST);
      drawGrains(matrix, camera, size, grains);
      gl.disable(gl.DEPTH_TEST);

      // seen through the plate, since a line dips below it wherever the plate does
      if (lines.vertices.length) {
        gl.enable(gl.BLEND);
        gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
        drawLines(matrix, size, lines);
        gl.disable(gl.BLEND);
      }
      gl.bindVertexArray(sandVertexArray);
    },
  };
}
