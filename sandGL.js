// GPU sand: moves and draws millions of grains with WebGL2, so a simulation can use as many grains as a real rig.
// Grain positions live in a floating-point texture. Each frame the CPU uploads the plate's transport fields once, a
// fragment shader advances every grain from them, and the grains are splatted into a density texture and coloured.

const QUAD_VS = `#version 300 es
void main() {
  vec2 corner = vec2(gl_VertexID == 1 ? 3.0 : -1.0, gl_VertexID == 2 ? 3.0 : -1.0);
  gl_Position = vec4(corner, 0.0, 1.0);
}`;

const FIELD_SAMPLER = `
uniform highp sampler2D uField; // RGBA32F: drift x, drift y, diffusivity, drive level
uniform int uG;
vec4 fieldAt(vec2 p) { // bilinear, clamped to the grid
  float cells = float(uG) / 2.0, end = float(uG) - 1.001;
  vec2 f = clamp((p + 1.0) * cells - 0.5, vec2(0.0), vec2(end));
  ivec2 c = ivec2(f);
  vec2 a = f - vec2(c);
  vec4 bottom = mix(texelFetch(uField, c, 0), texelFetch(uField, c + ivec2(1, 0), 0), a.x);
  vec4 top = mix(texelFetch(uField, c + ivec2(0, 1), 0), texelFetch(uField, c + ivec2(1, 1), 0), a.x);
  return mix(bottom, top, a.y);
}`;

/** GLSL: the bilinear value of a G×G texture at plate coordinates, clamped to the grid, as one or two channels. Needs uniform int uG. */
export const GRID_VALUE = `
vec2 gridPair(sampler2D grid, vec2 p) {
  float cells = float(uG) / 2.0, end = float(uG) - 1.001;
  vec2 f = clamp((p + 1.0) * cells - 0.5, vec2(0.0), vec2(end));
  ivec2 c = ivec2(f);
  vec2 a = f - vec2(c);
  vec2 bottom = mix(texelFetch(grid, c, 0).rg, texelFetch(grid, c + ivec2(1, 0), 0).rg, a.x);
  vec2 top = mix(texelFetch(grid, c + ivec2(0, 1), 0).rg, texelFetch(grid, c + ivec2(1, 1), 0).rg, a.x);
  return mix(bottom, top, a.y);
}
float gridValue(sampler2D grid, vec2 p) { return gridPair(grid, p).r; }`;

const STEP_FS = `#version 300 es
precision highp float; precision highp int; precision highp sampler2D;
uniform sampler2D uPos;
uniform float uDt;
uniform uint uFrame;
uniform bool uCircle;
uniform vec2 uHalf;
uniform float uMaxStep;
uniform sampler2D uCounts; // R32F: grains per cell
uniform float uCapacity;
uniform float uCrowdSpread;
out vec4 outPos;
${FIELD_SAMPLER}
uvec2 pcg2d(uvec2 v) { // Jarzynski & Olano 2020
  v = v * 1664525u + 1013904223u;
  v.x += v.y * 1664525u; v.y += v.x * 1664525u;
  v ^= v >> 16u;
  v.x += v.y * 1664525u; v.y += v.x * 1664525u;
  v ^= v >> 16u;
  return v;
}
vec2 gaussian(uvec2 seed) { // Box–Muller
  vec2 u = vec2(pcg2d(seed)) / 4294967296.0;
  float radius = sqrt(-2.0 * log(max(u.x, 1e-9)));
  return radius * vec2(cos(6.28318530718 * u.y), sin(6.28318530718 * u.y));
}
vec2 reflectInside(vec2 p, vec2 limit) {
  p = mix(p, 2.0 * limit - p, vec2(greaterThan(p, limit)));
  p = mix(p, -2.0 * limit - p, vec2(lessThan(p, -limit)));
  return clamp(p, -limit, limit);
}
void main() {
  ivec2 cell = ivec2(gl_FragCoord.xy);
  vec2 p = texelFetch(uPos, cell, 0).xy;
  vec4 f = fieldAt(p);
  // decided per grain, so a quiet zone has a smooth edge rather than the grid's staircase (activityAt in grains.js)
  float activity = clamp((f.a - 0.95) / 0.7, 0.0, 1.0);
  float cells = float(uG) / 2.0, end = float(uG) - 1.001;
  vec2 grid = clamp((p + 1.0) * cells - 0.5, vec2(0.0), vec2(end));
  ivec2 c = ivec2(grid);
  vec2 a = grid - vec2(c);
  // how overfull the cell is, as updateCrowding in grains.js
  vec4 corners = vec4(texelFetch(uCounts, c, 0).r, texelFetch(uCounts, c + ivec2(1, 0), 0).r,
    texelFetch(uCounts, c + ivec2(0, 1), 0).r, texelFetch(uCounts, c + ivec2(1, 1), 0).r);
  corners = clamp(corners / uCapacity - 1.0, 0.0, 4.0);
  float crowded = mix(mix(corners.x, corners.y, a.x), mix(corners.z, corners.w, a.x), a.y);
  if (activity > 0.0 || crowded > 0.0) {
    f.xyz *= activity;
    f.z += crowded * uCrowdSpread;
    vec2 noise = gaussian(uvec2(cell) + uFrame * uvec2(4099u, 7919u));
    p += clamp(f.xy * uDt, vec2(-uMaxStep), vec2(uMaxStep)) + sqrt(2.0 * max(f.z, 0.0) * uDt) * noise;
    if (uCircle) {
      float r = length(p);
      if (r > 0.995) p *= max(0.0, 1.99 - r) / r;
    } else {
      p = reflectInside(p, uHalf * 0.995);
    }
  }
  outPos = vec4(p, 0.0, 1.0);
}`;

const SPLAT_VS = `#version 300 es
precision highp float; precision highp int; precision highp sampler2D;
uniform sampler2D uPos;
uniform int uSide;
uniform float uSize;
uniform float uGrainSize;
void main() {
  vec2 p = texelFetch(uPos, ivec2(gl_VertexID % uSide, gl_VertexID / uSide), 0).xy;
  vec2 pixel = min(floor((p + 1.0) * uSize / 2.0), vec2(uSize - uGrainSize));
  gl_Position = vec4((pixel + uGrainSize / 2.0) / uSize * 2.0 - 1.0, 0.0, 1.0);
  gl_PointSize = uGrainSize;
}`;

const SPLAT_FS = `#version 300 es
precision highp float;
out vec4 outDensity;
void main() { outDensity = vec4(1.0); }`;

const SHOW_FS = `#version 300 es
precision highp float; precision highp int; precision highp sampler2D;
uniform sampler2D uDensity;
uniform sampler2D uBase;   // RGBA8 plate shading, top row first
uniform sampler2D uSand;   // 32 × 1 colours by grain count
uniform sampler2D uMap;    // 512 × 1 acceleration colours
uniform bool uShowMap;
uniform float uMapTop;
uniform float uSize;
out vec4 outColor;
uniform highp sampler2D uGamma; // R32F, the acceleration map
uniform int uG;
${GRID_VALUE}
void main() {
  // the canvas is bottom-up, the plate's coordinates and base image are top-down
  ivec2 pixel = ivec2(gl_FragCoord.x, uSize - gl_FragCoord.y);
  float density = texelFetch(uDensity, pixel, 0).r;
  if (density >= 0.5) {
    outColor = texelFetch(uSand, ivec2(min(31, int(density + 0.5)), 0), 0);
    return;
  }
  vec4 base = texelFetch(uBase, pixel, 0);
  if (uShowMap && base.a > 0.5) {
    float g = gridValue(uGamma, (vec2(pixel) + 0.5) / uSize * 2.0 - 1.0);
    int shade = (g < 1.0 ? 0 : 256) + min(255, int(g / uMapTop * 255.0));
    outColor = texelFetch(uMap, ivec2(shade, 0), 0);
    return;
  }
  outColor = vec4(base.rgb, 1.0);
}`;

export function compile(gl, vertexSource, fragmentSource) {
  const shader = (type, source) => {
    const s = gl.createShader(type);
    gl.shaderSource(s, source);
    gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
    return s;
  };
  const program = gl.createProgram();
  gl.attachShader(program, shader(gl.VERTEX_SHADER, vertexSource));
  gl.attachShader(program, shader(gl.FRAGMENT_SHADER, fragmentSource));
  gl.linkProgram(program);
  if (!gl.getProgramParameter(program, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(program));
  const uniforms = new Map();
  const count = gl.getProgramParameter(program, gl.ACTIVE_UNIFORMS);
  for (let i = 0; i < count; i++) {
    const { name } = gl.getActiveUniform(program, i);
    uniforms.set(name, gl.getUniformLocation(program, name));
  }
  return { program, uniforms };
}

export function makeTexture(gl, internalFormat, format, type, width, height, data = null) {
  const texture = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, texture);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  gl.texImage2D(gl.TEXTURE_2D, 0, internalFormat, width, height, 0, format, type, data);
  return texture;
}

function makeTarget(gl, texture) {
  const framebuffer = gl.createFramebuffer();
  gl.bindFramebuffer(gl.FRAMEBUFFER, framebuffer);
  gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, texture, 0);
  return framebuffer;
}

const bytesOf = pixels => new Uint8Array(pixels.buffer, pixels.byteOffset, pixels.byteLength);

/** Returns a GPU sand engine for the canvas, or null where the browser lacks what it needs (the CPU path then runs). */
export function createSandGL(canvas) {
  // tried on a spare canvas first: a canvas that has given out a WebGL context can never give the CPU path a 2D one
  const probe = document.createElement('canvas').getContext('webgl2');
  if (!probe || !probe.getExtension('EXT_color_buffer_float') || !probe.getExtension('EXT_float_blend')) return null;
  const gl = canvas.getContext('webgl2', { antialias: false, alpha: false, preserveDrawingBuffer: true });
  gl.getExtension('EXT_color_buffer_float');
  gl.getExtension('EXT_float_blend');

  const step = compile(gl, QUAD_VS, STEP_FS);
  const splat = compile(gl, SPLAT_VS, SPLAT_FS);
  const show = compile(gl, QUAD_VS, SHOW_FS);
  gl.bindVertexArray(gl.createVertexArray());

  const textures = { field: null, gamma: null, base: null, sand: null, map: null, density: null, counts: null, image: null };
  let countsTarget = null, imageTarget = null, capacity = 1;
  let G = 0, size = 0, side = 0, frame = 0, densityTarget = null;
  let positions = []; // [{ texture, target }, { texture, target }], swapped every step

  const bind = (unit, texture, location) => {
    gl.activeTexture(gl.TEXTURE0 + unit);
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.uniform1i(location, unit);
  };

  function replace(key, ...args) {
    if (textures[key]) gl.deleteTexture(textures[key]);
    textures[key] = makeTexture(gl, ...args);
  }

  return {
    gl,

    /** The finished image, the grain counts and the grains themselves, for views that build on them. */
    get image() { return textures.image; },
    get counts() { return textures.counts; },
    get positions() { return positions[0].texture; },
    get side() { return side; },

    resize(newSize) {
      size = newSize;
      replace('density', gl.R16F, gl.RED, gl.HALF_FLOAT, size, size);
      densityTarget = makeTarget(gl, textures.density);
      replace('image', gl.RGBA8, gl.RGBA, gl.UNSIGNED_BYTE, size, size);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
      imageTarget = makeTarget(gl, textures.image);
    },

    /** Upload grain positions in plate coordinates [-1, 1]. */
    setGrains(x, y) {
      for (const { texture } of positions) gl.deleteTexture(texture);
      side = Math.ceil(Math.sqrt(x.length));
      const data = new Float32Array(side * side * 4);
      // spare texels stay out of the draw call, so their contents do not matter
      for (let i = 0; i < x.length; i++) {
        data[4 * i] = x[i];
        data[4 * i + 1] = y[i];
      }
      positions = [data, null].map(pixels => {
        const texture = makeTexture(gl, gl.RGBA32F, gl.RGBA, gl.FLOAT, side, side, pixels);
        return { texture, target: makeTarget(gl, texture) };
      });
    },

    /** Upload the G×G fields: RGBA per cell = drift x, drift y (plate units/s), diffusivity (plate units²/s), level. */
    setFields(grid, data, gamma) {
      if (grid !== G) {
        G = grid;
        replace('field', gl.RGBA32F, gl.RGBA, gl.FLOAT, G, G, data);
        replace('gamma', gl.R32F, gl.RED, gl.FLOAT, G, G, gamma);
        replace('counts', gl.R32F, gl.RED, gl.FLOAT, G, G);
        countsTarget = makeTarget(gl, textures.counts);
        return;
      }
      gl.bindTexture(gl.TEXTURE_2D, textures.field);
      gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, G, G, gl.RGBA, gl.FLOAT, data);
      gl.bindTexture(gl.TEXTURE_2D, textures.gamma);
      gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, G, G, gl.RED, gl.FLOAT, gamma);
    },

    /** Colours as the packed RGBA values the CPU path uses. */
    setBase(pixels) {
      replace('base', gl.RGBA8, gl.RGBA, gl.UNSIGNED_BYTE, size, size, bytesOf(pixels));
    },
    setSandColors(pixels) {
      replace('sand', gl.RGBA8, gl.RGBA, gl.UNSIGNED_BYTE, pixels.length, 1, bytesOf(pixels));
    },
    setMapColors(pixels) {
      replace('map', gl.RGBA8, gl.RGBA, gl.UNSIGNED_BYTE, pixels.length, 1, bytesOf(pixels));
    },

    /** Count grains per cell, once per frame before stepping, to find overfull cells. */
    updateCrowding(grainCount, cellCapacity) {
      capacity = cellCapacity;
      gl.bindFramebuffer(gl.FRAMEBUFFER, countsTarget);
      gl.viewport(0, 0, G, G);
      gl.clearColor(0, 0, 0, 0);
      gl.clear(gl.COLOR_BUFFER_BIT);
      gl.enable(gl.BLEND);
      gl.blendFunc(gl.ONE, gl.ONE);
      gl.useProgram(splat.program);
      const u = splat.uniforms;
      bind(0, positions[0].texture, u.get('uPos'));
      gl.uniform1i(u.get('uSide'), side);
      gl.uniform1f(u.get('uSize'), G);
      gl.uniform1f(u.get('uGrainSize'), 1);
      gl.drawArrays(gl.POINTS, 0, grainCount);
    },

    /** Advance every grain by dt seconds through the uploaded fields. */
    step(dt, { circle, hx, hy }) {
      const [from, to] = positions;
      gl.bindFramebuffer(gl.FRAMEBUFFER, to.target);
      gl.viewport(0, 0, side, side);
      gl.disable(gl.BLEND);
      gl.useProgram(step.program);
      const u = step.uniforms;
      bind(0, from.texture, u.get('uPos'));
      bind(1, textures.field, u.get('uField'));
      bind(2, textures.counts, u.get('uCounts'));
      gl.uniform1f(u.get('uCapacity'), capacity);
      gl.uniform1f(u.get('uCrowdSpread'), (2 / G) ** 2); // crowdDiffusivity in grains.js
      gl.uniform1i(u.get('uG'), G);
      gl.uniform1f(u.get('uDt'), dt);
      gl.uniform1f(u.get('uMaxStep'), 1 / G); // half a grid cell, as in stepGrains
      gl.uniform1ui(u.get('uFrame'), ++frame);
      gl.uniform1i(u.get('uCircle'), circle ? 1 : 0);
      gl.uniform2f(u.get('uHalf'), hx, hy);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
      positions = [to, from];
    },

    /** Splat grains into the density texture and colour the canvas, or with offscreen colour only the plate, into the image texture. */
    draw(grainCount, { grainSize, showMap, mapTop, offscreen = false }) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, densityTarget);
      gl.viewport(0, 0, size, size);
      gl.clearColor(0, 0, 0, 0);
      gl.clear(gl.COLOR_BUFFER_BIT);
      if (!offscreen) {
        gl.enable(gl.BLEND);
        gl.blendFunc(gl.ONE, gl.ONE);
        gl.useProgram(splat.program);
        const u = splat.uniforms;
        bind(0, positions[0].texture, u.get('uPos'));
        gl.uniform1i(u.get('uSide'), side);
        gl.uniform1f(u.get('uSize'), size);
        gl.uniform1f(u.get('uGrainSize'), grainSize);
        gl.drawArrays(gl.POINTS, 0, grainCount);
      }

      gl.bindFramebuffer(gl.FRAMEBUFFER, offscreen ? imageTarget : null);
      gl.viewport(0, 0, size, size);
      gl.disable(gl.BLEND);
      gl.useProgram(show.program);
      const u = show.uniforms;
      bind(0, textures.density, u.get('uDensity'));
      bind(1, textures.base, u.get('uBase'));
      bind(2, textures.sand, u.get('uSand'));
      bind(3, textures.map, u.get('uMap'));
      bind(4, textures.gamma, u.get('uGamma'));
      gl.uniform1i(u.get('uG'), G);
      gl.uniform1i(u.get('uShowMap'), showMap ? 1 : 0);
      gl.uniform1f(u.get('uMapTop'), mapTop);
      gl.uniform1f(u.get('uSize'), size);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
      if (offscreen) { // mipmaps keep the rim and the acceleration map steady when the plate is seen at an angle
        gl.bindTexture(gl.TEXTURE_2D, textures.image);
        gl.generateMipmap(gl.TEXTURE_2D);
      }
    },

    /** Current grain positions, for statistics and tests. Slow: reads the texture back. */
    readGrains() {
      const data = new Float32Array(side * side * 4);
      gl.bindFramebuffer(gl.FRAMEBUFFER, positions[0].target);
      gl.readPixels(0, 0, side, side, gl.RGBA, gl.FLOAT, data);
      return data;
    },
  };
}
