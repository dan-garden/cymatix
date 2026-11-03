# Validation

The full evidence behind the summary in the [README](../README.md#how-accurate-is-it). Every paper is linked in the README's [references](../README.md#references).

- [Test files](#test-files)
- [Theory checks](#theory-checks)
- [Square plates](#square-plates)
- [Cantilever plates](#cantilever-plates)
- [Circular plates](#circular-plates)
- [Recovering an unknown rig from a photo](#recovering-an-unknown-rig-from-a-photo)
- [Leissa's book, table by table](#leissas-book-table-by-table)
- [Errors found in Leissa's book](#errors-found-in-leissas-book)
- [Book data not used](#book-data-not-used)
- [Physics added later](#physics-added-later)

## Test files

Run everything with `npm test`. Results on this page are in the tests unless it says otherwise.

| File | Tests | Covers |
|---|---|---|
| [physics.test.js](../physics.test.js) | 13 | Core maths and real experiments |
| [leissa.test.js](../leissa.test.js) | 6 | Leissa chapter 2, Tables 4.3–4.20 |
| [leissa-ch4b.test.js](../leissa-ch4b.test.js) | 6 | Leissa Tables 4.22–4.76 |
| [mindlin.test.js](../mindlin.test.js) | 13 | True Mindlin modes against exact theory and 3D elasticity |
| [fluid.test.js](../fluid.test.js) | 17 | Fluid added mass, radiation, damping and sand mass |
| [patch.test.js](../patch.test.js) | 12 | Finite clamps against exact annular plates and experiments |
| [nonlinear.test.js](../nonlinear.test.js) | 10 | Large-amplitude stiffening |
| [grains.test.js](../grains.test.js) | 14 | Sand and powder transport |
| [plateModel.test.js](../plateModel.test.js) | 13 | All of it working together |
| [plateSurface.test.js](../plateSurface.test.js) | 13 | The plate's height and each tone's shaking, against exact theory |
| [fieldLines.test.js](../fieldLines.test.js) | 7 | Field lines against the drift law and the real grain step |
| [orbitCamera.test.js](../orbitCamera.test.js) | 9 | The 3D camera, picking and zooming |
| [calibrate.test.js](../calibrate.test.js) | 4 | Fitting a plate to measurements |

Where a printed value disagrees, an independent exact solution (scipy) or a measurement decides. In every such case it sides with the app.

## Theory checks

| What | Checked against | Result |
|---|---|---|
| Bessel functions | Their power series and derivative identities | Pass |
| Circular-plate eigenvalues | Closed-form frequency equations, solved separately in scipy | Pass |
| Square-plate eigenvalues | Exact (simply supported) and published (clamped, free) values | Pass |
| Static centre deflection | Plate theory: clamped circle `Fa²/16πD`, simply supported square `0.0116 Fa²/D` ([Timoshenko & Woinowsky-Krieger](https://archive.org/details/theoryofplatessh0000timo_r3s0)) | Within 0.5 % and 1.5 % |
| Mindlin correction, simply supported square | [Hashemi & Arsanjani 2005](https://doi.org/10.1016/j.ijsolstr.2004.06.063), Table 9 | Within 2×10⁻⁵ |
| Cantilever eigenvalues | Leissa Tables 4.51–4.52 | Converged, at most 0.5 % below his 1969 upper bounds |
| Peak finder | A real Blackman-windowed analyser spectrum | Within 0.05 Hz at 5.4 Hz bins, level within 0.05 dB |

Three modelling results also matter:

- **Thin plates still use a per-mode correction,** below h/a = 0.005 on circles and h/L = 0.01 on rectangles. It is already 0.3–0.5 % off at 0.01 on free edges, about the size of a real plate's own uncertainty. Thicker plates use true Mindlin modes: see [below](#thick-plates).
- **Bolted plates:** without the high-mode term, bolted resonances drift up to 3 % high at the app's usual cut-off. With it, they match exact point-support and point-mass solutions to 0.03 %.
- **Long plates:** the basis has the same density along both sides, so long plates get more functions along their length (13 × 102 at b/a = 8).

## Square plates

### Shu et al. 2022

[Shu et al.](https://doi.org/10.3390/e24020215), Fig. 2.

- **Plate:** 320 × 320 × 1 mm aluminium, screwed to a shaker at its centre.
- **Material, from the paper:** E = 70 GPa, ρ = 2700 kg/m³, ν = 0.33.
- **Measured:** 25 sand-pattern resonances, 263–4817 Hz.

| Model | Fitted | Resonances, 200–5000 Hz | f1, f2, f3 error | f4–f25 median | RMS, all 25 |
|---|---|---|---|---|---|
| Rigid bolt | Nothing | 24 (measured: 25) | −20, −5, −7 % | 0.4 % | 4.6 % |
| Shaker, moving mass 20 g | Moving mass (16–25 g all fit) | 24 | −2.4, 0.0, +0.5 % | 0.4 % | 1.7 % |

- **The match is not chance.** Shifting the rigid model by ±3–10 % raises its f4–f25 median to 1.7–3.2 %.
- **Missing resonance.** One close pair near 2.2 kHz merges in the model.
- **Patterns.** Compared by eye with the photos at f1–f7, they match in topology: rings, diamonds, corner loops and the X at 697 Hz.

### Other rigs, scale-free

Not in the tests. One overall frequency factor is fitted, which absorbs unstated thickness or stiffness. "By chance" is the RMS a wrong model reaches, as a median over factors 0.7–1.3.

| Experiment | Mount model | Plate | Best factor | RMS | By chance |
|---|---|---|---|---|---|
| [Shu et al. 2022](https://doi.org/10.3390/e24020215) | Shaker, 20 g | 320 × 1 mm Al, centre screw | 1.007 | 1.7 % | 5.5 % |
| [D'Alessio 2021](https://doi.org/10.1007/s42452-020-04062-6) | Shaker, 20 g | 240 × 0.9 mm Al, PASCO driver rod | 0.88 | 1.8 % | 13 % |
| [Kovacic & Kanovic 2023](https://doi.org/10.3390/sym15091748) | Rigid bolt | 250 × 2 mm, LDS V408 shaker, material not stated | 0.82 | 1.7 % | 16 % |
| [Lin et al. 2024](https://doi.org/10.3390/sym16111460) | Off-centre drive | 320 × 1 mm Al | None | 3.6 % | 3.3 % |

- The order of resonances matches wherever the comparison can discriminate.
- The absolute level is only pinned down where the plate is fully specified, which means Shu et al.
- The D'Alessio and Kovacic rows use only 4–5 frequencies, so they are supporting evidence only.
- Lin's off-centre drive excites so many modes that nearest-peak matching is at chance level. It neither confirms nor refutes the model without matching each photographed pattern.

## Cantilever plates

From [Leissa 1969](https://ntrs.nasa.gov/citations/19700009156), pp. 82–86. These are measured frequency parameters `ωa²√(ρh/D)` of real cantilevers, with one clamped edge.

| Table | Plates | Modes | Measured vs model | In tests |
|---|---|---|---|---|
| 4.55 | Aluminium, a/b = ½, 1, 2, 5 | 14 | Mean −1.3 %, worst −4.2 % | Yes |
| 4.57 | 6061-T6, 7.5 × 7.5 × 0.125 in ([Craig, Plass & Caughfield](https://archive.org/details/DTIC_ADA307623)) | 4 | Mean −3.0 %, worst −4.1 % | Yes |
| 4.58 | Steel, 10 × 10 × 0.0627 in ([Gustafson, Stokey & Zorowski](https://doi.org/10.2514/8.3148)) | 5 | Mean −1.9 %, worst −3.6 % | Yes |
| 4.56 | Steel strips, a/b = 2–5, h/a = 0.02 ([Plunkett](https://doi.org/10.1243/jmes_jour_1963_005_020_02), [Wilson & Plunkett](https://archive.org/details/DTIC_ADA307623)), up to λ ≈ 290 | 45 | Median 1.0 %, mean −0.7 %, worst −5.4 % | Yes |
| 4.53 | Mild steel, 5.12 × 2.76 × 0.053 in, in Hz with handbook steel values | 22 | Mean −3.1 %, worst −5.8 % | No |

- All 68 nondimensional values (4.55–4.58) lie within −5.4 % and +2.6 %.
- Measurements run a little below the model, as expected. A real clamp is never perfectly rigid, which lowers measured frequencies, most of all the fundamental.
- For Table 4.56, the Mindlin correction halves the median error at a/b = 2, from 1.0 % to 0.5 % (chance level: 5.2 %).
- For Table 4.53, Leissa's own printed 1969 theory misses the same measurements by 0.2–8.9 %.

Also compared, not in the tests:

- **Table 4.47:** the square cantilever runs 4–6 % below the model.
- **Fig. 4.46, [Heiba](https://archive.org/details/DTIC_ADA307623):** ⅛ in steel cantilevers, 22 labelled modes, a/b = 1 to 0.4, handbook steel. They average 1.8 % below the model with 4 % scatter.
- **Fig. 4.45:** by eye, the app reproduces the four sand photographs of a square cantilever in mode order: a centre line (torsion), a line near the free edge, two lines along the sides, and a cross.

## Circular plates

- **[Strehlke 1855](https://doi.org/10.1002/andp.18551710805), via [Rayleigh](https://archive.org/details/theorysound07raylgoog) §220:**
  - Kirchhoff's calculated radii are reproduced within 2×10⁻⁴ r/a.
  - Strehlke's measured radii are reproduced within 0.003 r/a.
- **[Waller 1938](https://doi.org/10.1088/0959-5309/50/1/306), nodal circles, Leissa Table 2.9:**
  - 60 measured circles on free brass plates, 0–5 nodal diameters and 1–4 circles.
  - Median error 0.0014 r/a, worst 0.006.
- **[Waller 1938](https://doi.org/10.1088/0959-5309/50/1/306), frequencies, Leissa Table 2.8:**
  - 76 measured frequency ratios, up to 24 nodal diameters and 6 circles.
  - Median error 1.6 %.
- **[Escaler & De La Torre 2018](https://doi.org/10.1016/j.jfluidstructs.2018.07.017), centre-driven plate:**
  - 240 × 0.8 mm aluminium, driven through a bar, moving mass 0.08 kg as stated by the authors. Nothing is fitted.
  - Free-plate frequencies agree with the authors' own theory within 0.1 %.
  - All 6 measured resonances (70–4129 Hz) are reproduced, with exactly 6 in band: +1.8, +1.1, +0.5, +0.3, +0.1 and −0.4 %.
  - Measured nodal circles beyond r/a = 0.2 agree within 0.022 (mean 0.008).
  - The innermost circles lie on the 8 mm drive bar, which the point mount doesn't resolve.

## Recovering an unknown rig from a photo

This is not in the tests. A widely shared photo shows 12 square-plate patterns (345–5907 Hz) but gives no plate details.

1. **What the photo implies.** The centre dot and fourfold symmetry mean a free square bolted at its centre.
2. **What can be fitted.** Frequencies then depend only on `h·√(E/ρ)/L²` and the shaker's moving mass relative to the plate.
3. **The fit.** Those two values give 1.3–1.5 % RMS against the 12 frequencies (chance: 6.7 %), for a moving mass of 5–8 % of the plate.
4. **An example rig.** A 24 cm × 0.75 mm aluminium square on a shaker with about 7 g of moving mass. This is the "12-pattern photo plate" preset.
5. **The patterns.** At the photo's frequencies, 345, 2041, 3835, 4840 and 5201 Hz match closely, and most of the rest match in structure. The photo's 4444 Hz pattern is itself lopsided, which suggests an imperfection in the real plate.
6. **The alternative.** A rigid bolt fits the frequencies at 2.2 % but gives the wrong patterns, so the photos decide between them.

## Leissa's book, table by table

[Leissa, *Vibration of Plates*, NASA SP-160 (1969)](https://ntrs.nasa.gov/citations/19700009156), read from the [scanned pages](https://archive.org/details/DTIC_ADA307623). The tests check every table in chapters 2 and 4 that covers what the app models. Chapter 12's only isotropic table, Tomar's thick circular plate, is one Leissa himself calls "highly questionable", so it is not used.

Edges are listed as SS (simply supported), C (clamped) and F (free).

| Tables | Plate | Values | Agreement |
|---|---|---|---|
| 2.1, 2.2 | Clamped circle, up to 14 diameters and 9 circles | 64 + 88 | Median 0.001 %, radii median 0.0003 r/a |
| 2.3, 2.4 | Simply supported circle | 12 + 18 | 0.01 %. The table is labelled ν = 0.3 but holds the ν = 1/3 roots |
| 2.5, 2.6 | Free circle, ν = 0.33 and 0.25 | 74 + 38 | Median 0.13 % and 0.01 % |
| 2.7 | Free-circle nodal radii | 88 | 67 within 0.005 r/a. In the other 15 measured cases, Waller's measurement sides with the app |
| 2.12–2.15 | Circles held at the centre (free, SS, C), clamped circle with a central mass | 20 | Exact solutions matched to 0.03 % through the app's bolt and shaker models |
| 4.3–4.8 | SS-C-SS-C rectangles, many aspect ratios | 176 | Median 0.03 % |
| 4.10–4.20 | SS-C-SS-SS, SS-C-SS-F, SS-SS-SS-F, SS-F-SS-F rectangles | 275 | Median 0.01–0.2 %, exact Lévy solutions matched to 0.01 % |
| 4.22–4.31 | C-C-C-C: [Claassen & Thorne](https://doi.org/10.1115/1.3641676)'s grid (a/b = 1 to 0.02), [Young](https://doi.org/10.1115/1.4010175), [Iguchi](https://archive.org/details/DTIC_ADA307623), [Ödman](https://archive.org/details/DTIC_ADA307623), mode veering | 497 + 107 | 490 of 497 within 0.01 %. The other 7 are 5 digit slips plus 0.11–0.15 % |
| 4.25, 4.26, 4.32 | C-C-C-C rigorous bounds ([Bazley, Fox & Stadter](https://doi.org/10.1002/zamm.19670470307), [Aronszajn](https://archive.org/details/DTIC_ADA307623)), b/a = 1 to 8 | 110 | All inside |
| 4.33–4.44 | C-C-C-SS, C-C-SS-SS, C-C-F-F, C-F-C-F, C-F-SS-F | About 240 | To the printed figures. Beam-function Ritz upper bounds sit above ours |
| 4.45–4.50 | Cantilevers: [rigorous bounds](https://archive.org/details/DTIC_ADA307623) for a/b = 0.125 to 8, [Barton](https://doi.org/10.1115/1.4010265) and Young upper bounds, [Gontkevich](https://archive.org/details/DTIC_ADA307623) lower bounds | 70 + 20 | All inside |
| 4.61–4.76 | Free rectangles: [rigorous bounds](https://doi.org/10.1007/bf01601716) (b/a = 1 to 8), [Lemke](https://doi.org/10.1002/andp.19283911304) at four ν, [Iguchi](https://doi.org/10.1007/bf00535853) | 80 + 24 + 19 | All bounds hold. Lemke's upper bounds are 0.1–1.2 % above ours |
| 4.68 | [Waller](https://doi.org/10.1088/0959-5309/51/5/312)'s measured free square brass plate, ratios up to 147 | 89 | Median 0.9 %, worst 4.4 % |

## Errors found in Leissa's book

The comparison exposed these misprints and errors. Exact solutions, rigorous bounds, measurements or the trend within the table settle each one.

### Chapter 2

| Table | Printed | Should be |
|---|---|---|
| 2.5, (1,7) | 651.8 | About 553 |
| 2.7 | .925, .926, .993 | .898–.915, from both theory and Waller's measurement |
| 2.15 | "Precise" central-mass roots | Up to 0.34 % higher |

### Chapter 4, first half

| Table | Problem |
|---|---|
| 4.6 | λ* = 24.047 at a/b = 1.5 should be 25.047 |
| 4.7 | The a/b = 1 values contradict Table 4.3 |
| 4.9 | An older approximate source, up to 3.4 % off Table 4.6. Not used |
| 4.17 | Wrong at a/b = 0.6, 0.8 and 5, by up to 7 % |
| 4.19, 4.20 | SS-F-SS-F (m,1) values are wrong, as Leissa himself notes |

### Chapter 4, second half

| Table | Printed | Should be |
|---|---|---|
| 4.30 | 32.9981 | 32.7980 |
| 4.30 | 42.2530 | 42.3253 |
| 4.30 | 75.5403 | 75.4503 |
| 4.30 | 16.2699 | 22.3821 |
| 4.30 | 96.0908 | 96.2324 |
| 4.40 | λ* 18.93 | 17.93 |
| 4.42 | λ* 2.3 | 2.4 |
| 4.42, 4.43 | 61.2 | Too low: it sits below both neighbours |
| 4.67 | [Sigillito](https://archive.org/details/DTIC_ADA307623)'s "upper bound" 13.464 | The true value is 13.468, so this is not an upper bound |

Larger disagreements with named sources:

- **[Kanazawa & Kawai](https://archive.org/details/DTIC_ADA307623), Tables 4.28 and 4.35–4.37:** non-square values run a uniform 2.1–2.3 % high against Claassen & Thorne and the app.
- **[Dill & Pister](https://archive.org/details/DTIC_ADA307623) and Kanazawa & Kawai, C-C-C-SS at a/b = 0.5:** their value of 24.48 lies above a converged Ritz value of 24.14275, unchanged from 20 to 44 functions. A Ritz value bounds the exact one from above.
- **[Iguchi](https://doi.org/10.1007/bf00535853), Table 4.65:** two roots cannot exist (65.368 and 173.695, the second in a gap that Table 4.75's rigorous bounds prove empty). Two more values are above Table 4.67's rigorous upper bounds.

## Book data not used

| Table | Source | Why |
|---|---|---|
| 4.24 | Bolotin | Approximations |
| 4.47 | Theory row | A different plate theory |
| 4.59 | Walton | No plate size given |
| 4.62, 4.63 | Ritz single-term rows | Median 0.5 % off |
| 4.64, 4.72 | Ödman | Single-product estimates |
| 4.69–4.71 | Pavlik | Units unresolved, though ratios are within 1–3 % |
| 4.77 | Waller's rectangles | Medians 1.2–3.4 % with unidentified outliers |

## Physics added later

Each of these answers a limit the README used to list. "Before" is the model as it stood.

### Thick plates

True Mindlin modes replace the per-mode correction for thicker plates. Worst error over the first 12–15 modes:

| Reference | Plate | Before | After |
|---|---|---|---|
| Own Lévy exact solver, independent code | Five edge combinations, h/a = 0.1 | Up to 5.98 % | 5×10⁻⁶ % |
| Own Lévy exact solver, independent code | Same, h/a = 0.2 | Up to 13.1 % | 5×10⁻⁶ % |
| 3D elasticity (own Ritz solver) | Free square, h/a = 0.1, 0.2 | 3.3 %, 5.2 % | 0.60 %, 1.33 % |
| 3D elasticity | Clamped square, h/a = 0.1, 0.2 | 5.2 %, 10.5 % | 1.68 %, 2.86 % |
| 3D elasticity | Free disc, h/a = 0.1 | 0.92 % | 0.51 % |
| 3D elasticity | Clamped disc, h/a = 0.1 | 2.12 % | 0.99 % |
| Own polar Ritz solver, no Bessel functions | All edges, h/a = 0.1, 0.2 | Not modelled | 9×10⁻⁴ %, same root counts |

- The clamped 3D reference holds the whole face, which is stiffer than Mindlin's clamp, so part of those errors may be the boundary condition.
- Published values for n = 0 from Senjanović et al., *Arch. Mech.* 66, 389 (2014), [open copy via the journal](https://am.ippt.gov.pl/), agree within 0.002 in Λ. Their n ≥ 1 values use a different theory and were not used.
- Effect on thin experiments is small: Escaler moves by at most 0.003 % on axisymmetric modes.

**Waller's outer circles.** Leissa Table 2.8 gives her plate thickness: h/a = 0.0086–0.0156. Exact Mindlin and 3D move the outer circle by under 0.0002 r/a there, against the 0.012–0.021 gap. Closing it would need h/a of about 0.2, where frequencies fall by half or more. Poisson's ratio moves it by at most 0.0003. Her s = 5 and 6 outer radii (.956, .958) lie within 0.0034 and 0.0021 of Leissa's eq. 2.19, a formula that drops the I_n term, while her s ≤ 4 data follow exact theory. Leissa's printed Table 2.7 s = 5 outer radii are eq. 2.19 for all six n. The chance of both landing that close is about 3 %, a rough estimate. Whether her two entries were measured needs her paper, which could not be read.

### Fluid loading and damping

| Reference | Quantity | Before | After |
|---|---|---|---|
| Rigid disc, closed form | Added mass and rocking | Not modelled | 10⁻¹² |
| Kwak 1991, as printed in Escaler Table 3 | Free disc, six modes | Not modelled | 1.5×10⁻⁴ |
| Own vortex-ring and vortex-lattice solvers | Circles and rectangles | Not modelled | 3×10⁻⁴, 1 % on cantilevers |
| Escaler Table 5, their FEM | Free plate in water | Standard factor 1.7–11 % off | Within 3 % |
| Escaler Table 7, measured | Five modes in water | Standard factor 7.0 % RMS | 3.4 % RMS, all low |
| Lambourg et al. Fig. 10a, read from the figure | 21 decay rates | Best single ζ: factor 1.4 | Factor 1.17 |

What did not improve:

- **Shu et al., 25 resonances.** 20 g: 1.74 % in vacuum, 1.92 % with air. Moving mass refitted: 1.86 % at 15 g. With a free scale as well: 1.66 % in vacuum, 1.62 % with air.
- **Escaler in air, nothing fitted** (λ ≤ 25 / λ ≤ 45). Vacuum: 0.93 / 0.91 % RMS. Air: 0.61 / 0.91 %.
- **Lindholm's cantilevers in water,** read by eye from a secondary figure: the model runs 10–30 % low. Unresolved.
- **A second measured plate** ([Ege et al.](https://arxiv.org/abs/1210.3333), 2 mm aluminium): 0.22 % measured against 0.038 % modelled, so physical damping is an option, not the default. A rig also loses energy at the bolt, the shaker and the sand, which are not modelled.

### Finite clamps

| Reference | Quantity | Result |
|---|---|---|
| Own exact annular solver, b/a 0.02–0.1, three outer edges | λ² | Within 0.35 % |
| Same, b/a = 0.3 | λ² | Within 1 % |
| Leissa Tables 2.24, 2.30, 2.31 | λ², b/a 0.1–0.4 | Within 1.3 % |
| Own exact hub with 0.08 kg, Escaler's plate | Frequencies, nodal radii | Within 0.3 %, 0.003 r/a |
| Own exact hub with a suspended shaker | Four resonances | Within 0.3 % |
| r₀ → 0 against the point model | Coefficients | 10⁻¹⁵ on rectangles, 10⁻⁷ on a centred circle |

Compared with Escaler's measurements, nothing fitted:

| Quantity | Point | 4 mm clamp |
|---|---|---|
| Five innermost circles, RMS r/a | 0.033 | 0.0021 |
| Sixteen outer circles, RMS r/a | 0.0093 | 0.0071 |
| Six frequencies | +1.8 … −0.3 % | +3.8 … +6.0 % |
| The same, with air | +1.0 … −0.7 % | +3.1 … +5.3 % |

- The 4 mm radius assumes the paper's 8 mm bar is a diameter.
- Removing the mean offset leaves 0.7–0.8 % RMS either way, so frequencies do not choose between them. A plate about 5 % thinner or softer would cancel the offset. That is unverified.
- On Shu's plate no washer size helps: RMS rises from 1.45 % to 1.85 % at 4 mm, and a rigid 10 mm clamp still has f1 at −14.7 %.
- Exact theory shows a clamp shifts n = 0 modes by +1.1 to +2.7 % at r₀/a = 0.02 and +13 to +22 % at 0.1. Modes with n ≥ 2 move by at most 1.2 % up to r₀/a = 0.05.
- Some Leissa annular cells disagree with exact theory, for example free n = 1 printed 3.14 (exact 3.478), which may be OCR error, as the scan was not checked.

### Large-amplitude stiffening

| Reference | Quantity | Error |
|---|---|---|
| Yamaki 1961 via Leissa Table 12.4 | Circles, four cases | Within 0.5 % |
| Yamaki via Leissa Sec. 12.2, closed form | Simply supported square | Within 0.1 % |
| VK-Gong manual (Touzé group), ν = 0.38 | Three circular modes | Within 0.01 %, with one at −0.7 %, a likely slip in the manual |
| Own exact integrals, own displacement Ritz | Circles, rectangles | 0.05 %, 0.2 % |
| Own time integration of a forced Duffing mode | Amplitude and hysteresis, both sweep directions | Within 1.5 % |

The arbitrary warning at a deflection of 20 % of the thickness is replaced by `s = ⅜·κ_A·(a/ω²h)²/ζ`:

| s | Meaning |
|---|---|
| Below 0.1 | Negligible |
| 0.1–1.5 | Corrected, within 1 % at ζ ≤ 0.005 |
| 1.5–3 | Caution |
| Above 3 | The one-mode model has broken down |

- Nothing here has been checked against a measurement. The free-edge experiments of Thomas et al. could not be read.
- Bolted plates need coefficients for the pinned plate, which are not wired in, so they warn and do not correct.
- Stiffening is under 0.1 half-widths above about 500 Hz on 1 mm plates and about 1 kHz on 0.3 mm plates, even at 100 g.

### Sand and powder

| Reference | Quantity | Result |
|---|---|---|
| Zhou et al. 2016 Fig. 3b, digitised | Voltage for a 175 µm displacement, 55 notes | Factor 1.5–1.8 on average, 1.9–3 note to note. Shuffled notes: 5.7 |
| Zhou Fig. 3a, 3990 Hz | Displacement against note length | .010/.020/.040/.079/.159 mm, measured .010/.020/.044/.105/.151 |
| The same note with tilt only | Displacement | ×150 too small |
| Own polygon-grain simulation | Drift across Γ, e, μ, h, f, D | Within 20 % |
| Own boundary-layer solvers | Streaming coefficient | Rayleigh's ratio to 1.0000 |
| Faraday 1831 | Powder to antinodes, sand to nodes | Right direction, qualitative only |
| Loubet et al. 2007 | Lycopodium settling speed | 4.07 cm/s, measured 4.2 |

What the mechanism turned out to be:

- Friction against the top surface's in-plane motion drives sand. Tilt is undetectable: at Γ = 5, 1 kHz and h = 1 mm, friction gives −1.37 mm/s and tilt gives 0.001 ± 0.067 mm/s.
- Drift does not depend on grain size and is linear in thickness and in the amplitude gradient.
- Bounce scatter from irregular shape dominates each hop and sets the random spread.

What is assumed or unchecked:

- Formation time for real sand has no measured reference in anything read.
- Restitution (0.7) and friction (0.35) are from sand-grain simulations, not sand on metal.
- The random spread comes from a 2D simulation, varies by ×3 between runs and ×10 across friction, and is capped by air drag using a derived estimate.
- [Dorrestijn et al. 2007](https://doi.org/10.1103/physrevlett.98.026102) studied cantilevers in water, and rolling below Γ = 1 is [van Gerner et al. 2010](https://doi.org/10.1103/physreve.82.012301). Sand is irregular and does not roll, so that regime is not modelled.
- Whether powder goes to antinodes depends on drive, and there is no single crossover size. The direction agrees with the simulation in 67 of 75 runs, with mean error ×1.6 and ×3.7 on the runs that go to antinodes.
- With the air-drag cap on spread, a Faraday-style run (35 µm, 308 Hz, 100 g) ends with 37 % at the centre, 45 % at the rim and 17 % at the node. This is the lead's own check, against Faraday's description only.

### GPU and CPU

The same fixed field and the same 60× run on both engines give the same pattern by eye. The GPU engine's statistical test on a known field gives a settled width of 0.016 against about 0.014 from theory, with no grain leaving the plate.

### Numerical bug found on the way

The integral for I_n sums terms about e^(n²/2x) times larger than the answer. Raising the circle limit to λ = 120 exposed it as false modes: 7851 found, about 1970 expected. A power series now takes over, and 1873 modes result. A sample of 74 roots agrees with scipy to 9×10⁻¹³. The old limit of λ = 50 stayed below where it started, so earlier results were safe.

### 3D view

The view draws existing results, so it adds no physics of its own and cannot make the model more accurate.

| Check | Result | In tests |
|---|---|---|
| Slowly pushed clamped disc against `Fa²/16πD` (Timoshenko), in phase | Within 0.5 %, and 10⁻⁴ of that a quarter cycle later | Yes |
| Plate at resonance | Lags the force by a quarter cycle, with under 5 % in phase | Yes |
| Each tone on a simply supported plate against the exact Navier series, in size and phase, at 220, 660 and 1760 Hz | Worst 0.42 %, 0.26 % and 0.70 % of the tone's peak | Yes |
| Height after a full cycle, and against the amplitude | Same within 10⁻⁵, and never above it | Yes |
| Several tones together, each at its own phase | The sum of each alone within 10⁻⁵, and never past the peak deflection | Yes |
| Tones chosen by shaking, from a real saw wave's heights | Picks 220 and 660 Hz, with 660 Hz shaking hardest at a fifth of the height | Yes |
| Clamped rim, with air and true Mindlin modes | Under 2 % of peak one grid cell from the edge | Yes |
| Heap volume against the grains' volume | Conserved to 10⁻¹² | Yes |
| Sand lines on an exact (2,2) mode against the drift law `u ∝ −∇\|W\|` | Direction within 18° everywhere off the antinodes, amplitude always falling | Yes |
| Where sand lines end on that mode at 5 g | Within one step of the 0.95 g threshold | Yes |
| A grain moved by the real sand step, from partway along a line | Stays within one cell of it and stops at its end | Yes |
| Powder lines climbing the amplitude, on the same mode | 0 % at 5 g, 53 % at 100 g, 72 % at 200 g, each ending on an antinode, as [Faraday 1831](https://doi.org/10.1098/rstl.1831.0018) saw | Yes |
| Settled sand within one cell of where the plate stands almost still (nodes of the standing shape, moving under 10 % of the peak), on a free square at 982 Hz | 90.5 %, on places covering 24.5 % of the plate | Yes |
| The same on a free disc at 3072 Hz, and at 1000 Hz between modes | 78.9 % on 30.4 %, and 86.2 % on 26.9 % | No, one-off check |
| Highest field line on a clamped disc at its fundamental, seen side-on at its peak, against true height × magnification through the real camera | Pixel row predicted 298.3, 273.8 and 107.9, drawn 298, 273 and 107 (flat, ×200, ×1,000) | No, browser check |
| 3D from almost straight above at the 2D view's scale, with lines off, against the 2D image, same grains | Correlation 0.979 in 16-pixel blocks and 0.938 in 8-pixel blocks, rank correlation 0.987 and 0.980 | No, browser check |
| Zooming towards the pointer | The plate point under it stays there within 10⁻⁴, from any angle | Yes |
| Frame rate in 3D with field lines on, NVIDIA RTX 3080 | 52, 55 and 56 frames/s at 100,000, 1,000,000 and 3,000,000 grains, so grains were not the limit | No, browser check |
| Clicking the 3D plate at (0.4, 0.3) | Driver at r = 0.50, 37°, target 0.50, 36.9° | No, browser check |

The field lines move about the plate's own plane because the height is measured from the plate. A line therefore barely moves where the plate stands still, which is where the sand sits, from any angle.

Grains are drawn from the same positions the simulation moves, so they cannot disagree with the 2D view on where a grain is. The correlation above stays below 1 because the two views colour them differently: 2D brightens where grains overlap, and 3D draws each one as a lit ball. A grain's height within a pile is a display estimate, since grain contact is not simulated.

### Bugs found by the 3D view

**Plates sharing a solve.** Solved plates were reused by a key that left out the edge type, Poisson's ratio and aspect ratio. In the app, a plate with the same thickness ratio as an earlier one could silently take that plate's modes: the Clamped drum preset ran at 597 Hz, a free plate's frequency, instead of 1433 Hz. The 3D view exposed it, because the clamped rim visibly moved. The key now holds every input to the solve, and a test checks it. Published results were not affected, since the tests always solve fresh, but the earlier gallery and demo showed the wrong clamped plate and have been replaced.

**False motion at the rim.** Cells just outside the plate copy their nearest plate cell, but the slopes were still taken across the end of those copies, and the air-streaming term took them across the rim itself. Field lines showed it as an impossible top speed. Fixes in `transportFields`: streaming now uses the copied rim, and the outside cells take their plate cell's motion. A regression test covers it.

| Case | Fastest inside | Rim before | Rim after | Grains in the outer 5 % of the radius, before and after |
|---|---|---|---|---|
| Lycopodium, free disc, 308 Hz, 100 g | 11.8 mm/s | 1677 mm/s | 13.6 mm/s | 39.4 %, 43.0 % |
| Sand, free disc, 220 Hz, 10 g | 0.42 mm/s | 7.46 mm/s | 0.57 mm/s | 4.9 %, 5.0 % |
| Sand, free disc, 3072 Hz, 12 g | 0.16 mm/s | 0.86 mm/s | 0.17 mm/s | 2.2 %, 2.2 % |

Each step is capped at half a grid cell, so the false speed moved few grains. The published patterns barely change, and the gallery is regenerated anyway.

**Peak deflection with several tones.** The readout divided the total acceleration by the strongest tone's ω², which overstated it: the square-wave preset read 51.3 µm when its tones add up to at most 23.3 µm. It now adds the tones' real amplitudes and says "up to". A single tone is unchanged.

### Calibration

Fitting Shu's 25 resonances, from the stated plate with air, recovers a frequency scale of ×1.011 and a moving mass of 22 g, at 1.6 % RMS against 5.4 % by chance. The 12-pattern photo refits to ×1.009 and 7.0 g, at 1.28 % against 3.96 %. A preset's fit means little unless it beats chance.
