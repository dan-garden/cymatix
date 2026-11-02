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

## Test files

Run everything with `npm test`. Results on this page are in the tests unless it says otherwise.

| File | Tests | Covers |
|---|---|---|
| [physics.test.js](../physics.test.js) | 12 | Core maths and real experiments |
| [leissa.test.js](../leissa.test.js) | 6 | Leissa chapter 2, Tables 4.3–4.20 |
| [leissa-ch4b.test.js](../leissa-ch4b.test.js) | 6 | Leissa Tables 4.22–4.76 |

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

- **Thick plates, other edges:** at h/a = 0.1, thin-plate theory is 6–20 % high and the corrected value is 0.4–5 % off. Below h/L ≈ 0.02 the correction is negligible. This comparison is not in the tests.
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
