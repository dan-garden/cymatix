# Cymatix

Sound drives a simulated vibrating plate, and sand collects on its nodal lines (Chladni figures), in real time.
The source can be an audio/video file, the microphone, or a test tone (sine, square, sawtooth or triangle, any frequency).
Presets load a plate, driver and tone and clear the sand; "Reset all" returns every setting to its default, as does reloading.
Under Plate → Eigenmodes, clicking a mode plays its exact frequency (hover a circular mode for its nodal radii).

```bash
npm start   # serves on http://localhost:8137 (ES modules need http; the mic needs localhost or https)
npm test    # physics checks against theory and experiment
```

No build step and no dependencies: `index.html`, `app.js` (audio, sand, drawing), `physics.js` (all the maths).

## What is computed

- **Spectrum** – Web Audio FFT; the strongest peaks (parabolic interpolation) become sinusoidal drives at the driver,
  scaled so a 0 dBFS sine gives the "drive at 0 dBFS". Auto level instead sets the drive so the antinodes reach the target g.
- **Plate** – Kirchhoff thin plate, D = Eh³/12(1−ν²).
  Circle: exact modes [J_n(kr) + C·I_n(kr)]·cos nθ from the free, clamped or simply supported frequency equations.
  Rectangle: Rayleigh–Ritz on orthonormalised Legendre bases, each edge free, simply supported or clamped independently
  (e.g. a cantilever: one clamped edge, three free); the eigenproblem is solved by Householder tridiagonalisation + QL.
  The basis has the same density along both sides, so long plates get more functions along their length (b/a = 8: 13 × 102).
- **Thick plates** – each mode's frequency is corrected for shear deformation and rotary inertia with Mindlin's dispersion
  relation at the mode's bending wavenumber (κ² = π²/12). Exact for simply supported edges; for clamped/free edges it removes
  most of the thick-plate error (at h/a = 0.1: thin-plate theory is 6–20 % high, corrected 0.4–5 %). Below h/L ≈ 0.02 it is negligible.
- **Response** – steady-state modal superposition with modal damping ζ. G(x,x₀) = Σ φ_k(x)φ_k(x₀) / (m(ω_k² − ω² + 2iζω_kω)).
  - *Point force* F: W = F·G(x,x₀).
  - *Bolted to shaker* (a real Chladni rig): coil force F drives the shaker's moving mass m_a (armature, screw, washer), which is
    fixed to the plate; the plate receives F / (1 − m_a ω² G(x₀,x₀)). Resonances sit between the free plate's (m_a → 0) and the
    plate pinned at the bolt (m_a → ∞).
  - *Rigid bolt*: the m_a → ∞ limit, W = W₀·G(x,x₀)/G(x₀,x₀).
  - The modes above the cut-off ω_K still add flexibility at the bolt. Up there a plate's point receptance tends to an
    infinite plate's, so G(x₀,x₀) gets their sum in closed form, ln((ω_K+ω)/(ω_K−ω)) / (8π ω √(Dρh)). Without it the
    bolted resonances drift up to 3 % high at the app's usual cut-off; with it they match exact point-support and
    point-mass solutions to 0.03 %.

  Peaks at different frequencies add in power; Γ = plate acceleration / g. The readout gives the deflection as a share of
  the thickness and warns above 20 %, where the linear theory used here stops being accurate.
- **Sand** – grains move only where Γ > 1, drifting down ∇Γ with random scatter, so they settle where the plate is quieter than gravity.
  Where the sand ends up follows from the vibration field; how fast it gets there is a phenomenological model (mobility, scatter).

## Checked against theory

`physics.test.js`: Bessel values against their power series; circular-plate eigenvalues against the closed-form frequency
equations (independently solved in scipy); square-plate Ritz eigenvalues against exact (simply supported) and published
(clamped, free) values; static centre deflection against plate theory (circle Fa²/16πD, square 0.0116 Fa²/D);
the Mindlin correction against exact Mindlin solutions (Hashemi & Arsanjani 2005, Table 9: within 2×10⁻⁵);
cantilever eigenvalues against Leissa NASA SP-160 Tables 4.51–4.52 (converged, at most 0.5 % below his 1969 upper bounds);
the peak finder on a real Blackman-windowed analyser spectrum (frequency within 0.05 Hz at 5.4 Hz bins, level within 0.05 dB).

## Checked against experiment

All of these run in `npm test` except the scale-free table and Leissa's Table 4.53.

### Square plates

Shu et al., [*Entropy* 24, 215 (2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8870825/), Fig. 2: a 320 × 320 × 1 mm aluminium
square screwed to a shaker at its centre, with 25 measured sand-pattern resonances from 263 to 4817 Hz. Simulated with the
paper's material values (E 70 GPa, ρ 2700, ν 0.33):

| Model | Free parameters | Resonances 200–5000 Hz | f1–f3 error | f4–f25 median | RMS, all 25 |
|---|---|---|---|---|---|
| Rigid bolt | none | 24 (measured: 25) | −20, −5, −7 % | 0.4 % | 4.6 % |
| Shaker, moving mass 20 g | m_a (not given in the paper; 16–25 g all fit) | 24 | −2.4, 0.0, +0.5 % | 0.4 % | 1.7 % |

Shifting the rigid model by ±3–10 % raises its f4–f25 median to 1.7–3.2 %, so the agreement is not chance matching.
Patterns compared by eye with the photos at f1–f7 match in topology (rings, diamonds, corner loops, the X at 697 Hz).

Further rigs, compared scale-free (one overall frequency factor fitted, which absorbs unstated thickness or stiffness),
against the RMS a wrong model reaches by chance (median over factors 0.7–1.3):

| Experiment | Plate | Best factor | RMS | By chance |
|---|---|---|---|---|
| Shu et al. 2022 (shaker, 20 g) | 320 × 1 mm Al, centre screw | 1.007 | 1.7 % | 5.5 % |
| D'Alessio, *SN Appl. Sci.* 3:60 (2021) (shaker, 20 g) | 240 × 0.9 mm Al, PASCO driver rod | 0.88 | 1.8 % | 13 % |
| Kovacic & Kanovic, *Symmetry* 15:1748 (2023) (rigid bolt) | 250 × 2 mm Al, LDS V408 shaker, material not stated | 0.82 | 1.7 % | 16 % |
| Lin et al., *Symmetry* 16:1460 (2024), off-centre drive | 320 × 1 mm Al | — | 3.6 % | 3.3 % |

The sequence of resonances matches in every case where the comparison can discriminate; the absolute level is
pinned down only where the plate is fully specified (Shu). Lin's off-centre drive excites so many modes that
nearest-peak matching is at chance level: it neither confirms nor refutes the model without matching each photographed pattern.
The mount model per row (shaker vs rigid) follows the driver used; with 4–5 frequencies those two rows are supporting evidence only.

### Cantilever plates

Leissa, NASA SP-160 (1969), pp. 82–86, measured frequency parameters ωa²√(ρh/D) of real cantilevers (one clamped edge):

| Table | Plates | Modes | Measured vs model |
|---|---|---|---|
| 4.55 | aluminium, a/b = ½, 1, 2, 5 | 14 | mean −1.3 %, worst −4.2 % |
| 4.57 | 6061-T6, 7.5 × 7.5 × 0.125 in (Craig, Plass & Caughfield) | 4 | mean −3.0 %, worst −4.1 % |
| 4.58 | steel, 10 × 10 × 0.0627 in (Gustafson, Stokey & Zorowski) | 5 | mean −1.9 %, worst −3.6 % |
| 4.56 | steel strips, a/b = 2–5, h/a = 0.02 (Plunkett & Wilson), modes up to λ ≈ 290 | 45 | median 1.0 %, mean −0.7 %, worst −5.4 % |
| 4.53 | mild steel 5.12 × 2.76 × 0.053 in, in Hz with handbook steel values | 22 | mean −3.1 %, worst −5.8 % |

All 68 nondimensional values (4.55–4.58) lie within −5.4 … +2.6 % and run a little below the model, as expected: a real
clamp is never perfectly rigid, which lowers the measured frequencies, most of all the fundamental. For 4.56 the Mindlin
correction halves the median error at a/b = 2 (1.0 → 0.5 %; chance level 5.2 %). For 4.53 Leissa's printed 1969 theory
misses the same measurements by 0.2–8.9 %. The app also reproduces the four sand photographs of a square cantilever in
Leissa's Fig. 4.45 in mode order: centre line (torsion), line near the free edge, two lines along the sides, a cross.
Two more measured sets (README only): Table 4.47's square cantilever runs 4–6 % below the model, and Heiba's ⅛-in steel
cantilevers (Fig. 4.46, 22 labelled modes, a/b = 1 … 0.4, handbook steel) average 1.8 % below it with 4 % scatter.

### Circular plates

- **Nodal circles, Strehlke (1855)** via Rayleigh, *Theory of Sound* §220: the exact solver reproduces Kirchhoff's calculated
  radii to within 2×10⁻⁴ r/a, and Strehlke's measured ones to within 0.003 r/a.
- **Nodal circles, Waller (1938)** via Leissa, NASA SP-160, Table 2.9: 60 measured circles on free brass plates
  (0–5 nodal diameters, 1–4 circles): median error 0.0014 r/a, worst 0.006.
- **Frequencies, Waller (1938)**, Table 2.8: 76 measured frequency ratios (up to 24 nodal diameters, 6 circles): median error 1.6 %.
- **Centre-driven Chladni plate**, Escaler & De La Torre, *J. Fluids Struct.* (2018): 240 × 0.8 mm aluminium, driven through
  a bar, moving mass 0.08 kg as stated by the authors (nothing fitted here). Free-plate frequencies agree with the
  authors' own theory within 0.1 %. All 6 measured resonances (70–4129 Hz) are reproduced, exactly 6 in band, at
  +1.8, +1.1, +0.5, +0.3, +0.1, −0.4 %. Measured nodal circles beyond r/a = 0.2 agree within 0.022 (mean 0.008); the
  innermost ones lie on the 8 mm drive bar, which the point mount does not resolve.

### Recovering an unknown rig from a photo

A widely shared photo of 12 square-plate patterns (345–5907 Hz) gives no plate details. Its centre dot and fourfold
symmetry mean a free square bolted at the centre; frequencies then depend only on h·√(E/ρ)/L² and the shaker's moving
mass relative to the plate. Fitting those two to the 12 frequencies gives RMS 1.3–1.5 % (chance: 6.7 %) for a moving
mass of 5–8 % of the plate, e.g. a 24 cm × 0.75 mm aluminium square on a ~7 g shaker (the "12-pattern photo plate"
preset). Rendered at the photo's frequencies, 345, 2041, 3835, 4840 and 5201 Hz match the photographs closely and most
of the rest in structure; the photo's 4444 Hz pattern is itself lopsided (a real plate's imperfection). The rigid-bolt
alternative fits the frequencies at 2.2 % but gives the wrong patterns, so the photos decide between them.

## Checked against Leissa, *Vibration of Plates* (NASA SP-160, 1969)

`leissa.test.js` and `leissa-ch4b.test.js` check every table in chapters 2 and 4 that covers what the app models, read from the
[scanned pages](https://archive.org/details/DTIC_ADA307623) (chapter 12's only isotropic table, Tomar's thick circular
plate, is one Leissa himself calls "highly questionable"). Where a printed value disagrees, an independent exact solution
(scipy) or a measurement decides, and in every such case it sides with the app:

| Leissa | Plate | Values | Agreement |
|---|---|---|---|
| 2.1, 2.2 | clamped circle: λ² up to 14 diameters / 9 circles, nodal radii | 64 + 88 | median 0.001 %; radii median 0.0003 r/a |
| 2.3, 2.4 | simply supported circle | 12 + 18 | 0.01 % (the table is labelled ν = 0.3 but is the ν = 1/3 roots) |
| 2.5, 2.6 | free circle, ν = 0.33 and 0.25 | 74 + 38 | median 0.13 % and 0.01 % |
| 2.7 | free-circle nodal radii | 88 | 67 within 0.005 r/a; in the other 15 measured cases Waller's measurement sides with the app |
| 2.12–2.15 | free / simply supported / clamped circle held at the centre; clamped circle with a central mass | 20 | exact solutions matched to 0.03 % through the app's rigid-bolt and shaker models |
| 4.3–4.8 | SS-C-SS-C rectangles, many aspect ratios | 176 | median 0.03 % |
| 4.10–4.20 | SS-C-SS-SS, SS-C-SS-F, SS-SS-SS-F, SS-F-SS-F rectangles | 275 | median 0.01–0.2 %; exact Lévy solutions matched to 0.01 % |
| 4.22–4.31 | C-C-C-C: Claassen & Thorne's grid (a/b = 1 … 0.02), Young, Iguchi, Ödman, mode veering | 497 + 107 | 490 of 497 within 0.01 %; the other 7 are 5 digit slips plus 0.11–0.15 % |
| 4.25, 4.26, 4.32 | C-C-C-C rigorous lower/upper bounds (Bazley, Fox & Stadler; Aronszajn), b/a = 1 … 8 | 110 | all inside |
| 4.33–4.44 | C-C-C-SS, C-C-SS-SS, C-C-F-F, C-F-C-F (λ and λ* tables), C-F-SS-F | ~240 | to the printed figures; beam-function Ritz upper bounds sit above ours |
| 4.45–4.50 | cantilevers: rigorous bounds for a/b = 0.125 … 8, Barton/Young upper and Gontkevich lower bounds | 70 + 20 | all inside |
| 4.61–4.76 | completely free: rigorous bounds (b/a = 1 … 8), Lemke at four ν, Iguchi | 80 + 24 + 19 | all bounds hold; Lemke's upper bounds 0.1–1.2 % above ours |
| 4.68 | Waller's measured free square brass plate, frequency ratios up to 147 | 89 | median 0.9 %, worst 4.4 % |

Misprints or errors in the book that the comparison exposed: Table 2.5 (1,7) prints 651.8 for ≈ 553; Table 2.7 prints
.925/.926/.993 where theory and Waller's measurement both give .898–.915; Table 2.15's "precise" central-mass roots are up
to 0.34 % low; Table 4.6 prints λ* = 24.047 at a/b = 1.5 for 25.047; Table 4.7 at a/b = 1 contradicts Table 4.3;
Table 4.17 is wrong at a/b = 0.6, 0.8 and 5 (by up to 7 %); Table 4.9, an older approximate source, is up to 3.4 % off
Table 4.6 and is not used; the SS-F-SS-F (m,1) values of Tables 4.19–4.20 are wrong, as Leissa notes. In the second half
of chapter 4: digit slips in Table 4.30 (32.9981 for 32.7980, 42.2530 for 42.3253, 75.5403 for 75.4503, 16.2699 for
22.3821, 96.0908 for 96.2324), Table 4.40 (λ* 18.93 for 17.93), Table 4.42 (λ* 2.3 for 2.4) and Tables 4.42/4.43 (61.2,
below both neighbours); Kanazawa & Kawai's non-square values (Tables 4.28, 4.35–4.37) run a uniform 2.1–2.3 % high
against Claassen & Thorne and the app; Dill & Pister's and Kanazawa & Kawai's C-C-C-SS a/b = 0.5 value (24.48) lies above
a converged Ritz value (24.14275, unchanged from 20 to 44 functions), which bounds the exact one from above; Iguchi's
Table 4.65 has two roots that cannot exist (65.368; 173.695, in a gap Table 4.75's rigorous bounds prove empty) and two
values above Table 4.67's rigorous upper bounds; Sigillito's "upper bound" 13.464 (Table 4.67) is below the true 13.468.
Not used: Bolotin's approximations (4.24), a different plate theory (4.47's theory row), Walton (4.59, no plate size),
Ödman's single-product estimates (4.64, 4.72), the Ritz single-term rows of 4.62/4.63 (median 0.5 % off), Pavlik
(4.69–4.71, units unresolved; ratios within 1–3 %) and Waller's rectangles (4.77, medians 1.2–3.4 % with unidentified outliers).

## Known limits

- **Clamp size**: the bolt is a point. A real washer or bar clamps a small disk; that is finer than the modal basis resolves,
  and it shows only in the nodal circle closest to the clamp.
- **Air and sand loading** are not modelled. Air adds mass mainly to the lowest modes of thin, light plates
  (for 1 mm aluminium, up to a few % at the fundamental, under 1 % for patterns with many nodal lines).
- **Linear theory**: no large-amplitude stiffening (see the deflection readout), no Faraday waves, no fine powder
  gathering at antinodes from air currents.
- **Near a free edge**, measured outermost nodal circles of high-order modes (Waller, 5–6 circles) sit 0.012–0.021 r/a
  further out than thin-plate theory, while every inner circle matches within 0.004.
- **Rectangular modes** are valid up to the frequency shown under "Plate modes"; a thicker or smaller plate covers higher notes.
- **Inputs dominate in practice**: frequencies scale with h·√(E/ρ), so a 2 % thickness or stiffness error shifts every
  resonance by 1–2 %. Measure your plate.
