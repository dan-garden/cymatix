# Cymatix

**See the shape of sound.** Cymatix turns any sound into Chladni figures: sand patterns on a simulated vibrating plate, computed live in your browser from published plate physics.

![Cymatix in action: clicking presets makes sand gather into rings, squares, stars and grids](docs/demo.gif)

## Features

- **Three sound sources:** an audio or video file, your microphone, or a test tone (sine, square, sawtooth or triangle).
- **Any plate:** a circle or rectangle of any size and thickness, in aluminium, steel, brass, glass, acrylic or your own material.
- **Real mountings:** free, clamped or simply supported edges, set per edge on rectangles.
- **Real drives:** a point force, a shaker or a rigid bolt, with a finite clamp size if you know it.
- **Real surroundings:** air, vacuum or water.
- **Real grains:** sand, fine sand, salt or lycopodium powder, moving at speeds calculated in real units.
- **Millions of grains:** runs on the GPU with WebGL2, and falls back to the CPU.
- **Calibration:** enter resonances you measured, and it fits the model to your plate.
- **Ten presets:** including a plate from a [published experiment](https://doi.org/10.3390/e24020215).
- **Mode explorer:** click any natural mode to hear its exact frequency and see its pattern.
- **Lightweight:** plain JavaScript with no runtime dependencies. Vite is the only development tool.

## Quick start

You need Node 20.19+ or 22.12+.

```bash
git clone https://github.com/dan-garden/cymatix.git
cd cymatix
npm install
npm run dev
```

Then open <http://localhost:5173> and click a preset.

- The microphone only works on `localhost` or HTTPS.
- Click the plate to move the driver.
- **Reset all** restores every default, and so does reloading the page.
- `npm test` runs the physics checks.
- `npm run build` writes the site to `dist/`, and `npm run preview` serves that build.

## Deploying

Import the repo into [Vercel](https://vercel.com). It detects Vite, so no settings or config file are needed.

## Gallery

Every pattern below is a preset in the app, captured at ×100 time-lapse once the grains had settled.

![Nine sand patterns from the app's presets](docs/gallery.png)

## How it works

1. **Read the sound.** The browser's FFT finds the strongest frequency peaks.
2. **Solve the plate.** The app computes the plate's natural modes from its shape, material and edges.
3. **Shake the plate.** Each peak drives the plate at the driver point. Adding up the modes gives the vibration everywhere.
4. **Move the sand.** Grains hop where the plate accelerates faster than gravity. Friction against the plate's bending surface pushes them towards the nodal lines, where they come to rest.

## The science

Every formula and where it comes from. Each module is a separate file, listed under [project files](#project-files).

| Part | Model | Source |
|---|---|---|
| Plate bending | Thin-plate theory, `D = Eh³ / 12(1 − ν²)` | [Kirchhoff 1850](https://doi.org/10.1515/crll.1850.40.51), [Leissa 1969](https://ntrs.nasa.gov/citations/19700009156) |
| Circular plates | Exact modes, `[J_n(λρ) + C·I_n(λρ)]·cos nθ` | [Kirchhoff 1850](https://doi.org/10.1515/crll.1850.40.51), [Leissa 1969](https://ntrs.nasa.gov/citations/19700009156) ch. 2 |
| Bessel functions | Integral forms with the trapezoid rule, and a power series where the integral loses digits | [DLMF §10.9](https://dlmf.nist.gov/10.9), [DLMF §10.32](https://dlmf.nist.gov/10.32), [Trefethen & Weideman 2014](https://doi.org/10.1137/130932132) |
| Rectangular plates | Rayleigh–Ritz with orthogonal Legendre polynomials, each edge independent | [Ritz 1909](https://doi.org/10.1515/crll.1909.135.1), [Bhat 1985](https://doi.org/10.1016/s0022-460x(85)80109-7), [Dickinson & Di Blasio 1986](https://doi.org/10.1016/s0022-460x(86)80310-8) |
| Eigenvalue solver | Householder tridiagonalisation, then QL, in worker threads for large plates | [Martin et al. 1968](https://doi.org/10.1007/bf02161841), [Bowdler et al. 1968](https://doi.org/10.1007/bf02166681), [JAMA](https://math.nist.gov/javanumerics/jama/) |
| Thick plates | True Mindlin modes: three-field Ritz on rectangles, exact three-wavenumber Bessel on circles, `κ² = π²/12` | [Mindlin 1951](https://doi.org/10.1115/1.4010217), [Mindlin et al. 1956](https://doi.org/10.1115/1.4011349) |
| Forced vibration | Modal sum with damping, `G(x,x₀) = Σ φ_k(x)·φ_k(x₀) / m(ω_k² − ω² + 2iζω_k·ω)` | [Bishop & Johnson 1960](https://www.cambridge.org/9781107402454) |
| Shaker mount | Force reaching a plate that carries moving mass `m_a`, `F / (1 − m_a·ω²·G₀)`, where `G₀ = G(x₀,x₀)` | Derived here from receptance coupling, [Bishop & Johnson 1960](https://www.cambridge.org/9781107402454) |
| Rigid bolt | The limit of large `m_a`, `W = W₀·G(x,x₀) / G₀` | Same derivation, checked against [Leissa 1969](https://ntrs.nasa.gov/citations/19700009156) Tables 2.12–2.15 |
| Finite clamp | A rigid disk held by smooth ring forces, with its own translation and tilt | Derived here, checked against exact annular plates and [Leissa 1969](https://ntrs.nasa.gov/citations/19700009156) Tables 2.24–2.31 |
| Modes above the cut-off | Added to `G₀`: `ln((ω_K+ω)/(ω_K−ω)) / (8πω√(Dρh))`, with a near-field form for a finite clamp | Derived here from plate modal density, [Cremer et al. 2005](https://doi.org/10.1007/b137728), [Skudrzyk 1980](https://doi.org/10.1121/1.384169) |
| Surrounding fluid | Coupled added mass per mode for an unbaffled plate, `ω_fluid = ω / √(1 + β_k)` | Derived here, checked against [Kwak 1991](https://doi.org/10.1115/1.2897209) as printed by [Escaler & De La Torre 2018](https://doi.org/10.1016/j.jfluidstructs.2018.07.017) |
| Physical damping (optional) | Thermoelastic, viscous and radiation losses per mode, for aluminium | [Lambourg et al. 2001](https://doi.org/10.1121/1.1354201) |
| Large amplitudes | Von Kármán stiffening reduced to a cubic term per mode, solved by harmonic balance | [Yamaki 1961](https://doi.org/10.1002/zamm.19610411204) via Leissa, [Thomas & Bilbao 2008](https://doi.org/10.1016/j.jsv.2008.04.014) |
| Spectrum | Blackman-window FFT | [W3C Web Audio API](https://www.w3.org/TR/webaudio/), [Harris 1978](https://doi.org/10.1109/proc.1978.10837) |
| Peak finding | Parabolic interpolation between FFT bins | [Smith, *Spectral Audio Signal Processing*](https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html) |
| Several tones | Different frequencies add in power for the plate's field. Grain drift shares out by velocity amplitude | Modelling choice for the field, grain rule from simulation here |
| Sand threshold | Grains leave the plate near `Γ = acceleration / g > 1`. Their drift starts near Γ ≈ 0.95 and is full by about 1.65 | [Bougie & Duckert 2011](https://doi.org/10.1103/physreve.83.011303), simulation here |
| Sand motion | Friction against the top surface: `u = −Φ(Γ)·(h/2)·ω·∇A`, with `A` the deflection amplitude and `Φ` about 0.7 | Derived here from a bouncing-grain simulation, checked against [Zhou et al. 2016](https://doi.org/10.1038/ncomms12764) |
| Fine powder | Air drag over the vibrating plate, plus acoustic streaming, which carry light grains to antinodes under strong drive | Derived here, direction from [Faraday 1831](https://doi.org/10.1098/rstl.1831.0018) |
| Crowding | A cell holds about three layers of grains. Above that, grains jostle outwards until it no longer overflows | Modelling choice, not measured |
| Materials | Each preset names a grade and links its datasheet | [materials.js](materials.js) |

## How accurate is it?

`npm test` runs 105 tests. The tables below are the summary. [docs/validation.md](docs/validation.md) has the full results, including errors found in the reference book.

### Against real measurements

| Study | Plate | Compared | Result | In `npm test` |
|---|---|---|---|---|
| [Shu et al. 2022](https://doi.org/10.3390/e24020215) | 320 × 1 mm aluminium square on a shaker | 25 resonances | 1.7 % RMS in vacuum, 1.9 % with air (chance: 5.5 %) | Yes |
| [Escaler & De La Torre 2018](https://doi.org/10.1016/j.jfluidstructs.2018.07.017) | 240 × 0.8 mm aluminium circle, in air | 6 resonances | 0.9 % RMS with or without air, nothing fitted | Yes |
| [Escaler & De La Torre 2018](https://doi.org/10.1016/j.jfluidstructs.2018.07.017) | The same plate in water | 5 resonances | 3.4 % RMS, against 7.0 % for the standard added-mass factor | Yes |
| [Escaler & De La Torre 2018](https://doi.org/10.1016/j.jfluidstructs.2018.07.017) | The same plate, with a 4 mm clamp | 5 innermost nodal circles | 0.0021 r/a RMS, against 0.033 for a point | Yes |
| [Waller 1938](https://doi.org/10.1088/0959-5309/50/1/306) | Free brass circles | 76 frequency ratios | Median 1.6 % | Yes |
| [Waller 1938](https://doi.org/10.1088/0959-5309/50/1/306) | Free brass circles | 60 nodal circles | Median 0.0014 r/a | Yes |
| [Strehlke 1855](https://doi.org/10.1002/andp.18551710805) | Free circles | Nodal circles | Within 0.003 r/a | Yes |
| [Waller 1939](https://doi.org/10.1088/0959-5309/51/5/312) | Free brass square | 89 frequency ratios | Median 0.9 %, worst 4.4 % | Yes |
| [Leissa 1969](https://ntrs.nasa.gov/citations/19700009156) Tables 4.55–4.58 | Real cantilevers | 68 frequencies | −5.4 % to +2.6 % | Yes |
| [Zhou et al. 2016](https://doi.org/10.1038/ncomms12764) | Silicon plate with flattened solder balls | Voltage for a fixed grain displacement, 55 notes | Factor 1.5–1.8 on average, 1.9–3 note to note (shuffled notes: 5.7) | Yes |
| [Lambourg et al. 2001](https://doi.org/10.1121/1.1354201) | Free 4 mm aluminium plate | 21 decay rates | Factor 1.17 (best single damping value: 1.4) | Yes |
| [D'Alessio 2021](https://doi.org/10.1007/s42452-020-04062-6) | 240 × 0.9 mm aluminium square | 4–5 resonances | 1.8 % RMS after fitting one factor (chance: 13 %) | No |
| [Kovacic & Kanovic 2023](https://doi.org/10.3390/sym15091748) | 250 × 2 mm square | 4–5 resonances | 1.7 % RMS after fitting one factor (chance: 16 %) | No |
| [Lin et al. 2024](https://doi.org/10.3390/sym16111460) | 320 × 1 mm aluminium square, off-centre drive | Resonances | 3.6 % RMS (chance: 3.3 %), so **inconclusive** | No |

How to read this:

- **"Chance"** is the error a wrong model reaches by accident. A result only means something if it beats chance.
- **Shu et al.:** one value is fitted, the shaker's moving mass, because the paper doesn't give it: 20 g in vacuum, 15 g with air. Without any fit, a rigid bolt gives 4.6 % RMS.
- **D'Alessio and Kovacic:** the fitted factor absorbs unknown thickness or stiffness. With so few frequencies, these rows are supporting evidence only.
- **Zhou et al.:** the grains are flattened solder balls on silicon, not sand. This is the only absolute speed check, and it is not a check of sand.
- **Decay rates:** a second free plate ([Ege et al.](https://arxiv.org/abs/1210.3333)) measured 5.7 times more damping than the model, so physical damping stays an option.
- **Patterns:** the Shu et al. patterns for f1–f7 were compared with the photos by eye. They match in topology, but this is not a numerical test.

### Against theory

| Reference | Compared | Result |
|---|---|---|
| [Leissa 1969](https://ntrs.nasa.gov/citations/19700009156) ch. 2 and 4 | Every table covering plates the app models | Typical median 0.001 % to 0.2 % |
| Rigorous upper and lower bounds, listed in Leissa | Clamped, cantilever and free rectangles | All inside the bounds |
| [Hashemi & Arsanjani 2005](https://doi.org/10.1016/j.ijsolstr.2004.06.063) | Exact thick-plate solutions | Within 2×10⁻⁵ |
| Own 3D elasticity solver | Thick free and clamped plates at h/a = 0.1 | 0.5–1.7 % worst case, against 0.9–5.2 % for the old correction |
| Own exact Mindlin solver, independent code | Five edge combinations, first 12 modes | Within 5×10⁻⁶ % |
| Own exact annular solver | Finite clamps, r₀/a up to 0.1 | Within 0.35 % |
| Exact point-support and point-mass solutions, Leissa Tables 2.12–2.15 | Bolted plate models | Within 0.03 % |
| Own vortex-lattice and vortex-ring solvers | Fluid added mass | Within 0.03 % on smooth modes, 1 % on a cantilever |
| [Yamaki 1961](https://doi.org/10.1002/zamm.19610411204) via Leissa, time integration | Stiffening coefficients and backbone curves | Within 0.5 % and 1.5 % |
| [Timoshenko & Woinowsky-Krieger 1959](https://archive.org/details/theoryofplatessh0000timo_r3s0) | Static centre deflection of a clamped circle and a simply supported square | Within 0.5 % and 1.5 % |

## Known limits

What the model still cannot claim:

- **Sand speed is calculated, but only partly checked.** It comes from mechanics in real units, with no tuned values. The only absolute check is on solder balls, not sand. Restitution and friction are typical values, and the random spread is untested and uncertain by about ×3.
- **Powder is the least certain part.** It reaches antinodes under strong drive, as [Faraday 1831](https://doi.org/10.1098/rstl.1831.0018) saw. Its speed is unvalidated, and its direction is unreliable below 30 µm or below Γ ≈ 100.
- **Grain contact is approximate.** Crowding only limits how many grains fit in a cell, so line widths at large sand masses are not validated.
- **Stiffening covers point-force mounts only.** It matches exact theory, but no measurement. On a bolt or shaker the app warns instead, and past its range it says the model has broken down.
- **Granular waves are not simulated.** They need a deep layer at low frequency ([Melo et al. 1995](https://doi.org/10.1103/physrevlett.75.3838)), outside normal Chladni conditions. The app warns if you get close.
- **Air cannot be confirmed from air data.** It lowers frequencies by at most 1.6 % on a 320 mm × 1 mm square, and under 1 % above about 150 Hz. That is no larger than the uncertainty in a real plate's properties, so it is included because it is physically there, not because it improves the fit.
- **Sand mass is not applied.** Resting grains shift frequencies, but they collect where the plate barely moves: 10 g on the Escaler plate gives −5.2 % spread evenly and under −1 % on the nodal lines. Bouncing grains add less, as about 1/Γ ([Ferreyra et al. 2021](https://doi.org/10.1007/s10035-021-01090-w), abstract only).
- **A finite clamp is an input, not a discovery.** It shifts axisymmetric modes by 1–22 % as r₀/a goes from 0.02 to 0.1. On Escaler's plate a 4 mm clamp predicts the inner circles but moves frequencies to +3 to +5 %, which would need a plate about 5 % thinner or softer. That is unverified, so the default stays a point.
- **Thick plates are Mindlin, not 3D.** The remaining error against 3D elasticity is up to 1.7 % at h/a = 0.1 and grows with thickness.
- **Two of Waller's nodal circles don't fit.** Her outermost circles at 5–6 circles sit 0.012–0.021 r/a beyond theory. Exact Mindlin and 3D rule out thickness: they move these by under 0.0002 r/a. The two values lie within 0.0034 of an old approximate formula (Leissa eq. 2.19), so they may have been computed rather than measured. Only her paper can settle that, and it could not be read.
- **Very large thin plates still have a frequency limit.** A 1 m × 1 mm square now reaches 6.8 kHz, up from 900 Hz, after about a minute in worker threads.
- **Bolted resonances depend slightly on the mode count.** Escaler's sixth resonance moves by 0.6 % as more modes are added, even with the high-mode term.
- **Your inputs still matter most.** A 2 % thickness or stiffness error shifts every resonance by 1–2 %. Measure your plate, or use the calibration tool.
- **Material values are typical for each grade.** A real sheet can differ by a few percent. Acrylic's value is a short-term static one.
- **Quick-pick tones make no health claim.** The Solfeggio, 432 Hz and 440 Hz options are only convenient frequencies. The pattern depends on the plate, not on any special frequency.

## Project files

- [index.html](index.html): the page and its controls.
- [app.js](app.js): audio, drawing and the main loop.
- [physics.js](physics.js): thin-plate modes and the driven plate.
- [plateModel.js](plateModel.js): the full plate, with each effect below optional.
- [mindlin.js](mindlin.js): true thick-plate modes.
- [fluid.js](fluid.js): fluid added mass, radiation, physical damping and sand mass.
- [patch.js](patch.js): finite clamps.
- [nonlinear.js](nonlinear.js): large-amplitude stiffening.
- [grains.js](grains.js): sand and powder transport.
- [sandGL.js](sandGL.js): the GPU grain engine.
- [calibrate.js](calibrate.js): fitting a plate to measurements.
- [materials.js](materials.js): material presets with sources.
- [solvePlan.js](solvePlan.js), [solverPool.js](solverPool.js), [solver.worker.js](solver.worker.js): solving off the main thread.
- `*.test.js` and [testkit.js](testkit.js): the checks against theory, measurements and Leissa's book.

## References

### Plate theory

- Kirchhoff (1850). [Über das Gleichgewicht und die Bewegung einer elastischen Scheibe](https://doi.org/10.1515/crll.1850.40.51). *J. reine angew. Math.* 40, 51–88.
- Rayleigh (1894). [The Theory of Sound, vol. 1, 2nd ed.](https://archive.org/details/theorysound07raylgoog) Macmillan.
- Ritz (1909). [Über eine neue Methode zur Lösung gewisser Variationsprobleme der mathematischen Physik](https://doi.org/10.1515/crll.1909.135.1). *J. reine angew. Math.* 135, 1–61.
- Mindlin (1951). [Influence of rotatory inertia and shear on flexural motions of isotropic, elastic plates](https://doi.org/10.1115/1.4010217). *J. Appl. Mech.* 18, 31–38.
- Mindlin, Schacknow & Deresiewicz (1956). [Flexural vibrations of rectangular plates](https://doi.org/10.1115/1.4011349). *J. Appl. Mech.* 23, 430–436.
- Chu & Herrmann (1956). [Influence of large amplitudes on free flexural vibrations of rectangular elastic plates](https://doi.org/10.1115/1.4011396). *J. Appl. Mech.* 23, 532–540.
- Timoshenko & Woinowsky-Krieger (1959). [Theory of Plates and Shells, 2nd ed.](https://archive.org/details/theoryofplatessh0000timo_r3s0) McGraw-Hill.
- Bishop & Johnson (1960). [The Mechanics of Vibration](https://www.cambridge.org/9781107402454). Cambridge University Press.
- Leissa (1969). [Vibration of Plates, NASA SP-160](https://ntrs.nasa.gov/citations/19700009156). Also as a [full scan](https://archive.org/details/DTIC_ADA307623).
- Thomas & Bilbao (2008). [Geometrically nonlinear flexural vibrations of plates: in-plane boundary conditions and some symmetry properties](https://doi.org/10.1016/j.jsv.2008.04.014). *J. Sound Vib.* 315, 569–590.
- Yamaki (1961). [Influence of large amplitudes on flexural vibrations of elastic plates](https://doi.org/10.1002/zamm.19610411204). *ZAMM* 41, 501–510. Read through Leissa Sec. 12.2.
- Kwak (1991). [Vibration of circular plates in contact with water](https://doi.org/10.1115/1.2897209). *J. Appl. Mech.* 58, 480–483. Values as printed by Escaler & De La Torre.
- Lambourg, Chaigne & Matignon (2001). [Time-domain simulation of damped impacted plates. II. Numerical model and results](https://doi.org/10.1121/1.1354201). *J. Acoust. Soc. Am.* 109, 1433–1447.
- Skudrzyk (1980). [The mean-value method of predicting the dynamic response of complex vibrators](https://doi.org/10.1121/1.384169). *J. Acoust. Soc. Am.* 67, 1105–1135.
- Bhat (1985). [Natural frequencies of rectangular plates using characteristic orthogonal polynomials in Rayleigh-Ritz method](https://doi.org/10.1016/s0022-460x(85)80109-7). *J. Sound Vib.* 102, 493–499.
- Dickinson & Di Blasio (1986). [On the use of orthogonal polynomials in the Rayleigh-Ritz method for the study of the flexural vibration and buckling of isotropic and orthotropic rectangular plates](https://doi.org/10.1016/s0022-460x(86)80310-8). *J. Sound Vib.* 108, 51–62.
- Hashemi & Arsanjani (2005). [Exact characteristic equations for some of classical boundary conditions of vibrating moderately thick rectangular plates](https://doi.org/10.1016/j.ijsolstr.2004.06.063). *Int. J. Solids Struct.* 42, 819–853.
- Cremer, Heckl & Petersson (2005). [Structure-Borne Sound, 3rd ed.](https://doi.org/10.1007/b137728) Springer.

### Experiments

- Chladni (1787). [Entdeckungen über die Theorie des Klanges](https://www.digitale-sammlungen.de/en/view/bsb10908815). Weidmanns Erben und Reich.
- Strehlke (1855). [Ueber die Schwingungen homogener elastischer Scheiben](https://doi.org/10.1002/andp.18551710805). *Ann. Phys.* 171, 577–602. The radii are quoted in Rayleigh §220.
- Waller (1938). [Vibrations of free circular plates. Part 1: Normal modes](https://doi.org/10.1088/0959-5309/50/1/306). *Proc. Phys. Soc.* 50, 70–76. Leissa Ref. 2.21.
- Waller (1939). [Vibrations of free square plates: part I. Normal vibrating modes](https://doi.org/10.1088/0959-5309/51/5/312). *Proc. Phys. Soc.* 51, 831–844. Leissa Ref. 4.119.
- Ege, Boncompagne, Laulagnet & Guyader (2012). [Experimental estimations of viscoelastic properties of multilayer damped plates in broad-band frequency range](https://arxiv.org/abs/1210.3333). arXiv preprint, used only for its measured aluminium loss factor.
- Escaler & De La Torre (2018). [Axisymmetric vibrations of a circular Chladni plate in air and fully submerged in water](https://doi.org/10.1016/j.jfluidstructs.2018.07.017). *J. Fluids Struct.* 82, 432–445.
- D'Alessio (2021). [Forced free vibrations of a square plate](https://doi.org/10.1007/s42452-020-04062-6). *SN Appl. Sci.* 3, 60.
- Shu, Tseng, Lai, Yu, Huang & Chen (2022). [Exploring the origin of maximum entropy states relevant to resonant modes in modern Chladni plates](https://doi.org/10.3390/e24020215). *Entropy* 24, 215.
- Kovacic & Kanovic (2023). [Chladni plate in anechoic chamber: symmetry in vibrational and acoustic response](https://doi.org/10.3390/sym15091748). *Symmetry* 15, 1748.
- Lin, Hsu, Su, Liang & Chen (2024). [Exploring the multiplication of resonant modes in off-center-driven Chladni plates from maximum entropy states](https://doi.org/10.3390/sym16111460). *Symmetry* 16, 1460.

### Sand and granular motion

- Melo, Umbanhowar & Swinney (1995). [Hexagons, kinks, and disorder in oscillated granular layers](https://doi.org/10.1103/physrevlett.75.3838). *Phys. Rev. Lett.* 75, 3838–3841.
- Umbanhowar & Swinney (2000). [Wavelength scaling and square/stripe and grain mobility transitions in vertically oscillated granular layers](https://doi.org/10.1016/s0378-4371(00)00432-5). *Physica A* 288, 344–362.
- Dorrestijn et al. (2007). [Chladni figures revisited based on nanomechanics](https://doi.org/10.1103/physrevlett.98.026102). *Phys. Rev. Lett.* 98, 026102.
- van Gerner, van der Hoef, van der Meer & van der Weele (2010). [Inversion of Chladni patterns by tuning the vibrational acceleration](https://doi.org/10.1103/physreve.82.012301). *Phys. Rev. E* 82, 012301. Abstract only.
- Ferreyra, Baldini, Pugnaloni & Job (2021). [Effect of lateral confinement on the apparent mass of granular dampers](https://doi.org/10.1007/s10035-021-01090-w). *Granul. Matter* 23, 45. Abstract only.
- Loubet, Jarosz, Saint-Jean & Huber (2007). [A method for measuring the settling velocity distribution of large biotic particles](https://doi.org/10.1007/s10453-007-9054-2). *Aerobiologia* 23, 159–169. Abstract only, for lycopodium's settling speed.
- Faraday (1831). [On a peculiar class of acoustical figures, and on certain forms assumed by groups of particles upon vibrating elastic surfaces](https://doi.org/10.1098/rstl.1831.0018). *Phil. Trans. R. Soc.* 121, 299–340.
- Bougie & Duckert (2011). [Continuum simulations of shocks and patterns in vertically oscillated granular layers](https://doi.org/10.1103/physreve.83.011303). *Phys. Rev. E* 83, 011303. [Free copy on arXiv](https://arxiv.org/abs/1009.2742).
- Arango & Reyes (2016). [Stochastic models for Chladni figures](https://doi.org/10.1017/s0013091515000139). *Proc. Edinb. Math. Soc.* 59, 287–300.
- Zhou, Sariola, Latifi & Liimatainen (2016). [Controlling the motion of multiple objects on a Chladni plate](https://doi.org/10.1038/ncomms12764). *Nat. Commun.* 7, 12764.

### Numerics and signal processing

- Martin, Reinsch & Wilkinson (1968). [Householder's tridiagonalization of a symmetric matrix](https://doi.org/10.1007/bf02161841). *Numer. Math.* 11, 181–195.
- Bowdler, Martin, Reinsch & Wilkinson (1968). [The QR and QL algorithms for symmetric matrices](https://doi.org/10.1007/bf02166681). *Numer. Math.* 11, 293–306.
- Harris (1978). [On the use of windows for harmonic analysis with the discrete Fourier transform](https://doi.org/10.1109/proc.1978.10837). *Proc. IEEE* 66, 51–83.
- Trefethen & Weideman (2014). [The exponentially convergent trapezoidal rule](https://doi.org/10.1137/130932132). *SIAM Rev.* 56, 385–458.
- Hicklin, Moler, Webb et al. [JAMA: A Java Matrix Package](https://math.nist.gov/javanumerics/jama/). NIST and MathWorks.
- NIST. [Digital Library of Mathematical Functions, §10.9](https://dlmf.nist.gov/10.9) and [§10.32](https://dlmf.nist.gov/10.32).
- Smith. [Quadratic interpolation of spectral peaks](https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html), in *Spectral Audio Signal Processing*.
- W3C. [Web Audio API](https://www.w3.org/TR/webaudio/).

### Data read through Leissa's book

These values were read from [Leissa 1969](https://ntrs.nasa.gov/citations/19700009156), not from the originals. "Ref." is the reference number in the book. Reports with no online copy link to the book's listing.

<details>
<summary>Show all 21 sources</summary>

- Lemke (1928). [Experimentelle Untersuchungen zur W. Ritzschen Theorie der Transversalschwingungen quadratischer Platten](https://doi.org/10.1002/andp.19283911304). *Ann. Phys.* 391, 717–750. Ref. 4.110.
- Iguchi (1938). Die Eigenwertprobleme für die elastische rechteckige Platte. *Mem. Fac. Eng., Hokkaido Univ.*, 305–372. [Ref. 4.9](https://archive.org/details/DTIC_ADA307623).
- Aronszajn (1950). The Rayleigh-Ritz and the Weinstein methods for approximation of eigenvalues, part III. Tech. Rept. 3, Oklahoma A&M College. [Ref. 4.49](https://archive.org/details/DTIC_ADA307623).
- Young (1950). [Vibration of rectangular plates by the Ritz method](https://doi.org/10.1115/1.4010175). *J. Appl. Mech.* 17, 448–453. Ref. 4.47.
- Barton (1951). [Vibration of rectangular and skew cantilever plates](https://doi.org/10.1115/1.4010265). *J. Appl. Mech.* 18, 129–134. Ref. 4.73.
- Kanazawa & Kawai (1952). On the lateral vibration of anisotropic rectangular plates. *Proc. 2nd Jap. Natl. Congr. Appl. Mech.*, 333–338. [Ref. 4.16](https://archive.org/details/DTIC_ADA307623).
- Iguchi (1953). [Die Eigenschwingungen und Klangfiguren der vierseitig freien rechteckigen Platte](https://doi.org/10.1007/bf00535853). *Ingenieur-Archiv* 21, 303–322. Ref. 4.113.
- Wilson & Plunkett (1953). Vibration of cantilever plates with rectangular and wedge-shaped cross-sections. Report DF 53GL17, General Electric. [Ref. 4.85](https://archive.org/details/DTIC_ADA307623).
- Gustafson, Stokey & Zorowski (1954). [The effect of tip removal on the natural vibrations of uniform cantilevered triangular plates](https://doi.org/10.2514/8.3148). *J. Aeronaut. Sci.* 21, 621–633. Ref. 4.87.
- Heiba (1954). Vibration characteristics of a cantilever plate with sweptback leading edge. Report 82, Cranfield College of Aeronautics. [Ref. 4.83](https://archive.org/details/DTIC_ADA307623).
- Ödman (1955). Studies of boundary value problems, part II: characteristic functions of rectangular plates. Proc. NR 24, Swedish Cement and Concrete Research Institute. [Ref. 4.13](https://archive.org/details/DTIC_ADA307623).
- Dill & Pister (1958). Vibration of rectangular plates and plate systems. *Proc. 3rd U.S. Natl. Congr. Appl. Mech.*, 123–132. [Ref. 4.24](https://archive.org/details/DTIC_ADA307623).
- Gontkevich (1960). The lower bounds of natural frequencies of a plate in flexural vibrations. *Prikl. Mekh.* 6, 346–350. [Ref. 4.55](https://archive.org/details/DTIC_ADA307623).
- Claassen & Thorne (1961). [Vibrations of thin rectangular isotropic plates](https://doi.org/10.1115/1.3641676). *J. Appl. Mech.* 28, 304–305. Refs. 4.35 and 4.36.
- Plunkett (1963). [Natural frequencies of uniform and non-uniform rectangular cantilever plates](https://doi.org/10.1243/jmes_jour_1963_005_020_02). *J. Mech. Eng. Sci.* 5, 146–156. Ref. 4.84.
- Craig, Plass & Caughfield (1964). Experimental determination of frequencies and mode shapes of cantilever and hub-pin plates. Report DRL-518, University of Texas. [Ref. 4.86](https://archive.org/details/DTIC_ADA307623).
- Gontkevich (1964). Natural Vibrations of Plates and Shells. Naukova Dumka. [Ref. 4.51](https://archive.org/details/DTIC_ADA307623).
- Bazley, Fox & Stadter (1965). Upper and lower bounds for frequencies of rectangular cantilever plates. Tech. Memo TG-705, Johns Hopkins APL. [Ref. 4.75](https://archive.org/details/DTIC_ADA307623).
- Sigillito (1965). Improved upper bounds for frequencies of rectangular free and cantilever plates. Eng. Memo EM-4012, Johns Hopkins APL. [Ref. 4.76](https://archive.org/details/DTIC_ADA307623).
- Bazley, Fox & Stadter (1967). [Upper and lower bounds for the frequencies of rectangular clamped plates](https://doi.org/10.1002/zamm.19670470307). *ZAMM* 47, 191–198. Journal version of Ref. 4.58.
- Bazley, Fox & Stadter (1967). [Upper and lower bounds for the frequencies of rectangular free plates](https://doi.org/10.1007/bf01601716). *ZAMP* 18, 445–460. Journal version of Ref. 4.118.

</details>
