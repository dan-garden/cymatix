// Plate material presets: E (GPa), ρ (kg/m³), ν, each read from the linked document. These are typical values for the
// named grade. A real sheet can differ by a few percent, so use the calibration tool when a measurement is available.
export const MATERIALS = {
  aluminum: {
    name: 'Aluminium 6061-T6',
    E: 68.9,
    rho: 2700,
    nu: 0.33,
    source: 'ASM Aerospace Specification Metals (ν estimated there)',
    url: 'https://www.aerospacemetals.com/wp-content/uploads/2023/06/Aluminum-6061-T6-6061-T651.pdf',
  },
  steel: {
    name: 'Mild steel 1020',
    E: 210,
    rho: 7800,
    nu: 0.29,
    source: 'Roylance, Mechanics of Materials, appendix',
    url: 'https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Mechanics_of_Materials_(Roylance)/07%3A_Appendices/7.01%3A_Material_Properties',
  },
  brass: {
    name: 'Cartridge brass C26000',
    E: 110.3,
    rho: 8530,
    nu: 0.33,
    source: 'Fisk Alloy C26000 datasheet (ν from its E and G)',
    url: 'https://fiskalloy.com/pdf/c26000/',
  },
  glass: {
    name: 'Soda-lime float glass',
    E: 72,
    rho: 2500,
    nu: 0.22,
    source: 'Pilkington Technical Bulletin ATS-129',
    url: 'https://www.pilkington.com/-/media/pilkington/site-content/usa/window-manufacturers/technical-bulletins/ats-129---properties-of-soda-lime-silica-float-glass.pdf',
  },
  acrylic: {
    name: 'Cast acrylic (PLEXIGLAS GS)',
    E: 3.3,
    rho: 1190,
    nu: 0.37,
    source: 'Röhm PLEXIGLAS GS/XT Technical Information 211-1 (short-term value)',
    url: 'https://igepa.ee/wp-content/uploads/2017/06/211-1-PLEXIGLAS-GS-XT-en.pdf',
  },
};
