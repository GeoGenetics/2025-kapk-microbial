# Supplementary Tables & Extended Data Figures — Manuscript Integration

## WHERE TO ADD SUPPLEMENTARY TABLE CITATIONS IN MAIN TEXT

---

### Supplementary Table S4 — Taxonomy and damage (full table)
**Currently commented out in R script — must be uncommented.**

Add to **Para 45** (Community Structure section):
> "We identified a diverse set of bacterial, archaeal and viral taxa in the damage-supported depositional component **(Supplementary Table S4)**."
✅ Already cited — just ensure the table is regenerated.

Add to **Para 48** (archaeal assemblage):
> "The archaeal assemblage included several methanogenic lineages characteristic of wetland ecosystems, including hydrogenotrophic and acetoclastic methanogens as well as taxa capable of methane oxidation **(Fig. 4A; Supplementary Table S4)**."
✅ Already cited.

Add to **Para 52** (carbon/nitrogen cycling):
> "Beyond methane cycling, the reconstructed community included taxa involved in carbon fermentation and nitrogen transformations **(Supplementary Table S4)**."
✅ Already cited.

---

### Supplementary Table S12 — Custom KEGG module definitions
**Table does not exist yet. Needs new sup_table_4.xlsx with module definitions from user_modules/ on denbi-h-micro.**

Add to **Methods — Carbon metabolism** (already drafted):
> "We expanded the standard KEGG modules collection using anvi-setup-user-modules to include custom modules **(Supplementary Table S12**; defined as Woodcroft2018 and Borrel2023 in the module class) designed to explore carbon metabolism in permafrost."
✅ Already in methods text — table needs to be created.

---

### Supplementary Tables S13–S16 — KEGG and CAZy DART results
**Need to be added to Results — Carbon metabolism section (Para 63–65 area). Currently NO citation of KEGG/CAZy result tables in main text.**

Suggested insertion after the functional heatmap description:
> "The complete KEGG module coverage and damage statistics across all samples are provided in **Supplementary Table S13**; modules with authenticated ancient damage (posterior ≥ 0.7) are listed in **Supplementary Table S14**. CAZy family assignments are provided in **Supplementary Tables S15–S16**."

---

### Supplementary Table S17 — Viral references (DART authenticated)
Add to **Para 80** (virome section):

Current text:
> "The amino acid searches also provided information about the functional potential of the viruses we recruited **(Supplementary Table S17)**."
✅ Already cited — table needs regeneration with DART output.

Also add to **Para 54** (virome abundance):
> "We detected 1,352 authenticated viral references across all samples **(Supplementary Table S17)**."
⚠️ **Not currently cited here — add citation.**

---

### Supplementary Table S18 — IMGVR annotation
**Not cited anywhere in main text.**

Add to **Para 81** (ecosystem origins of virome):
> "We leveraged the ecosystem information associated with each reference in the IMG/VR data to gain insights into the environmental origins of the references **(Fig. 6E; Supplementary Table S18)**."

---

### Supplementary Tables S19–S23 — Simulations and benchmarks
Add to **Methods — Benchmarking** (if section exists) or Supplementary Information:
> "Parameter optimisation and benchmark results are provided in **Supplementary Tables S19–S23**."
✅ Already cited as S18–S22 — renumber to S19–S23 after adding new S12.

---

## EXTENDED DATA FIGURES — WHAT IS NEEDED

### Existing (ED Figs 1–9) — keep as-is
| ED Fig | Content | Status |
|--------|---------|--------|
| ED Fig. 1 | Taxonomic database workflow | ✅ keep |
| ED Fig. 2 | Taxonomic + functional profiling workflow | ⚠️ **update orange panel to show DART instead of xFilter** |
| ED Fig. 3 | Extraction controls | ✅ keep |
| ED Fig. 4 | Damage thresholds and sample selection | ✅ keep |
| ED Fig. 5 | Extraction bias and bloom assessment | ✅ keep |
| ED Fig. 6 | Read alignment saturation | ✅ keep |
| ED Fig. 7 | Taxonomic classification assessment | ✅ keep |
| ED Fig. 8 | Taxonomic abundance assessment | ✅ keep |
| ED Fig. 9 | Damage estimate assessment | ✅ keep |

### New Extended Data Figures needed

**ED Fig. 10 — DART damage authentication validation**
- Content: damage distribution across KEGG/CAZy/viral references; bimodal distribution showing 0.5 vs ≥0.7 posterior; DART AUC-ROC on synthetic benchmark
- Available figures: `results/functional_agp/fig_damage_distribution.pdf`
- Caption draft: "Distribution of DART posterior damage probability across KEGG, CAZy and viral protein references in all Kap København samples. The bimodal distribution reflects unauthenticated (posterior ≈ 0.5) and authentically ancient (posterior ≥ 0.7) references. The dashed line indicates the authentication threshold applied throughout this study."

**ED Fig. 11 — Key enzyme dotplot**
- Content: dotplot of key enzymes/KOs across samples with damage authentication overlay
- Available figures: `results/functional_agp/fig_key_enzyme_dotplot.pdf`
- Caption draft: "Coverage and ancient DNA authentication of key metabolic enzymes across Kap København samples. Each dot represents a KO detected in a sample; dot size reflects average protein coverage (log₂ scale) and colour indicates DART posterior damage probability. Only enzymes with authenticated ancient damage (posterior ≥ 0.7) are shown."

**ED Fig. 12 — Methanoflorentales BEAST2 molecular clock**
- Content: phylogenetic tree with node ages, rate estimates, calibration
- Available: `analysis/beast2/` directory (untracked)
- Caption draft: "Bayesian molecular clock analysis of Methanoflorentales. Divergence times estimated under a relaxed lognormal clock (preferred over strict clock; 2 × ln BF = +1,914) using the two-million-year Kap København constraint. Node bars show 95% HPD intervals."

---

## SUMMARY OF ACTIONS NEEDED

### In R script (14-summary-tables.R)
- [ ] Uncomment S4 (taxonomy and damage)
- [ ] Add new sup_table_4: custom KEGG module definitions (S12)
- [ ] Renumber all S-references +1 from S13 onwards
- [ ] Regenerate sup_table_5 (KEGG/CAZy) and sup_table_6 (Viral) with DART outputs
- [ ] Verify kegg_module_damage.tsv + cazy_family_damage.tsv exist locally

### In main text
- [ ] Para 54: add **(Supplementary Table S17)** to viral reference count sentence
- [ ] Para 63–65 area: add citations to S13–S16 for KEGG/CAZy results
- [ ] Para 81: add **(Supplementary Table S18)** for IMGVR ecosystem annotation
- [ ] Add ED Fig. 10 citation where DART authentication is described (Methods or Results virome/functional)
- [ ] Add ED Fig. 11 citation in carbon metabolism Results section
- [ ] Add ED Fig. 12 citation in Methanoflorentales molecular clock paragraph (Para 76)

### In Extended Data
- [ ] Update ED Fig. 2 functional workflow panel to show DART
- [ ] Create ED Fig. 10 (DART validation)
- [ ] Create ED Fig. 11 (KEGG enrichment)
- [ ] Create ED Fig. 12 (BEAST2 clock) — if BEAST2 results are finalised
