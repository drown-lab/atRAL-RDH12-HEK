# Manuscript revision to-do — atRAL / RDH12 multiomics

Audit of `atRAL_manuscript_draft_V11.docx` against the underlying data in this repository.
Every number below was re-derived from the source tables; the file and column are named so each
item can be checked independently.

Priority key: **P1** = factual error, must fix. **P2** = statistical framing, must address.
**P3** = qualification or softening. **P4** = housekeeping.

---

## 1. Must-fix factual errors (P1)

| # | Location | Problem | Correct value / action |
|---|---|---|---|
| 1.1 | Results, ferroptosis section; Discussion, ferrostatin paragraph | RDH12 + ferrostatin-1 IC50 given as "~285 µM" | Fitted value is **351.5 µM (95% CI 330.8–372.2)**. Source: `Cell_Viability/figures/atRAL_IC50_estimates.csv`. 285 appears nowhere in the data or the per-replicate fits. |
| 1.2 | Results, ferroptosis section vs Discussion | The two sections interpret the ferrostatin experiment in **opposite** directions. Results say ferrostatin protection shows ferroptosis is a critical executioner and that RDH12 acts independently of lipid peroxidation. Discussion says ferrostatin failed to protect WT cells and that ferroptosis emerges only once clearance is intact. | Data support the **Discussion** reading. Rewrite the Results passage; see §5.2. |
| 1.3 | Abstract | "Acute atRAL exposure triggered antioxidant defense programs by upregulating NRF2 and HMOX1" | NFE2L2 is **not significant at either acute dose** (adjusted p = 0.31 at 100 µM, 0.26 at 200 µM). HMOX1 rises only at 100 µM (adjusted p < 0.01) and is flat at 200 µM (adjusted p = 0.54). NRF2 induction is a **recovery** finding. See §5.1. |
| 1.4 | Figure 3E legend | States "fold change ±1.5 and BH-adjusted p < 0.01" | The generating script uses **±1.4-fold and adjusted p < 0.1** (`Proteomic_Rscripts/Figures/HeatMap_Fig1_horizontal_maintext.R`, lines 27–28). Several plotted proteins, for example PFDN4 at adjusted p 0.038 and PFDN5 at 0.060, fail the stated cutoff. The Figure 5E legend is already correct. |
| 1.5 | Discussion, glycerophospholipid paragraph | "we did not observe changes in cytosolic or calcium-independent phospholipase A2" | **PLA2G4A, cytosolic phospholipase A2 alpha, increases significantly at 200 µM** (log2FC 0.67, adjusted p 0.010) and is quantified in three of three replicates in both groups. PLA2G15 also rises (0.42, adjusted p 0.030). PLA2G6 was never quantified. See §5.5. |
| 1.6 | Discussion, glycerophospholipid paragraph | "LPCAT3 … was not among the altered proteins" | LPCAT3 has **no row in the results table**; it was never quantified. Absence of change cannot be distinguished from absence of measurement. |
| 1.7 | Results, recovery genotype section | Control-enriched proteins "associated with … extracellular matrix organization, cell adhesion …" | Neither term is enriched. Extracellular matrix organization appears in **none** of the 186 raw GO terms; cell-adhesion terms sit at adjusted p of 0.166 or worse. Heme biosynthesis is marginal at 0.052. Source: `GOanalysis/output/clusterProfiler_batch/RDH12_200_atRAL5hr.24h_recvr_vs_GFP_200_atRAL5hr.24h_recvr_down/`. See §5.7. |
| 1.8 | Results, recovery lipid section | "the exclusive detection of PC(18:0/Aze) … in recovery samples but not in the acute dataset" | The name **does not exist as a standalone identification**. It is the unselected half of `PC(20:0/8:0(COOH))_PC(18:0/Aze)`; the pipeline's `picked_candidate` column selects the former, and "Aze" has **zero** hits in `Recovery_DEA_results_v1_corrected.csv`. The feature is present in all 18 recovery samples **and in the blank**, with a max-to-blank ratio of 4.18. Its absence from the acute data is a **panel difference**: precursor 694 is not in the acute transition list, whose nearest members are 688, 690 and 692. See §5.8. |
| 1.9 | Results, Figure 2D text | RDH12 IC50 "95% CI: 213.3 to 226.9" | Lower bound is **213.4** (213.396). Truncation rather than rounding. |
| 1.10 | Results and Methods, viability | "a broad concentration range (0-400 µM atRAL)" | **There is no 0 µM measurement.** Vehicle arms used 15 concentrations from 25 to 400 µM, with 375 absent. Ferrostatin arms used **only 7**, from 50 to 400 µM. The 100% baseline is a fixed model parameter, not data. Source: `Cell_Viability/data/atRAL_cell_viability.csv`. |

---

## 2. Statistical framing that must be addressed (P2)

### 2.1 Lipidomics: unadjusted p-values presented as significance

The acute volcano and class-count scripts filter on the **raw** p-value column under a variable
named `padj_cutoff`:

```r
# Lipidomics/RScripts/Acute/Acute_Volcanos_stnd.R
padj_cutoff = 0.1,
lfc_cutoff  = log2(1.3),
...
padj_col  <- paste0(contrast, "_p.val")   # reads the RAW p column
```

The same pattern is in `Lipidomics/RScripts/Recovery_Count_SigUpDownperclass.R`. Figure axes are
labelled as p-value, so the panels are internally honest, but the Results text calls these marks
significant.

After Benjamini-Hochberg correction, using `Acute_DEA_results_v1_corrected.csv` with 247 features:

| Acute contrast | adjusted p ≤ 0.05 | adjusted p ≤ 0.10 | Figure rule: raw p ≤ 0.1 and 1.3-fold |
|---|---|---|---|
| 100 µM vs vehicle | **0** | **0** | 12 |
| 200 µM vs vehicle | 2 | 10 | 16 |
| 200 vs 100 µM | 1 | 1 | 29 |

**Action.** Choose one convention and apply it consistently.

- *Recommended.* Keep the figures, relabel them explicitly as unadjusted p, and restate the acute
  lipid results as trends. Name the two adjusted-significant species. See §5.4.
- *Alternative.* Re-run the volcano scripts against `_p.adj` and accept that the 100 µM panel will
  be empty.

Note that the results table's own `*_significant` column uses a **third** rule, adjusted p ≤ 0.05
and at least 1.5-fold, set in `Lipidomics/RScripts/Acute/4b_FixBHpadjust.R`. Pick one and state it
in Methods.

### 2.2 Lipidomics: the composition tests already exist and are null

`Lipidomics/Figures/Acute/DBE_composition/Acute_DBE_composition_fisher_tests_vs_vehicle.csv`
contains a per-class Fisher exact test against vehicle for all 26 class-by-dose cells.
**Every BH-adjusted p is exactly 1.0**, and the smallest unadjusted p in the file is 0.476.
Stating the double-bond compositional shifts without this caveat contradicts an analysis that
already lives in this repository. See §5.6.

### 2.3 Proteomics: detection depth differs across groups and biases findings in one direction

The acute vehicle group is the least deeply sampled group in the dataset.

| Group | Quantified in all 3 | Undetected in all 3 |
|---|---|---|
| Acute RDH12 vehicle | **5,984** | **893** |
| Acute RDH12 100 µM | 7,514 | 194 |
| All other groups | 6,358 to 7,219 | 195 to 418 |

From `Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv`:

- **360 of the 535 increases at 100 µM (67%)** and **451 of 664 at 200 µM (68%)** involve proteins
  undetected in every vehicle replicate. Their fold changes are imputation output.
- Only about 26 to 27 percent of reported increases rest on a fully quantified vehicle baseline.
- **Decreases are largely real.** Only 4 to 6 percent arise from non-detection in the treated group.

**Action.** Add the disclosure paragraph in §5.3, and add the `missclass_*` columns to the SI
protein tables so readers can filter.

### 2.4 The chaperonin increase is dose-limited

All eight CCT subunits rise at 100 µM, with log2 changes of 0.74 to 1.33 and adjusted p at or
below 0.026. At 200 µM **every one is flat**, with log2 changes of -0.03 to 0.55 and adjusted p at
or above 0.19. Tubulins follow the same pattern; the sole exception, TUBA4A, is undetected in
vehicle. The prefoldin decrease is solid at both doses on clean baselines. See §5.9.

---

## 3. Qualifications and softening (P3)

| # | Item | Finding | Suggested handling |
|---|---|---|---|
| 3.1 | POLD4 "log2FC = -9.87" | Undetected in all three replicates at 200 µM. Not a measured 900-fold loss. | State it as loss of detection, or drop the number. |
| 3.2 | MOCS3, PAN2, ZDHHC18 as "largest increase" | All three **undetected in all three vehicle replicates**. MOCS3 ranks eighth at 200 µM, not first; PAN2 and ZDHHC18 fall outside the top ten. | Reframe as proteins detected only after treatment. Ranking by imputed fold change is not meaningful. |
| 3.3 | PRKCE "increased abundance is notable" | Undetected in all three vehicle replicates; log2 changes of 4.25 and 4.51 are imputation output. | Remove, or flag explicitly. |
| 3.4 | PNPLA6 | Only significant acute contrast is 200 versus vehicle (adjusted p 0.016), on a one-of-three vehicle baseline. Both 100 versus vehicle and 200 versus 100 are null. Its acute vehicle value is the lowest across all nine groups, and untreated vehicle cells rise by 1.18 in log2 (adjusted p 0.0054) between acute and recovery with no atRAL present. | Present as a possibility, not an observation. Handled in §5.5. |
| 3.5 | Recovery ferroptosis GO enrichment | Real at adjusted p 8.0e-4, but driven by **exactly three genes out of 24**: HMOX1, NFE2L2, SLC39A7. NFE2L2 and SLC39A7 are **one-of-three versus zero-of-three** between genotypes. | State the gene count and the detection caveat. |
| 3.6 | Control-associated recovery proteins | CCN1, CCN2, FN1, WNT5A, ERBB2, FGFR2, SPRY2 and UROS are **all undetected in RDH12 cells**. ATF3, ABCB10, SLC11A2 and IBA57 are clean three-of-three. | Lead with the four well-measured proteins. |
| 3.7 | "four ceramide and dihydroceramide species" at 100 µM | The count is exact, but it is **three ceramides and one dihydroceramide**, with two C16 species and one C20; the C20 is the dihydroceramide. The fourth is a C24 hydroxy species not mentioned. None survives correction, lowest adjusted p 0.248, and two further ceramides with larger effects missed on p. | See §5.4. |
| 3.8 | PE depletion | Direction is strong, 24 of 25 species negative, but **no PE reaches adjusted p below 0.10**; the best is 0.067. | See §5.4. |
| 3.9 | LPC 16:0 | log2FC 0.499, adjusted p **0.120**. Not significant under any adjusted rule. LPC 18:1 at 0.0045 and PC 32:0 at 0.017 do hold. | See §5.4. |
| 3.10 | "increased number of DG species" | Detection counts run 10 in vehicle, 17 at 100 µM, 13 at 200 µM. **Non-monotonic**, and a species count rather than abundance. | Say "at 100 µM" and label it a detection count. |
| 3.11 | "rapid, dose-dependent collapse" for control cells | Fitted Hill slopes are -3.15 for control vehicle and **-6.49 for RDH12 vehicle**. The RDH12 curve is the steepest of the four. | Avoid implying control cells die more abruptly. |

---

## 4. Figure panels and the scripts that build them

Regenerate only what the corresponding text change requires. All paths are repo-relative.
Prepend `C:\Program Files\R\R-4.4.1\bin` to PATH, since `Rscript` is not on PATH. Running a script
from the repo root drops a stray `Rplots.pdf`; delete it afterwards.

| Manuscript panel | Script | Output |
|---|---|---|
| Fig 2D, Fig 6B, IC50 and ferrostatin | `Cell_Viability/R/01_atRAL_dose_response.R` | `Cell_Viability/figures/atRAL_dose_response_curves.{pdf,png}`, `atRAL_IC50_by_group.{pdf,png}`, `atRAL_IC50_estimates.csv`, `atRAL_IC50_anova.txt` |
| Fig 3A, 3B, acute proteome volcanoes | `Proteomic_Rscripts/Create_Volcanos_allContrasts_toDirctry.R` | `Proteomic_Figs/VolcanoPlots/StandardVolcano2/` |
| Fig 3A, 3B alternative with missingness shapes | `Proteomic_Rscripts/Create_Volcanos_wMissingess_allContrasts_toDirctry.R` | `Proteomic_Figs/VolcanoPlots/wMissingessInfo/` |
| **Fig 3E, acute pathway heatmap** | `Proteomic_Rscripts/Figures/HeatMap_Fig1_horizontal_maintext.R`; thresholds at lines 27–28, gene lists at lines 64–97 | `Proteomic_Figs/enriched/Fig1_horizontal_maintext_heatmap_centered_sigOnly.{pdf,svg}` plus a manual copy in `Proteomic_Figs/enriched/1_Maintext_version/` |
| Fig 4A, lipid class counts | `Lipidomics/RScripts/Acute/Acute_DistinctPrecursors_LipidClass_Figure.R` | `Lipidomics/Figures/Acute/` |
| **Fig 4B, 4C, acute lipid volcanoes** | `Lipidomics/RScripts/Acute/Acute_Volcanos_stnd.R`; the raw-p issue is at the `padj_col` assignment | `Lipidomics/Figures/Acute/VolcanoPlots/StandardVolcano/` |
| Fig 4D, 4E, class direction counts | `Lipidomics/RScripts/Acute/Acute_Count_SigUpDownperclass.R` | `Lipidomics/Figures/Acute/Lipid_class_direction_counts/` |
| Fig 4F, acute targeted lipid heatmap | `Lipidomics/RScripts/Acute/HeatMap_KeyLipids2.R`; lipid list at line 61, patterns at line 81 | `Lipidomics/Figures/Acute/Acute_targeted_heatmap/` |
| **Fig 4G, 4H, double-bond composition** | `Lipidomics/RScripts/Acute/Acute_DBE_Composition_Figure.R` | `Lipidomics/Figures/Acute/DBE_composition/`, which also holds the null Fisher test table |
| Fig 5A–C, recovery volcanoes | Same proteome volcano scripts as Figure 3 | `Proteomic_Figs/VolcanoPlots/StandardVolcano2/` |
| Fig 5E, recovery pathway heatmap | `Proteomic_Rscripts/Figures/HeatMap_Fig4wFerrop.R`; thresholds at lines 41–42 | `Proteomic_Figs/enriched/Fig4a_pathway_heatmap_centered_noSigAnno_WTfarRight.{pdf,svg}` plus copy in `1_Maintext_version/` |
| Fig 5G, recovery class direction | `Lipidomics/RScripts/Recovery_Count_SigUpDownperclass.R`; same raw-p issue | `Lipidomics/Figures/Recovery/Lipid_class_direction_counts/` |
| Fig 5H, recovery lipid heatmap | `Lipidomics/RScripts/Rvr_HeatMap.R` | `Lipidomics/Figures/Recovery/Recovery_targeted_heatmap/` |
| Fig 6A, ferroptosis heatmap | `Proteomic_Rscripts/Figures/HeatMap_Ferroptosis_maintext_Fig.R` | `Proteomic_Figs/ferroptosis_heatmap/` |
| Fig S5, S7, GO dotplots | Scripts under `GOanalysis/` | `GOanalysis/output/clusterProfiler_batch/<contrast>_{up,down}/` |

**Note on the main-text folder.** `Proteomic_Figs/enriched/1_Maintext_version/` holds hand-copied
duplicates. The scripts write to the parent `enriched/` directory, so a re-run does **not** update
the main-text copies. Either change `out_dir` in the two heatmap scripts or copy the files after
each run.

---

## 5. Revised passages

Drop-in replacements. µM is used throughout; the draft currently mixes "µM" and "uM" in the
Discussion.

### 5.1 Abstract, antioxidant sentence

> Acute atRAL exposure engaged the supply arm of the KEAP1–NRF2 program, with increased HMOX1,
> SLC7A11, GCLC, MGST2 and MGST3 at 100 µM, while the terminal glutathione peroxidase arm
> declined; NRF2 itself rose only during recovery. Acute exposure was accompanied by membrane
> remodeling, indicated by increased glycerophospholipid enzyme abundance and by directional
> increases in lysophosphatidylcholine and ceramide species.

### 5.2 Results, ferrostatin-1 experiment

Replaces the closing passage of the ferroptosis section.

> To test functionally whether ferroptosis contributes to atRAL-induced death, we measured the
> atRAL IC50 in WT and RDH12-expressing cells after pre-treatment with 30 µM ferrostatin-1.
> Ferrostatin-1 did not protect WT cells: their IC50 was 133.7 µM (95% CI 127.0–140.4) without
> and 123.5 µM (95% CI 113.6–133.5) with the inhibitor, a 0.92-fold change indistinguishable from
> no effect (Tukey p = 0.53). In RDH12-expressing cells the same treatment raised the IC50 from
> 220.2 µM (95% CI 213.4–226.9) to 351.5 µM (95% CI 330.8–372.2), a 1.60-fold increase
> (Tukey p = 1.7 × 10⁻⁴). The genotype-by-ferrostatin interaction was significant
> (F(1,8) = 46.3, p = 1.4 × 10⁻⁴), and the separation between genotypes widened from 1.65-fold to
> approximately 2.85-fold (Figure 6B). Radical trapping therefore protects only in the RDH12
> background, indicating that lipid peroxyl radical propagation becomes a rate-limiting
> contributor to death only once aldehyde clearance is intact.

This also resolves item 1.2. The Discussion paragraph on ferrostatin already argues this reading
and needs only its two IC50 values and the fold change updated.

### 5.3 Results or Methods, detection depth disclosure

New paragraph.

> Detection depth differed across sample groups, and the acute vehicle group was the least deeply
> sampled: 5,984 protein groups were quantified in all three replicates, compared with 6,358 to
> 7,514 in the other eight groups, and 893 proteins were undetected in all three acute vehicle
> replicates, compared with 194 to 418 elsewhere. As a consequence, 360 of the 535 increases at
> 100 µM and 451 of the 664 increases at 200 µM involve proteins not detected in any vehicle
> replicate, so their reported fold changes are determined by the imputation model rather than by
> measurement. Decreases are far less affected, with only 4 to 6 percent arising from
> non-detection in the treated group. Fold-change magnitudes for proteins with a missing
> reference condition should therefore be read as directional. The missingness class of every
> protein in every group is reported in SI Table X.

### 5.4 Results, acute lipidome

Replaces the opening of that section.

> In the acute RDH12-expressing lipidome, changes were subtle. The acute panel comprised 247
> features and treatment effects were small, so few species survived Benjamini-Hochberg
> correction: none at 100 µM, and two at 200 µM, namely LPC 18:1 (log2FC 1.09, adjusted p 0.0045)
> and PC 32:0 (log2FC 0.86, adjusted p 0.017), with a further eight species at an adjusted p below
> 0.1. The volcano plots in Figure 4B, C are plotted on unadjusted p-values and should be read as
> showing trends rather than adjusted significance. LPC 16:0 rose in the same direction without
> reaching significance after adjustment (log2FC 0.50, adjusted p 0.12). Phosphatidylethanolamine
> species shifted consistently downward at 200 µM, with 24 of 25 detected species showing a
> negative fold change, although none reached significance after adjustment (lowest adjusted p
> 0.067). This consistent direction is notable because atRAL condenses with the primary amine of
> phosphatidylethanolamine to form N-retinylidene-PE.
>
> At 100 µM atRAL, four sphingolipid species increased relative to vehicle, comprising three
> ceramides and one dihydroceramide and including two C16 species and one C20 species
> (Figure 4B, F). None reached significance after correction (lowest adjusted p 0.25), and all
> four returned to vehicle levels at 200 µM, so we report this as a dose-specific trend rather
> than an established change.

### 5.5 Discussion, glycerophospholipid enzymes

Replaces the whole paragraph.

> The proteome showed a coincident increase in enzymes of glycerophospholipid synthesis and
> remodeling, though the strength of evidence differs among them because the acute vehicle group
> was the least deeply sampled group in the dataset. Several observations rest on proteins
> quantified in every replicate of both groups compared. LPCAT1, a reacylase favoring saturated
> and monounsaturated acyl-CoA donors, increased at both doses; the cytosolic calcium-dependent
> phospholipase A2 PLA2G4A increased at 200 µM, as did the lysosomal phospholipase A2 PLA2G15,
> indicating that the deacylation arm of the Lands cycle was engaged rather than silent; and
> PTDSS1, which supports phosphatidylserine synthesis by base exchange, increased at both doses
> from a vehicle baseline quantified in two of three replicates. A second and larger set of
> enzymes also scored as increased, comprising GPAT4 and AGPAT3/4/5 for phosphatidic acid, CHKA
> and ETNK1 in the Kennedy pathway, and DGKE linking diacylglycerol back to phosphatidic acid.
> Each of these was undetected in all three vehicle replicates, so the magnitude of the reported
> change is set by imputation and should be read as directional rather than quantitative. Two
> enzymes central to the interpretation could not be evaluated at all: LPCAT3, which incorporates
> arachidonate and promotes ferroptosis through the ACSL4 axis, and the calcium-independent
> phospholipase A2 PLA2G6 were not detected in any sample, so absence of change cannot be
> distinguished from absence of measurement. ACSL4 itself was quantified throughout and changed
> by less than 1.3-fold at both doses. Read together, the coincidence of elevated LPC with
> elevated LPCAT1 is consistent with reacylation skewed toward less oxidizable acyl chains, and
> the concurrent rise in PC 32:0 would be compatible with LPCAT1-mediated reacylation of LPC 16:0
> with palmitoyl-CoA, though our MRM-based annotation does not allow us to demonstrate that route
> definitively. PNPLA6 was also scored as increased at 200 µM, but on a vehicle baseline
> quantified in only one of three replicates, and neither the 100 µM comparison nor the 200 versus
> 100 µM comparison reached significance. In retinal pigment epithelial cells PNPLA6 acts as a
> phospholipase B that deacylates PC and then hydrolyzes the resulting LPC to
> glycerophosphocholine, which re-enters PC synthesis through the Kennedy pathway; if the increase
> is genuine it would therefore favor enhanced engagement of this PC regeneration loop over net
> LPC production, but we do not rest any conclusion on it. These data are consistent with an
> attempted membrane-remodeling response, but they do not establish that repair outpaces damage.

### 5.6 Results, double-bond composition

Replaces the closing sentences of that paragraph.

> From this analysis, we observed that DGs had a higher proportion of PUFA-containing precursors
> at 100 µM atRAL and TGs at 200 µM, that ceramide precursors shifted toward saturated species,
> and that PC composition remained comparatively stable relative to vehicle (Figure 4G, H). These
> are descriptive proportions of distinct detected precursors and are not abundance-weighted.
> Per-class Fisher exact tests against vehicle were not significant for any lipid class at either
> dose (all BH-adjusted p = 1.0; smallest unadjusted p = 0.48), and several classes rest on very
> few precursors, with the ceramide shift resting on zero, two and one saturated species across
> the three conditions. We therefore present these compositional observations as descriptive and
> hypothesis-generating rather than as statistically supported findings.

### 5.7 Results, control-enriched recovery processes

> In contrast, proteins more abundant in WT cells than in RDH12-expressing cells were enriched for
> regulation of the ERK1 and ERK2 cascade, fibroblast growth factor receptor signaling, and
> morphogenesis-related processes, with heme biosynthesis approaching significance (adjusted p =
> 0.052) (Figure S7D). Individual matricellular and extracellular matrix proteins, including CCN1,
> CCN2 and FN1, were also more abundant in WT cells, although extracellular matrix organization
> and cell adhesion were not themselves enriched terms. Among these WT-associated proteins, ATF3,
> ABCB10, SLC11A2 and IBA57 were quantified in all replicates of both genotypes, whereas CCN1,
> CCN2, FN1, WNT5A, ERBB2, FGFR2, SPRY2 and UROS were undetected in RDH12-expressing cells, so
> their fold changes reflect loss of detection.

### 5.8 Results, oxidized phosphatidylcholine species

> We also searched the recovery dataset for putatively oxidized PC species using the LIPID MAPS
> database. A feature at MRM transition m/z 694.5→184 received an ambiguous annotation between
> PC(20:0/8:0(COOH)) and PC(18:0/Aze); our annotation pipeline selected the former, and we refer
> to it under that name. It was more abundant in RDH12-expressing cells than in WT cells after
> recovery from 200 µM atRAL (log2FC 0.79, adjusted p 0.0098) and more abundant than in
> vehicle-treated RDH12 cells (log2FC 0.96, adjusted p 0.014). This transition is not present in
> the acute panel, so the two timepoints cannot be compared for this species.
> PC(16:0/8:0(COOH)), at m/z 652.5→184, was detected in both datasets and was not significantly
> altered between conditions (Figure 5H). These observations show that putatively oxidized PC
> species remain detectable after atRAL removal.

### 5.9 Results and Discussion, chaperonin and prefoldin

Results:

> Another feature of the heatmap was the increased abundance of multiple subunits of the
> chaperonin-containing TCP-1 complex (CCT/TRiC), along with several tubulin proteins, at 100 µM
> atRAL (Figure 3E). This increase was not sustained at 200 µM, where CCT subunit abundance
> returned to vehicle levels. In contrast, several subunits of the prefoldin complex were
> decreased at both doses.

Discussion, replacing the final sentence of the acute-proteome paragraph:

> The divergence within the cytoskeletal folding pathway fits the same logic, with the
> qualification that it is dose-limited: CCT/TRiC subunits and tubulins rise at 100 µM and return
> to baseline at 200 µM, whereas prefoldin subunits fall at both doses. Folding capacity for
> existing clients is transiently reinforced while delivery of newly synthesized actin and
> tubulin, which depends on prefoldin, is curtailed along with translation.

### 5.10 Figure 3E legend

Replace "that met a significance threshold of fold change ±1.5 and BH-adjusted p < 0.01" with:

> that met a significance threshold of fold change ±1.4 and BH-adjusted p < 0.1

### 5.11 Methods, viability assay

> Cells were exposed to atRAL across 15 concentrations from 25 to 400 µM in the vehicle arms and 7
> concentrations from 50 to 400 µM in the ferrostatin-1 arms, with n = 3 biological replicates per
> group. Dose-response curves were fitted with a four-parameter log-logistic model (`drc::LL.4`)
> with the upper and lower asymptotes fixed at 100% and 0%. No zero-dose viability measurement was
> collected, so the 100% baseline is a model constraint rather than an observation. IC50 values are
> pooled fits across replicates with delta-method 95% confidence intervals; per-replicate IC50
> values were compared by two-factor ANOVA on log10(IC50) with Tukey HSD.

---

## 6. Repository housekeeping (P4)

| # | Item |
|---|---|
| 6.1 | The 5 h atRAL exposure and 24 h ferrostatin pre-treatment durations are **not recorded anywhere** in `Cell_Viability/`. Add them to the data file or a README so the Methods statement is traceable. |
| 6.2 | Rename `padj_cutoff` and `padj_col` in `Acute_Volcanos_stnd.R` and `Recovery_Count_SigUpDownperclass.R` to reflect that they read `_p.val`, or switch them to `_p.adj`. The current naming will mislead the next reader. |
| 6.3 | Add the `missclass_*` columns to `SI_Table_Protein_DEA.csv` so reviewers can filter imputation-driven fold changes. |
| 6.4 | `Lipidomics/RScripts/Acute/HeatMap_KeyLipids.R` and `HeatMap_KeyLipids2.R` write the **same output filenames**. Delete or rename the superseded first version. |
| 6.5 | One lipid name has an unclosed parenthesis, `Cer(d18:1/20:0(OH)`, in the `name` column. Harmless to the analysis but worth fixing before the SI tables ship. |
| 6.6 | Consider pointing the two main-text heatmap scripts directly at `Proteomic_Figs/enriched/1_Maintext_version/` so re-runs update the figures actually used. |
| 6.7 | The acute vehicle group's low detection depth deserves a sentence in the limitations paragraph, independent of the disclosure in §5.3. |

---

## 7. Claims that checked out, no action needed

- 8,049 protein groups quantified.
- 535 and 1,850 at 100 µM, 664 and 823 at 200 µM. The text says "less than" and "greater than"
  where the code used inclusive bounds; strictly these are adjusted p ≤ 0.01 and at least 2-fold.
- 28 increased and 49 decreased after 200 µM recovery.
- HIF1AN log2 changes of -1.4 and -1.3.
- Retinol nearly tripled, 8.90 over 3.00 being 2.97, and retinoic acid suppressed by over 87%,
  actually 87.6%.
- WT IC50 133.7 µM (95% CI 127.0–140.4) and RDH12 IC50 220.2 µM; the ratio is 1.65-fold.
- Over 99% of protein groups shared at baseline, 8,328 of 8,353 being 99.7%, Pearson r = 0.99,
  and RDH12 detected only in RDH12-expressing cells.
- No induction of canonical unfolded protein response chaperones. In fact HSPA5, ATF4, ATF6,
  DNAJB9, CALR and CANX all decrease, which is a stronger result than the claim made.
- Prefoldin decrease: all six subunits down, five at adjusted p below 0.05, on clean baselines.
- N-glycan machinery: STT3A, STT3B, RPN1, RPN2 and UGGT1 are solid. ALG1, ALG11 and DPAGT1 are
  imputation-driven.
- Acute GO enrichment for transmembrane transport, cellular lipid metabolic process and protein
  N-linked glycosylation.
- Recovery microtubule enrichment in RDH12 cells, with adjusted p down to 2.9e-14.
- **Recovery genotype lipid comparison**, the strongest result in the paper and currently
  undersold: ether-PC mean log2FC +0.469 with 24 of 25 adjusted-significant hits going up; DG with
  48 of 49 species down and 24 adjusted-significant; PC at +0.195. This survives proper adjustment
  at 0.05 with no raw-p rule needed.
