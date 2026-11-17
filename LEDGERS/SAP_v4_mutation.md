# SAP v4 — Study 4: seeded faults in a Voronoi implementation (confirmatory)

Manuscript 261114SW (H. Tajima). Written 2026-11-17 before any mutant was generated or run and before any Study 4 input was generated. This plan is published in a public GitHub repository so that its time of publication is recorded by a third party; Studies 1–3 were hashed only locally.

## 1. Question
Studies 1–3 measured what oracle-free tests miss on faults that happened to exist in five libraries, two of them used outside their primary domain. Study 4 asks the same question on controlled, in-domain faults: of the seeded faults that an exact reference detects without any exception, how many do the oracle-free tests of the paper detect?

## 2. Subject
`study4/code/hp_voronoi.py`: a 40-line Voronoi implementation by half-plane clipping in double precision (pure Python), written for this study. First-order mutants are generated deterministically by `study4/code/mutate.py` with the operators AOR, ROR, CRP, LCR, NEG and SDL (defined in the file's docstring); mutants identical in source to the original or to an earlier mutant are removed.

## 3. Inputs
Candidate base inputs: the generators of Study 1 (`code/vor_common.py`) for classes G1, G2, G3, G4, G5 and G9, with n ∈ {12, 40} and new seeds 300–303 (48 base inputs), each with five follow-ups: translation (MR-T), 90° rotation (MR-R), scaling by 2^10 (MR-S), permutation (MR-P) and reflection (MR-F). G6 (exact duplicates) is excluded because the subject's duplicate convention differs from the reference's; G7 and G8 are excluded because the unmutated subject exceeds τ on them (pilot, Section 8). A base input is **valid** if, for the unmutated subject, all six runs complete, pass the exact reference at τ = 1e-9, raise no definition-check flag, differ from Shapely by at most τ, and all five metamorphic distances are at most τ. Only valid base inputs are used for mutants. This rule is applied after this plan is published.

## 4. Oracles and kills (per mutant, over all runs of all valid base inputs)
- **Non-silent mutant:** any run raises an exception, exceeds 5 s, or returns a list of the wrong length.
- **Reference kill:** a completed run has a missing cell, a misassigned cell or a normalised symmetric difference above τ = 1e-9 against the exact rational reference (`vor_common.reference_cells`, `compare_to_reference`).
- **Silent reference-killed mutant (SRK):** not non-silent and reference-killed. All recalls below use the SRK mutants as denominator (unit: mutant).
- **Definition checks (DEF):** P1–P3 of Study 1 (completeness, area partition at τ, nearest-site test in exact arithmetic) flag at least one completed run.
- **Metamorphic relations:** the distance between the mutant's base cells and its follow-up cells mapped back exceeds τ for at least one family; exact MRs are MR-R, MR-S, MR-P, MR-F; all MRs add MR-T.
- **Differential test (DIFF):** at least one completed run differs from Shapely's cells (ordered mode) by more than τ.

## 5. Predictions and decision rules (each at α = 0.01; Bonferroni over five)
- **P1** Definition checks miss some SRK mutants: supported if the one-sided upper 99% Clopper–Pearson limit of DEF recall is below 0.95.
- **P2** Definition checks detect more SRK mutants than the exact MRs: supported if the exact one-sided McNemar test (discordant pairs, binomial) gives p < 0.01.
- **P3** The union of all oracle-free tests of the paper (DEF ∪ all MRs) misses at least 5% of SRK mutants: supported if the one-sided lower 99% Clopper–Pearson limit of the miss rate is at least 0.05.
- **P4** A trusted second implementation is nearly an oracle: supported if the one-sided lower 99% Clopper–Pearson limit of DIFF recall is at least 0.90.
- **P5** What the oracle-free tests miss is small: among SRK mutants, the worst reference error of those missed by DEF ∪ all MRs is smaller than of those detected (one-sided Mann–Whitney, p < 0.01; not supported if fewer than three are missed).

Rationale: Studies 1–3 found that definition checks catch missing and misassigned cells but not small geometric errors, that exact metamorphic relations are blind to faults that commute with the transformation, and that Shapely agreed with the exact reference to 1e-12 on inputs without repeats.

## 6. Exploratory (not confirmatory)
Per-operator and per-MR breakdowns; counts of non-silent and not-killed mutants; any analysis not listed in Section 5.

## 7. Code and procedure (SHA-256 at publication)
- `study4/code/hp_voronoi.py`: `5bf98fb9f51ad1c8de64c73133fe7fe213b0f7c69fac3c2b8f64724c43f096d3`
- `study4/code/mutate.py`: `28ada294410f1c41b20170a4590ba0446d8d40bff69de1753e559790576190d6`
- `study4/code/s4_common.py`: `2e7dd35f68baacfd621a9590a85f6100f9618ff96738bbc0ad803354b3960e27`
- `study4/code/run_study4.py`: `848dac7ccf4a0879e0eab7f79b5885754478c35ac951495034a2341e6d353f23`
- `study4/code/analyze_study4.py`: `65dbadaef8ef71c41cb01bacd5680acc8a15700538bdaf9836c72c2c0f7689a9`
- `study4/code/pilot_unmutated.py`: `4aa92b82ac55213944057d3433a6e68fe191d68ac0ec015996b75d59b5f94c42`
- `code/vor_common.py`: `dcf4572cdbc80a015342123e9248ae6d567af19bbe99f46349710248b84e04fb`
- `code/vor_jobs.py`: `dbc647cd8c8e6ec154c13142ee0eba0672eeea663119893291685518a56b4542`

Run: `python3 study4/code/run_study4.py` then `python3 study4/code/analyze_study4.py` (outputs in `study4/results/`). Any deviation is logged in the manuscript's deviation log before the affected analysis is run.

## 8. What was done before publication
The unmutated subject was run on Study 1 base inputs with seeds 0–1 (`pilot_unmutated.py`): it passed at τ on all classes except G7 (0/4) and G8 (0/4), which is why they are excluded. The harness was tested once on the Study 1 input G1, n = 12, seed 0, with the unmutated subject and with one hand-written change (`fu <= 0` → `fu < 0`); those outputs were discarded. The mutant generator was tested on an unrelated five-line toy function (not on `hp_voronoi.py`), which exposed and fixed a defect in locating nested operators. No mutant of `hp_voronoi.py` was generated, and no input with seeds 300–303 was generated.
