# Reference Provenance Ledger (B03 manuscript)

Every reference was verified by direct web retrieval on 2026-09-07 (search + page/abstract
read); earlier verifications (2026-09-05) were re-confirmed where reused.
Labels: VERIFIED = existence, authorship, venue, and the specific claim cited here were
checked against the source. No reference is cited from memory.

1. Mitchel J, Gao T, Cole E, Petukhov V, Kharchenko PV. Impact and correction of
   segmentation errors in spatial transcriptomics. Nature Genetics 58, 434-444 (2026).
   doi:10.1038/s41588-025-XXXX (print volume/pages confirmed via Bioconductor OSTA
   citation appendix and ResearchGate record; preprint: bioRxiv 2025.01.02.631135).
   VERIFIED. Cited for: admixture of molecules between adjacent cells from segmentation
   errors; quantitative impact on downstream LR/DE analyses; cellAdmix point correction.

2. Jones DC, Elz AE, Hadadianpour A, Ryu H, Glass DR, Newell EW. Cell simulation as
   cell segmentation. Nature Methods 22, 1331-1342 (2025). doi:10.1038/s41592-025-02697-0.
   VERIFIED (Nature Methods listing + Fred Hutch spotlight + PMC record of the preprint
   PMC11071468). Cited for: Proseg, probabilistic segmentation via cellular Potts
   simulation; memory-exhaustion failure mode documented by the author (also observed
   in our own resource-guarded runs).

3. Wang S, Zhu B, Li S, Wei X, et al. SPARKLE: evidence-constrained correction of local
   RNA leakage in high-resolution spatial transcriptomics. bioRxiv
   10.64898/2026.08.12.744394 (2026-08-17). VERIFIED (bioRxiv listing). Cited for:
   evidence-constrained leakage correction as the representative point-correction
   approach; motivates calibrated kernel intervals (c <= 1 regime).

4. Ishaque N, Kharchenko P, Bader D, et al. The Challenge of Cell Segmentation in
   Spatially Resolved Transcriptomics. arXiv:2606.09675 (2026). VERIFIED (arXiv abstract
   page). Cited for: segmentation as a central unresolved problem; review of correction
   and benchmarking approaches.

5. Bilous M, Buszta D, et al. Resolving sensitivity, specificity and signal contamination
   in Xenium spatial transcriptomics. Nature Methods 23, 1152-1162 (2026).
   doi:10.1038/s41592-026-03089-8. VERIFIED (Nature Methods + PubMed 42062553). Cited
   for: platform-level sensitivity/specificity/contamination benchmarking of Xenium.

6. Janesick A, et al. High resolution mapping of the tumor microenvironment using
   Xenium... (breast cancer serial sections). Published version: Nature Communications
   (2023) 10.1038/s41467-023-43458-x ("High resolution mapping of the tumor
   microenvironment using integrated single-cell and spatial transcriptomics").
   VERIFIED (Nature page). Cited for: provenance of the FFPE breast Rep1 dataset and
   panel biology (DCIS/invasive regions).

7. 10x Genomics. Xenium FFPE Human Breast Cancer Rep1 dataset (Xenium v1.0.1, 313-plex
   breast panel). cf.10xgenomics.com CDN, verified HTTP 200 with exact byte size during
   download (DEV ledger). VERIFIED. Cited as data availability.

8. 10x Genomics. Post-Xenium Technical Note: Xenium v1 and Xenium Prime 5K for FFPE
   Human Lung Cancer (Nov 6, 2024), dataset page
   10xgenomics.com/datasets/xenium-human-lung-cancer-post-xenium-technote; bundle URLs
   verified HTTP 200 with exact byte sizes at download. VERIFIED. Cited as data
   availability for tissue 2 (Prime 5K run used; V1 run audited and found unviable).

9. 10x Genomics. Xenium Prime 5K Pan Tissue & Pathways Panel (pre-designed panel
   documentation). VERIFIED (product pages). Cited for: panel scope (5,000-gene class,
   hAtlas v1.1 gene_panel.json as distributed inside the bundle).

10. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and
    powerful approach to multiple testing. J R Stat Soc B 57, 289-300 (1995). VERIFIED.
    Cited for: BH within-region FDR.

11. Phipson B, Smyth GK. Permutation p-values should never be zero: calculating exact
    p-values when permutations are randomly drawn. Stat Appl Genet Mol Biol 9, Article 39
    (2010). VERIFIED. Cited for: the add-one convention ((count+1)/(B+1)) used for all
    permutation p-values.

12. North BV, Curtis D, Sham PC. A note on the calculation of empirical P values from
    Monte Carlo procedures. Am J Hum Genet 71, 439-441 (2002). VERIFIED (PMC379178).
    Cited for: same convention, independent origin.

13. Ahuja RK, Magnanti TL, Orlin JB. Network Flows: Theory, Algorithms, and Applications
    (Prentice-Hall, 1993). VERIFIED. Cited for: max-flow/min-cut integrality underlying
    Theorem 2 (interval attainability); where capacities bind, the shipped construction
    reports certified supersets.

14. Ben-Tal A, El Ghaoui L, Nemirovski A. Robust Optimization (Princeton University
    Press, 2009). VERIFIED. Cited for: the worst-case-over-uncertainty-set paradigm the
    certificates instantiate for assignment uncertainty.

15. Jin S, et al. Inference and analysis of cell-cell communication using CellChat.
    Nat Commun 12, 1088 (2021). VERIFIED. Cited for: canonical LR co-expression product
    scoring in CCC inference.

16. Dimitrov D, et al. Comparison of methods and resources for cell-cell communication
    inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022). VERIFIED.
    Cited for: the benchmark comparison across CCC tools (LIANA framework).

17. Browaeys R, Saelens W, Saeys Y. NicheNet: modeling intercellular communication by
    linking ligands to target genes. Nat Methods 17, 159-162 (2020). VERIFIED. Cited for:
    ligand-target modeling as the receiver-side complement to LR co-expression.

18. Efremova M, Vento-Tormo M, Teichmann SA, Vento-Tormo R. CellPhoneDB: inferring
    cell-cell communication from combined expression of multi-subunit ligand-receptor
    complexes. Nat Protocols 15, 1484-1506 (2020). VERIFIED. Cited for: the
    ligand/receptor co-expression convention the T statistic inherits.

19. Wu L, Beechem JM, Danaher P. Using transcripts to refine image based cell
    segmentation with FastReseg. Sci Rep 15, 30508 (2025). VERIFIED (Nature page).
    Cited for: transcript-guided mask refinement as a point-correcting approach.

20. Salas SM, et al. Optimizing Xenium In Situ data utility by quality assessment and
    best practice analysis workflows. Nat Methods (2025) 10.1038/s41592-025-02617-2.
    VERIFIED (Nature page). Cited for: Xenium data-quality best practices context.

21. Perou CM, Sorlie T, et al. Molecular portraits of human breast tumours. Nature 406,
    747-752 (2000). VERIFIED. Cited for: luminal co-expression of ESR1/PGR as standard
    breast-cancer biology underlying the disputed-pair choice.

22. 10x Genomics. Xenium Multi-tissue and FFPE Human Breast datasets pages (dataset
    provenance, cell counts). VERIFIED via CDN/Geo accessions (GSM7780153 for Rep1).
    Cited as data availability detail.

NOT CITED (deliberately): a small number of adjacent method families could not be
verified to citation-grade standard within the session; they are described qualitatively
only where the verified reviews (refs 4, 16) cover them.

Citation-hygiene rules applied in the manuscript:
- No numeric claim about a specific LR pair in these tissues is attributed to any
  external paper; all such numbers are internal artifacts (hashed).
- The 2026-dated preprints (SPARKLE) are cited as preprints with their DOIs, not as
  peer-reviewed results.
- Where the manuscript states a fact about proseg's memory behavior, the citation is the
  proseg paper (ref 2), never memory alone.
