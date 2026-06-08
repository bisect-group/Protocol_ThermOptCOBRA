# Protocol for Thermodynamically Optimal Construction and Analysis of Metabolic Networks using ThermOptCOBRA

This directory contains the implementations and tools supporting the protocol to construct and analyze thermodynamically optimal genome-scale metabolic networks.

## Directory Structure and Descriptions

All software packages and algorithms are located inside the `Protocol_ThermOptCOBRA` directory:

* **[`ThermOptCOBRA`](ThermOptCOBRA)**:
    The core suite of MATLAB-based functions for incorporating thermodynamic constraints into metabolic network analysis. It contains:
  * `ThermOptEnumerator`: For enumerating all thermodynamically infeasible cycles (TICs) in a model.
  * `ThermOptCC`: For consistency checking and identifying thermodynamically blocked reactions.
  * `ThermOptiCS`: For constructing context-specific models (CSMs) under thermodynamic constraints.
  * `ThermOptFlux`: For post-processing flux distributions to remove cycle-related fluxes.
  * *Associated Paper:* Kumar & Bhatt (2025), *iScience*.

* **[`Localgini`](Localgini)**:
    Contains the `GiniReactionImportance.m` implementation. This tool quantifies gene expression variability across samples using Gini coefficients to establish thresholding and extract core reactions for context-specific model reconstruction.
  * *Associated Paper:* Kumar & Bhatt (2025), *npj Systems Biology and Applications*.

* **[`spectraCC`](spectraCC)**:
    Contains the flux consistency checking modules (`spectraCC.m`, `forwardcc.m`, `reverse.m`) belonging to the SPECTRA framework. SPECTRA is a generalist method to reconstruct metabolic networks from multi-omics data at a large scale.
  * *Associated Paper:* Kumar et al. (2026), *bioRxiv*.

* **[`looplessFluxSampler`](looplessFluxSampler)**:
    An efficient toolbox utilizing the Adaptive Direction Sampling on a Box (ADSB) algorithm to sample the non-convex loopless and mass-balanced flux solution space of metabolic models.
  * *Associated Paper:* Saa et al. (2024), *BMC Bioinformatics*.

---

## Citations and References

If you use these algorithms or resources in your research, please cite the respective papers:

1. **ThermOptCOBRA**
    * Kumar, S. P., & Bhatt, N. P. (2025). ThermOptCobra: Thermodynamically optimal construction and analysis of metabolic networks for reliable phenotype predictions. *iScience*, 28(8), 113005.
    * [Journal Link / DOI](https://www.cell.com/iscience/fulltext/S2589-0042(25)01266-0)

2. **Localgini**
    * Kumar, S. P., & Bhatt, N. P. (2025). Modelling reliable metabolic phenotypes by analysing the context-specific transcriptomics data. *npj Systems Biology and Applications*, 11(1), 23.
    * [Journal Link / DOI](https://www.nature.com/articles/s41540-025-00617-8)

3. **SPECTRA (spectraCC)**
    * Kumar, S. P., Sridhar, S., Alsmadi, N., Mahadevan, R., & Bhatt, N. P. (2026). Generalist method to reconstruct metabolic networks from multi-omics data at large-scale. *bioRxiv* preprint.
    * [Preprint Link / DOI](https://www.biorxiv.org/content/10.64898/2026.04.02.716249v1.full)

4. **LooplessFluxSampler**
    * Saa, P. A., Zapararte, S., Drovandi, C. C., & Nielsen, L. K. (2024). LooplessFluxSampler: An efficient algorithm for sampling the loopless flux solution space of metabolic models. *BMC Bioinformatics*, 25, 12.
    * [Journal Link / DOI](https://link.springer.com/article/10.1186/s12859-023-05616-2)
