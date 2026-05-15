# An Efficient Adjoint-Free Melnikov–Lyapunov Diagnostic for Transition Detection in Slow–Fast Chemical Reaction Networks
DOI: https://doi.org/10.5281/zenodo.20196726

## Authors
- **Sri Venkata Durga Sudarsan Madhyannapu**
- **Pradheep Kumar S.**

## Affiliations
1. Freshmen Engineering Department, NRI Institute of Technology (Autonomous), Pothavarappadu, Agiripalli, Vijayawada, 521212, Andhra Pradesh, India
2. Research Scholar, Jawaharlal Nehru Technological University Kakinada, Andhra Pradesh, India
3. School of Basic Sciences, SRM University AP, Neerukonda, Mangalagiri Mandal, Guntur, 522240, Andhra Pradesh, India

## Manuscript Information
| | |
|---|---|
| **Journal** | Communications in Nonlinear Science and Numerical Simulation (Elsevier), ISSN: 1007-5704 |
| **Manuscript ID** | CNSNS-D-26-00848 |
| **Status** | Under Review, 2026 |

## Overview
MATLAB implementations of the Adjoint-Free Melnikov–Lyapunov diagnostic for detecting bifurcation transitions in slow–fast chemical reaction networks (CRNs). Validated on 7 biochemical models (2D to 10D), achieving 235–612× speedup vs. MATCONT with average error 5.7%.

## Performance Summary

| Model | Dim | ε | Speedup | Error |
|---|---|---|---|---|
| Substrate Inhibition Oscillator | 2 | 0.005 | 235× | 7.2% |
| Enzymatic Feedback Network | 3 | 0.008 | 303× | 5.1% |
| Extended Glycolysis | 7 | 0.0008 | 446× | 6.3% |
| MAPK Cascade | 8 | 0.0005 | 468× | 4.8% |
| Mammalian Circadian Rhythm | 10 | 0.0002 | 612× | 5.9% |
| Modified Oregonator | 3 | 0.0004 | 341× | 8.1% |
| Brusselator | 2 | 0.001 | — | 2.8% |
| **Average** | | | **401×** | **5.7%** |

## Repository Structure
```
├── README.md
├── CITATION.txt
├── brusselator_2D_melnikov.m
├── glycolysis_7D_melnikov.m
├── mapk_cascade_8D_melnikov.m
├── circadian_10D_melnikov.m
├── melnikov_substrate.m
├── melnikov_enzymatic.m
├── melnikov_oregonator_final.m
├── brute_substrate.m
├── brute_enzymatic.m
├── brute_oregonator_final.m
├── benchmark_substrate_clean.m
├── benchmark_enzymatic_clean.m
├── benchmark_oregonator_final.m
├── substrate_rhs.m
├── enzymatic_rhs.m
├── oregonator_rhs.m
└── glycolysis_phase_portraits.m
```

## Quick Start
```matlab
results = brusselator_2D_melnikov();     % 2D  — fastest
results = glycolysis_7D_melnikov();      % 7D
results = mapk_cascade_8D_melnikov();    % 8D
results = circadian_10D_melnikov();      % 10D — extreme stiffness
```

## Citation
```bibtex
@article{sudarsan2026melnikov,
  title   = {An Efficient Adjoint-Free Melnikov--Lyapunov Diagnostic for
             Transition Detection in Slow--Fast Chemical Reaction Networks},
  author  = {Madhyannapu, Sri Venkata Durga Sudarsan and Pradheep Kumar, S.},
  journal = {Communications in Nonlinear Science and Numerical Simulation},
  publisher = {Elsevier},
  issn    = {1007-5704},
  year    = {2026},
  note    = {Under Review, Manuscript ID: CNSNS-D-26-00848}
}
```

## License
MIT License — provided for academic research purposes only.
