Reproducibility files for
"Domination Number and A_alpha-Spectral Properties of Bipartite Graphs"

FILES
-----
reproduce_section4.py
    Python script reproducing the computational consistency checks in Section 4.

requirements.txt
    Minimal Python package requirements.

RUN
---
python -m pip install -r requirements.txt
python reproduce_section4.py

EXPECTED MANUSCRIPT CHECKPOINTS
-------------------------------
Exhaustive Graph Atlas:
- connected nonisomorphic bipartite graphs, n=2,...,7:
  1, 1, 3, 5, 17, 44 (total 71)
- gamma=2 counts for n=4,...,7:
  2, 4, 15, 31 (total 52)
- alpha values:
  0, 0.1, 0.3, 0.5, 0.51, 0.7, 0.9, 0.99
- 568 spectral tests
- no violation (tolerance 1e-9)

Random experiment:
- NumPy default RNG seed 20260926
- n = 10, 12, 14, 16
- balanced bipartitions
- edge probabilities = 0.15, 0.30, 0.50, 0.70
- 40 generated graphs for each (n, probability)
- disconnected realizations discarded
- manuscript retained counts:
  n=10: 78
  n=12: 80
  n=14: 89
  n=16: 97
  total: 344
- 2752 spectral tests
- no violation (tolerance 1e-9)

Together the two experiments contain 3320 graph-parameter tests.

NOTE
----
The computations are supplementary numerical consistency checks. They are not
used to prove the mathematical results.
