# Stage 07 Blocker - TYPE6 off-axis D3039 orientation mismatch

Author: J.C. Vaught

The post-matrix `LAW12 + TYPE6/SOL_ORTH` D3039 decks now start, run, write animation output, convert through `anim_to_vtk`, and post-process with PyVista. The remaining blocker is numerical, not starter parsing.

Observed results from `results/timeseries.csv`:

- Run 7A, 0 deg: `E_FEM = 170.22 GPa` vs `E1 = 171.40 GPa`, error `0.69%`, PASS.
- Run 7B, 45 deg: `E_FEM = 171.15 GPa` vs analytic `Ex(45) = 13.28 GPa`, error `1189%`, FAIL.

The 45 deg value is E1-like, indicating the verified `TYPE6` solid property path is not applying the in-plane material-frame rotation in the way required by the stage 07 off-axis D3039 gate. I also probed alternate `TYPE6` `Ip`, `Iorth`, `Phi`, skew-frame, and reference-vector combinations; none recovered the 45 deg plane-stress transformed modulus within the 2% tolerance.

Verdict: `FAIL` due to a real off-axis modulus mismatch after the matrix substitution, not due to `ERROR 3047` or another starter parse blocker.
