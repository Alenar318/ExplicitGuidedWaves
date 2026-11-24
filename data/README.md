# Reference data

Inputs used only for comparison; nothing in `src/` depends on them except
`reference_signals`, which reads `bem/` and `comsol_implicit/`.

## `comsol_implicit/`

Implicit COMSOL solution of the notched pipe, L = 1 m (validated reference of the
paper). Copied unchanged from `BEM_HIBRIDO_AXI_3D_JULIA/datos/comsol/`.

- `entalla_T01_L1_40kHz_x{0.25,0.75}_th{0,90,180}.txt`: displacement (m) on the outer
  surface at z = 0.25 and 0.75 m, theta = 0, 90, 180 deg. Columns: time (s), u_r,
  u_theta, u_z; 1 us sampling, 0-500 us. COMSOL 6.0, header lines start with `%`.
- `resumen_comsol_entalla_T01_L1_40kHz.txt`: COMSOL's summary of that run (mesh,
  solver, time).

## `bem/`

Hybrid and full 3D BEM signals of the same case, inverted with Durbin's method from the
Laplace-domain solutions of the BEM project (`hibrido_entalla_N53.jls`,
`completo_entalla_N53.jls`, N = 53 values of s, T0 = 4.0e-4 s). Exported on 22 Nov 2026
with `bem/export_bem_notch.jl` (kept here as the record of how they were made; it ran
from the v0.1 project and calls `BEMAxi.durbin` of the BEM project).

- `bem_hibrido_entalla_x{0.25,0.75}_th{0,90,180}.txt`, `bem_completo_entalla_...`:
  displacement (m). Columns: time (s), u_r, u_theta, u_x (axial); 1 us sampling.

## `reference_costs.csv`

Run times and peak memory of the other codes on the same laptop, with their source.
COMSOL explicit: `Comparacion_FEM_explicito/README.md` (memory from the COMSOL run log);
implicit COMSOL: notched-pipe section of the paper (`articulo_hibrido_axi_3D_CAS_v07.tex`).
BEM: two sets of rows, those of the v07 paper (an earlier version of the BEM code) and
those measured in cost mode on 23 Nov 2026 with the current code, which the revised
paper uses; `02_notched_pipe.ipynb` shows only the current ones. Empty fields: not measured.
