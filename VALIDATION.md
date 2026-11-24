# Validation and cost measurements

All runs on the same laptop as the COMSOL and BEM runs of the paper, Julia 1.11.6,
Gridap 0.18.12. The numbers below were produced by v0.1 of this code (project
`Explicit_FEM_Gridap_UGW`, 21-22 Nov 2026). A one-off regression check on 22 Nov 2026
ran v0.1 and v0.2 side by side on two small cases (intact pipe, nc = 24; notched pipe,
nc = 8): same number of steps, v_theta equal to 1e-9 relative (intact; the node
coordinates differ by round-off) and u_theta identical (notched). That script needed
the v0.1 project next to this one and is not part of the package. The notebooks in `notebooks/examples/`
regenerate the tables from `results/`.

Error measure, as in the paper: NRMSE = RMS of the difference in % of the range
(max - min) of the reference signal.

## 1. Intact pipe against the closed form

Case of the COMSOL intact-pipe cost series: L = 0.5 m, 30 kHz 5-cycle burst, 1 MPa,
t_max = 400 us with the end reflections, v_theta on the outer surface at z = 0.10,
0.25 and 0.40 m, against the method-of-images solution (`analytic_surface_velocity`).
One cell through the wall.

### Which direction limits order 2 (`run_intact_pipe.jl directions`)

| series | nc | nz | h_z (mm) | dofs | NRMSE % (0.10 / 0.25 / 0.40) |
|---|---|---|---|---|---|
| A | 12 | 100 | 5.00 | 45 225 | 0.034 / 0.043 / 0.044 |
| A | 24 | 100 | 5.00 | 88 641 | 0.029 / 0.037 / 0.038 |
| A | 48 | 100 | 5.00 | 175 473 | 0.0064 / 0.0080 / 0.0083 |
| A | 96 | 100 | 5.00 | 349 137 | 0.00069 / 0.00085 / 0.00091 |
| B | 48 | 50 | 10.00 | 88 173 | 0.025 / 0.032 / 0.032 |
| B | 48 | 200 | 2.50 | 350 073 | 0.0080 / 0.0101 / 0.0104 |
| B | 48 | 400 | 1.25 | 699 273 | 0.0041 / 0.0052 / 0.0053 |

Refining the circumference (A) lowers the error by an order of magnitude from 48 to
96 sectors; refining the axis (B) beyond 5 mm does little. The circumference limits
the accuracy, and straight-sided cells are enough: curved cells are not needed to
reach COMSOL's accuracy.

### At matched accuracy with COMSOL (L = 0.5 m)

| code | mesh | dofs | NRMSE % | loop time | threads |
|---|---|---|---|---|---|
| COMSOL, quartic DG (Elastic Waves, Time Explicit) | hx = 20 mm | 337 500 | 1.0e-3 to 2.3e-3 | 2 491 s | 8 |
| this code, order 2, lumped mass | nc = 96, h_z = 5 mm | 349 137 | 6.9e-4 to 9.1e-4 | 202 s | 8 |

The stable steps are similar (4.9e-8 s against 4.2e-8 s); the difference is the cost
per step. Most of COMSOL's explicit cost comes from its implementation (velocity-stress
DG with quartic elements), not from the explicit method.

### Cost against pipe length (`run_intact_pipe.jl length`, 8 threads)

| L (m) | dofs | steps | loop (s) | NRMSE % (L/5, L/2, 4L/5) | COMSOL (s) | COMSOL / this |
|---|---|---|---|---|---|---|
| 0.5 | 349 137 | 8 201 | 202 | 6.9e-4 / 8.5e-4 / 9.1e-4 | 2 491 | 12.3 |
| 1.0 | 696 537 | 16 402 | 868 | 9.5e-4 / 1.07e-3 / 1.48e-3 | 11 107 | 12.8 |
| 2.0 | 1 391 337 | 32 803 | 3 736 | 1.22e-3 / 1.40e-3 / 2.27e-3 | 37 489 | 10.0 |

Doubling L multiplies the loop time by 4.3 (about L^2.1): dofs and steps both grow
with L (t_max = 400 us L / 0.5 m), and the cost per dof-step rises slowly as the
vectors outgrow the cache. 8 threads give only about 2x over one thread: the sparse
product is bound by memory bandwidth (~19 GB/s on one thread, ~33 GB/s on eight).

## 2. Notched pipe

L = 1 m, 40 kHz burst, t_max = 500 us, order 2, two cells through the wall, 5 mm along
the pipe and 1.25 mm across the notch. The notch is cut exactly (cells removed and
volume checked before every run).

| mesh | dofs | steps | assembly | loop (8 threads) |
|---|---|---|---|---|
| nc = 48, 1.25 mm | 598 755 | 7 326 | 16 s | 547 s |
| nc = 96, 1.25 mm | 1 191 315 | 10 528 | 28 s | 1 540 s |
| nc = 48, 0.625 mm | 619 731 | 14 383 | 19 s | 1 185 s |

**Mesh convergence.** The three meshes differ by at most 0.02 % in u_theta and 3 ns in
delay; nc = 48 with 1.25 mm is converged.

**Against the hybrid BEM** (u_theta NRMSE in % of the BEM range; window: the scattered
wave only, in % of its own range):

| model | z = 0.75 m (echo) | z = 0.25 m (transmitted) | delay |
|---|---|---|---|
| full 3D BEM | 0.009-0.011 % | 0.009-0.015 % | 1-3 ns |
| this code, nc = 48 | 0.23-0.25 % (window 0.45-0.76 %) | 0.054-0.062 % (window 0.07-0.08 %) | 1-5 ns |
| implicit COMSOL | 0.49-0.67 % (window 1.8-2.7 %) | 1.87-1.90 % (window 1.8 %) | +0.13-0.16 us / +0.44 us |

This code and the BEM agree to 0.06 % downstream with no measurable delay; the
implicit COMSOL lags both by the same amount, so its 1.9 % is its own numerical
dispersion. The 0.25 % residual in the echo does not change with nc or with the axial
size across the notch: it is either the BEM (the hybrid and the full 3D model share
the notch discretization) or this code's far-field mesh acting on the scattered
flexural modes. Not separated.

**Longer pipe, L = 3 m.** Notch and sensors moved 1 m along z, so the signals must be
those of L = 1 m delayed by 319.45 us until end reflections arrive. 1 762 755 dofs,
12 012 steps, loop 2 722 s, peak memory 7.9 GB (about 4.5 GB per million dofs).

| comparison (0-385 us of 1 m time) | echo | transmitted | delay |
|---|---|---|---|
| 3 m against 1 m, same mesh | 0.06 % | 0.06 % | 1 ns |
| 3 m against hybrid BEM at 1 m | 0.10-0.12 % | 0.07-0.08 % | 1-6 ns |

One more metre of propagation adds 0.06 %: the numerical dispersion is negligible.

## 3. Cost against the hybrid BEM, and its caveats

Measured, same laptop, 8 threads, notched pipe at L = 1 m: this code 547 s, hybrid BEM
4 216 s (53 values of s x 80 s, 8.3 GB), full 3D BEM 7 764 s (9.6 GB), implicit COMSOL
23 084 s (20 cores).

Long pipe, measured 23 Nov 2026 with the BEM script's cost mode (3 of the 206 values of
s): the hybrid BEM at L = 10 m (notch and sensors 4 m farther along the pipe, t_max =
2 ms) needs 61.5 s per value of s and 5.8 GB, i.e. about 12 700 s (3.5 h) for the
transient. Re-measured the same way at L = 1 m: 61.1 s per value of s and 5.4 GB, i.e. about
3 240 s (54 min) for the 53 values; the 70.3 min (80 s per value, 8.3 GB) of the paper
were measured with an earlier version of the BEM code. The cost per value of s is the
same at 1 and 10 m, as the method predicts, and at 1 m this code is 5.3-5.9 times faster
than the hybrid BEM (547-607 s). Full 3D BEM, same mode: 87.1 s per value of s and
9.5 GB (about 77 min), so the hybrid BEM is 1.43 times faster with 43 % less memory at 1 m
(the paper had 129.4 min, 147 s per value, measured with the earlier code). This code at 10 m: about 5.8 M dofs, 26 GB and 6 h,
extrapolated from 1 and 3 m; it does not fit in the laptop.

**Run-time variability.** v0.2 reproduces every v0.1 result to round-off (signals) and
the same NRMSE, but the loop times measured on the evening of 22 Nov were 1.1-1.7 times
those of the morning on the same machine: notched pipe 607 s (547 s), intact pipe 352 /
1 359 / 5 707 s at 0.5 / 1 / 2 m (202 / 868 / 3 736 s). The regression check showed v0.1
and v0.2 equally fast in the same session, so the difference is the state of the laptop
(load, temperature, power management). The tables above use the morning series, run back
to back at a constant rate; timings for the paper should come from one such session.

Two things must enter any comparison at other lengths:

- **The number of Laplace values grows with the observation time**, N ~ 1.6 f_max t_max.
  The cost per value of s of the hybrid BEM hardly depends on L (8.4-9.6 s for the
  intact pipe from 0.5 to 8 m), but its total cost grows linearly with t_max. An
  earlier comparison in this project multiplied the cost per value by a fixed N = 32 at
  every length, while the FEM ran t_max proportional to L; that overstated the
  advantage of the BEM in long pipes. With N proportional to t_max, the intact-pipe
  ratio FEM / hybrid BEM is 0.75, 1.6 and 3.3 at 0.5, 1 and 2 m (measured FEM,
  per-value BEM cost from the paper).
- **Memory.** This code needs about 4.5 GB per million dofs (stiffness assembled), so
  the laptop (about 12 GB free) holds a notched pipe of about 4.5 m. The hybrid BEM
  needs about 8.5 GB at any length. A matrix-free explicit code would need much less,
  which a fair comparison should state.

## 4. Errors found during the validation

1. `num_nodes` and `symmetric_part` are not exported at the top level of Gridap 0.18;
   `symmetric_gradient` is.
2. The load was first applied on `tag_26`, which is the outer surface, not the end:
   Gridap numbers the faces from the last parametric coordinate to the first, so the
   end z = L is `tag_22`. The first complete run came out 3.7 times too large.
3. A `CartesianDiscreteModel` with a curvilinear map integrates the faces normal to the
   axis as if every cell had the width of its inner edge (constant Jacobian), short by
   2 ri / (ri + ro) = 0.929 with one radial layer. Volumes and lateral faces are exact,
   so nothing looks wrong: the load is 7.5 % short and the NRMSE stops at a
   mesh-independent 1.3 %. Fixed by converting to an `UnstructuredDiscreteModel`.
4. The analytic peak velocity used a factor 0.5 too many; the burst peak is now sampled.
5. The lumped mass was first obtained by summing the rows of the assembled consistent
   mass, which did not fit in memory at L = 2 m. It is now assembled directly as a
   vector (identical to round-off, `test/runtests.jl`).
6. The power iteration for the time step had no fixed seed, so a resumed run could get
   a slightly different step. Seeded.
7. Windows power management slowed a one-thread run by up to 2.8x partway through; cost
   runs are made with the laptop plugged in, in maximum-performance mode, and the rate
   is checked in the checkpoint log.
