# Stage 14 - Compression After Impact (ASTM D7137/D7137M)

Author. J.C. Vaught
Date. 2026-04-29
Status. Specification, pre-run. Chained from Stage 13 (Low-velocity impact, ASTM D7136).

---

## 1. Stage summary

This stage closes the LVI to CAI loop on the same coupon as Stage 13. A
$100 \times 150 \times 4$ mm 24-ply $[0/+45/-45/90]_{3s}$ T800-SC-24K /
P2362W laminate that has just absorbed a drop-weight impact (Stage 13,
ASTM D7136 \cite{ASTM_D7136}) is loaded edgewise in compression inside
an ASTM D7137/D7137M-17 anti-buckling fixture \cite{ASTM_D7137} until
catastrophic load drop. The damage state from Stage 13 (intra-laminar
Hashin damage variables, plastic and damage history on the
\texttt{/PROP/TYPE14} solid bricks, eroded element flags, and the
inter-laminar cohesive-zone delamination pattern) is preserved exactly
by reading Stage 13's final state file (\texttt{rootname\_0001.sta},
written by \texttt{/STATE} cards at the end of the LVI engine run) and
seeding the Stage 14 starter via \texttt{/INIBRI/...} initial-brick
keywords plus an explicit \texttt{/RUN/<rootname>/2} engine restart.
The chaining mode is explicit-explicit, in line with audit cross-cut D
\cite{OpenRadiossEndToEndAudit2026}, which records that explicit-explicit
state-file restart is documented (Spring-back tutorial,
\cite{OpenRadiossConfluenceSpringback2024}) while the
implicit-after-explicit direction is not. The CAI run is therefore
quasi-static explicit, with end-shortening velocity slow enough that
the kinetic energy stays below five percent of the internal energy
throughout the loading ramp; this is checked from the global energy
time-history file (\texttt{T01}). The pass criterion is the residual
compressive strength, defined as the peak edge-load divided by the
gross cross-section before catastrophic load drop, within ten percent
of the canonical Lopes and Camanho 2009 reference \cite{Lopes2009LVIPart2}
for a comparable quasi-isotropic dispersed-stacking laminate at the
same impact energy. The Soutis and Curtis 1996 sublaminate-buckling
closed form \cite{SoutisCurtis1996} provides the corroborating
analytical envelope. The UofSC AP-PLY companion paper of Kodagali et al.
2024 \cite{Kodagali2024LowVelocityImpact} provides the architecture-
specific cross-check on the same T800-SC-24K / P2362W material card.

---

## 2. Geometry

The specimen is the unmodified Stage 13 coupon. Mesh continuity is a
hard requirement of OpenRadioss restart, so no remeshing, no node
renumbering, and no element renumbering is permitted between the two
stages. The Stage 14 starter deck includes the Stage 13 final mesh by
referencing the same \texttt{Include} files and the same node, brick,
and part identifiers.

### 2.1 Dimension table

| Symbol | Quantity | Value | Source |
|---|---|---|---|
| $L$ | Specimen length (load axis) | $150.0$ mm | ASTM D7137 \cite{ASTM_D7137} sec.~8 |
| $W$ | Specimen width | $100.0$ mm | ASTM D7137 \cite{ASTM_D7137} sec.~8 |
| $t$ | Specimen thickness | $4.0$ mm | Kodagali 2023 \cite{Kodagali2023MesoArchitectured} |
| $n_{\text{ply}}$ | Number of plies | $24$ | $[0/+45/-45/90]_{3s}$ |
| $t_{\text{ply}}$ | Cured ply thickness | $\approx 0.167$ mm | $t/n_{\text{ply}}$ |
| $A_g$ | Gross compressive cross-section | $W \cdot t = 400.0$ mm$^2$ | derived |
| $L_{\text{free}}$ | Free length between knife edges | $115.0$ mm | ASTM D7137 \cite{ASTM_D7137} Fig.~1 |

Note that the ASTM D7137 standard specifies $4 \times 6$ inch
($\approx 101.6 \times 152.4$ mm) as the canonical coupon size; the
master plan brief uses the close metric value $100 \times 150$ mm in
keeping with the SI convention used throughout the project. The
deviation is below two percent and is well inside the standard's
geometric tolerance.

### 2.2 Coordinate system

The same Stage 13 right-handed Cartesian frame is reused. The $x$-axis
is the long (load) axis, the $y$-axis is the short transverse in-plane
axis, the $z$-axis is the through-thickness axis. The impact face is
$z = +t/2$, the back face is $z = -t/2$. The compressive load axis is
$\pm x$.

### 2.3 ASCII sketch of the D7137 fixture footprint

```
   side support                                   side support
   (knife edge,                                   (knife edge,
    y = +W/2)                                      y = -W/2)
        v                                                v
        |                                                |
   +---------+--------------------------------+---------+
   |         |    upper end-load platen       |         |
   |         |  (+x face, uniform u_x = -d)   |         |
   |         +--------------------------------+         |
   |                                                    |
   |         (free length L_free = 115 mm,              |
   |          free in y on the long sides)              |
   |                                                    |
   |         +--------------------------------+         |
   |         |    lower end-load platen       |         |
   |         |  (-x face, fixed u_x = 0)      |         |
   +---------+--------------------------------+---------+
        ^                                                ^
        |                                                |
   knife-edge contact line                          knife-edge contact line
   (z = -t/2 and z = +t/2,                          (mirror image at y = -W/2)
    y = +W/2,
    only u_z = 0 enforced)
```

Through-thickness extent ($z$) is the 4 mm laminate; the fixture knife
edges contact only the surface nodes on $z = \pm t/2$, leaving the
in-plane $u_y$ DOF free along the long edges so that Poisson contraction
is not over-constrained.

---

## 3. Mesh

The mesh is identical to Stage 13. No re-meshing. No re-numbering. No
change of element formulation. The same \texttt{/PROP/TYPE14} solid
brick property and the same per-ply orientation are used; only the
boundary conditions, the loading function, and the analysis-control
cards change between Stage 13 and Stage 14.

### 3.1 Element budget (inherited from Stage 13)

| Region | Element type | Element count | Through-thickness layers |
|---|---|---|---|
| Laminate plies | HEXA8, \texttt{/PROP/TYPE14} | $\approx 24 \times N_{xy}$ | one element per ply |
| Inter-laminar cohesive layer | HEXA8 zero-thickness, \texttt{/MAT/LAW83} or \texttt{/MAT/LAW117} | $23 \times N_{xy}$ | between adjacent plies |

$N_{xy}$ is the in-plane brick count from Stage 13. Stage 14 adds zero
new elements; it adds only knife-edge \texttt{/RBODY} masters and an
end-loading platen represented as a kinematic \texttt{/RBODY} on the
$\pm x$ end faces.

### 3.2 Why mesh continuity is mandatory

Restart in OpenRadioss is implemented as a node-by-node and
element-by-element state read. The state file
\texttt{rootname\_0001.sta} written by Stage 13's \texttt{/STATE/...}
cards stores per-element stress tensors, per-element history variables
(damage, plastic strain), per-element deletion flags, and per-node
positions and velocities, all keyed by integer ID. If Stage 14 changes
any ID, adds an element, deletes an element, or alters element
formulation, the state read fails or, worse, silently maps the wrong
state onto the wrong element. The restart literature
\cite{OpenRadiossConfluenceSpringback2024} and the audit
\cite{OpenRadiossEndToEndAudit2026} both flag this as the single
hardest-to-debug class of restart failure.

---

## 4. Boundary conditions and loading (anti-buckling fixture approximation)

The ASTM D7137 fixture is a four-piece steel anti-buckling jig.
The two long sides of the coupon are pressed between knife-edge
supports that prevent out-of-plane motion ($u_z = 0$ along
$y = \pm W/2$ surface-to-surface contact lines) but leave the in-plane
$u_x$ and $u_y$ DOFs free. The two short ends of the coupon sit
between flat steel platens that are constrained to remain parallel
during loading. The bottom platen is fixed in $x$; the top platen is
displaced in $-x$ to compress the coupon. ASTM D7137
\cite{ASTM_D7137} Figure~1 is the canonical sketch.

In OpenRadioss with solid elements only (the master-plan hard
constraint), the fixture is approximated as follows.

### 4.1 Knife-edge approximation on the long sides

On each long side ($y = +W/2$ and $y = -W/2$), the line of surface
nodes shared by the impact face and the back face is restrained in the
through-thickness direction with a \texttt{/BCS} card that sets
$u_z = 0$ and leaves $u_x$ and $u_y$ free. The knife edges are
modelled as kinematic restraints rather than as deformable steel
parts; the steel of the real fixture is at least one hundred times
stiffer than the through-thickness laminate stiffness, so the
restraint approximation is a sound lower-bound on the buckling load
\cite{ASTM_D7137}. The line restraint covers only the interior of the
free length $L_{\text{free}} = 115$ mm; nodes at the very ends of the
long sides (within the platen footprint) are released, so that the
transition between knife-edge support and end-platen clamping is
geometrically clean.

### 4.2 End-platen approximation on the short ends

On each short end ($x = \pm L/2$), all surface nodes on that face are
collected into a kinematic \texttt{/RBODY} master. The bottom platen
master node ($x = -L/2$) is fixed in all six DOFs via a \texttt{/BCS}
card. The top platen master node ($x = +L/2$) is fixed in
$u_y, u_z, \theta_x, \theta_y, \theta_z$ and is driven in $-u_x$ by an
\texttt{/IMPDISP} imposed-displacement card, ramped linearly from zero
to the target end shortening. The platen \texttt{/RBODY} replaces the
steel platen of the physical fixture by an idealized rigid plane; this
matches ASTM D7137 \cite{ASTM_D7137} which specifies that the platens
must remain parallel to within $0.05$ mm out-of-parallelism.

### 4.3 Loading rate (kinetic-energy budget)

ASTM D7137 \cite{ASTM_D7137} specifies a quasi-static rate of
$1.25$ mm/min (head displacement). For an explicit dynamic solver
this is impractically slow; the explicit time step is on the order of
$\Delta t \sim 10^{-7}$ s for a 4 mm laminate, giving a real-time
ramp of order $10^{12}$ steps which is infeasible. The standard
explicit quasi-static workaround is to compress at a velocity
$v_{\text{end}}$ that is slow enough on the natural-period scale of
the coupon. The first compressive natural period of an undamaged
$100 \times 150 \times 4$ mm IM7-class laminate is $T_1 \approx 1.5$ ms
(estimated from a slender-plate Euler-buckling argument with
$E_x \approx 70$ GPa). A safe explicit ramp imposes the full
expected end shortening (about $0.5$ mm at the experimentally measured
CAI strength of $\approx 200$ MPa) over $50 T_1 \approx 75$ ms, giving
$v_{\text{end}} \approx 0.5/0.075 = 6.7$ mm/s. This is approximately
$300\times$ faster than the ASTM rate but $200\times$ slower than the
laminate's first natural period, well inside the quasi-static envelope.

The quasi-static condition is then verified post-hoc on the actual
run by checking, from the global \texttt{T01} energy history, that

$$
\text{KE}(t) \le 0.05 \cdot \text{IE}(t)
$$

where KE is the global kinetic energy and IE is the global internal
(strain plus damage dissipation) energy, throughout the loading ramp
\cite{ASTM_D7137,Lopes2009LVIPart2}. If this is violated the run is
re-launched with $v_{\text{end}}$ halved.

### 4.4 Initial conditions (this is the restart)

At $t = 0$ of Stage 14, all nodal velocities, all element stresses,
and all element history variables are inherited from Stage 13's final
state. The Stage 14 starter does not zero any state; the
\texttt{/INIBRI/...} cards (Section 6.3 below) are the mechanism by
which Stage 13's stress, plastic strain, damage, and erosion-flag
fields are mapped onto Stage 14's identical mesh. The end-platen
\texttt{/IMPDISP} ramp begins at $t = 0$ from zero displacement.

---

## 5. Material card (inherited from Stage 13, not redefined)

The composite material law is identical to Stage 13: \texttt{/MAT/LAW25}
with \texttt{Iform = 1} (CRASURV formulation), with the \texttt{/FAIL/HASHIN}
failure card and \texttt{Ifail} flag set so that elements are eroded
when the criterion is met. The inter-laminar cohesive layer uses
\texttt{/MAT/LAW83} (or \texttt{/MAT/LAW117} on a recent build) with
the mode-I and mode-II fracture-energy values from Kodagali 2023
\cite{Kodagali2023MesoArchitectured}. Both are confirmed solid-element
compatible by the Altair composite material introduction page
\cite{AltairRadiossCompositeIntro,AltairRadiossLAW25}.

The Stage 14 starter does not redefine these material cards. It
\texttt{\#include}s the same Stage 13 material file. Redefining a
material card with even subtly different numerical values silently
corrupts the restart because the integration of the stored history
variable (e.g. damage $D \in [0, 1]$) requires the same yield surface
and the same softening modulus that produced it.

| Material | Card | Source |
|---|---|---|
| T800-SC-24K / P2362W ply (LAW25 CRASURV) | \texttt{/MAT/LAW25} \cite{AltairRadiossLAW25} | Kodagali 2023 \cite{Kodagali2023MesoArchitectured} |
| Hashin failure with erosion | \texttt{/FAIL/HASHIN}, \texttt{Ifail = 2} \cite{AltairRadiossFailHashin} | Stage 13 inherited |
| Inter-laminar cohesive (mode I/II) | \texttt{/MAT/LAW83} or \texttt{/MAT/LAW117} | Stage 13 inherited |
| Solid property | \texttt{/PROP/TYPE14} \cite{AltairRadiossPropType14} | Stage 13 inherited |

---

## 6. OpenRadioss deck skeleton (restart cards)

The OpenRadioss restart workflow is a two-stage starter-engine pair
per run, with the restart linkage between runs implemented through
state files written by \texttt{/STATE} cards in Stage 13's engine
deck and read by \texttt{/INIBRI/...} cards in Stage 14's starter
deck, plus a \texttt{/RUN} card in Stage 14's engine deck that
advances the run counter. Card names below are confirmed against
the Altair Radioss reference manual pages cited in
\cite{OpenRadiossEndToEndAudit2026} and the Spring-back restart
tutorial \cite{OpenRadiossConfluenceSpringback2024}.

### 6.1 Files produced by Stage 13 (the upstream stage)

Stage 13's engine deck, at the end of the LVI run (just before the
impactor rebounds clear), writes:

- \texttt{rootname\_0001.sta} -- the binary state file; produced by
  the \texttt{/STATE/BRICK/FULL} card requested at the final time
  step. This file holds, per brick, the Cauchy stress tensor, the
  internal-energy density, the strain tensor, the integration-point
  history variables (damage indicators, plastic strain), and the
  element-deletion flag.
- \texttt{rootname\_0001.rst} -- the binary restart file (engine-side
  restart format), the runtime-state companion of the \texttt{.sta}
  file used internally by the engine to resume kinematics.
- \texttt{T01} time-history -- written for the validation check.

### 6.2 Stage 14 starter deck (\texttt{rootname2\_0000.rad})

The starter deck for Stage 14 is structurally equivalent to Stage 13's
starter, with the following block-level changes.

```
#include 'common_geom.inc'         // mesh, parts, properties, materials
#include 'common_mat.inc'          // LAW25, LAW83 cards (unchanged)

/INIBRI/STRS_F/<part_id>           // initial brick stress tensor
   <reads rootname_0001.sta>
/INIBRI/STRA_F/<part_id>           // initial brick strain tensor
   <reads rootname_0001.sta>
/INIBRI/EPSP/<part_id>             // initial plastic strain
   <reads rootname_0001.sta>
/INIBRI/AUX/<part_id>              // initial auxiliary variables
   <reads rootname_0001.sta>       // damage variables, history flags
/INIBRI/THICK/<part_id>            // initial cohesive thickness state

/INIVEL/...                        // initial nodal velocity (from .sta)

// Boundary conditions for CAI compression (replace LVI BCs)
/RBODY/...   id=1001               // bottom-platen rigid body, x = -L/2
/RBODY/...   id=1002               // top-platen rigid body, x = +L/2
/BCS/...                           // u_z = 0 on long-side knife-edge nodes
/BCS/...                           // bottom-platen master fully fixed
/BCS/...                           // top-platen master constrained except u_x

// Loading
/IMPDISP/<top-platen-master>/<func_id>/<dir = -X>/<scale>
/FUNCT/<func_id>                   // ramp 0 to d_target over t_ramp = 75 ms

// Output requests
/TH/RBODY/...   id=1002            // top-platen reaction force history
/TH/PART/...                       // global energy history (KE, IE)
/STATE/BRICK/FULL                  // (optional) write a Stage 14 final state
                                   // for any downstream chain
```

The Stage 14 starter is *not* an "in-place" restart -- it is a fresh
starter run that re-initializes the mesh and then overlays Stage 13's
state via \texttt{/INIBRI/...} reads from the Stage 13
\texttt{.sta} file. This is the documented OpenRadioss pattern
\cite{OpenRadiossConfluenceSpringback2024}.

### 6.3 Stage 14 engine deck (\texttt{rootname2\_0001.rad})

```
/RUN/rootname2/1                   // declare engine run number 1 of rootname2
/STOP/EMAX/1.05                    // stop on energy or load drop
/DT/BRICK/CST/0.667e-7             // explicit time step (Stage 13 inherited)
/PRINT/-100                        // print every 100 cycles
/ANIM/DT                           // animation cadence
/TH/...                            // time-history cadence
/STATE/DT                          // optional final-state write at end of run
```

### 6.4 Card-by-card audit of the restart cards

| Card | Function | Source |
|---|---|---|
| \texttt{/STATE/BRICK/FULL} (Stage 13 engine) | Write final brick stress, strain, history at end of LVI run | \cite{OpenRadiossConfluenceSpringback2024,OpenRadiossEndToEndAudit2026} |
| \texttt{/INIBRI/STRS\_F} (Stage 14 starter) | Initial Cauchy stress tensor on every brick of a part | \cite{OpenRadiossConfluenceSpringback2024} |
| \texttt{/INIBRI/STRA\_F} | Initial total-strain tensor on every brick | \cite{OpenRadiossConfluenceSpringback2024} |
| \texttt{/INIBRI/EPSP} | Initial equivalent plastic strain | \cite{OpenRadiossConfluenceSpringback2024} |
| \texttt{/INIBRI/AUX} | Initial auxiliary history (damage flags, Hashin damage indicators, eroded-element flag) | \cite{OpenRadiossConfluenceSpringback2024} |
| \texttt{/INIBRI/THICK} | Initial cohesive layer thickness (for cohesive solid layer) | \cite{OpenRadiossConfluenceSpringback2024} |
| \texttt{/INIVEL} | Initial nodal velocity (read from \texttt{.sta}) | \cite{OpenRadiossConfluenceSpringback2024} |
| \texttt{/RUN/rootname2/1} (Stage 14 engine) | Engine run-number directive; maps to the Stage 14 deck identity | \cite{OpenRadiossInstall2026} |

### 6.5 Restart-format quirks (documented in this run)

Three quirks must be honored or the restart silently corrupts state.

First quirk. Element ID-set, node ID-set, and part ID-set must
match Stage 13's exactly. The \texttt{.sta} file is keyed on integer
IDs, not on coordinates. Renumbering during \texttt{inp2rad}
conversion or during a \texttt{/REORDER/.../} card is a footgun.

Second quirk. \texttt{/INIBRI/AUX} maps the auxiliary-variable
array index by index, and the index layout is law-specific. For
\texttt{/MAT/LAW25} the auxiliary array stores Hashin damage
variables, plastic strain history, and the eroded-element flag in a
fixed slot ordering that the \texttt{.sta} writer encodes; the reader
must use the same law on the same property or the slot mapping is
silently wrong. The cited Spring-back tutorial warns about this for
LAW2 \cite{OpenRadiossConfluenceSpringback2024}; the same warning
applies to LAW25.

Third quirk. The Spring-back tutorial documents the
explicit-then-implicit case (forming-then-springback). The
audit \cite{OpenRadiossEndToEndAudit2026} cross-cut D notes that the
reverse direction (implicit-then-explicit) and the implicit
restart from an explicit state file are *not* documented. This stage
honors that note by chaining explicit-explicit only. The CAI is
quasi-static in the physical experiment but is implemented here as a
slow-velocity explicit ramp with the kinetic-energy check of
Section 4.3.

### 6.6 Fallback if state-file restart proves problematic

If the \texttt{/INIBRI/...} reads either fail (missing \texttt{.sta}
field) or silently corrupt state (wrong eroded-element mapping), the
fallback is a single combined LVI-then-CAI deck. In that fallback
the Stage 13 LVI run is re-executed in the Stage 14 deck (no
restart), and immediately after the impactor has rebounded clear
(measured by the impactor-laminate contact force returning to zero on
the \texttt{T01} curve), the boundary conditions are switched in the
same engine deck via a sensor-triggered \texttt{/SENSOR/...} card
that activates the CAI \texttt{/IMPDISP} on the top platen and
de-activates the impactor \texttt{/RBODY}. The fallback wastes the
LVI compute time on every CAI parameter sweep but avoids the
restart-format failure modes entirely. The audit
\cite{OpenRadiossEndToEndAudit2026} explicitly notes that this
single-deck fallback is the safer path when state-file restart is
not yet known to work for the chosen material law.

---

## 7. Reference solution

The canonical FEM benchmark for CAI on a quasi-isotropic dispersed
stacking laminate after a $6.7$ J/mm impact is Lopes and Camanho 2009
\cite{Lopes2009LVIPart2}. Their Part II numerical-simulations paper
reports residual compressive strength figures for the same family of
quasi-isotropic stacking sequences as the master-plan brief
$[0/+45/-45/90]_{3s}$. The reported residual compressive strength
across the dispersed-stacking sequences they simulate falls in the
band $\sigma_{\text{CAI}} = 200 \pm 30$ MPa for IM7/8552, which is
the closest-relative material to the master-plan T800-SC-24K /
P2362W card. Lopes and Camanho 2009 Figure~6 (load--displacement)
and Table~3 (residual strength values per stacking sequence) are
the comparison targets for this stage \cite{Lopes2009LVIPart2}.

The Soutis and Curtis 1996 sublaminate-buckling closed form
\cite{SoutisCurtis1996} provides the analytical envelope. Their
prediction depends on the projected damage diameter $d_d$ from
Stage 13. For a circa-$30$ mm post-impact damage diameter (typical
for the 6.7 J/mm impact energy class) on the 4 mm 24-ply laminate,
their semi-empirical formula gives
$\sigma_{\text{CAI}} \approx 180\text{--}220$ MPa, in agreement
with the Lopes-Camanho band.

The UofSC AP-PLY companion paper of Kodagali et al. 2024
\cite{Kodagali2024LowVelocityImpact} reports CAI residual strength
on the matching T800-SC-24K / P2362W material for hybrid
pseudo-woven (PW) and baseline (non-PW) coupons. Their baseline
24-ply data is the numerical match for this stage. Even though
Kodagali's coupons are the AP-PLY architecture, not a flat tape
laminate, the baseline (non-PW) tape-laminate datapoint they report
is directly comparable.

---

## 8. Validation success criterion

The pass criterion is one quantitative measurement on the load-time
record from \texttt{T01}.

Define the residual compressive strength as

$$
\sigma_{\text{CAI}} = \frac{F_{\text{peak}}}{A_g}
$$

where $F_{\text{peak}}$ is the peak compressive reaction force on the
top-platen \texttt{/RBODY} (taken before the catastrophic load drop
defined as $F$ falling below $0.6 F_{\text{peak}}$ within
$\le 0.5$ ms of the peak), and $A_g = W \cdot t = 400$ mm$^2$ is the
gross undamaged cross-section (per the ASTM D7137 convention of
dividing by gross, not net, area \cite{ASTM_D7137}).

Pass criterion. $\sigma_{\text{CAI}}$ from the FEM run within ten
percent of the Lopes and Camanho 2009 reference value
\cite{Lopes2009LVIPart2} for the matching impact energy and stacking
sequence

$$
\frac{|\sigma_{\text{CAI}}^{\text{FEM}} - \sigma_{\text{CAI}}^{\text{ref}}|}{\sigma_{\text{CAI}}^{\text{ref}}} \le 0.10
$$

with the reference value taken as the centerline of the
$200 \pm 30$ MPa band.

Quasi-static check. The kinetic-to-internal energy ratio from the
\texttt{T01} global-energy history must satisfy
$\text{KE}(t) / \text{IE}(t) \le 0.05$ throughout the loading ramp,
per Section 4.3. Failure of this check invalidates the residual
strength reading and triggers a re-run with $v_{\text{end}}/2$.

Sanity checks. The buckling mode at peak load must be local
sublaminate buckling at the impact site (visible as the back-face
$z$-displacement field on \texttt{/ANIM} output before the load drop),
not global Euler buckling of the coupon. Global buckling indicates
either an under-stiffened knife-edge approximation (Section 4.1)
or an unrealistically low residual strength prediction.

---

## 9. Toolchain runner

The runner.py companion script implements the full Stage 14 chain.
The script lives at
\texttt{tests/stage\_14\_CAI\_D7137/runner.py}.

### 9.1 Top-level steps performed by runner.py

1. Locate Stage 13's output. Search
   \texttt{tests/stage\_13\_LVI\_D7136/out/} for the most recent
   completed LVI run; identify
   \texttt{rootname13\_0001.sta},
   \texttt{rootname13\_0001.rst}, and the time-history
   \texttt{T01}; verify all three exist and that the LVI run
   exited cleanly (final cycle is past the impactor-rebound
   contact-force-zero crossing).

2. Render the Stage 14 starter deck from a Jinja2 template.
   The template ingests Stage 13's \texttt{rootname13} so that
   \texttt{/INIBRI/...} cards point to the correct \texttt{.sta}
   file by name and so that the Stage 14 mesh include files
   \texttt{\#include} the same Stage 13 mesh blocks.

3. Render the Stage 14 engine deck from a Jinja2 template,
   parameterized by $v_{\text{end}}$, ramp time $t_{\text{ramp}}$,
   target end shortening $d_{\text{target}}$, animation cadence
   $\Delta t_{\text{anim}}$, and the time-history cadence
   $\Delta t_{\text{TH}}$.

4. Stage the deck plus all \texttt{\#include}d files plus the
   Stage 13 \texttt{.sta} into a Lima-shared working directory
   on the macOS host. The working directory layout matches
   what the OpenRadioss starter expects: all auxiliary files sit
   next to the rootname-based deck.

5. Execute starter inside Lima Apptainer
   \cite{OpenRadiossInstall2026,OpenRadiossDiscussion2125}:

   ```bash
   limactl shell apptainer -- \
       /OpenRadioss/exec/starter_linuxa64 \
       -i rootname14_0000.rad -nt $NTHREADS
   ```

6. Execute engine inside Lima Apptainer:

   ```bash
   limactl shell apptainer -- \
       /OpenRadioss/exec/engine_linuxa64 \
       -i rootname14_0001.rad -nt $NTHREADS
   ```

7. Parse the resulting \texttt{T01} time-history with
   \texttt{Vortex-Radioss} \cite{VortexRadioss2024} or by direct
   Fortran-binary read; extract the top-platen reaction force time
   series $F(t)$ and the global energy histories KE$(t)$ and IE$(t)$.

8. Compute $F_{\text{peak}}$ as the maximum compressive reaction
   force on the top-platen \texttt{/RBODY} that occurs before the
   load-drop criterion ($F$ falls below $0.6 F_{\text{peak}}$ in
   $\le 0.5$ ms). Compute
   $\sigma_{\text{CAI}} = F_{\text{peak}} / A_g$.

9. Verify the quasi-static condition KE$(t) / $IE$(t) \le 0.05$
   throughout the loading ramp. If violated, halve $v_{\text{end}}$
   and return to step 3.

10. Convert \texttt{rootname14\_0001.anim} to VTKHDF using the
    Kitware converter \cite{KitwareOpenRadiossVTKHDF2024} for
    diagnostic visualization.

11. Write a CSV record of $F(t)$, KE$(t)$, IE$(t)$ for
    downstream Typst + CeTZ plotting per the project plotting
    convention.

12. Compare $\sigma_{\text{CAI}}^{\text{FEM}}$ against the
    Lopes-Camanho band; emit pass/fail; emit a summary log.

### 9.2 Inputs to runner.py

| Argument | Default | Purpose |
|---|---|---|
| \texttt{--stage13-out} | \texttt{tests/stage\_13\_LVI\_D7136/out} | Path to Stage 13 output directory |
| \texttt{--rootname13} | \texttt{lvi\_d7136} | Stage 13 deck rootname |
| \texttt{--rootname14} | \texttt{cai\_d7137} | Stage 14 deck rootname |
| \texttt{--vend} | $6.7$ mm/s | End-platen velocity |
| \texttt{--t-ramp} | $75$ ms | Loading-ramp duration |
| \texttt{--d-target} | $1.0$ mm | Target end shortening (over-shoot for safety) |
| \texttt{--nthreads} | $8$ | OpenRadioss thread count |
| \texttt{--workdir} | \texttt{out/} | Lima-shared working directory |
| \texttt{--ke-ie-thresh} | $0.05$ | Quasi-static KE/IE threshold |

### 9.3 Outputs from runner.py

| File | Format | Contents |
|---|---|---|
| \texttt{out/cai\_d7137\_0001.anim} | OpenRadioss animation | Per-cycle field output |
| \texttt{out/cai\_d7137\_0001.sta} | OpenRadioss state | Final state for any further chain |
| \texttt{out/T01} | OpenRadioss time-history | Reaction-force, energy histories |
| \texttt{out/cai\_load\_disp.csv} | CSV | $F(t)$, $u_x^{\text{platen}}(t)$ |
| \texttt{out/cai\_energy.csv} | CSV | KE$(t)$, IE$(t)$, ratio |
| \texttt{out/cai\_summary.json} | JSON | $F_{\text{peak}}$, $\sigma_{\text{CAI}}$, pass/fail, audit metadata |
| \texttt{out/cai\_d7137.vtkhdf} | VTKHDF | For ParaView / PyVista |

---

## 10. Risks and unknowns

1. State-file restart card-name drift. The exact spelling of the
   \texttt{/INIBRI/...} sub-cards (\texttt{STRS\_F}, \texttt{STRA\_F},
   \texttt{EPSP}, \texttt{AUX}, \texttt{THICK}) is taken from the
   Spring-back tutorial \cite{OpenRadiossConfluenceSpringback2024} and
   from the audit cross-cut D \cite{OpenRadiossEndToEndAudit2026}. The
   Altair Radioss reference manual does have an \texttt{/INIBRI}
   landing page with the same spelling; verify on the actual install
   before running. If the card spelling is different in the build
   we link, the runner emits a clear error from the starter log.

2. \texttt{/INIBRI/AUX} slot ordering. The auxiliary-variable layout
   is law-specific. For \texttt{/MAT/LAW25} this is documented in the
   LAW25 reference page \cite{AltairRadiossLAW25} but the precise
   slot ordering between Hashin-fiber-tension, Hashin-fiber-compression,
   Hashin-matrix-tension, Hashin-matrix-compression, plastic strain,
   and erosion flag is not in the same place. If the slot ordering is
   wrong, the restart silently corrupts the damage state. Validate
   by re-running Stage 13 to its end with a clean state-file write,
   immediately running Stage 14 with zero applied load, and comparing
   the Stage 14 starter-instant stress and damage fields with the
   Stage 13 final-cycle fields. They must agree to machine precision.

3. Mesh-continuity drift. Any change to the Stage 13 deck (mesh
   refinement, ID renumber, deck reformat with \texttt{inp2rad}) will
   silently corrupt the restart. The runner pins Stage 13's
   \texttt{rootname13.rad} and its include files by content hash and
   refuses to run if Stage 13's deck has changed since the
   \texttt{.sta} was written.

4. Eroded element bookkeeping. Stage 13 produces some number of
   eroded brick elements (Hashin failure with \texttt{Ifail = 2}).
   The eroded-element flag is one of the auxiliary variables. On the
   Stage 14 side, the starter must accept eroded elements as deleted
   from the kinematic mesh before \texttt{/INIVEL} reads node
   velocities. This requires the \texttt{/INIBRI/AUX} read to come
   before any contact interface re-definition. The Spring-back
   tutorial \cite{OpenRadiossConfluenceSpringback2024} uses this
   ordering; honor it.

5. Quasi-static velocity tuning. The first ramp at $v_{\text{end}} =
   6.7$ mm/s might fail the KE/IE $\le 0.05$ check, in which case the
   runner halves and re-runs. Each re-run is approximately the full
   75 ms simulation time. Budget several restarts.

6. Buckling mode misfire. If the knife-edge \texttt{/BCS} restraint
   on the long sides is too lax (covering too few node lines), the
   coupon globally Euler-buckles before sublaminate buckling at the
   impact site dominates, and the predicted $\sigma_{\text{CAI}}$ is
   far below the Lopes-Camanho band. The fix is a denser
   knife-edge node set covering the full $L_{\text{free}}$ free
   length on both long sides.

7. Cohesive-layer state mapping. \texttt{/MAT/LAW83} or
   \texttt{/MAT/LAW117} cohesive bricks store mode-I and mode-II
   damage variables. \texttt{/INIBRI/AUX} must map these correctly.
   If the Stage 13 deck used \texttt{LAW83} and the Stage 14 deck
   resolves to \texttt{LAW117} (or vice versa) the slot mapping
   silently corrupts. The runner pins the cohesive-law card across
   the two stages.

8. Single combined LVI+CAI fallback. If the restart-format quirks of
   risks 1, 2, or 7 cannot be controlled, the fallback is the
   single-deck combined LVI-then-CAI run described in Section 6.6.
   The fallback is documented and supported by the audit
   \cite{OpenRadiossEndToEndAudit2026}, but doubles the per-CAI-
   parameter compute cost. Treat as a hot backup, not a primary
   path.

9. Material card identity across stages. Even a one-decimal-place
   change in any of the LAW25 stiffness or strength values between
   Stage 13 and Stage 14 corrupts the restart. The runner enforces
   identity by sharing a single \texttt{common\_mat.inc} file
   between both decks.

10. ASTM D7137 metric-vs-imperial dimension drift. The standard is
    written in $4 \times 6$ inch units; the master plan brief uses
    $100 \times 150$ mm. The two-percent dimensional drift is below
    the standard's geometric tolerance \cite{ASTM_D7137} and is
    ignored, but the gross-area $A_g$ used in the residual strength
    calculation is the metric value, not the imperial one. The
    runner outputs $A_g$ explicitly in its summary JSON to make
    this convention auditable.

---

## References (BibTeX keys)

The keys below resolve in the project's combined bibliography
(\texttt{references/*.bib}).

- \texttt{ASTM\_D7137} -- ASTM D7137/D7137M-17 \cite{ASTM_D7137} (CAI standard)
- \texttt{ASTM\_D7136} -- ASTM D7136/D7136M-20 \cite{ASTM_D7136} (LVI standard, upstream)
- \texttt{Lopes2009LVIPart2} -- Lopes and Camanho 2009 Part II \cite{Lopes2009LVIPart2} (canonical CAI-after-LVI numerical benchmark)
- \texttt{SoutisCurtis1996} -- Soutis and Curtis 1996 \cite{SoutisCurtis1996} (closed-form sublaminate buckling)
- \texttt{Kodagali2023MesoArchitectured} -- Kodagali 2023 PhD \cite{Kodagali2023MesoArchitectured} (UofSC material card and laminate)
- \texttt{Kodagali2024LowVelocityImpact} -- Kodagali et al. 2024 \cite{Kodagali2024LowVelocityImpact} (UofSC AP-PLY LVI+CAI cross-check)
- \texttt{OpenRadiossConfluenceSpringback2024} -- Spring-back restart tutorial \cite{OpenRadiossConfluenceSpringback2024}
- \texttt{OpenRadiossEndToEndAudit2026} -- in-project audit \cite{OpenRadiossEndToEndAudit2026}
- \texttt{OpenRadiossInstall2026} -- INSTALL.md / HOWTO.md \cite{OpenRadiossInstall2026}
- \texttt{OpenRadiossDiscussion2125} -- Lima + Apptainer macOS path \cite{OpenRadiossDiscussion2125}
- \texttt{AltairRadiossLAW25} -- LAW25 (CRASURV) reference \cite{AltairRadiossLAW25}
- \texttt{AltairRadiossFailHashin} -- /FAIL/HASHIN reference \cite{AltairRadiossFailHashin}
- \texttt{AltairRadiossPropType14} -- /PROP/TYPE14 solid property \cite{AltairRadiossPropType14}
- \texttt{AltairRadiossCompositeIntro} -- composite material introduction \cite{AltairRadiossCompositeIntro}
- \texttt{KitwareOpenRadiossVTKHDF2024} -- VTKHDF converter \cite{KitwareOpenRadiossVTKHDF2024}
- \texttt{VortexRadioss2024} -- Vortex-Radioss Python reader \cite{VortexRadioss2024}
