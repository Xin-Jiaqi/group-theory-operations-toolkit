# G2: The same layer group can forbid every macroscopic \(E_\parallel\to S_z\) tensor yet allow local spin patterns — a Wyckoff-orbit test

**2026-10-22. Scientific status: SOURCE-GROUNDED GROUP-ALGEBRA PASS; REAL MAGNET/FE PHYSICS NOT ESTABLISHED; NOVELTY NOT ESTABLISHED.**  
**Executable:** [p3_site_orbit_addressability_gate.py](p3_site_orbit_addressability_gate.py) (NumPy only; tested local independent copy, deterministic assertions). Code is in **Draft PR #32**, not merged. No ab initio job, no spin-dynamics/relaxation, no guaranteed Kubo coefficient, no new magnetic layer-group catalog.

## 1. Why these explicit site orbits, not random character tables

The [AFLOW two-dimensional crystallography Wyckoff table](https://aflowlib.org/prototype-encyclopedia/Tutorials/Two_Dimensional_Periodic_Systems_Part4/) lists plane group **p3 #13** positions:
- General **\(3d\)**: \((x,y),(-y,x-y),(-x+y,-x)\);
- Special **\(1a\)**: \((0,0)\), **\(1b\)**: \((1/3,2/3)\), **\(1c\)**: \((2/3,1/3)\).

The corresponding crystallographic **layer group is \(p3\) #65** in [International Tables Vol. E](https://onlinelibrary.wiley.com/iucr/itc/Ea/ch4o1v0001/sgtable4o1o065/). For 2D hexagonal lattice vectors \(\mathbf a_1=(1,0),\mathbf a_2=(-1/2,\sqrt3/2)\), its threefold operation is
\[
R^{\mathrm{frac}}_3=\begin{pmatrix}0&-1\\1&-1\end{pmatrix},\qquad
R^{\mathrm{cart}}_3=\begin{pmatrix}-1/2&-\sqrt3/2\\\sqrt3/2&-1/2\end{pmatrix}.
\]
Directly apply \(R_3\) to each atomic position, match modulo **integer lattice translations** and construct the site-permutation matrix \(P(g)\). For an arbitrary reference \((0.137,0.224)\in3d\),
\[
P_{3d}(C_3)=\begin{pmatrix}0&0&1\\1&0&0\\0&1&0\end{pmatrix},\qquad
P_{1a+1b+1c}(C_3)=I_3.
\]
The latter sites are each fixed **modulo a lattice translation**, not as identical points in Cartesian space.

Both constructions have six *synthetic* axial local-spin sites (three top and three bottom), with a collinear layer-A-type order \((+,+,+;-,-,-)\hat z\) preserving \(C_3\). We stipulate a valid **only-\(C_3\)** decorated structural/magnetic environment; extra inversion, mirrors, half translations or antiunitary operations are **not** asserted. In particular, three chemically identical special-center sites with identical decorations might possess an extra smaller lattice translation: such an accidental enlarged symmetry must be eliminated by physically inequivalent ligand/site environments before these are called actual \(p3\) material realizations. This is **not a MnPS3 source model**.

## 2. Physical local observable vs uniform electric field

Input \(V_E=\mathbb R^2\) is a **uniform in-plane polar field** \((E_x,E_y)\) transforming as \(R_3\). Output is local, site-resolved *axial* electronic \(S^z_i\), whose axial \(z\) component is preserved by the **proper** \(C_3\) rotation and whose six site labels transform by \(P_T\oplus P_B\). The symmetry-allowed local response matrix \(K\in\mathbb R^{6\times2}\) obeys
\[
Q_z(C_3)K=K R_3,\qquad
Q_z=\operatorname{diag}(P(C_3),P(C_3)).
\]

Use two **independent** checks:
1. Build the full real nullspace of \(QK-KR\) by singular value decomposition.
2. Evaluate the **character inner product**
\[
\dim\mathrm{Hom}_{C_3}(V_E,V_Q)=\tfrac13\sum_{g=e,C_3,C_3^2}\chi_Q(g)\chi_E(g).
\]
For a generic \(3d\) layer, \(\chi_P=(3,0,0)\). For three fixed-center sites \(1a+1b+1c\), \(\chi_P=(3,3,3)\). The in-plane electric vector has \(\chi_E=(2,-1,-1)\).

### Verifiable table, with two layers in the specified model

| Output representation | Generic \(3d\) per layer | Three fixed special positions per layer |
| --- | ---: | ---: |
| Six local longitudinal \(S^z_i\) observables | **4** independent \(E_{x,y}\)-to-site coefficients | **0** |
| Three-component layer-odd longitudinal \((S^z_{Ti}-S^z_{Bi})/2\) pattern | **2** | **0** |
| Six site-local **transverse** \((S^x_i,S^y_i)\), with proper axial spin rotation included | **8** | **12** |
| One-layer three-site transverse channel | **4** | **6** |
| Global \(\sum_i S_i^z\) or summed layer-odd \(S^z\) | **0** | **0** |

The transverse entry is key: \(Q_\perp=P_{\mathrm{site}}\otimes R_{3}\), **not** just the site-permutation matrix. It is a representation of potential local transverse magnetic fluctuations, **NOT** a dynamically diagonalized magnon or Néel collective eigenmode.

Every symmetry-allowed \(3d\) longitudinal response has vanishing total-site sum and separate vanishing sum within each layer; hence normal **macroscopic \(S_z\leftarrow E_\parallel\)** and bulk layer-odd \(S_z\) are both **forbidden** even when a **local compensated site spin-density texture** is allowed. One explicit normalized allowed pattern from a basis matrix under \(E_x\) is approximately
\[
\delta\langle S^z_{T,i}\rangle\propto (-1.00000,\ 0.44868,\ 0.55132),
\]
whose sum is zero. These values are **arbitrary normalized representation-basis coefficients**, not physical signal strengths; the exact ratios are not forced by \(C_3\), only the allowed two-dimensional subspace.

**Important physics limitation:** A longitudinal \(\delta S_z\) is NOT automatically a driven *transverse magnon normal mode*. The transverse subspace dimension is separately tabulated, but physical magnetic-mode frequencies and oscillator strengths require an accepted exchange/an\-isotropy dynamical matrix and full electronic Kubo vertices; these are not in the experiment.

## 3. A stringent no-go: this is NOT an ordinary ±FE switching demonstration

The two crystal/magnetic decorations have **the same point/layer-group type** but occupy **different Wyckoff orbits** (\(3d\) versus \(1a\oplus1b\oplus1c\)). They are **not** known or computed to be equal-energy, opposite-polarization, or connected by a parent **full SOC/magnetic** operation. Their local-response Hom spaces need not have the same dimension precisely because we have **not** established symmetry-equivalent daughter states.

Indeed, a pure parent symmetry element conjugates site stabilizers and preserves their orbit multiplicities. It **cannot** map an orbit of three symmetry-equivalent \(3d\) sites into three individually fixed \(1a/1b/1c\) sites as an exact symmetry-related FE partner. Any actual transition between these decorations involves a more substantial site-orbit reorganization and may break \(C_3\) along the path. This defeats the attractive but invalid story "**ordinary FE switching changes 4 local response channels to zero**." A site-orbit change **can** be physically relevant across genuinely inequivalent structural/metastable phases; its existence and electric switchability require separate physical evidence.

**Related negative control:** merely relabeling the three \(3d\) sites with a permutation \(V\) sends \(Q\to VQV^{-1}\) and \(K\to VK\). Nullspace dimension and physical local response do **not** change. The script asserts this.

## 4. Close established predecessors: no premature priority claim

- [Watanabe & Yanase, *Group-theoretical classification of multipole order: Emergent responses and candidate materials*, PRB **98**, 245129 (2018)](https://doi.org/10.1103/PhysRevB.98.245129): **site/sublattice magnetic multipoles and electric response** already classified.
- [*Hidden physical effects in noncentrosymmetric crystals*, PRB **102**, 165143 (2020)](https://doi.org/10.1103/PhysRevB.102.165143): globally forbidden, local-sector-allowed response symmetry already systematic.
- [Tan et al., *Hidden radial spin polarization mediated by chiral sublattice symmetry*, PRB **113**, 085115 (2026)](https://doi.org/10.1103/26wm-zhv5): hidden site/sublattice spin textures, symmetry protection and real-material calculations.
- [Li et al., *Generation and Annihilation of Local Magnetic Moments by an Electric Field*, Nano Lett. **25**, 12762 (2025)](https://doi.org/10.1021/acs.nanolett.5c03407): electrically switched trimer local magnetic moments via breathing mode, a direct warning against claiming switched 3-sublattice magnetic readout.
- [Urru & Spaldin, *Magnetic octupole tensor decomposition and second-order magnetoelectric effect*, Annals of Physics **447**, 168964 (2022)](https://doi.org/10.1016/j.aop.2022.168964): locally active, macroscopically compensated magnetic responses and multipoles.

**Consequence:** the explicit \(p3\) Hom/character result is **standard induced-representation physics**, not a publishable new effect on its own. Its value is a genuine software/physics **acceptance gate for full layer-group theory**, distinguishing site-resolved from merely point-group permitted global responses.

## 5. Highest-value NEXT gate (G3), before any full 80/528 sweep

1. Find an actual **source-supported structural switching pair** that is *not a full relativistic magnetic-domain symmetry pair*, **keeping the same atomic inventory**, and extract its ordered site permutations; or prove a no-go for a proposed pair. A hypothetical \(3d\to 1a+1b+1c\) FE partner is NOT an acceptable shortcut.
2. Calculate **actual transverse magnetic-mode representations** in both states using symmetry-compatible exchange/anisotropy matrices and compare physical electric-drive modes, not just local \(S_z\) density labels.
3. For antiunitary relations, perform the **source/measurement operator-swapped Onsager** test rather than same-index AC parity. Check contact terms for an actual field-coupled Kubo computation.
4. Keep **one exact complete-domain partner** and **one pure site-relabel** as negative controls, so any new claimed channel is clearly intrinsic rather than a coordinate/sublattice gauge artifact.
5. Compare direct mechanisms with 2018 PRB multipoles, 2020 hidden physical effects, 2025 trimer electro-magnetism, 2026 HSP and 2026 magnon/SOT group theories **at the same observable level**. Stop if the alleged effect is merely another hidden-spin site multipole.

**Scientific decision:** G2 **PASS as a precise mathematics/structure-vs-site-orbit demonstration**, but **FAIL as proof of FE-controlled new magnetic physics**. Extension to exhaustive groups without a genuine switching source pair, a physical mode and a distinct measurable response would mostly replicate existing representation theory.
