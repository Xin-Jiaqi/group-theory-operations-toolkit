# FE domain switching and magnetic-mode *addressability*: a symmetry no-go before any 80/528 catalog

**22 October 2026. EXPERIMENTAL GROUP-THEORY RESEARCH; not a production feature, a full layer-group classification, or a new-physics claim.**  
**Execution:** \`python research/fe_domain_mode_addressability_20261022/domain_covariance_gate.py\` (NumPy only). The standalone locally executed finite-dimensional gates pass. This branch is intended as a **draft PR**, so repository review and CI decide integration.

## 1. Why another complete symmetry table is not sufficient

Important verified earlier work:

1. Zhu et al., [*Sliding Ferroelectric Control of Unconventional Magnetism in Stacked Bilayers*, PRL 135, 056801 (2025)](https://doi.org/10.1103/dmzg-ck2t), already **systematically survey 80 layer groups** using magnetic structural connection operations and switched spin splitting/AHE.
2. Du et al., [*Stacking theory for bilayer two-dimensional magnets*, PRB 113, 214455 (2026)](https://doi.org/10.1103/g471-xrkl) enumerate **448 collinear spin layer groups** and solve a general magnetic stacking theory. **Do not confuse these with the 528 magnetic layer-group types** in the crystallographic Litvin tables; they are different symmetry categories.
3. Sun et al., [*Unified symmetry framework for spin–ferroelectric coupling*, Nat. Commun. 17, 3101 (2026)](https://doi.org/10.1038/s41467-026-69635-2) classify three FE-spin-band switching types (including real bilayer MnPS3 magnetic connections).
4. Chen et al., [*Unconventional magnons in collinear magnets dictated by spin space groups*, Nature 640, 349 (2025)](https://doi.org/10.1038/s41586-025-08715-7) comprehensively classify unconventional spin-space-group magnons. [Liu et al., *Oriented spin space groups*, Nature 652 (2026)](https://doi.org/10.1038/s41586-026-10401-1) and [SSG projective irrep theory, arXiv:2609.21289](https://arxiv.org/abs/2609.21289) block a general claim of newly introducing magnetic spin-group representations.
5. Xie et al., [*General Theory of Chiral Splitting of Magnons in Two-Dimensional Magnets*, arXiv:2601.15031](https://arxiv.org/abs/2601.15031) already propose an electric-field chiral magnon classification using **464 collinear spin layer groups** (PREPRINT, category conventions may differ from Du's 448; **do not conflate counts or assert one is wrong**). Ni et al., [*Electric-Switchable Chiral Magnons in PT-Symmetric Antiferromagnets*, arXiv:2601.16052](https://arxiv.org/abs/2601.16052) identify magnetic layer groups permitting electric magnon splitting (PREPRINT).
6. Liu et al., [*Symmetry-Dictated Switching of Antiferromagnetic Magnon Transport in 2D Multiferroics*, Nano Lett. 26, 11017 (2026)](https://doi.org/10.1021/acs.nanolett.6c02164) already use FE to switch magnon Berry curvature and thermal Hall response via exchange/DMI asymmetry.
7. [Tang & Wan, PRB 105, 155156 (2022)](https://doi.org/10.1103/PhysRevB.105.155156) classify band nodes of **528 magnetic layer groups**; the [IUCr/Litvin tables](https://www.iucr.org/what-we-do/publications/books/Magnetic-Group-Tables) are the established symmetry data.

**Verdict:** an allowed/forbidden magnetic response tensor table for 80/528 groups, generic FE chirality switching or the usual magnon irrep catalog is not a sufficient independently new paper. We seek a **testable distinction between physically inequivalent switching and laboratory-coordinate relabeling**.

## 2. Exact domain-covariance no-go (unitary full-Hamiltonian partners)

Let the real in-plane electric-field space \(V_E\) carry an orthogonal representation \(D_E\), and normalized *actual magnetic-mode coordinate* space \(V_Q\) carry an orthogonal representation \(D_Q\). A static real linear generalized-force coupling \(K:V_E\to V_Q\) is allowed by an **individual state's** stabilizer \(H\) if

\[
D_Q(h) K = K D_E(h),\quad h\in H,
\qquad \dim\mathrm{Hom}_H(V_E,V_Q)
\]
is its number of independent symmetry-allowed coefficients.

**Conditional theorem (ordinary exact unitary domain relation only).** If a parent **full magnetic Hamiltonian**, including the magnetic order and SOC whenever present, exactly relates domain \(A\) to \(B\) by \(c\), then in properly mapped orthonormal mode and field coordinates

\[
K_B = D_Q(c)K_A D_E(c)^{-1},\qquad H_B=cH_Ac^{-1}.
\]

Hence

\[
\operatorname{rank}K_B=\operatorname{rank}K_A,\quad
s_i(K_B)=s_i(K_A),\quad
\dim\mathrm{Hom}_{H_B}=\dim\mathrm{Hom}_{H_A}.
\]

These are **elementary consequences of orthogonal/unitary similarity of domains**, not a newly discovered physical law. They give a powerful research **NO-GO**: a fixed-laboratory polarization can turn a particular *component* OFF→ON upon rotating the domain, but this cannot be advertised as creation of a genuinely independent magnetic excitation channel unless the actual mode, lab geometry and complete state transformation are treated correctly.

**Important limits:** this is a **static real/uni­tary statement**. An antiunitary map can mix complex conjugation, time-reversal parity, domain magnetization and retarded/advanced causal prescriptions. Do **not** apply the above finite-frequency equality or singular-value argument automatically to a generic retarded \(\chi(\omega)\) under antiunitary partner switching. The source MnPS3 Eq.13 \([E\parallel M_zM_y]\) is a **nonrelativistic spin-space** operation and by itself does not imply a full **SOC** unitary connection.

### A concrete falsification of "FE makes a new bright mode"

For a mirror-even scalar magnetic mode and mirror \(m_x\), a polar in-plane electric field has \(D_E(m_x)=\mathrm{diag}(-1,+1)\). The invariant response is \(K_A=(0,\kappa)\): \(E_x\) dark, \(E_y\) bright. Rotate both complete states by \(C_{4z}\), with \(D_E(C_{4z})=\begin{pmatrix}0&-1\\1&0\end{pmatrix}\). The partner has \(K_B=(-\kappa,0)\): a fixed-lab \(E_x\) appears bright, yet \(s_1(K_A)=s_1(K_B)=|\kappa|\), both have rank one, and the field rotated together with the domain remains equally bright.

The standalone [\`domain_covariance_gate.py\`](domain_covariance_gate.py) evaluates \(\mathrm{Hom}\) dimensions, exact partner stabilizer covariance, two-dimensional rank/singular-value invariance, and this fixed-lab OFF→ON counterexample. All are **synthetic representation matrices**, NOT full-layer-group or MnPS3 site results.

## 3. Three genuinely different symmetry questions, ordered by scientific value

**Class 0 — symmetry-trivial partner:** complete SOC magnetic states are exact unitary partners. The physical modal addressability pattern changes only by covariance under corresponding coordinates, although a fixed device's lab-axis selection can still be useful engineering. **Not intrinsic-new-channel physics**.

**Class 1 — structural partner but NOT full magnetic partner:** stacking operation connects ionic structures or *nonrelativistic* spin-space models but does **not** map the full SOC magnetic Hamiltonian (true magnetization axes and local-moment site map fail). The null statement above does not apply. Symmetry permits a different **physical** mode-driving matrix, but neither its magnitude nor actual nonzero value follows from symmetry. The detailed SOC obstruction must be tested, not assumed to create an interesting effect.

**Class 2 — nonconjugate magnetic isotropy / distinct physical order:** different magnetic stabilizers, occupied site orbits, magnetic orientation/texture or magnetic-mode representations. The allowed map space \(\mathrm{Hom}_H(V_E,V_Q)\) can change dimension, providing a **robust group-theoretic difference** in the number of independent allowed couplings. However this does **not** prove FE energetic degeneracy, switchability, finite oscillator strength or a detectable AFM resonance.

**Potential route:** classify **switching orbits and genuine connection classes** modulo source and destination stabilizers (candidate double-coset \(H_B\backslash G/H_A\) when both embed consistently in one parent group), followed by source-correct site permutation and physical-mode induced representations. A double-coset labels *symmetry equivalence classes of operations*, not energetic switching paths. This generic group fact alone is **not an originality claim**.

A no-material first benchmark is a **state pair with an actual stabilizer change** vs the symmetry-trivial partner in the test above. Restrict initially to **unitary SOC groups**, where causal antiunitary response subtleties cannot be swept under the rug.

## 4. Where full **layer** groups (rather than point groups) become essential

For a truly **uniform** field and \(\Gamma\)-point excitations, many factorized translations disappear from the external-field representation; ordinary bulk tensor dimensions often reduce to the magnetic point/co-group. The full layer-group information is still essential for **actual Wyckoff/site permutations, layer exchange and induced local-mode representations**, which the point-group tensor table misses.

At **finite in-plane wave vector** \(q\), nonsymmorphic \((R|\mathbf t)\) acts with Bloch translation phase, and a generic point-group table is demonstrably inadequate. For a glide \(g\) with \(g^2=T_a\),

\[
D_q(g)^2=e^{-iq a}\,I,
\]
so at \(q=\pi/a\) the glide eigenvalues are \(\pm i\), not the naive mirror \(\pm1\). The test script verifies this for a two-site exchange representation.

This is **standard nonsymmorphic projective representation theory**. Its possible research use is systematic **electric/phonon/near-field addressability** of site- or layer-resolved magnons at selected \(q\), in the presence of FE switching. A spatially uniform free-space electric drive generally has \(q\approx0\) and CANNOT by itself excite a zone-edge magnetic mode in a translationally invariant crystal. A plausible finite-\(q\) experiment requires a justified patterned field, phonon coupling, near-field probe, defect momentum mixing, or supercell folding.

## 5. What to build from the existing toolkit, and what not to promise

The existing \`group-theory-operations-toolkit\` already stores **80 ordinary LG, 528 magnetic LG, 122 MPGs**, standard Seitz operations and FE stacking classes, ordinary LEE and multiple response-tensor invariant bases. This research branch **does not modify production APIs or magnetic-group tables**.

**Phase A (low cost, ~pure group algebra):**
- Freeze two actual magnetic states and their site/mode bases; obtain individual stabilizers, physical magnetic connecting operations, and the FE connection graph from structure/origin-corrected Seitz operations.
- Compute and compare \(\mathrm{Hom}_{H_A}(V_E,V_Q)\) vs \(\mathrm{Hom}_{H_B}(V_E,V_Q)\), along with **covariance versus nonconjugacy**; reject any claimed "new modes" that result from rotated coordinates.
- Use at least **one nonconjugate** example with a source-validated ordered Wyckoff orbit and one symmetry-equivalent negative control. Enumerating group names without site-orbit/magnetic-mode data does not generate an honest all-groups magnon catalog.

**Phase B (only if A yields an authentic discriminator):**
- Add finite-\(q\) little groups and nonsymmorphic/projective phases; compare the few cases where phase factors change the actual *field-to-mode* selection rather than only magnon band connectivity.
- Implement antiunitary co-representation and causal-frequency transformation carefully; no naïve time-odd sign assignment to complex retarded susceptibilities.
- Extend from representative site orbit/stacking cases to a **properly evidence-qualified** 80 LG/528 MLG mode-addressability table, using magnetic site/projector and mode info and clear no-go/unknown labels.

**Phase C (small real model/DFT only if earlier phases survive peer comparison):** construct one observable with a genuine lab-feasible drive and an experimentally distinguishable spectral signature. Then assign a single real material. No material workflow launches from this note.

**No-go to stop entirely:** if all exact FE partners are mere coordinate-covariant tensor changes; all apparent channel additions are identical to published point-group/SOT predictions; or finite-q field/mode coupling has no feasible experimental momentum source, **archive a rigorous negative result rather than create a new tensor name**.

## 6. Paper-level bar and active compute boundary

A publishable *general* theory would need a **nontrivial new physically testable classification invariant or selection principle**, exhaustive reproducible evidence for the stated class, correct magnetic SOC vs nonrelativistic SSG and antiunitary handling, plus an observable beyond existing PRL stacking theory, 2025 Nature SSG magnon catalogs, 2026 PRB spin-layer-group stacking theory and 2026 field/FE-magnon classifications.

**Current verdict:** the verified **domain-covariance no-go is useful but textbook group theory**. The group-level scientific opening is **not yet proved**; Class 1 or Class 2 physics, possibly finite \(q\), is a research search strategy subject to strict novelty checks.

No MnPS3 DFT jobs or R1 restart authorized. The original Xiang PRL reproduction is untouched. For real magnetic site mapping, see the previously approved read-only B+D R2-P0 task in \`quantum-geometry-reproductions\`; full magnetic mode projectors have not been returned yet.


### Important additional publication status: the full spin-layer-group *representation theory* is also now available as a preprint

[Zhang et al., *Spin layer groups and their corepresentations*, arXiv:2605.25484v2 (29 May 2026)](https://arxiv.org/html/2605.25484v2) already tabulates **40,087** inequivalent spin layer groups under its **supercell-index cutoff \(I_k\le12\)**: **448 collinear, 6,083 coplanar-noncollinear, 33,556 noncoplanar**, and constructs little-group corepresentations. This is a PREPRINT, not verified here as a journal publication. It makes the suggestion "we will newly classify all magnetic/spin layer-group irreps" particularly weak. Our possible distinction must be an **actual electric-drive-to-site/mode addressability *functional* under FE switching**, with a physically available excitation probe, not a spin layer-group table or representation construction. The **464 count** used in [arXiv:2601.15031](https://arxiv.org/abs/2601.15031) is a different reported enumeration that requires definition comparison; this note does **not** resolve the discrepancy.


## 7. NEW: independent antiunitary finite-frequency response gate (2026-10-22)

**Correction to any naive extension of Sec. 2:** The unitary full-domain rule \(K_B=D_QKD_E^{-1}\) is a legitimate static representation-covariance result. For an **ANTIUNITARY** domain connection \(\Theta=UK\) the correct finite-frequency retarded relation is *not*, in general, same-index \(\chi_{A B}^{(B)}(\omega)=\eta_A\eta_B\chi_{A B}^{(A)}(\omega)\) or \(\eta_A\eta_B\chi_{A B}^{(A)}(\omega)^*\).

With *strict definitions*: the partner Hamiltonian \(H_B=U H_A^*U^\dagger\), two Hermitian measurement/source operators \(A,B\) of definite antiunitary parities \(\Theta A\Theta^{-1}=\eta_A A\), \(\Theta B\Theta^{-1}=\eta_B B\) (\(\eta=\pm1\)); the Lehmann-retarded response with \(e^{-i\omega t}\) and identical causal \(+i\eta\) **satisfies the Onsager-Casimir cross response**
\[
\boxed{\chi_{A B}^{(B)}(\omega)=
\eta_A\eta_B\chi_{B A}^{(A)}(\omega).}
\]
For tensor components and nontrivial group rotations, transform both physical operator bases and their indices appropriately. This is a known **generalized Onsager relation**, NOT a discovery or a new classification invariant. It does **not** state that the *forward* cross susceptibility for the two magnetic domains is the negative of the same forward cross susceptibility in one domain. At finite \(\omega\) the forward and reciprocal functions can differ; broadening/contact protocols must match for both.

[**Independent two-band and physical \(\Theta^2=-1\) four-band Lehmann test**](antiunitary_kubo_onsager_gate.py) was separately executed on the same complete spectral Kubo sum and demonstrates:
- spinless \(\Theta=K\), \(H_A=0.5\sigma_x+0.8\sigma_y+0.7\sigma_z\), \(A=\sigma_z,B=\sigma_x\), at \(\omega=.6,\eta=.03\): \(\chi_{AB}^A=0.22188817+0.15989348i\) whereas \(\chi_{AB}^B=\chi_{BA}^A=0.23993095-0.15667204i\).
- spinful \(\Theta=(I_{\rm orbital}\otimes i\sigma_y)K\), \(A=I\otimes\sigma_z/2\) (odd), \(B=\tau_x\otimes I\) (even): at \(\omega=.24\), \(\chi_{AB}^A=0.03240507-0.00152205i\), but \(\chi_{AB}^B=-\chi_{BA}^A=-0.03199982-0.00165070i\).
- Reciprocity residuals \(\sim10^{-17}\); **naive same-index sign fails** by nonzero finite-\(\omega\) difference. The example has random Hermitian \(4\times4\) matrices as a purely algebraic verification (not a material, and no causality contact term added).

**Consequences for group-wide G1:** distinguish (a) exact unitary domain covariance of **same type of response**, (b) antiunitary **Onsager-swapped** cross correlations with appropriate time parity, (c) structural or nonrelativistic-spin connections that are not symmetries of the *relativistic* SOC Hamiltonian, and (d) different stabilizers/site orbits. Full antiunitary corepresentation catalog and experimentally normalized electric-field coupling are still OPEN. The primary source MnPS3 path1 Eq.13 \([E\parallel M_zM_y]\) is a **nonrelativistic spin-space** statement, not by itself a valid physical \(\Theta\) for this Kubo identity.


## 8. Falsification-first gate: claim-level NO-GO does NOT delete scientific assets

The project now has a public [B+D idea elimination protocol and 20-claim register](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/blob/main/papers/xiang-2025-spin-quantum-geometry/research/brainstorm-20261021/docs/BD_IDEA_ELIMINATION_GATES_20261022.md) and a [20-entry machine-readable registry](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/blob/main/papers/xiang-2025-spin-quantum-geometry/research/brainstorm-20261021/data/BD_IDEA_REGISTRY_20261022.json).

Research-only [idea_elimination_gate.py](idea_elimination_gate.py) provides a conservative claim-screening utility (Python stdlib only; self-test mode supplied) with [NO-GO example](example_BD08_NO_GO.json) and [OPEN but currently BLOCKED example](example_BD12_BLOCKED.json). Its possible statuses are NO_GO, BLOCKED, PHYSICS_CANDIDATE and **none** is a novelty certificate. It cannot certify that user-filled FE-state mappings, papers, same-gauge vertices, or original material outputs are true. This is an **unmerged PR branch**, not a released group-toolkit API; do not silently invoke it in production.

Closest newly identified **direct precedence** is particularly strong:
- [González-Hernández et al., *Non-relativistic torque and Edelstein effect in non-collinear magnets*, Nat Commun 15, 7663 (2024)](https://doi.org/10.1038/s41467-024-51565-6) already analyzes **local magnetic-sublattice projected Edelstein response and exchange torque**, including insulating systems, with spin groups. A mere site-resolved spin torque or spin-group-vs-MPG classification is not new.
- [*Local Magnetoelectric Effects as Predictors of Surface Magnetic Order*, PRX 15, 021094 (2025)](https://doi.org/10.1103/7brd-lynv) predicts surface magnetic modifications from local bulk magnetic multipoles; simple bulk-forbidden interface-local allowed magnetism is not new.
- [*Shift photocurrent vortices from topological polarization textures*, npj 2D Mater Appl 9, 43 (2025)](https://doi.org/10.1038/s41699-025-00552-6) derives moiré polar-vortex-induced **shift-current vortices with non-Abelian Berry geometry**; generic vortex-based optical stories are already occupied.
- [Su & Lin, *Topological sliding moiré heterostructure*, PRB 101, 041113(R) (2020)](https://doi.org/10.1103/PhysRevB.101.041113) establishes sliding synthetic-coordinate Chern/topological charge pumping. An arbitrary path-holonomy/topological sliding pump cannot be claimed without a sharply different observable.
- [Ke, Liu & Liu, *Superlubric Motion of Wavelike Domain Walls in Sliding Ferroelectrics*, PRL 135, 046201 (2025)](https://doi.org/10.1103/jhlq-2dd7) warns that **rigid uniform layer sliding is often not the experimentally relevant switching path**; wall motion and Born-charge physics matter.

**Next physical test, not a discovery:** use **actual BSF FE domain/connection evidence** and magnetic-site/projector operators to determine whether *driven real FE domain-wall switching* has a measurable compensated layer/magnetic-mode signal not predicted by already known local SOT, magnetic multipoles and wall electrodynamics. A wall-local allowed tensor, alone, fails the originality test. No VASP/HPC authority is implied.


## 9. G3 high-return decision: do not run domain-wall DFT, screen registry-gradient representations first

The [B+D source/ROI report](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/blob/main/papers/xiang-2025-spin-quantum-geometry/research/brainstorm-20261021/docs/BD_G3_DOMAIN_DFT_COST_AND_REGISTRY_GRADIENT_20261022.md) explicitly **defers magnetic FE domain-wall DFT**: the source-derived MnPS3 R1 primitive-cell SCF already failed under 2 h walltime, and larger wall supercells have more structural, magnetic, electrostatic and response-vertex degrees of freedom. [Nasir/Atkinson PRB 113, 035425 (2026)](https://doi.org/10.1103/y87v-34d8) have FE-wall Edelstein precedent; [Dong et al. RRL 2024](https://doi.org/10.1002/pssr.202400062) magnetic sliding bilayer spin/orbital Edelstein; [Yarmohammadi et al. PRB 113, 184403 (2026)](https://doi.org/10.1103/4txy-gy8f) slow-phonon Edelstein; [WTe2 coherent sliding/SHG Nature Physics online 21 Oct 2026](https://doi.org/10.1038/s41567-026-03491-1) direct sliding-mode dynamics. Plain wall, sliding-phonon and local spin response claims are **not new**.

New, limited group-theoretic [\`registry_gradient_gate.py\`](registry_gradient_gate.py) checks (for an **assumed physical representation of relative displacement U**):
\[
K_0\in\mathrm{Hom}_{H_*}(E,O),\qquad
\Gamma=\partial_\delta K\in\mathrm{Hom}_{H_*}(U\otimes E,O).
\]
For a mirror \(m_x\), axial \(S_z\) is odd, \(E_x,u_x\) odd and \(E_y,u_y\) even. \(E_y\to S_z\) is forbidden at the symmetric reference but a first-order \(u_xE_y\to S_z\) is symmetry-allowed; equally \(u_yE_x\to S_z\) is allowed. This is a **simple tensor-product/Neumann result**, not a newly discovered phonon effect or a nonadiabatic response. **If a space operation swaps top/bottom layers, U need not transform as a simple polar field:** verify explicit Seitz + layer indexing before running across 80 LG.

**G3.1 next non-DFT evidence:** use two actual 80-layer-group BSF-generated *source-supported* cases with proper layer-exchange and local-magnetic mode bases; compare rank/induced operators and candidate \(\Gamma\) to the closest already-published microscopic phonon/Edelstein response. If result is only ordinary tensor activation or standard SOT, preserve as METHOD_PASS / NO_GO_NOVELTY. Do not claim a full 80/528 census from this one mirror proof. No HPC or DFT job authorized.


## 10. Physics-first correction: direct competitors vs Introduction and a one-parameter seed (22 Oct 2026)

User-requested [core scientific reset and paper-by-claim reference roles](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/blob/main/papers/xiang-2025-spin-quantum-geometry/research/brainstorm-20261021/docs/BD_CORE_PROBLEM_REFERENCE_ROLES_AND_SMALL_SEED_20261022.md) replaces over-eager topic-wide "NO-GO" heuristics. Every NO_GO rejects ONLY the specific general originality claim, not the possibility of a different source/measurement/frequency-specific effect. A paper is a **DIRECT competitor** only relative to the SAME observable, source/experimental drive, response regime and overlapping conclusion (or a demonstrably inclusive theorem). Related methods and materials are legitimate **Introduction context**, not automatic project rejection.

Two standalone, **locally executed and corrected** symmetry gates now give an actually focused small next calculation:
- [physical_layer_exchange_gate.py](physical_layer_exchange_gate.py): for layer-swapping operations, \(\mathbf u=\mathbf r_T-\mathbf r_B\) transforms as \(D_u(g)=-R_\parallel(g)\) and layer-odd spin adds a permutation sign to the axial-vector representation. Under \`mz\`, naive \(\mathbf u\to R_\parallel\mathbf u\) would wrongly forbid allowed first-derivative components.
- [stacking_first_seed_gate.py](stacking_first_seed_gate.py): in an **IDEAL** bilayer with \`C3z\`, vertical \`mx\`, layer-swapping horizontal \`mz\`, a uniform in-plane charge-current-driven nonequilibrium layer-odd \(S_-^z=S_T^z-S_B^z\) has no zero-shift linear response; its first order in slide is uniquely
\[
\delta S_-^z=\gamma\,(u_x j_y-u_y j_x)+O(u^2j).
\]
A full \(12\times4\) linear constraint SVD has rank 3 and nullity 1, with the \(\pm\) normalized antisymmetric matrix as its only basis. This is **STANDARD** equivariant tensor algebra, and **does not** show \(\gamma\ne0\), a real material bilayer symmetry, FE switchability, or a novel spin accumulation. Unlike the 2026 PRL's vertical-electric-field-activated response, this **relative-shift parameter derivative** may be an analytically useful target, but the primary PRL SM must be reviewed before treating it as an independent contribution.

**Proposed first microscopic output, not yet computed:** one controlled \(\chi_{S_-^z,j_y}(u_x)\) curve in a physical SOC/transport model, with odd-slide/mirror/SOC-off and fixed-current-vs-fixed-field checks. Use the existing 2026 LEE source's minimal model *only if* it actually permits the required registry deformation and layer-operator embedding. No new VASP/HPC jobs.


## 11. Test against the *published* LEE minimal Hamiltonian (22 Oct 2026)

Added [\`lee_published_model_baseline.py\`](lee_published_model_baseline.py), a dimensionless Eq. (2)-style Boltzmann sanity test of Zhou et al., PRL 137, 076401 (2026), Eq. (3). With \(t_\perp=1,\alpha=0.5,m=1,\mu=0,T=0.02\), 161 values of \(\Delta/t_\perp\) satisfy:
- layer-odd in-plane \(\chi_{S_-^y,j_x}\) even in \(\Delta\);
- global in-plane \(\chi_{S_+^y,j_x}\) odd in \(\Delta\);
- \(\chi_{S_-^z,j_y}=0\) exactly in the published minimal model, because the spin-helicity eigenstates are in-plane and \(\langle\tau_z\sigma_z\rangle=0\).

The numerical parity residuals were \(2.22\times10^{-16}\). This **is not a material fit or novel response**. Crucially, it shows why the proposed slide-induced \(S_z\) curve cannot be computed by simply changing scalar \(t_\perp\) or \(\Delta\) in Eq. (3). Need a source-backed spin-z SOC/valley term AND microscopic registry derivative \(\partial_uH\); do not add a free spinful hopping merely to force nonzero \(\gamma\). See [source model audit](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/blob/main/papers/xiang-2025-spin-quantum-geometry/research/brainstorm-20261021/docs/BD25_LEE_MODEL_AUDIT_20261022.md). No DFT/HPC jobs.


## 12. Source-backed geometry prerequisite for a registry-dependent SOC model

Added [\`aa_mos2_registry_geometry_gate.py\`](aa_mos2_registry_geometry_gate.py) using the AA bilayer MoS\(_2\) coordinates and interlayer spacing from Venkateswarlu et al., PRB 102, 081103(R) (2020). The script was run locally and passes 450 interlayer pair/lattice-translation checks for \(u\) and \(-u\); the maximum distance covariance error is \(1.776\times10^{-15}\) Å. It also recovers the source's \(3.676\) Å interfacial S–S plane spacing.

This is **geometry-only**, not a tight-binding Hamiltonian, SOC, band structure or response. The initial test correctly failed on a bug in the assertion that upper/lower S planes are the same chemical species; after changing both to chemical label S and retaining their distinct mirror indices, the test passed. Next is the full orbital Slater–Koster matrix covariance, including orbital parity/rotation, before SOC/Kubo. No DFT/HPC jobs.


## 13. BD27–BD28: a genuine *general-theory seed* from configuration-space response topology (2026-10-22)

**Research report:** [BD27 stacking-torus zero/index and parity theory](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/blob/main/papers/xiang-2025-spin-quantum-geometry/research/brainstorm-20261021/docs/BD27_REGISTRY_TORUS_ZERO_INDEX_GATE_20261022.md). The key methodological correction remains: a specific earlier material example **does not automatically pre-empt** a later family-wide theory. Any priority comparison must match the level and actual proposition; [BSF PRL 130, 146801 (2023)](https://doi.org/10.1103/PhysRevLett.130.146801) and [LEE PRL 137, 076401 (2026)](https://doi.org/10.1103/qhvs-pfkq) are **direct foundations for high-symmetry registry geometry and general layer-spin effects**, but their SM needs to be compared to **the exact response zero/index rule** before deciding novelty. The 2025 [shift-current vortex paper](https://doi.org/10.1038/s41699-025-00552-6) treats another source/output and is an Introduction/nearest-texture result for this claim.

**All tests below were executed on standalone NumPy scripts locally; no DFT or electronic current–spin susceptibility was computed:**

| Script | Physically defined symmetry comparison | Result |
|---|---|---|
| [stacking_response_torus_gate.py](stacking_response_torus_gate.py) | Hex C3 + vertical mirror + horizontal layer swap; current \(\to S_-^z\) real response field over \(\mathbb T^2\) | **Six forced zero registries** G/A/B/C/M/N; Fourier witness indices +1 at G/M/N and −1 at A/B/C; sum zero |
| [registry_fixedpoints_exact.py](registry_fixedpoints_exact.py) | **Exact rational** solve \((I-R)\delta\in\Lambda\) | C3 has 3 fixed registries; C4 has 2; layer-swap inversion has 4; unions hex 6 and square 4 |
| [stacking_square_response_torus_gate.py](stacking_square_response_torus_gate.py) | Square C4+mirror+layer swap, same current \(\to S_-^z\) | **Four forced zero registries**, witness indices +1,+1,−1,−1 |
| [registry_layer_even_spin_gate.py](registry_layer_even_spin_gate.py) | Hex current \(\to S_+^z\) instead of \(S_-^z\) | Three C3-forced zero registries; **quadratic** onset and generic G winding **−2**; M/N +1, total zero in witness |
| [layer_parity_activation_gate.py](layer_parity_activation_gate.py) | Same hex family, layer even vs odd; square C4 vs horizontal layer swap | Hex \(S_-^z\) first allowed \(i\gamma u\) with index +1; \(S_+^z\) first allowed \(\beta\bar u^2\), index −2 if nonzero. Square \(C_{2z}\) versus \(m_z\) requires **\(S_+^z\) response identically zero for all registry translations** |

A successful **method gate** is distinct from a physics discovery. Poincaré–Hopf on a smooth real vector field over the torus gives \(\sum\mathrm{ind}=0\); a simple C3-fixed zero has +1 index; if these are the *only* six zeros, half-translation A/B/C must be −1 saddles. But an actual spin response can vanish everywhere or have additional zeros, complex AC needs separate treatment, and \(\gamma,\beta\) remain entirely unknown.

**Research novelty test:** an S2 classification is worthwhile only if response-specific *operator parity, global zeros and activation/index rules* create a nontrivial, measurably distinct prediction beyond a restatement of existing BSF registry-point tables, LEE static tensor selection, and standard symmetry/index theorems. Start with one source-backed microscopic model and compare the original LEE SM equation-by-equation. The next physical task is still **not** a 528-MLG scan or domain-wall DFT.


## 14. Direct BiDB geometric validation resolves the "ideal symmetry, no real material" blocker

The user supplied the full actual BiDB ASE/SQLite dataset; [BD29 database analysis](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/blob/main/papers/xiang-2025-spin-quantum-geometry/research/brainstorm-20261021/docs/BD29_BIDB_DB_DIRECT_SCREEN_20261022.md) validated 11,184 systems and uses exact species ratios (not BiDB's abstract ABC stoichiometry). Full source database is **not public in this repository**.

**New read-only independent geometry code:** [bidb_geometry_layer_swap_gate.py](bidb_geometry_layer_swap_gate.py), invoked with the user's own database path and explicit IDs; validates source atomic species and position permuations modulo in-plane lattice under C3z, a vertical mirror and *pure horizontal layer-exchanging mz*. Executed locally with these results:

| BiDB id | Material | Structural C3 / vertical mirror / horizontal layer swap | Maximum matching residual |
| ---: | --- | --- | ---: |
| **980** | **Pd2Te2** | 8/8 PASS / 8/8 PASS / **8/8 PASS** | **0.001291 Å** |
| **975** | **Pd2Se2** | 8/8 PASS / 8/8 PASS / **8/8 PASS** | **0.001044 Å** |
| 2903 | NBa2 | all 6/6 PASS | 0 |
| 3555 | ScSe2 | all 6/6 PASS | 0 |
| 2684/2685 | WSe2 | C3, mv PASS; **mz FAIL** | fails first site ~1.916 Å |
| 8055 | VTeS | C3, mv PASS; **mz FAIL** | ~0.673 Å |
| 7964 | VSeTe | C3, mv PASS; **mz FAIL** | ~2.005 Å |

**First microscopic test family is now Pd2Te2 id980** (Stable p-6m2 / D3h, 8 atoms, BiDB nonmagnetic/metallic flags) with Pd2Se2 id975 as control, precisely because its physical layer exchange matches the original BD27 hypothesis. The Pd2Te2 parent monolayer is **0.134 eV/atom above the hull**, so it is a **theoretical demonstrator**, not a confirmed feasible material. Need Pd/Te orbital symmetry/SOC TB, fixed-current spin operator and causal response; the group construction alone proves neither nonzero spin accumulation nor a new effect.

Magnetic extension: VSeTe (8 Stable FM bilayers), VTeS (2 Stable FM), then AgVP2Se6 (17 Stable but 20-atom low-symmetry), with the critical warning that BiDB interlayer FM/AFM fields do **not** certify intralayer Néel/spiral or the SOC magnetic group.

No new DFT, HPC or external access to the source database has been initiated. No source structures or full BiDB contents are committed to GitHub; selected research notes and self-contained symmetry test code only.


## 15. First nonzero geometry-derived spin-current model after BiDB data, 2026-10-22

[BD29 microresponse report](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/blob/main/papers/xiang-2025-spin-quantum-geometry/research/brainstorm-20261021/docs/BD29_PD2TE2_TEP_SOC_MICROSCOPIC_RESPONSE_20261022.md) provides a **real-BiDB-geometry-based but UNFITTED 24-band toy**, not a Pd2Te2 material prediction. New local-data-dependent [pd2te2_te_p_soc_boltzmann_gate.py](pd2te2_te_p_soc_boltzmann_gate.py) reads **BiDB id 980** (id 975 as a composition control) via a user-supplied SQLite file path; no database bytes are committed. It selects four real Te sites × three p orbitals × two spin components; Pd-d is omitted. Real bond geometry enters SK p-p hopping, atomic \(\lambda{\bf L}\cdot{\bf S}\) supplies SOC, no arbitrary spin-dependent hopping was added. Its hopping, onsite energies, chemical potential and SOC amplitudes are chosen illustrative parameters, **not DFT fitted**.

A local equivalent implementation passed full 24×24 C3, vertical mirror and layer-swapping mirror matrix covariance (≤\(5.6\times10^{-16}\) relative error), time reversal on H and velocities, u→−u layer odd/even spin parity, atomic SOC-off and interlayer-hopping-OFF controls. A raw BiDB lattice had a ~0.003 Å exact-hex deviation, so the model **explicitly** projects to exact hex; do not pretend the raw coordinates were perfectly symmetric. Degenerate k-point traps in diagonal Boltzmann formula were rejected by shifted, full-C3v-orbit integration; the script is only a **model Boltzmann gate**, not a fully gauge-independent causal Kubo package.

**First result:** at \(u_x=0.10\) Å, \(\chi_{S_-^z,j_y}\) in normalized toy units is −0.006525 (12×12), −0.006511 (16×16), −0.006507 (20×20), −0.006506 (24×24). It reverses under u and vanishes for SOC=0; all interlayer hopping OFF drives it to ~\(10^{-16}\), confirming both SOC and interlayer hybridization are necessary in this model. Layer-even Sz current response is even in u and much smaller, but its amplitude is not converged to claim a real β or winding index. The reported coefficients have no experimental units/calibration. [Full 11-point computed CSV](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/blob/main/papers/xiang-2025-spin-quantum-geometry/research/brainstorm-20261021/data/BD29_PD2TE2_TEP_SOC_BOLTZMANN_20261022.csv).

**Science novelty unresolved:** ordinary sliding-induced Edelstein signals have prior art; the specific broad **S2/S3 response-field / global registry zero and operator parity** framework still requires a complete LEE SM equation-level comparison and a source-parametrized spinful material Hamiltonian. A specific published material example cannot by itself defeat a broader theory, but standard group covariance cannot be rebranded as a new physical theorem. No DFT/WSL/Slurm/Japan HPC work was authorized by this code.


## 16. Important correction: only the local even index is conditionally −2

The exact source-geometry toy script was rerun locally against the user's BiDB SQLite. The baseline \(r=0.125\) Å contour at \(\mu=1.5\) eV, SOC \(=0.48\) eV gives odd +1/even −2 at \(N_k=12,16\). But stress-testing \(\mu\) and SOC shows **finite-radius even winding is not universal**: at the same radius, \(\mu=1.5\) eV and SOC \(=0.72\) eV give even +4 at \(N_k=8,12\). At smaller \(r=0.02\)–0.05 Å the tested parameter sets recover the conditional local −2 index, with a very small minimum even-field norm. This is consistent with extra response zeros crossing the finite loop, not a violation of the local Taylor rule. The layer-odd +1 remained stable in nine tested parameter pairs. Full results: [BD30 audit and data](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/blob/main/papers/xiang-2025-spin-quantum-geometry/research/brainstorm-20261021/docs/BD30_LEE_SCOPE_AND_MODEL_ROBUSTNESS_AUDIT_20261022.md).

**Direct-prior-art gate:** published LEE already covers generic bilayer layer-spin response and 80 layer-group stacking routes. Its reference list cites a same-group submitted companion titled *Control of the Rashba-Edelstein effect via sliding ferroelectricity in bilayer systems*. Until that text is reviewed, continuous-sliding novelty is unresolved. The public repos should not claim a new universal finite-radius winding.

The local group-theory note [BD31](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/blob/main/papers/xiang-2025-spin-quantum-geometry/research/brainstorm-20261021/docs/BD31_LOCAL_EQUIVARIANT_RESPONSE_FIELD_RULE_20261022.md) derives \(F_-=i\gamma z+O(|z|^3)\) and \(F_+=\beta\bar z^2+O(|z|^4)\), conditional on nonzero leading coefficients. This is an application of standard equivariance, not a new mathematical theorem by itself. A read-only [parameter/radius winding script](pd2te2_registry_winding_scan.py) is now included.
