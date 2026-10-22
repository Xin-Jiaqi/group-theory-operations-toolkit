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
