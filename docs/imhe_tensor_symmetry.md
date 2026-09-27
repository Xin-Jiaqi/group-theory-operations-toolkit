# Intrinsic magnetoelectric Hall tensor symmetry

This document gives an independent symmetry derivation for the planar IMHE response discussed by Das and Agarwal, arXiv:2604.20249v1.

The response is

$$
j_a=\chi_{ab;zz}E_b\mathcal E B,
\qquad a,b\in\{x,y\}.
$$

Define

$$
X=
\begin{pmatrix}
\chi_{xx;zz}&\chi_{xy;zz}\\
\chi_{yx;zz}&\chi_{yy;zz}
\end{pmatrix},
$$

so that

$$
\mathbf j_\parallel=X\mathbf E_\parallel\mathcal E B.
$$

The microscopic transport derivation gives two coefficient subspaces:

$$
X_{\rm in}
=
H
\begin{pmatrix}
0&-1\\
1&0
\end{pmatrix},
$$

for intrinsic IMHE, and

$$
X_{\rm ex}
=
\begin{pmatrix}
A&C\\
C&D
\end{pmatrix},
$$

for the extrinsic sector.

For a layer-normal-preserving operation

$$
R=
\begin{pmatrix}
R_\parallel&0\\
0&s_z
\end{pmatrix},
$$

the gate field is a polar z component and the magnetic field is an axial z component, so

$$
\mathcal E'B'=(\det R)\mathcal E B.
$$

Neumann's principle gives

$$
\boxed{
R_\parallel X=(\det R)XR_\parallel
}.
$$

For a magnetic operation $(R,\theta)$, with $\theta=1$ when time reversal is included, the condition becomes

$$
\boxed{
R_\parallel X
=
\eta_T^\theta(\det R)XR_\parallel
}.
$$

The two sectors have

$$
\eta_T^{\rm in}=+1,
\qquad
\eta_T^{\rm ex}=-1.
$$

This directly explains why pure time reversal can retain intrinsic IMHE but eliminates the extrinsic coefficient.

## Key analytic checks

For inversion $P$ and horizontal mirror $M_z$, both sectors vanish.

For a vertical mirror $M_x$,

$$
X_{\rm in}
=
H
\begin{pmatrix}
0&-1\\
1&0
\end{pmatrix},
$$

while

$$
X_{\rm ex}
=
\begin{pmatrix}
0&C\\
C&0
\end{pmatrix}.
$$

Thus $M_x$ forbids $\chi^{\rm ex}_{xx;zz}$ but not the entire extrinsic tensor.

For $C_{2x}$, intrinsic IMHE is forbidden and the extrinsic tensor is diagonal.

For $C_{3z}$,

$$
X_{\rm ex}=AI,
$$

so

$$
\chi^{\rm ex}_{yx;zz}=0,
\qquad
\chi^{\rm ex}_{xx;zz}
=
\chi^{\rm ex}_{yy;zz}.
$$

The intrinsic Hall component remains allowed. This independently establishes the clean transverse separation used for the R5G discussion without relying on the missing Supplemental Material.

For $C_{3z}T$, intrinsic IMHE remains allowed while the symmetric extrinsic tensor is forced to zero.

For $C_{4z}T$,

$$
X_{\rm ex}
=
\begin{pmatrix}
A&C\\
C&-A
\end{pmatrix},
$$

so the longitudinal extrinsic component may survive.

For $PT$, the independent transformation law gives

$$
X_{\rm in}=0,
$$

but allows the symmetric extrinsic tensor. This disagrees with the extrinsic PT entry in arXiv:2604.20249v1 Table I and is therefore retained as an explicit validation discrepancy rather than adjusted to match the paper.

The repeated $C_{2z}T$, $C_{3z}T$ and $C_{6z}T$ labels in Table I are mutually incompatible as printed. The solver retains its result as an independent check; the missing SM S2 is required before assigning the discrepancy to a specific published row.

## Scope

The fixed component set $\chi_{ab;zz}$ is closed only for operations that preserve the chosen z axis up to sign. Therefore:

- the implementation supports 106 of the 122 standard magnetic point-group embeddings;
- the 16 magnetic point groups with cubic parent symmetry are rejected because some operations rotate z into x or y;
- all 528 magnetic layer groups are supported, because the layer normal is intrinsic to the group definition.

A complete all-122 treatment would require defining a full 3D parent response tensor before selecting the zz field components.

## API

For magnetic point groups:

```python
from group_theory_operations import magnetic_imhe_tensor_basis

intrinsic = magnetic_imhe_tensor_basis("3", "intrinsic")
extrinsic = magnetic_imhe_tensor_basis("3", "extrinsic")
```

For magnetic layer groups:

```python
from group_theory_operations import magnetic_layer_imhe_tensor_basis

result = magnetic_layer_imhe_tensor_basis("6.5.25", "intrinsic")
```

The solver uses the repository's existing canonical operation matrices and explicit time-reversal labels. Published tables are treated as cross-checks rather than hard-coded ground truth.
