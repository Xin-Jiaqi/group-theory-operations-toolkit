# Edelstein response tensor

`response_tensor_basis(point_group, "edelstein")` solves the homogeneous
spatial selection rule for the linear relation

$$
S_i=\chi_{ij}j_j,
$$

where current $\mathbf j$ is a polar vector and spin magnetization
$\mathbf S$ is an axial vector.  For every point operation $R$, it uses

$$
\chi=\det(R)R\chi R^T.
$$

Thus an inversion-containing point group has no ordinary homogeneous
Edelstein tensor.  The result is a 3 by 3 real tensor basis; its rows are
spin components $(S_x,S_y,S_z)$ and its columns are current components
$(j_x,j_y,j_z)$.  For a two-dimensional in-plane transport problem, read
the first two columns.

For example, `4mm` returns the one-dimensional Rashba basis

$$
\begin{pmatrix}
0&-1&0\\
1&0&0\\
0&0&0
\end{pmatrix},
$$

which gives $\mathbf S_\parallel\propto\hat z\times\mathbf j$.
The alias `rashba-edelstein` is accepted as well.  This API reports only
the ordinary point-group-allowed response of one homogeneous system; it
does not by itself classify layer-resolved or electric-field-induced LEE.
