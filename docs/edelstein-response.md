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
the ordinary point-group-allowed response of one homogeneous system.

## Bilayer layer Edelstein effect

`bilayer_edelstein_point_groups(monolayer_point_group, bilayer_point_group)`
adds the layer degree of freedom for standard point-group embeddings.  It
partitions the bilayer operations into layer-preserving $R_B^+$ and
layer-exchanging $R_B^-$ subsets by their action on the layer normal.  For a
layer-exchange operation it uses

$$
\chi^{L'}=\det(R_B^-)R_B^-\chi^L(R_B^-)^T.
$$

The returned component-wise classification is `type_i`, `type_ii`, `both`,
or `forbidden`.  Type-I components already exist locally and obey
$\chi^{L'}=-\chi^L$ at zero field.  Type-II components are forbidden in both
the isolated monolayer and zero-field global bilayer, but become allowed when
a vertical field removes $R_B^-$ while retaining $R_B^+$.

For arbitrary, nonstandard coordinate embeddings, use
`bilayer_edelstein_analysis(monolayer_operations, intralayer_operations,
interlayer_operations)`.  The result is symmetry allowance only: it does not
calculate layer-resolved transport weights, interlayer tunnelling, or a
response magnitude.
