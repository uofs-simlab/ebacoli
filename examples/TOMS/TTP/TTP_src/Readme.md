[./TenTusscher.f](./TenTusscher.f) is generated from the CellML repository.

Jacobian of this model is evaluated using the functions `f_d_COLX.f` which have been created through automatic differentiation. These represent the `X`th column of the Jacobian.

[./TTP\_ebacoli.f](./TTP_ebacoli.f) has two purposes:

1.  It combines the cell model with the tissue equation, and
2.  It wraps the CellML function evaluation and AD-generated Jacobian routines in a way callable from eBACOLI.
