# Draft PR scope

Included:

- continuous Spalding wall-law algebra;
- robust tangential wall-shear evaluation;
- unit and wide-range numerical tests;
- CI validation gates;
- weak-form and staged-integration documentation.

Excluded:

- any change to BC 530/532/901;
- any Step-1 traction assembly;
- SA transported-variable wall-function boundary condition;
- thermal wall treatment;
- QCR2000.

This separation is intentional so the algebraic foundation can be reviewed and validated independently.
