# SA Wall-Treatment CI

The foundation branch adds `.github/workflows/wall-treatment-foundation.yml` with four independent gates:

1. Linux algebraic wall-law verification with `-Wall -Wextra -Werror -pedantic`, AddressSanitizer and UndefinedBehaviorSanitizer.
2. Windows MSVC algebraic verification with `/W4 /WX`.
3. Existing serial/OpenMP solver build plus the pre-existing CTest suite.
4. Existing MPI/OpenMP/PETSc distributed solver build.

The algebraic jobs also run a wall-Reynolds-number sweep from `Re_y=1e-12` through `1e10` and require reconstruction of `Re_y=u+ y+` to relative error no larger than `2e-12`.

This CI is a foundation gate only. It does not constitute physical validation of a wall-modelled RANS calculation.
