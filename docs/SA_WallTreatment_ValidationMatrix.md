# SA Wall-Treatment Validation Matrix

This matrix separates what is verified in the algebraic foundation from what remains unverified before production activation.

| Item | Current status | Acceptance criterion |
|---|---|---|
| Tangential projection | Implemented/tested | `u_t.n = 0` to round-off |
| Face-normal orientation | Implemented/tested | `n -> -n` leaves predicted shear unchanged |
| Spalding wall-unit identity | Implemented/tested | `u+ y+ = U_t y/nu` within solver tolerance |
| Viscous asymptote | Implemented/tested | `tau_w -> mu U_t/y` for `y+ << 1` |
| Dissipation | Implemented/tested | `t_w.u_t <= 0` |
| Zero tangential speed | Implemented/tested | exactly zero modeled shear |
| Invalid input handling | Implemented/tested | fail fast on invalid geometry/properties |
| Wide `Re_y` robustness | Implemented/tested | finite solution and identity error <= `2e-12` for `1e-12..1e10` sweep |
| Linux sanitizers | CI gate | ASan/UBSan clean |
| Windows MSVC warnings | CI gate | `/W4 /WX` clean |
| Existing serial solver | CI regression gate | build + current CTest pass |
| Existing MPI/PETSc solver | CI regression gate | distributed build passes |
| TRI3 wall-load sign | Not implemented | synthetic FE traction test |
| Wall-face inventory | Not implemented | exact owned-face count and valid adjacent fluid tet |
| Wall sample definition | Not implemented | verified distance/sample on skewed TET4 and corners |
| Normal-only wall constraint | Not implemented | no penetration without strong tangential no-slip |
| Step-1 wall traction | Not implemented | replace, not duplicate, natural wall flux |
| MPI wall residual | Not implemented | P1/P2/P40 rank-independent wall load |
| SA wall scalar treatment | Not implemented | reference flat-plate agreement |
| Thermal wall model | Not implemented | heated-channel `Nu`, `T+`, heat balance |
| QCR2000 | Not implemented | square-duct secondary-flow benchmark |
| First-wall helium CHT | Not started | mesh/y+ independence and global energy balance |

No item marked "Not implemented" is implied to be correct by the present foundation branch.
