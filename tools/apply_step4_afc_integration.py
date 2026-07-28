#!/usr/bin/env python3
"""Apply the validated thermal AFC path to the solver call sites.

This temporary deterministic patcher exists because the validation branch is
being exercised directly on Sunbird. It refuses to continue unless every source
snippet exactly matches the audited branch state.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def replace_once(path: Path, old: str, new: str) -> None:
    text = path.read_text(encoding="utf-8")

    if new in text:
        print(f"already integrated: {path.relative_to(ROOT)}")
        return

    count = text.count(old)
    if count != 1:
        raise RuntimeError(
            f"expected exactly one audited match in {path}, found {count}"
        )

    path.write_text(text.replace(old, new, 1), encoding="utf-8")
    print(f"updated: {path.relative_to(ROOT)}")


def main() -> None:
    cmake = ROOT / "CMakeLists.txt"
    steps = ROOT / "src/solver/Steps.cpp"
    distributed = ROOT / "src/solver/SolverDistributedProductionLoop.cpp"

    replace_once(
        cmake,
        """    src/assembly/MomentumAssembly.cpp
    src/assembly/EnergyAssembly.cpp

    src/assembly/SpalartAllmarasAssembly.cpp
""",
        """    src/assembly/MomentumAssembly.cpp
    src/assembly/EnergyAssembly.cpp
    src/assembly/EnergyElementMatrix.cpp
    src/assembly/ThermalAfc.cpp

    src/assembly/SpalartAllmarasAssembly.cpp
""",
    )

    replace_once(
        cmake,
        """    src/parallel/PartitionMetadata.cpp
    src/parallel/HaloExchange.cpp
)
""",
        """    src/parallel/PartitionMetadata.cpp
    src/parallel/HaloExchange.cpp
    src/parallel/HaloExchangeExtrema.cpp
)
""",
    )

    replace_once(
        steps,
        """#include \"cbs/assembly/EnergyAssembly.hpp\"
#include \"cbs/assembly/MomentumAssembly.hpp\"
""",
        """#include \"cbs/assembly/EnergyAssembly.hpp\"
#include \"cbs/assembly/ThermalAfc.hpp\"
#include \"cbs/assembly/MomentumAssembly.hpp\"
""",
    )

    replace_once(
        steps,
        """#include <chrono>
#include <cmath>
""",
        """#include <chrono>
#include <cmath>
#include <cstdlib>
""",
    )

    replace_once(
        steps,
        """        EnergyAssembly::assembleStep4Rhs(s);

#pragma omp parallel for schedule(static)
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            s.temperature(ip) =
                s.temperature1(ip) + s.rhs1(ip) * s.elcoe2p(ip);
        }

        Boundary::applyTemperature(s);
""",
        """        EnergyAssembly::assembleStep4Rhs(s);

        const char* thermal_afc = std::getenv(\"CBS3D_THERMAL_AFC\");
        const bool afc_enabled =
            thermal_afc != nullptr
            && thermal_afc[0] != '\\0'
            && std::string(thermal_afc) != \"0\";

        if (afc_enabled)
        {
            ThermalAfc::applySerial(s);
            Boundary::applyTemperature(s);
            return;
        }

#pragma omp parallel for schedule(static)
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            s.temperature(ip) =
                s.temperature1(ip) + s.rhs1(ip) * s.elcoe2p(ip);
        }

        Boundary::applyTemperature(s);
""",
    )

    replace_once(
        distributed,
        """#include \"cbs/assembly/EnergyAssembly.hpp\"
#include \"cbs/assembly/MomentumAssembly.hpp\"
""",
        """#include \"cbs/assembly/EnergyAssembly.hpp\"
#include \"cbs/assembly/ThermalAfc.hpp\"
#include \"cbs/assembly/MomentumAssembly.hpp\"
""",
    )

    replace_once(
        distributed,
        """            // CBS Step 4: thermal assembly over fluid and solid elements.
            if (s_.cfg.temp_calc > 0)
            {
                {
                    ScopedHeatFluxSuppression suppress_heat_flux(s_);
                    EnergyAssembly::assembleStep4Rhs(s_);
                }

                HaloExchange::sumGhostContributionsToOwners(
                    s_.rhs1,
                    s_.partition_metadata,
                    MPI_COMM_WORLD);

                for (const Int ip : s_.owned_nodes)
                {
                    s_.rhs1(ip) +=
                        thermal_boundary.heat_flux_load(ip);

                    s_.temperature(ip) =
                        s_.temperature1(ip) +
                        s_.rhs1(ip) * s_.elcoe2p(ip);
                }

                apply_owned_temperature_constraints(
                    s_,
                    thermal_boundary);

                HaloExchange::broadcastOwnedToGhosts(
                    s_.temperature,
                    s_.partition_metadata,
                    MPI_COMM_WORLD);
            }
""",
        """            // CBS Step 4: thermal assembly over fluid and solid elements.
            if (s_.cfg.temp_calc > 0)
            {
                {
                    ScopedHeatFluxSuppression suppress_heat_flux(s_);
                    EnergyAssembly::assembleStep4Rhs(s_);
                }

                if (environment_flag_enabled(\"CBS3D_THERMAL_AFC\"))
                {
                    ThermalAfc::applyDistributed(
                        s_,
                        thermal_boundary.fixed,
                        thermal_boundary.value,
                        thermal_boundary.heat_flux_load,
                        MPI_COMM_WORLD);
                }
                else
                {
                    HaloExchange::sumGhostContributionsToOwners(
                        s_.rhs1,
                        s_.partition_metadata,
                        MPI_COMM_WORLD);

                    for (const Int ip : s_.owned_nodes)
                    {
                        s_.rhs1(ip) +=
                            thermal_boundary.heat_flux_load(ip);

                        s_.temperature(ip) =
                            s_.temperature1(ip) +
                            s_.rhs1(ip) * s_.elcoe2p(ip);
                    }

                    apply_owned_temperature_constraints(
                        s_,
                        thermal_boundary);

                    HaloExchange::broadcastOwnedToGhosts(
                        s_.temperature,
                        s_.partition_metadata,
                        MPI_COMM_WORLD);
                }
            }
""",
    )

    print("Step-4 AFC integration applied successfully.")


if __name__ == "__main__":
    main()
