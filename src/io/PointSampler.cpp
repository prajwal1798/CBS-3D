//=============================================================================
// CBS3D++_SI
//
// FILE      : PointSampler.cpp
// PURPOSE   : PROBE TRACES AND REGULAR-GRID SAMPLES, SEE PointSampler.hpp.
//=============================================================================

#include "cbs/io/PointSampler.hpp"
#include "cbs/assembly/HeatPulse.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

namespace fs = std::filesystem;

namespace cbs
{
    namespace
    {
        // ================================
        // ONE SET OF POINTS AND THE PART OF IT THIS RANK EVALUATES
        // ================================
        struct PointSet
        {
            bool enabled = false;
            std::size_t npoints = 0;
            std::vector<std::array<Real, 3>> xyz;      // ALL POINTS (EVERY RANK)
            std::vector<std::size_t> mine;             // POINT INDICES OWNED HERE
            std::vector<Int> element;                  // CONTAINING LOCAL ELEMENT
            std::vector<std::array<Real, 4>> weight;   // BARYCENTRIC WEIGHTS
            std::vector<int> found;                    // 1 IF SOME RANK OWNS IT (RANK 0)
            Real dt = 0.0;                             // TIME INTERVAL (TRANSIENT)
            Int every = 0;                             // ITERATION INTERVAL (STEADY)
            Real next_time = 0.0;
        };

        PointSet probes;
        PointSet grid;
        std::vector<std::string> probe_names;
        int grid_n[3] = {0, 0, 0};
        Real grid_box[6] = {0, 0, 0, 0, 0, 0};

        const char* const field_names[PointSampler::nfields] =
            {"T", "u", "v", "w", "p", "wx", "wy", "wz"};

        std::string env(const char* name, const std::string& fallback = "")
        {
            const char* v = std::getenv(name);
            return (v == nullptr) ? fallback : std::string(v);
        }

        fs::path out_dir(const std::string& case_name)
        {
            return fs::path("output") / case_name;
        }

        // ================================
        // MPI HELPERS (NO-OPS IN A SERIAL BUILD)
        // ================================
        bool distributed(const CBSStateSI& s)
        {
#ifdef CBS3D_USE_MPI
            return s.mpi_enabled;
#else
            (void)s;
            return false;
#endif
        }

        void allreduce_min_int(const CBSStateSI& s, std::vector<int>& v)
        {
#ifdef CBS3D_USE_MPI
            if (distributed(s) && !v.empty())
            {
                MPI_Allreduce(MPI_IN_PLACE, v.data(), static_cast<int>(v.size()), MPI_INT, MPI_MIN,
                              MPI_COMM_WORLD);
            }
#else
            (void)s;
            (void)v;
#endif
        }

        template <class T>
        void reduce_sum_to_root(const CBSStateSI& s, std::vector<T>& v)
        {
#ifdef CBS3D_USE_MPI
            if (distributed(s) && !v.empty())
            {
                const MPI_Datatype type = std::is_same<T, float>::value ? MPI_FLOAT
                                        : std::is_same<T, double>::value ? MPI_DOUBLE : MPI_INT;
                if (s.mpi_rank == 0)
                {
                    MPI_Reduce(MPI_IN_PLACE, v.data(), static_cast<int>(v.size()), type, MPI_SUM, 0,
                               MPI_COMM_WORLD);
                }
                else
                {
                    MPI_Reduce(v.data(), nullptr, static_cast<int>(v.size()), type, MPI_SUM, 0,
                               MPI_COMM_WORLD);
                }
            }
#else
            (void)s;
            (void)v;
#endif
        }

        void allreduce_box(const CBSStateSI& s, Real box[6])
        {
#ifdef CBS3D_USE_MPI
            if (distributed(s))
            {
                double lo[3] = {box[0], box[2], box[4]};
                double hi[3] = {box[1], box[3], box[5]};
                MPI_Allreduce(MPI_IN_PLACE, lo, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, hi, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                box[0] = lo[0]; box[2] = lo[1]; box[4] = lo[2];
                box[1] = hi[0]; box[3] = hi[1]; box[5] = hi[2];
            }
#else
            (void)s;
            (void)box;
#endif
        }

        // ================================
        // TETRAHEDRON GEOMETRY
        // ================================
        void element_nodes(const CBSStateSI& s, Int ie, Real x[4][3])
        {
            for (int a = 0; a < 4; ++a)
            {
                const Int node = s.intma(a + 1, ie);

                for (int d = 0; d < 3; ++d)
                {
                    x[a][d] = s.coord(d + 1, node);
                }
            }
        }

        // BARYCENTRIC WEIGHTS OF p; FALSE IF p IS OUTSIDE (TOLERANCE tol)
        bool barycentric(const Real x[4][3], const Real p[3], Real tol, std::array<Real, 4>& w)
        {
            const Real a[3] = {x[1][0] - x[0][0], x[1][1] - x[0][1], x[1][2] - x[0][2]};
            const Real b[3] = {x[2][0] - x[0][0], x[2][1] - x[0][1], x[2][2] - x[0][2]};
            const Real c[3] = {x[3][0] - x[0][0], x[3][1] - x[0][1], x[3][2] - x[0][2]};
            const Real r[3] = {p[0] - x[0][0], p[1] - x[0][1], p[2] - x[0][2]};

            const Real det = a[0] * (b[1] * c[2] - b[2] * c[1])
                           - b[0] * (a[1] * c[2] - a[2] * c[1])
                           + c[0] * (a[1] * b[2] - a[2] * b[1]);

            if (std::abs(det) < 1e-300)
            {
                return false;
            }

            // CRAMER'S RULE FOR r = l1 a + l2 b + l3 c
            const Real l1 = (r[0] * (b[1] * c[2] - b[2] * c[1])
                           - b[0] * (r[1] * c[2] - r[2] * c[1])
                           + c[0] * (r[1] * b[2] - r[2] * b[1])) / det;
            const Real l2 = (a[0] * (r[1] * c[2] - r[2] * c[1])
                           - r[0] * (a[1] * c[2] - a[2] * c[1])
                           + c[0] * (a[1] * r[2] - a[2] * r[1])) / det;
            const Real l3 = (a[0] * (b[1] * r[2] - b[2] * r[1])
                           - b[0] * (a[1] * r[2] - a[2] * r[1])
                           + r[0] * (a[1] * b[2] - a[2] * b[1])) / det;
            const Real l0 = 1.0 - l1 - l2 - l3;

            if (l0 < -tol || l1 < -tol || l2 < -tol || l3 < -tol)
            {
                return false;
            }

            w = {l0, l1, l2, l3};
            return true;
        }

        Real dN(const CBSStateSI& s, Int ie, Int dim, Int a)
        {
            const Int index = (ie - 1) * s.cfg.ndim * s.cfg.nep + (dim - 1) * s.cfg.nep + a;
            return s.dNkdx(index);
        }

        // ================================
        // FIND THE OWNED ELEMENT OF EVERY POINT; KEEP THE POINTS THIS RANK
        // WINS (LOWEST RANK AMONG THOSE THAT FOUND IT)
        // ================================
        void locate(const CBSStateSI& s, PointSet& set, bool regular)
        {
            const std::size_t n = set.npoints;
            std::vector<Int> local_element(n, 0);
            std::vector<std::array<Real, 4>> local_weight(n);

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                Real x[4][3];
                element_nodes(s, ie, x);
                Real lo[3];
                Real hi[3];

                for (int d = 0; d < 3; ++d)
                {
                    lo[d] = std::min(std::min(x[0][d], x[1][d]), std::min(x[2][d], x[3][d]));
                    hi[d] = std::max(std::max(x[0][d], x[1][d]), std::max(x[2][d], x[3][d]));
                }

                const Real size = std::max(hi[0] - lo[0], std::max(hi[1] - lo[1], hi[2] - lo[2]));
                const Real tol = 1e-9;

                auto test = [&](std::size_t ip) {
                    if (local_element[ip] != 0)
                    {
                        return;
                    }
                    const Real p[3] = {set.xyz[ip][0], set.xyz[ip][1], set.xyz[ip][2]};
                    for (int d = 0; d < 3; ++d)
                    {
                        if (p[d] < lo[d] - tol * size || p[d] > hi[d] + tol * size)
                        {
                            return;
                        }
                    }
                    std::array<Real, 4> w;
                    if (barycentric(x, p, tol, w))
                    {
                        local_element[ip] = ie;
                        local_weight[ip] = w;
                    }
                };

                if (!regular)
                {
                    for (std::size_t ip = 0; ip < n; ++ip)
                    {
                        test(ip);
                    }
                    continue;
                }

                // REGULAR GRID: ONLY THE CELL CENTRES INSIDE THE ELEMENT BOX
                int range[3][2];
                bool empty = false;

                for (int d = 0; d < 3; ++d)
                {
                    const Real h = (grid_box[2 * d + 1] - grid_box[2 * d]) / grid_n[d];
                    const Real g0 = grid_box[2 * d];
                    range[d][0] = std::max(0, static_cast<int>(std::ceil((lo[d] - g0) / h - 0.5 - 1e-9)));
                    range[d][1] = std::min(grid_n[d] - 1,
                                           static_cast<int>(std::floor((hi[d] - g0) / h - 0.5 + 1e-9)));
                    empty = empty || range[d][0] > range[d][1];
                }

                if (empty)
                {
                    continue;
                }

                for (int k = range[2][0]; k <= range[2][1]; ++k)
                {
                    for (int j = range[1][0]; j <= range[1][1]; ++j)
                    {
                        for (int i = range[0][0]; i <= range[0][1]; ++i)
                        {
                            test(static_cast<std::size_t>(i) +
                                 static_cast<std::size_t>(grid_n[0]) *
                                     (static_cast<std::size_t>(j) +
                                      static_cast<std::size_t>(grid_n[1]) * static_cast<std::size_t>(k)));
                        }
                    }
                }
            }

            const int none = std::numeric_limits<int>::max();
            std::vector<int> owner(n, none);

            for (std::size_t ip = 0; ip < n; ++ip)
            {
                if (local_element[ip] != 0)
                {
                    owner[ip] = static_cast<int>(s.mpi_rank);
                }
            }

            allreduce_min_int(s, owner);

            set.mine.clear();
            set.element.clear();
            set.weight.clear();
            set.found.assign(n, 0);

            for (std::size_t ip = 0; ip < n; ++ip)
            {
                set.found[ip] = owner[ip] != none ? 1 : 0;

                if (owner[ip] == static_cast<int>(s.mpi_rank))
                {
                    set.mine.push_back(ip);
                    set.element.push_back(local_element[ip]);
                    set.weight.push_back(local_weight[ip]);
                }
            }
        }

        // ================================
        // THE 8 FIELDS AT ONE OWNED POINT
        // ================================
        void evaluate(const CBSStateSI& s, Int ie, const std::array<Real, 4>& w, Real out[8])
        {
            for (int f = 0; f < 8; ++f)
            {
                out[f] = 0.0;
            }

            for (int a = 0; a < 4; ++a)
            {
                const Int node = s.intma(a + 1, ie);
                out[0] += w[a] * (s.cfg.temp_calc > 0 ? s.temperature(node) : 0.0);
                out[1] += w[a] * s.unkno(1, node);
                out[2] += w[a] * s.unkno(2, node);
                out[3] += w[a] * s.unkno(3, node);
                out[4] += w[a] * s.pres(node);
            }

            // CURL OF THE P1 VELOCITY: omega_i = eps_ijk du_k/dx_j
            Real g[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};  // g[k][j] = du_k/dx_j

            for (int a = 1; a <= 4; ++a)
            {
                const Int node = s.intma(a, ie);

                for (int k = 0; k < 3; ++k)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        g[k][j] += dN(s, ie, j + 1, a) * s.unkno(k + 1, node);
                    }
                }
            }

            out[5] = g[2][1] - g[1][2];
            out[6] = g[0][2] - g[2][0];
            out[7] = g[1][0] - g[0][1];
        }

        bool due(const CBSStateSI& s, PointSet& set, Int iteration, bool first)
        {
            if (!set.enabled)
            {
                return false;
            }

            if (s.cfg.transient_on > 0)
            {
                if (first)
                {
                    return s.cfg.rtime <= 0.0;
                }

                if (s.cfg.rtime + 1e-9 * set.dt >= set.next_time)
                {
                    while (set.next_time <= s.cfg.rtime + 1e-9 * set.dt)
                    {
                        set.next_time += set.dt;
                    }
                    return true;
                }

                return false;
            }

            return first ? (iteration == 0) : (set.every > 0 && iteration % set.every == 0);
        }

        void write_probes(const CBSStateSI& s, const std::string& case_name, Int iteration)
        {
            std::vector<double> buf(probes.npoints * 8, 0.0);

            for (std::size_t m = 0; m < probes.mine.size(); ++m)
            {
                Real v[8];
                evaluate(s, probes.element[m], probes.weight[m], v);

                for (int f = 0; f < 8; ++f)
                {
                    buf[probes.mine[m] * 8 + static_cast<std::size_t>(f)] = v[f];
                }
            }

            reduce_sum_to_root(s, buf);

            if (s.mpi_rank != 0)
            {
                return;
            }

            const fs::path path = out_dir(case_name) / (case_name + "_probes.csv");
            std::ofstream out(path, std::ios::app);
            out << iteration << ',' << std::setprecision(10) << s.cfg.rtime << ',' << s.heat_source_scale;
            out << std::setprecision(9);

            for (std::size_t ip = 0; ip < probes.npoints; ++ip)
            {
                for (int f = 0; f < 8; ++f)
                {
                    out << ',';
                    if (probes.found[ip])
                    {
                        out << buf[ip * 8 + static_cast<std::size_t>(f)];
                    }
                    else
                    {
                        out << "nan";
                    }
                }
            }

            out << '\n';
        }

        void write_grid(const CBSStateSI& s, const std::string& case_name, Int iteration)
        {
            const std::size_t n = grid.npoints;
            std::vector<float> buf(n * 8, 0.0f);

            for (std::size_t m = 0; m < grid.mine.size(); ++m)
            {
                Real v[8];
                evaluate(s, grid.element[m], grid.weight[m], v);

                for (int f = 0; f < 8; ++f)
                {
                    buf[static_cast<std::size_t>(f) * n + grid.mine[m]] = static_cast<float>(v[f]);
                }
            }

            reduce_sum_to_root(s, buf);

            if (s.mpi_rank != 0)
            {
                return;
            }

            std::ostringstream name;
            name << "sample_" << std::setw(8) << std::setfill('0') << iteration << ".f32";
            const fs::path dir = out_dir(case_name) / "fno_grid";
            std::ofstream bin(dir / name.str(), std::ios::binary | std::ios::trunc);
            bin.write(reinterpret_cast<const char*>(buf.data()),
                      static_cast<std::streamsize>(buf.size() * sizeof(float)));

            std::ofstream index(dir / "index.csv", std::ios::app);
            index << name.str() << ',' << iteration << ',' << std::setprecision(10) << s.cfg.rtime << ','
                  << s.heat_source_scale << '\n';
        }

        void read_probe_file(const std::string& path, std::vector<std::string>& names,
                             std::vector<std::array<Real, 3>>& xyz)
        {
            std::ifstream in(path);

            if (!in)
            {
                throw std::runtime_error("PointSampler - cannot open CBS3D_PROBE_FILE " + path);
            }

            std::string line;

            while (std::getline(in, line))
            {
                const std::size_t first = line.find_first_not_of(" \t\r");

                if (first == std::string::npos || line[first] == '#')
                {
                    continue;
                }

                std::istringstream row(line);
                std::string name;
                Real x = 0.0;
                Real y = 0.0;
                Real z = 0.0;

                if (!(row >> name >> x >> y >> z))
                {
                    throw std::runtime_error("PointSampler - bad probe line: " + line);
                }

                names.push_back(name);
                xyz.push_back({x, y, z});
            }
        }
    }

    // =========================================================================
    void PointSampler::initialise(CBSStateSI& s, const std::string& case_name)
    {
        probes = PointSet();
        grid = PointSet();
        probe_names.clear();

        const std::string probe_file = env("CBS3D_PROBE_FILE");
        const std::string grid_dims = env("CBS3D_GRID_DIMS");

        if (probe_file.empty() && grid_dims.empty())
        {
            return;
        }

        if (s.mpi_rank == 0)
        {
            fs::create_directories(out_dir(case_name));
        }

        // ================================
        // PROBES (EVERY RANK READS THE SMALL TEXT FILE)
        // ================================
        if (!probe_file.empty())
        {
            read_probe_file(probe_file, probe_names, probes.xyz);
            probes.npoints = probes.xyz.size();
            probes.enabled = probes.npoints > 0;
            probes.dt = std::atof(env("CBS3D_PROBE_DT", "0.01").c_str());
            probes.every = std::atoi(env("CBS3D_PROBE_EVERY", "100").c_str());
            locate(s, probes, false);

            if (s.mpi_rank == 0)
            {
                std::cout << "PROBES      : " << probes.npoints << " from " << probe_file << ", every "
                          << (s.cfg.transient_on > 0 ? std::to_string(probes.dt) + " s"
                                                     : std::to_string(probes.every) + " iterations")
                          << "\n";

                for (std::size_t ip = 0; ip < probes.npoints; ++ip)
                {
                    std::cout << "  " << std::setw(16) << std::left << probe_names[ip] << std::right
                              << " (" << probes.xyz[ip][0] << ", " << probes.xyz[ip][1] << ", "
                              << probes.xyz[ip][2] << ")" << (probes.found[ip] ? "" : "  NOT IN THE MESH")
                              << "\n";
                }

                const fs::path path = out_dir(case_name) / (case_name + "_probes.csv");

                if (!fs::exists(path))
                {
                    std::ofstream out(path);
                    out << "iteration,time,heat_scale";

                    for (const std::string& name : probe_names)
                    {
                        for (int f = 0; f < 8; ++f)
                        {
                            out << ',' << name << '_' << field_names[f];
                        }
                    }

                    out << '\n';
                }
            }
        }

        // ================================
        // REGULAR GRID
        // ================================
        if (!grid_dims.empty())
        {
            std::istringstream dims(grid_dims);

            if (!(dims >> grid_n[0] >> grid_n[1] >> grid_n[2]) || grid_n[0] < 1 || grid_n[1] < 1 ||
                grid_n[2] < 1)
            {
                throw std::runtime_error("PointSampler - CBS3D_GRID_DIMS must be 'nx ny nz' (all >= 1)");
            }

            const std::string bounds = env("CBS3D_GRID_BOUNDS");

            if (!bounds.empty())
            {
                std::istringstream b(bounds);

                for (int d = 0; d < 6; ++d)
                {
                    if (!(b >> grid_box[d]))
                    {
                        throw std::runtime_error("PointSampler - CBS3D_GRID_BOUNDS must be 'x0 x1 y0 y1 z0 z1'");
                    }
                }
            }
            else
            {
                const Real big = std::numeric_limits<Real>::max();
                Real box[6] = {big, -big, big, -big, big, -big};

                for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                {
                    for (int a = 1; a <= 4; ++a)
                    {
                        const Int node = s.intma(a, ie);

                        for (int d = 0; d < 3; ++d)
                        {
                            box[2 * d] = std::min(box[2 * d], s.coord(d + 1, node));
                            box[2 * d + 1] = std::max(box[2 * d + 1], s.coord(d + 1, node));
                        }
                    }
                }

                allreduce_box(s, box);
                std::copy(box, box + 6, grid_box);
            }

            grid.npoints = static_cast<std::size_t>(grid_n[0]) * static_cast<std::size_t>(grid_n[1]) *
                           static_cast<std::size_t>(grid_n[2]);
            grid.xyz.resize(grid.npoints);

            for (int k = 0; k < grid_n[2]; ++k)
            {
                for (int j = 0; j < grid_n[1]; ++j)
                {
                    for (int i = 0; i < grid_n[0]; ++i)
                    {
                        const int ijk[3] = {i, j, k};
                        std::array<Real, 3> p;

                        for (int d = 0; d < 3; ++d)
                        {
                            const Real h = (grid_box[2 * d + 1] - grid_box[2 * d]) / grid_n[d];
                            p[d] = grid_box[2 * d] + (ijk[d] + 0.5) * h;
                        }

                        grid.xyz[static_cast<std::size_t>(i) +
                                 static_cast<std::size_t>(grid_n[0]) *
                                     (static_cast<std::size_t>(j) +
                                      static_cast<std::size_t>(grid_n[1]) * static_cast<std::size_t>(k))] = p;
                    }
                }
            }

            grid.enabled = true;
            grid.dt = std::atof(env("CBS3D_GRID_DT", "0.1").c_str());
            grid.every = std::atoi(env("CBS3D_GRID_EVERY", "1000").c_str());
            locate(s, grid, true);

            // MATERIAL MASK ON RANK 0: OWNER CONTRIBUTES (MATERIAL + 1)
            std::vector<int> mask(grid.npoints, 0);

            for (std::size_t m = 0; m < grid.mine.size(); ++m)
            {
                mask[grid.mine[m]] = static_cast<int>(s.mat_elem(grid.element[m])) + 1;
            }

            reduce_sum_to_root(s, mask);

            if (s.mpi_rank == 0)
            {
                const fs::path dir = out_dir(case_name) / "fno_grid";
                fs::create_directories(dir);
                std::vector<std::int8_t> m8(grid.npoints);
                std::size_t inside = 0;

                for (std::size_t ip = 0; ip < grid.npoints; ++ip)
                {
                    m8[ip] = static_cast<std::int8_t>(mask[ip] - 1);
                    inside += mask[ip] > 0 ? 1 : 0;
                }

                std::ofstream mfile(dir / "mask.i8", std::ios::binary | std::ios::trunc);
                mfile.write(reinterpret_cast<const char*>(m8.data()), static_cast<std::streamsize>(m8.size()));

                std::ofstream json(dir / "grid.json", std::ios::trunc);
                json << std::setprecision(12) << "{\n"
                     << "  \"dims\": [" << grid_n[0] << ", " << grid_n[1] << ", " << grid_n[2] << "],\n"
                     << "  \"bounds\": [" << grid_box[0] << ", " << grid_box[1] << ", " << grid_box[2] << ", "
                     << grid_box[3] << ", " << grid_box[4] << ", " << grid_box[5] << "],\n"
                     << "  \"points\": \"cell centres: x_i = x0 + (i + 0.5)(x1 - x0)/nx\",\n"
                     << "  \"fields\": [\"T\", \"u\", \"v\", \"w\", \"p\", \"omega_x\", \"omega_y\", \"omega_z\"],\n"
                     << "  \"units\": [\"K\", \"m/s\", \"m/s\", \"m/s\", \"solver pressure\", \"1/s\", \"1/s\", \"1/s\"],\n"
                     << "  \"dtype\": \"float32 little-endian\",\n"
                     << "  \"layout\": \"[field][k][j][i], i fastest\",\n"
                     << "  \"mask\": \"mask.i8: int8 per point, -1 outside the mesh, else material id (0 water)\",\n"
                     << "  \"index\": \"index.csv: file,iteration,time,heat_scale\",\n"
                     << "  \"points_inside\": " << inside << "\n}\n";

                if (!fs::exists(dir / "index.csv"))
                {
                    std::ofstream index(dir / "index.csv");
                    index << "file,iteration,time,heat_scale\n";
                }

                std::cout << "GRID SAMPLE : " << grid_n[0] << " x " << grid_n[1] << " x " << grid_n[2]
                          << " (" << inside << " of " << grid.npoints << " points inside the mesh), every "
                          << (s.cfg.transient_on > 0 ? std::to_string(grid.dt) + " s"
                                                     : std::to_string(grid.every) + " iterations")
                          << ", " << std::setprecision(3)
                          << grid.npoints * 8 * 4 / 1.0e6 << " MB per sample -> " << dir.string() << "\n";
            }
        }

        // ================================
        // FIRST SAMPLE TIMES; THE INITIAL STATE OF A FRESH RUN IS SAMPLED NOW
        // ================================
        for (PointSet* set : {&probes, &grid})
        {
            if (set->enabled && s.cfg.transient_on > 0)
            {
                if (!(set->dt > 0.0))
                {
                    throw std::runtime_error("PointSampler - CBS3D_PROBE_DT / CBS3D_GRID_DT must be > 0");
                }

                set->next_time = (std::floor(s.cfg.rtime / set->dt + 1e-9) + 1.0) * set->dt;
            }
        }

        const Int start = s.cfg.iiter_total;

        if (due(s, probes, start, true))
        {
            write_probes(s, case_name, start);
        }

        if (due(s, grid, start, true))
        {
            write_grid(s, case_name, start);
        }
    }

    // =========================================================================
    void PointSampler::sample(CBSStateSI& s, const std::string& case_name, Int iteration)
    {
        if (due(s, probes, iteration, false))
        {
            write_probes(s, case_name, iteration);
        }

        if (due(s, grid, iteration, false))
        {
            write_grid(s, case_name, iteration);
        }
    }
}
