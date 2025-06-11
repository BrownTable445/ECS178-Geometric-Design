#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/boundary_loop.h>

#include <stdexcept>
#ifndef TRIANGLE_SURF_HH
#define TRIANGLE_SURF_HH


struct TriangleSurface {
    using VXd = Eigen::VectorXd;
    using SpMat = Eigen::SparseMatrix<double>;

    TriangleSurface() {
        load(std::string(MESH_PATH) + "/armadillo.obj");
    }

    SpMat laplacian() const;

    // Per-vertex differential geometry properties
    Eigen::MatrixXd normals() const;
    VXd        voronoiAreas() const;
    VXd  gaussianCurvatures() const;
    VXd      meanCurvatures() const;
    VXd             kappa_1() const;
    VXd             kappa_2() const;

    void load(const std::string &path) {
        if (path.size() == 0) throw std::runtime_error("Empty path");

        Eigen::MatrixX3d V_new;
        Eigen::MatrixX3i F_new;
        bool success = igl::read_triangle_mesh(path, V_new, F_new);

        if (!success) throw std::runtime_error("Load failed");

        std::cout << "Loaded vertices, faces: " << V_new.rows() << ", " << F_new.rows() << std::endl;
        setMesh(V_new, F_new);
    }

    void save(const std::string &path) const {
        if (path.size() == 0) throw std::runtime_error("Empty path");

        bool success = igl::write_triangle_mesh(path, V, F);
        if (!success) throw std::runtime_error("write failed");
    }

    void setMesh(const Eigen::MatrixX3d &V_new, const Eigen::MatrixX3i &F_new) {
        V = V_new;
        F = F_new;

        std::vector<std::vector<int>> boundaryLoop;
        igl::boundary_loop(F, boundaryLoop);
        const int nv = V.rows();
        m_onBoundary.assign(nv, false);
        for (const auto &l : boundaryLoop) {
            for (int i : l) m_onBoundary[i] = true;
        }
    }

    template<typename Derived>
    void zeroOutBoundaryValues(Eigen::MatrixBase<Derived> &val) const {
        for (size_t i = 0; i < val.rows(); ++i)
            if (m_onBoundary[i]) val.row(i).setZero();
    }

    Eigen::MatrixX3d V;
    Eigen::MatrixX3i F;
private:
    std::vector<bool> m_onBoundary;
};

#endif /* end of include guard: TRIANGLE_SURF_HH */
