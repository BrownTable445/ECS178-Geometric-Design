#include "triangle_surf.hh"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

auto TriangleSurface::laplacian() const -> SpMat {
    SpMat L;
    igl::cotmatrix(V, F, L);
    return L;
}

Eigen::VectorXd TriangleSurface::voronoiAreas() const {
    SpMat M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    return M.diagonal();
}

Eigen::MatrixXd TriangleSurface::normals() const {
    // TODO 3.1: construct the per-vertex area-weighted normal vector field.
    Eigen::MatrixX3d normal_per_triangle(F.rows(), 3);
    Eigen::MatrixX3d normal_per_vertex = Eigen::MatrixX3d::Zero(V.rows(), 3);
    for (size_t i = 0; i < F.rows(); ++i) {
        Eigen::Vector3d v1 = V.row(F.row(i)[0]) - V.row(F.row(i)[1]);
        Eigen::Vector3d v2 = V.row(F.row(i)[0]) - V.row(F.row(i)[2]);
        normal_per_triangle.row(i) = v1.cross(v2) / (2 * v1.cross(v2).norm());
    }
    for (size_t i = 0; i < F.rows(); ++i) {
        for (size_t j = 0; j < 3; ++j) {
            size_t current_vertex = F.row(i)[j];
            normal_per_vertex.row(current_vertex) += normal_per_triangle.row(i);
        }
    }
    for (size_t i = 0; i < normal_per_vertex.rows(); ++i) {
        normal_per_vertex.row(i).normalize();
    }
    return normal_per_vertex;
}

Eigen::VectorXd TriangleSurface::gaussianCurvatures() const {
    // TODO 3.2: calculate the vector of per-vertex discrete Gaussian
    // curvatures using the angle defect formula.
    // Be sure to zero out the boundary values using `zeroOutBoundaryValues`
    VXd K = VXd::Zero(V.rows()), voronoi_area = voronoiAreas();
    for (size_t i = 0; i < F.rows(); ++i) {
        Eigen::Vector3d v1 = V.row(F.row(i)[0]) - V.row(F.row(i)[1]);
        Eigen::Vector3d v2 = V.row(F.row(i)[0]) - V.row(F.row(i)[2]);
        K[F.row(i)[0]] += atan2(v1.cross(v2).norm(), v1.dot(v2));

        v1 = V.row(F.row(i)[1]) - V.row(F.row(i)[0]);
        v2 = V.row(F.row(i)[1]) - V.row(F.row(i)[2]);
        K[F.row(i)[1]] += atan2(v1.cross(v2).norm(), v1.dot(v2));

        v1 = V.row(F.row(i)[2]) - V.row(F.row(i)[0]);
        v2 = V.row(F.row(i)[2]) - V.row(F.row(i)[1]);
        K[F.row(i)[2]] += atan2(v1.cross(v2).norm(), v1.dot(v2));
    }
    for (size_t i = 0; i < K.rows(); ++i) {
        K[i] = (2 * M_PI - K[i]) / abs(voronoi_area[i]);
    }
    zeroOutBoundaryValues(K);
    return K;
    //return VXd::Zero(V.rows());
}

Eigen::VectorXd TriangleSurface::meanCurvatures() const {
    // TODO 3.3: calculate the vector of per-vertex discrete mean
    // curvatures using the discrete Laplace-Beltrami operator.
    // Be sure to zero out the boundary values using `zeroOutBoundaryValues`
    VXd h(V.rows());
    VXd A = voronoiAreas();
    SpMat L = laplacian();
    Eigen::MatrixX3d n = normals(), left_hand_side = L * V / 2;
    for (size_t i = 0; i < left_hand_side.rows(); ++i) {
        left_hand_side.row(i) /= abs(A[i]);
        h[i] = left_hand_side.row(i).dot(n.row(i)); //linear combintion
    }
    zeroOutBoundaryValues(h);
    return h;
    //return VXd::Zero(V.rows());
}

Eigen::VectorXd TriangleSurface::kappa_1() const {
    // TODO 3.4: calculate the vector of per-vertex first principal curvatures.
    Eigen::VectorXd K = gaussianCurvatures(), H = meanCurvatures(), k1(K.rows());
    for (size_t i = 0; i < k1.rows(); ++i) {
        double discriminant = H[i] * H[i] - K[i];
        if (discriminant < 0) {
            discriminant = 0;
        }
        k1[i] = H[i] + sqrt(discriminant);
    }
    return k1;
    //return VXd::Zero(V.rows());
}

Eigen::VectorXd TriangleSurface::kappa_2() const {
    // TODO 3.4: calculate the vector of per-vertex second principal curvatures.
    Eigen::VectorXd K = gaussianCurvatures(), H = meanCurvatures(), k2(K.rows());
    for (size_t i = 0; i < k2.rows(); ++i) {
        double discriminant = H[i] * H[i] - K[i];
        if (discriminant < 0) {
            discriminant = 0;
        }
        k2[i] = H[i] - sqrt(discriminant);
    }
    return k2;
    //return VXd::Zero(V.rows());
}
