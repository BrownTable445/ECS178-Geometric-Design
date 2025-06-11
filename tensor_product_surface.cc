#include "tensor_product_surface.hh"
#include "spline.hh"
#include <array>

template<size_t Dimension>
Eigen::Matrix<float, Dimension + 1, 1> toHomogeneousCoordinates(const Eigen::Matrix<float, Dimension, 1> &x, float w) {
    Eigen::Matrix<float, Dimension + 1, 1> result;
    // TODO 1.1: convert rational B-spline control point data (x, w) into a single
    // vector [w x; w] in homogeneous coordinates.
    result.head(Dimension) = x;
    result[Dimension] = 1;
    result *= w;
    return result;
}

template<int Dimension>
Eigen::Matrix<float, Dimension, 1> fromHomogeneousCoordinates(const Eigen::Matrix<float, Dimension + 1, 1> &x) {
    // TODO 1.1: convert point `x` expressed in homogeneous coordinates into
    // its Euclidean coordinates representation.
    return x.head(Dimension) / x[Dimension];
}

template<size_t Dimension>
auto TensorProductSurface_T<Dimension>::evalPtHomogCoord(double u, double v) const -> HVec {
    size_t I_u = u_spline.findSegmentContaining(u);
    size_t I_v = v_spline.findSegmentContaining(v);
    size_t L_u = u_spline.numControlPts() - 1;
    size_t L_v = v_spline.numControlPts() - 1;
    size_t n_u = u_spline.degree();
    size_t n_v = v_spline.degree();

    // TODO 1.1: evaluate the NURBS surface point's homogeneous coordinates xtilde(u, v)

    Eigen::MatrixXf evaluated_row_points(n_u + 1, Dimension + 1);
    for (size_t r = 0; r < n_u + 1; ++r) {
        Eigen::MatrixXf row_spline_control_points(n_v + 1, Dimension + 1);
        auto temp = controlPts.row(I_u - (n_u - 1) + r);
        for (size_t c = 0; c < n_v + 1; ++c) {
            row_spline_control_points.row(c) = toHomogeneousCoordinates(temp[I_v - (n_v - 1) + c], weights(r, I_v - (n_v - 1) + c));
            v_spline.setControlPt(I_v - (n_v - 1) + c, row_spline_control_points.row(c));
        }
        evaluated_row_points.row(r) = v_spline.evalPt(SplineEvalMethod::DE_BOOR, v);
        u_spline.setControlPt(I_u - (n_u - 1) + r, evaluated_row_points.row(r));
    }
    HVec final_point = u_spline.evalPt(SplineEvalMethod::DE_BOOR, u);
    return final_point;
}

template<size_t Dimension>
auto TensorProductSurface_T<Dimension>::evalPt(double u, double v) const -> Vec {

    // TODO 1.1: evaluate the NURBS surface point's Euclidean coordinates x(u, v)
    return fromHomogeneousCoordinates<Dimension>(evalPtHomogCoord(u, v));
}

template<size_t Dimension>
void TensorProductSurface_T<Dimension>::eval(size_t resolution, MX2f &U, MXDf &V, Eigen::MatrixX3i &F) const {
    //resolution = 3 ;
    const size_t numCellsU = resolution - 1;
    const size_t numCellsV = resolution - 1;
    const size_t   numVtxU = resolution;
    const size_t   numVtxV = resolution;

    const auto &cellIdx = [=](size_t a, size_t b) { return a * numCellsV + b; };
    const auto & vtxIdx = [=](size_t a, size_t b) { return a *   numVtxV + b; };

    float u_min = u_spline.domainStart();
    float u_max = u_spline.domainEnd();
    float v_min = v_spline.domainStart();
    float v_max = v_spline.domainEnd();

    // TODO 1.2: sample the NURBS surface on a triangulated regular grid in the
    // UV domain, obtaining triangle mesh (V, F). Also record the UV
    // coordinates (parameter domain) of each vertex in `U`.
    U.resize((resolution + 1) * (resolution + 1), 2);
    V.resize((resolution + 1) * (resolution + 1), 3);
    F.resize(resolution * resolution * 2, 3);

    for (size_t a = 0; a <= resolution; ++a) {
        for (size_t b = 0; b <= resolution; ++b) {
            float alpha = (float)a / resolution;
            float beta = (float)b / resolution;
            float u = (1 - alpha) * u_min + alpha * u_max;
            float v = (1 - beta) * v_min + beta * v_max;
            U.row(a * (resolution + 1) + b) = Eigen::Vector2f(u, v);
            //leave F indexing for next for loop
            V.row(a * (resolution + 1) + b) = evalPt(u, v);
        }
    }
    size_t i = 0, counter = 0;
    while (i <= (resolution + 1) * (resolution + 1) - 1 - (resolution + 2)) {
        F.row(counter * 2) = Eigen::Vector3i(i, i + resolution + 2, i + 1);
        F.row(counter * 2 + 1) = Eigen::Vector3i(i, i + resolution + 1, i + resolution + 2);
        ++i; //move up
        if (i >= resolution && (i - resolution) % (resolution + 1) == 0) {
            ++i;
        }
        ++counter;
    }
}

template<size_t Dimension>
void TensorProductSurface_T<Dimension>::insertVKnot(double vbar) {
    // TODO 1.3: insert the knot `vbar` into each of the "v curves", creating a
    // new column of data points.
    ControlPtGrid newControlPts = controlPts;
    WeightGrid newWeights = weights;

    size_t rows = newControlPts.rows();
    size_t cols = newControlPts.cols();

    newControlPts.conservativeResize(rows, cols + 1);
    newWeights.conservativeResize(rows, cols + 1);
    for (size_t i = 0; i < rows; ++i) { 
        auto current_row = getControlPtRow(i);
        Eigen::MatrixXf current_row_control_points(cols, Dimension + 1);
        for (size_t j = 0; j < cols; ++j) {
            current_row_control_points.row(j).head(Dimension) = current_row[j];
            current_row_control_points.row(j)[Dimension] = weights(i, j);
        }
        HSpline v_spline_copy = v_spline;
        v_spline_copy.setControlPts(current_row_control_points);
        v_spline_copy.insertKnot(vbar);
        for (size_t j = 0; j < cols + 1; ++j) {
            newControlPts(i, j) = v_spline_copy.controlPt(j).head(Dimension);
            newWeights(i, j) = v_spline_copy.controlPts(j, Dimension);
        }
        if (i == rows - 1) {
            v_spline = v_spline_copy;
        }
    }

    controlPts = newControlPts;
    weights = newWeights;
}

template<size_t Dimension>
void TensorProductSurface_T<Dimension>::insertUKnot(double ubar) {
    // TODO 1.3: insert the knot `ubar` into each of the "u curves", creating a
    // new row of data points.
    ControlPtGrid newControlPts = controlPts;
    WeightGrid newWeights = weights;

    size_t rows = newControlPts.rows();
    size_t cols = newControlPts.cols();
    newControlPts.conservativeResize(rows + 1, cols);
    newWeights.conservativeResize(rows + 1, cols);
    for (size_t i = 0; i < cols; ++i) {
        auto current_col = getControlPtCol(i);
        Eigen::MatrixXf current_col_control_points(rows, Dimension + 1);
        for (size_t j = 0; j < rows; ++j) {
            current_col_control_points.row(j).head(Dimension) = current_col[j];
            current_col_control_points.row(j)[Dimension] = weights(j, i);
        }
        HSpline u_spline_copy = u_spline;
        u_spline_copy.setControlPts(current_col_control_points);
        u_spline_copy.insertKnot(ubar);
        for (size_t j = 0; j < rows + 1; ++j) {
            newControlPts(j, i) = u_spline_copy.controlPt(j).head(Dimension);
            newWeights(j, i) = u_spline_copy.controlPts(j, Dimension);
        }
        if (i == cols - 1) {
            u_spline = u_spline_copy;
        }
    }

    controlPts = newControlPts;
    weights = newWeights;
}

template<size_t Dimension>
std::shared_ptr<TensorProductSurface_T<Dimension>> ruled_surface(const NURBS<Dimension> &spline_1, const NURBS<Dimension> &spline_2) {
    auto result_ptr = std::make_shared<TensorProductSurface_T<Dimension>>();
    auto &result = *result_ptr;

    const size_t n = spline_1.degree();
    if (n != spline_2.degree()) throw std::runtime_error("Splines must be of the same degree");

    result.controlPts.resize(spline_1.numControlPts(), 2);
    result.weights   .resize(spline_1.numControlPts(), 2);
    result.getControlPtColFlattened(0) = spline_1.getControlPts();
    result.getControlPtColFlattened(1) = spline_2.getControlPts();
    result.weights.col(0) = spline_1.weights();
    result.weights.col(1) = spline_2.weights();
    result.u_spline = spline_1.homogeneousSpline();

    // Convert v_spline into a simple spline/Bézier curve with just two control
    // points and knots.
    result.v_spline.degree() = 1;
    result.v_spline.controlPts.conservativeResize(2, Dimension + 1);
    result.v_spline.knots.setLinSpaced(2, 0, 1);

    return result_ptr;
}

template<size_t Dimension>
std::shared_ptr<TensorProductSurface_T<Dimension>> surface_of_revolution(const NURBS<Dimension> &spline, size_t numSegments) {
    NURBS<Dimension> circ = NURBSCircle<Dimension>(1.0, numSegments);
    auto result_ptr = std::make_shared<TensorProductSurface_T<Dimension>>();
    auto &result = *result_ptr;

    result.controlPts.resize(circ.numControlPts(), spline.numControlPts());
    result.weights   .resize(circ.numControlPts(), spline.numControlPts());

    for (size_t i = 0; i < circ.numControlPts(); ++i) {
        for (size_t j = 0; j < spline.numControlPts(); ++j) {
            float r = spline.controlPt(j)[0];
            // Use the circ control point to rotate (and scale) around the y axis.
            result.controlPts(i, j) << r * circ  .controlPt(i)[1],
                                           spline.controlPt(j)[1],
                                       r * circ  .controlPt(i)[0];
            result.weights(i, j) = circ.weight(i) * spline.weight(j);
        }
    }

    result.u_spline = circ.homogeneousSpline();
    result.v_spline = spline.homogeneousSpline();

    return result_ptr;
}

template<size_t Dimension>
auto TensorProductSurface_T<Dimension>::partial_xtilde(double u, double v, size_t diff_u, size_t diff_v) const -> HVec {
    // TODO 2.2: compute the partial derivative of the homogeneous coordinate
    // vector "xtilde" `diff_u` times with respect to `u` and `diff_v` times
    // with respect to `v` (applying `diff_u + diff_v` derivatives in total).
    size_t rows = controlPts.rows();
    size_t cols = controlPts.cols();
    Eigen::MatrixXf evaluated_row_points_derivative(rows, Dimension + 1);

    for (size_t r = 0; r < rows; ++r) {
        Eigen::MatrixXf row_spline_control_points(cols, Dimension + 1);
        auto temp = controlPts.row(r);
        for (size_t c = 0; c < cols; ++c) {
            row_spline_control_points.row(c).head(Dimension) = temp[c];
            row_spline_control_points.row(c)[Dimension] = weights(r, c);
        }
        v_spline.setControlPts(row_spline_control_points);
        v_spline.inferKnots();
        HSpline v_spline_copy = v_spline;
        for (size_t derivative_count = 0; derivative_count < diff_v; ++derivative_count) {
            v_spline_copy = derivative(v_spline_copy);
        }
        evaluated_row_points_derivative.row(r) = v_spline_copy.evalPt(SplineEvalMethod::DE_BOOR, v);
    }
    u_spline.setControlPts(evaluated_row_points_derivative);
    u_spline.inferKnots();
    HSpline u_spline_copy = u_spline;
    for (size_t derivative_count = 0; derivative_count < diff_u; ++derivative_count) {
        u_spline_copy = derivative(u_spline_copy);
    }
    Eigen::VectorXf final_point = u_spline_copy.evalPt(SplineEvalMethod::DE_BOOR, u);
    return final_point;
    //return HVec::Zero();
}

template<size_t Dimension>
Eigen::Matrix<float, Dimension, 2> TensorProductSurface_T<Dimension>::jacobian(double u, double v) const {
    HVec xtilde = evalPtHomogCoord(u, v);

    Eigen::Matrix<float, Dimension + 1, 2> J_xtilde;
    // TODO 2.3: Compute the Jacobian of `xtilde` in `J_xtilde`.
    J_xtilde.setZero();

    J_xtilde.col(0) = partial_xtilde(u, v, 1, 0);
    J_xtilde.col(1) = partial_xtilde(u, v, 0, 1);

    return J_xtilde.template topRows<Dimension>() / xtilde[Dimension]
            - xtilde.template head<Dimension>() * (J_xtilde.template bottomRows<1>() / std::pow(xtilde[Dimension], 2));
}

template<size_t Dimension>
auto TensorProductSurface_T<Dimension>::normal(double u, double v) const -> Vec {
    // TODO 2.3: Compute the unit surface normal at (u, v)
    Eigen::Matrix<float, Dimension, 2> J = jacobian(u, v);
    Vec normal = J.col(0).cross(J.col(1));
    return normal.normalized();
    //return Vec::Zero;
}

template<size_t Dimension>
auto TensorProductSurface_T<Dimension>::firstFundamentalForm(double u, double v) const -> M2d {
    // TODO 2.3: Calculate the first fundamental form at (u, v)
    //return M2d::Identity();
    Eigen::Matrix<float, Dimension, 2> J = jacobian(u, v);
    return J.transpose() * J;
}

template<size_t Dimension>
auto TensorProductSurface_T<Dimension>::secondFundamentalForm(double u, double v) const -> M2d {
    // Quantities needed in the chain rule terms for the perspective divide.
    HVec xtilde = evalPtHomogCoord(u, v);

    Eigen::Matrix<float, Dimension + 1, 2> J_xtilde;
    // TODO 2.4: copy your J_xtilde computation from 2.3 here.
    J_xtilde.col(0) = partial_xtilde(u, v, 1, 0);
    J_xtilde.col(1) = partial_xtilde(u, v, 0, 1);

    std::array<std::array<HVec, 2>, 2> d2xtilde;
    // TODO 2.4: calculate the "vector-valued matrix" of second partial derivatives of `xtilde`
    d2xtilde[0][0] = partial_xtilde(u, v, 2, 0);
    d2xtilde[0][1] = partial_xtilde(u, v, 1, 1);
    d2xtilde[1][0] = partial_xtilde(u, v, 1, 1);
    d2xtilde[1][1] = partial_xtilde(u, v, 0, 2);

    // The following applies the chain rule to obtain the matrix of second partial derivatives
    // of Euclidean surface point `x`.
    std::array<std::array<Vec, 2>, 2> d2x;
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            // Convert derivatives of the surface point in homogeneous
            // coordinates into derivatives of the Euclidean coordinates.
            d2x[i][j] = d2xtilde[i][j].template head<Dimension>() / xtilde[Dimension]
                      -  (J_xtilde(Dimension, i) * J_xtilde.col(j).template head<Dimension>()
                        + J_xtilde(Dimension, j) * J_xtilde.col(i).template head<Dimension>()
                        + xtilde.template head<Dimension>() * d2xtilde[i][j][Dimension]) / std::pow(xtilde[Dimension], 2)
                      + xtilde.template head<Dimension>() * (2 * (J_xtilde(Dimension, i) * J_xtilde(Dimension, j)) / std::pow(xtilde[Dimension], 3))
                    ;
        }
    }
    d2x[0][1] = d2x[1][0]; // The loop above only computes the lower triangle (exploiting symmetry)

    // TODO 2.4: evaluate the second fundamental form.
    M2d second_fundamental_form;
    second_fundamental_form(0, 0) = d2x[0][0].transpose() * normal(u, v);
    second_fundamental_form(0, 1) = d2x[0][1].transpose() * normal(u, v);
    second_fundamental_form(1, 0) = d2x[1][0].transpose() * normal(u, v);
    second_fundamental_form(1, 1) = d2x[1][1].transpose() * normal(u, v);
    return second_fundamental_form;
    //return M2d::Identity();
}

template<size_t Dimension>
auto TensorProductSurface_T<Dimension>::principalCurvatures(double u, double v) const -> V2d {
    // TODO: 2.4: solve a generalized eigenvalue problem for the principal curvatures.
    M2d FFF = firstFundamentalForm(u, v), SFF = secondFundamentalForm(u, v);
    Eigen::GeneralizedSelfAdjointEigenSolver<M2d> solver(SFF, FFF);
    V2d K = solver.eigenvalues();
    float k2 = K[0];
    K[0] = K[1];
    K[1] = k2;
    return K;
}

template<size_t Dimension>
auto TensorProductSurface_T<Dimension>::gaussianAndMeanCurvature(double u, double v) const -> V2d {
    // TODO: 2.4: evaluate (K, H) using the principal curvatures.
    V2d principal_curvatures = principalCurvatures(u, v);
    return V2d(principal_curvatures[0] * principal_curvatures[1], (principal_curvatures[0] + principal_curvatures[1]) / 2);
}

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation of all class and function templates and functions.
////////////////////////////////////////////////////////////////////////////////
template struct TensorProductSurface_T<3>;

template std::shared_ptr<TensorProductSurface> ruled_surface<3>(const NURBS<3> &spline_1, const NURBS<3> &spline_2);
template std::shared_ptr<TensorProductSurface> surface_of_revolution<3>(const NURBS<3> &spline, size_t numSegments);
