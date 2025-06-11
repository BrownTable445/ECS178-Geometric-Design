#include "spline.hh"
#include <iostream>
#include <unordered_map>
#include <unordered_set>

template<size_t Dimension>
double SplineBase<Dimension>::bsplineBasisFunction(size_t n, size_t i, double u) const {
    // TODO 1.3: Evaluate the B-Spline basis function N_i^n(u) using recursion.
    if (n == 0) {
        if (u >= knots[clamp(i - 1, (size_t)0, numKnots() - 1)] && u < knots[clamp(i, (size_t)0, numKnots() - 1)])
            return 1;
        return 0;
    }
    double l = bsplineBasisFunction(n - 1, i, u), r = bsplineBasisFunction(n - 1, i + 1, u);
    double sum = 0;
    if (l > 0)
        sum += (u - knots[clamp(i - 1, (size_t)0, numKnots() - 1)]) / (knots[clamp(n - 1 + i, (size_t)0, numKnots() - 1)] - knots[clamp(i - 1, (size_t)0, numKnots() - 1)]) * l;
    if (r > 0)
        sum += (knots[clamp(n + i, (size_t)0, numKnots() - 1)] - u) / (knots[clamp(n + i, (size_t)0, numKnots() - 1)] - knots[clamp(i, (size_t)0, numKnots() - 1)]) * r;
    return sum;
}

template<size_t Dimension>
size_t Spline_T<Dimension>::findSegmentContaining(double u) const {
    size_t I = std::distance(knots.data(),
        std::upper_bound(knots.data(), knots.data() + knots.size(), u)) - 1;
    const size_t L = numControlPts() - 1;
    return std::min(I, L - 1); // Clamp to the last curve segment [U_(L - 1), u_L]; this is needed when evaluating exactly at the upper knot bound.
}

template<size_t Dimension>
auto Spline_T<Dimension>::getControlPtsInfluencingSegment(size_t I) const ->std::vector<Vec> {
    const size_t n = degree();
    // TODO 1.1: extract the de Boor control points influencing the segment [u_I, u_{I + 1}].
    std::vector<Vec> v(n + 1, Vec::Zero());
    for (int i = (int)I - (n - 1), j = 0; i <= (int)I + 1; ++i, ++j) {
        v[j] = controlPt(i);
    }
    return v;
}

template<size_t Dimension>
typename Spline_T<Dimension>::Vec Spline_T<Dimension>::evalPt(EvalMethod method, double u) const {
    const size_t n = degree();

    size_t I = findSegmentContaining(u);
    if (n == 0) {
        if (knots.size() == 0)
            return controlPt(0);
        else if (u < knots[0])
            return controlPt(0);
        return controlPt(I + 1);
    }
    if (method == EvalMethod::DE_BOOR) {
        // TODO 1.1: implement the triangular scheme of the de Boor algorithm.
        d.resize(n + 1);
        d[0].resize(n + 1);
        d[0] = getControlPtsInfluencingSegment(I);
        for (size_t j = 1; j <= n; ++j) {
            d[j].resize(n + 1 - j);
            for (size_t i = 0; i <= n - j; ++i) {
                float alpha = (u - knots[I - n + j + i]) / (knots[I + 1 + i] - knots[I - n + j + i]);
                d[j][i] = (1 - alpha) * d[j - 1][i] + alpha * d[j - 1][i + 1];
            }
        }
        return d[n][0];
    }
    if (method == EvalMethod::BASIS) {
        // TODO 1.3: evaluate the B-Spline using the summation formula and your
        // `bsplineBasisFunction` implementation.
        Vec sum = Vec::Zero();
        for (size_t i = 0; i < numControlPts(); ++i)
            sum += controlPt(i) * bsplineBasisFunction(n, i, u);
        return sum;
    }
    throw std::runtime_error("Unimplemented method");
}

template<size_t Dimension>
void SplineBase<Dimension>::eval(size_t resolution, EvalMethod method, std::vector<Vec>& result) const {
    //result.assign(resolution, Vec::Zero());

    // TODO 1.1: evaluate the curve at `resolution` equispaced points between
    // `domainStart()` and `domainEnd()` by calling `evalPt`.
    result.clear();
    float epsilon = 0.00000001;
    for (size_t i = 0; i <= resolution; ++i) {
        result.push_back(evalPt(method, domainStart() + epsilon + (domainEnd() - epsilon - domainStart()) / resolution * i));
    }
}


template<size_t Dimension>
void Spline_T<Dimension>::visualizeDeBoor(double u, MXDf& V, Eigen::MatrixX2i& E) const {
    // TODO 1.2: Visualization of the generations of the de Boor algorithm.
    size_t n = degree();
    V.resize((n + 1) * (n + 2) / 2, 3);
    E.resize(n * (n + 1) / 2, 2);
    int vertex_index = 0, edge_index = 0;
    for (size_t i = 0; i < n + 1; ++i) {
        for (size_t j = 0; j < d[i].size(); ++j) {
            V.row(vertex_index) = d[i][j];
            if (j > 0) {
                E.row(edge_index) = Eigen::Vector2i(vertex_index, vertex_index - 1);
                ++edge_index;
            }
            ++vertex_index;
        }
    }
}

template<size_t Dimension>
void Spline_T<Dimension>::inferKnots() {
    const int n = degree();
    const int L = numControlPts() - 1;
    const int K = L + n - 1;
    //knots.setLinSpaced(K + 1, 0, 1);

    knots.resize(K + 1);


    // TODO 1.4.1 - 1.4.3: fill in `knots` by inferring the parameter values
    // using various approaches depending on the `paramType` and `repeatedEndKnots`
    // member variables.
    std::vector<float> distances(K + 1);
    for (int i = 0; i < K + 1; ++i) {
        float sum = 0;
        for (int j = i - (n - 1); j <= i; ++j) {//not j <= i + 1
            if (j >= n - 2 && j - (n - 2) + 1 <= L) {
                float distance = (controlPt(j - (n - 2)) - controlPt(j - (n - 2) + 1)).norm();
                if (paramType == ParametrizationType::CENTRIPETAL)
                    distance = sqrt(distance);
                else if (paramType == ParametrizationType::UNIFORM)
                    distance = 1;
                sum += distance;
            }
            else if (repeatedEndKnots) {
                sum = 0;
                break;
            }
        }
        distances[i] = sum;
    }
    float total_length = 0;
    for (size_t i = 0; i < K + 1; ++i)
        total_length += distances[i];
    knots = Eigen::VectorXf(K + 1);
    for (size_t i = 0; i < K + 1; ++i) {
        knots[i] = distances[i] / total_length;
        if (i > 0)
            knots[i] += knots[i - 1];
    }
    Base::validateKnots(knots);
}


template<size_t Dimension>
void Spline_T<Dimension>::insertKnot(double ubar) {
    // TODO 2.1: Calculate the new, expanded `knots` and `controlPts`
    // produced by inserting the knot `ubar`.
    //VXf knots_new = knots;
    //MXDf controlPts_new = controlPts;

    size_t n = degree();
    evalPt(EvalMethod::DE_BOOR, ubar);
    size_t I = findSegmentContaining(ubar);
    controlPts.conservativeResize(controlPts.rows() + 1, controlPts.cols());
    for (size_t i = controlPts.rows() - 1; i > I + 1; --i)
        controlPts.row(i) = controlPts.row(i - 1);
    for (size_t i = 0; i < d[1].size(); ++i)
        controlPts.row(I - n + 2 + i) = d[1][i];
    knots.conservativeResize(knots.size() + 1);
    for (size_t i = knots.size() - 1; i > I;--i)
        knots[i] = knots[i - 1];
    knots[I + 1] = ubar;

    // Apply the new curve data.
    //knots = knots_new;
    //controlPts = controlPts_new;

    Base::validateSizes();
}

template<size_t Dimension>
std::vector<BezierCurve> splineToBezierSegments(Spline_T<Dimension> spline /* copy is intentional */) {
    auto& knots = spline.knots;
    const size_t n = spline.degree();

    // TODO 2.2: extract the Bézier curve representation of each polynomial piece
    // of the B-Spline `spline`. This should be done by inserting each knot up
    // to multiplicity `m` and then extracting the de Boor points influencing
    // each nonempty curve segment.
    size_t L = spline.numControlPts() - 1;
    std::unordered_map<float, size_t> u_m;
    std::unordered_set<float> unique_knots;
    for (size_t i = 0; i < knots.size(); ++i) {
        ++u_m[knots[i]];
        if (i >= n - 1 && i <= L)
            unique_knots.insert(knots[i]);
    }

    for (auto i : unique_knots) {
        for (size_t j = u_m[i]; j < n; ++j) {
            spline.insertKnot(i);
        }
    }

    L = spline.numControlPts() - 1;

    std::vector<BezierCurve> result;
    for (size_t i = n - 1; i < L; ++i) {
        std::vector<Eigen::Vector3f> v(n + 1);
        for (size_t j = 0; j <= n; ++j) {
            v[j] = spline.controlPt(i - (n - 1) + j);
        }
        result.emplace_back(v);
    }

    return result;
}

////////////////////////////////////////////////////////////////////////////////
// NURBS functionality (using homogeneous coordinates)
////////////////////////////////////////////////////////////////////////////////
template<size_t Dimension>
typename NURBS<Dimension>::Vec NURBS<Dimension>::controlPt(size_t i) const {
    // TODO 4: get the Euclidean coordinates corresponding to the homogenizedf
    // coordinate vector `m_homogeneousSpline.controlPt(i)`.
    return (m_homogeneousSpline.controlPt(i) / weight(i)).head(Dimension);
}

template<size_t Dimension>
typename NURBS<Dimension>::MXDf NURBS<Dimension>::getControlPts() const {
    // TODO 4: get the Euclidean coordinates corresponding to all of the
    // homogenized coordinate vectors stored in the rows of
    // `m_homogeneousSpline.getControlPts()`.
    size_t num_control_pts = m_homogeneousSpline.numControlPts();
    Eigen::MatrixXf euclidean_control_pts = Eigen::MatrixXf::Zero(num_control_pts, Dimension);
    for (size_t i = 0; i < num_control_pts; ++i) {
        euclidean_control_pts.row(i) = controlPt(i);
    }
    return euclidean_control_pts;
}

template<size_t Dimension>
void NURBS<Dimension>::setControlPt(size_t i, const Vec& v) {
    // TODO 4: set the Euclidean coordinates of control point `i`
    // while leaving its weight unchanged.
    Eigen::VectorXf homogeneous_v = Eigen::VectorXf::Zero(Dimension + 1);
    homogeneous_v.head(Dimension) = v;
    homogeneous_v[Dimension] = weight(i);
    m_homogeneousSpline.setControlPt(i, homogeneous_v);
}

template<size_t Dimension>
void NURBS<Dimension>::setWeight(size_t i, float w) {
    // TODO 4: set the weight of control point `i` while leaving its Euclidean
    // coordinates unchanged.
    m_homogeneousSpline.controlPts(i) /= weight(i);
    m_homogeneousSpline.controlPts(i) *= w;
}

template<size_t Dimension>
typename NURBS<Dimension>::Vec NURBS<Dimension>::evalPt(EvalMethod method, double u) const {
    // TODO 4: evaluate s(u) by first evaluating a point in homogenized
    // coordinates using `m_homogeneousSpline.evalPt` and then doing a
    // perspective divide.
    Eigen::VectorXf evaluated_pt = m_homogeneousSpline.evalPt(method, u);
    return (evaluated_pt / evaluated_pt[Dimension]).head(Dimension);
}

// Calculate the derivative of a spline.
template<size_t Dimension>
Spline_T<Dimension> derivative(const Spline_T<Dimension> &s) {
    using VXf  = typename Spline_T<Dimension>::VXf;
    using MXDf = typename Spline_T<Dimension>::MXDf;

    // HW5 TODO 2.1: calculate the one-lower degree spline that is the derivative of `s`.
    // Note: if you've copied this method into your old `spline.cc`, be sure also to copy
    // the two explicit template instantiations of this method at the bottom of the file.
    size_t n = s.degree();
    VXf  u = s.knots;
    MXDf d = s.getControlPts();

    if (n == 0) {
        d.setZero();
        return s;
    }

    size_t K = u.size() - 1;
    size_t L = d.rows() - 1;
    MXDf d2(L, Dimension);
    for (size_t i = 0; i < L; ++i) {
        d2.row(i) = n * (d.row(i + 1) - d.row(i)) / (u[i + n] - u[i]);
    }

    return Spline_T<Dimension>(n - 1, u.segment(1, K - 1), d2);
}

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation of all class and function templates and functions.
////////////////////////////////////////////////////////////////////////////////
template struct SplineBase<3>;
template struct Spline_T<3>;
template struct SplineBase<4>;
template struct Spline_T<4>;
template struct NURBS<3>;
template std::vector<BezierCurve> splineToBezierSegments<3>(Spline_T<3> spline);

template Spline_T<3> derivative(const Spline_T<3> &s);
template Spline_T<4> derivative(const Spline_T<4> &s);
