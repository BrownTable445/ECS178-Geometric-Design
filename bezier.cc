#include "bezier.hh"
#include <iostream>

Eigen::Vector3f deCasteljauRecursive(std::vector<Eigen::Vector3f>& b0, double t, size_t j, size_t i) {
    if (j == 0)
        return b0[i];
    return (1 - t) * deCasteljauRecursive(b0, t, j - 1, i) + t * deCasteljauRecursive(b0, t, j - 1, i + 1);
}

long long factorial(size_t n) {
    if (n <= 1)
        return 1;
    return n * factorial(n - 1);
}

double getBernsteinBasis(double t, size_t j, size_t i) {
    return factorial(j) / (factorial(j - i) * factorial(i)) * pow(t, i) * pow(1 - t, j - i);
}

Eigen::Vector3f BezierCurve::evalPt(EvalMethod method, double t) const {
    if (method == EvalMethod::DE_CASTELJAU) {
        // TODO: Problem 1.1
        // Implement de Casteljau's algorithm, filling in the point table
        // stored in member variable `b`.
        for (int j = 1; j < b[0].size(); ++j) {
            b[j].resize(b[0].size() - j);
            for (int i = 0; i < b[0].size() - j; ++i) {
                b[j][i] = (1 - t) * b[j - 1][i] + t * b[j - 1][i + 1];
            }
        }
        return Eigen::Vector3f(b[b[0].size() - 1][0]);
    }
    if (method == EvalMethod::BERNSTEIN) {
        // TODO: Problem 2.1
        // Evaluate the Bézier curve using the Bernstein basis.
        Eigen::Vector3f accumulator = Eigen::Vector3f::Zero();
        for (size_t i = 0; i <= degree(); ++i) {
            accumulator += b[0][i].transpose() * getBernsteinBasis(t, degree(), i);
        }
        return accumulator;
    }
    if (method == EvalMethod::HORNER) {
        // TODO: Problem 2.2
        // Evaluate the Bézier curve efficiently using the Horner approach.
        size_t n = degree(); //n is the maximum INDEX in b[0]
        double s = 1;
        Eigen::Vector3f accumulator = b[0].back();
        for (int i = n - 1; i >= 0; --i) {
            s *= (i + 1.0) / (n - i) * (1 - t);
            accumulator = (accumulator * t) + b[0][i] * s;
        }
        return accumulator;
    }
    throw std::runtime_error("Unimplemented method");
}

void BezierCurve::eval(size_t resolution, EvalMethod method, std::vector<Eigen::Vector3f>& result) const {
    // TODO: Problem 1.2
    // Evaluate the discrete polyline approximation to this curve in `result`,
    // using `resolution` points evenly spaced over parameter interval [0, 1].
    std::vector<Eigen::Vector3f> evaluated_points(resolution);
    for (size_t i = 0; i < resolution; ++i) {
        evaluated_points[i] = evalPt(method, 1.0 / (resolution - 1) * i);
    }
    result = evaluated_points;
}

void BezierCurve::visualizeDeCasteljau(double t, Eigen::MatrixX3f& V, Eigen::MatrixX2i& E) const {
    // TODO: Problem 1.3
    // Fill `V` and `E` with an indexed edge set representation of
    // the lines/points generated during the de Casteljau process.
    size_t n = degree();
    unsigned int total_points = n * (n + 1) / 2, total_edges = (n - 1) * n / 2;
    V.resize(total_points, 3);
    E.resize(total_edges, 2);
    evalPt(EvalMethod::DE_CASTELJAU, t);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n - i; ++j) {
            size_t current_index = total_points - (n - i) * (n - i + 1) / 2 + j;
            V.row(current_index) = b[i + 1][j];
            if (j > 0) {
                E.row(total_edges - (n - i - 1) * (n - i) / 2 + j - 1) = Eigen::Vector2i(current_index - 1, current_index);
            }
        }
    }
}

// Get the tangent from the last evaluation
Eigen::Vector3f BezierCurve::probeCurve(double t, Eigen::Vector3f& tangent, Eigen::Vector3f& normal) const {
    // TODO: Problem 4.1
    // Evaluate the position, unit tangent, and unit normal of the curve at
    // parameter `t`.
    Eigen::Vector3f position = evalPt(EvalMethod::DE_CASTELJAU, t);
    tangent = (b[b[0].size() - 2][1] - b[b[0].size() - 2][0]).normalized();
    Eigen::Vector3f c_2nd_derivative = degree() * (degree() - 1) * (b[degree() - 2][2] - 2 * b[degree() - 2][1] + b[degree() - 2][0]);
    normal = (c_2nd_derivative - tangent * tangent.dot(c_2nd_derivative)).normalized();
    return position;
}

void BezierCurve::setDegree(size_t new_degree) {
    if (new_degree < 1) throw std::runtime_error("Invalid degree");
    size_t npts = new_degree + 1;

    //const auto &oldPts = getControlPts();
    auto oldPts = b[0];
    std::vector<Eigen::Vector3f> newPts = oldPts;


    if (new_degree > degree()) {
        while (newPts.size() < npts) {
            newPts.emplace_back(oldPts.back());
            for (size_t i = 1; i < oldPts.size(); ++i) {
                float alpha = i / float(oldPts.size());
                newPts[i] = alpha * oldPts[i - 1] + (1 - alpha) * oldPts[i];
            }
            oldPts = newPts;
        }
    }
    else if (new_degree < degree()) {
        while (oldPts.size() > npts) {
            size_t n = oldPts.size() - 2; //the degree we want to downsize to
            Eigen::MatrixXf M(n + 2, n + 1);
            M.setZero();
            M(0, 0) = 1;
            M(n + 1, n) = 1;
            for (size_t j = 1; j < n + 1; ++j) {
                float alpha = j / (n + 1.0f);
                M(j, j - 1) = alpha;
                M(j, j) = 1 - alpha;
            }
            Eigen::MatrixXf A = M.transpose() * M;
            Eigen::MatrixX3f x_tilde(n + 2, 3);
            for (size_t i = 0; i < n + 2; ++i)
                x_tilde.row(i) = oldPts[i];
            Eigen::MatrixX3f x(n + 1, 3);
            x = A.llt().solve(M.transpose() * x_tilde);
            oldPts.resize(n + 1);
            for (size_t i = 0; i < n + 1; ++i) {
                oldPts[i] = x.row(i);
            }
        }
        newPts = oldPts;
    }

    setControlPts(newPts);
}

void BezierCurve::subdivide(double t, BezierCurve& c1, BezierCurve& c2) const {
    // TODO: Problem 3.2
    // Split this curve at the parameter `t` into curves `c` and `c2`.

    evalPt(EvalMethod::DE_CASTELJAU, t);
    std::vector<Eigen::Vector3f> left_control_points(b[0].size()), right_control_points(b[0].size());
    for (int i = 0; i < b[0].size(); ++i) {
        left_control_points[i] = b[i][0];
        right_control_points[i] = b[degree() - i][i];
    }

    c1 = BezierCurve(left_control_points);
    c2 = BezierCurve(right_control_points);
}

BoundingBox BezierCurve::boundingBox() const {
    //BoundingBox bb;
    // TODO: Problem 4.2
    // Compute the bounding box of the control points for this curve.
    Eigen::Vector3f min = b[0][0], max = b[0][0];
    for (auto& i : b[0]) {
        if (i(0, 0) < min(0, 0))
            min(0, 0) = i(0, 0);
        if (i(0, 0) > max(0, 0))
            max(0, 0) = i(0, 0);
        if (i(1, 0) < min(1, 0))
            min(1, 0) = i(1, 0);
        if (i(1, 0) > max(1, 0))
            max(1, 0) = i(1, 0);
        if (i(2, 0) < min(2, 0))
            min(2, 0) = i(2, 0);
        if (i(2, 0) > max(2, 0))
            max(2, 0) = i(2, 0);
    }
    return BoundingBox(min, max);
}

bool BoundingBox::overlaps(const BoundingBox& b) const {
    // TODO: Problem 4.2
    // Return whether this bounding box overlaps the other bounding box `b`.
    float this_x_length = this->maxCorner(0, 0) - this->minCorner(0, 0);
    float this_y_length = this->maxCorner(1, 0) - this->minCorner(1, 0);
    float this_z_length = this->maxCorner(2, 0) - this->minCorner(2, 0);
    float other_x_length = b.maxCorner(0, 0) - b.minCorner(0, 0);
    float other_y_length = b.maxCorner(1, 0) - b.minCorner(1, 0);
    float other_z_length = b.maxCorner(2, 0) - b.minCorner(2, 0);
    bool x_collision = this->minCorner(0, 0) + this_x_length >= b.minCorner(0, 0) && b.minCorner(0, 0) + other_x_length >= this->minCorner(0, 0);
    bool y_collision = this->minCorner(1, 0) + this_y_length >= b.minCorner(1, 0) && b.minCorner(1, 0) + other_y_length >= this->minCorner(1, 0);
    bool z_collision = this->minCorner(2, 0) + this_z_length >= b.minCorner(2, 0) && b.minCorner(2, 0) + other_z_length >= this->minCorner(2, 0);
    return x_collision && y_collision && z_collision;
}

void getIntersections(const BezierCurve& c1, const BezierCurve& c2, std::vector<Eigen::Vector3f>& result, float tol) {
    // TODO: Problem 4.2
    // *Append* approximations of *all* intersections between curves `c1` and
    // `c2` into the result collection `result`.
    if (c1.boundingBox().overlaps(c2.boundingBox())) {
        BoundingBox temp = c1.boundingBox().unionWith(c2.boundingBox());
        if (temp.diameter() < tol)
            result.emplace_back(temp.center());
        else {
            BezierCurve c1_p1, c1_p2, c2_p1, c2_p2;
            c1.subdivide(0.5, c1_p1, c1_p2);
            c2.subdivide(0.5, c2_p1, c2_p2);
            getIntersections(c1_p1, c2_p1, result, tol);
            getIntersections(c1_p1, c2_p2, result, tol);
            getIntersections(c1_p2, c2_p1, result, tol);
            getIntersections(c1_p2, c2_p2, result, tol);
        }
    }
}

void mergeDuplicateIntersections(std::vector<Eigen::Vector3f>& ipoints, float tol) {
    // TODO: Bonus
    // Merge all redundant, nearly equal copies of intersection points
    // appearing in `ipoints`.
    for (int i = 0; i < ipoints.size(); ++i) {
        for (int j = 0; j < ipoints.size(); ++j) {
            if (i != j && (ipoints[i] - ipoints[j]).norm() < 2 * tol) {
                ipoints.emplace_back((ipoints[i] + ipoints[j]) / 2);
                ipoints.erase(ipoints.begin() + i);
                ipoints.erase(ipoints.begin() + j);
                i = 0;
                break;
            }
        }
    }
}

void BezierCurve::reduceDegree() {
    // TODO: HW3 Problem 2.1
    // Implement one step of degree reduction, computing the control points of
    // the one-degree-lower Bézier curve whose *elevation* best approximates
    // the current control points.
    Eigen::MatrixXf X_tilde = pointsToMatrixRows(b[0]);
    setControlPts(pointsFromMatrixRows(X_tilde));
}

void BezierCurve::drag(float t, const Eigen::Vector3f& p) {
    // TODO: HW3 Problem 2.2
    // Calculate and apply the minimum-norm displacement to the current control points
    // such that the curve passes through `p` at parameter `t` (i.e., `c(t) = p`).
    Eigen::VectorXf M_tilde(degree() - 1);
    for (size_t i = 1; i < degree(); ++i) {
        M_tilde[i - 1] = getBernsteinBasis(t, degree(), i);
    }
    Eigen::Vector3f c_t = evalPt(EvalMethod::HORNER, t);
    Eigen::MatrixX3f d(degree() - 1, 3);
    //d = (M_tilde.transpose() / M_tilde.squaredNorm()) * (p - c_t);
    for (size_t i = 0; i < degree() - 1; ++i) {
        d.row(i) = (M_tilde[i] / M_tilde.squaredNorm()) * (p - c_t);
    }
    for (size_t i = 1; i < degree(); ++i) {
        b[0][i] += d.row(i - 1).transpose();
    }
}

void BezierCurve::approximate(const Eigen::VectorXd& t_values, const std::vector<Eigen::Vector3f>& pts, double smoothingWeight) {
    size_t m = t_values.size(), n = degree();
    if (pts.size() != m) throw std::runtime_error("Time and point size mismatch");

    // TODO: HW3 Problem 3.1
    // Solve for the Bézier control point positions that represent the best-fit
    // polynomial curve for the data {(t_values[i], pts[i])}.
    Eigen::MatrixXf B = pointsToMatrixRows(b[0]);

    Eigen::MatrixXf V(m, n + 1);

    for (size_t r = 0; r < m; ++r) {
        for (size_t j = 0; j < n + 1; ++j) {
            V(r, j) = getBernsteinBasis(t_values[r], n, j);
        }
    }

    //Eigen::MatrixXf A = V.transpose() * V;
    Eigen::MatrixX3f P(pts.size(), 3);

    for (size_t i = 0; i < pts.size(); ++i) {
        P.row(i) = pts[i];
    }

    //B = A.llt().solve(V.transpose() * P);

    //setControlPts(pointsFromMatrixRows(B));

    // TODO: HW3 Problem 3.2
    // Incorporate a smoothing regularization term into the least-squares approximation
    // using the weight `smoothingWeight`.

    Eigen::MatrixXf S(n - 1, n + 1);

    S.setZero();

    for (size_t i = 0; i < n - 1; ++i) {
        S(i, i) = 1;
        S(i, i + 1) = -2;
        S(i, i + 2) = 1;
    }

    Eigen::MatrixXf A = V.transpose() * V + smoothingWeight * S.transpose() * S;

    B = A.llt().solve(V.transpose() * P);

    setControlPts(pointsFromMatrixRows(B));
}
