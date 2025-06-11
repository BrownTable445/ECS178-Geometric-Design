#ifndef BEZIER_HH
#define BEZIER_HH

#include <vector>
#include <Eigen/Dense>

inline static Eigen::MatrixXf pointsToMatrixRows(const std::vector<Eigen::Vector3f> pts) {
    Eigen::MatrixXf result(pts.size(), 3);
    for (size_t i = 0; i < pts.size(); ++i)
        result.row(i) = pts[i].transpose();
    return result;
}

inline static std::vector<Eigen::Vector3f> pointsFromMatrixRows(const Eigen::MatrixXf &B) {
    if (B.cols() != 3) throw std::runtime_error("Unexpected shape of B");

    std::vector<Eigen::Vector3f> result(B.rows());
    for (size_t i = 0; i < result.size(); ++i)
        result[i] = B.row(i).transpose();
    return result;
}


struct BoundingBox {
    // Default constructor creates an empty bounding box
    BoundingBox() : minCorner(Eigen::Vector3f::Constant(std::numeric_limits<float>::max())),
                    maxCorner(Eigen::Vector3f::Constant(std::numeric_limits<float>::min())) { }


    BoundingBox(const Eigen::Vector3f &minCorner, const Eigen::Vector3f &maxCorner)
        : minCorner(minCorner), maxCorner(maxCorner) { }

    Eigen::Vector3f minCorner = Eigen::Vector3f::Constant(std::numeric_limits<float>::max()),
                    maxCorner = Eigen::Vector3f::Constant(std::numeric_limits<float>::min());

    bool overlaps(const BoundingBox &b) const;
    float diameter()     const { return (maxCorner - minCorner).norm(); }
    Eigen::Vector3f center() const { return 0.5 * (minCorner + maxCorner); }

    BoundingBox unionWith(const BoundingBox &b) const { return BoundingBox{minCorner.cwiseMin(b.minCorner), maxCorner.cwiseMax(b.maxCorner)}; }
};

inline std::ostream &operator<<(std::ostream &os, BoundingBox bb) {
    os << "[(" << bb.minCorner.transpose() << "), (" << bb.maxCorner.transpose() <<  ")]";
    return os;
}

struct BezierCurve {
    enum class EvalMethod : int { DE_CASTELJAU, BERNSTEIN, HORNER };

    BezierCurve(Eigen::Vector3f translation = Eigen::Vector3f::Zero()) {
        std::vector<Eigen::Vector3f> defaultControlPts = {
            Eigen::Vector3f(-1.0, 0.0, 0.0) + translation,
            Eigen::Vector3f(0.0,  1.0, 0.0) + translation,
            Eigen::Vector3f(1.0,  0.0, 0.0) + translation
        };

        setControlPts(defaultControlPts);
    }

    BezierCurve(const std::vector<Eigen::Vector3f> &controlPts) {
        setControlPts(controlPts);
    }

    void setControlPts(const std::vector<Eigen::Vector3f> &controlPts) {
        size_t n = controlPts.size() - 1;
        b.resize(n + 1);
        b[0] = controlPts;
    }

    const std::vector<Eigen::Vector3f> &getControlPts() const {
        checkPolygon();
        return b[0];
    }

    std::vector<Eigen::Vector3f> &getControlPts() {
        checkPolygon();
        return b[0];
    }

    void checkPolygon() const {
        if (b.size() < 1) throw std::runtime_error("Control polygon was not set");
    }

    Eigen::Vector3f evalPt(EvalMethod method, double t) const;
    void            eval  (size_t resolution, EvalMethod method, std::vector<Eigen::Vector3f> &result) const;

    // Get the unit tangent and normal
    Eigen::Vector3f probeCurve(double t, Eigen::Vector3f &tangent, Eigen::Vector3f &normal) const;

    void visualizeDeCasteljau(double t, Eigen::MatrixX3f &V, Eigen::MatrixX2i &E) const;

    size_t degree() const { checkPolygon(); return b[0].size() - 1; }

    void setDegree(size_t degree);
    void reduceDegree();

    void subdivide(double t, BezierCurve &c1, BezierCurve &c2) const;

    // Drag the curve point c(t) to p by applying a minimal change to the control polygon.
    void drag(float t, const Eigen::Vector3f &p);

    void approximate(const Eigen::VectorXd &t_values, const std::vector<Eigen::Vector3f> &pts, double smoothingWeight);

    BoundingBox boundingBox() const;
    // Single member variable:
    // Triangular table of evaluation points b^j_i used in de Casteljau algorithm
    mutable std::vector<std::vector<Eigen::Vector3f>> b; // triangular table of evaluation points b^j_i used in de Casteljau algorithm
};

void getIntersections(const BezierCurve &c1, const BezierCurve &c2, std::vector<Eigen::Vector3f> &result, float tol);

void mergeDuplicateIntersections(std::vector<Eigen::Vector3f> &ipoints, float tol);

#endif /* end of include guard: BEZIER_HH */
