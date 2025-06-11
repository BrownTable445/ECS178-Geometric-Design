////////////////////////////////////////////////////////////////////////////////
// spline.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  B-spline interpolation and manipulation code for ECS 178.
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
//  Company:  University of California, Davis
//  Created:  11/13/2022 20:27:51
*///////////////////////////////////////////////////////////////////////////////
#ifndef SPLINE_HH
#define SPLINE_HH

#include <memory>
#include "bezier.hh"
#include <iostream>

template<typename T>
T clamp(T x, T a, T b) {
    return std::max(a, std::min(b, x));
}

enum class SplineEvalMethod : int { DE_BOOR, BASIS };
enum class SplineParametrizationType : int { UNIFORM, CHORDLEN, CENTRIPETAL };

template<size_t Dimension>
struct SplineBase {
    using EvalMethod = SplineEvalMethod;
    using ParametrizationType = SplineParametrizationType;

    using Vec = Eigen::Matrix<float, Dimension, 1>;
    using VXf     = Eigen::VectorXf;
    using MXDf    = Eigen::Matrix<float, Eigen::Dynamic, Dimension>;

    virtual size_t numControlPts() const = 0;

    size_t numKnots()   const { return getKnots().size(); }
    // Curve is defined over the parameter interval [u_n - 1, u_L]
    float domainStart() const { return getKnots()[degree() - 1]; }
    float domainEnd  () const { return getKnots()[numControlPts() - 1]; }

    virtual Vec controlPt(size_t    i) const = 0;
    virtual void setControlPt(size_t i, const Vec &p) = 0;

    virtual void setControlPts(const MXDf &P) = 0;
    virtual MXDf getControlPts() const = 0;

    void make2D() {
        if (Dimension < 3) return;
        for (size_t i = 0; i < numControlPts(); ++i)
            controlPt(i)[2] = 0;
    }

    // Must be overridden in the NURBS case
    virtual const VXf &getKnots()  const { return knots; }
    virtual       VXf &getKnots()        { return knots; }
    virtual size_t        degree() const { return m_degree; }
    virtual size_t       &degree()       { return m_degree; }

    virtual Vec evalPt(EvalMethod method, double u) const = 0;
    void eval(size_t resolution, EvalMethod method, std::vector<Vec> &result) const;

    virtual void inferKnots() = 0;
    virtual void insertKnot(double ubar) = 0;

    double bsplineBasisFunction(size_t n, size_t i, double u) const;

    static void validateKnots(const VXf &knots) {
        // Sanity check: validate ordering
        for (size_t i = 1; i < knots.size(); ++i) {
            if (knots[i] < knots[i - 1])
                throw std::runtime_error("Knots are not sorted.");
        }
    }

    void validateSizes() const {
        size_t K = numKnots() - 1;
        size_t L = numControlPts() - 1;

        if (K != L + degree() - 1) {
            std::cerr << "Mismatched control point and knot sizes:" << std::endl
                << "    L = " << L << ", K = " << K << ", n = " << degree() << std::endl;
            throw std::runtime_error("Size mismatch");
        }
    }

    virtual ~SplineBase() { }

    // Both ordinary B-Splines and NURBS share the same knot/degree representation.
    VXf  knots;
    ParametrizationType paramType = ParametrizationType::UNIFORM;
    bool repeatedEndKnots = true;
protected:
    size_t m_degree = 3;
};

template<size_t Dimension>
struct Spline_T : SplineBase<Dimension> {
    using Base = SplineBase<Dimension>;
    using Vec     = typename Base::Vec;
    using MXDf    = typename Base::MXDf;
    using VXf  = Eigen::VectorXf;

    using EvalMethod          = typename Base::EvalMethod;
    using ParametrizationType = typename Base::ParametrizationType;

    Spline_T(const Vec &t = Vec::Zero()) {
        controlPts.setZero(5, Dimension);
        controlPts.col(0).setLinSpaced(-2, 2);
        controlPts.col(1) = VXf::LinSpaced(5, 0.0, 2 * M_PI).array().sin().matrix();
        controlPts.rowwise() += t.transpose();
        inferKnots();
    }

    Spline_T(size_t n, const VXf &u, const MXDf &d) {
        Base::m_degree = n;
        knots = u;
        controlPts = d;
        Base::validateSizes();
    }

    size_t numControlPts() const override { return controlPts.rows();  }

    void setControlPts(const MXDf &P) override { controlPts = P; }
    Vec controlPt(size_t i) const override { return controlPts.row(i).transpose(); }
    void setControlPt(size_t i, const Vec &v) override { controlPts.row(i) = v.transpose(); }

    MXDf getControlPts() const override { return controlPts; }
    std::vector<Vec> getControlPtsInfluencingSegment(size_t I) const;

    void inferKnots() override;
    void insertKnot(double ubar) override;

    size_t findSegmentContaining(double u) const;
    Vec evalPt(EvalMethod method, double u) const override;

    void visualizeDeBoor(double u, MXDf &V, Eigen::MatrixX2i &E) const;

    // Triangular table of evaluation points d^j_i for the de Boor algorithm.
    mutable std::vector<std::vector<Vec>> d; 

    void info() const {
        std::cout << "Degree " << degree() << " spline with" << std::endl
                  << "\t" << numKnots() << " knots: " << knots.transpose() << std::endl
                  << "\t" << numControlPts() << " controlPts: " << std::endl << controlPts << std::endl;
    }

    // Parameter values of the spline
    MXDf controlPts;
    using Base::numKnots;
    using Base::knots;
    using Base::degree;
    using Base::paramType;
    using Base::repeatedEndKnots;
    using Base::bsplineBasisFunction;
};

template<size_t Dimension>
std::vector<BezierCurve> splineToBezierSegments(Spline_T<Dimension> spline);

using Spline = Spline_T<3>; // Non-rational B-Spline in 3D.

// NURBS curves can be implemented using ordinary splines in homogeneous
// coordinates.
template<size_t Dimension>
struct NURBS : public SplineBase<Dimension> {
    using Base = SplineBase<Dimension>;
    using Vec     = typename Base::Vec;
    using VXf     = typename Base::VXf;
    using MXDf    = typename Base::MXDf;

    using EvalMethod          = typename Base::EvalMethod;
    using ParametrizationType = typename Base::ParametrizationType;

    NURBS() {
        m_homogeneousSpline.controlPts.template rightCols<1>().setOnes();
    }

    NURBS(size_t degree, const VXf &knots, const MXDf &controlPts, const VXf &weights) {
        m_homogeneousSpline.degree() = degree;
        m_homogeneousSpline.controlPts.resize(controlPts.rows(), Dimension + 1);
        m_homogeneousSpline.controlPts.template leftCols<Dimension>() = weights.asDiagonal() * controlPts;
        m_homogeneousSpline.controlPts.template rightCols<1>() = weights;
        m_homogeneousSpline.knots = knots;
    }

    // Promotion from non-rational B-Spline
    NURBS(const Spline_T<Dimension> &ordinarySpline) {
        Base::paramType = ordinarySpline.paramType;
        Base::repeatedEndKnots = ordinarySpline.repeatedEndKnots;

        m_homogeneousSpline.controlPts.resize(ordinarySpline.numControlPts(), Dimension + 1);

        m_homogeneousSpline.controlPts.resize(ordinarySpline.numControlPts(), Dimension + 1);
        m_homogeneousSpline.controlPts.template leftCols<Dimension>() = ordinarySpline.getControlPts();
        m_homogeneousSpline.controlPts.template rightCols<1>().setOnes();
        m_homogeneousSpline.knots = ordinarySpline.knots;
        m_homogeneousSpline.degree() = ordinarySpline.degree();

        Base::validateSizes();
    }

    using HSpline = Spline_T<Dimension + 1>;

    void setHomogeneousControlPts(const typename HSpline::MXDf &Pw) {
        m_homogeneousSpline.setControlPts(Pw);
    }

    void setControlPts(const MXDf &P) override {
        throw std::runtime_error("Dangerous; avoid.");
    }

    Vec     controlPt(size_t i)         const override;
    void setControlPt(size_t i, const Vec &v) override;

    virtual const VXf &getKnots() const override { return m_homogeneousSpline.getKnots(); }
    virtual       VXf &getKnots()       override { return m_homogeneousSpline.getKnots(); }

    size_t  degree() const override { return m_homogeneousSpline.degree(); }
    size_t &degree()       override { return m_homogeneousSpline.degree(); }

    float   weight(size_t i         ) const { return m_homogeneousSpline.controlPts(i, Dimension); }
    void setWeight(size_t i, float w);

    auto weights() const { return m_homogeneousSpline.controlPts.template rightCols<1>(); }

    MXDf getControlPts() const override;

    Vec evalPt(EvalMethod method, double u) const override;

    void inferKnots() override { m_homogeneousSpline.inferKnots(); }
    void insertKnot(double ubar) override { m_homogeneousSpline.insertKnot(ubar); }

    void visualizeDeBoor(double u, MXDf &V, Eigen::MatrixX2i &E) const;

    size_t numControlPts() const override { return m_homogeneousSpline.numControlPts();  }

    const HSpline &homogeneousSpline() const { return m_homogeneousSpline; }

private:
    using Base::knots;
    HSpline m_homogeneousSpline;
};

// Get a NURBS version of `genericSpline`, promoting to a new copy if necessary.
template<size_t Dimension>
std::shared_ptr<NURBS<Dimension>> NURBSForSpline(const std::shared_ptr<SplineBase<Dimension>> &genericSpline) {
    auto nurbsSpline = std::dynamic_pointer_cast<NURBS<Dimension>>(genericSpline);
    if (!nurbsSpline) {
        auto ordinarySpline = std::dynamic_pointer_cast<Spline_T<Dimension>>(genericSpline);
        if (!ordinarySpline) throw std::logic_error("Generic spline must be either ordinary or NURBS");
        nurbsSpline = std::make_shared<NURBS<Dimension>>(*ordinarySpline);
    }

    return nurbsSpline;
}

template<size_t Dimension>
NURBS<Dimension> NURBSCircle(float radius, size_t numSegments, size_t keptSegments = -1) {
    if (numSegments < 3) throw std::runtime_error("Segments must be >= 3");
    keptSegments = std::min(keptSegments, numSegments);

    // http://julianpanetta.com/teaching/GMSlides/17-NURBS-deck.html#/how-to-represent-a-perfect-circle/7
    float theta = M_PI / numSegments; // half the subtended angle
    float w_1 = std::cos(theta);
    Eigen::Rotation2D<float> R(2 * theta);

    using MXDf = typename NURBS<Dimension>::MXDf;
    MXDf controlPts;
    controlPts.setZero(1 + keptSegments * 2, Dimension); // Two control points contributed by each segment (after accounting for sharing)
    Eigen::VectorXf weights(controlPts.rows());
    Eigen::VectorXf knots(2 * (keptSegments + 1)); // Two knots contributed by each segment (after account for sharing)

    Eigen::Vector2f b0(radius, 0),
                    b2(radius * cos(2 * theta), radius * sin(2 * theta));

    Eigen::Vector2f vPerp(b2[1] - b0[1], b0[0] - b2[0]); // 90 degree **cw** rotation
    Eigen::Vector2f b1 = 0.5 * ((b0 + b2) + vPerp * tan(theta));

    // Rotate to make a 4-segment circle starts in the bottom-right quadrant
    // (so that when `keptSegments = 2` we get the right semicircle that can
    //  be revolved around the y axis).
    b0 = R.inverse() * b0;
    b1 = R.inverse() * b1;
    b2 = R.inverse() * b2;

    auto controlPtsXY = controlPts.template leftCols<2>();

    controlPtsXY.row(0) = b0.transpose();
    controlPtsXY.row(1) = b1.transpose();
    controlPtsXY.row(2) = b2.transpose();
    weights[0] = 1;
    weights[1] = w_1;
    weights[2] = 1;
    knots[0] = 0;
    knots[1] = 0;
    knots[2] = 1.0 / numSegments;
    knots[3] = 1.0 / numSegments;

    for (size_t i = 1; i < keptSegments; ++i) {
        controlPtsXY.row(2 * i + 1) = R * controlPtsXY.row(2 * (i - 1) + 1).transpose();
        controlPtsXY.row(2 * i + 2) = R * controlPtsXY.row(2 * (i - 1) + 2).transpose();
        weights[2 * i + 1] = w_1;
        weights[2 * i + 2] = 1;

        knots[2 * (i + 1) + 0] = (1.0 / numSegments) + knots[2 * i + 0];
        knots[2 * (i + 1) + 1] = (1.0 / numSegments) + knots[2 * i + 1];
    }

    NURBS<Dimension> result(2, knots, controlPts, weights);
    result.validateSizes();
    result.validateKnots(knots);
    result.homogeneousSpline().info();
    return result;
}

// Calculate the derivative of a spline.
template<size_t Dimension>
Spline_T<Dimension> derivative(const Spline_T<Dimension> &s);

#endif /* end of include guard: SPLINE_HH */
