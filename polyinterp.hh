#ifndef POLYINTERP_HH
#define POLYINTERP_HH

#include <vector>
#include <Eigen/Dense>
#include "bezier.hh"

struct PolyInterp {
    enum class EvalMethod : int { AITKEN, LAGRANGE, LAGRANGE_BARYCENTRIC };
    enum class ParametrizationType : int { UNIFORM, CHORDLEN, XCOORD };

    PolyInterp(Eigen::Vector3f translation = Eigen::Vector3f::Zero()) {
        std::vector<Eigen::Vector3f> defaultInterpPts = {
            Eigen::Vector3f(-1.0, 0.0, 0.0) + translation,
            Eigen::Vector3f(0.0,  1.0, 0.0) + translation,
            Eigen::Vector3f(1.0,  0.0, 0.0) + translation
        };

        setInterpPts(defaultInterpPts);
    }

    Eigen::Vector3f evalPt(EvalMethod method, double t) const;
    void            eval  (size_t resolution, EvalMethod method, std::vector<Eigen::Vector3f> &result) const;

    void setInterpPts(const std::vector<Eigen::Vector3f> &interpPts) {
        size_t n = interpPts.size() - 1;
        p.resize(n + 1);
        p[0] = interpPts;
        m_samplePointsUpdated();
    }

    const std::vector<Eigen::Vector3f> &getInterpPts() const { return p[0]; }

    ParametrizationType getParamType() const { return m_paramType; }
    void setParamType(ParametrizationType ptype) { m_paramType = ptype; m_samplePointsUpdated(); }

    void visualizeAitken(double t, Eigen::MatrixX3f &V, Eigen::MatrixX2i &E) const;

    size_t degree() const { return p[0].size() - 1; }

    // Triangular table of evaluation points p^j_i used in Aitken's method
    mutable std::vector<std::vector<Eigen::Vector3f>> p; 
    // Parameter values at which the curve should pass through the points.
    Eigen::VectorXd t_values, w;

private:
    void m_samplePointsUpdated();
    ParametrizationType m_paramType = ParametrizationType::UNIFORM;
};

#endif /* end of include guard: POLYINTERP_HH */
