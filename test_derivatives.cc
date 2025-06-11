#include "spline.hh"
#include <algorithm>
#include <numeric>
#include <cstdlib>

using Spline = Spline_T<3>;

int main(int argc, const char *argv[]) {
    size_t numTests = 1000;
    float maxError = 0;
    for (size_t i = 0; i < numTests; ++i) {
        size_t n = 1 + rand() % 5;
        size_t numPolynomialPieces = 1 + rand() % 5;
        size_t L = numPolynomialPieces + n - 1;
        size_t K = L + n - 1;

        // Generate a rand polyline.
        Spline_T<3>::MXDf controlPts;
        controlPts.setZero(L + 1, 3);
        controlPts.col(0).setLinSpaced(0, 1);
        controlPts.col(1).setLinSpaced(0, 1);
        controlPts += 0.5 * Spline_T<3>::MXDf::Random(L + 1, 3);

        // Generate knots using rand spacings in [1/4, 3/4]
        Eigen::VectorXf knots = (0.5 + Eigen::VectorXf::Random(K + 1).array() / 4).matrix();
        knots[0] = 0;
        std::partial_sum(knots.data(), knots.data() + knots.size(), knots.data());
        knots /= knots[knots.size() - 1];

        Spline s(n, knots, controlPts);

        float u = 0;

        float fd_eps = 0.5e-3;

        // Generate a rand evaluation point such that the central finite
        // difference stencil is within bounds and does not straddle a knot.
        do {
            float alpha = rand() / float(RAND_MAX);
            u = (1 - alpha) * (s.domainStart() + fd_eps) + alpha * (s.domainEnd() - fd_eps);
        } while (s.findSegmentContaining(u + fd_eps) != s.findSegmentContaining(u - fd_eps));

        SplineEvalMethod m = SplineEvalMethod::BASIS;

        // Analytical derivative
        Spline s_prime = derivative(s);
        Spline::Vec an = s_prime.evalPt(m, u);

        // Centered finite difference approximation.
        Spline::Vec fd = (s.evalPt(m, u + fd_eps) - s.evalPt(m, u - fd_eps)) / (2 * fd_eps);
        float relError = (fd - an).norm() / fd.norm();
        std::cout << relError << " (" << fd.transpose() << " vs " << an.transpose() << ")" << std::endl;
        maxError = std::max(maxError, relError);
    }
    std::cout << "Maximum relative error: " << maxError << std::endl;

    return 0;
}
