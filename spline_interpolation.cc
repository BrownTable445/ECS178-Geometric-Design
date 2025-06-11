#include "spline_interpolation.hh"
#include "TriDiagonalSystem.hh"
#include <algorithm>

std::vector<Eigen::Vector3f> solveInterplationSystem(const std::vector<float>& alpha,
    const std::vector<float>& beta,
    const std::vector<float>& gamma,
    const Eigen::Vector3f& leftEndConditionConstraintRow,
    const Eigen::Vector3f& rightEndConditionConstraintRow,
    const std::vector<Eigen::Vector3f>& w_x,
    bool verbose) {
    const size_t m = alpha.size() + 2;
    if ((m != beta.size() + 2) || (m != gamma.size() + 2) || (m != w_x.size() + 2))
        throw std::runtime_error("All diagonals should be the same length");

    std::vector<float> a(m - 1), d(m), c(m - 1);
    std::vector<Eigen::Vector3f> rhs(m);

    // The linear system is unfortunately not tri-diagonal as originally
    // formulated because the natural spline conditions in the first and last
    // rows contain three nonzero entries (rather than 2). However, because we
    // enforce a parametrization with zero-length padding, the second and
    // second-to-last rows have a single nonzero in their first/last columns.
    // Therefore, we obtain a tri-diagonal system by exchanging the first and
    // last equation pairs.
    //
    // [ d_0 c_0                                        ]
    // [ a_0 d_1 c_1                                    ]
    // [     a_1 d_2 c_2                                ]
    // [             . . .                              ]
    // [               . . .                            ]
    // [                   a_{m - 3} d_{m - 2} c_{m - 2}]
    // [                             a_{m - 2} d_{m - 1}]

    // alpha[1] up to alpha_[-1]
    for (size_t i = 1; i < m - 3; ++i) {
        a[i] = alpha[i];
        d[i + 1] = beta[i];
        c[i + 1] = gamma[i];
        rhs[i + 1] = w_x[i];
    }

    // After swapping, the first equation "[d[0], c[0] ... ] = rhs[0]" comes from alpha[0]
    d[0] = alpha[0];
    c[0] = beta[0];
    rhs[0] = w_x[0];
    if (gamma[0] != 0) throw std::runtime_error("gamma_0 should be exactly zero for zero-length padding parametrizations");

    // The second equation comes from the left end condition
    a[0] = leftEndConditionConstraintRow[0];
    d[1] = leftEndConditionConstraintRow[1];
    c[1] = leftEndConditionConstraintRow[2];
    rhs[1].setZero();

    // The second-to-last equation (m - 2) comes from the right end condition
    a[m - 3] = rightEndConditionConstraintRow[0];
    d[m - 2] = rightEndConditionConstraintRow[1];
    c[m - 2] = rightEndConditionConstraintRow[2];
    rhs[m - 2].setZero();

    // The last equation looks like [     a.back(), d.back() ] = [     0 gamma.back()]
    a[m - 2] = beta.back();
    d[m - 1] = gamma.back();
    if (alpha.back() != 0) throw std::runtime_error("gamma_0 should be exactly zero for zero-length padding parametrizations");
    rhs.back() = w_x.back();

    if (verbose) {
        std::cout << "diagonal: " << Eigen::Map<Eigen::RowVectorXf>(d.data(), d.size()) << std::endl;
        std::cout << "subdiagonal: " << Eigen::Map<Eigen::RowVectorXf>(a.data(), a.size()) << std::endl;
        std::cout << "supdiagonal: " << Eigen::Map<Eigen::RowVectorXf>(c.data(), c.size()) << std::endl;

        std::cout << "rhs:" << std::endl;
        for (size_t i = 0; i < rhs.size(); ++i)
            std::cout << "\t" << rhs[i].transpose() << std::endl;
    }

    TriDiagonalSystem<float> sys(a, d, c);
    return sys.solve(rhs);
}

Spline naturalCubicSplineInterpolant(std::vector<Eigen::Vector3f>& dataPts, Spline::ParametrizationType paramType) {
    bool verbose = true;
    if (verbose) std::cout << "Generating interpolant for " << dataPts.size() << " data points" << std::endl;
    Spline result;

    const size_t n = 3;
    const size_t P = dataPts.size() - 1;
    const size_t L = P + 2;

    // TODO 3.1: infer the knot sequence from the data point spacings
    // and the requested parameterization type. You should always
    // use lengths of zero for the first two and last two parameter intervals
    // (i.e., pretend that `repeatedEndKnots = true` was passed).
    auto& knots = result.knots;
    //knots.setLinSpaced(P + 2 * (n - 1) + 1, 0, 1);
    knots.conservativeResize(P + 2 * (n - 1) + 1);
    if (paramType == Spline::ParametrizationType::UNIFORM) {
        for (size_t i = 2; i <= L; ++i) {
            float temp = 1.0f / (L - 2 + 1) * (i + 1 - 2);
            knots[i] = 1.0f / (L - 2 + 1) * (i + 1 - 2);
        }
        knots[0] = knots[2];
        knots[1] = knots[2];
        knots[L + 1] = knots[L];
        knots[L + 2] = knots[L];
        //std::cout << knots << std::endl;
    }
    else {
        knots[0] = 0;
        knots[1] = 0;
        for (size_t i = 2; i < L; ++i) {
            knots[i] = (dataPts[i - 2] + dataPts[i + 1 - 2]).norm() + knots[i - 1];
            if (paramType == Spline::ParametrizationType::CENTRIPETAL)
                knots[i] = sqrt(knots[i]);
        }
        knots[0] = knots[2];
        knots[1] = knots[2];
        knots[L + 1] = knots[L];
        knots[L + 2] = knots[L];
    }

    // TODO 3.2: solve for the de Boor points such that `dataPts` are
    // interpolated. The tri-diagonal system set up and solved here
    // essentially uses the identity matrix (up to a permutation
    // of the first two and last two rows).
    // Update the variables `alpha`, `beta`, `gamma`, w_x`,
    // `leftEndConditionConstraintRow`, and `rightEndConditionConstraintRow`
    // appropriately.
    std::vector<float> alpha(P + 1, 0.0f),
        beta(P + 1, 1.0f),
        gamma(P + 1, 0.0f);
    std::vector<Eigen::Vector3f> w_x(dataPts.size(), Eigen::Vector3f::Zero());
    Eigen::Vector3f leftEndConditionConstraintRow(1, 0, 0),
        rightEndConditionConstraintRow(0, 0, 1);

    for (size_t i = 0; i <= P; ++i) {
        float delta_i = /*i == 0 || i == 1 || i == L || i == L + 1 ? 0 :*/ knots[i] - knots[i + 1];
        float delta_i_plus_1 = /*i + 1 == 0 || i + 1 == 1 || i + 1 == L || i + 1 == L + 1 ? 0 :*/ knots[i + 1] - knots[i + 2];
        float delta_i_plus_2 = /*i + 2 == 0 || i + 2 == 1 || i + 2 == L || i + 2 == L + 1 ? 0 :*/ knots[i + 2] - knots[i + 3];
        float delta_i_plus_3 = /*i + 3 == 0 || i + 3 == 1 || i + 3 == L || i + 3 == L + 1 ? 0 :*/ knots[i + 3] - knots[i + 4];
        alpha[i] = delta_i_plus_2 * delta_i_plus_2 / (delta_i + delta_i_plus_1 + delta_i_plus_2);
        beta[i] = delta_i_plus_2 * (delta_i + delta_i_plus_1) / (delta_i + delta_i_plus_1 + delta_i_plus_2) + delta_i_plus_1 * (delta_i_plus_2 + delta_i_plus_3) / (delta_i_plus_1 + delta_i_plus_2 + delta_i_plus_3);
        gamma[i] = delta_i_plus_1 * delta_i_plus_1 / (delta_i_plus_1 + delta_i_plus_2 + delta_i_plus_3);
        w_x[i] = (delta_i_plus_1 + delta_i_plus_2) * dataPts[i];
    }
    leftEndConditionConstraintRow[1] = (knots[3] - knots[4]) / (knots[2] - knots[3] + knots[3] - knots[4]) - 2;
    leftEndConditionConstraintRow[2] = (knots[2] - knots[3]) / (knots[2] - knots[3] + knots[3] - knots[4]);
    rightEndConditionConstraintRow[0] = (knots[L - 1] - knots[L]) / (knots[L - 2] - knots[L - 1] + knots[L - 1] - knots[L]);
    rightEndConditionConstraintRow[1] = (knots[L - 2] - knots[L - 1]) / (knots[L - 2] - knots[L - 1] + knots[L - 1] - knots[L]) - 2;

    if (verbose) {
        std::cout << "leftEndConditionConstraintRow: " << leftEndConditionConstraintRow.transpose() << std::endl;
        std::cout << "rightEndConditionConstraintRow: " << rightEndConditionConstraintRow.transpose() << std::endl;
    }

    result.controlPts = pointsToMatrixRows(solveInterplationSystem(
        alpha, beta, gamma,
        leftEndConditionConstraintRow,
        rightEndConditionConstraintRow, w_x, verbose));

    if (verbose) {
        std::cout << "Computed ";
        result.info();
    }

    return result;
}
