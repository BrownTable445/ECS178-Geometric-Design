////////////////////////////////////////////////////////////////////////////////
// tensor_product_surface.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Tensor product NURBS surface.
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
//  Company:  University of California, Davis
//  Created:  12/03/2022 17:28:49
*///////////////////////////////////////////////////////////////////////////////
#ifndef TENSOR_PRODUCT_SURFACE_HH
#define TENSOR_PRODUCT_SURFACE_HH

#include "spline.hh"

template<size_t Dimension>
struct TensorProductSurface_T {
    using Scalar  = float;
    using Vec     = Eigen::Matrix<Scalar, Dimension, 1>;
    using M2d     = Eigen::Matrix<Scalar, 2, 2>;
    using V2d     = Eigen::Matrix<Scalar, 2, 1>;
    using Mat     = Eigen::Matrix<Scalar, Dimension, Dimension>;
    using VXf     = Eigen::VectorXf;
    using MXDf    = Eigen::Matrix<Scalar, Eigen::Dynamic, Dimension>;
    using MX2f    = Eigen::Matrix<Scalar, Eigen::Dynamic, 2>;
    using HSpline = Spline_T<Dimension + 1>;
    using HVec    = typename HSpline::Vec;

    TensorProductSurface_T(const Vec &t = Vec::Zero()) {
        const size_t gridSize = 5;
        controlPts.resize(gridSize, gridSize);
        weights.setOnes(gridSize, gridSize);

        for (size_t i = 0; i < gridSize; ++i) {
            float x = 2 * (i / (gridSize - 1.0f)) - 1;
            for (size_t j = 0; j < gridSize; ++j) {
                float z = 2 * (j / (gridSize - 1.0f)) - 1;
                controlPts(i, j) << x, 0.5 * x * z, -z;
            }
        }

        // Initialize the u- and v-spline curves with sensible values.
        u_spline.controlPts.setZero(gridSize, Dimension + 1);
        v_spline.controlPts.setZero(gridSize, Dimension + 1);
        u_spline.controlPts.rightCols(1).setOnes();
        v_spline.controlPts.rightCols(1).setOnes();

        u_spline.inferKnots();
        v_spline.inferKnots();
    }

    HVec evalPtHomogCoord(double u, double v) const;
    Vec evalPt(double u, double v) const;

    const Vec &controlPt(size_t i, size_t j) const { return controlPts(i, j); }
    void    setControlPt(size_t i, size_t j, const Vec &p) { controlPts(i, j) = p; }
    void eval(size_t resolution, MX2f &U, MXDf &V, Eigen::MatrixX3i &F) const;

    void insertUKnot(double ubar);
    void insertVKnot(double vbar);

    using  FlatPtType = Eigen::Map<      Eigen::Matrix<Scalar, Eigen::Dynamic, Dimension, Eigen::RowMajor>>;
    using CFlatPtType = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Dimension, Eigen::RowMajor>>;
    CFlatPtType getControlPtsFlattened() const { return CFlatPtType((const Scalar *) controlPts.data(), controlPts.size(), Dimension); }
     FlatPtType getControlPtsFlattened()       { return  FlatPtType(      (Scalar *) controlPts.data(), controlPts.size(), Dimension); }

    // Note: grid point columns are stored contiguously
    CFlatPtType getControlPtColFlattened(int j) const { return CFlatPtType((const Scalar *) controlPts.col(j).data(), controlPts.rows(), Dimension); }
     FlatPtType getControlPtColFlattened(int j)       { return  FlatPtType(      (Scalar *) controlPts.col(j).data(), controlPts.rows(), Dimension); }

    void unflattenIndex(size_t flatIdx, size_t &i, size_t &j) const {
        const size_t rs = controlPts.rows();
        i = flatIdx % rs;
        j = flatIdx / rs;
    }

    std::vector<Vec> getControlPtRow(int i) const {
        if (i >= controlPts.rows()) throw std::runtime_error("Index out of bounds");
        std::vector<Vec> result(controlPts.cols());
        for (size_t j = 0; j < controlPts.cols(); ++j)
            result[j] = controlPts(i, j);
        return result;
    }

    std::vector<Vec> getControlPtCol(int j) const {
        if (j >= controlPts.cols()) throw std::runtime_error("Index out of bounds");
        std::vector<Vec> result(controlPts.rows());
        for (size_t i = 0; i < controlPts.rows(); ++i)
            result[i] = controlPts(i, j);
        return result;
    }

    // Differential geometry quantities
    HVec partial_xtilde(double u, double v, size_t diff_u, size_t diff_v) const;
    Eigen::Matrix<Scalar, Dimension    , 2> jacobian(double u, double v) const;
    Vec normal(double u, double v) const;
    M2d firstFundamentalForm(double u, double v) const;
    M2d secondFundamentalForm(double u, double v) const;

    V2d principalCurvatures(double u, double v) const;
    V2d gaussianAndMeanCurvature(double u, double v) const;
    double    meanCurvature(double u, double v) const;

    using WeightGrid = Eigen::MatrixXf;
    using ControlPtGrid = Eigen::Matrix<Vec, Eigen::Dynamic, Eigen::Dynamic>;
    WeightGrid weights;
    ControlPtGrid controlPts;

    // Mutable to allow modification of the control points in the `const`
    // `evalPt` routine...
    // It would be better to implement the spline curve evaluation routines to
    // operate on external data sources.
    mutable HSpline u_spline, v_spline;
};

using TensorProductSurface = TensorProductSurface_T<3>;

template<size_t Dimension>
std::shared_ptr<TensorProductSurface_T<Dimension>> ruled_surface(const NURBS<Dimension> &spline_1, const NURBS<Dimension> &spline_2);

template<size_t Dimension>
std::shared_ptr<TensorProductSurface_T<Dimension>> surface_of_revolution(const NURBS<Dimension> &spline, size_t numSegments);

#endif /* end of include guard: TENSOR_PRODUCT_SURFACE_HH */
