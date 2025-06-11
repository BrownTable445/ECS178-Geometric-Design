////////////////////////////////////////////////////////////////////////////////
// KnotViewer.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Viewer plugin for visualizing and manipulating the knot sequence.
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
//  Company:  University of California, Davis
//  Created:  11/16/2022 20:39:03
*///////////////////////////////////////////////////////////////////////////////
#ifndef KNOTVIEWER_HH
#define KNOTVIEWER_HH
#include <Eigen/Dense>
#include <igl/ortho.h>
#include <igl/unproject.h>
#include <algorithm>
#include "LinePlotter.hh"
#include "spline.hh"
#include "EditorBase.hh"

struct KnotViewer : public LinePlotter {
    using Base = LinePlotter;
    using VXf = Eigen::VectorXf;
    using SB = SplineBase<3>;

    const float tickHeight = 0.075f;
    const float repeatedKnotSpacing = 5e-3;

    virtual void init(igl::opengl::glfw::Viewer *_viewer) override {
        LinePlotter::init(_viewer);

        update();
    }

    size_t       numKnots() const { return m_spline->numKnots(); }
    float    knot(size_t i) const { return m_spline->getKnots()[i]; }
    float   &knot(size_t i)       { return m_spline->getKnots()[i]; }
    const float *knotData() const { return m_spline->getKnots().data(); }
          float *knotData()       { return m_spline->getKnots().data(); }

    void attachEditor(EditorBase *editor) { m_splineEditor = editor; }
    bool disabled() { return m_spline == nullptr; }
    void disable() { m_spline.reset(); clear(); }
    void enable(const std::shared_ptr<SB> &s) {
        m_spline = s;
        update();
    }

    void update();

    virtual Eigen::Matrix4f modelViewProjection() const override {
        Eigen::Matrix4f P;
        // Position the plot window in the lower right.
        igl::ortho(-0.95, 1.05, -0.25, 3.75, -1, 1, P);
        return P;
    }

    Eigen::Vector3f unprojectPt(int screen_x, int screen_y) const {
        const auto &c = viewer->core();
        Eigen::Vector3f win(screen_x, c.viewport(3) - screen_y, 0.0f);
        return igl::unproject(win, Eigen::Matrix4f::Identity().eval(), modelViewProjection(), c.viewport);
    }

    float knotCoordinateForPt(const Eigen::Vector3f &pt) const {
        return (1 - pt[0]) * knot(0) + pt[0] * knot(numKnots() - 1);
    }

    // The index of the first knot value greater than or equal to `u`.
    // (This is the "lower bound" in the sense of std::lower_bound.)
    size_t knotLowerBoundIndex(float u) const {
        return std::distance(knotData(), std::lower_bound(knotData(), knotData() + numKnots(), u));
    }

    virtual bool mouse_down(int button, int modifier) override;

    bool isDragging() const { return drag_idx < numKnots(); }

    virtual bool mouse_up(int button, int modifier) override {
        if (disabled() || !isDragging()) return false;
        drag_idx = NONE;
        return true;
    }

    virtual bool mouse_move(int mouse_x, int mouse_y) override;

    virtual void post_resize(int w, int h) override { Base::post_resize(w, h); update(); }

private:
    static constexpr size_t NONE = std::numeric_limits<size_t>::max();
    size_t drag_idx = NONE;
    std::shared_ptr<SB> m_spline;
    using LinePlotter::addEdges;
    using LinePlotter::addPolyline;
    EditorBase * m_splineEditor;
};

#endif /* end of include guard: KNOTVIEWER_HH */
