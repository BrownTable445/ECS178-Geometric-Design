////////////////////////////////////////////////////////////////////////////////
// PolyapproxEditor.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Polynomial approximation user interface for ECS 178.
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
//  Company:  University of California, Davis
//  Created:  10/28/2022 18:14:57
*///////////////////////////////////////////////////////////////////////////////
#ifndef POLYAPPROXEDITOR_HH
#define POLYAPPROXEDITOR_HH

#include "EditorBase.hh"
#include "BezierEditor.hh"
#include "PolyinterpEditor.hh"

struct PolyapproxEditor : public PolyinterpEditor {
    PolyapproxEditor(AssignmentGUI *gui, DataSource &ds, PolyInterp &pi)
        : PolyinterpEditor(gui, ds, pi) {
        bezierCurve.setDegree(4);
    }
    BezierCurve bezierCurve;
    float smoothingWeight = 0;

    virtual void updateView(IGLViewer &v) override {
        auto &plotter        = *dynamic_cast<LinePlotter         *>(v.plugins[0]);
        auto &vectorRenderer = *dynamic_cast<VectorFieldRenderer *>(v.plugins[1]);

        vectorRenderer.clear();
        plotter.clear();
        auto &vdata = v.data();
        vdata.clear();
        float dpiScale = static_cast<Viewer &>(v).getDPIScale();

        bezierCurve.approximate(polynomialInterpolant.t_values, polynomialInterpolant.getInterpPts(), smoothingWeight);

        BezierEditor::drawBezierCurve(v, resolution, BezierCurve::EvalMethod::HORNER, bezierCurve, false);

        vdata.add_points(pointsToMatrixRows(polynomialInterpolant.getInterpPts()).cast<double>(), /* C = */ BezierEditor::cageColor(true).head<3>().transpose().cast<double>());
        vdata.point_size = 10.0f * dpiScale;
    }

    virtual bool draw_menu(IGLViewer &viewer) override {
        bool handled = false, changed = false;
        if (ImGui::InputInt("Curve Eval Resolution", &resolution, 1)) {
            resolution = std::min(std::max(resolution, 2), 1000);
            changed = true;
        }
        int editModeInt = int(editMode);
        if (ImGui::Combo("Dimension", &editModeInt, "2D\0003D\0")) {
            editMode = EditingMode(editModeInt);
            applyEditingMode(viewer);
            if (editMode == EditingMode::TWO_D) {
                // When switching to 2D editing mode, we need to project the design onto the Y = 0 axis.
                auto allPts = polynomialInterpolant.getInterpPts(); // we must edit and then set a copy of all points for proper updating...
                for (auto &p : allPts)
                    p[2] = 0.0;
                polynomialInterpolant.setInterpPts(allPts);
                changed = true;
            }
            handled = true;
        }

        if (dataSource.draw_menu(viewer)) {
            polynomialInterpolant.setInterpPts(dataSource.generate());
            changed = true;
        }

        int ptypeInt = int(polynomialInterpolant.getParamType());
        if (ImGui::Combo("Parametrization Type", &ptypeInt, "Uniform\0Chord Length\0XCOORD\0")) {
            polynomialInterpolant.setParamType(PolyInterp::ParametrizationType(ptypeInt));
            changed = true;
        }

        int degree = bezierCurve.degree();
        if (ImGui::InputInt("Curve Degree", &degree, 1)) {
            degree = std::min(std::max(degree, 1), 25);
            bezierCurve.setDegree(degree);
            changed = true;
        }

        changed |= logarithmicSliderWithZero("Smoothing Regularization", smoothingWeight, 1e-5, 1e4);

        if (changed) updateView(viewer);
        return changed || handled;
    }
};

#endif /* end of include guard: POLYAPPROXEDITOR_HH */
