////////////////////////////////////////////////////////////////////////////////
// SplineinterpEditor.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Spline interpolation user interface for ECS 178.
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
//  Company:  University of California, Davis
//  Created:  10/28/2022 10:34:11
*///////////////////////////////////////////////////////////////////////////////
#ifndef SPLINEINTERPEDITOR_HH
#define SPLINEINTERPEDITOR_HH

#include "EditorBase.hh"
#include "BezierEditor.hh"
#include "PolyinterpEditor.hh"
#include "SplineEditor.hh"
#include "spline.hh"
#include "spline_interpolation.hh"

struct SplineInterpEditor : public EditorBase {
    SplineInterpEditor(AssignmentGUI *gui, DataSource &ds) : EditorBase(gui), dataSource(ds) {
        interpolant = std::make_shared<Spline>();
        interpPts = ds.generate();
        recomputeInterpolant();
    }
    DataSource &dataSource;

    Spline::ParametrizationType paramType = Spline::ParametrizationType::UNIFORM;

    std::shared_ptr<Spline> interpolant;
    std::vector<Eigen::Vector3f> interpPts;

    struct Selection {
        int controlPt = 0;
        bool dragging = false;
        bool dragAll = false;
    };

    Selection selection;
    int resolution = 100;
    bool showCage = true;

    virtual void updateView(IGLViewer &v) override {
        auto &plotter        = *dynamic_cast<LinePlotter         *>(v.plugins[0]);
        auto &vectorRenderer = *dynamic_cast<VectorFieldRenderer *>(v.plugins[1]);
        vectorRenderer.clear();
        plotter.clear();
        auto &vdata = v.data();
        vdata.clear();
        float dpiScale = static_cast<Viewer &>(v).getDPIScale();

        Eigen::Vector4f cageColor = Eigen::Vector4f(0.0f, 0.0f, 0.0f, 1.0f);
        Eigen::Vector4f curveColor = Eigen::Vector4f(0.3, 0.6f, 1.0f, 1.0f);

        vdata.add_points(pointsToMatrixRows(interpPts).cast<double>(), /* C = */ BezierEditor::cageColor(true).head<3>().transpose().cast<double>());
        SplineEditor::drawSpline(v, resolution, Spline::EvalMethod::DE_BOOR, *interpolant, /* selected */ false, /* showCage */ showCage);

        SplineEditor::knotViewer(v).enable(interpolant);

        vdata.point_size = 10.0f * dpiScale;
    }

    virtual bool callback_key_pressed(IGLViewer &viewer, unsigned int key, int modifiers) override {
        key = std::toupper(key);

        bool handled = false;
        switch (key) {
            case '=':
            case '+':
            case '-':
            case '_':
                resolution = std::min(std::max(resolution + ((key == '=' || key == '+') ? 1 : -1), 2), 1000);
                handled = true;
        }
        if (handled) updateView(viewer);
        return handled;
    }

    virtual bool callback_mouse_down(IGLViewer &viewer, int button, int modifier) override {
        if (button == int(IGLViewer::MouseButton::Right))
            return false;

        MouseRay(viewer, viewer.down_mouse_x, viewer.down_mouse_y)
           .select(1, [&](size_t ci) -> const std::vector<Eigen::Vector3f> & { return interpPts; },
                  [&](size_t /* ci */, size_t pi) {
                      selection.controlPt = pi;
                      selection.dragging  = true;
                      selection.dragAll = modifier == GLFW_MOD_SHIFT;
                  });

        return selection.dragging;
    }

    void recomputeInterpolant() {
        *interpolant = naturalCubicSplineInterpolant(interpPts, paramType);
    }

    virtual bool callback_mouse_move(IGLViewer &viewer, int mouse_x, int mouse_y) override {
        if (selection.dragging) {
            Eigen::Vector3f &cp = interpPts[selection.controlPt];
            Eigen::Vector3f newPt = dragPoint(viewer, cp, mouse_x, mouse_y, editMode == EditingMode::TWO_D);

            if (selection.dragAll) {
                Eigen::Vector3f diff = newPt - cp;
                for (auto &p : interpPts)
                    p += diff;
            }
            else cp = newPt;

            recomputeInterpolant();

            updateView(viewer);
            return true;
        }
        return false;
    }

    virtual bool callback_mouse_up(IGLViewer &viewer, int button, int modifier) override {
        if (selection.dragging) {
            selection.dragging = false;
            return true;
        }

        return false;
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
                for (auto &p : interpPts)
                    p[2] = 0.0;
                recomputeInterpolant();
                
                changed = true;
            }
            handled = true;
        }

        changed |= ImGui::Checkbox("Show Control Polyline", &showCage);

        if (dataSource.draw_menu(viewer)) {
            interpPts = dataSource.generate();
            recomputeInterpolant();
            changed = true;
        }

        int ptypeInt = int(paramType);
        if (ImGui::Combo("Parametrization Type", &ptypeInt, "Uniform\0Chord Length\0Centripetal\0")) {
            paramType = Spline::ParametrizationType(ptypeInt);
            recomputeInterpolant();
            changed = true;
        }

        if (ImGui::Button("Convert to Bézier Curves")) {
            auto result = splineToBezierSegments(*interpolant);
            BezierEditor *bezierEditor = dynamic_cast<BezierEditor *>(gui->getEditor("bezier"));
            if (bezierEditor == nullptr) throw std::logic_error("Failed to access Bézier editor");
            bezierEditor->bezierCurves = result;
            bezierEditor->editMode = editMode;
            gui->switchEditor("bezier", /* ignorePrevCameraSettings = */ true);
        }

        if (changed) updateView(viewer);
        return changed || handled;
    }

protected:
    virtual void m_enterImpl(IGLViewer &v) override {
        SplineEditor::knotViewer(v).attachEditor(this);
        SplineEditor::knotViewer(v).enable(interpolant);
    }

    virtual void m_leaveImpl(IGLViewer &v) override {
        SplineEditor::knotViewer(v).disable();
    }
};

#endif /* end of include guard: SPLINEINTERPEDITOR_HH */
