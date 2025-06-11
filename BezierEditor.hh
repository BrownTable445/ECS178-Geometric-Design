////////////////////////////////////////////////////////////////////////////////
// BezierEditor.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Bézier curve editor user interface for ECS 178.
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
//  Company:  University of California, Davis
//  Created:  10/28/2022 10:34:11
*///////////////////////////////////////////////////////////////////////////////
#ifndef BEZIEREDITOR_HH
#define BEZIEREDITOR_HH

#include <igl/project.h>
#include <igl/unproject.h>
#include "EditorBase.hh"

struct BezierEditor : public EditorBase {
    BezierEditor(AssignmentGUI *gui) :
        EditorBase(gui), bezierCurves(1) { }

    std::vector<BezierCurve> bezierCurves;

    enum class EvalProbeVisualization {
        NONE = 0, DE_CASTELJAU = 1, FRAME = 2
    };

    int   resolution      = 100;
    float eval_t          = 0.5;
    BezierCurve::EvalMethod method = BezierCurve::EvalMethod::DE_CASTELJAU;

    struct Selection {
        int curve = 0;
        int controlPt = 0;
        float clicked_t = 0;
        bool dragging = false;
        bool dragAll = false;
    };

    Selection selection;
    EvalProbeVisualization visMode = EvalProbeVisualization::NONE;

    float intersectionTolerance = 1e-4;
    std::vector<Eigen::Vector3f> intersectionPts;

    static Eigen::Vector4f cageColor(bool selected) { return selected ? Eigen::Vector4f(0.0f, 0.0f, 0.0f, 1.0f)
                                                                      : Eigen::Vector4f(0.6f, 0.6f, 0.6f, 1.0f); }

    static Eigen::Vector4f curveColor(bool selected) { return selected ? Eigen::Vector4f(0.3, 0.6f, 1.0f, 1.0f)
                                                                       : Eigen::Vector4f(0.7f, 0.8, 1.0f, 1.0f); }
    static void drawBezierCurve(IGLViewer &v, int resolution, BezierCurve::EvalMethod method, BezierCurve &c, bool selected) {
        auto &plotter        = *dynamic_cast<LinePlotter         *>(v.plugins[0]);
        float dpiScale = static_cast<Viewer &>(v).getDPIScale();

        // Draw the control polygon
        plotter.addPolyline(c.getControlPts(), /* closed = */ false, /* thickness = */ 2.0f * dpiScale, /* color = */ cageColor(selected));
        v.data().add_points(pointsToMatrixRows(c.getControlPts()).cast<double>(), /* C = */ cageColor(selected).head<3>().transpose().cast<double>());

        // Evaluate and draw the Bezier curve
        static std::vector<Eigen::Vector3f> evaluatedCurve; // re-used across evaluations/updates
        c.eval(resolution, method, evaluatedCurve);
        plotter.addPolyline(evaluatedCurve, /* closed = */ false, /* thickness = */ 5.0f * dpiScale, /* color = */ curveColor(selected));
    }

    void updateView(IGLViewer &v) override {
        auto &plotter        = *dynamic_cast<LinePlotter         *>(v.plugins[0]);
        auto &vectorRenderer = *dynamic_cast<VectorFieldRenderer *>(v.plugins[1]);
        vectorRenderer.clear();
        plotter.clear();
        auto &vdata = v.data();
        vdata.clear();
        float dpiScale = static_cast<Viewer &>(v).getDPIScale();

        for (size_t i = 0; i < bezierCurves.size(); ++i) {
            const auto &c = bezierCurves[i];
            drawBezierCurve(v, resolution, method, bezierCurves[i], i == size_t(selection.curve));
        }

        if (size_t(selection.curve) < bezierCurves.size()) {
            // Visualize quantities for the selected curve.
            const auto &c = bezierCurves[selection.curve];
            // Evaluate at the probe point
            Eigen::Vector3f p = c.evalPt(method, eval_t);
            vdata.add_points(p.transpose().cast<double>(), /* C = */ curveColor(true).head<3>().transpose().cast<double>());

            if (visMode == EvalProbeVisualization::FRAME) {
                Eigen::Vector3f tangent, normal;
                p = c.probeCurve(eval_t, tangent, normal);
                VectorFieldRenderer::PtArray frenetPts(3, 3), frenetVectors(3, 3);
                VectorFieldRenderer::ColorArray colors(3, 4);
                frenetPts << p.transpose(),
                             p.transpose(),
                             p.transpose();
                frenetVectors << tangent.transpose(),
                                  normal.transpose(),
                                  tangent.cross(normal).normalized().transpose();
                colors << 1, 0, 0, 1,
                          0, 1, 0, 1,
                          0, 0, 1, 1;
                vectorRenderer.setField(frenetPts, frenetVectors, colors);
            }
            if (visMode == EvalProbeVisualization::DE_CASTELJAU) {
                Eigen::MatrixX3f V;
                Eigen::MatrixX2i E;
                bezierCurves[selection.curve].visualizeDeCasteljau(eval_t, V, E);
                plotter.addEdges(V, E, 1.0f * dpiScale, Eigen::Vector4f(0.5f, 0.5f, 0.5f, 1.0f));
                vdata.add_points(V.cast<double>(), /* C = */ Eigen::RowVector3d(0.5f, 0.5f, 0.5f));
            }
        }

        vdata.add_points(pointsToMatrixRows(intersectionPts).cast<double>(), Eigen::RowVector3d(0.3, 1.0f, 0.3f));

        vdata.point_size = 10.0f * dpiScale;
    }

    bool callback_key_pressed(IGLViewer &viewer, unsigned int key, int modifiers) override {
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

    bool callback_mouse_down(IGLViewer &viewer, int button, int modifier) override {
        if (button == int(IGLViewer::MouseButton::Right))
            return false;

        int oldSel = selection.curve;
        constexpr size_t res = 100; // resolution used for Bezier curve click detection
        MouseRay(viewer, viewer.down_mouse_x, viewer.down_mouse_y).select(
            // First detect clicks on the control points.
            bezierCurves.size(),
            [&](size_t bi) { return bezierCurves[bi].getControlPts(); },
            [&](size_t ci, size_t pi) {
                selection.curve     = ci;
                selection.controlPt = pi;
                selection.dragging  = true;
                selection.dragAll = modifier == GLFW_MOD_SHIFT;
            },
            // Next, detect clicks on the curve itself.
            bezierCurves.size(),
            [&](size_t bi) {
                static std::vector<Eigen::Vector3f> evaluatedCurve; // re-used
                bezierCurves[bi].eval(res, BezierCurve::EvalMethod::HORNER, evaluatedCurve);
                return evaluatedCurve;
            },
            [&](size_t ci, size_t pi) {
                selection.curve     = ci;
                selection.controlPt = -1;
                selection.dragging  = true;
                selection.clicked_t = float(pi) / (res - 1);
                selection.dragAll = modifier == GLFW_MOD_SHIFT;
            }
        );

        if (selection.curve != oldSel) updateView(viewer);

        return selection.dragging;
    }

    bool callback_mouse_move(IGLViewer &viewer, int mouse_x, int mouse_y) override {
        if (selection.dragging) {
            auto &c = bezierCurves[selection.curve];
            if (selection.controlPt >= 0) {
                // Dragging on control point.
                Eigen::Vector3f &cp = c.getControlPts()[selection.controlPt];

                Eigen::Vector3f newPt = dragPoint(viewer, cp, mouse_x, mouse_y, editMode == EditingMode::TWO_D);
                if (selection.dragAll) {
                    Eigen::Vector3f diff = newPt - cp;
                    for (auto &p : c.getControlPts())
                        p += diff;
                }
                else cp = newPt;
            }
            else {
                // Dragging a point on the curve directly.
                Eigen::Vector3f origPt = c.evalPt(BezierCurve::EvalMethod::HORNER, selection.clicked_t);
                Eigen::Vector3f  newPt = dragPoint(viewer, origPt, mouse_x, mouse_y, editMode == EditingMode::TWO_D);

                if (selection.dragAll) {
                    Eigen::Vector3f diff = newPt - origPt;
                    for (auto &p : c.getControlPts())
                        p += diff;
                }
                else c.drag(selection.clicked_t, newPt);
            }

            updateView(viewer);
            return true;
        }
        return false;
    }

    bool callback_mouse_up(IGLViewer &viewer, int button, int modifier) override {
        if (selection.dragging) {
            selection.dragging = false;
            return true;
        }

        return false;
    };

    bool draw_menu(IGLViewer &viewer) override {
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
                for (auto &c : bezierCurves) {
                    for (auto &p : c.getControlPts())
                        p[2] = 0.0;
                }
                changed = true;
            }
            handled = true;
        }

        if (ImGui::Button("Add Curve", ImVec2(ImGui::GetContentRegionAvail().x*0.5f, 0.0f))) {
            Eigen::Vector3f t = Eigen::Vector3f::Zero();
            t.topRows(2).setRandom();
            bezierCurves.emplace_back(t);
            selection.curve = bezierCurves.size() - 1;
            changed = true;
        }
        ImGui::SameLine();

        auto conditionallyDisabledButton = [](bool disabled, const char *label, const ImVec2 &size = ImVec2(0,0)) {
            if (disabled)
                ImGui::BeginDisabled();
            bool result = ImGui::Button(label, size);
            if (disabled) ImGui::EndDisabled();
            return result;
        };

        if (conditionallyDisabledButton(bezierCurves.size() <= 1, "Delete Curve", ImVec2(ImGui::GetContentRegionAvail().x, 0.0f))) {
            bezierCurves.erase(bezierCurves.begin() + selection.curve);
            selection.curve = std::min(size_t(selection.curve), bezierCurves.size() - 1);
            changed = true;
        }

        int degree = bezierCurves[selection.curve].degree();
        if (ImGui::InputInt("Curve Degree", &degree, 1)) {
            degree = std::min(std::max(degree, 1), 12);
            bezierCurves[selection.curve].setDegree(degree);
            changed = true;
        }

        int evalMethodInt = int(method);
        if (ImGui::Combo("Eval Method", &evalMethodInt, "De Casteljau\0Bernstein\0Horner\0")) {
            method = BezierCurve::EvalMethod(evalMethodInt);
            changed = true;
        }

        if (ImGui::SliderFloat("Eval t", &eval_t, /* vmin = */ 0.0f, /* vmax = */ 1.0f))
            changed = true;

        int visModeInt = int(visMode);
        if (ImGui::Combo("Eval Pt Visualization", &visModeInt, "Point Only\0De Casteljau\0Frenet Frame\0")) {
            visMode = EvalProbeVisualization(visModeInt);
            changed = true;
        }

        if (ImGui::Button("Subdivide")) {
            auto c = bezierCurves[selection.curve];
            bezierCurves.emplace_back();
            c.subdivide(eval_t, bezierCurves[selection.curve], bezierCurves.back());

            changed = true;
        }

        // ImGUI's logarithmic sliders somehow *still* seem to be broken even
        // when passing the flags `ImGuiSliderFlags_Logarithmic | ImGuiSliderFlags_NoRoundToFormat` :(
        // Use the workaround suggested in "https://github.com/ocornut/imgui/issues/642#issuecomment-541240419"
        auto logarithmicSlider = [](const char *label, float &val, float vmin, float vmax) {
            char val_str[64];
            snprintf(val_str, 64, "%.3e", val);
            val = log(val);
            bool ret = ImGui::SliderFloat(label, &val, log(vmin), log(vmax), val_str);
            val = exp(val);
            return ret;
        };
        handled |= logarithmicSlider("Intersection Tolerance", intersectionTolerance, 1e-6, 1e-1);

        if (conditionallyDisabledButton(bezierCurves.size() <= 1, "Calculate Intersections")) {
            intersectionPts.clear();
            // Perform intersection tests against all pairs of curves.
            for (size_t i = 0; i < bezierCurves.size(); ++i) {
                for (size_t j = i + 1; j < bezierCurves.size(); ++j) {
                    getIntersections(bezierCurves[i], bezierCurves[j], intersectionPts, intersectionTolerance);
                }
            }

            std::cout << "Found " << intersectionPts.size() << " approximate intersection point(s)";

            mergeDuplicateIntersections(intersectionPts, intersectionTolerance);

            std::cout << " which merged to " << intersectionPts.size() << " intersection(s)" << std::endl;
            changed = true;
        }

        if (changed) updateView(viewer);
        return changed || handled;

    }
};

#endif /* end of include guard: BEZIEREDITOR_HH */
