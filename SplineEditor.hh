////////////////////////////////////////////////////////////////////////////////
// SplineEditor.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Spline interpolation and manipulation user interface for ECS 178.
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
//  Company:  University of California, Davis
//  Created:  11/13/2022 19:12:20
*///////////////////////////////////////////////////////////////////////////////
#ifndef SPLINEEDITOR_HH
#define SPLINEEDITOR_HH

#include "EditorBase.hh"
#include "BezierEditor.hh"
#include "SurfaceEditor.hh"
#include "PolyinterpEditor.hh"
#include "spline.hh"
#include "KnotViewer.hh"
#include "AssignmentGUI.hh"

struct SplineEditor : public EditorBase {
    using SB = SplineBase<3>;
    std::vector<std::shared_ptr<SB>> splines;

    SplineEditor(AssignmentGUI *gui) : EditorBase(gui), splines(1, std::make_shared<Spline>()) { }

    using Selection = BezierEditor::Selection;
    Selection selection;
    int   resolution      = 100;
    float eval_u          = 0.5;
    Spline::EvalMethod method = Spline::EvalMethod::DE_BOOR;
    bool showCage         = true;
    bool showDeBoor       = false;

    bool recomputeParametrizationOnDrag = true;
    int  circleSegments = 3, revSegments = 4, keptSegments = 3;

          SB &selectedSpline()       { return *splines[selection.curve]; }
    const SB &selectedSpline() const { return *splines[selection.curve]; }

    static KnotViewer &knotViewer(IGLViewer &v) {
        return *dynamic_cast<KnotViewer*>(v.plugins[2]);
    }

    std::shared_ptr<SB> selectedSplinePtr() {
        if (selection.curve > splines.size()) throw std::runtime_error("No valid spline selection");
        return splines[selection.curve];
    }

    virtual void updateView(IGLViewer &v) override {
        auto &plotter        = *dynamic_cast<LinePlotter         *>(v.plugins[0]);
        auto &vectorRenderer = *dynamic_cast<VectorFieldRenderer *>(v.plugins[1]);
        vectorRenderer.clear();
        plotter.clear();
        auto &vdata = v.data();
        vdata.clear();
        float dpiScale = static_cast<Viewer &>(v).getDPIScale();

        for (size_t i = 0; i < splines.size(); ++i)
            drawSpline(v, resolution, method, *splines[i], i == size_t(selection.curve), showCage);

        knotViewer(v).enable(selectedSplinePtr());

        // Evaluate at the probe point
        const auto &s = selectedSpline();
        eval_u = clamp(eval_u, s.domainStart(), s.domainEnd());
        Eigen::Vector3f p = s.evalPt(method, eval_u);
        vdata.add_points(p.transpose().cast<double>(), /* C = */ BezierEditor::curveColor(true).head<3>().transpose().cast<double>());

        auto ordinarySpline = std::dynamic_pointer_cast<Spline>(selectedSplinePtr());
        if (showDeBoor && ordinarySpline) {
            Eigen::MatrixX3f V;
            Eigen::MatrixX2i E;
            ordinarySpline->visualizeDeBoor(eval_u, V, E);
            plotter.addEdges(V, E, 1.0f * dpiScale, Eigen::Vector4f(0.5f, 0.5f, 0.5f, 1.0f));
            vdata.add_points(V.cast<double>(), /* C = */ Eigen::RowVector3d(0.5f, 0.5f, 0.5f));
        }

        vdata.point_size = 10.0f * dpiScale;
    }

    static void drawSpline(IGLViewer &v, int resolution, Spline::EvalMethod method, SB &c, bool selected, bool showCage) {
        auto &plotter  = *dynamic_cast<LinePlotter *>(v.plugins[0]);
        float dpiScale = static_cast<Viewer &>(v).getDPIScale();

        if (showCage) {
            // Draw the control polygon
            const auto &d = c.getControlPts();
            plotter.addPolyline(pointsFromMatrixRows(d), /* closed = */ false, /* thickness = */ 2.0f * dpiScale, /* color = */ BezierEditor::cageColor(selected));
            v.data().add_points(d.cast<double>(), /* C = */ BezierEditor::cageColor(selected).head<3>().transpose().cast<double>());
        }

        // Evaluate and draw the spline curve
        static std::vector<Eigen::Vector3f> evaluatedCurve; // re-used across evaluations/updates
        c.eval(resolution, method, evaluatedCurve);
        plotter.addPolyline(evaluatedCurve, /* closed = */ false, /* thickness = */ 5.0f * dpiScale, /* color = */ BezierEditor::curveColor(selected));
    }

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
                // When switching to 2D editing mode, we need to project the design onto the Z = 0 axis.
                for (auto &c : splines)
                    c->make2D();
                changed = true;
            }
            handled = true;
        }

        changed |= ImGui::Checkbox("Show Control Polyline", &showCage);

        ImGui::Separator();
        ImGui::Text("Evaluation");

        int evalMethodInt = int(method);
        if (ImGui::Combo("Eval Method", &evalMethodInt, "De Boor\0Basis Functions\0")) {
            method = Spline::EvalMethod(evalMethodInt);
            changed = true;
        }

        changed |= ImGui::Checkbox("de Boor Visualization", &showDeBoor);

        if (ImGui::SliderFloat("Eval u", &eval_u, /* vmin = */ selectedSpline().domainStart(), /* vmax = */ selectedSpline().domainEnd())) {
            changed = true;
            knotViewer(viewer).update();
        }

        SB &s = selectedSpline();
        if (ImGui::Button("Write Evaluated Basis Functions")) {
            std::ofstream outFile("basis_functions.txt");
            if (!outFile.is_open()) throw std::runtime_error("Failed to open output file");

            float start = s.domainStart() + 1e-6;
            float   end = s.domainEnd() - 1e-6;
            for (size_t i = 0; i < resolution; ++i) {
                double u = start + (end - start) * i / (resolution - 1);
                outFile << u;
                for (size_t l = 0; l < s.numControlPts(); ++l)
                    outFile << "\t" << s.bsplineBasisFunction(s.degree(), l, u);
                outFile << std::endl;
            }
        }

        ImGui::Separator();
        ImGui::Text("Operations");

        if (ImGui::Button("Add Curve", ImVec2(ImGui::GetContentRegionAvail().x*0.5f, 0.0f))) {
            Eigen::Vector3f t = Eigen::Vector3f::Zero();
            t.topRows(2).setRandom();
            splines.emplace_back(std::make_shared<Spline>(t));
            selection.curve = splines.size() - 1;
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

        if (conditionallyDisabledButton(splines.size() <= 1, "Delete Curve", ImVec2(ImGui::GetContentRegionAvail().x, 0.0f))) {
            splines.erase(splines.begin() + selection.curve);
            selection.curve = std::min(size_t(selection.curve), splines.size() - 1);
            changed = true;
        }

        if (ImGui::Button("Insert Knot")) {
            selectedSpline().insertKnot(eval_u);
            changed = true;
        }

        ImGui::SameLine();

        auto ordinarySpline = std::dynamic_pointer_cast<Spline>(selectedSplinePtr());
        if (conditionallyDisabledButton(ordinarySpline == nullptr, "Convert to Bézier Curves")) {
            auto result = splineToBezierSegments(*ordinarySpline);
            BezierEditor *bezierEditor = dynamic_cast<BezierEditor *>(gui->getEditor("bezier"));
            if (bezierEditor == nullptr) throw std::logic_error("Failed to access Bézier editor");
            bezierEditor->bezierCurves = result;
            bezierEditor->editMode = editMode;
            gui->switchEditor("bezier", /* ignorePrevCameraSettings = */ true);
        }

        // Spline num pts
        {
            int npts = s.numControlPts();
            ImGui::InputInt("Num Pts", &npts, 1);
            npts = std::max<size_t>(npts, 2);
            if (npts != s.numControlPts()) {
                const auto &oldPts = s.getControlPts();
                size_t N = oldPts.rows() - 1; // old number of edges
                Spline::MXDf newPts(npts, 3);

                for (size_t i = 0; i < npts - 1; ++i) {
                    // Sample the old control polyline
                    double t = i / double(npts - 1);

                    size_t j = std::floor(t * N);
                    double alpha = N * t - j;

                    newPts.row(i) = (1 - alpha) * oldPts.row(j)
                                  +      alpha  * oldPts.row(j + 1);
                }
                newPts.row(npts - 1) = oldPts.row(N);
                s.setControlPts(newPts);
                s.inferKnots();
                changed = true;
            }
        }

        // Spline degree
        {
            int degree = s.degree();
            ImGui::InputInt("Spline Degree", &degree, 1);
            // Maximum degree is numPts - 1
            degree = clamp<size_t>(degree, 1, s.numControlPts() - 1);
            if (degree != s.degree()) {
                s.degree() = degree;
                s.inferKnots();
                changed = true;
            }
        }

        ImGui::Separator();
        ImGui::Text("Surface Construction");

        if (conditionallyDisabledButton((splines.size() != 2) || (splines[0]->degree() != splines[1]->degree()), "Ruled Surface")) {
            auto result = ruled_surface(*NURBSForSpline(splines[0]), *NURBSForSpline(splines[1]));
            SurfaceEditor *surfEditor = dynamic_cast<SurfaceEditor *>(gui->getEditor("surface"));
            if (surfEditor == nullptr) throw std::logic_error("Failed to access surface editor");
            surfEditor->appendSurface(result);
            gui->switchEditor("surface");
        }

        ImGui::SameLine();

        if (ImGui::Button("Surface of Revolution")) {
            auto result = surface_of_revolution(*NURBSForSpline(selectedSplinePtr()), revSegments);
            SurfaceEditor *surfEditor = dynamic_cast<SurfaceEditor *>(gui->getEditor("surface"));
            if (surfEditor == nullptr) throw std::logic_error("Failed to access surface editor");
            surfEditor->appendSurface(result);
            gui->switchEditor("surface");
        }
        if (ImGui::InputInt("Revolution Segments", &revSegments, 1)) {
            revSegments = std::max(3, revSegments);
        }

        {
            ImGui::Separator();
            ImGui::Text("Parametrization Settings");
            bool paramChanged = false;

            int ptypeInt = int(s.paramType);
            if (ImGui::Combo("Type", &ptypeInt, "Uniform\0Chord Length\0Centripetal\0")) {
                s.paramType = Spline::ParametrizationType(ptypeInt);
                paramChanged = true;
            }

            paramChanged |= ImGui::Checkbox("Repeated End Knots", &s.repeatedEndKnots);
            ImGui::Checkbox("Recompute on Drag", &recomputeParametrizationOnDrag);

            if (paramChanged) {
                s.inferKnots();
                changed = true;
            }
        }

        {
            ImGui::Separator();
            ImGui::Text("Rational Splines (NURBS)");
            bool rational = !ordinarySpline;
            if (ImGui::Checkbox("Rational", &rational)) {
                if (rational) { // conversion to rational requested
                    splines[selection.curve] = std::make_shared<NURBS<3>>(*ordinarySpline);
                }
                else { // conversion from rational requested
                    auto sr = std::dynamic_pointer_cast<NURBS<3>>(selectedSplinePtr());
                    if (!sr) throw std::runtime_error("Curve is not rational");
                    auto snew = std::make_shared<Spline>();

                    // Conversion from NURBS loses weights...
                    snew->setControlPts(sr->getControlPts());
                    snew->knots            = sr->homogeneousSpline().knots;
                    snew->repeatedEndKnots = sr->repeatedEndKnots;
                    snew->paramType        = sr->paramType;
                    snew->degree()         = sr->degree();

                    splines[selection.curve] = snew;
                }
                changed = true;

                // we swapped out the currently selected spline...
                knotViewer(viewer).enable(selectedSplinePtr());
            }

            // WARNING: reference `s` will be dangling here if a rationality
            // conversion was performed.

            if (rational) {
                auto sr = std::dynamic_pointer_cast<NURBS<3>>(selectedSplinePtr());
                if (!sr) throw std::runtime_error("Curve is not rational");
                for (size_t i = 0; i < sr->numControlPts(); ++i) {
                    float w = sr->weight(i);
                    if (ImGui::SliderFloat(("w_" + std::to_string(i)).c_str(), &w, /* vmin = */ 1e-3, /* vmax = */ 10)) {
                        changed = true;
                        sr->setWeight(i, w);
                    }
                }

                if (ImGui::Button("Make Circle")) {
                    auto c = std::make_shared<NURBS<3>>();
                    *c = NURBSCircle<3>(1.0, circleSegments, keptSegments);
                    c->validateSizes();

                    splines[selection.curve] = c;
                    // we swapped out the currently selected spline...
                    knotViewer(viewer).enable(selectedSplinePtr());
                    changed = true;
                }
                if (ImGui::InputInt("Circle Segments", &circleSegments, 1)) {
                    circleSegments = std::max(3, circleSegments);
                    keptSegments = circleSegments;
                }
                if (ImGui::InputInt("Kept Segments", &keptSegments, 1)) {
                    circleSegments = std::max(std::max(1, keptSegments), circleSegments);
                }
            }
        }

        if (changed) updateView(viewer);
        return changed || handled;
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
        MouseRay(viewer, viewer.down_mouse_x, viewer.down_mouse_y)
            .select(splines.size(),
                [&](size_t si) { return splines[si]->getControlPts(); },
            [&](size_t ci, size_t pi) {
                selection.curve     = ci;
                selection.controlPt = pi;
                selection.dragging  = true;
                selection.dragAll = modifier == GLFW_MOD_SHIFT;
            }
        );

        if (selection.curve != oldSel) updateView(viewer);

        return selection.dragging;
    }

    bool callback_mouse_move(IGLViewer &viewer, int mouse_x, int mouse_y) override {
        if (selection.dragging) {
            if (selection.controlPt >= 0) {
                // Dragging on control point.
                auto &s = selectedSpline();
                auto cp = s.controlPt(selection.controlPt);
                Eigen::Vector3f new_pt = dragPoint(viewer, cp, mouse_x, mouse_y, editMode == EditingMode::TWO_D);

                if (selection.dragAll) {
                    for (size_t i = 0; i < s.numControlPts(); ++i)
                        s.setControlPt(i, s.controlPt(i) + new_pt - cp);
                }
                else {
                    s.setControlPt(selection.controlPt, new_pt);
                    // update parametrization on shape change?
                    if (recomputeParametrizationOnDrag)
                        s.inferKnots();
                }
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

protected:
    virtual void m_enterImpl(IGLViewer &v) override {
        knotViewer(v).attachEditor(this);
        knotViewer(v).enable(selectedSplinePtr());
    }

    virtual void m_leaveImpl(IGLViewer &v) override {
        knotViewer(v).disable();
    }
};

#endif /* end of include guard: SPLINEEDITOR_HH */
