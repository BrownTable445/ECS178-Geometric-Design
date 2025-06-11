////////////////////////////////////////////////////////////////////////////////
// SurfaceEditor.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Surface manipulation user interface for ECS 178.
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
//  Company:  University of California, Davis
//  Created:  11/13/2022 19:12:20
*///////////////////////////////////////////////////////////////////////////////
#ifndef SURFACEEDITOR_HH
#define SURFACEEDITOR_HH

#include "EditorBase.hh"
#include "TriSurfEditor.hh"
#include "spline.hh"
#include "KnotViewer.hh"
#include "AssignmentGUI.hh"
#include "tensor_product_surface.hh"
#include "UVDomainWidget.hh"
#include "CurvatureVisualization.hh"

struct SurfaceEditor : public EditorBase {
    using Surf = TensorProductSurface_T<3>;
    std::vector<std::shared_ptr<Surf>> surfaces;

    CurvatureVis curvatureVis;

    UVDomainWidget domainWidget;

    SurfaceEditor(AssignmentGUI *gui) : EditorBase(gui), surfaces(1, std::make_shared<Surf>()) {
        editMode = EditingMode::THREE_D;
        concatenatedMesh.reevalMeshes(resolution, surfaces);
    }

    struct Selection {
        int surf = 0;
        int controlPt = 0; // "flattened" linear index
        bool dragging = false;
        bool dragAll = false;
    };

    Selection selection;
    int   resolution = 30;
    float eval_u     = 0.5;
    float eval_v     = 0.5;
    bool showCage    = true;
    bool showFrame   = true;

    // Concatenation of all sampled surfaces into a single mesh. Also
    // implements caching to avoid re-evaluation when changing the query point/visualization.
    struct ConcatenatedMesh {
        std::vector<Eigen::MatrixX2f> surf_U;
        std::vector<Eigen::MatrixX3f> surf_V;
        std::vector<Eigen::MatrixX3i> surf_F;

        size_t numMeshes() const {
            const size_t nm = surf_V.size();
            if (nm != surf_F.size()) throw std::logic_error("V, F mismatch");
            return nm;
        }

        void reevalMeshes(size_t resolution, const std::vector<std::shared_ptr<Surf>> &s) {
            surf_U.resize(s.size());
            surf_V.resize(s.size());
            surf_F.resize(s.size());
            for (size_t i = 0; i < s.size(); ++i)
                reevalMesh(resolution, i, *s[i]);
            markConcatenationDirty();
        }

        void reevalMesh(size_t resolution, size_t i, const Surf &s) {
            if (i >= numMeshes()) throw std::runtime_error("Mesh index out of bounds");
            s.eval(resolution, surf_U[i], surf_V[i], surf_F[i]);
            markConcatenationDirty();
        }

        void updateConcatenation(Selection &selection, const std::vector<std::shared_ptr<Surf>> &surfaces,
                                 const CurvatureVis &curvatureVis) {
            const size_t nm = numMeshes();
            int nv = 0, nf = 0;
            for (size_t i = 0; i < nm; ++i) {
                nv += surf_V[i].rows();
                nf += surf_F[i].rows();
            }
            V_concat.resize(nv, 3);
            C_concat.resize(nv, 3);
            F_concat.resize(nf, 3);

            if (curvatureVis.type == CurvatureVis::Type::None || selection.surf >= nm) {
                // No curvature to visualize.
                scalarField.resize(0);
            }

            size_t v_offset = 0;
            size_t f_offset = 0;
            for (size_t i = 0; i < nm; ++i) {
                const auto &U = surf_U[i];
                const auto &V = surf_V[i];
                const auto &F = surf_F[i];

                V_concat.middleRows(v_offset, V.rows()) = V.cast<double>();
                F_concat.middleRows(f_offset, F.rows()) = F.array() + v_offset;

                if (selection.surf == i && curvatureVis.type != CurvatureVis::Type::None) {
                    scalarField.resize(U.rows());
                    Eigen::VectorXd &curvature = scalarField;
                    const auto &s = *surfaces[i];

                    for (size_t vi = 0; vi < U.rows(); ++vi) {
                        if (curvatureVis.type == CurvatureVis::Type::Gauss)
                            curvature[vi] = s.gaussianAndMeanCurvature(U(vi, 0), U(vi, 1))[0];
                        else if (curvatureVis.type == CurvatureVis::Type::Mean)
                            curvature[vi] = s.gaussianAndMeanCurvature(U(vi, 0), U(vi, 1))[1];
                        else if (curvatureVis.type == CurvatureVis::Type::Kappa1)
                            curvature[vi] = s.principalCurvatures(U(vi, 0), U(vi, 1))[0];
                        else
                            curvature[vi] = s.principalCurvatures(U(vi, 0), U(vi, 1))[1];
                    }

                    C_concat.middleRows(v_offset, V.rows()) = curvatureVis.colorsForCurvature(curvature);
                }
                else {
                    C_concat.middleRows(v_offset, V.rows()).rowwise() =
                        BezierEditor::curveColor(selection.surf == i).head<3>().cast<double>().transpose();
                }

                v_offset += V.rows();
                f_offset += F.rows();
            }
            m_concatenationDirty = false;
        }

        bool concatenationDirty() const { return m_concatenationDirty; }
        void markConcatenationDirty() { m_concatenationDirty = true; }

        Eigen::MatrixXd V_concat;
        Eigen::MatrixXi F_concat;
        Eigen::MatrixXd C_concat;

        Eigen::VectorXd scalarField;

    private:
        bool m_concatenationDirty = true;
    };

    ConcatenatedMesh concatenatedMesh;

          Surf &selectedSurf()       { return *surfaces[selection.surf]; }
    const Surf &selectedSurf() const { return *surfaces[selection.surf]; }

    static KnotViewer &knotViewer(IGLViewer &v) {
        return *dynamic_cast<KnotViewer*>(v.plugins[2]);
    }

    std::shared_ptr<Surf> selectedSurfPtr() {
        if (selection.surf > surfaces.size()) throw std::runtime_error("No valid surface selection");
        return surfaces[selection.surf];
    }

    virtual void updateView(IGLViewer &v) override {
        auto &plotter        = *dynamic_cast<LinePlotter         *>(v.plugins[0]);
        auto &vectorRenderer = *dynamic_cast<VectorFieldRenderer *>(v.plugins[1]);
        vectorRenderer.clear();
        plotter.clear();
        auto &vdata = v.data();

        // Evaluate and draw the concatenated surfaces
        if (concatenatedMesh.concatenationDirty()) {
            concatenatedMesh.updateConcatenation(selection, surfaces, curvatureVis);
            vdata.clear();
            vdata.set_mesh(concatenatedMesh.V_concat, concatenatedMesh.F_concat);
            vdata.set_colors(concatenatedMesh.C_concat);
        }

        float dpiScale = static_cast<Viewer &>(v).getDPIScale();

        vdata.clear_points();
        for (size_t i = 0; i < surfaces.size(); ++i)
            drawSurfCage(v, *surfaces[i], i == size_t(selection.surf), showCage);

        // Evaluate at the probe point
        const auto &s = selectedSurf();
        eval_u = clamp(eval_u, s.u_spline.domainStart(), s.u_spline.domainEnd());
        eval_v = clamp(eval_v, s.v_spline.domainStart(), s.v_spline.domainEnd());
        Eigen::Vector3f p = s.evalPt(eval_u, eval_v);
        vdata.add_points(p.transpose().cast<double>(), /* C = */ BezierEditor::curveColor(true).head<3>().transpose().cast<double>());

        // Visualize tangent and normal vectors.
        if (showFrame) {
            auto J = s.jacobian(eval_u, eval_v);

            VectorFieldRenderer::PtArray frenetPts(3, 3), frenetVectors(3, 3);
            VectorFieldRenderer::ColorArray colors(3, 4);
            frenetPts << p.transpose(),
                         p.transpose(),
                         p.transpose();
            frenetVectors << J.col(0).normalized().transpose(),
                             J.col(1).normalized().transpose(),
                             s.normal(eval_u, eval_v).transpose();
            colors << 1, 0, 0, 1,
                      0, 1, 0, 1,
                      0, 0, 1, 1;
            vectorRenderer.setField(frenetPts, frenetVectors, colors);
        }

        vdata.point_size = 10.0f * dpiScale;
    }

    void drawSurfCage(IGLViewer &v, Surf &s, bool selected, bool showCage) {
        auto &plotter  = *dynamic_cast<LinePlotter *>(v.plugins[0]);
        float dpiScale = static_cast<Viewer &>(v).getDPIScale();

        if (showCage) {
            // Draw the control polygon
            const auto &d = s.getControlPtsFlattened();

            auto &cp = s.controlPts;
            for (size_t i = 0; i < cp.rows(); ++i)
                plotter.addPolyline(s.getControlPtRow(i), /* closed = */ false, /* thickness = */ 2.0f * dpiScale, /* color = */ BezierEditor::cageColor(selected));
            for (size_t j = 0; j < cp.cols(); ++j)
                plotter.addPolyline(s.getControlPtCol(j), /* closed = */ false, /* thickness = */ 2.0f * dpiScale, /* color = */ BezierEditor::cageColor(selected));

            Eigen::MatrixXd pts = s.getControlPtsFlattened().cast<double>();
            Eigen::MatrixXd ptColors;
            if (selected == false)
                ptColors = BezierEditor::cageColor(selected).head<3>().transpose().cast<double>();
            else {
                ptColors.resize(pts.rows(), 3);
                ptColors.rowwise() = BezierEditor::cageColor(selected).head<3>().transpose().cast<double>();
                ptColors.row(selection.controlPt) << 1, 0, 0;
            }
            v.data().add_points(pts, ptColors);
        }
    }

    void appendSurface(const std::shared_ptr<Surf> &s) {
        surfaces.emplace_back(s);
        selection.surf = surfaces.size() - 1;
        selection.controlPt = 0;
        concatenatedMesh.reevalMeshes(resolution, surfaces);
    }

    bool draw_menu(IGLViewer &viewer) override {
        bool handled = false, changed = false;
        if (ImGui::InputInt("Sampling Resolution", &resolution, 1)) {
            resolution = std::min(std::max(resolution, 2), 1000);
            concatenatedMesh.reevalMeshes(resolution, surfaces);
        }

        bool showWireframe = viewer.data().show_lines;
        if (ImGui::Checkbox("Show Wireframe", &showWireframe))
            viewer.data().show_lines = showWireframe;

        changed |= ImGui::Checkbox("Show Control Polyline", &showCage);

        ImGui::Separator();

        {
            bool changedKnots, changedEvalPt;
            domainWidget.run(selectedSurf(), eval_u, eval_v, changedKnots, changedEvalPt);
            if (changedKnots)
                concatenatedMesh.reevalMesh(resolution, selection.surf, selectedSurf());
            changed |= changedEvalPt;
        }

        ImGui::Separator();
        ImGui::Text("Operations");

        if (ImGui::Button("Add Surf", ImVec2(ImGui::GetContentRegionAvail().x*0.5f, 0.0f))) {
            Eigen::Vector3f t = Eigen::Vector3f::Zero();
            t.topRows(2).setRandom();
            appendSurface(std::make_shared<Surf>(t));
        }
        ImGui::SameLine();

        auto conditionallyDisabledButton = [](bool disabled, const char *label, const ImVec2 &size = ImVec2(0,0)) {
            if (disabled)
                ImGui::BeginDisabled();
            bool result = ImGui::Button(label, size);
            if (disabled) ImGui::EndDisabled();
            return result;
        };

        if (conditionallyDisabledButton(surfaces.size() <= 1, "Delete Surface", ImVec2(ImGui::GetContentRegionAvail().x, 0.0f))) {
            surfaces.erase(surfaces.begin() + selection.surf);
            selection.surf = std::min(size_t(selection.surf), surfaces.size() - 1);
            selection.controlPt = 0;
            concatenatedMesh.reevalMeshes(resolution, surfaces);
        }

        if (ImGui::Button("Insert U Knot")) {
            selectedSurf().insertUKnot(eval_u);
            concatenatedMesh.reevalMesh(resolution, selection.surf, selectedSurf());
        }

        ImGui::SameLine();

        if (ImGui::Button("Insert V Knot")) {
            selectedSurf().insertVKnot(eval_v);
            concatenatedMesh.reevalMesh(resolution, selection.surf, selectedSurf());
        }

        auto &s = selectedSurf();
        {
            ImGui::Separator();
            ImGui::Text("Parametrization Settings");

            if (ImGui::Checkbox("Repeated End Knots", &s.u_spline.repeatedEndKnots)) {
                s.v_spline.repeatedEndKnots = s.u_spline.repeatedEndKnots;
                s.u_spline.inferKnots();
                s.v_spline.inferKnots();
                concatenatedMesh.reevalMesh(resolution, selection.surf, selectedSurf());
            }
        }

        {
            size_t i, j;
            s.unflattenIndex(selection.controlPt, i, j);
            if (ImGui::SliderFloat(("Weight w_" + std::to_string(i) + "," + std::to_string(j)).c_str(), &s.weights(i, j), /* vmin = */ 1e-3, /* vmax = */ 10))
                concatenatedMesh.reevalMesh(resolution, selection.surf, selectedSurf());
        }

        ImGui::Separator();
        ImGui::Text("Differential Geometry");
        changed |= ImGui::Checkbox("Show Eval Pt. Frame", &showFrame);
        if (curvatureVis.widget(concatenatedMesh.scalarField)) {
            concatenatedMesh.markConcatenationDirty(); // needed for color update
        }

        auto kappa = s.principalCurvatures(eval_u, eval_v);
        auto curvature = s.gaussianAndMeanCurvature(eval_u, eval_v);

        ImGui::Text("Curvatures at Eval Probe: (%f, %f)", eval_u, eval_v);
        ImGui::Text("\t\tkappa_1: %0.4f", kappa[0]);
        ImGui::Text("\t\tkappa_2: %0.4f", kappa[1]);
        ImGui::Text("\t\tK: %0.4f, H: %0.4f", curvature[0], curvature[1]);

        ImGui::Separator();
        if (ImGui::Button("Copy to Triangle Mesh Editor")) {
            TriSurfEditor *trisurfEditor = dynamic_cast<TriSurfEditor *>(gui->getEditor("trisurf"));
            if (trisurfEditor == nullptr) throw std::logic_error("Failed to access TriSurfEditor");
            trisurfEditor->surface.setMesh(concatenatedMesh.surf_V[selection.surf].cast<double>(),
                                           concatenatedMesh.surf_F[selection.surf]);
            gui->switchEditor("trisurf", /* ignorePrevCameraSettings = */ true);
        }

        changed |= concatenatedMesh.concatenationDirty();
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
        if (!showCage) return false; // Only drag control pts if they're shown.

        MouseRay ray(viewer, viewer.down_mouse_x, viewer.down_mouse_y);

        float selectionRadius = 0.05f;
        float minDist = selectionRadius;

        // First detect clicks on the control points.
        selection.dragging = false;
        int oldSel = selection.surf;
        int oldSelPt = selection.controlPt;
        MouseRay(viewer, viewer.down_mouse_x, viewer.down_mouse_y).select(
            surfaces.size(),
            [&](size_t si) { return surfaces[si]->getControlPtsFlattened(); },
            [&](size_t si, size_t pi) {
                selection.surf      = si;
                selection.controlPt = pi;
                selection.dragging  = true;
                selection.dragAll   = modifier == GLFW_MOD_SHIFT;
            }
        );

        if ((selection.surf != oldSel) || (selection.controlPt != oldSelPt)) {
            concatenatedMesh.markConcatenationDirty(); // needed for color update
            updateView(viewer);
        }

        return selection.dragging;
    }

    bool callback_mouse_move(IGLViewer &viewer, int mouse_x, int mouse_y) override {
        if (selection.dragging) {
            auto draggedPt = [&](Eigen::Vector3f pt) {
                const auto &c = viewer.core();
                // Determine screen-space depth of control point (to be preserved by drag)
                Eigen::Vector3f win = igl::project(pt, c.view, c.proj, c.viewport);

                win.head<2>() << mouse_x, c.viewport(3) - mouse_y; // GLFW's y coordinate is flipped vs unproject_ray's
                pt = igl::unproject(win, c.view, c.proj, c.viewport);

                return pt;
            };

            if (selection.controlPt >= 0) {
                // Dragging on control point.
                auto cp = selectedSurf().getControlPtsFlattened().row(selection.controlPt);
                Eigen::Vector3f new_pt = draggedPt(cp);
                if (selection.dragAll)
                    selectedSurf().getControlPtsFlattened().rowwise() += (new_pt.transpose() - cp).eval();
                else cp = new_pt;

                concatenatedMesh.reevalMesh(resolution, selection.surf, selectedSurf());
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
        domainWidget.reset();
        concatenatedMesh.markConcatenationDirty();
    }
};

#endif /* end of include guard: SURFACEEDITOR_HH */
