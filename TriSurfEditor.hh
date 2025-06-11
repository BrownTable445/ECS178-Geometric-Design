////////////////////////////////////////////////////////////////////////////////
// TriSurfEditor.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  A basic interface for inspecting triangle meshes for ECS178.
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
//  Company:  University of California, Davis
//  Created:  03/13/2023 09:57:33
*///////////////////////////////////////////////////////////////////////////////
#ifndef TRISURFEDITOR_HH
#define TRISURFEDITOR_HH

#include "EditorBase.hh"
#include "AssignmentGUI.hh"
#include "BezierEditor.hh"
#include "triangle_surf.hh"
#include <igl/file_dialog_open.h>
#include <igl/file_dialog_save.h>
#include "CurvatureVisualization.hh"

struct TriSurfEditor : public EditorBase {
    using Surf = TriangleSurface;
    Surf surface;

    CurvatureVis curvatureVis;

    TriSurfEditor(AssignmentGUI *gui) : EditorBase(gui) {
        editMode = EditingMode::THREE_D;
    }

    struct Selection {
        int controlPt = 0; // "flattened" linear index
        bool dragging = false;
        bool dragAll = false;
    };

    Selection selection;
    bool showPts    = true;
    bool showNormals   = true;
    Eigen::VectorXd curvatures;

    virtual void updateView(IGLViewer &v) override {
        auto &plotter        = *dynamic_cast<LinePlotter         *>(v.plugins[0]);
        auto &vectorRenderer = *dynamic_cast<VectorFieldRenderer *>(v.plugins[1]);
        vectorRenderer.clear();
        plotter.clear();
        auto &vdata = v.data();

        float dpiScale = static_cast<Viewer &>(v).getDPIScale();

        vdata.clear();
        vdata.set_mesh(surface.V, surface.F);

        if (curvatureVis.type == CurvatureVis::Type::Gauss) {
            curvatures = surface.gaussianCurvatures();
        }
        else if (curvatureVis.type == CurvatureVis::Type::Mean) {
            curvatures = surface.meanCurvatures();
        }
        else if (curvatureVis.type == CurvatureVis::Type::Kappa1) {
            curvatures = surface.kappa_1();
        }
        else if (curvatureVis.type == CurvatureVis::Type::Kappa2) {
            curvatures = surface.kappa_2();
        }
        else {
            curvatures.resize(0);
            vdata.set_colors(BezierEditor::curveColor(true).cast<double>().transpose());
        }
        if (curvatures.size()) {
            vdata.set_colors(curvatureVis.colorsForCurvature(curvatures));
        }
        
        if (showNormals) {
            VectorFieldRenderer::PtArray pts(surface.V.cast<float>()), vecs(surface.V.rows(), 3);
            vecs = 0.25 * surface.normals().cast<float>();
            VectorFieldRenderer::ColorArray colors(surface.V.rows(), 4);
            colors.col(0).setOnes();
            colors.middleCols(1, 2).setZero();
            colors.col(3).setOnes();
            vectorRenderer.setField(pts, vecs, colors);
        }

        if (showPts) {
            v.data().add_points(surface.V, Eigen::RowVector3d(0, 0, 0));
        }

        vdata.point_size = 4.0f * dpiScale;
    }

    bool draw_menu(IGLViewer &viewer) override {
        bool handled = false, changed = false;
        bool showWireframe = viewer.data().show_lines;
        if (ImGui::Checkbox("Show Wireframe", &showWireframe))
            viewer.data().show_lines = showWireframe;

        changed |= ImGui::Checkbox("Enable Editing", &showPts);

        if (ImGui::Button("Load Mesh")) {
            std::string path = igl::file_dialog_open();
            try {
                surface.load(path);
                changed = true;
            }
            catch (...) {
                std::cout << "Failed to open surface: " << path << std::endl;
            }
        }

        ImGui::SameLine();

        if (ImGui::Button("Save Mesh")) {
            std::string path = igl::file_dialog_save();
            try {
                surface.save(path);
            }
            catch (...) {
                std::cout << "Failed to save surface to path: " << path << std::endl;
            }
        }

        ImGui::Separator();
        ImGui::Text("Differential Geometry");
        changed |= ImGui::Checkbox("Show Vertex Normals", &showNormals);

        changed |= curvatureVis.widget(curvatures);

        if (changed) updateView(viewer);
        return changed || handled;
    }

    bool callback_key_pressed(IGLViewer &viewer, unsigned int key, int modifiers) override {
        return false;
    }

    bool callback_mouse_down(IGLViewer &viewer, int button, int modifier) override {
        if (button == int(IGLViewer::MouseButton::Right))
            return false;
        if (!showPts) return false; // Only drag vertices if they're shown.

        int oldSelPt = selection.controlPt;
        MouseRay(viewer, viewer.down_mouse_x, viewer.down_mouse_y).select(
            1, [&](size_t /* si */) { return surface.V; },
            [&](size_t si, size_t pi) {
                selection.controlPt = pi;
                selection.dragging  = true;
                selection.dragAll = modifier == GLFW_MOD_SHIFT;
            }
        );

        if (selection.controlPt != oldSelPt) {
            updateView(viewer);
        }

        return selection.dragging;
    }

    bool callback_mouse_move(IGLViewer &viewer, int mouse_x, int mouse_y) override {
        if (selection.dragging) {
            auto draggedPt = [&](Eigen::Vector3d pt) {
                const auto &c = viewer.core();
                // Determine screen-space depth of control point (to be preserved by drag)
                Eigen::Vector3f win = igl::project(pt.cast<float>().eval(), c.view, c.proj, c.viewport);

                win.head<2>() << mouse_x, c.viewport(3) - mouse_y; // GLFW's y coordinate is flipped vs unproject_ray's
                pt = igl::unproject(win, c.view, c.proj, c.viewport).cast<double>();

                return pt;
            };

            if (selection.controlPt >= 0) {
                // Dragging on control point.
                auto cp = surface.V.row(selection.controlPt);
                Eigen::Vector3d new_pt = draggedPt(cp);
                if (selection.dragAll)
                    surface.V.rowwise() += (new_pt.transpose() - cp).eval();
                else cp = new_pt;
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
};

#endif /* end of include guard: TRISURFEDITOR_HH */
