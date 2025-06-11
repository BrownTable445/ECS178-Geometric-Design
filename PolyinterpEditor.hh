////////////////////////////////////////////////////////////////////////////////
// PolyinterpEditor.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Polynomial interpolation user interface for ECS 178.
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
//  Company:  University of California, Davis
//  Created:  10/28/2022 10:34:11
*///////////////////////////////////////////////////////////////////////////////
#ifndef POLYINTERPEDITOR_HH
#define POLYINTERPEDITOR_HH

#include "EditorBase.hh"
#include "BezierEditor.hh"
#include "polyinterp.hh"

struct DataSource {
    DataSource() { std::cout << "Constructed data source" << std::endl; }
    int selectedFunction = 0;
    int numDataPts = 3;
    float noiseMagnitude = 1e-16;
    std::vector<std::function<Eigen::Vector3f(double)>> demoFunctions = {{
        [](double t) { double x = 2 * t - 1; return Eigen::Vector3f(x, 1 - x * x, 0.0); },
        [](double t) { double x = 2 * t - 1; return Eigen::Vector3f(x, 1.0 / (1.0 + 25 * (x * x)), 0.0); },
        [](double t) { return Eigen::Vector3f(2 * cos(M_PI / 2 * t), sin(M_PI / 2 * t), 0.0); },
        [](double t) { return Eigen::Vector3f(2 * cos(2 * M_PI * t), 2 * sin(2 * M_PI * t), 0.0); }
    }};

    enum class SampleSpacing {
        UNIFORM = 0, CHEBYSHEV = 1
    };

    SampleSpacing sampleSpacing = SampleSpacing::UNIFORM;

    bool draw_menu(IGLViewer &viewer) {
        bool functionSettingsChanged = false;
        if (ImGui::BeginCombo("Point Data Source", ("Function " + std::to_string(selectedFunction)).c_str())) {
            for (size_t i = 0; i < demoFunctions.size(); ++i) {
                bool selected = (i == selectedFunction);
                if (ImGui::Selectable(("Function " + std::to_string(i)).c_str(), &selected)) {
                    selectedFunction = i;
                    functionSettingsChanged = true;
                }
                if (selected) ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        if (ImGui::InputInt("Num Data Points", &numDataPts, 1)) {
            numDataPts = std::min(std::max(numDataPts, 2), 100);
            functionSettingsChanged = true;
        }

        int ssInt = int(sampleSpacing);
        functionSettingsChanged |= ImGui::Combo("Sample Spacing", &ssInt, "Uniform\0Chebyshev\0");
        sampleSpacing = SampleSpacing(ssInt);

        functionSettingsChanged |= logarithmicSliderWithZero("Data Noise", noiseMagnitude, 1e-16, 1e-1);
        functionSettingsChanged |= ImGui::Button("Regenerate Data Points");

        return functionSettingsChanged;
    }

    std::vector<Eigen::Vector3f> generate() {
        std::vector<Eigen::Vector3f> pts(numDataPts);
        for (size_t i = 0; i < numDataPts; ++i) {
            double t;
            if (sampleSpacing == SampleSpacing::UNIFORM)
                t = double(i) / (numDataPts - 1);
            else if (sampleSpacing == SampleSpacing::CHEBYSHEV) {
                t = 0.5 + 0.5 * cos(((2 * i + 1) * M_PI) / (2 * numDataPts));
            }
            else throw std::runtime_error("Unknown sample spacing");
            pts[i] = demoFunctions[selectedFunction](t);
            pts[i].head<2>() += Eigen::Vector2f::Random() * noiseMagnitude;
        }
        return pts;
    }
};

struct PolyinterpEditor : public EditorBase {
    PolyinterpEditor(AssignmentGUI *gui, DataSource &ds, PolyInterp &pi)
        : EditorBase(gui), dataSource(ds), polynomialInterpolant(pi) {
        polynomialInterpolant.setInterpPts(dataSource.generate());
    }

    PolyInterp &polynomialInterpolant;

    enum class EvalProbeVisualization {
        NONE = 0, AITKEN = 1
    };

    bool  showAitken = false;
    int   resolution = 100;
    float eval_t     = 0.5;
    PolyInterp::EvalMethod method = PolyInterp::EvalMethod::AITKEN;

    DataSource &dataSource;

    struct Selection {
        int controlPt = 0;
        bool dragging = false;
        bool dragAll = false;
    };

    Selection selection;
    EvalProbeVisualization visMode = EvalProbeVisualization::NONE;

    virtual void updateView(IGLViewer &v) override {
        auto &plotter        = *dynamic_cast<LinePlotter         *>(v.plugins[0]);
        auto &vectorRenderer = *dynamic_cast<VectorFieldRenderer *>(v.plugins[1]);

        vectorRenderer.clear();
        plotter.clear();
        auto &vdata = v.data();
        vdata.clear();
        float dpiScale = static_cast<Viewer &>(v).getDPIScale();

        static std::vector<Eigen::Vector3f> evaluatedCurve; // re-used across updates

        Eigen::Vector4f cageColor = Eigen::Vector4f(0.0f, 0.0f, 0.0f, 1.0f);
        Eigen::Vector4f curveColor = Eigen::Vector4f(0.3, 0.6f, 1.0f, 1.0f);

        // Evaluate and draw the interpolating curve
        polynomialInterpolant.eval(resolution, method, evaluatedCurve);
        plotter.addPolyline(evaluatedCurve, /* closed = */ false, /* thickness = */ 5.0f * dpiScale, /* color = */ curveColor);

        // Evaluate at the probe point
        Eigen::Vector3f p = polynomialInterpolant.evalPt(method, eval_t);

        vdata.add_points(pointsToMatrixRows(polynomialInterpolant.getInterpPts()).cast<double>(), /* C = */ BezierEditor::cageColor(true).head<3>().transpose().cast<double>());

        vdata.add_points(p.transpose().cast<double>(), /* C = */ curveColor.head<3>().transpose().cast<double>());

        vectorRenderer.clear();
        if (visMode == EvalProbeVisualization::AITKEN) {
            Eigen::MatrixX3f V;
            Eigen::MatrixX2i E;
            polynomialInterpolant.visualizeAitken(eval_t, V, E);
            plotter.addEdges(V, E, 1.0f * dpiScale, Eigen::Vector4f(0.5f, 0.5f, 0.5f, 1.0f));
            vdata.add_points(V.cast<double>(), /* C = */ Eigen::RowVector3d(0.5f, 0.5f, 0.5f));
        }

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
           .select(1, [&](size_t ci) { return polynomialInterpolant.getInterpPts(); },
                  [&](size_t /* ci */, size_t pi) {
                      selection.controlPt = pi;
                      selection.dragging  = true;
                      selection.dragAll = modifier == GLFW_MOD_SHIFT;
                  });

        return selection.dragging;
    }

    virtual bool callback_mouse_move(IGLViewer &viewer, int mouse_x, int mouse_y) override {
        if (selection.dragging) {
            auto allPts = polynomialInterpolant.getInterpPts(); // we must edit and then set a copy of all points for proper updating...
            Eigen::Vector3f &cp = allPts[selection.controlPt];
            Eigen::Vector3f newPt = dragPoint(viewer, cp, mouse_x, mouse_y, editMode == EditingMode::TWO_D);

            if (selection.dragAll) {
                Eigen::Vector3f diff = newPt - cp;
                for (auto &p : allPts)
                    p += diff;
            }
            else cp = newPt;

            polynomialInterpolant.setInterpPts(allPts);

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

        int evalMethodInt = int(method);
        if (ImGui::Combo("Evaluation Method", &evalMethodInt, "Aitken\0Lagrange\0Lagrange Barycentric\0")) {
            method = PolyInterp::EvalMethod(evalMethodInt);
            changed = true;
        }

        if (ImGui::SliderFloat("Eval t", &eval_t, /* vmin = */ 0.0f, /* vmax = */ 1.0f))
            changed = true;

        int visModeInt = int(visMode);
        if (ImGui::Combo("Eval Pt Visualization", &visModeInt, "Point Only\0Aitken\0")) {
            visMode = EvalProbeVisualization(visModeInt);
            changed = true;
        }

        if (changed) updateView(viewer);
        return changed || handled;
    }
};

#endif /* end of include guard: POLYINTERPEDITOR_HH */
