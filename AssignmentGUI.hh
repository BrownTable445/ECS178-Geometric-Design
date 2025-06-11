////////////////////////////////////////////////////////////////////////////////
// AssignmentGUI.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Wrap a viewer and the associated editing modes available for an assignment.
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
//  Company:  University of California, Davis
//  Created:  11/18/2022 22:06:28
*///////////////////////////////////////////////////////////////////////////////
#ifndef ASSIGNMENTGUI_HH
#define ASSIGNMENTGUI_HH
#include "Viewer.hh"
#include <igl/project.h>
#include <imgui.h>
#include <cmath>
#include <memory>
#include <algorithm>
#include "EditorBase.hh"
#include "LinePlotter.hh"
#include "KnotViewer.hh"
#include "VectorFieldRenderer.hh"

struct AssignmentGUI {
    AssignmentGUI(const std::string &windowTitle,
                  const std::function<void(AssignmentGUI &)> &registerEditors,
                  const std::string &initialEditor)
        : editorName(initialEditor),
          viewer(windowTitle, std::vector<IGLViewerPlugin *>{&plotter, &vectorRenderer, &knotViewer})
    {
        registerEditors(*this);

        // Ensure that each editor is displayed with the default camera settings
        // when switched to the first time.
        for (auto &e : editors) e.second->cacheCameraSettings(viewer);

        editor = editors.at(initialEditor).get();
        knotViewer.attachEditor(editor);
        vectorRenderer.arrowSizePx_x *= viewer.getDPIScale();
        viewer.core().background_color << 0.9f, 0.9f, 0.9f, 1.0f;
        editor->enter(viewer);
        editor->updateView(viewer);

        viewer.menu.callback_draw_viewer_menu = [&]() {
            bool handled = false;

            if (ImGui::BeginCombo("Editor", editorGuiNames.at(editorName).c_str())) {
                for (auto &e : editorGuiNames) {
                    bool selected = e.first == editorName;
                    if (ImGui::Selectable(e.second.c_str(), &selected)) {
                        switchEditor(e.first);
                        handled = true;
                    }
                    if (selected) ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }

            return handled || editor->draw_menu(viewer);
        };

        // Route event callbacks to the current editor.
        viewer.callback_key_down    = [&](IGLViewer &viewer, int a, int b) -> bool { return editor->callback_key_pressed(viewer, a, b); };
        viewer.callback_key_pressed = [&](IGLViewer &viewer, int a, int b) -> bool { return editor->callback_key_pressed(viewer, a, b); };
        viewer.callback_mouse_down  = [&](IGLViewer &viewer, int a, int b) -> bool { return editor->callback_mouse_down (viewer, a, b); };
        viewer.callback_mouse_move  = [&](IGLViewer &viewer, int a, int b) -> bool { return editor->callback_mouse_move (viewer, a, b); };
        viewer.callback_mouse_up    = [&](IGLViewer &viewer, int a, int b) -> bool { return editor->callback_mouse_up   (viewer, a, b); };
    }

    void switchEditor(const std::string &name, bool ignorePrevCameraSettings = false) {
        EditorBase *oldEditor = editor;
        editorName = name;
        editor = getEditor(editorName);
        editor->updateView(viewer);

        // Notify old editor and new editor of the switch.
        if (oldEditor != editor) {
            oldEditor->leave(viewer);
            editor->enter(viewer, ignorePrevCameraSettings);
        }
    }

    EditorBase *getEditor(const std::string &editorName) {
        return editors.at(editorName).get();
    }

    void run() {
        try {
            viewer.run();
        }
        catch (std::runtime_error &e) {
            std::cout << "Exiting due to error: " << e.what() << std::endl;
        }
    }

    EditorBase *editor;
    std::string editorName;

    std::map<std::string, std::string> editorGuiNames;
    std::map<std::string, std::unique_ptr<EditorBase>> editors;

    LinePlotter plotter;
    VectorFieldRenderer vectorRenderer;
    KnotViewer knotViewer;

    Viewer viewer;
};

#endif /* end of include guard: ASSIGNMENTGUI_HH */
