#ifndef EDITORBASE_HH
#define EDITORBASE_HH

#include "LinePlotter.hh"
#include "VectorFieldRenderer.hh"
#include "bezier.hh"
#include "imgui.h"
#include <memory>

inline bool logarithmicSliderWithZero(const char *label, float &val, float vmin, float vmax) {
    char val_str[64];
    snprintf(val_str, 64, "%.3e", val);
    val = std::max(val, vmin);
    val = log(val);
    bool ret = ImGui::SliderFloat(label, &val, log(vmin), log(vmax), val_str);
    val = exp(val);
    if (val < vmin * (1 + 1e-6)) val = 0.0f;
    return ret;
}

struct AssignmentGUI;

struct EditorBase {
    enum class EditingMode {
        TWO_D = 0, THREE_D = 1
    };
    EditingMode editMode = EditingMode::TWO_D;

    EditorBase(AssignmentGUI *gui) : gui(gui) { }

    virtual void updateView(IGLViewer &v) = 0;
    virtual bool callback_key_pressed(IGLViewer &viewer, unsigned int key, int modifiers) = 0;
    virtual bool callback_mouse_down(IGLViewer &viewer, int button, int modifier) = 0;
    virtual bool callback_mouse_move(IGLViewer &viewer, int mouse_x, int mouse_y) = 0;
    virtual bool callback_mouse_up(IGLViewer &viewer, int button, int modifier) = 0;
    virtual bool draw_menu(IGLViewer &viewer) = 0;

    // Struct for saving the libigl viewer's camera setting to enable restoring view
    // upon switching between editors.
    struct CameraSettings {
        CameraSettings(const IGLViewer &viewer) {
            const auto &c = viewer.core();
            camera_base_zoom        = c.camera_base_zoom;
            camera_zoom             = c.camera_zoom;
            orthographic            = c.orthographic;
            camera_base_translation = c.camera_base_translation;
            camera_translation      = c.camera_translation;
            camera_eye              = c.camera_eye;
            camera_up               = c.camera_up;
            camera_center           = c.camera_center;
            camera_view_angle       = c.camera_view_angle;
            camera_dnear            = c.camera_dnear;
            camera_dfar             = c.camera_dfar;
            trackball_angle         = c.trackball_angle;
        }

        void apply(IGLViewer &viewer) const {
            auto &c = viewer.core();
            c.camera_base_zoom        = camera_base_zoom;
            c.camera_zoom             = camera_zoom;
            c.orthographic            = orthographic;
            c.camera_base_translation = camera_base_translation;
            c.camera_translation      = camera_translation;
            c.camera_eye              = camera_eye;
            c.camera_up               = camera_up;
            c.camera_center           = camera_center;
            c.camera_view_angle       = camera_view_angle;
            c.camera_dnear            = camera_dnear;
            c.camera_dfar             = camera_dfar;
            c.trackball_angle         = trackball_angle;
        }

        float camera_base_zoom;
        float camera_zoom;
        bool orthographic;
        Eigen::Vector3f camera_base_translation;
        Eigen::Vector3f camera_translation;
        Eigen::Vector3f camera_eye;
        Eigen::Vector3f camera_up;
        Eigen::Vector3f camera_center;
        Eigen::Quaternionf trackball_angle;
        float camera_view_angle;
        float camera_dnear;
        float camera_dfar;
    };

    // Notification of the GUI switching to/from this editor.
    void enter(IGLViewer &viewer, bool ignorePrevCameraSettings = false) {
        m_enterImpl(viewer);
        if (m_prevCameraSettings && !ignorePrevCameraSettings)
            m_prevCameraSettings->apply(viewer);
        m_prevCameraSettings.reset();
        applyEditingMode(viewer);
    }

    void cacheCameraSettings(const IGLViewer &viewer) {
        m_prevCameraSettings = std::make_unique<CameraSettings>(viewer);
    }

    void leave(IGLViewer &viewer) {
        m_leaveImpl(viewer);
        cacheCameraSettings(viewer);
    }

    void applyEditingMode(IGLViewer &viewer) {
        if (editMode == EditingMode::TWO_D) {
            viewer.core().trackball_angle = Eigen::Quaternionf::Identity();
            viewer.core().set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_NO_ROTATION);
        }
        else {
            viewer.core().set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_TRACKBALL);
        }
    }

    virtual ~EditorBase() { }

protected:
    virtual void m_enterImpl(IGLViewer &/* viewer */) { }
    virtual void m_leaveImpl(IGLViewer &/* viewer */) { }

    AssignmentGUI *gui;
    std::unique_ptr<CameraSettings> m_prevCameraSettings;
};

#endif /* end of include guard: EDITORBASE_HH */
