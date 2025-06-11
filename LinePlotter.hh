////////////////////////////////////////////////////////////////////////////////
// LinePlotter.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Thick line rendering in OpenGL 3.3.
//  GL_LINE_SMOOTH + glLineWidth is no longer supported on many drivers. We use
//  an improved version of the approach originally suggested in:
//  https://stackoverflow.com/questions/60440682/drawing-a-line-in-modern-opengl
//
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
//  Company:  University of California, Davis
//  Created:  10/09/2022 17:06:09
*///////////////////////////////////////////////////////////////////////////////
#ifndef LINEPLOTTER_HH
#define LINEPLOTTER_HH

#include <igl/opengl/create_shader_program.h>
#include <igl/opengl/destroy_shader_program.h>
#include <igl/opengl/bind_vertex_attrib_array.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/ViewerPlugin.h>

#include <list>

inline std::string readFileToString(const std::string &path) {
    std::ifstream inFile(path);
    if (!inFile.is_open()) throw std::runtime_error("Couldn't open input file " + path);
    std::stringstream ss;
    ss << inFile.rdbuf();
    return ss.str();
};


struct LinePlotter : public igl::opengl::glfw::ViewerPlugin {
    // Points stored coordinate-contiguously in the columns (differing from `libigl` conventions)
    using VArray = Eigen::Matrix<float, 4, Eigen::Dynamic, Eigen::ColMajor>;

    LinePlotter() { viewer = nullptr; }

    struct LineGeometry {
        LineGeometry(const VArray &varray, bool polyline,
                     float thickness = 10.0f,
                     const Eigen::Vector4f &color = Eigen::Vector4f(1.0f, 0.0f, 0.0f, 1.0f))
            : color(color),
              thickness(thickness),
              contiguous(polyline),
              m_primitiveType(polyline ? Type::LINE_STRIP : Type::LINE) {

            m_numEdges = polyline ? varray.cols() - 2 - 1 : varray.cols() / 2;
            glGenTextures(1, &m_vertexCoordTexture);
            glGenBuffers(1,  &m_vertexCoordBuffer);

            glBindBuffer(GL_TEXTURE_BUFFER, m_vertexCoordBuffer);
            glBufferData(GL_TEXTURE_BUFFER, varray.size() * sizeof(float), varray.data(), GL_DYNAMIC_DRAW);

            glBindTexture(GL_TEXTURE_BUFFER, m_vertexCoordTexture);
            glTexBuffer  (GL_TEXTURE_BUFFER, GL_RGBA32F, m_vertexCoordBuffer);
        }

        void draw(GLuint shader) const {
            glUseProgram(shader);
            glUniform1f(glGetUniformLocation(shader, "u_thickness"),     thickness);
            glUniform1i(glGetUniformLocation(shader, "u_lineStripMode"), (m_primitiveType == Type::LINE_STRIP));
            glUniform1i(glGetUniformLocation(shader, "u_contiguous"),    contiguous);
            glUniform4fv(glGetUniformLocation(shader, "u_color"), 1, color.data());

            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_BUFFER, m_vertexCoordTexture);
            glDrawArrays(GL_TRIANGLES, 0, 6 * m_numEdges);
        }

        LineGeometry &operator=(const LineGeometry &) = delete; // copying would be problematic (breaks RAII resource management).

        ~LineGeometry() {
            glDeleteTextures(1, &m_vertexCoordTexture);
            glDeleteBuffers(1,  &m_vertexCoordBuffer);
        }

        Eigen::Vector4f color;
        float thickness;
        bool contiguous;

    private:
        size_t m_numEdges;
        enum class Type { LINE, LINE_STRIP };
        Type m_primitiveType;
        GLuint m_vertexCoordBuffer, m_vertexCoordTexture;
    };

    void clear() {
        m_lineGeometry.clear();
    }

    void addEdges(const Eigen::MatrixX3f &V, const Eigen::MatrixX2i &E,
                  float thickness = 5.0f, const Eigen::Vector4f &color = Eigen::Vector4f(0.0, 0.0, 0.0, 1.0)) {
        if (E.rows() < 1) return;
        VArray data(4, 2 * E.rows());
        for (int i = 0; i < E.rows(); ++i) {
            if ((E(i, 0) >= V.rows()) || (E(i, 1) >= V.rows()))
                throw std::runtime_error("Edge endpoint index out of bounds");
            data.col(2 * i + 0) << V.row(E(i, 0)).transpose(), 1.0;
            data.col(2 * i + 1) << V.row(E(i, 1)).transpose(), 1.0;
        }

        m_lineGeometry.emplace_back(data, false, thickness, color);
    }

    void addPolyline(const std::vector<Eigen::Vector3f> &pts, bool closed = false,
                     float thickness = 5.0f, const Eigen::Vector4f &color = Eigen::Vector4f(0.0, 0.0, 0.0, 1.0)) {
        if (pts.size() < 2) throw std::runtime_error("Polyline must have at least one edge");

        // Polylines (LINE_STRIP) must be padded by an additional vertex at their
        // beginnings and ends for proper miters. Closed polylines
        // additionally need the first vertex to be repeated at the end.
        VArray data(4, pts.size() + closed + 2); 

        // Add padding at the beginning
        if (closed) data.col(0) << pts.back(), 1.0;
        else        data.col(0) << 2 * pts[0] - pts[1], 1.0;

        for (int i = 0; i < pts.size(); ++i)
            data.col(i + 1) << pts[i], 1.0;

        if (closed) {
            data.col(pts.size() + 1) << pts[0], 1.0; // Repeat first vertex to close polyline
            data.col(pts.size() + 2) << pts[1], 1.0; // End padding
        }
        else {
            data.col(pts.size() + 1) << 2 * pts.back() - pts[pts.size() - 2], 1.0; // End padding
        }

        m_lineGeometry.emplace_back(data, true, thickness, color);
    }

    virtual void init(igl::opengl::glfw::Viewer *_viewer) override {
        assert(_viewer);
        ViewerPlugin::init(_viewer);

        igl::opengl::create_shader_program(readFileToString(SHADER_PATH "/line.vert"),
                                           readFileToString(SHADER_PATH "/line.frag"), {}, m_shader);

        glGenVertexArrays(1, &m_vao);
        glBindVertexArray(m_vao);

        // Initialize the resolution.
        post_resize(viewer->core().viewport[2], viewer->core().viewport[3]);
    }

    virtual void post_resize(int w, int h) override { m_resolution << w, h; }

    virtual Eigen::Matrix4f modelViewProjection() const {
        const Eigen::Matrix4f &proj = viewer->core().proj;
        const Eigen::Matrix4f &view = viewer->core().view;
        return proj * view;
    }

    virtual bool post_draw() override {
        glBindVertexArray(m_vao);
        glUseProgram(m_shader);
        
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // for debugging
        glUniformMatrix4fv(glGetUniformLocation(m_shader, "u_mvp"), /* count */ 1, /* transpose */ GL_FALSE, modelViewProjection().data());

        glUniform2f(glGetUniformLocation(m_shader, "u_resolution"), m_resolution[0], m_resolution[1]);
        glUniform1i(glGetUniformLocation(m_shader, "u_vcoord_tex"), 0);

        glBindVertexArray(m_vao);

        for (const auto &lg : m_lineGeometry)
            lg.draw(m_shader);

        return false;
    };

    virtual ~LinePlotter() {
        glDeleteShader(m_shader);
        glDeleteVertexArrays(1, &m_vao);
    }

private:
    Eigen::Vector2f m_resolution;
    GLuint m_shader, m_vao;
    std::list<LineGeometry> m_lineGeometry;
};
    
#endif /* end of include guard: LINEPLOTTER_HH */
