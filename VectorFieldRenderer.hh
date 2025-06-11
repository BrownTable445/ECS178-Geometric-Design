////////////////////////////////////////////////////////////////////////////////
// VectorFieldRenderer.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Visualize a vector field with instanced cylindrical arrow geometry.
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
//  Company:  University of California, Davis
//  Created:  10/10/2022 11:38:04
*///////////////////////////////////////////////////////////////////////////////
#ifndef VECTORFIELDRENDERER_HH
#define VECTORFIELDRENDERER_HH
#include "LinePlotter.hh"
#include "Viewer.hh"
#include "GLErrors.hh"
#include <cmath>

struct VectorFieldRenderer : public igl::opengl::glfw::ViewerPlugin {
    using    PtArray = Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>;
    using ColorArray = Eigen::Matrix<float, Eigen::Dynamic, 4, Eigen::RowMajor>;
    using   IdxArray = Eigen::Matrix<int  , Eigen::Dynamic, 3, Eigen::RowMajor>;

    VectorFieldRenderer() { viewer = nullptr; }

    virtual void init(igl::opengl::glfw::Viewer *_viewer) override {
        assert(_viewer);
        ViewerPlugin::init(_viewer);

        igl::opengl::create_shader_program(readFileToString(SHADER_PATH "/vector_field.vert"),
                                           readFileToString(SHADER_PATH "/vector_field.frag"), {}, m_shader);
        glGenVertexArrays(1, &m_vao);
        glBindVertexArray(m_vao);

        PtArray V, N;
        IdxArray F;
        arrowGeometry(0.4, 0.12, 0.025, /* numSegments = */ 20, V, N, F);
        m_numTrisPerArrow = F.rows();

        glGenBuffers(3, m_arrowGeomBuffers.data()); // V,  N, F
        glBindBuffer(GL_ARRAY_BUFFER, m_arrowGeomBuffers[0]);
        glBufferData(GL_ARRAY_BUFFER, V.size() * sizeof(float), V.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(0, V.cols() /* components per vertex */, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(0);

        glBindBuffer(GL_ARRAY_BUFFER, m_arrowGeomBuffers[1]);
        glBufferData(GL_ARRAY_BUFFER, N.size() * sizeof(float), N.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(1, N.cols() /* components per vertex */, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(1);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_arrowGeomBuffers[2]);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, F.size() * sizeof(int), F.data(), GL_STATIC_DRAW);

        glGenBuffers(3, m_instancedAttributeBuffers.data()); // pos, vec, color
        glCheckErrorWarn();

    }

    virtual bool post_draw() override {
        glUseProgram(m_shader);
        glBindVertexArray(m_vao);

        const Eigen::Matrix4f &proj = viewer->core().proj;
        const Eigen::Matrix4f &view = viewer->core().view;
        m_targetDepth = -view(2, 3) / view(3, 3);

        glUniformMatrix4fv(glGetUniformLocation(m_shader, "u_modelViewMatrix" ), /* count */ 1, /* transpose */ GL_FALSE, view.data());
        glUniformMatrix4fv(glGetUniformLocation(m_shader, "u_projectionMatrix"), /* count */ 1, /* transpose */ GL_FALSE, proj.data());
        // WARNING: doing the `transpose()` with Eigen below won't work unless we explicitly
        // evaluate to an Eigen::Matrix3f since `.transpose()` evaluates to a **row-major** temporary.
        glUniformMatrix3fv(glGetUniformLocation(m_shader, "u_normalMatrix"    ), /* count */ 1, /* transpose */ GL_TRUE, view.topLeftCorner<3, 3>().inverse().eval().data());

        glUniform1f(glGetUniformLocation(m_shader, "u_arrowAlignment"         ), int(alignment) / 2.0f);
        glUniform1f(glGetUniformLocation(m_shader, "u_arrowRelativeScreenSize"), arrowSizePx_x / viewer->core().viewport[2]);
        glUniform1f(glGetUniformLocation(m_shader, "u_targetDepth"            ), m_targetDepth);

        glUniform3f(glGetUniformLocation(m_shader, "u_lightEyePos"),      0.0f, 1.0f, 5.0f); // behind and above viewer
        glUniform3f(glGetUniformLocation(m_shader, "u_diffuseIntensity"), 0.6f, 0.6f, 0.6f);
        glUniform3f(glGetUniformLocation(m_shader, "u_ambientIntensity"), 0.5f, 0.5f, 0.5f);

        glDrawElementsInstanced(GL_TRIANGLES, 3 * m_numTrisPerArrow, GL_UNSIGNED_INT, NULL, m_numArrows);
        glCheckErrorWarn();
        
        return false;
    }

    void clear() {
        m_numArrows = 0;
    }

    void setField(const PtArray &P, const PtArray &V, const ColorArray &C) {
        m_numArrows = P.rows();
        if ((V.rows() != m_numArrows) || (C.rows() != m_numArrows)) throw std::runtime_error("Input size mismatch");

        glBindVertexArray(m_vao);
        glBindBuffer(GL_ARRAY_BUFFER, m_instancedAttributeBuffers[0]);
        glBufferData(GL_ARRAY_BUFFER, P.size() * sizeof(float), P.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(2, P.cols() /* components per vertex */, GL_FLOAT, GL_FALSE, 0, 0);
        glVertexAttribDivisor(2, 1); // instanced

        glBindBuffer(GL_ARRAY_BUFFER, m_instancedAttributeBuffers[1]);
        glBufferData(GL_ARRAY_BUFFER, V.size() * sizeof(float), V.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(3, V.cols() /* components per vertex */, GL_FLOAT, GL_FALSE, 0, 0);
        glVertexAttribDivisor(3, 1); // instanced

        glBindBuffer(GL_ARRAY_BUFFER, m_instancedAttributeBuffers[2]);
        glBufferData(GL_ARRAY_BUFFER, C.size() * sizeof(float), C.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(4, C.cols() /* components per vertex */, GL_FLOAT, GL_FALSE, 0, 0);
        glVertexAttribDivisor(4, 1); // instanced

        glEnableVertexAttribArray(2);
        glEnableVertexAttribArray(3);
        glEnableVertexAttribArray(4);

        glCheckErrorWarn();
    }

    // Which part of the arrow glyph to position at the arrow location.
    enum class GlyphAlignment : int {
        TAIL = 0, CENTER = 1, TIP = 2
    };

    GlyphAlignment alignment = GlyphAlignment::TAIL;

    inline void arrowGeometry(float headHeight, float headRadius, float shaftRadius, int ns, PtArray &V, PtArray &N, IdxArray &F) {
        PtArray profile(5, 3);
        profile << 0, 0, 0,
                   0, shaftRadius, 0,
                   1 - headHeight, shaftRadius, 0,
                   1 - headHeight, headRadius, 0,
                   1, 0, 0;
        Eigen::Vector3f coneNormal(headRadius, headHeight, 0);
        coneNormal.normalize();

        PtArray normals(5, 3);

        normals << -1, 0, 0,
                    0, 1, 0,
                    0, 1, 0,
                    coneNormal.transpose(),
                    coneNormal.transpose();
        float a = 2 * M_PI / ns; // Rotation angle increment
        Eigen::Matrix3f R;       // Rotation around the arrow (x) axis
        R << 1,      0,       0,
             0, cos(a), -sin(a),
             0, sin(a),  cos(a);

        size_t npp = profile.rows();
        // Each strip of the revolution consists of one triangle at each end + two triangles to triangulate each inner quad
        size_t nst = 2 * (npp - 1);
        IdxArray stripTris(nst, 3);
        stripTris.row(0) << 0, 1, npp + 1;                          // Tri from first end
        stripTris.row(1) << npp + (npp - 2), npp - 2, npp - 1;      // Tri from last end
        for (size_t i = 1; i < npp - 1; ++i) {
            stripTris.row(2 * i    ) << i,     i + 1,         npp + i;
            stripTris.row(2 * i + 1) << i + 1, npp + (i + 1), npp + i;
        }

        V.resize(ns * npp, 3);
        N.resize(ns * npp, 3);
        F.resize(ns * nst, 3);

        for (size_t i = 0; i < ns; ++i) {
            size_t vstart = i * npp;
            V.middleRows(vstart, npp) = profile;
            N.middleRows(vstart, npp) = normals;

            F.middleRows(i * nst, nst) = (stripTris.array() + vstart).unaryExpr([&](const int j) -> int { return j % V.rows(); }).matrix();
            profile = profile * R;
            normals = normals * R;
        }
    }

    ~VectorFieldRenderer() {
        glDeleteShader(m_shader);
        glDeleteVertexArrays(1, &m_vao);
        glDeleteBuffers(3, m_arrowGeomBuffers.data());
        glDeleteBuffers(3, m_instancedAttributeBuffers.data());
    }

    float arrowSizePx_x = 72.0f;

private:
    GLuint m_shader, m_vao;
    std::array<GLuint, 3> m_instancedAttributeBuffers;
    std::array<GLuint, 3> m_arrowGeomBuffers;
    float m_targetDepth = 5.0f; // based
    size_t m_numTrisPerArrow, m_numArrows = 0;
};

#endif /* end of include guard: VECTORFIELDRENDERER_HH */
