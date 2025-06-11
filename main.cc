////////////////////////////////////////////////////////////////////////////////
// main.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
// Main executable/user interface for ECS 178 Assignment 5
// Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
// Company:  University of California, Davis
// Created:  10/09/2022 17:06:09
*///////////////////////////////////////////////////////////////////////////////
#include "Viewer.hh"
#include <igl/project.h>
#include <cmath>
#include <memory>
#include <algorithm>
#include "LinePlotter.hh"
#include "KnotViewer.hh"
#include "VectorFieldRenderer.hh"

#include "AssignmentGUI.hh"
#include "BezierEditor.hh"
#include "PolyinterpEditor.hh"
#include "PolyapproxEditor.hh"
#include "SplineEditor.hh"
#include "SplineInterpEditor.hh"
#include "SurfaceEditor.hh"
#include "TriSurfEditor.hh"

std::string str_tolower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), static_cast<int (*)(int)>(::tolower));
    return s;
}

int main(int argc, char *argv[]) {
    std::string editorName = "surface"; // default mode
    if (argc == 2) { editorName = str_tolower(argv[1]); }

    // Use a shared DataSource and PolyInterp object for approximation and
    // interpolation editors.
    DataSource dataSource;
    PolyInterp polyInterp;

    AssignmentGUI gui("ECS178 Homework 5",
            [&](AssignmentGUI &gui) {
                gui.editors.emplace("bezier",        std::make_unique<      BezierEditor>(&gui));
                gui.editors.emplace("polyinterp",    std::make_unique<  PolyinterpEditor>(&gui, dataSource, polyInterp));
                gui.editors.emplace("polyapprox",    std::make_unique<  PolyapproxEditor>(&gui, dataSource, polyInterp));
                gui.editors.emplace("spline",        std::make_unique<      SplineEditor>(&gui));
                gui.editors.emplace("spline_interp", std::make_unique<SplineInterpEditor>(&gui, dataSource));
                gui.editors.emplace("surface",       std::make_unique<     SurfaceEditor>(&gui));
                gui.editors.emplace("trisurf",       std::make_unique<     TriSurfEditor>(&gui));

                gui.editorGuiNames = {{"bezier",        "Bézier"},
                                      {"polyinterp",    "Polynomial Interp"},
                                      {"polyapprox",    "Polynomial Approx"},
                                      {"spline",        "Spline"},
                                      {"spline_interp", "Spline Interpolation"},
                                      {"surface",       "NURBS Surface"},
                                      {"trisurf",       "Triangle Mesh"}};
            },
            editorName);
    gui.run();

    return 0;
}
