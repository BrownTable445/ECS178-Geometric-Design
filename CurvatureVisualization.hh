////////////////////////////////////////////////////////////////////////////////
// CurvatureVisualization.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  GUI support for controlling the visualization of curvature quantities.
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
//  Company:  University of California, Davis
//  Created:  03/13/2023 15:29:47
*///////////////////////////////////////////////////////////////////////////////
#ifndef CURVATUREVISUALIZATION_HH
#define CURVATUREVISUALIZATION_HH

#include<igl/colormap.h>
#include<Eigen/Dense>

struct CurvatureVis {
    enum class Type { None, Gauss, Mean, Kappa1, Kappa2 };
    Type type = Type::None;
    bool autoscale = true;
    float colormapMin = -1, colormapMax = 1;


    Eigen::MatrixX3d colorsForCurvature(const Eigen::VectorXd &curvature) const {
        Eigen::MatrixX3d colors(curvature.rows(), 3);
        if (autoscale)
            igl::colormap(igl::ColorMapType::COLOR_MAP_TYPE_VIRIDIS, curvature, true, colors);
        else
            igl::colormap(igl::ColorMapType::COLOR_MAP_TYPE_VIRIDIS, curvature, colormapMin, colormapMax, colors);
        return colors;
    }

    bool widget(const Eigen::VectorXd &curvatures) {
        int curvatureVisInt = int(type);
        bool changed = false;
        if (ImGui::Combo("Curvature Visualization", &curvatureVisInt, "None\0Gaussian\0Mean\0Kappa 1\0Kappa 2\0")) {
            type = Type(curvatureVisInt);
            changed = true;
        }
        changed |= ImGui::Checkbox("Autoscale Colormap", &autoscale);
        if (curvatures.rows() == 0) return changed;

        if (autoscale) {
            colormapMin = curvatures.minCoeff();
            colormapMax = curvatures.maxCoeff();
            ImGui::Text("Colormap Min: %f", colormapMin);
            ImGui::Text("Colormap Max: %f", colormapMax);
        }
        else {
            changed |=ImGui::InputFloat("Colormap Min", &colormapMin);
            changed |=ImGui::InputFloat("Colormap Max", &colormapMax);
        }
        return changed;
    }
};


#endif /* end of include guard: CURVATUREVISUALIZATION_HH */
