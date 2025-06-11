#include "KnotViewer.hh"
#include "SplineEditor.hh"

void KnotViewer::update() {
    clear();
    if (disabled()) return;

    const size_t n = numKnots();
    Eigen::MatrixX3f V(2 * (1 + n), 3);
    Eigen::MatrixX2i E(1 + n, 2);

    // Plot the U axis
    V.topRows<2>() << 0, 0, 0,
                      1, 0, 0;
    E.row(0) << 0, 1;

    // Plot a tick mark for each knot, clustering the repeated knots.
    auto alpha_for_u = [&](float u) { return (u - knot(0)) / (knot(n - 1) - knot(0)); };
    size_t multiplicity = 0;
    for (size_t i = 0; i < n; i += multiplicity) {
        float u = knot(i);
        multiplicity = 0;
        while ((i + multiplicity < n) && (knot(i + multiplicity) == u)) ++multiplicity;

        float alpha = alpha_for_u(u);
        for (size_t j = i; j < i + multiplicity; ++j) {
            // Center the repeated knot tickmarks around the knot value.
            float alpha_spaced = alpha + repeatedKnotSpacing * ((j - i) - 0.5 * (multiplicity - 1));
            V.row(2 * (j + 1) + 0) << alpha_spaced,  tickHeight, 0.0f;
            V.row(2 * (j + 1) + 1) << alpha_spaced, -tickHeight, 0.0f;
            E.row(1 + j) << 2 * (j + 1), 2 * (j + 1) + 1;
        }
    }

    float dpiScale = static_cast<Viewer &>(*viewer).getDPIScale();

    // Curve evaluation location marker.
    auto splineEditor = dynamic_cast<SplineEditor *>(m_splineEditor);
    if (splineEditor) {
        const size_t numSegments = 20;
        const float aspectRatio = (2.0f * viewer->core().viewport[2]) / viewer->core().viewport[3]; // account for nonuniform stretch by modelViewProjection() and window;
        std::vector<Eigen::Vector3f> markerPts;
        for (size_t i = 0; i < numSegments; ++i) {
            double theta = (2 * M_PI * i) / numSegments;
            markerPts.emplace_back(alpha_for_u(splineEditor->eval_u) + 0.5 * (tickHeight / aspectRatio) * cos(theta), 0.5 * tickHeight * sin(theta), 0.0f);
        }
        addPolyline(markerPts, true, 4 * dpiScale, BezierEditor::curveColor(true));
    }

    addEdges(V, E, 2.5 * dpiScale);
}

bool KnotViewer::mouse_down(int button, int modifier) {
    // Only allow manipulation if the m_splineEditor is an ordinary SplineEditor (not in interpolation mode)
    if (dynamic_cast<SplineEditor *>(m_splineEditor) == nullptr) return false;
    if (disabled() || (numKnots() == 0)) return false;

    Eigen::Vector3f pt = unprojectPt(viewer->down_mouse_x, viewer->down_mouse_y);

    if (pt[1] < -tickHeight || pt[1] > tickHeight) return false; // outside knot bar

    // Search for the closest knot to the click point, omitting the first and last
    float u = knotCoordinateForPt(pt);
    const size_t n = numKnots();

    // Do a binary search for the closest knot, omitting the end knots (which we
    // don't want to drag).
    size_t i = clamp<size_t>(knotLowerBoundIndex(u), 1, n - 2);

    float dist = std::abs(u - knot(i));
    // Bear in mind the nearest knot could be the one behind `u`...
    if ((i > 1) && std::abs(u - knot(i - 1)) < dist) {
        dist = std::abs(u - knot(i - 1));
        --i;
    }

    if (dist / (knot(n - 1) - knot(0)) < 1e-2)
        drag_idx = i;

    return true;
}

bool KnotViewer::mouse_move(int mouse_x, int mouse_y) {
    if (disabled() || !isDragging()) return false;

    size_t n = numKnots();
    float u = clamp(knotCoordinateForPt(unprojectPt(mouse_x, mouse_y)), knot(0), knot(n - 1));
    knot(drag_idx) = u;

    // Maintain ordering (this may change the index of the knot we're dragging!)
    std::sort(knotData(), knotData() + numKnots());

    // Locate the knot we're dragging (or one of its copies) in the sorted sequence.
    drag_idx = clamp<size_t>(knotLowerBoundIndex(u), 1, n - 2);
    if (knot(drag_idx) != u) throw std::runtime_error("Lost track of knot while dragging");

    // Notify the spline viewer of the curve update...
    m_splineEditor->updateView(*viewer);
    update();

    return true;
}
