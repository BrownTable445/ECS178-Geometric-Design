#define IMGUI_DEFINE_MATH_OPERATORS
#include <imgui_internal.h>
#include <imgui.h>
#include "tensor_product_surface.hh"
#include "UVDomainWidget.hh"

void UVDomainWidget::run(Surf &surf, float &eval_u, float &eval_v, bool &changedKnots, bool &changedEvalPt) {
    ImGui::Text("UV Domain");
    ImVec2 canvas_pos = ImGui::GetCursorScreenPos();
    ImVec2 canvas_size = ImGui::GetContentRegionAvail();

    // Leave padding for the vertical scrollbar even when it is not shown.
    // This prevents a "feedback loop" where the emergence of a scrollbar
    // shrinks the widget vertically, which triggers removal of the the scrollbar,
    // expanding content vertically and requiring a scrollbar...
    if (ImGui::GetCurrentWindow()->ScrollbarSizes.x == 0) {
        canvas_size.x -= ImGui::GetStyle().ScrollbarSize;
    }

    constexpr float knotRadius   = 4.0f;
    constexpr float evalPtRadius = 5.0f;
    constexpr float labelMargin  = 4.0f;
    constexpr float  knotMargin  = 2.0f;
    constexpr float knotLineThickness = 2.0f;

    ImVec2 multiplicityLabelSize = ImGui::CalcTextSize("5");


    ImU32 const frameCol    = ImGui::GetColorU32(ImGuiCol_FrameBg);
    ImS32 const colorBlue	= ImGui::GetColorU32(ImVec4( 70.0f/255.0f, 102.0f/255.0f, 230.0f/255.0f, 1.0f));
    ImS32 const colorOrange	= ImGui::GetColorU32(ImVec4(255.0f/255.0f, 128.0f/255.0f,  64.0f/255.0f, 1.0f));

    const ImVec2 margins(multiplicityLabelSize.x + labelMargin + knotMargin + 2 * knotRadius,
                         multiplicityLabelSize.y + labelMargin + knotMargin + 2 * knotRadius);

    ImDrawList &drawlist = *ImGui::GetWindowDrawList();

    // Get bounding box of the UV domain.
    auto &u_knots = surf.u_spline.knots;
    auto &v_knots = surf.v_spline.knots;
    ImRect uvDomain(u_knots.minCoeff(), v_knots.minCoeff(),
                    u_knots.maxCoeff(), v_knots.maxCoeff());

    // Determine the vertical canvas size from the horizontal size, margins, and aspect ratio of the domain rect.
    ImVec2 uv_vis_size(canvas_size.x - margins.x, 0);
    uv_vis_size.y = uv_vis_size.x * uvDomain.GetHeight() / uvDomain.GetWidth();
    canvas_size.y = uv_vis_size.y + margins.y;
    ImGui::InvisibleButton("canvas", canvas_size);

    // Determine region on screen containing the full padded uv domain
    ImRect uvDomainFrame(canvas_pos.x + margins.x, canvas_pos.y,
                         canvas_pos.x + canvas_size.x,
                         canvas_pos.y + canvas_size.y - margins.y);

    auto screen_coord_for_u_param = [&](float u) {
        float alpha = (u - uvDomain.Min.x) / uvDomain.GetWidth();
        return (1 - alpha) * uvDomainFrame.Min.x + alpha * uvDomainFrame.Max.x;
    };

    auto screen_coord_for_v_param = [&](float v) {
        float alpha = (v - uvDomain.Min.y) / uvDomain.GetHeight();
        // Note: v values should increase from the bottom (`Max.y`) to the top (`Min.y`)!
        return (1 - alpha) * uvDomainFrame.Max.y + alpha * uvDomainFrame.Min.y;
    };

    auto screen_pt_for_uv = [&](float u, float v) {
        return ImVec2(screen_coord_for_u_param(u), screen_coord_for_v_param(v));
    };

    // Draw the curve's domain of definition
    ImGui::RenderFrame(screen_pt_for_uv(surf.u_spline.domainStart(), surf.v_spline.domainStart()),
                       screen_pt_for_uv(surf.u_spline.domainEnd()  , surf.v_spline.domainEnd()),
                       ImGui::GetColorU32(ImGuiCol_FrameBg), false, 0.0f);

    auto u_marker_pt_center = [&](float u) {
        return ImVec2(screen_coord_for_u_param(u), canvas_pos.y + uv_vis_size.y + knotRadius + knotMargin);
    };
    auto v_marker_pt_center = [&](float v) {
        return ImVec2(canvas_pos.x + multiplicityLabelSize.x + labelMargin + knotRadius, screen_coord_for_v_param(v));
    };

    auto sqNorm = [](ImVec2 v) { return v.x * v.x + v.y * v.y; };

    // Draw faint vertical/horizontal lines at the knots
    for (int i = 0; i < u_knots.size(); ++i) {
        drawlist.AddLine(screen_pt_for_uv(u_knots[i], uvDomain.Min.y),
                         screen_pt_for_uv(u_knots[i], uvDomain.Max.y), ImGui::GetColorU32(ImGuiCol_Border), knotLineThickness);
    }
    for (int j = 0; j < v_knots.size(); ++j) {
        drawlist.AddLine(screen_pt_for_uv(uvDomain.Min.x, v_knots[j]),
                         screen_pt_for_uv(uvDomain.Max.x, v_knots[j]), ImGui::GetColorU32(ImGuiCol_Border), knotLineThickness);
    }

    // Draw u grid lines and knots as little circles with multiplicity text below them.
    for (int i = 0; i < u_knots.size(); /* incremented inside */) {
        float u = u_knots[i];
        ImVec2 p = u_marker_pt_center(u);
        drawlist.AddCircleFilled(p, knotRadius, colorOrange, 4);

        int m = 1;
        while ((i + m < u_knots.size()) && (u_knots(i + m) == u)) ++m;
        if (m > 1) {
            char buffer[16];
            ImFormatString(buffer, IM_ARRAYSIZE(buffer), "%d", m);
            drawlist.AddText(p + ImVec2(-multiplicityLabelSize.x / 2 + 1, 3.0f), colorOrange, buffer);
        }
        i += m;
    }

    // Draw v grid lines and knots as little circles with multiplicity text to the left of them.
    for (int j = 0; j < v_knots.size(); /* incremented inside */) {
        float v = v_knots[j];
        ImVec2 p = v_marker_pt_center(v);
        drawlist.AddCircleFilled(p, knotRadius, colorOrange, 4);

        int m = 1;
        while ((j + m < v_knots.size()) && (v_knots(j + m) == v)) ++m;
        if (m > 1) {
            char buffer[16];
            ImFormatString(buffer, IM_ARRAYSIZE(buffer), "%d", m);
            drawlist.AddText(p + ImVec2(-multiplicityLabelSize.x - knotRadius, -multiplicityLabelSize.y / 2), colorOrange, buffer);
        }
        j += m;
    }

    bool adding_preview = false;

    if (ImGui::IsMouseReleased(ImGuiMouseButton_Left))
        state = State::DRAGGING_NONE;

    const ImVec2 mp = ImGui::GetIO().MousePos;
    if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
        if (uvDomainFrame.Contains(mp))
            state = State::DRAGGING_EVALPT;
        else {
            // Check if we're dragging a knot
            // Find the closest u or v knot, excluding the first and last
            float sqdist = knotRadius * knotRadius;
            int closest_u = -1, closest_v = -1;
            for (size_t i = 1; i < u_knots.rows() - 1; ++i) {
                float d = sqNorm(u_marker_pt_center(u_knots[i]) - mp);
                if (d < sqdist) {
                    sqdist = d;
                    closest_u = i;
                }
            }
            for (size_t j = 1; j < v_knots.rows() - 1; ++j) {
                float d = sqNorm(v_marker_pt_center(v_knots[j]) - mp);
                if (d < sqdist) {
                    sqdist = d;
                    closest_v = j;
                }
            }

            if      (closest_v > -1) { state = State::DRAGGING_VKNOT; drag_idx = closest_v;  }
            else if (closest_u > -1) { state = State::DRAGGING_UKNOT; drag_idx = closest_u;  }
        }
    }
    if (state == State::DRAGGING_EVALPT) {
        float u_frac = (mp.x - uvDomainFrame.Min.x) / uvDomainFrame.GetWidth();
        float v_frac = (uvDomainFrame.Max.y - mp.y) / uvDomainFrame.GetHeight();
        eval_u = (1 - u_frac) * uvDomain.Min.x + u_frac * uvDomain.Max.x;
        eval_v = (1 - v_frac) * uvDomain.Min.y + v_frac * uvDomain.Max.y;
    }

    if (state == State::DRAGGING_UKNOT) {
        float alpha = clamp((mp.x - uvDomainFrame.Min.x) / uvDomainFrame.GetWidth(), 0.0f, 1.0f);
        float u = u_knots[drag_idx] = (1 - alpha) * uvDomain.Min.x + alpha * uvDomain.Max.x;
        std::sort(u_knots.data(), u_knots.data() + u_knots.size());
        drag_idx = 1 + std::distance(u_knots.data() + 1, std::find(u_knots.data(), u_knots.data() + u_knots.size(), u));
        if (drag_idx > u_knots.size()) throw std::runtime_error("Lost track of u knot while dragging");
    }

    if (state == State::DRAGGING_VKNOT) {
        float alpha = clamp((uvDomainFrame.Max.y - mp.y) / uvDomainFrame.GetHeight(), 0.0f, 1.0f);
        float v = v_knots[drag_idx] = (1 - alpha) * uvDomain.Min.y + alpha * uvDomain.Max.y;
        std::sort(v_knots.data(), v_knots.data() + v_knots.size());
        drag_idx = 1 + std::distance(v_knots.data() + 1, std::find(v_knots.data(), v_knots.data() + v_knots.size(), v));
        if (drag_idx > v_knots.size()) throw std::runtime_error("Lost track of v knot while dragging");
    }

    // Make sure the evaluation point is maintained within bounds (also when knots are dragged)
    eval_u = clamp(eval_u, surf.u_spline.domainStart(), surf.u_spline.domainEnd());
    eval_v = clamp(eval_v, surf.v_spline.domainStart(), surf.v_spline.domainEnd());

    changedKnots  = (state == State::DRAGGING_UKNOT) || (state == State::DRAGGING_VKNOT);
    changedEvalPt = state == State::DRAGGING_EVALPT;

    // Draw evaluation point
    drawlist.AddCircleFilled(screen_pt_for_uv(eval_u, eval_v), knotRadius, colorBlue, 16);
}
