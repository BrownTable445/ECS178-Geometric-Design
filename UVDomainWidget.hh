#ifndef UVDOMAINWIDGET_HH
#define UVDOMAINWIDGET_HH

#include "tensor_product_surface.hh"

struct UVDomainWidget {
    // Interaction state
    enum class State {
        DRAGGING_NONE, DRAGGING_EVALPT, DRAGGING_UKNOT, DRAGGING_VKNOT
    };

    State state = State::DRAGGING_NONE;
    int drag_idx = -1;

    using Surf = TensorProductSurface_T<3>;
    void run(Surf &surf, float &eval_u, float &eval_v, bool &changedKnots, bool &changedEvalPt);

    void reset() { state = State::DRAGGING_NONE; drag_idx = -1; }
};

#endif /* end of include guard: UVDOMAINWIDGET_HH */
