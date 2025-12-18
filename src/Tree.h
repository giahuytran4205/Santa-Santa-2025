#ifndef TREE_H
#define TREE_H
#include "Geometry.h"

CompositeShape createTree() {
    CompositeShape tree;

    // --- PART 1: TOP TIER ---
    tree.addPart({
        {0.0f, 0.8f},
        {0.125f, 0.5f},
        {0.0625f, 0.5f},
        {-0.0625f, 0.5f},
        {-0.125f, 0.5f}
    });

    // --- PART 2: MIDDLE TIER ---
    tree.addPart({
        {0.0625f, 0.5f},
        {0.2f, 0.25f},
        {0.1f, 0.25f},
        {-0.1f, 0.25f},
        {-0.2f, 0.25f},
        {-0.0625f, 0.5f}
    });

    // --- PART 3: BOTTOM FOLIAGE ---
    tree.addPart({
        {0.1f, 0.25f},
        {0.35f, 0.0f},
        {0.075f, 0.0f},
        {-0.075f, 0.0f},
        {-0.35f, 0.0f},
        {-0.1f, 0.25f}
    });

    // --- PART 4: TRUNK ---
    tree.addPart({
        {0.075f, 0.0f},
        {0.075f, -0.2f},
        {-0.075f, -0.2f},
        {-0.075f, 0.0f}
    });

    return tree;
}

#endif // TREE_H