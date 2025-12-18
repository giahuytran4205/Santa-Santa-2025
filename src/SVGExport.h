#ifndef SVG_EXPORT_H
#define SVG_EXPORT_H
#include <string>
#include "Sparrow.h"

void exportToSVG(const std::string& filename, const SparrowSolver& solver, double scale = 50.0);

#endif