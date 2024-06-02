/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2024 Edoardo Fusa <edoardo@fusa.science>
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#pragma once

#include <geos/export.h>

#include <geos/index/kdtree/Metric.h>

namespace geos {
namespace index { // geos::index
namespace kdtree { // geos::index::kdtree


/**
 * A class that represents a measure of a generic quantity. Its value is associated with a {@link Metric}.
 */
class GEOS_DLL Measure {
public:
    double value;
    const Metric& metric;

    Measure(double value_in, const Metric& metric_in) : value(value_in), metric(metric_in) {}
};


} // namespace geos::index::kdtree
} // namespace geos::index
} // namespace geos