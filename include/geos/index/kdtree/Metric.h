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
#include <geos/geom/Coordinate.h>
#include <geos/util/IllegalArgumentException.h>

using geos::geom::CoordinateXY;

namespace geos {
namespace index { // geos::index
namespace kdtree { // geos::index::kdtree


class GEOS_DLL Metric {
public:
    enum class Type {
        L2,
        L2Squared
    };

    Metric(const Metric&) = delete;
    Metric& operator=(const Metric&) = delete;

    virtual ~Metric() = default;

    virtual inline double distance(const CoordinateXY& c0, const CoordinateXY& c1) const = 0;

    /**
     * Converts a distance from a Metric to the current Metric.
     */
    virtual inline double convertFrom(double distance, const Metric& from) const = 0;

    /**
     * Converts a distance to a Metric from the current Metric.
     */
    virtual inline double convertTo(double distance, const Metric& to) const = 0;

    virtual Type getType() const { return type; }

private:
    Type type;

protected:
    explicit Metric(const Type& t) : type(t) {}
};


class GEOS_DLL L2Squared : public Metric {
public:

    inline double distance(const CoordinateXY& c0, const CoordinateXY& c1) const override {
        double dx = c0.x - c1.x;
        double dy = c0.y - c1.y;
        return dx * dx + dy * dy;
    }

    inline double convertFrom(double distance, const Metric& from) const override {
        if (from.getType() == Type::L2) {
            return distance * distance;
        } else if (from.getType() == Type::L2Squared) {
            return distance;
        } else {
            throw util::IllegalArgumentException("Unknown Metric type");
        }
    }

    inline double convertTo(double distance, const Metric& to) const override {
        if (to.getType() == Type::L2) {
            return std::sqrt(distance);
        } else if (to.getType() == Type::L2Squared) {
            return distance;
        } else {
            throw util::IllegalArgumentException("Unknown Metric type");
        }
    }

    static Metric& getMetric() {
        static L2Squared metric;
        return metric;
    }

private:
    L2Squared() : Metric(Type::L2Squared) {}
};


class GEOS_DLL L2 : public Metric {
public:
    inline double distance(const CoordinateXY& c0, const CoordinateXY& c1) const override {
        double dx = c0.x - c1.x;
        double dy = c0.y - c1.y;
        return std::sqrt(dx * dx + dy * dy);
    }

    inline double convertFrom(double distance, const Metric& from) const override {
        if (from.getType() == Type::L2) {
            return distance;
        } else if (from.getType() == Type::L2Squared) {
            return std::sqrt(distance);
        } else {
            throw util::IllegalArgumentException("Unknown Metric type");
        }
    }

    inline double convertTo(double distance, const Metric& to) const override {
        if (to.getType() == Type::L2) {
            return distance;
        } else if (to.getType() == Type::L2Squared) {
            return distance * distance;
        } else {
            throw util::IllegalArgumentException("Unknown Metric type");
        }
    }

    static Metric& getMetric() {
        static L2 metric;
        return metric;
    }

private:
    L2() : Metric(Type::L2) {}
};


} // namespace geos::index::kdtree
} // namespace geos::index
} // namespace geos
