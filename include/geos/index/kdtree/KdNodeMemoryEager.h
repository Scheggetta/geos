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

#include <memory>

#include <geos/export.h>
#include <geos/geom/Coordinate.h> // for composition

using geos::geom::CoordinateXY;

namespace geos {
namespace index { // geos::index
namespace kdtree { // geos::index::kdtree


/**
 * A class that represents no payload for a {@link KdNodeMemoryEager}.
 */
class GEOS_DLL NoPayload {
public:
    // Delete all the constructors
    NoPayload() = delete;
    NoPayload(const NoPayload&) = delete;
    NoPayload(NoPayload&&) = delete;
    NoPayload& operator=(const NoPayload&) = delete;
    NoPayload& operator=(NoPayload&&) = delete;
};


/**
 * A node of a {@link KdTreeMemoryEager}, which represents one or more points in the same location.
 */
template <typename T>  // T is the payload type
class GEOS_DLL KdNodeMemoryEager {
private:
    using TPtr = std::unique_ptr<T>;

    CoordinateXY p;
    TPtr payload;  // The payload is managed uniquely by the node
    KdNodeMemoryEager* left;
    KdNodeMemoryEager* right;
    std::size_t count;

public:

    KdNodeMemoryEager(double x, double y, TPtr payload_in) :
        p(x, y),
        payload(std::move(payload_in)),
        left(nullptr),
        right(nullptr),
        count(1) {}

    KdNodeMemoryEager(const CoordinateXY& p_in, TPtr payload_in) :
        p(p_in),
        payload(std::move(payload_in)),
        left(nullptr),
        right(nullptr),
        count(1) {}

    inline double getX() const { return p.x; }
    inline double getY() const { return p.y; }
    inline const CoordinateXY& getCoordinate() const { return p; }
    inline T* getData() const { return payload.get(); }
    inline KdNodeMemoryEager* getLeft() const { return left; }
    inline KdNodeMemoryEager* getRight() const { return right; }
    inline void increment() { count++; }
    inline std::size_t getCount() const { return count; }
    inline bool isRepeated() const { return count > 1; }
    inline void setLeft(KdNodeMemoryEager* p_left) { left = p_left; }
    inline void setRight(KdNodeMemoryEager* p_right) { right = p_right; }

};


} // namespace geos::index::kdtree
} // namespace geos::index
} // namespace geos
