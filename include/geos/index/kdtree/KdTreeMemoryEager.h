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
#include <geos/geom/Envelope.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/prep/PreparedGeometryFactory.h>
#include <geos/geom/prep/PreparedPolygon.h>

#include <geos/index/kdtree/KdNodeVisitorMemoryEager.h>
#include <geos/index/kdtree/KdNodeMemoryEager.h>
#include <geos/index/kdtree/Metric.h>
#include <geos/index/kdtree/Measure.h>

#include <memory>
#include <vector>
#include <string>
#include <algorithm>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4251) // warning C4251: needs to have dll-interface to be used by clients of class
#endif

using geos::geom::Geometry;
using geos::geom::Envelope;
using geos::geom::prep::PreparedGeometry;
using geos::geom::prep::PreparedPolygon;
using geos::geom::prep::PreparedGeometryFactory;

namespace geos {
namespace index { // geos::index
namespace kdtree { // geos::index::kdtree


template <typename T>
class GEOS_DLL KdTreeMemoryEager {

private:
    using TPtr = std::unique_ptr<T>;

    std::deque<KdNodeMemoryEager<T>> nodes;  // Deque to ensure pointers are not invalidated
    KdNodeMemoryEager<T>* root;
    std::size_t numberOfNodes;
    double toleranceMetric;
    const Metric& metric;
    const geom::GeometryFactory::Ptr factory = geom::GeometryFactory::create();
    

    KdNodeMemoryEager<T>* findClosestNode(const CoordinateXY& p) {
        ClosestNodeVisitor<T> visitor(p, metric);
        Geometry::Ptr circle = factory->createPoint(p)->
                buffer(metric.convertTo(toleranceMetric, L2::getMetric()));
        PreparedPolygon prepCircle(circle.get());

        query(prepCircle, visitor);
        if (visitor.noResultFound()) {
            return nullptr;
        } else {
            return visitor.getNodeAt(0);
        }
    }


    /**
    * Create a node on a locally managed vector
    */
    inline KdNodeMemoryEager<T>* createNode(const CoordinateXY& p, TPtr payload) {
        nodes.emplace_back(p, std::move(payload));
        return &nodes.back();
    }


    KdNodeMemoryEager<T>* insertExact(const CoordinateXY& p, TPtr payload) {
        KdNodeMemoryEager<T>* currentNode = root;
        KdNodeMemoryEager<T>* leafNode = root;
        bool isOddLevel = true;
        bool isLessThan = true;

        /**
        * Traverse the tree, first cutting the plane left-right (by X ordinate)
        * then top-bottom (by Y ordinate)
        */

        while (currentNode != nullptr) {
            // test if point is already a node (not strictly necessary because already checked)
            // bool isInTolerance = metric.distance(p, currentNode->getCoordinate()) <= toleranceMetric;

            // check if point is already in tree (up to tolerance) and if so simply
            // return existing node
            // if (isInTolerance) {
            //     currentNode->increment();
            //     return currentNode;
            // }

            if (isOddLevel) {
                isLessThan = p.x < currentNode->getX();
            } else {
                isLessThan = p.y < currentNode->getY();
            }
            leafNode = currentNode;
            if (isLessThan) {
                currentNode = currentNode->getLeft();
            } else {
                currentNode = currentNode->getRight();
            }

            isOddLevel = not isOddLevel;
        }

        // no node found, add new leaf node to tree
        numberOfNodes++;
        KdNodeMemoryEager<T>* node = createNode(p, std::move(payload));
        if (isLessThan) {
            leafNode->setLeft(node);
        } else {
            leafNode->setRight(node);
        }
        return node;
    }


    void queryNode(KdNodeMemoryEager<T>* currentNode,
                   const PreparedGeometry& prepGeom,
                   bool odd,
                   KdNodeVisitorMemoryEager<T>& visitor) {
        // NOTE: `visitor` can have a different metric than the tree. If that is the case, it doesn't matter for the
        //       purpose of this function since the two classes are independent.

        // Non-recursive formulation of in-order traversal from
        // http://web.cs.wpi.edu/~cs2005/common/iterative.inorder
        // Otherwise we may blow up the stack
        // See https://github.com/qgis/QGIS/issues/45226

        using Pair = std::pair<KdNodeMemoryEager<T>*, bool>;
        std::stack<Pair> activeNodes;

        const Envelope* env = prepGeom.getGeometry().getEnvelopeInternal();

        while (true) {
            if (currentNode != nullptr) {
                double min;
                double discriminant;

                if (odd) {
                    min = env->getMinX();
                    discriminant = currentNode->getX();
                } else {
                    min = env->getMinY();
                    discriminant = currentNode->getY();
                }
                bool searchLeft = min < discriminant;

                activeNodes.emplace(currentNode, odd);

                // search is computed via in-order traversal
                KdNodeMemoryEager<T>* leftNode = nullptr;
                if (searchLeft) {
                    leftNode = currentNode->getLeft();
                }
                if (leftNode) {
                    currentNode = leftNode;
                    odd = not odd;
                } else {
                    currentNode = nullptr;
                }

            } else if (not activeNodes.empty()) {
                currentNode = activeNodes.top().first;
                odd = activeNodes.top().second;
                activeNodes.pop();

                if (prepGeom.contains(factory->createPoint(currentNode->getCoordinate()).get())) {
                    visitor.visit(currentNode);
                    if (visitor.stopSearch()) {
                        break;
                    }
                }

                double max;
                double discriminant;

                if (odd) {
                    max = env->getMaxX();
                    discriminant = currentNode->getX();
                } else {
                    max = env->getMaxY();
                    discriminant = currentNode->getY();
                }
                bool searchRight = discriminant <= max;

                if (searchRight) {
                    currentNode = currentNode->getRight();
                    if (currentNode) {
                        odd = not odd;
                    }
                } else {
                    currentNode = nullptr;
                }
            } else {
                break;
            }
        }
    }


    KdNodeMemoryEager<T>* queryNodePoint(KdNodeMemoryEager<T>* currentNode,
                              const CoordinateXY& queryPt,
                              bool odd) {
        while (currentNode != nullptr) {
            if (currentNode->getCoordinate().equals2D(queryPt)) {
                return currentNode;
            }

            double ord;
            double discriminant;
            if (odd) {
                ord = queryPt.x;
                discriminant = currentNode->getX();
            }
            else {
                ord = queryPt.y;
                discriminant = currentNode->getY();
            }

            bool searchLeft = (ord < discriminant);
            odd = not odd;
            if (searchLeft) {
                currentNode = currentNode->getLeft();
            }
            else {
                currentNode = currentNode->getRight();
            }
        }
        return nullptr;
    }

public:

    /**
    * Converts a collection of {@link KdNodeMemoryEager}s to an vector of {@link geom::CoordinateXY}s.
    *
    * @param kdnodes a collection of nodes
    * @return an vector of the coordinates represented by the nodes
    */
    inline static std::unique_ptr<std::vector<geom::CoordinateXY>> toCoordinates(std::vector<KdNodeMemoryEager<T>*>& kdnodes) {
        return toCoordinates(kdnodes, false);
    }


    /**
    * Converts a collection of {@link KdNodeMemoryEager}s
    * to an vector of {@link geom::CoordinateXY}s,
    * specifying whether repeated nodes should be represented
    * by multiple coordinates.
    *
    * @param kdnodes a collection of nodes
    * @param includeRepeated true if repeated nodes should
    *   be included multiple times
    * @return an vector of the coordinates represented by the nodes
    */
    static std::unique_ptr<std::vector<geom::CoordinateXY>> toCoordinates(std::vector<KdNodeMemoryEager<T>*>& kdnodes,
                                                                          bool includeRepeated) {

        std::unique_ptr<std::vector<CoordinateXY>> coord(new std::vector<CoordinateXY>);
        for (auto node: kdnodes) {
            std::size_t count = includeRepeated ? node->getCount() : 1;
            for (std::size_t i = 0; i < count; i++) {
                coord->emplace_back(node->getCoordinate());
            }
        }
        if (not includeRepeated) {
            // Remove duplicate Coordinates from coordList
            coord->erase(std::unique(coord->begin(), coord->end()), coord->end());
        }
        return coord;
    }


    explicit KdTreeMemoryEager(const Metric& metric_in = L2Squared::getMetric()) :
        nodes(),
        root(nullptr),
        numberOfNodes(0),
        toleranceMetric(0.0),
        metric(metric_in)
        {};


    explicit KdTreeMemoryEager(const Measure& measure) :
        nodes(),
        root(nullptr),
        numberOfNodes(0),
        toleranceMetric(measure.value),
        metric(measure.metric) {

        if (measure.value < 0.0) {
            throw util::IllegalArgumentException("The value of the Measure obj must be non-negative");
        }
    };


    inline size_t size() const { return numberOfNodes; }

    inline bool isEmpty() const { return root == nullptr; }


    /**
    * Inserts a new point in the kd-tree.
    */
    inline KdNodeMemoryEager<T>* insert(const geom::CoordinateXY& p) {
        return insert(p, nullptr);
    }


    KdNodeMemoryEager<T>* insert(const geom::CoordinateXY& p, TPtr payload) {
        if (root == nullptr) {
            root = createNode(p, std::move(payload));
            numberOfNodes++;
            return root;
        }

        /**
        * Check if the point is already in the tree, up to toleranceMetric.
        * If toleranceMetric is zero, this phase of the insertion can be skipped.
        */
        if (toleranceMetric > 0) {
            KdNodeMemoryEager<T>* closestNode = findClosestNode(p);
            if (closestNode != nullptr) {
                // point already in index - increment counter
                closestNode->increment();
                return closestNode;
            }
        }

        return insertExact(p, std::move(payload));
    }


    /**
    * Performs a search of the points in the index and visits the ones in the given {@link PreparedGeometry}
    */
    const typename KdNodeVisitorMemoryEager<T>::QueryResults& query(const PreparedGeometry& prepGeom,
                                                                    KdNodeVisitorMemoryEager<T>& visitor,
                                                                    bool ordered = true) {
        queryNode(root, prepGeom, true, visitor);
        return visitor.getResults(ordered);
    }

    /**
    * Performs a search of the points in the index and visits the ones in the given {@link Envelope}
    */
    const typename KdNodeVisitorMemoryEager<T>::QueryResults& query(const Envelope& queryEnv,
                                                                    KdNodeVisitorMemoryEager<T>& visitor,
                                                                    bool ordered = true) {
        Geometry::Ptr rect = factory->toGeometry(&queryEnv);
        PreparedPolygon prepRect(rect.get());

        queryNode(root, prepRect, true, visitor);
        return visitor.getResults(ordered);
    }

    /**
     * Performs a search of the points in the index and returns the ones in the given {@link Envelope}.
     * Legacy method, use {@link #query(const Envelope&, KdNodeVisitorMemoryEager)} instead.
     */
    std::unique_ptr<std::vector<KdNodeMemoryEager<T>*>> query(const Envelope& queryEnv, bool ordered = true) {
        // The distances of the results are wrt the centre of the query envelope
        std::unique_ptr<std::vector<KdNodeMemoryEager<T>*>> result(new std::vector<KdNodeMemoryEager<T>*>);

        CoordinateXY centre;
        queryEnv.centre(centre);
        AccumulatingVisitor<T> visitor(centre, metric);

        const auto& queryResults = query(queryEnv, visitor, ordered);  // queryResults can be used as long as visitor is alive
        for (const auto& qr : queryResults) {
            result->push_back(qr.first);
        }
        return result;
    }

    /**
    * Searches for a given point in the index and returns its node if found.
    */
    KdNodeMemoryEager<T>* query(const geom::CoordinateXY& queryPt) {
        return queryNodePoint(root, queryPt, true);
    }

    /**
     * Performs a radius search of the points in the index and visits the ones in the given radius
     */
    const typename KdNodeVisitorMemoryEager<T>::QueryResults& radius_search(const Measure& radiusMeasure,
                                                                            AccumulatingVisitor<T>& visitor,
                                                                            bool ordered = true) {

        double radius = radiusMeasure.metric.convertTo(radiusMeasure.value, L2::getMetric());

        // `radius` must be non-squared and positive and the returned distances will have the tree's metric
        if (radius <= 0) {
            throw util::IllegalArgumentException("Radius must be positive");
        }

        Geometry::Ptr circle = factory->createPoint(visitor.getQueryPoint())->buffer(radius);
        PreparedPolygon prepCircle(circle.get());

        return query(prepCircle, visitor, ordered);
    }

};

} // namespace geos::index::kdtree
} // namespace geos::index
} // namespace geos

#ifdef _MSC_VER
#pragma warning(pop)
#endif
