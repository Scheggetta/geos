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

#include <algorithm>

#include <geos/export.h>

#include <geos/index/kdtree/KdNodeMemoryEager.h>
#include <geos/index/kdtree/Metric.h>

namespace geos {
namespace index { // geos::index
namespace kdtree { // geos::index::kdtree


template <typename T>
class GEOS_DLL KdNodeVisitorMemoryEager {
public:
    // Declare type as non-copyable
    KdNodeVisitorMemoryEager(const KdNodeVisitorMemoryEager&) = delete;
    KdNodeVisitorMemoryEager& operator=(const KdNodeVisitorMemoryEager&) = delete;

    explicit KdNodeVisitorMemoryEager(const CoordinateXY& p_in) : p(p_in), resultList() {};
    KdNodeVisitorMemoryEager() : p(CoordinateXY::getNull()), resultList() {};

    virtual ~KdNodeVisitorMemoryEager() = default;

    using QueryResults = std::vector<std::pair<KdNodeMemoryEager<T>*, double>>;

    virtual void visit(KdNodeMemoryEager<T>* node) = 0;
    virtual inline bool stopSearch() const = 0;

    inline bool hasValidQueryPoint() const { return not p.isNull(); }
    inline const CoordinateXY& getQueryPoint() const { return p; }

    inline bool noResultFound() const { return resultList.empty(); }
    inline std::size_t size() const { return resultList.size(); }

    inline const QueryResults& getResults(bool ordered) {
        if (ordered) {
            std::sort(resultList.begin(), resultList.end(), [](const auto& a, const auto& b) {
                return a.second < b.second;
            });
        }
        return resultList;
    }

    inline const T& getPayloadAt(std::size_t i) const { return *(getNodeAt(i)->getData()); }

    inline void add(KdNodeMemoryEager<T>* node, double distance) { resultList.emplace_back(node, distance); }

    void setAt(std::size_t i, KdNodeMemoryEager<T>* node, double distance) {
        resultList[i].first = node;
        resultList[i].second = distance;
    }

    inline void clear() { resultList.clear(); }

    inline const std::pair<KdNodeMemoryEager<T>*, double>& getAt(std::size_t i) const { return resultList[i]; }

    inline KdNodeMemoryEager<T>* getNodeAt(std::size_t i) const { return resultList[i].first; }

    inline double getDistanceAt(std::size_t i) const { return resultList[i].second; }

private:
    CoordinateXY p;
    QueryResults resultList;  // node, metric distance

protected:
    void initializeResults() { resultList = QueryResults(0); }

    inline void setQueryPoint(const CoordinateXY& point) { p = point; }
};


/**
 * ClosestNodeVisitor retains the closest node found during a search.
 */
template <typename T>
class ClosestNodeVisitor : public KdNodeVisitorMemoryEager<T> {
public:
    explicit ClosestNodeVisitor(const CoordinateXY& p, const Metric& metric_in) :
        KdNodeVisitorMemoryEager<T>(p),
        closestDistance(std::numeric_limits<double>::max()),
        metric(metric_in) {
        this->initializeResults();  // to avoid maybe-uninitialized warning
    }

    inline void visit(KdNodeMemoryEager<T>* node) override {
        double distance = metric.distance(this->getQueryPoint(), node->getCoordinate());
        if (this->noResultFound()) {
            this->add(node, distance);
            closestDistance = distance;
        } else if (distance < closestDistance) {
            this->setAt(0, node, distance);
            closestDistance = distance;
        }
    }

    inline bool stopSearch() const override { return false; }

    inline void clearResults() {
        this->clear();
        closestDistance = std::numeric_limits<double>::max();
    }

    inline void clearResults(const CoordinateXY& newPt) {
        clearResults();
        this->setQueryPoint(newPt);
    }

    inline KdNodeMemoryEager<T>* getNode() const { return this->getNodeAt(0); }

    inline double getDistance() const { return this->getDistanceAt(0); }

private:
    double closestDistance;
    const Metric& metric;
};


/**
 * AccumulatingVisitor stores all nodes visited during a search.
 */
template <typename T>
class AccumulatingVisitor : public KdNodeVisitorMemoryEager<T> {
public:

    explicit AccumulatingVisitor(const CoordinateXY& p, const Metric& metric_in) :
        KdNodeVisitorMemoryEager<T>(p),
        metric(metric_in) {
        this->initializeResults();  // to avoid maybe-uninitialized warning
    }

    inline void visit(KdNodeMemoryEager<T>* node) override {
        double distance = metric.distance(this->getQueryPoint(), node->getCoordinate());
        this->add(node, distance);
    }

    inline bool stopSearch() const override { return false; }

    inline void clearResults() {
        this->clear();
    }

    inline void clearResults(const CoordinateXY& newPt) {
        clearResults();
        this->setQueryPoint(newPt);
    }

private:
    const Metric& metric;
};


} // namespace geos::index::kdtree
} // namespace geos::index
} // namespace geos
