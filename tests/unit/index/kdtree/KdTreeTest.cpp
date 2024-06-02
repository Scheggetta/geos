#include <tut/tut.hpp>
// geos
#include <geos/index/kdtree/KdTree.h>
#include <geos/index/kdtree/KdTreeMemoryEager.h>
#include <geos/geom/Envelope.h>
#include <geos/io/WKTReader.h>

using namespace geos::index::kdtree;
using namespace geos::geom;

namespace tut {

// dummy data, not used
struct test_kdtree_data {

    geos::io::WKTReader reader_;
    GeometryFactory::Ptr factory = GeometryFactory::create();

    template<typename T>
    std::unique_ptr<std::vector<KdNodeMemoryEager<T>*>> queryResultsToNodeVector(
            const typename KdNodeVisitorMemoryEager<T>::QueryResults& queryResults) {

        std::unique_ptr<std::vector<KdNodeMemoryEager<T>*>> result(new std::vector<KdNodeMemoryEager<T>*>);
        for (const auto& qr : queryResults) {
            result->push_back(qr.first);
        }
        return result;
    }

    template <typename T>
    void testQueryGeom_(std::string& wktInput, Measure tolerance, const CoordinateXY& queryPt,
                        const PreparedGeometry& prepGeom, std::string& wktExpected, bool includeRepeated) {

        KdTreeMemoryEager<T> index(tolerance);

        // Read input and fill tree with it
        auto geo = reader_.read(wktInput);
        std::unique_ptr<CoordinateSequence> coords = geo->getCoordinates();
        for (std::size_t i = 0; i < coords->size(); i++) {
            index.insert(coords->getAt<CoordinateXY>(i));
        }

        // Read expected output into vector of coordinates
        auto geoExpected = reader_.read(wktExpected);
        auto csExpected = geoExpected->getCoordinates();
        std::vector<CoordinateXY> expectedCoord;
        csExpected->toVector(expectedCoord);

        // Read tree into vector of coordinates
        AccumulatingVisitor<T> visitor(queryPt, L2Squared::getMetric());

        typename KdNodeVisitorMemoryEager<T>::QueryResults results;

        results = index.query(prepGeom, visitor);

        std::unique_ptr<std::vector<KdNodeMemoryEager<T>*>> result_ = queryResultsToNodeVector<T>(results);

        std::unique_ptr<std::vector<CoordinateXY>> result = KdTreeMemoryEager<T>::toCoordinates(*result_, includeRepeated);

        std::sort(result->begin(), result->end());
        std::sort(expectedCoord.begin(), expectedCoord.end());

        ensure(result->size() == expectedCoord.size());
        ensure(*result == expectedCoord);
    }

    template<typename T>
    void testQueryGeom(std::string& wktInput, Measure tolerance, const CoordinateXY& queryPt,
                       const PreparedGeometry& prepGeom, std::string& wktExpected) {
        testQueryGeom_<T>(wktInput, tolerance, queryPt, prepGeom, wktExpected, false);
    }

    template<typename T>
    void testQueryGeomRepeated(std::string& wktInput, Measure tolerance, const CoordinateXY& queryPt,
                               const PreparedGeometry& prepGeom, std::string& wktExpected) {
        testQueryGeom_<T>(wktInput, tolerance, queryPt, prepGeom, wktExpected, true);
    }

    void testQuery(std::string& wktInput, double tolerance, const Envelope& queryEnv, std::string& wktExpected, bool includeRepeated) {
        KdTree index(tolerance);
        // Read input and fill tree with it
        auto geo = reader_.read(wktInput);
        std::unique_ptr<CoordinateSequence> coords = geo->getCoordinates();
        for (std::size_t i = 0; i < coords->size(); i++) {
            index.insert(coords->getAt(i));
        }
        // Read expected output into vector of coordinates
        auto geoExpected = reader_.read(wktExpected);
        auto csExpected = geoExpected->getCoordinates();
        std::vector<Coordinate> expectedCoord;
        csExpected->toVector(expectedCoord);
        // Read tree into vector of coordinates
        std::unique_ptr<std::vector<Coordinate>> result = KdTree::toCoordinates(*(index.query(queryEnv)), includeRepeated);

        std::sort(result->begin(), result->end());
        std::sort(expectedCoord.begin(), expectedCoord.end());

        ensure("Result count not equal to expected count", result->size() == expectedCoord.size());
        ensure("Expected result coordinates not found", *result == expectedCoord);
    }

    void testQuery(std::string& wktInput, double tolerance, const Envelope& queryEnv, std::string& wktExpected) {
        testQuery(wktInput, tolerance, queryEnv, wktExpected, false);
    }

    void testQueryRepeated(std::string& wktInput, double tolerance, const Envelope& queryEnv, std::string& wktExpected) {
        testQuery(wktInput, tolerance, queryEnv, wktExpected, true);
    }

};

using group = test_group<test_kdtree_data>;
using object = group::object;

group test_kdtree_group("geos::index::kdtree::KdTree");

//
// testSinglePoint
//
template<>
template<>
void object::test<1> ()
{
    KdTree index(.001);

    KdNode* node1 = index.insert(Coordinate(1, 1));
    KdNode* node2 = index.insert(Coordinate(1, 1));

    ensure("Inserting 2 identical points should create one node", node1 == node2);

    Envelope queryEnv(0, 10, 0, 10);
    std::unique_ptr<std::vector<KdNode*>> result = index.query(queryEnv);

    ensure("query should return 1 result", result->size() == 1);

    KdNode* node = result->at(0);
    ensure("node should have two entries", node->getCount() == 2);
    ensure("node should be repeated", node->isRepeated());
}

//
// testMultiplePoint
//
template<>
template<>
void object::test<2> ()
{
    std::string wkt_in = "MULTIPOINT ((1 1), (2 2))";
    double tolerance = 0.0;
    Envelope env(0, 10, 0, 10);
    std::string wkt_out = "MULTIPOINT ((1 1), (2 2))";
    testQuery(wkt_in, tolerance, env, wkt_out);

    CoordinateXY centroid;
    Geometry::Ptr geom = factory->toGeometry(&env);
    geom->getCentroid(centroid);
    auto prepGeom = PreparedGeometryFactory::prepare(geom.get());
    if (prepGeom == nullptr) { throw std::runtime_error("prepGeom is nullptr"); }
    Measure toleranceMeasure(tolerance * tolerance, L2Squared::getMetric());
    testQueryGeom<int>(wkt_in, toleranceMeasure, centroid, *prepGeom, wkt_out);
}

//
// testSubset
//
template<>
template<>
void object::test<3> ()
{
    std::string wkt_in = "MULTIPOINT ( (1 1), (2 2), (3 3), (4 4) )";
    double tolerance = 0.0;
    Envelope env(1.5, 3.4, 1.5, 3.5);
    std::string wkt_out = "MULTIPOINT ( (2 2), (3 3) )";
    testQuery(wkt_in, tolerance, env, wkt_out);

    CoordinateXY centroid;
    Geometry::Ptr geom = factory->toGeometry(&env);
    geom->getCentroid(centroid);
    auto prepGeom = PreparedGeometryFactory::prepare(geom.get());
    if (prepGeom == nullptr) { throw std::runtime_error("prepGeom is nullptr"); }
    Measure toleranceMeasure(tolerance * tolerance, L2Squared::getMetric());
    testQueryGeom<int>(wkt_in, toleranceMeasure, centroid, *prepGeom, wkt_out);
}

//
// testToleranceFailure
//
template<>
template<>
void object::test<4> ()
{
    std::string wkt_in = "MULTIPOINT ( (0 0), (-.1 1), (.1 1) )";
    double tolerance = 1.0;
    Envelope env(-9, 9, -9, 9);
    std::string wkt_out = "MULTIPOINT ( (0 0), (-.1 1) )";
    testQuery(wkt_in, tolerance, env, wkt_out);

    CoordinateXY centroid;
    Geometry::Ptr geom = factory->toGeometry(&env);
    geom->getCentroid(centroid);
    auto prepGeom = PreparedGeometryFactory::prepare(geom.get());
    if (prepGeom == nullptr) { throw std::runtime_error("prepGeom is nullptr"); }
    Measure toleranceMeasure(tolerance * tolerance, L2Squared::getMetric());
    testQueryGeom<int>(wkt_in, toleranceMeasure, centroid, *prepGeom, wkt_out);
}

//
// testTolerance2
//
template<>
template<>
void object::test<5> ()
{
    std::string wkt_in = "MULTIPOINT ((10 60), (20 60), (30 60), (30 63))";
    double tolerance = 9.0;
    Envelope env(0,99, 0, 99);
    std::string wkt_out = "MULTIPOINT ((10 60), (20 60), (30 60))";
    testQuery(wkt_in, tolerance, env, wkt_out);

    CoordinateXY centroid;
    Geometry::Ptr geom = factory->toGeometry(&env);
    geom->getCentroid(centroid);
    auto prepGeom = PreparedGeometryFactory::prepare(geom.get());
    if (prepGeom == nullptr) { throw std::runtime_error("prepGeom is nullptr"); }
    Measure toleranceMeasure(tolerance * tolerance, L2Squared::getMetric());
    testQueryGeom<int>(wkt_in, toleranceMeasure, centroid, *prepGeom, wkt_out);
}

//
// testTolerance2_perturbedY
//
template<>
template<>
void object::test<6> ()
{
    std::string wkt_in = "MULTIPOINT ((10 60), (20 61), (30 60), (30 63))";
    double tolerance = 9.0;
    Envelope env(0,99, 0, 99);
    std::string wkt_out = "MULTIPOINT ((10 60), (20 61), (30 60))";
    testQuery(wkt_in, tolerance, env, wkt_out);

    CoordinateXY centroid;
    Geometry::Ptr geom = factory->toGeometry(&env);
    geom->getCentroid(centroid);
    auto prepGeom = PreparedGeometryFactory::prepare(geom.get());
    if (prepGeom == nullptr) { throw std::runtime_error("prepGeom is nullptr"); }
    Measure toleranceMeasure(tolerance * tolerance, L2Squared::getMetric());
    testQueryGeom<int>(wkt_in, toleranceMeasure, centroid, *prepGeom, wkt_out);
}

//
// testSnapToNearest
//
template<>
template<>
void object::test<7> ()
{
    std::string wkt_in = "MULTIPOINT ( (10 60), (20 60), (16 60))";
    double tolerance = 5.0;
    Envelope env(0,99, 0, 99);
    std::string wkt_out = "MULTIPOINT ( (10 60), (20 60), (20 60))";
    testQueryRepeated(wkt_in, tolerance, env, wkt_out);

    CoordinateXY centroid;
    Geometry::Ptr geom = factory->toGeometry(&env);
    geom->getCentroid(centroid);
    auto prepGeom = PreparedGeometryFactory::prepare(geom.get());
    if (prepGeom == nullptr) { throw std::runtime_error("prepGeom is nullptr"); }
    Measure toleranceMeasure(tolerance * tolerance, L2Squared::getMetric());
    testQueryGeomRepeated<int>(wkt_in, toleranceMeasure, centroid, *prepGeom, wkt_out);
}



//
// testSinglePointWithPayload
//
template<>
template<>
void object::test<8> ()
{
    Measure tolerance(.001, L2::getMetric());
    KdTreeMemoryEager<int> index(tolerance);
    KdNodeMemoryEager<int>* node1 = index.insert(Coordinate(1, 1), std::make_unique<int>(42));
    KdNodeMemoryEager<int>* node2 = index.insert(Coordinate(1, 1), std::make_unique<int>(43));

    ensure(node1 == node2);

    Envelope queryEnv(0, 10, 0, 10);
    std::unique_ptr<std::vector<KdNodeMemoryEager<int>*>> result = index.query(queryEnv);

    ensure(result->size() == 1);

    auto* node = (KdNodeMemoryEager<int>*)(*result)[0];
    ensure(node->getCount() == 2);
    ensure(node->isRepeated());
    ensure(*(node->getData()) == 42);
}


//
// testPayloadWithEnvelope
//
template<>
template<>
void object::test<9> ()
{
    Measure tolerance(.001, L2::getMetric());
    KdTreeMemoryEager<int> index(tolerance);
    KdNodeMemoryEager<int>* node1 = index.insert(Coordinate(1, 1), std::make_unique<int>(-12));
    KdNodeMemoryEager<int>* node2 = index.insert(Coordinate(1, 1), std::make_unique<int>(100));

    ensure(node1 == node2);

    Envelope queryEnv(0, 10, 0, 10);
    std::unique_ptr<std::vector<KdNodeMemoryEager<int>*>> result = index.query(queryEnv);

    ensure(result->size() == 1);

    KdNodeMemoryEager<int>* node = result->at(0);

    ensure(node->getCount() == 2);
    ensure(node->isRepeated());
    ensure(*(node->getData()) == -12);
}


//
// testPayloadWithPreparedGeometry
//
template<>
template<>
void object::test<10> ()
{
    Measure tolerance(.001, L2::getMetric());
    KdTreeMemoryEager<int> index(tolerance);
    KdNodeMemoryEager<int>* node1 = index.insert(Coordinate(1, 1), std::make_unique<int>(-12));
    KdNodeMemoryEager<int>* node2 = index.insert(Coordinate(1, 1), std::make_unique<int>(100));

    ensure(node1 == node2);

    Envelope queryEnv(0, 10, 0, 10);
    Geometry::Ptr geom = factory->toGeometry(&queryEnv);
    PreparedPolygon prepPol(geom.get());

    ClosestNodeVisitor<int> visitor(CoordinateXY(0.5, 0.5), L2Squared::getMetric());
    const KdNodeVisitorMemoryEager<int>::QueryResults& results = index.query(prepPol, visitor);

    ensure(results.size() == 1);

    KdNodeMemoryEager<int>* node = visitor.getNodeAt(0);

    ensure(node->getCount() == 2);
    ensure(node->isRepeated());
    ensure(*(node->getData()) == -12);
}


//
// testRadiusSearch
//
template<>
template<>
void object::test<11> ()
{
    KdTreeMemoryEager<int> index;
    index.insert(CoordinateXY(-2.0, 0.0), std::make_unique<int>(0));
    index.insert(CoordinateXY(0.5, 0.0), std::make_unique<int>(1));
    index.insert(CoordinateXY(0.45, 0.57), std::make_unique<int>(2));
    index.insert(CoordinateXY(-0.499, 0.5), std::make_unique<int>(3));

    CoordinateXY queryPt(0.5, 0.5);
    // even though the radius is 1.0 and has been defined with L2,
    // the index of the tree is L2Squared (default value) so the visitor will use L2Squared
    Measure radiusMeasure(1.0, L2::getMetric());
    AccumulatingVisitor<int> visitor(queryPt, L2Squared::getMetric());
    bool ordered = true;

    const KdNodeVisitorMemoryEager<int>::QueryResults& results = index.radius_search(radiusMeasure, visitor, ordered);

    ensure(results.size() == 3);
    ensure(visitor.getPayloadAt(0) == 2);
    ensure(visitor.getPayloadAt(1) == 1);
    ensure(visitor.getPayloadAt(2) == 3);
}


//
// testRadiusSearchNoPayload
//
template<>
template<>
void object::test<12> ()
{
    KdTreeMemoryEager<NoPayload> index;
    index.insert(CoordinateXY(-2.0, 0.0));
    index.insert(CoordinateXY(0.5, 0.0));
    index.insert(CoordinateXY(0.45, 0.57));
    index.insert(CoordinateXY(-0.499, 0.5));

    CoordinateXY queryPt(0.5, 0.5);
    // even though the radius is 1.0 and has been defined with L2,
    // the index of the tree is L2Squared (default value) so the visitor will use L2Squared
    Measure radiusMeasure(1.0, L2::getMetric());
    AccumulatingVisitor<NoPayload> visitor(queryPt, L2Squared::getMetric());
    bool ordered = true;

    const KdNodeVisitorMemoryEager<NoPayload>::QueryResults& results = index.radius_search(radiusMeasure, visitor, ordered);

    ensure(results.size() == 3);
    ensure(visitor.getNodeAt(0) == index.query(CoordinateXY(0.45, 0.57)));
    ensure(visitor.getNodeAt(1) == index.query(CoordinateXY(0.5, 0.0)));
    ensure(visitor.getNodeAt(2) == index.query(CoordinateXY(-0.499, 0.5)));
}


} // namespace tut
