#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <cmath>
#include <chrono>
#include <vector>
#include <iostream>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

using std::vector;
using std::cerr;

using boost::geometry::model::point;
using boost::geometry::model::box;
using boost::geometry::cs::cartesian;
using boost::geometry::index::rtree;
using boost::geometry::index::quadratic;
using boost::geometry::index::intersects;


void test_pack_method(){
    typedef point<float, 2, cartesian> point;

    vector<point> points;

    auto t1 = std::chrono::high_resolution_clock::now();
    for (size_t i=0; i<400; i++){
        for (size_t j=0; j<40; j++) {
            points.emplace_back(point {i+0.5f,j+0.5f});
        }
    }

    // Construct using "pack" method (bulk loading)
    rtree<point, quadratic<30> > tree(points);

    auto t2 = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    cerr << "Completed in: " << elapsed.count() << " milliseconds" << '\n';

    box q(point{2,2},point{4,4});

    std::vector<point> result;
    tree.query(intersects(q), std::back_inserter(result));

    cerr << result.size() << '\n';
}


void test_dynamic_method(){
    typedef point<float, 2, cartesian> point;

    vector<point> points;

    // Construct using "pack" method (bulk loading)
    rtree<point, quadratic<30> > tree;

    auto t1 = std::chrono::high_resolution_clock::now();
    for (size_t i=0; i<400; i++){
        for (size_t j=0; j<400; j++) {
            tree.insert(point {i+0.5f,j+0.5f});
        }
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    cerr << "Completed in: " << elapsed.count() << " milliseconds" << '\n';

    box q(point{2,2},point{4,4});

    std::vector<point> result;
    tree.query(intersects(q), std::back_inserter(result));

    cerr << result.size() << '\n';
}

int main()
{
    test_pack_method();
    test_dynamic_method();
    return 0;
}
