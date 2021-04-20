#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

using boost::geometry::wkt;
using boost::geometry::model::point;
using boost::geometry::model::box;
using boost::geometry::cs::cartesian;
using boost::geometry::index::rtree;
using boost::geometry::index::quadratic;
using boost::geometry::index::intersects;


#include <cmath>
#include <chrono>
#include <vector>
#include <iostream>
#include <utility>

using std::vector;
using std::cerr;
using std::make_pair;
using std::pair;


void test_pack_method(){
    typedef point<float, 2, cartesian> point;

    vector<point> points;

    auto t1 = std::chrono::high_resolution_clock::now();
    for (size_t i=0; i<800; i++){
        for (size_t j=0; j<80; j++) {
            points.emplace_back(point {i+0.5f,j+0.5f});
        }
    }

    // Construct using "pack" method (bulk loading)
    rtree<point, quadratic<30> > tree(points);

    auto t2 = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    cerr << "Construction completed in: " << elapsed.count() << " milliseconds" << '\n';

    std::vector<point> result;

    for (size_t i=0; i<800; i++){
        for (size_t j=0; j<80; j++) {
            box q(point{i,j},point{i+10,j+10});

            tree.query(intersects(q), std::back_inserter(result));
        }
    }

    auto t3 = std::chrono::high_resolution_clock::now();
    auto elapsed2 = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
    cerr << "Queries completed in: " << elapsed2.count() << " milliseconds" << '\n';
}


void test_dynamic_method(){
    typedef point<float, 2, cartesian> point;

    vector<point> points;

    // Construct using "pack" method (bulk loading)
    rtree<point, quadratic<30> > tree;

    auto t1 = std::chrono::high_resolution_clock::now();
    for (size_t i=0; i<800; i++){
        for (size_t j=0; j<80; j++) {
            tree.insert(point {i+0.5f,j+0.5f});
        }
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    cerr << "Construction completed in: " << elapsed.count() << " milliseconds" << '\n';

    std::vector<point> result;

    for (size_t i=0; i<800; i++){
        for (size_t j=0; j<80; j++) {
            box q(point{i,j},point{i+10,j+10});

            tree.query(intersects(q), std::back_inserter(result));
        }
    }

    auto t3 = std::chrono::high_resolution_clock::now();
    auto elapsed2 = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
    cerr << "Queries completed in: " << elapsed2.count() << " milliseconds" << '\n';
}


void test_object_storage(){
    typedef point<size_t, 2, cartesian> xy_point;
    typedef pair<xy_point, size_t> xy_object;

    vector<xy_object> points;

    // Construct using "pack" method (bulk loading)
    rtree<xy_object, quadratic<30> > tree;

    for (size_t i=0; i<40; i++){
        for (size_t j=0; j<40; j++) {
            tree.insert(make_pair(xy_point {i,j}, i*j));
        }
    }

    cerr << '\n';
    for (size_t i=0; i<40; i++){
        for (size_t j=0; j<40; j++) {
            std::vector<xy_object> result;

            box q(xy_point{i,j},xy_point{i+2,j+2});

            tree.query(intersects(q), std::back_inserter(result));

//            for (const auto& item: result){
//                cerr << i << ' ' << j << " (" << item.first.get<0>() << ',' << item.first.get<1>() << ")=" << item.second << '\n';
//            }
//            cerr << '\n';
        }
    }

}



int main()
{
    test_pack_method();
    test_dynamic_method();
    test_object_storage();
    return 0;
}
