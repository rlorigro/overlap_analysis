#include "boost/icl/split_interval_map.hpp"
#include "boost/icl/interval_map.hpp"
#include "boost/icl/interval.hpp"

using boost::icl::partial_enricher;
using boost::icl::total_enricher;
using boost::icl::partial_absorber;
using boost::icl::total_absorber;
using boost::icl::split_interval_map;
using boost::icl::interval_map;
using boost::icl::interval;

#include <iostream>
#include <set>

using std::cerr;
using std::set;


int main(){
    cerr << "TESTING interval_map <uint32_t, set<uint32_t>, partial_enricher>";
    {
        interval_map <uint32_t, set<uint32_t>, partial_enricher> m;
        uint32_t count = 0;

        {
            auto a = interval<uint32_t>::right_open(0, 10);
            set<uint32_t> b = {++count};
            m.add({a, b});
        }
        {
            auto a = interval<uint32_t>::right_open(5, 15);
            set<uint32_t> b = {++count};
            m.add({a, b});
        }
        {
            auto a = interval<uint32_t>::right_open(10, 20);
            set<uint32_t> b = {++count};
            m.add({a, b});
        }

        {
            uint32_t i = 1;
            cerr << "finding: " << i << '\n';
            auto iter = m.find(i);
            cerr << iter->first << " -> ";
            for (auto& item: iter->second){
                cerr << item << " " ;
            }
            cerr << '\n';
        }
        {
            uint32_t i = 6;
            cerr << "finding: " << i << '\n';
            auto iter = m.find(i);
            cerr << iter->first << " -> ";
            for (auto& item: iter->second){
                cerr << item << " " ;
            }
            cerr << '\n';
        }
        {
            uint32_t i = 11;
            cerr << "finding: " << i << '\n';
            auto iter = m.find(i);
            cerr << iter->first << " -> ";
            for (auto& item: iter->second){
                cerr << item << " " ;
            }
            cerr << '\n';
        }
        {
            uint32_t i = 16;
            cerr << "finding: " << i << '\n';
            auto iter = m.find(i);
            cerr << iter->first << " -> ";
            for (auto& item: iter->second){
                cerr << item << " " ;
            }
            cerr << '\n';
        }
        cerr << '\n';

        for (auto& item: m){
            cerr << item.first << " -> ";
            for (auto& item2: item.second){
                cerr << item2 << " " ;
            }
            cerr << '\n';
        }

    }

    cerr << "TESTING split_interval_map <uint32_t, set<uint32_t>, partial_enricher>";

    {
        split_interval_map <uint32_t, set<uint32_t>, partial_enricher> m;
        uint32_t count = 0;

        {
            auto a = interval<uint32_t>::right_open(0, 10);
            set<uint32_t> b = {++count};
            m.add({a, b});
        }
        {
            auto a = interval<uint32_t>::right_open(5, 15);
            set<uint32_t> b = {++count};
            m.add({a, b});
        }
        {
            auto a = interval<uint32_t>::right_open(10, 20);
            set<uint32_t> b = {++count};
            m.add({a, b});
        }

        {
            uint32_t i = 1;
            cerr << "finding: " << i << '\n';
            auto iter = m.find(i);
            cerr << iter->first << " -> ";
            for (auto& item: iter->second){
                cerr << item << " " ;
            }
            cerr << '\n';
        }
        {
            uint32_t i = 6;
            cerr << "finding: " << i << '\n';
            auto iter = m.find(i);
            cerr << iter->first << " -> ";
            for (auto& item: iter->second){
                cerr << item << " " ;
            }
            cerr << '\n';
        }
        {
            uint32_t i = 11;
            cerr << "finding: " << i << '\n';
            auto iter = m.find(i);
            cerr << iter->first << " -> ";
            for (auto& item: iter->second){
                cerr << item << " " ;
            }
            cerr << '\n';
        }
        {
            uint32_t i = 16;
            cerr << "finding: " << i << '\n';
            auto iter = m.find(i);
            cerr << iter->first << " -> ";
            for (auto& item: iter->second){
                cerr << item << " " ;
            }
            cerr << '\n';
        }
        cerr << '\n';

        for (auto& item: m){
            cerr << item.first << " -> ";
            for (auto& item2: item.second){
                cerr << item2 << " " ;
            }
            cerr << '\n';
        }

    }

    return 0;
}
