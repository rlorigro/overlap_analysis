#include "Kernel.hpp"
#include <iostream>

using std::cerr;


int main(){
    GaussianSmoother smoother(4,1);

    smoother.print_kernel(cerr);

    vector<size_t> x = {0,0,0,0,100,0,0,0,0};
    vector<size_t> y(x.size());

    smoother.smooth(x,y);

    for (auto item: x){
        cerr << item << ' ';
    }
    cerr << '\n';

    for (auto item: y){
        cerr << item << ' ';
    }
    cerr << '\n';

    GaussianSmoother smoother2(5,2);
    smoother2.print_kernel(cerr);

    return 0;
}