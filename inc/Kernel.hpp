#ifndef OVERLAP_ANALYSIS_KERNEL_HPP
#define OVERLAP_ANALYSIS_KERNEL_HPP

#include <stdexcept>
#include <iostream>
#include <vector>

using std::runtime_error;
using std::ostream;
using std::vector;
using std::cerr;
using std::max;


class GaussianSmoother {
private:
    vector<double> kernel;
    const int64_t width;
    const double sigma;

public:
    GaussianSmoother(int64_t width, double sigma);
    void generate_kernel();
    template<class T>
    void smooth(vector<T>& x, vector<T>& y);
    void print_kernel(ostream& o);
};


template <class T> void GaussianSmoother::smooth(vector<T>& x, vector<T>& y){
    if (x.size() != y.size()){
        throw runtime_error("ERROR: output vector must be same size as input");
    }

    for (int64_t i=0; i<int64_t(x.size()); i++) {
        size_t kernel_index = 0;

        for (int64_t j=-width; j<=width; j++){
            if (i+j > 0 and i+j < int64_t(x.size())){
                y[i] += kernel[kernel_index]*x[i+j];
            }
            kernel_index++;
        }
    }
}


#endif //OVERLAP_ANALYSIS_KERNEL_HPP
