#include "Kernel.hpp"

#include <cmath>


void GaussianSmoother::print_kernel(ostream& o) {
    size_t index = 0;
    for (int64_t i=-width; i<=width; i++){
        cerr << i << '\t' << kernel[index] << '\n';
        index++;
    }
}


GaussianSmoother::GaussianSmoother(int64_t width, double sigma):
        kernel(2*width+1),
        width(width),
        sigma(sigma)
{
    generate_kernel();
}


void GaussianSmoother::generate_kernel(){
    if (kernel.size() % 2 != 1){
        throw runtime_error("ERROR: kernel size must be odd number");
    }

    double e = 2.71828;
    double pi = 3.14159;

    double a;
    double b;

    a = 1/(sigma*sqrt(double(2)*pi));

    size_t index = 0;
    for (int64_t i=-width; i<=width; i++){
        b = (-0.5)*pow(double(i)/sigma, 2.0);

        kernel[index] = a*pow(e,b);

        index++;
    }
}
