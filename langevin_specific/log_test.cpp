#include <iostream>
#include <cmath>
#include <chrono>

const int X_SIZE = 5;

int enzyme_dup;
int enzyme_out;
int enzyme_const;

void __enzyme_autodiff(...);

// function featuring logarithm of product
double f1(double *x) {
    double y = 1.0;
    for (int i = 0; i < X_SIZE; i++) {
        y *= x[i];
    }
    return log(y);
}
void nabla_f1(double *x, double *dx) {
    __enzyme_autodiff(
        (void*) f1,
        enzyme_dup, x, dx
    );
}

// same function, but simplified
double f2(double *x) {
    double y = 0;
    for (int i = 0; i < X_SIZE; i++) {
        y += log(x[i]);
    }
    return y;
}
void nabla_f2(double *x, double *dx) {
    __enzyme_autodiff(
        (void*) f2,
        enzyme_dup, x, dx
    );
}

int main() {
    double x[X_SIZE];
    for (int i = 0; i < X_SIZE; i++) {
        x[i] = i + 1;
    }
    double dx[X_SIZE];

    // ==== differentiate f1 ====
    {
        std::cout << "Differentiate `log(x_1 * ... * x_" << X_SIZE << ")`" << std::endl;

        for (int i = 0; i < X_SIZE; i++) {
            dx[i] = 0.0;
        }

        // measure time
        auto start = std::chrono::high_resolution_clock::now();
        nabla_f1(x, dx);
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Elapsed time: " << (end - start).count() << "ms" << std::endl;

        std::cout << "Gradient: ( ";
        for (int i = 0; i < X_SIZE; i++) {
            std::cout << dx[i] << " ";
        }
        std::cout << ")" << std::endl;
    }

    std::cout << std::endl;
    
    // ==== differentiate f2 ==== 
    {
        std::cout << "Differentiate `log(x_1) + ... + log(x_" << X_SIZE << ")`" << std::endl;

        for (int i = 0; i < X_SIZE; i++) {
            dx[i] = 0.0;
        }

        auto start = std::chrono::high_resolution_clock::now();
        nabla_f2(x, dx);
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Elapsed time: " << (end - start).count() << "ms" << std::endl;

        std::cout << "Gradient: ( ";
        for (int i = 0; i < sizeof(dx) / sizeof(dx[0]); i++) {
            std::cout << dx[i] << " ";
        }
        std::cout << ")" << std::endl;    
    }
}