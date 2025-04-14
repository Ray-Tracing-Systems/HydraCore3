#include "../neural.h"
#include <iostream>
#include <vector>

int main()
{
    std::vector<float> A{1, 0, 0, 0,
                         0, 2, 0, 0,
                         0, 0, 1, 0,
                         0, 0, 0, 1};

    std::vector<float> B{1, 2,
                         3, 4,
                         5, 6,
                         7, 8};

    float out[4 * 2];
    nn::Matmul(A.data(), B.data(), out, 4, 4, 2);

    for(int i = 0; i < 4 * 2; ++i) {
        std::cout << out[i] << " ";
    }
    std::cout << std::endl;
}
