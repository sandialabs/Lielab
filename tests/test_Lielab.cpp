#define CATCH_CONFIG

#include <catch2/catch_session.hpp>

#include <Eigen/Core>
#include <Lielab.hpp>

#include <iostream>

int main(int argc, char* argv[])
{
    std::cout << "Testing Lielab " << Lielab::VERSION << std::endl;
    std::cout << "SIMD instruction sets in use: " << Lielab::get_simd_info() << std::endl;

    int result = Catch::Session().run(argc, argv);
    return result;
}
