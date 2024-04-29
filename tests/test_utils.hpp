#ifndef _LIELAB_TEST_H
#define _LIELAB_TEST_H

#include <iostream>

#include <Eigen/Core>
#include <Lielab.hpp>

// Tolerances
double TOL_FINE = 1e-6;
double TOL_COARSE = 1e-3;

// Types of numbers
int an_int = 2;
float a_float = 2.0;
double a_double = 2.0;
std::complex<int> an_imag_int(0, 2);
std::complex<float> an_imag_float(0.0, 2.0);
std::complex<double> an_imag_double(0.0, 2.0);

void assert_matrix(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2)
{
    const size_t nrows = mat1.rows();
    const size_t ncols = mat1.cols();

    CHECK(nrows == mat2.rows());
    CHECK(ncols == mat2.cols());

    for (int ii = 0; ii < nrows; ii++)
    {
        for (int jj = 0; jj < ncols; jj++)
        {
            CHECK(std::abs(mat1(ii,jj) - mat2(ii,jj)) < TOL_FINE);
        }
    }
}

void assert_matrix(Eigen::MatrixXcd mat1, Eigen::MatrixXcd mat2)
{
    const size_t nrows = mat1.rows();
    const size_t ncols = mat1.cols();

    CHECK(nrows == mat2.rows());
    CHECK(ncols == mat2.cols());

    for (int ii = 0; ii < nrows; ii++)
    {
        for (int jj = 0; jj < ncols; jj++)
        {
            CHECK(std::abs(mat1(ii,jj) - mat2(ii,jj)) < TOL_FINE);
        }
    }
}

template<typename L>
void assert_domain(L mat1, L mat2)
{
    assert_matrix(mat1.get_matrix(), mat2.get_matrix());
}

#endif
