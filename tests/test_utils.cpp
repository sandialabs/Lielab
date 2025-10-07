#include <iostream>

#include <catch2/catch_all.hpp>

#include <Eigen/Core>
#include <Lielab.hpp>
#include "test_utils.hpp"

void assert_matrix(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2)
{
    const ptrdiff_t nrows = mat1.rows();
    const ptrdiff_t ncols = mat1.cols();

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

void assert_matrix(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat1, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat2)
{
    const ptrdiff_t nrows = mat1.rows();
    const ptrdiff_t ncols = mat1.cols();

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
    const ptrdiff_t nrows = mat1.rows();
    const ptrdiff_t ncols = mat1.cols();

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

void assert_matrix(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat1, Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat2)
{
    const ptrdiff_t nrows = mat1.rows();
    const ptrdiff_t ncols = mat1.cols();

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
