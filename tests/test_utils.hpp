#ifndef LIELAB_TEST_HPP
#define LIELAB_TEST_HPP

#include <iostream>

#include <Eigen/Core>
#include <Lielab.hpp>

// Tolerances
constexpr double TOL_FINE = 1e-6;
constexpr double TOL_COARSE = 1e-3;

// Types of numbers
constexpr int an_int = 2;
constexpr float a_float = 2.0;
constexpr double a_double = 2.0;
constexpr std::complex<int> an_imag_int(0, 2);
constexpr std::complex<float> an_imag_float(0.0, 2.0);
constexpr std::complex<double> an_imag_double(0.0, 2.0);

void assert_matrix(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2);
void assert_matrix(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat1, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat2);
void assert_matrix(Eigen::MatrixXcd mat1, Eigen::MatrixXcd mat2);
void assert_matrix(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat1, Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat2);

template<typename L>
void assert_domain(L mat1, L mat2);


template <typename T>
void is_algebra(const std::vector<T> & basis);

template <typename T>
void is_liealgebra(const std::vector<T> & basis);

template <typename T>
void is_group(const std::vector<T> & elements, const T & identity);


#include "test_utils.tpp"

#endif
