#ifndef LIELAB_UTILS_EIGENTOOLS_HPP
#define LIELAB_UTILS_EIGENTOOLS_HPP

#include <Eigen/Core>

#include <initializer_list>
#include <vector>

namespace Lielab::utils
{

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> concatenate(const std::vector<Eigen::Matrix<F, Eigen::Dynamic, 1>>& vlist);

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> concatenate(std::initializer_list<Eigen::Matrix<F, Eigen::Dynamic, 1>> vlist);

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> arange(const F n0, const F nf);

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> arange(const F nf);

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> repeat(const Eigen::Matrix<F, Eigen::Dynamic, 1>& vec, const ptrdiff_t n);

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> repeat(std::initializer_list<F> vec, const ptrdiff_t n);

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> tile(const Eigen::Matrix<F, Eigen::Dynamic, 1>& vec, const ptrdiff_t n);

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> tile(std::initializer_list<F> vec, const ptrdiff_t n);

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> linspace(const F start, const F stop, const ptrdiff_t sz);

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> column_stack(const std::vector<Eigen::Matrix<F, Eigen::Dynamic, 1>>& vlist);

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> column_stack(std::initializer_list<Eigen::Matrix<F, Eigen::Dynamic, 1>> vlist);

Eigen::MatrixXd vertical_stack(const std::vector<Eigen::MatrixXd>& mlist);
Eigen::MatrixXd vertical_stack(std::initializer_list<Eigen::MatrixXd> mlist);

}

#include "eigentools.tpp"

#endif
