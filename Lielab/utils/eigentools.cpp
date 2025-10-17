#include "eigentools.hpp"

namespace Lielab::utils
{

template Eigen::VectorXi concatenate<int>(const std::vector<Eigen::VectorXi>& vlist);
template Eigen::VectorXd concatenate<double>(const std::vector<Eigen::VectorXd>& vlist);
template Eigen::VectorXi concatenate<int>(std::initializer_list<Eigen::VectorXi> vlist);
template Eigen::VectorXd concatenate<double>(std::initializer_list<Eigen::VectorXd> vlist);

template Eigen::VectorXi arange<int>(const int n0, const int nf);
template Eigen::VectorXd arange<double>(const double n0, const double nf);
template Eigen::VectorXi arange<int>(const int nf);
template Eigen::VectorXd arange<double>(const double nf);

template Eigen::VectorXi repeat<int>(const Eigen::VectorXi& vec, const ptrdiff_t n);
template Eigen::VectorXd repeat<double>(const Eigen::VectorXd& vec, const ptrdiff_t n);
template Eigen::VectorXi repeat<int>(std::initializer_list<int> vec, const ptrdiff_t n);
template Eigen::VectorXd repeat<double>(std::initializer_list<double> vec, const ptrdiff_t n);

template Eigen::VectorXi tile<int>(const Eigen::VectorXi& vec, const ptrdiff_t n);
template Eigen::VectorXd tile<double>(const Eigen::VectorXd& vec, const ptrdiff_t n);
template Eigen::VectorXi tile<int>(std::initializer_list<int> vec, const ptrdiff_t n);
template Eigen::VectorXd tile<double>(std::initializer_list<double> vec, const ptrdiff_t n);

template Eigen::VectorXd linspace<double>(const double start, const double stop, const ptrdiff_t sz);
template Eigen::VectorXd logspace<double>(const double start, const double stop, const ptrdiff_t sz);

template Eigen::MatrixXi column_stack<int>(const std::vector<Eigen::VectorXi>& vlist);
template Eigen::MatrixXd column_stack<double>(const std::vector<Eigen::VectorXd>& vlist);
template Eigen::MatrixXi column_stack<int>(std::initializer_list<Eigen::VectorXi> vlist);
template Eigen::MatrixXd column_stack<double>(std::initializer_list<Eigen::VectorXd> vlist);

template Eigen::MatrixXd vertical_stack<double>(const std::vector<Eigen::MatrixXd>& mlist);
template Eigen::MatrixXd vertical_stack<double>(std::initializer_list<Eigen::MatrixXd> mlist);

Eigen::VectorXd to_VectorXd(std::initializer_list<double> vec)
{
    const size_t sz = vec.size();

    Eigen::VectorXd out(sz);
    size_t ii = 0;
    for (double val : vec)
    {
        out(ii) = val;
        ii++;
    }

    return out;
}

template Eigen::MatrixXd linear_interpolate(const Eigen::VectorXd& x_interp, const Eigen::VectorXd& x, const Eigen::MatrixXd& y);

}
