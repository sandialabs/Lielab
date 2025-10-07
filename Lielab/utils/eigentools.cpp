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

template Eigen::MatrixXi column_stack<int>(const std::vector<Eigen::VectorXi>& vlist);
template Eigen::MatrixXd column_stack<double>(const std::vector<Eigen::VectorXd>& vlist);
template Eigen::MatrixXi column_stack<int>(std::initializer_list<Eigen::VectorXi> vlist);
template Eigen::MatrixXd column_stack<double>(std::initializer_list<Eigen::VectorXd> vlist);

Eigen::MatrixXd vertical_stack(const std::vector<Eigen::MatrixXd>& mlist)
{
    if (mlist.size() == 0) return Eigen::MatrixXd::Zero(0, 0);

    size_t height = mlist[0].rows();
    size_t width = mlist[0].cols();

    for (size_t ii = 1; ii < mlist.size(); ii++)
    {
        height += mlist[ii].rows();
        if (mlist[ii].cols() < width) width = mlist[ii].cols();
    }

    Eigen::MatrixXd out(height, width);

    size_t rowind = 0;

    for (size_t ii = 0; ii < mlist.size(); ii++)
    {
        const size_t nrows = mlist[ii].rows();
        out.block(rowind, 0, nrows, width) = mlist[ii].block(0, 0, nrows, width);
        rowind += nrows;
    }

    return out;
}

Eigen::MatrixXd vertical_stack(std::initializer_list<Eigen::MatrixXd> mlist)
{
    const std::vector<Eigen::MatrixXd> mlist_vec = mlist;
    return vertical_stack(mlist_vec);
}

}
