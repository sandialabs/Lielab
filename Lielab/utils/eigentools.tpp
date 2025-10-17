#ifndef LIELAB_UTILS_EIGENTOOLS_TPP
#define LIELAB_UTILS_EIGENTOOLS_TPP

#include "eigentools.hpp"

#include <Eigen/Core>

namespace Lielab::utils
{

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> concatenate(const std::vector<Eigen::Matrix<F, Eigen::Dynamic, 1>>& vlist)
{
    using vec_t = Eigen::Matrix<F, Eigen::Dynamic, 1>;

    if (vlist.size() == 0)
    {
        return vec_t::Zero(0);
    }

    size_t sz = 0;
    for (size_t ii = 0; ii < vlist.size(); ii++)
    {
        sz += vlist[ii].size();
    }

    vec_t out = vec_t::Zero(sz);
    size_t ind0 = 0;

    for (size_t ii = 0; ii < vlist.size(); ii++)
    {
        const size_t si = vlist[ii].size();
        out(Eigen::seqN(ind0, si)) = vlist[ii];
        ind0 += si;
    }

    return out;
}

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> concatenate(std::initializer_list<Eigen::Matrix<F, Eigen::Dynamic, 1>> vlist)
{
    const std::vector<Eigen::Matrix<F, Eigen::Dynamic, 1>> vlist_vec = vlist;
    return concatenate<F>(vlist_vec);
}

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> arange(const F n0, const F nf)
{
    using vec_t = Eigen::Matrix<F, Eigen::Dynamic, 1>;

    if (n0 >= nf) return vec_t::Zero(0);

    const ptrdiff_t sz = static_cast<ptrdiff_t>(std::ceil(nf - n0));

    vec_t out = vec_t::Zero(sz);

    for (ptrdiff_t ii = 0; ii < sz; ii++)
    {
        out(ii) = static_cast<F>(ii) + n0;
    }
    return out;
}

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> arange(const F nf)
{
    return arange<F>(F(0), nf);
}

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> repeat(const Eigen::Matrix<F, Eigen::Dynamic, 1>& vec, const ptrdiff_t n)
{
    using vec_t = Eigen::Matrix<F, Eigen::Dynamic, 1>;

    if (n <= 0) return vec_t::Zero(0);

    const ptrdiff_t sz = vec.size();

    if (sz == 0) return vec_t::Zero(0);

    vec_t out = vec_t::Zero(sz*n);
    for (ptrdiff_t ii = 0; ii < sz; ii++)
    {
        for (ptrdiff_t jj = 0; jj < n; jj++)
        {
            out(ii*n + jj) = vec(ii);
        }
    }
    return out;
}

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> repeat(std::initializer_list<F> vlist, const ptrdiff_t n)
{
    using vec_t = Eigen::Matrix<F, Eigen::Dynamic, 1>;
    
    // TODO: This initialization is very slow
    vec_t vlist_vec = vec_t::Zero(vlist.size());
    size_t ii = 0;
    for (F val : vlist)
    {
        vlist_vec(ii) = val;
        ii++;
    }

    return repeat<F>(vlist_vec, n);
}

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> tile(const Eigen::Matrix<F, Eigen::Dynamic, 1>& vec, const ptrdiff_t n)
{
    using vec_t = Eigen::Matrix<F, Eigen::Dynamic, 1>;

    if (n <= 0) return vec_t::Zero(0);

    const ptrdiff_t sz = vec.size();

    if (sz == 0) return vec_t::Zero(0);

    vec_t out = vec_t::Zero(sz*n);
    for (ptrdiff_t ii = 0; ii < n; ii++)
    {
        for (ptrdiff_t jj = 0; jj < sz; jj++)
        {
            out(ii*sz + jj) = vec(jj);
        }
    }
    return out;
}

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> tile(std::initializer_list<F> vlist, const ptrdiff_t n)
{
    using vec_t = Eigen::Matrix<F, Eigen::Dynamic, 1>;
    
    // TODO: This initialization is very slow
    vec_t vlist_vec = vec_t::Zero(vlist.size());
    ptrdiff_t ii = 0;
    for (F val : vlist)
    {
        vlist_vec(ii) = val;
        ii++;
    }

    return tile<F>(vlist_vec, n);
}

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> linspace(const F start, const F stop, const ptrdiff_t sz)
{
    using vec_t = Eigen::Matrix<F, Eigen::Dynamic, 1>;

    if (sz <= 0) return vec_t::Zero(0);

    vec_t out = vec_t::Zero(sz);

    out(0) = start;

    if (sz == 1) return out;

    const F m = (stop - start)/static_cast<F>(sz - 1);

    for (ptrdiff_t ii = 1; ii < sz - 1; ii++)
    {
        out(ii) = start + m*static_cast<double>(ii);
    }

    out(sz - 1) = stop;

    return out;
}

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, 1> logspace(const F start, const F stop, const ptrdiff_t sz)
{
    using vec_t = Eigen::Matrix<F, Eigen::Dynamic, 1>;

    if (sz <= 0) return vec_t::Zero(0);

    vec_t powspace = linspace<F>(start, stop, sz);
    vec_t out = vec_t::Zero(sz);

    for (ptrdiff_t ii = 0; ii < sz; ii++)
    {
        out(ii) = std::pow(F(10.0), powspace(ii));
    }

    return out;
}

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> column_stack(const std::vector<Eigen::Matrix<F, Eigen::Dynamic, 1>>& vlist)
{
    using vec_t = Eigen::Matrix<F, Eigen::Dynamic, 1>;
    using mat_t = Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>;

    mat_t out = mat_t::Zero(0, 0);
    if (vlist.size() == 0) return out;

    ptrdiff_t sz = vlist[0].size();

    for (const vec_t& vec : vlist)
    {
        const ptrdiff_t vz = vec.size();
        if (vz < sz) sz = vz;
    }

    out = mat_t::Zero(sz, vlist.size());

    ptrdiff_t ii = 0;
    for (const vec_t& vec : vlist)
    {
        out.col(ii) = vec(Eigen::seqN(0, sz));
        ii++;
    }

    return out;
}

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> column_stack(std::initializer_list<Eigen::Matrix<F, Eigen::Dynamic, 1>> vlist)
{
    const std::vector<Eigen::Matrix<F, Eigen::Dynamic, 1>> vlist_vec = vlist;
    return column_stack<F>(vlist_vec);
}

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> vertical_stack(const std::vector<Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>>& mlist)
{
    using mat_t = Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>;

    if (mlist.size() == 0) return mat_t::Zero(0, 0);

    ptrdiff_t height = mlist[0].rows();
    ptrdiff_t width = mlist[0].cols();

    for (ptrdiff_t ii = 1; ii < static_cast<ptrdiff_t>(mlist.size()); ii++)
    {
        height += mlist[ii].rows();
        if (mlist[ii].cols() < width) width = mlist[ii].cols();
    }

    mat_t out(height, width);

    size_t rowind = 0;

    for (size_t ii = 0; ii < mlist.size(); ii++)
    {
        const size_t nrows = mlist[ii].rows();
        out.block(rowind, 0, nrows, width) = mlist[ii].block(0, 0, nrows, width);
        rowind += nrows;
    }

    return out;
}

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> vertical_stack(std::initializer_list<Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>> mlist)
{
    const std::vector<Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>> mlist_vec = mlist;
    return vertical_stack<F>(mlist_vec);
}

template <typename F>
Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> linear_interpolate(const Eigen::Matrix<F, Eigen::Dynamic, 1>& x_interp, const Eigen::Matrix<F, Eigen::Dynamic, 1>& x, const Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>& y)
{
    // using Eigen::placeholders::last;

    const ptrdiff_t n_interp = x_interp.size();
    const ptrdiff_t n_cols = y.cols();
    const ptrdiff_t n_data = y.rows();

    // TODO: Check n_data == x.size()

    const ptrdiff_t last = n_data - 1;

    // TODO: Check last != 0, or check data has at least 2 points.

    Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> out(n_interp, n_cols);

    for (ptrdiff_t ii = 0; ii < n_interp; ii++)
    {
        if (x_interp(ii) < x(0))
        {
            out.row(ii) = (x_interp(ii) - x(0))*(y.row(1) - y.row(0))/(x(1) - x(0)) + y.row(0);
        }
        else if (x_interp(ii) >= x(last))
        {
            out.row(ii) = (x_interp(ii) - x(last - 1))*(y.row(last) - y.row(last-1))/(x(last) - x(last-1)) + y.row(last-1);
        }
        else
        {
            for (size_t jj = 0; jj < n_data - 1; jj++)
            {
                if (x_interp(ii) >= x(jj) && x_interp(ii) < x(jj+1))
                {
                    out.row(ii) = (x_interp(ii) - x(jj))*(y.row(jj+1) - y.row(jj))/(x(jj+1) - x(jj)) + y.row(jj);
                }
            }
        }
    }

    return out;
}

}

#endif
