#ifndef _LIELAB_DOMAIN_HALIE_HPP
#define _LIELAB_DOMAIN_HALIE_HPP

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <variant>

#include "abstract.hpp"
#include "utils.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>

// lielab dependencies
#include "constants.hpp"

using lielab::abstract::Real;
using lielab::abstract::Imaginary;

namespace lielab
{

namespace domain
{
class halie
{
    public:
    static constexpr size_t INDEX_gl = 0;
    static constexpr size_t INDEX_rn = 1;
    static constexpr size_t INDEX_so = 2;
    static constexpr size_t INDEX_sp = 3;
    static constexpr size_t INDEX_su = 4;

    typedef std::variant<lielab::domain::gl,
                         lielab::domain::rn,
                         lielab::domain::so,
                         lielab::domain::sp,
                         lielab::domain::su> TYPES;

    std::vector<TYPES> space;

    halie()
    {

    }

    halie(std::initializer_list<TYPES> others)
    {
        /*!
        * Initializer list constructor for halie
        *
        * Enables construction like:
        *     lielab::domain::rn R(3);
        *     lielab::domain::so O(3);
        *     lielab::domain::halie M{R, O};
        */

        this->space = std::vector<TYPES>{std::move(others)};
    }

    halie(const std::vector<TYPES> & others)
    {
        /*!
        * Vector list constructor for halie
        *
        * Not needed for C++, but enables construction in Python like:
        *     R = lielab.domain.rn(3)
        *     O = lielab.domain.so(3)
        *     M = lielab.domain.halie([R, O])
        */

        this->space = others;
    }

    size_t get_dimension() const
    {
        size_t dim = 0;

        for (int ii = 0; ii < space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_gl)
            {
                dim += std::get<lielab::domain::gl>(space[ii]).get_dimension();
            }
            else if (ind == INDEX_rn)
            {
                dim += std::get<lielab::domain::rn>(space[ii]).get_dimension();
            }
            else if (ind == INDEX_so)
            {
                dim += std::get<lielab::domain::so>(space[ii]).get_dimension();
            }
            else if (ind == INDEX_sp)
            {
                dim += std::get<lielab::domain::sp>(space[ii]).get_dimension();
            }
            else if (ind == INDEX_su)
            {
                dim += std::get<lielab::domain::su>(space[ii]).get_dimension();
            }
        }

        return dim;
    }

    Eigen::VectorXi get_shape() const
    {
        Eigen::VectorXi out(space.size());

        for (int ii = 0; ii < space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_gl)
            {
                out(ii) = static_cast<int>(std::get<lielab::domain::gl>(space[ii]).shape);
            }
            else if (ind == INDEX_rn)
            {
                out(ii) = static_cast<int>(std::get<lielab::domain::rn>(space[ii]).shape);
            }
            else if (ind == INDEX_so)
            {
                out(ii) = static_cast<int>(std::get<lielab::domain::so>(space[ii]).shape);
            }
            else if (ind == INDEX_sp)
            {
                out(ii) = static_cast<int>(std::get<lielab::domain::sp>(space[ii]).shape);
            }
            else if (ind == INDEX_su)
            {
                out(ii) = static_cast<int>(std::get<lielab::domain::su>(space[ii]).shape);
            }
        }

        return out;
    }

    Eigen::VectorXd get_vector() const
    {
        std::vector<Eigen::VectorXd> vectors = std::vector<Eigen::VectorXd>(this->space.size());

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_gl)
            {
                vectors[ii] = std::get<lielab::domain::gl>(this->space[ii]).get_vector();
            }
            else if (ind == INDEX_rn)
            {
                vectors[ii] = std::get<lielab::domain::rn>(this->space[ii]).get_vector();
            }
            else if (ind == INDEX_so)
            {
                vectors[ii] = std::get<lielab::domain::so>(this->space[ii]).get_vector();
            }
            else if (ind == INDEX_sp)
            {
                vectors[ii] = std::get<lielab::domain::sp>(this->space[ii]).get_vector();
            }
            else if (ind == INDEX_su)
            {
                vectors[ii] = std::get<lielab::domain::su>(this->space[ii]).get_vector();
            }
        }
        return lielab::utils::concatenate(vectors);
    }

    void set_vector(const Eigen::VectorXd & vec)
    {
        size_t jj = 0;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_gl)
            {
                const size_t dim = std::get<lielab::domain::gl>(this->space[ii]).get_dimension();
                std::get<lielab::domain::gl>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
                jj += dim;
            }
            else if (ind == INDEX_rn)
            {
                const size_t dim = std::get<lielab::domain::rn>(this->space[ii]).get_dimension();
                std::get<lielab::domain::rn>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
                jj += dim;
            }
            else if (ind == INDEX_so)
            {
                const size_t dim = std::get<lielab::domain::so>(this->space[ii]).get_dimension();
                std::get<lielab::domain::so>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
                jj += dim;
            }
            else if (ind == INDEX_sp)
            {
                const size_t dim = std::get<lielab::domain::sp>(this->space[ii]).get_dimension();
                std::get<lielab::domain::sp>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
                jj += dim;
            }
            else if (ind == INDEX_su)
            {
                const size_t dim = std::get<lielab::domain::su>(this->space[ii]).get_dimension();
                std::get<lielab::domain::su>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
                jj += dim;
            }
        }
    }

    halie operator+(const halie & other) const
    {
        assert(this->space.size() == other.space.size());

        halie out;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_gl)
            {
                out.space.push_back(std::get<lielab::domain::gl>(this->space[ii]) + std::get<lielab::domain::gl>(other.space[ii]));
            }
            else if (ind == INDEX_rn)
            {
                out.space.push_back(std::get<lielab::domain::rn>(this->space[ii]) + std::get<lielab::domain::rn>(other.space[ii]));
            }
            else if (ind == INDEX_so)
            {
                out.space.push_back(std::get<lielab::domain::so>(this->space[ii]) + std::get<lielab::domain::so>(other.space[ii]));
            }
            else if (ind == INDEX_sp)
            {
                out.space.push_back(std::get<lielab::domain::sp>(this->space[ii]) + std::get<lielab::domain::sp>(other.space[ii]));
            }
            else if (ind == INDEX_su)
            {
                out.space.push_back(std::get<lielab::domain::su>(this->space[ii]) + std::get<lielab::domain::su>(other.space[ii]));
            }
        }

        return out;
    }

    halie & operator+=(const halie & other)
    {
        assert(this->space.size() == other.space.size());

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_gl)
            {
                std::get<lielab::domain::gl>(this->space[ii]) += std::get<lielab::domain::gl>(other.space[ii]);
            }
            else if (ind == INDEX_rn)
            {
                std::get<lielab::domain::rn>(this->space[ii]) += std::get<lielab::domain::rn>(other.space[ii]);
            }
            else if (ind == INDEX_so)
            {
                std::get<lielab::domain::so>(this->space[ii]) += std::get<lielab::domain::so>(other.space[ii]);
            }
            else if (ind == INDEX_sp)
            {
                std::get<lielab::domain::sp>(this->space[ii]) += std::get<lielab::domain::sp>(other.space[ii]);
            }
            else if (ind == INDEX_su)
            {
                std::get<lielab::domain::su>(this->space[ii]) += std::get<lielab::domain::su>(other.space[ii]);
            }
        }

        return *this;
    }

    halie operator-(const halie & other) const
    {
        halie out;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_gl)
            {
                out.space.push_back(std::get<lielab::domain::gl>(this->space[ii]) - std::get<lielab::domain::gl>(other.space[ii]));
            }
            else if (ind == INDEX_rn)
            {
                out.space.push_back(std::get<lielab::domain::rn>(this->space[ii]) - std::get<lielab::domain::rn>(other.space[ii]));
            }
            else if (ind == INDEX_so)
            {
                out.space.push_back(std::get<lielab::domain::so>(this->space[ii]) - std::get<lielab::domain::so>(other.space[ii]));
            }
            else if (ind == INDEX_sp)
            {
                out.space.push_back(std::get<lielab::domain::sp>(this->space[ii]) - std::get<lielab::domain::sp>(other.space[ii]));
            }
            else if (ind == INDEX_su)
            {
                out.space.push_back(std::get<lielab::domain::su>(this->space[ii]) - std::get<lielab::domain::su>(other.space[ii]));
            }
        }

        return out;
    }

    halie & operator-=(const halie & other)
    {
        assert(this->space.size() == other.space.size());

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_gl)
            {
                std::get<lielab::domain::gl>(this->space[ii]) -= std::get<lielab::domain::gl>(other.space[ii]);
            }
            else if (ind == INDEX_rn)
            {
                std::get<lielab::domain::rn>(this->space[ii]) -= std::get<lielab::domain::rn>(other.space[ii]);
            }
            else if (ind == INDEX_so)
            {
                std::get<lielab::domain::so>(this->space[ii]) -= std::get<lielab::domain::so>(other.space[ii]);
            }
            else if (ind == INDEX_sp)
            {
                std::get<lielab::domain::sp>(this->space[ii]) -= std::get<lielab::domain::sp>(other.space[ii]);
            }
            else if (ind == INDEX_su)
            {
                std::get<lielab::domain::su>(this->space[ii]) -= std::get<lielab::domain::su>(other.space[ii]);
            }
        }

        return *this;
    }

    halie operator-() const
    {
        halie out;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_gl)
            {
                out.space.push_back(-std::get<lielab::domain::gl>(this->space[ii]));
            }
            else if (ind == INDEX_rn)
            {
                out.space.push_back(-std::get<lielab::domain::rn>(this->space[ii]));
            }
            else if (ind == INDEX_so)
            {
                out.space.push_back(-std::get<lielab::domain::so>(this->space[ii]));
            }
            else if (ind == INDEX_sp)
            {
                out.space.push_back(-std::get<lielab::domain::sp>(this->space[ii]));
            }
            else if (ind == INDEX_su)
            {
                out.space.push_back(-std::get<lielab::domain::su>(this->space[ii]));
            }
        }

        return out;
    }

    halie operator*(const Real auto other) const
    {
        halie out;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_gl)
            {
                out.space.push_back(std::get<lielab::domain::gl>(this->space[ii]) * other);
            }
            else if (ind == INDEX_rn)
            {
                out.space.push_back(std::get<lielab::domain::rn>(this->space[ii]) * other);
            }
            else if (ind == INDEX_so)
            {
                out.space.push_back(std::get<lielab::domain::so>(this->space[ii]) * other);
            }
            else if (ind == INDEX_sp)
            {
                out.space.push_back(std::get<lielab::domain::sp>(this->space[ii]) * other);
            }
            else if (ind == INDEX_su)
            {
                out.space.push_back(std::get<lielab::domain::su>(this->space[ii]) * other);
            }
        }

        return out;
    }

    halie & operator*=(const Real auto other)
    {
        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_gl)
            {
                std::get<lielab::domain::gl>(this->space[ii]) *= other;
            }
            else if (ind == INDEX_rn)
            {
                std::get<lielab::domain::rn>(this->space[ii]) *= other;
            }
            else if (ind == INDEX_so)
            {
                std::get<lielab::domain::so>(this->space[ii]) *= other;
            }
            else if (ind == INDEX_sp)
            {
                std::get<lielab::domain::sp>(this->space[ii]) *= other;
            }
            else if (ind == INDEX_su)
            {
                std::get<lielab::domain::su>(this->space[ii]) *= other;
            }
        }

        return *this;
    }

    halie operator*(const halie & other) const
    {
        halie out;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_gl)
            {
                out.space.push_back(std::get<lielab::domain::gl>(this->space[ii]) * std::get<lielab::domain::gl>(other.space[ii]));
            }
            else if (ind == INDEX_rn)
            {
                out.space.push_back(std::get<lielab::domain::rn>(this->space[ii]) * std::get<lielab::domain::rn>(other.space[ii]));
            }
            else if (ind == INDEX_so)
            {
                out.space.push_back(std::get<lielab::domain::so>(this->space[ii]) * std::get<lielab::domain::so>(other.space[ii]));
            }
            else if (ind == INDEX_sp)
            {
                out.space.push_back(std::get<lielab::domain::sp>(this->space[ii]) * std::get<lielab::domain::sp>(other.space[ii]));
            }
            else if (ind == INDEX_su)
            {
                out.space.push_back(std::get<lielab::domain::su>(this->space[ii]) * std::get<lielab::domain::su>(other.space[ii]));
            }
        }

        return out;
    }

    halie & operator*=(const halie & other)
    {
        assert(this->space.size() == other.space.size());

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_gl)
            {
                std::get<lielab::domain::gl>(this->space[ii]) *= std::get<lielab::domain::gl>(other.space[ii]);
            }
            else if (ind == INDEX_rn)
            {
                std::get<lielab::domain::rn>(this->space[ii]) *= std::get<lielab::domain::rn>(other.space[ii]);
            }
            else if (ind == INDEX_so)
            {
                std::get<lielab::domain::so>(this->space[ii]) *= std::get<lielab::domain::so>(other.space[ii]);
            }
            else if (ind == INDEX_sp)
            {
                std::get<lielab::domain::sp>(this->space[ii]) *= std::get<lielab::domain::sp>(other.space[ii]);
            }
            else if (ind == INDEX_su)
            {
                std::get<lielab::domain::su>(this->space[ii]) *= std::get<lielab::domain::su>(other.space[ii]);
            }
        }

        return *this;
    }

    halie operator/(const Real auto other) const
    {
        halie out;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_gl)
            {
                out.space.push_back(std::get<lielab::domain::gl>(this->space[ii]) / other);
            }
            else if (ind == INDEX_rn)
            {
                out.space.push_back(std::get<lielab::domain::rn>(this->space[ii]) / other);
            }
            else if (ind == INDEX_so)
            {
                out.space.push_back(std::get<lielab::domain::so>(this->space[ii]) / other);
            }
            else if (ind == INDEX_sp)
            {
                out.space.push_back(std::get<lielab::domain::sp>(this->space[ii]) / other);
            }
            else if (ind == INDEX_su)
            {
                out.space.push_back(std::get<lielab::domain::su>(this->space[ii]) / other);
            }
        }

        return out;
    }

    halie & operator/=(const Real auto other)
    {
        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_gl)
            {
                std::get<lielab::domain::gl>(this->space[ii]) /= other;
            }
            else if (ind == INDEX_rn)
            {
                std::get<lielab::domain::rn>(this->space[ii]) /= other;
            }
            else if (ind == INDEX_so)
            {
                std::get<lielab::domain::so>(this->space[ii]) /= other;
            }
            else if (ind == INDEX_sp)
            {
                std::get<lielab::domain::sp>(this->space[ii]) /= other;
            }
            else if (ind == INDEX_su)
            {
                std::get<lielab::domain::su>(this->space[ii]) /= other;
            }
        }

        return *this;
    }

    std::string to_string() const
    {
        std::string out = "";
        const size_t sz = this->space.size();
        Eigen::VectorXi shapes = this->get_shape();

        for (int ii = 0; ii < sz; ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_gl)
            {
                out += "gl(" + std::to_string(shapes(ii)) + ")";
            }
            else if (ind == INDEX_rn)
            {
                out += "rn(" + std::to_string(shapes(ii)) + ")";
            }
            else if (ind == INDEX_so)
            {
                out += "so(" + std::to_string(shapes(ii)) + ")";
            }
            else if (ind == INDEX_sp)
            {
                out += "sp(" + std::to_string(shapes(ii)) + ")";
            }
            else if (ind == INDEX_su)
            {
                out += "su(" + std::to_string(shapes(ii)) + ")";
            }

            if (ii < sz-1)
            {
                out += " x ";
            }
        }

        return out;
    }

    friend std::ostream & operator<<(std::ostream & os, const halie & other);
};

std::ostream & operator<<(std::ostream & os, const halie & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */

    os << other.to_string();
    return os;
}

halie operator*(const Real auto other, const halie & rhs)
{
    return rhs*other;
}
}
}

#endif
