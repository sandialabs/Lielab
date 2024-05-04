#ifndef _LIELAB_DOMAIN_COMPOSITEALGEBRA_HPP
#define _LIELAB_DOMAIN_COMPOSITEALGEBRA_HPP

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <variant>

#include "../abstract.hpp"
#include "../utils.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>

// Lielab dependencies
#include "../constants.hpp"

using Lielab::abstract::Real;
using Lielab::abstract::Imaginary;

namespace Lielab
{

namespace domain
{
class CompositeAlgebra
{
    public:
    static constexpr size_t INDEX_cn  = 0;
    static constexpr size_t INDEX_gl  = 1;
    static constexpr size_t INDEX_glc = 2;
    static constexpr size_t INDEX_rn  = 3;
    static constexpr size_t INDEX_se  = 4;
    static constexpr size_t INDEX_so  = 5;
    static constexpr size_t INDEX_sp  = 6;
    static constexpr size_t INDEX_su  = 7;

    typedef std::variant<Lielab::domain::cn,
                         Lielab::domain::gl,
                         Lielab::domain::glc,
                         Lielab::domain::rn,
                         Lielab::domain::se,
                         Lielab::domain::so,
                         Lielab::domain::sp,
                         Lielab::domain::su> TYPES;

    std::vector<TYPES> space;

    CompositeAlgebra()
    {

    }

    CompositeAlgebra(std::initializer_list<TYPES> others)
    {
        /*!
        * Initializer list constructor for CompositeAlgebra
        *
        * Enables construction like:
        *     Lielab::domain::rn R(3);
        *     Lielab::domain::so O(3);
        *     Lielab::domain::CompositeAlgebra M{R, O};
        */

        this->space = std::vector<TYPES>{std::move(others)};
    }

    CompositeAlgebra(const std::vector<TYPES> & others)
    {
        /*!
        * Vector list constructor for CompositeAlgebra
        *
        * Not needed for C++, but enables construction in Python like:
        *     R = lielab.domain.rn(3)
        *     O = lielab.domain.so(3)
        *     M = lielab.domain.CompositeAlgebra([R, O])
        */

        this->space = others;
    }

    size_t get_dimension() const
    {
        size_t dim = 0;

        for (int ii = 0; ii < space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_cn)
            {
                dim += std::get<Lielab::domain::cn>(space[ii]).get_dimension();
            }
            else if (ind == INDEX_gl)
            {
                dim += std::get<Lielab::domain::gl>(space[ii]).get_dimension();
            }
            else if (ind == INDEX_glc)
            {
                dim += std::get<Lielab::domain::glc>(space[ii]).get_dimension();
            }
            else if (ind == INDEX_rn)
            {
                dim += std::get<Lielab::domain::rn>(space[ii]).get_dimension();
            }
            else if (ind == INDEX_se)
            {
                dim += std::get<Lielab::domain::se>(space[ii]).get_dimension();
            }
            else if (ind == INDEX_so)
            {
                dim += std::get<Lielab::domain::so>(space[ii]).get_dimension();
            }
            else if (ind == INDEX_sp)
            {
                dim += std::get<Lielab::domain::sp>(space[ii]).get_dimension();
            }
            else if (ind == INDEX_su)
            {
                dim += std::get<Lielab::domain::su>(space[ii]).get_dimension();
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
            if (ind == INDEX_cn)
            {
                out(ii) = static_cast<int>(std::get<Lielab::domain::cn>(space[ii]).shape);
            }
            else if (ind == INDEX_gl)
            {
                out(ii) = static_cast<int>(std::get<Lielab::domain::gl>(space[ii]).shape);
            }
            else if (ind == INDEX_glc)
            {
                out(ii) = static_cast<int>(std::get<Lielab::domain::glc>(space[ii]).shape);
            }
            else if (ind == INDEX_rn)
            {
                out(ii) = static_cast<int>(std::get<Lielab::domain::rn>(space[ii]).shape);
            }
            else if (ind == INDEX_se)
            {
                out(ii) = static_cast<int>(std::get<Lielab::domain::se>(space[ii]).shape);
            }
            else if (ind == INDEX_so)
            {
                out(ii) = static_cast<int>(std::get<Lielab::domain::so>(space[ii]).shape);
            }
            else if (ind == INDEX_sp)
            {
                out(ii) = static_cast<int>(std::get<Lielab::domain::sp>(space[ii]).shape);
            }
            else if (ind == INDEX_su)
            {
                out(ii) = static_cast<int>(std::get<Lielab::domain::su>(space[ii]).shape);
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
            if (ind == INDEX_cn)
            {
                vectors[ii] = std::get<Lielab::domain::cn>(this->space[ii]).get_vector();
            }
            else if (ind == INDEX_gl)
            {
                vectors[ii] = std::get<Lielab::domain::gl>(this->space[ii]).get_vector();
            }
            else if (ind == INDEX_glc)
            {
                vectors[ii] = std::get<Lielab::domain::glc>(this->space[ii]).get_vector();
            }
            else if (ind == INDEX_rn)
            {
                vectors[ii] = std::get<Lielab::domain::rn>(this->space[ii]).get_vector();
            }
            else if (ind == INDEX_se)
            {
                vectors[ii] = std::get<Lielab::domain::se>(this->space[ii]).get_vector();
            }
            else if (ind == INDEX_so)
            {
                vectors[ii] = std::get<Lielab::domain::so>(this->space[ii]).get_vector();
            }
            else if (ind == INDEX_sp)
            {
                vectors[ii] = std::get<Lielab::domain::sp>(this->space[ii]).get_vector();
            }
            else if (ind == INDEX_su)
            {
                vectors[ii] = std::get<Lielab::domain::su>(this->space[ii]).get_vector();
            }
        }
        return Lielab::utils::concatenate(vectors);
    }

    void set_vector(const Eigen::VectorXd & vec)
    {
        size_t jj = 0;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_cn)
            {
                const size_t dim = std::get<Lielab::domain::cn>(this->space[ii]).get_dimension();
                std::get<Lielab::domain::cn>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
                jj += dim;
            }
            else if (ind == INDEX_gl)
            {
                const size_t dim = std::get<Lielab::domain::gl>(this->space[ii]).get_dimension();
                std::get<Lielab::domain::gl>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
                jj += dim;
            }
            else if (ind == INDEX_glc)
            {
                const size_t dim = std::get<Lielab::domain::glc>(this->space[ii]).get_dimension();
                std::get<Lielab::domain::glc>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
                jj += dim;
            }
            else if (ind == INDEX_rn)
            {
                const size_t dim = std::get<Lielab::domain::rn>(this->space[ii]).get_dimension();
                std::get<Lielab::domain::rn>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
                jj += dim;
            }
            else if (ind == INDEX_se)
            {
                const size_t dim = std::get<Lielab::domain::se>(this->space[ii]).get_dimension();
                std::get<Lielab::domain::se>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
                jj += dim;
            }
            else if (ind == INDEX_so)
            {
                const size_t dim = std::get<Lielab::domain::so>(this->space[ii]).get_dimension();
                std::get<Lielab::domain::so>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
                jj += dim;
            }
            else if (ind == INDEX_sp)
            {
                const size_t dim = std::get<Lielab::domain::sp>(this->space[ii]).get_dimension();
                std::get<Lielab::domain::sp>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
                jj += dim;
            }
            else if (ind == INDEX_su)
            {
                const size_t dim = std::get<Lielab::domain::su>(this->space[ii]).get_dimension();
                std::get<Lielab::domain::su>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
                jj += dim;
            }
        }
    }

    CompositeAlgebra operator+(const CompositeAlgebra & other) const
    {
        assert(this->space.size() == other.space.size());

        CompositeAlgebra out;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_cn)
            {
                out.space.push_back(std::get<Lielab::domain::cn>(this->space[ii]) + std::get<Lielab::domain::cn>(other.space[ii]));
            }
            else if (ind == INDEX_gl)
            {
                out.space.push_back(std::get<Lielab::domain::gl>(this->space[ii]) + std::get<Lielab::domain::gl>(other.space[ii]));
            }
            else if (ind == INDEX_glc)
            {
                out.space.push_back(std::get<Lielab::domain::glc>(this->space[ii]) + std::get<Lielab::domain::glc>(other.space[ii]));
            }
            else if (ind == INDEX_rn)
            {
                out.space.push_back(std::get<Lielab::domain::rn>(this->space[ii]) + std::get<Lielab::domain::rn>(other.space[ii]));
            }
            else if (ind == INDEX_se)
            {
                out.space.push_back(std::get<Lielab::domain::se>(this->space[ii]) + std::get<Lielab::domain::se>(other.space[ii]));
            }
            else if (ind == INDEX_so)
            {
                out.space.push_back(std::get<Lielab::domain::so>(this->space[ii]) + std::get<Lielab::domain::so>(other.space[ii]));
            }
            else if (ind == INDEX_sp)
            {
                out.space.push_back(std::get<Lielab::domain::sp>(this->space[ii]) + std::get<Lielab::domain::sp>(other.space[ii]));
            }
            else if (ind == INDEX_su)
            {
                out.space.push_back(std::get<Lielab::domain::su>(this->space[ii]) + std::get<Lielab::domain::su>(other.space[ii]));
            }
        }

        return out;
    }

    CompositeAlgebra & operator+=(const CompositeAlgebra & other)
    {
        assert(this->space.size() == other.space.size());

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_cn)
            {
                std::get<Lielab::domain::cn>(this->space[ii]) += std::get<Lielab::domain::cn>(other.space[ii]);
            }
            else if (ind == INDEX_gl)
            {
                std::get<Lielab::domain::gl>(this->space[ii]) += std::get<Lielab::domain::gl>(other.space[ii]);
            }
            else if (ind == INDEX_glc)
            {
                std::get<Lielab::domain::glc>(this->space[ii]) += std::get<Lielab::domain::glc>(other.space[ii]);
            }
            else if (ind == INDEX_rn)
            {
                std::get<Lielab::domain::rn>(this->space[ii]) += std::get<Lielab::domain::rn>(other.space[ii]);
            }
            else if (ind == INDEX_se)
            {
                std::get<Lielab::domain::se>(this->space[ii]) += std::get<Lielab::domain::se>(other.space[ii]);
            }
            else if (ind == INDEX_so)
            {
                std::get<Lielab::domain::so>(this->space[ii]) += std::get<Lielab::domain::so>(other.space[ii]);
            }
            else if (ind == INDEX_sp)
            {
                std::get<Lielab::domain::sp>(this->space[ii]) += std::get<Lielab::domain::sp>(other.space[ii]);
            }
            else if (ind == INDEX_su)
            {
                std::get<Lielab::domain::su>(this->space[ii]) += std::get<Lielab::domain::su>(other.space[ii]);
            }
        }

        return *this;
    }

    CompositeAlgebra operator-(const CompositeAlgebra & other) const
    {
        CompositeAlgebra out;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_cn)
            {
                out.space.push_back(std::get<Lielab::domain::cn>(this->space[ii]) - std::get<Lielab::domain::cn>(other.space[ii]));
            }
            else if (ind == INDEX_gl)
            {
                out.space.push_back(std::get<Lielab::domain::gl>(this->space[ii]) - std::get<Lielab::domain::gl>(other.space[ii]));
            }
            else if (ind == INDEX_glc)
            {
                out.space.push_back(std::get<Lielab::domain::glc>(this->space[ii]) - std::get<Lielab::domain::glc>(other.space[ii]));
            }
            else if (ind == INDEX_rn)
            {
                out.space.push_back(std::get<Lielab::domain::rn>(this->space[ii]) - std::get<Lielab::domain::rn>(other.space[ii]));
            }
            else if (ind == INDEX_se)
            {
                out.space.push_back(std::get<Lielab::domain::se>(this->space[ii]) - std::get<Lielab::domain::se>(other.space[ii]));
            }
            else if (ind == INDEX_so)
            {
                out.space.push_back(std::get<Lielab::domain::so>(this->space[ii]) - std::get<Lielab::domain::so>(other.space[ii]));
            }
            else if (ind == INDEX_sp)
            {
                out.space.push_back(std::get<Lielab::domain::sp>(this->space[ii]) - std::get<Lielab::domain::sp>(other.space[ii]));
            }
            else if (ind == INDEX_su)
            {
                out.space.push_back(std::get<Lielab::domain::su>(this->space[ii]) - std::get<Lielab::domain::su>(other.space[ii]));
            }
        }

        return out;
    }

    CompositeAlgebra & operator-=(const CompositeAlgebra & other)
    {
        assert(this->space.size() == other.space.size());

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_cn)
            {
                std::get<Lielab::domain::cn>(this->space[ii]) -= std::get<Lielab::domain::cn>(other.space[ii]);
            }
            else if (ind == INDEX_gl)
            {
                std::get<Lielab::domain::gl>(this->space[ii]) -= std::get<Lielab::domain::gl>(other.space[ii]);
            }
            else if (ind == INDEX_glc)
            {
                std::get<Lielab::domain::glc>(this->space[ii]) -= std::get<Lielab::domain::glc>(other.space[ii]);
            }
            else if (ind == INDEX_rn)
            {
                std::get<Lielab::domain::rn>(this->space[ii]) -= std::get<Lielab::domain::rn>(other.space[ii]);
            }
            else if (ind == INDEX_se)
            {
                std::get<Lielab::domain::se>(this->space[ii]) -= std::get<Lielab::domain::se>(other.space[ii]);
            }
            else if (ind == INDEX_so)
            {
                std::get<Lielab::domain::so>(this->space[ii]) -= std::get<Lielab::domain::so>(other.space[ii]);
            }
            else if (ind == INDEX_sp)
            {
                std::get<Lielab::domain::sp>(this->space[ii]) -= std::get<Lielab::domain::sp>(other.space[ii]);
            }
            else if (ind == INDEX_su)
            {
                std::get<Lielab::domain::su>(this->space[ii]) -= std::get<Lielab::domain::su>(other.space[ii]);
            }
        }

        return *this;
    }

    CompositeAlgebra operator-() const
    {
        CompositeAlgebra out;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_cn)
            {
                out.space.push_back(-std::get<Lielab::domain::cn>(this->space[ii]));
            }
            else if (ind == INDEX_gl)
            {
                out.space.push_back(-std::get<Lielab::domain::gl>(this->space[ii]));
            }
            else if (ind == INDEX_glc)
            {
                out.space.push_back(-std::get<Lielab::domain::glc>(this->space[ii]));
            }
            else if (ind == INDEX_rn)
            {
                out.space.push_back(-std::get<Lielab::domain::rn>(this->space[ii]));
            }
            else if (ind == INDEX_se)
            {
                out.space.push_back(-std::get<Lielab::domain::se>(this->space[ii]));
            }
            else if (ind == INDEX_so)
            {
                out.space.push_back(-std::get<Lielab::domain::so>(this->space[ii]));
            }
            else if (ind == INDEX_sp)
            {
                out.space.push_back(-std::get<Lielab::domain::sp>(this->space[ii]));
            }
            else if (ind == INDEX_su)
            {
                out.space.push_back(-std::get<Lielab::domain::su>(this->space[ii]));
            }
        }

        return out;
    }

    CompositeAlgebra operator*(const Real auto other) const
    {
        CompositeAlgebra out;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_cn)
            {
                out.space.push_back(std::get<Lielab::domain::cn>(this->space[ii]) * other);
            }
            else if (ind == INDEX_gl)
            {
                out.space.push_back(std::get<Lielab::domain::gl>(this->space[ii]) * other);
            }
            else if (ind == INDEX_glc)
            {
                out.space.push_back(std::get<Lielab::domain::glc>(this->space[ii]) * other);
            }
            else if (ind == INDEX_rn)
            {
                out.space.push_back(std::get<Lielab::domain::rn>(this->space[ii]) * other);
            }
            else if (ind == INDEX_se)
            {
                out.space.push_back(std::get<Lielab::domain::se>(this->space[ii]) * other);
            }
            else if (ind == INDEX_so)
            {
                out.space.push_back(std::get<Lielab::domain::so>(this->space[ii]) * other);
            }
            else if (ind == INDEX_sp)
            {
                out.space.push_back(std::get<Lielab::domain::sp>(this->space[ii]) * other);
            }
            else if (ind == INDEX_su)
            {
                out.space.push_back(std::get<Lielab::domain::su>(this->space[ii]) * other);
            }
        }

        return out;
    }

    CompositeAlgebra & operator*=(const Real auto other)
    {
        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_cn)
            {
                std::get<Lielab::domain::cn>(this->space[ii]) *= other;
            }
            else if (ind == INDEX_gl)
            {
                std::get<Lielab::domain::gl>(this->space[ii]) *= other;
            }
            else if (ind == INDEX_glc)
            {
                std::get<Lielab::domain::glc>(this->space[ii]) *= other;
            }
            else if (ind == INDEX_rn)
            {
                std::get<Lielab::domain::rn>(this->space[ii]) *= other;
            }
            else if (ind == INDEX_se)
            {
                std::get<Lielab::domain::se>(this->space[ii]) *= other;
            }
            else if (ind == INDEX_so)
            {
                std::get<Lielab::domain::so>(this->space[ii]) *= other;
            }
            else if (ind == INDEX_sp)
            {
                std::get<Lielab::domain::sp>(this->space[ii]) *= other;
            }
            else if (ind == INDEX_su)
            {
                std::get<Lielab::domain::su>(this->space[ii]) *= other;
            }
        }

        return *this;
    }

    CompositeAlgebra operator*(const CompositeAlgebra & other) const
    {
        CompositeAlgebra out;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_cn)
            {
                out.space.push_back(std::get<Lielab::domain::cn>(this->space[ii]) * std::get<Lielab::domain::cn>(other.space[ii]));
            }
            else if (ind == INDEX_gl)
            {
                out.space.push_back(std::get<Lielab::domain::gl>(this->space[ii]) * std::get<Lielab::domain::gl>(other.space[ii]));
            }
            else if (ind == INDEX_glc)
            {
                out.space.push_back(std::get<Lielab::domain::glc>(this->space[ii]) * std::get<Lielab::domain::glc>(other.space[ii]));
            }
            else if (ind == INDEX_rn)
            {
                out.space.push_back(std::get<Lielab::domain::rn>(this->space[ii]) * std::get<Lielab::domain::rn>(other.space[ii]));
            }
            else if (ind == INDEX_se)
            {
                out.space.push_back(std::get<Lielab::domain::se>(this->space[ii]) * std::get<Lielab::domain::se>(other.space[ii]));
            }
            else if (ind == INDEX_so)
            {
                out.space.push_back(std::get<Lielab::domain::so>(this->space[ii]) * std::get<Lielab::domain::so>(other.space[ii]));
            }
            else if (ind == INDEX_sp)
            {
                out.space.push_back(std::get<Lielab::domain::sp>(this->space[ii]) * std::get<Lielab::domain::sp>(other.space[ii]));
            }
            else if (ind == INDEX_su)
            {
                out.space.push_back(std::get<Lielab::domain::su>(this->space[ii]) * std::get<Lielab::domain::su>(other.space[ii]));
            }
        }

        return out;
    }

    CompositeAlgebra & operator*=(const CompositeAlgebra & other)
    {
        assert(this->space.size() == other.space.size());

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_cn)
            {
                std::get<Lielab::domain::cn>(this->space[ii]) *= std::get<Lielab::domain::cn>(other.space[ii]);
            }
            if (ind == INDEX_gl)
            {
                std::get<Lielab::domain::gl>(this->space[ii]) *= std::get<Lielab::domain::gl>(other.space[ii]);
            }
            if (ind == INDEX_glc)
            {
                std::get<Lielab::domain::glc>(this->space[ii]) *= std::get<Lielab::domain::glc>(other.space[ii]);
            }
            else if (ind == INDEX_rn)
            {
                std::get<Lielab::domain::rn>(this->space[ii]) *= std::get<Lielab::domain::rn>(other.space[ii]);
            }
            else if (ind == INDEX_se)
            {
                std::get<Lielab::domain::se>(this->space[ii]) *= std::get<Lielab::domain::se>(other.space[ii]);
            }
            else if (ind == INDEX_so)
            {
                std::get<Lielab::domain::so>(this->space[ii]) *= std::get<Lielab::domain::so>(other.space[ii]);
            }
            else if (ind == INDEX_sp)
            {
                std::get<Lielab::domain::sp>(this->space[ii]) *= std::get<Lielab::domain::sp>(other.space[ii]);
            }
            else if (ind == INDEX_su)
            {
                std::get<Lielab::domain::su>(this->space[ii]) *= std::get<Lielab::domain::su>(other.space[ii]);
            }
        }

        return *this;
    }

    CompositeAlgebra operator/(const Real auto other) const
    {
        CompositeAlgebra out;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_cn)
            {
                out.space.push_back(std::get<Lielab::domain::cn>(this->space[ii]) / other);
            }
            else if (ind == INDEX_gl)
            {
                out.space.push_back(std::get<Lielab::domain::gl>(this->space[ii]) / other);
            }
            else if (ind == INDEX_glc)
            {
                out.space.push_back(std::get<Lielab::domain::glc>(this->space[ii]) / other);
            }
            else if (ind == INDEX_rn)
            {
                out.space.push_back(std::get<Lielab::domain::rn>(this->space[ii]) / other);
            }
            else if (ind == INDEX_se)
            {
                out.space.push_back(std::get<Lielab::domain::se>(this->space[ii]) / other);
            }
            else if (ind == INDEX_so)
            {
                out.space.push_back(std::get<Lielab::domain::so>(this->space[ii]) / other);
            }
            else if (ind == INDEX_sp)
            {
                out.space.push_back(std::get<Lielab::domain::sp>(this->space[ii]) / other);
            }
            else if (ind == INDEX_su)
            {
                out.space.push_back(std::get<Lielab::domain::su>(this->space[ii]) / other);
            }
        }

        return out;
    }

    CompositeAlgebra & operator/=(const Real auto other)
    {
        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_cn)
            {
                std::get<Lielab::domain::cn>(this->space[ii]) /= other;
            }
            else if (ind == INDEX_gl)
            {
                std::get<Lielab::domain::gl>(this->space[ii]) /= other;
            }
            else if (ind == INDEX_glc)
            {
                std::get<Lielab::domain::glc>(this->space[ii]) /= other;
            }
            else if (ind == INDEX_rn)
            {
                std::get<Lielab::domain::rn>(this->space[ii]) /= other;
            }
            else if (ind == INDEX_se)
            {
                std::get<Lielab::domain::se>(this->space[ii]) /= other;
            }
            else if (ind == INDEX_so)
            {
                std::get<Lielab::domain::so>(this->space[ii]) /= other;
            }
            else if (ind == INDEX_sp)
            {
                std::get<Lielab::domain::sp>(this->space[ii]) /= other;
            }
            else if (ind == INDEX_su)
            {
                std::get<Lielab::domain::su>(this->space[ii]) /= other;
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
            if (ind == INDEX_cn)
            {
                out += "cn(" + std::to_string(shapes(ii)) + ")";
            }
            else if (ind == INDEX_gl)
            {
                out += "gl(" + std::to_string(shapes(ii)) + ")";
            }
            else if (ind == INDEX_glc)
            {
                out += "glc(" + std::to_string(shapes(ii)) + ")";
            }
            else if (ind == INDEX_rn)
            {
                out += "rn(" + std::to_string(shapes(ii)) + ")";
            }
            else if (ind == INDEX_se)
            {
                out += "se(" + std::to_string(shapes(ii)) + ")";
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

    friend std::ostream & operator<<(std::ostream & os, const CompositeAlgebra & other);
};

std::ostream & operator<<(std::ostream & os, const CompositeAlgebra & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */

    os << other.to_string();
    return os;
}

CompositeAlgebra operator*(const Real auto other, const CompositeAlgebra & rhs)
{
    return rhs*other;
}
}
}

#endif
