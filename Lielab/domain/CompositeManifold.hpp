#ifndef _LIELAB_DOMAIN_COMPOSITEMANIFOLD_HPP
#define _LIELAB_DOMAIN_COMPOSITEMANIFOLD_HPP

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

// lielab dependencies
#include "constants.hpp"

using Lielab::abstract::Real;
using Lielab::abstract::Imaginary;

namespace Lielab
{
namespace domain
{

class CompositeManifold
{
    public:

    static constexpr size_t INDEX_CN  = 0;
    static constexpr size_t INDEX_GL  = 1;
    static constexpr size_t INDEX_GLC = 2;
    static constexpr size_t INDEX_RN  = 3;
    static constexpr size_t INDEX_SE  = 4;
    static constexpr size_t INDEX_SO  = 5;
    static constexpr size_t INDEX_SP  = 6;
    static constexpr size_t INDEX_SU  = 7;
    static constexpr size_t INDEX_cn  = 8;
    static constexpr size_t INDEX_gl  = 9;
    static constexpr size_t INDEX_glc = 10;
    static constexpr size_t INDEX_rn  = 11;
    static constexpr size_t INDEX_se  = 12;
    static constexpr size_t INDEX_so  = 13;
    static constexpr size_t INDEX_sp  = 14;
    static constexpr size_t INDEX_su  = 15;

    typedef std::variant<Lielab::domain::CN,
                         Lielab::domain::GL,
                         Lielab::domain::GLC,
                         Lielab::domain::RN,
                         Lielab::domain::SE,
                         Lielab::domain::SO,
                         Lielab::domain::SP,
                         Lielab::domain::SU,
                         Lielab::domain::cn,
                         Lielab::domain::gl,
                         Lielab::domain::glc,
                         Lielab::domain::rn,
                         Lielab::domain::se,
                         Lielab::domain::so,
                         Lielab::domain::sp,
                         Lielab::domain::su> TYPES;

    std::vector<TYPES> space;

    CompositeManifold()
    {
        
    }

    CompositeManifold(std::initializer_list<TYPES> others)
    {
        /*!
        * Initializer list constructor for CompositeManifold
        *
        * Enables construction like:
        *     Lielab::domain::RN R(3);
        *     Lielab::domain::SO O(3);
        *     Lielab::domain::CompositeManifold M{R, O};
        */

        this->space = std::vector<TYPES>{std::move(others)};
    }

    CompositeManifold(const std::vector<TYPES> & others)
    {
        /*!
        * Vector list constructor for CompositeManifold
        *
        * Not needed for C++, but enables construction in Python like:
        *     R = lielab.domain.RN(3)
        *     O = lielab.domain.SO(3)
        *     M = lielab.domain.CompositeManifold([R, O])
        */

        this->space = others;
    }

    Eigen::VectorXi get_shape() const
    {
        Eigen::VectorXi out(space.size());

        for (int ii = 0; ii < space.size(); ii++)
        {
            const size_t ind = space[ii].index();
            if (ind == INDEX_CN)
            {
                out(ii) = static_cast<int>(std::get<Lielab::domain::CN>(space[ii]).shape);
            }
            else if (ind == INDEX_GL)
            {
                out(ii) = static_cast<int>(std::get<Lielab::domain::GL>(space[ii]).shape);
            }
            else if (ind == INDEX_GLC)
            {
                out(ii) = static_cast<int>(std::get<Lielab::domain::GLC>(space[ii]).shape);
            }
            else if (ind == INDEX_RN)
            {
                out(ii) = static_cast<int>(std::get<Lielab::domain::RN>(space[ii]).shape);
            }
            else if (ind == INDEX_SE)
            {
                out(ii) = static_cast<int>(std::get<Lielab::domain::SE>(space[ii]).shape);
            }
            else if (ind == INDEX_SO)
            {
                out(ii) = static_cast<int>(std::get<Lielab::domain::SO>(space[ii]).shape);
            }
            else if (ind == INDEX_SP)
            {
                out(ii) = static_cast<int>(std::get<Lielab::domain::SP>(space[ii]).shape);
            }
            else if (ind == INDEX_SU)
            {
                out(ii) = static_cast<int>(std::get<Lielab::domain::SU>(space[ii]).shape);
            }
            else if (ind == INDEX_cn)
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

    Eigen::VectorXd serialize() const
    {

        if (space.size() == 0)
        {
            return Eigen::VectorXd(0);
        }
        
        std::vector<Eigen::VectorXd> serials;
        const size_t sz = this->space.size();

        for (int ii = 0; ii < sz; ii++)
        {
            const size_t ind = space[ii].index();

            if (ind == INDEX_CN)
            {
                serials.push_back(std::get<Lielab::domain::CN>(space[ii]).serialize());
            }
            else if (ind == INDEX_GL)
            {
                serials.push_back(std::get<Lielab::domain::GL>(space[ii]).serialize());
            }
            else if (ind == INDEX_GLC)
            {
                serials.push_back(std::get<Lielab::domain::GLC>(space[ii]).serialize());
            }
            else if (ind == INDEX_RN)
            {
                serials.push_back(std::get<Lielab::domain::RN>(space[ii]).serialize());
            }
            else if (ind == INDEX_SE)
            {
                serials.push_back(std::get<Lielab::domain::SE>(space[ii]).serialize());
            }
            else if (ind == INDEX_SO)
            {
                serials.push_back(std::get<Lielab::domain::SO>(space[ii]).serialize());
            }
            else if (ind == INDEX_SP)
            {
                serials.push_back(std::get<Lielab::domain::SP>(space[ii]).serialize());
            }
            else if (ind == INDEX_SU)
            {
                serials.push_back(std::get<Lielab::domain::SU>(space[ii]).serialize());
            }
            else if (ind == INDEX_cn)
            {
                serials.push_back(std::get<Lielab::domain::cn>(space[ii]).get_vector());
            }
            else if (ind == INDEX_gl)
            {
                serials.push_back(std::get<Lielab::domain::gl>(space[ii]).get_vector());
            }
            else if (ind == INDEX_glc)
            {
                serials.push_back(std::get<Lielab::domain::glc>(space[ii]).get_vector());
            }
            else if (ind == INDEX_rn)
            {
                serials.push_back(std::get<Lielab::domain::rn>(space[ii]).get_vector());
            }
            else if (ind == INDEX_se)
            {
                serials.push_back(std::get<Lielab::domain::se>(space[ii]).get_vector());
            }
            else if (ind == INDEX_so)
            {
                serials.push_back(std::get<Lielab::domain::so>(space[ii]).get_vector());
            }
            else if (ind == INDEX_sp)
            {
                serials.push_back(std::get<Lielab::domain::sp>(space[ii]).get_vector());
            }
            else if (ind == INDEX_su)
            {
                serials.push_back(std::get<Lielab::domain::su>(space[ii]).get_vector());
            }
        }
        
        Eigen::VectorXd out = Lielab::utils::concatenate(serials);
        return out;
    }

    void unserialize(const Eigen::VectorXd &vec)
    {
        ptrdiff_t jj = 0;
        ptrdiff_t sz = 0;

        for (auto & M : this->space)
        {
            const size_t ind = M.index();
            jj += sz;

            if (ind == INDEX_CN)
            {
                sz = std::get<Lielab::domain::CN>(M).get_size();
                std::get<Lielab::domain::CN>(M).unserialize(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_GL)
            {
                sz = std::get<Lielab::domain::GL>(M).get_size();
                std::get<Lielab::domain::GL>(M).unserialize(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_GLC)
            {
                sz = std::get<Lielab::domain::GLC>(M).get_size();
                std::get<Lielab::domain::GLC>(M).unserialize(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_RN)
            {
                sz = std::get<Lielab::domain::RN>(M).get_size();
                std::get<Lielab::domain::RN>(M).unserialize(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_SE)
            {
                sz = std::get<Lielab::domain::SE>(M).get_size();
                std::get<Lielab::domain::SE>(M).unserialize(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_SO)
            {
                sz = std::get<Lielab::domain::SO>(M).get_size();
                std::get<Lielab::domain::SO>(M).unserialize(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_SP)
            {
                sz = std::get<Lielab::domain::SP>(M).get_size();
                std::get<Lielab::domain::SP>(M).unserialize(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_SU)
            {
                sz = std::get<Lielab::domain::SU>(M).get_size();
                std::get<Lielab::domain::SU>(M).unserialize(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_cn)
            {
                sz = std::get<Lielab::domain::cn>(M).get_dimension();
                std::get<Lielab::domain::cn>(M).set_vector(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_gl)
            {
                sz = std::get<Lielab::domain::gl>(M).get_dimension();
                std::get<Lielab::domain::gl>(M).set_vector(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_glc)
            {
                sz = std::get<Lielab::domain::glc>(M).get_dimension();
                std::get<Lielab::domain::glc>(M).set_vector(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_rn)
            {
                sz = std::get<Lielab::domain::rn>(M).get_dimension();
                std::get<Lielab::domain::rn>(M).set_vector(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_se)
            {
                sz = std::get<Lielab::domain::se>(M).get_dimension();
                std::get<Lielab::domain::se>(M).set_vector(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_so)
            {
                sz = std::get<Lielab::domain::so>(M).get_dimension();
                std::get<Lielab::domain::so>(M).set_vector(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_sp)
            {
                sz = std::get<Lielab::domain::sp>(M).get_dimension();
                std::get<Lielab::domain::sp>(M).set_vector(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_su)
            {
                sz = std::get<Lielab::domain::su>(M).get_dimension();
                std::get<Lielab::domain::su>(M).set_vector(vec(Eigen::seqN(jj, sz)));
            }
        }
    }

    CompositeManifold operator*(const CompositeManifold & other) const
    {
        CompositeManifold out;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_CN)
            {
                out.space.push_back(std::get<Lielab::domain::CN>(this->space[ii]) * std::get<Lielab::domain::CN>(other.space[ii]));
            }
            else if (ind == INDEX_GL)
            {
                out.space.push_back(std::get<Lielab::domain::GL>(this->space[ii]) * std::get<Lielab::domain::GL>(other.space[ii]));
            }
            else if (ind == INDEX_GLC)
            {
                out.space.push_back(std::get<Lielab::domain::GLC>(this->space[ii]) * std::get<Lielab::domain::GLC>(other.space[ii]));
            }
            else if (ind == INDEX_RN)
            {
                out.space.push_back(std::get<Lielab::domain::RN>(this->space[ii]) * std::get<Lielab::domain::RN>(other.space[ii]));
            }
            else if (ind == INDEX_SE)
            {
                out.space.push_back(std::get<Lielab::domain::SE>(this->space[ii]) * std::get<Lielab::domain::SE>(other.space[ii]));
            }
            else if (ind == INDEX_SO)
            {
                out.space.push_back(std::get<Lielab::domain::SO>(this->space[ii]) * std::get<Lielab::domain::SO>(other.space[ii]));
            }
            else if (ind == INDEX_SP)
            {
                out.space.push_back(std::get<Lielab::domain::SP>(this->space[ii]) * std::get<Lielab::domain::SP>(other.space[ii]));
            }
            else if (ind == INDEX_SU)
            {
                out.space.push_back(std::get<Lielab::domain::SU>(this->space[ii]) * std::get<Lielab::domain::SU>(other.space[ii]));
            }
            else if (ind == INDEX_cn)
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

    CompositeManifold & operator*=(const CompositeManifold & other)
    {
        assert(this->space.size() == other.space.size());

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_CN)
            {
                std::get<Lielab::domain::CN>(this->space[ii]) *= std::get<Lielab::domain::CN>(other.space[ii]);
            }
            else if (ind == INDEX_GL)
            {
                std::get<Lielab::domain::GL>(this->space[ii]) *= std::get<Lielab::domain::GL>(other.space[ii]);
            }
            else if (ind == INDEX_GLC)
            {
                std::get<Lielab::domain::GLC>(this->space[ii]) *= std::get<Lielab::domain::GLC>(other.space[ii]);
            }
            else if (ind == INDEX_RN)
            {
                std::get<Lielab::domain::RN>(this->space[ii]) *= std::get<Lielab::domain::RN>(other.space[ii]);
            }
            else if (ind == INDEX_SE)
            {
                std::get<Lielab::domain::SE>(this->space[ii]) *= std::get<Lielab::domain::SE>(other.space[ii]);
            }
            else if (ind == INDEX_SO)
            {
                std::get<Lielab::domain::SO>(this->space[ii]) *= std::get<Lielab::domain::SO>(other.space[ii]);
            }
            else if (ind == INDEX_SP)
            {
                std::get<Lielab::domain::SP>(this->space[ii]) *= std::get<Lielab::domain::SP>(other.space[ii]);
            }
            else if (ind == INDEX_SU)
            {
                std::get<Lielab::domain::SU>(this->space[ii]) *= std::get<Lielab::domain::SU>(other.space[ii]);
            }
            else if (ind == INDEX_cn)
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

    CompositeManifold inverse() const
    {
        CompositeManifold out;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_CN)
            {
                out.space.push_back(std::get<Lielab::domain::CN>(this->space[ii]).inverse());
            }
            if (ind == INDEX_GL)
            {
                out.space.push_back(std::get<Lielab::domain::GL>(this->space[ii]).inverse());
            }
            if (ind == INDEX_GLC)
            {
                out.space.push_back(std::get<Lielab::domain::GLC>(this->space[ii]).inverse());
            }
            else if (ind == INDEX_RN)
            {
                out.space.push_back(std::get<Lielab::domain::RN>(this->space[ii]).inverse());
            }
            else if (ind == INDEX_SE)
            {
                out.space.push_back(std::get<Lielab::domain::SE>(this->space[ii]).inverse());
            }
            else if (ind == INDEX_SO)
            {
                out.space.push_back(std::get<Lielab::domain::SO>(this->space[ii]).inverse());
            }
            else if (ind == INDEX_SP)
            {
                out.space.push_back(std::get<Lielab::domain::SP>(this->space[ii]).inverse());
            }
            else if (ind == INDEX_SU)
            {
                out.space.push_back(std::get<Lielab::domain::SU>(this->space[ii]).inverse());
            }
            else if (ind == INDEX_cn)
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

    std::string to_string() const
    {
        std::string out = "";
        const size_t sz = this->space.size();
        Eigen::VectorXi shapes = this->get_shape();

        for (int ii = 0; ii < sz; ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_CN)
            {
                out += "CN(" + std::to_string(shapes(ii)) + ")";
            }
            else if (ind == INDEX_GL)
            {
                out += "GL(" + std::to_string(shapes(ii)) + ")";
            }
            else if (ind == INDEX_GLC)
            {
                out += "GLC(" + std::to_string(shapes(ii)) + ")";
            }
            else if (ind == INDEX_RN)
            {
                out += "RN(" + std::to_string(shapes(ii)) + ")";
            }
            else if (ind == INDEX_SE)
            {
                out += "SE(" + std::to_string(shapes(ii)) + ")";
            }
            else if (ind == INDEX_SO)
            {
                out += "SO(" + std::to_string(shapes(ii)) + ")";
            }
            else if (ind == INDEX_SP)
            {
                out += "SP(" + std::to_string(shapes(ii)) + ")";
            }
            else if (ind == INDEX_SU)
            {
                out += "SU(" + std::to_string(shapes(ii)) + ")";
            }
            else if (ind == INDEX_cn)
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

    friend std::ostream & operator<<(std::ostream & os, const CompositeManifold & other);
};

std::ostream & operator<<(std::ostream & os, const CompositeManifold & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */

    os << other.to_string();
    return os;
}

}
}

#endif
