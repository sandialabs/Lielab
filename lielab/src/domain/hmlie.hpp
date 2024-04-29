#ifndef _LIELAB_DOMAIN_HMLIE_HPP
#define _LIELAB_DOMAIN_HMLIE_HPP

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

class hmlie
{
    public:

    static constexpr size_t INDEX_GL = 0;
    static constexpr size_t INDEX_RN = 1;
    static constexpr size_t INDEX_SE = 2;
    static constexpr size_t INDEX_SO = 3;
    static constexpr size_t INDEX_SP = 4;
    static constexpr size_t INDEX_SU = 5;
    static constexpr size_t INDEX_gl = 6;
    static constexpr size_t INDEX_rn = 7;
    static constexpr size_t INDEX_se = 8;
    static constexpr size_t INDEX_so = 9;
    static constexpr size_t INDEX_sp = 10;
    static constexpr size_t INDEX_su = 11;

    typedef std::variant<lielab::domain::GL,
                         lielab::domain::RN,
                         lielab::domain::SE,
                         lielab::domain::SO,
                         lielab::domain::SP,
                         lielab::domain::SU,
                         lielab::domain::gl,
                         lielab::domain::rn,
                         lielab::domain::se,
                         lielab::domain::so,
                         lielab::domain::sp,
                         lielab::domain::su> TYPES;

    std::vector<TYPES> space;

    hmlie()
    {
        
    }

    hmlie(std::initializer_list<TYPES> others)
    {
        /*!
        * Initializer list constructor for hmlie
        *
        * Enables construction like:
        *     lielab::domain::RN R(3);
        *     lielab::domain::SO O(3);
        *     lielab::domain::hmlie M{R, O};
        */

        this->space = std::vector<TYPES>{std::move(others)};
    }

    hmlie(const std::vector<TYPES> & others)
    {
        /*!
        * Vector list constructor for hmlie
        *
        * Not needed for C++, but enables construction in Python like:
        *     R = lielab.domain.RN(3)
        *     O = lielab.domain.SO(3)
        *     M = lielab.domain.hmlie([R, O])
        */

        this->space = others;
    }

    Eigen::VectorXi get_shape() const
    {
        Eigen::VectorXi out(space.size());

        for (int ii = 0; ii < space.size(); ii++)
        {
            const size_t ind = space[ii].index();
            if (ind == INDEX_GL)
            {
                out(ii) = static_cast<int>(std::get<lielab::domain::GL>(space[ii]).shape);
            }
            else if (ind == INDEX_RN)
            {
                out(ii) = static_cast<int>(std::get<lielab::domain::RN>(space[ii]).shape);
            }
            else if (ind == INDEX_SE)
            {
                out(ii) = static_cast<int>(std::get<lielab::domain::SE>(space[ii]).shape);
            }
            else if (ind == INDEX_SO)
            {
                out(ii) = static_cast<int>(std::get<lielab::domain::SO>(space[ii]).shape);
            }
            else if (ind == INDEX_SP)
            {
                out(ii) = static_cast<int>(std::get<lielab::domain::SP>(space[ii]).shape);
            }
            else if (ind == INDEX_SU)
            {
                out(ii) = static_cast<int>(std::get<lielab::domain::SU>(space[ii]).shape);
            }
            else if (ind == INDEX_gl)
            {
                out(ii) = static_cast<int>(std::get<lielab::domain::gl>(space[ii]).shape);
            }
            else if (ind == INDEX_rn)
            {
                out(ii) = static_cast<int>(std::get<lielab::domain::rn>(space[ii]).shape);
            }
            else if (ind == INDEX_se)
            {
                out(ii) = static_cast<int>(std::get<lielab::domain::se>(space[ii]).shape);
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

            if (ind == INDEX_GL)
            {
                serials.push_back(std::get<lielab::domain::GL>(space[ii]).serialize());
            }
            else if (ind == INDEX_RN)
            {
                serials.push_back(std::get<lielab::domain::RN>(space[ii]).serialize());
            }
            else if (ind == INDEX_SE)
            {
                serials.push_back(std::get<lielab::domain::SE>(space[ii]).serialize());
            }
            else if (ind == INDEX_SO)
            {
                serials.push_back(std::get<lielab::domain::SO>(space[ii]).serialize());
            }
            else if (ind == INDEX_SP)
            {
                serials.push_back(std::get<lielab::domain::SP>(space[ii]).serialize());
            }
            else if (ind == INDEX_SU)
            {
                serials.push_back(std::get<lielab::domain::SU>(space[ii]).serialize());
            }
            else if (ind == INDEX_gl)
            {
                serials.push_back(std::get<lielab::domain::gl>(space[ii]).get_vector());
            }
            else if (ind == INDEX_rn)
            {
                serials.push_back(std::get<lielab::domain::rn>(space[ii]).get_vector());
            }
            else if (ind == INDEX_se)
            {
                serials.push_back(std::get<lielab::domain::se>(space[ii]).get_vector());
            }
            else if (ind == INDEX_so)
            {
                serials.push_back(std::get<lielab::domain::so>(space[ii]).get_vector());
            }
            else if (ind == INDEX_sp)
            {
                serials.push_back(std::get<lielab::domain::sp>(space[ii]).get_vector());
            }
            else if (ind == INDEX_su)
            {
                serials.push_back(std::get<lielab::domain::su>(space[ii]).get_vector());
            }
        }
        
        Eigen::VectorXd out = lielab::utils::concatenate(serials);
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

            if (ind == INDEX_GL)
            {
                sz = std::get<lielab::domain::GL>(M).get_size();
                std::get<lielab::domain::GL>(M).unserialize(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_RN)
            {
                sz = std::get<lielab::domain::RN>(M).get_size();
                std::get<lielab::domain::RN>(M).unserialize(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_SE)
            {
                sz = std::get<lielab::domain::SE>(M).get_size();
                std::get<lielab::domain::SE>(M).unserialize(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_SO)
            {
                sz = std::get<lielab::domain::SO>(M).get_size();
                std::get<lielab::domain::SO>(M).unserialize(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_SP)
            {
                sz = std::get<lielab::domain::SP>(M).get_size();
                std::get<lielab::domain::SP>(M).unserialize(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_SU)
            {
                sz = std::get<lielab::domain::SU>(M).get_size();
                std::get<lielab::domain::SU>(M).unserialize(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_gl)
            {
                sz = std::get<lielab::domain::gl>(M).get_dimension();
                std::get<lielab::domain::gl>(M).set_vector(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_rn)
            {
                sz = std::get<lielab::domain::rn>(M).get_dimension();
                std::get<lielab::domain::rn>(M).set_vector(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_se)
            {
                sz = std::get<lielab::domain::se>(M).get_dimension();
                std::get<lielab::domain::se>(M).set_vector(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_so)
            {
                sz = std::get<lielab::domain::so>(M).get_dimension();
                std::get<lielab::domain::so>(M).set_vector(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_sp)
            {
                sz = std::get<lielab::domain::sp>(M).get_dimension();
                std::get<lielab::domain::sp>(M).set_vector(vec(Eigen::seqN(jj, sz)));
            }
            else if (ind == INDEX_su)
            {
                sz = std::get<lielab::domain::su>(M).get_dimension();
                std::get<lielab::domain::su>(M).set_vector(vec(Eigen::seqN(jj, sz)));
            }
        }
    }

    hmlie operator*(const hmlie & other) const
    {
        hmlie out;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_GL)
            {
                out.space.push_back(std::get<lielab::domain::GL>(this->space[ii]) * std::get<lielab::domain::GL>(other.space[ii]));
            }
            else if (ind == INDEX_RN)
            {
                out.space.push_back(std::get<lielab::domain::RN>(this->space[ii]) * std::get<lielab::domain::RN>(other.space[ii]));
            }
            else if (ind == INDEX_SE)
            {
                out.space.push_back(std::get<lielab::domain::SE>(this->space[ii]) * std::get<lielab::domain::SE>(other.space[ii]));
            }
            else if (ind == INDEX_SO)
            {
                out.space.push_back(std::get<lielab::domain::SO>(this->space[ii]) * std::get<lielab::domain::SO>(other.space[ii]));
            }
            else if (ind == INDEX_SP)
            {
                out.space.push_back(std::get<lielab::domain::SP>(this->space[ii]) * std::get<lielab::domain::SP>(other.space[ii]));
            }
            else if (ind == INDEX_SU)
            {
                out.space.push_back(std::get<lielab::domain::SU>(this->space[ii]) * std::get<lielab::domain::SU>(other.space[ii]));
            }
            else if (ind == INDEX_gl)
            {
                out.space.push_back(std::get<lielab::domain::gl>(this->space[ii]) + std::get<lielab::domain::gl>(other.space[ii]));
            }
            else if (ind == INDEX_rn)
            {
                out.space.push_back(std::get<lielab::domain::rn>(this->space[ii]) + std::get<lielab::domain::rn>(other.space[ii]));
            }
            else if (ind == INDEX_se)
            {
                out.space.push_back(std::get<lielab::domain::se>(this->space[ii]) + std::get<lielab::domain::se>(other.space[ii]));
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

    hmlie & operator*=(const hmlie & other)
    {
        assert(this->space.size() == other.space.size());

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_GL)
            {
                std::get<lielab::domain::GL>(this->space[ii]) *= std::get<lielab::domain::GL>(other.space[ii]);
            }
            else if (ind == INDEX_RN)
            {
                std::get<lielab::domain::RN>(this->space[ii]) *= std::get<lielab::domain::RN>(other.space[ii]);
            }
            else if (ind == INDEX_SE)
            {
                std::get<lielab::domain::SE>(this->space[ii]) *= std::get<lielab::domain::SE>(other.space[ii]);
            }
            else if (ind == INDEX_SO)
            {
                std::get<lielab::domain::SO>(this->space[ii]) *= std::get<lielab::domain::SO>(other.space[ii]);
            }
            else if (ind == INDEX_SP)
            {
                std::get<lielab::domain::SP>(this->space[ii]) *= std::get<lielab::domain::SP>(other.space[ii]);
            }
            else if (ind == INDEX_SU)
            {
                std::get<lielab::domain::SU>(this->space[ii]) *= std::get<lielab::domain::SU>(other.space[ii]);
            }
            else if (ind == INDEX_gl)
            {
                std::get<lielab::domain::gl>(this->space[ii]) += std::get<lielab::domain::gl>(other.space[ii]);
            }
            else if (ind == INDEX_rn)
            {
                std::get<lielab::domain::rn>(this->space[ii]) += std::get<lielab::domain::rn>(other.space[ii]);
            }
            else if (ind == INDEX_se)
            {
                std::get<lielab::domain::se>(this->space[ii]) += std::get<lielab::domain::se>(other.space[ii]);
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

    hmlie inverse() const
    {
        hmlie out;

        for (int ii = 0; ii < this->space.size(); ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_GL)
            {
                out.space.push_back(std::get<lielab::domain::GL>(this->space[ii]).inverse());
            }
            else if (ind == INDEX_RN)
            {
                out.space.push_back(std::get<lielab::domain::RN>(this->space[ii]).inverse());
            }
            else if (ind == INDEX_SE)
            {
                out.space.push_back(std::get<lielab::domain::SE>(this->space[ii]).inverse());
            }
            else if (ind == INDEX_SO)
            {
                out.space.push_back(std::get<lielab::domain::SO>(this->space[ii]).inverse());
            }
            else if (ind == INDEX_SP)
            {
                out.space.push_back(std::get<lielab::domain::SP>(this->space[ii]).inverse());
            }
            else if (ind == INDEX_SU)
            {
                out.space.push_back(std::get<lielab::domain::SU>(this->space[ii]).inverse());
            }
            else if (ind == INDEX_gl)
            {
                out.space.push_back(-std::get<lielab::domain::gl>(this->space[ii]));
            }
            else if (ind == INDEX_rn)
            {
                out.space.push_back(-std::get<lielab::domain::rn>(this->space[ii]));
            }
            else if (ind == INDEX_se)
            {
                out.space.push_back(-std::get<lielab::domain::se>(this->space[ii]));
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

    std::string to_string() const
    {
        std::string out = "";
        const size_t sz = this->space.size();
        Eigen::VectorXi shapes = this->get_shape();

        for (int ii = 0; ii < sz; ii++)
        {
            const size_t ind = this->space[ii].index();
            if (ind == INDEX_GL)
            {
                out += "GL(" + std::to_string(shapes(ii)) + ")";
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
            else if (ind == INDEX_gl)
            {
                out += "gl(" + std::to_string(shapes(ii)) + ")";
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

    friend std::ostream & operator<<(std::ostream & os, const hmlie & other);
};

std::ostream & operator<<(std::ostream & os, const hmlie & other)
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
