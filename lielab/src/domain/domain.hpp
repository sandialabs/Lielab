#ifndef _LIELAB_DOMAIN_HPP
#define _LIELAB_DOMAIN_HPP

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

struct IndexError : public std::exception
{
   std::string s;
   IndexError(std::string ss) : s(ss) {}
   ~IndexError() throw () {}
   const char* what() const throw() { return s.c_str(); }
};

struct NotImplementedError : public std::exception
{
   std::string s;
   NotImplementedError(std::string ss) : s(ss) {}
   ~NotImplementedError() throw () {}
   const char* what() const throw() { return s.c_str(); }
};

struct SizeError : public std::exception
{
   std::string s;
   SizeError(std::string ss) : s(ss) {}
   ~SizeError() throw () {}
   const char* what() const throw() { return s.c_str(); }
};

namespace lielab
{
    /*!
    * The domain namespace. Houses all functions and classes related to domains available in lielab.
    * 
    * The naming convention used within the domain namespace is different than the standard naming
    * convention used elsewhere. Typically, classes are in PascalCase and functions are in lowercase.
    * Somewhat deviating from this standard, lie algebras are in lowercase and Lie Groups are in
    * UpperCase. This distinction is made to follow mathematical convections. For example, 
    * su is the special unitary lie algebra and SU is the associated Special Unitary group.
    */

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

        class hmlie
        {
            public:

            static constexpr size_t INDEX_GL = 0;
            static constexpr size_t INDEX_RN = 1;
            static constexpr size_t INDEX_SO = 2;
            static constexpr size_t INDEX_SP = 3;
            static constexpr size_t INDEX_SU = 4;

            typedef std::variant<lielab::domain::GL,
                                 lielab::domain::RN,
                                 lielab::domain::SO,
                                 lielab::domain::SP,
                                 lielab::domain::SU> TYPES;

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
                }
                
                Eigen::VectorXd out = lielab::utils::concatenate(serials);
                return out;
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
