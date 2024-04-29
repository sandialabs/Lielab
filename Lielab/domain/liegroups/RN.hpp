#ifndef _LIELAB_DOMAIN_RN_HPP
#define _LIELAB_DOMAIN_RN_HPP

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab
{
    namespace domain
    {
        class RN
        {
            /*!
            * The RN class.
            */
            public:
            static constexpr bool abelian = true;
            size_t shape = 0;

            typedef Eigen::VectorXd Base;
            Base _data;

            // Initialization methods

            RN() : RN(0)
            {
                /*! \f{equation*}{ () \rightarrow RN \f}
                * 
                * Empty initialization function. Enables instantiation like:
                * 
                *     Lielab::domain::RN x, y, z;
                * 
                */

            }

            RN(const size_t shape)
            {
                /*! \f{equation*}{ (\mathbb{Z}) \rightarrow RN \f}
                *
                * Constructor instantiating an \f$RN\f$ object.
                * 
                * Enables instantiation like:
                * 
                *     Lielab::domain::RN x(2), y(3), z(4);
                * 
                * @param[in] shape The shape of the data matrix.
                */
                
                this->shape = shape;
                const size_t L = (shape == 0) ? 0 : shape - 1;
                this->_data = Eigen::VectorXd::Zero(L);
            }

            RN(std::initializer_list<double> other)
            {
                /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow RN \f}
                *
                * Constructor instantiating an \f$RN\f$ object from an \f$n \times 1\f$
                * real vector.
                *
                * @param[in] other The object to instantiate from as a real matrix.
                */

                this->shape = other.size() + 1;
                this->_data = Eigen::VectorXd{std::move(other)};
            }

            template<typename OtherDerived>
            RN(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& RN \\ (\mathbb{R}^{n \times 1}) &\rightarrow& RN \f}
                *
                * Constructor instantiating an \f$RN\f$ object from either an
                * \f$n \times n\f$ real matrix or \f$n \times 1\f$ real vector.
                *
                * @param[in] other The object to instantiate from as a real matrix.
                */

                if (other.rows() == other.cols() && other.rows() != 1)
                {
                    this->shape = other.rows();
                    this->_data = Eigen::VectorXd::Zero(shape-1);
                    for (int ii = 0; ii < shape-1; ii++)
                    {
                        _data(ii) = other(ii, shape-1);
                    }
                }
                else if (other.cols() == 1)
                {
                    this->shape = other.rows()+1;
                    this->_data = other;
                }
                else
                {
                    throw Errorx("Input data to RN malformed.");
                }
            }

            template<typename OtherDerived>
            RN & operator=(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{eqnarray*}{ RN &:= \mathbb{R}^{n \times n} \\ &:= \mathbb{R}^{n \times 1} \f}
                * 
                * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `RN`.
                */

                if (other.rows() == other.cols() && other.rows() != 1)
                {   
                    this->shape = other.rows();
                    this->_data = Eigen::VectorXd::Zero(shape-1);
                    for (int ii = 0; ii < other.rows()-1; ii++)
                    {
                        _data(ii) = other(ii, shape-1);
                    }
                }
                else if (other.cols() == 1)
                {
                    shape = other.rows()+1;
                    _data = other;
                }
                return *this;
            }

            static Eigen::MatrixXd project(const Eigen::MatrixXd & other)
            {
                /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in RN \f}
                *
                * Projects a matrix suitable for data.
                */

                if (other.rows() != other.cols())
                {
                    throw Errorx("Size of the matrix must be square.");
                }

                Eigen::MatrixXd out = other;

                for (int ii = 0; ii < out.rows(); ii++)
                {
                    for (int jj = 0; jj < out.cols(); jj++)
                    {
                        if (jj != other.cols()-1)
                        {
                            out(ii,jj) = 0.0;
                        }
                    }
                }

                out += Eigen::MatrixXd::Identity(other.rows(), other.cols());

                out(out.rows()-1, out.cols()-1) = 0.0;
                return out;
            }

            size_t get_dimension() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
                * 
                * Gets the dimension of the group.
                */

                return static_cast<size_t>(this->shape - 1);
            }

            size_t get_size() const
            {
                /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
                 *
                 * Gets the size of the data representation.
                 */

                return static_cast<size_t>(this->shape - 1);
            }

            Eigen::MatrixXd get_matrix() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times n} \f}
                * 
                * Returns a matrix representation.
                * 
                * Formerly called "get_ados_representation()".
                * 
                * Ado, Igor D. "Note on the representation of finite continuous groups by
                *               means of linear substitutions, Izv. Fiz." Mat. Obsch.(Kazan)
                *               7.1 (1935): 935.
                * 
                * Ado, Igor D. "The representation of Lie algebras by matrices." Uspekhi
                *               Matematicheskikh Nauk 2.6 (1947): 159-173.
                */

                Eigen::MatrixXd out = Eigen::MatrixXd::Identity(shape, shape);

                for (int ii = 0; ii < shape-1; ii++)
                {
                    out(ii, shape-1) = _data(ii);
                }

                return out;
            }

            RN inverse() const
            {
                /*! \f{equation*}{ (RN) \rightarrow RN \f}
                * 
                * Returns the inverse.
                */

                return -this->_data;
            }

            // Data representation

            Eigen::VectorXd serialize() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
                * 
                * Returns a serialized representation.
                */

                return _data;
            }

            void unserialize(const Eigen::VectorXd &vec)
            {
                /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow () \f}
                * 
                * Sets the RN object from a serialized vector.
                */

                this->_data = vec(Eigen::seqN(0, this->shape-1));
            }

            double operator()(const size_t index) const
            {
                /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{R} \f}
                *
                * Gets a value in the column vector representation.
                */

                return _data(index);
            }

            double & operator()(const size_t index)
            {
                /*! \f{equation*}{ RN(\mathbb{Z}) := \mathbb{R} \f}
                *
                * Assignment of a value in the column vector representation.
                */

                return _data(index);
            }

            double operator()(const size_t index1, const size_t index2) const
            {
                /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
                *
                * Gets a value in the square matrix representation.
                */

                Eigen::MatrixXd mat = this->get_matrix();
                return mat(index1, index2);
            }

            size_t rows() const
            {
                // TODO: Delete this method
                return shape;
            }

            size_t cols() const
            {
                // TODO: Delete this method
                return shape;
            }

            RN operator*(const RN & other) const
            {
                /*! \f{equation*}{ (RN, RN) \rightarrow RN \f}
                *
                * Group product.
                */

                assert(this->shape == other.shape);
                Base out = _data + other._data;
                return out;
            }

            RN & operator*=(const RN & other)
            {
                /*! \f{equation*}{ (RN, RN) \rightarrow RN \f}
                *
                * In place group product.
                */

                assert(this->shape == other.shape);
                this->_data += other._data;
                return *this;
            }

            friend std::ostream & operator<<(std::ostream& os, const RN & other);
        };
        
        std::ostream & operator<<(std::ostream& os, const RN & other)
        {
            /*!
            * Overloads the "<<" stream insertion operator.
            */

            os << static_cast<const Eigen::VectorXd>(other._data);
            return os;
        }
    }
}

#endif
