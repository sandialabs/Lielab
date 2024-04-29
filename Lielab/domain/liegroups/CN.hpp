#ifndef _LIELAB_DOMAIN_CN_HPP
#define _LIELAB_DOMAIN_CN_HPP

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab
{
    namespace domain
    {
        class CN
        {
            /*!
            * The CN class.
            */
            public:
            static constexpr bool abelian = true;
            size_t shape = 0;

            typedef Eigen::VectorXcd Base;
            Base _data;

            // Initialization methods

            CN() : CN(0)
            {
                /*! \f{equation*}{ () \rightarrow CN \f}
                * 
                * Empty initialization function. Enables instantiation like:
                * 
                *     Lielab::domain::CN x, y, z;
                * 
                */

            }

            CN(const size_t shape)
            {
                /*! \f{equation*}{ (\mathbb{Z}) \rightarrow CN \f}
                *
                * Constructor instantiating an \f$CN\f$ object.
                * 
                * Enables instantiation like:
                * 
                *     Lielab::domain::CN x(2), y(3), z(4);
                * 
                * @param[in] shape The shape of the data matrix.
                */
                
                this->shape = shape;
                const size_t L = (shape == 0) ? 0 : shape - 1;
                this->_data = Eigen::VectorXcd::Zero(L);
            }

            CN(std::initializer_list<std::complex<double>> other)
            {
                /*! \f{equation*}{ (\mathbb{C}^{n \times 1}) \rightarrow CN \f}
                *
                * Constructor instantiating an \f$CN\f$ object from an \f$n \times 1\f$
                * imaginary vector.
                *
                * @param[in] other The object to instantiate from as an imaginary matrix.
                */

                this->_data = Eigen::VectorXcd{std::move(other)};
                this->shape = this->_data.size() + 1;
            }

            template<typename OtherDerived>
            CN(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{eqnarray*}{(\mathbb{C}^{n \times n}) &\rightarrow& CN \\ (\mathbb{C}^{n \times 1}) &\rightarrow& CN \f}
                *
                * Constructor instantiating an \f$CN\f$ object from either an
                * \f$n \times n\f$ imaginary matrix or \f$n \times 1\f$ imaginary vector.
                *
                * @param[in] other The object to instantiate from as a real matrix.
                */

                if (other.rows() == other.cols() && other.rows() != 1)
                {
                    this->shape = other.rows();
                    this->_data = Eigen::VectorXcd::Zero(shape-1);
                    for (int ii = 0; ii < shape-1; ii++)
                    {
                        _data(ii) = other(ii, shape-1);
                    }
                }
                else if (other.cols() == 1)
                {
                    this->_data = other;
                    this->shape = other.rows()+1;
                }
                else
                {
                    throw Errorx("Input data to CN malformed.");
                }
            }

            template<typename OtherDerived>
            CN & operator=(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{eqnarray*}{ CN &:= \mathbb{R}^{n \times n} \\ &:= \mathbb{C}^{n \times 1} \f}
                * 
                * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `CN`.
                */

                if (other.rows() == other.cols() && other.rows() != 1)
                {   
                    this->shape = other.rows();
                    this->_data = Eigen::VectorXcd::Zero(shape-1);
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

            static Eigen::MatrixXcd project(const Eigen::MatrixXcd & other)
            {
                /*! \f{equation*}{ (\mathbb{C}^{n \times n}) \rightarrow \mathbb{C}^{n \times n} \in CN \f}
                *
                * Projects a matrix suitable for data.
                */

                if (other.rows() != other.cols())
                {
                    throw Errorx("Size of the matrix must be square.");
                }

                Eigen::MatrixXcd out = other;

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

                out += Eigen::MatrixXcd::Identity(other.rows(), other.cols());

                out(out.rows()-1, out.cols()-1) = 0.0;
                return out;
            }

            size_t get_dimension() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
                * 
                * Gets the dimension of the group.
                */

                return static_cast<size_t>(2*(this->shape - 1));
            }

            size_t get_size() const
            {
                /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
                 *
                 * Gets the size of the data representation.
                 */

                return static_cast<size_t>(2*(this->shape - 1));
            }

            Eigen::MatrixXcd get_ados_representation() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times n} \f}
                * 
                * Returns a matrix representation.
                */

                Eigen::MatrixXcd out = Eigen::MatrixXcd::Identity(shape, shape);

                for (int ii = 0; ii < shape-1; ii++)
                {
                    out(ii, shape-1) = _data(ii);
                }

                return out;
            }

            CN inverse() const
            {
                /*! \f{equation*}{ (CN) \rightarrow CN \f}
                * 
                * Returns the inverse.
                */

                return -this->_data;
            }

            // Data representation

            Eigen::VectorXd serialize() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times 1} \f}
                * 
                * Returns a serialized representation.
                */

                const size_t dim = this->get_dimension();
                Eigen::VectorXd out = Eigen::VectorXd::Zero(dim);

                ptrdiff_t kk = 0;
                for (ptrdiff_t ii = 0; ii < dim/2; ii++)
                {
                    out(kk) = std::real(this->_data(ii));
                    kk += 1;
                    out(kk) = std::imag(this->_data(ii));
                    kk += 1;
                }

                return out;
            }

            void unserialize(const Eigen::VectorXd &vector)
            {
                /*! \f{equation*}{ (\mathbb{C}^{n \times 1}) \rightarrow () \f}
                * 
                * Sets the CN object from a serialized vector.
                */
                
                const size_t dim = vector.size();
                Eigen::VectorXcd temp = Eigen::VectorXcd::Zero(dim/2);

                ptrdiff_t kk = 0;
                for (ptrdiff_t ii = 0; ii < static_cast<ptrdiff_t>(dim/2); ii++)
                {
                    temp(ii) = std::complex<double>(vector(kk), vector(kk+1));
                    kk += 2;
                }

                this->_data = temp;
                this->shape = this->_data.size() + 1;
            }

            std::complex<double> operator()(const size_t index) const
            {
                /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{C} \f}
                *
                * Gets a value in the column vector representation.
                */

                return _data(index);
            }

            std::complex<double> & operator()(const size_t index)
            {
                /*! \f{equation*}{ CN(\mathbb{Z}) := \mathbb{C} \f}
                *
                * Assignment of a value in the column vector representation.
                */

                return _data(index);
            }

            std::complex<double> operator()(const size_t index1, const size_t index2) const
            {
                /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{C} \f}
                *
                * Gets a value in the square matrix representation.
                */

                Eigen::MatrixXcd mat = get_ados_representation();
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

            CN operator*(const CN & other) const
            {
                /*! \f{equation*}{ (CN, CN) \rightarrow CN \f}
                *
                * Group product.
                */

                assert(this->shape == other.shape);
                Base out = _data + other._data;
                return out;
            }

            CN & operator*=(const CN & other)
            {
                /*! \f{equation*}{ (CN, CN) \rightarrow CN \f}
                *
                * In place group product.
                */

                assert(this->shape == other.shape);
                this->_data += other._data;
                return *this;
            }

            friend std::ostream & operator<<(std::ostream& os, const CN & other);
        };
        
        std::ostream & operator<<(std::ostream& os, const CN & other)
        {
            /*!
            * Overloads the "<<" stream insertion operator.
            */

            os << static_cast<const Eigen::VectorXcd>(other._data);
            return os;
        }
    }
}

#endif
