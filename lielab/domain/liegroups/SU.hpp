#ifndef _LIELAB_DOMAIN_SU_HPP
#define _LIELAB_DOMAIN_SU_HPP

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab
{
    namespace domain
    {
        class SU
        {
            /*!
            * The SU class.
            */
            public:
            static constexpr bool abelian = false;
            size_t shape = 0;

            Eigen::MatrixXcd _data;
            
            // Initialization methods

            SU() : SU(0)
            {
                /*! \f{equation*}{ () \rightarrow SU \f}
                * 
                * Empty initialization function. Enables instantiation like:
                * 
                *     Lielab::domain::SU x, y, z;
                * 
                */

            }

            SU(const size_t shape)
            {
                /*! \f{equation*}{ (\mathbb{Z}) \rightarrow SU \f}
                *
                * Constructor instantiating an \f$SU\f$ object.
                * 
                * Enables instantiation like:
                * 
                *     Lielab::domain::SU x(2), y(4), z(6);
                * 
                * @param[in] shape The shape of the data matrix.
                */

                this->_data = Eigen::MatrixXcd::Identity(shape, shape);
                this->shape = shape;
            }

            static SU Quaternion()
            {
                /*! \f{equation*}{ () \rightarrow SU \f}
                 *
                 * Constructor instantiating an Identity Quaternion as an \f$SU\f$ object.
                 * 
                 * Enables instatiation like:
                 * 
                 *     Lielab::domain::SU Quaternion0 = Lielab::domain::SU::Quaternion();
                 * 
                 * @param[out] quaternion An SU object representing Identity Quaternion.
                 */

                return SU(2);
            }

            static SU Quaternion(const double e0, const double e1, const double e2, const double e3)
            {
                /*! \f{equation*}{ (\mathbb{R}^4) \rightarrow SU \f}
                 *
                 * Constructor instantiating a Quaternion as an \f$SU\f$ object.
                 * 
                 * Enables instatiation like:
                 * 
                 *     Lielab::domain::SU Quaternion0 = Lielab::domain::SU::Quaternion(1.0, 0.0, 0.0, 0.0);
                 * 
                 * @param[out] quaternion An SU object representing the Quaternion.
                 */

                constexpr std::complex<double> i(0.0,1.0);
                SU qout = SU(2);

                qout._data(0,0) = e0 + e1*i;
                qout._data(1,1) = e0 - e1*i;
                qout._data(0,1) = -e2 + e3*i;
                qout._data(1,0) = e2 + e3*i;

                return qout;
            }

            template<typename OtherDerived>
            SU(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{equation}{(\mathbb{C}^{n \times n}) \rightarrow SU \f}
                *
                * Constructor instantiating an \f$SU\f$ object from an
                * \f$n \times n\f$ imaginary matrix.
                *
                * @param[in] other The object to instantiate from as an imaginary matrix.
                */

                if (other.rows() != other.cols())
                {
                    throw Errorx("Size of the matrix must be square.");
                }

                this->_data = Eigen::MatrixXcd(other);
                this->shape = other.rows();
            }

            template<typename OtherDerived>
            SU & operator=(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{equation}{ SU := \mathbb{C}^{n \times n} \f}
                * 
                * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `SU`.
                */

                if (other.rows() != other.cols())
                {
                    throw Errorx("Size of the matrix must be square.");
                }
                
                this->_data = Eigen::MatrixXcd(other);
                this->shape = other.rows();
                return *this;
            }

            // TODO: Project

            size_t get_dimension() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
                * 
                * Gets the dimension of the group.
                */

                return this->shape * this->shape - 1;
            }

            size_t get_size() const
            {
                /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
                 *
                 * Gets the size of the data representation.
                 */

                return static_cast<size_t>(2*std::pow(this->shape, 2));
            }

            Eigen::MatrixXcd get_ados_representation() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times n} \f}
                * 
                * Returns a matrix representation.
                */

                return this->_data;
            }

            SU inverse() const
            {
                /*! \f{equation*}{ () \rightarrow SU \f}
                * 
                * Returns the inverse.
                */

                return this->_data.inverse();
            }

            Eigen::VectorXd serialize() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times 1} \f}
                * 
                * Returns a serialized representation.
                */

                const size_t shape2 = this->shape*this->shape;
                Eigen::VectorXd out(2*shape2);

                for (size_t jj = 0; jj < this->shape; jj++)
                {
                    for (size_t ii = 0; ii < this->shape; ii++)
                    {
                        out(jj*this->shape + ii) = std::real(this->_data(ii, jj));
                    }
                }

                for (size_t jj = 0; jj < this->shape; jj++)
                {
                    for (size_t ii = 0; ii < this->shape; ii++)
                    {
                        out(shape2 + jj*this->shape + ii) = std::imag(this->_data(ii, jj));
                    }
                }

                return out;
            }

            //    if (this->shape == 2)
            //    {
                // TODO: This is the old serialization
            //        Eigen::VectorXd out(4);
            //        out(0) = (std::real(this->_data(0,0)) + std::real(this->_data(1,1))) / 2.0;
            //        out(1) = (std::imag(this->_data(0,0)) - std::imag(this->_data(1,1))) / 2.0;
            //        out(2) = (std::real(this->_data(1,0)) - std::real(this->_data(0,1))) / 2.0;
            //        out(3) = (std::imag(this->_data(1,0)) + std::imag(this->_data(0,1))) / 2.0;
            //        return out;
            //    }
                
                // TODO: General real serialization strategy
                // Eigen::VectorXcd out(this->shape * this->shape);

                // for (size_t jj = 0; jj < this->shape; jj++)
                // {
                //     for (size_t ii = 0; ii < this->shape; ii++)
                //     {
                //         out(jj*this->shape + ii) = this->_data(ii, jj);
                //     }
                // }

                // return out;
            // }

            void unserialize(const Eigen::VectorXd &vec)
            {
                /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow SU \f}
                * 
                * Sets the SP object from a serialized vector.
                */

                const size_t _size = static_cast<size_t>(this->shape*this->shape);
                Eigen::MatrixXd vreal = vec(Eigen::seqN(0, _size)).reshaped(this->shape, this->shape);
                Eigen::MatrixXd vimag = vec(Eigen::seqN(_size, _size)).reshaped(this->shape, this->shape);
                this->_data = vreal.cast<std::complex<double>>() + std::complex<double>(0.0, 1.0)*(vimag.cast<std::complex<double>>());

                // This is the old unserialization
                // if (vec.size() == 4)
                // {
                //     SU out(2);
                //     out._data(0,0) = std::complex<double>(vec(0), vec(1));
                //     out._data(0,1) = std::complex<double>(-vec(2), vec(3));
                //     out._data(1,0) = std::complex<double>(vec(2), vec(3));
                //     out._data(1,1) = std::complex<double>(vec(0), -vec(1));
                //     return out;
                // }

                // TODO: General real unserialization strategy
            }

            std::complex<double> operator()(const size_t index1, const size_t index2) const
            {
                /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{C} \f}
                *
                * Gets a value in the square matrix representation.
                */

                return this->_data(index1, index2);
            }

            std::complex<double> & operator()(const size_t index1, const size_t index2)
            {
                /*! \f{equation*}{ SU(\mathbb{Z}, \mathbb{Z}) := \mathbb{C} \f}
                *
                * Assignment of a value in the square matrix representation.
                */

                return this->_data(index1, index2);
            }

            SU adjoint() const
            {
                // TODO: Delete this?
                return _data.adjoint();
            }

            size_t rows() const
            {
                // TODO: delete this method
                return this->_data.rows();
            }

            size_t cols() const
            {
                // TODO: delete this method
                return this->_data.cols();
            }

            SU operator*(const SU & other) const
            {
                /*! \f{equation*}{ (SU, SU) \rightarrow SU \f}
                *
                * Group product.
                */

                assert(this->shape == other.shape);
                return this->_data * other._data;
            }

            SU & operator*=(const SU & other)
            {
                /*! \f{equation*}{ (SU, SU) \rightarrow SU \f}
                *
                * In place group product.
                */

                assert(this->shape == other.shape);
                this->_data *= other._data;
                return *this;
            }

            friend std::ostream & operator<<(std::ostream& os, const SU & other);
        };

        std::ostream & operator<<(std::ostream& os, const SU & other)
        {
            /*!
            * Overloads the "<<" stream insertion operator.
            */

            os << static_cast<const Eigen::MatrixXcd>(other._data);
            return os;
        }
    }
}

#endif
