#ifndef _LIELAB_DOMAIN_GLC_HPP
#define _LIELAB_DOMAIN_GLC_HPP

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab
{
    namespace domain
    {
        class GLC
        {
            /*!
            * The GLC class.
            */
            public:
            static constexpr bool abelian = false;
            size_t shape = 0;

            typedef Eigen::MatrixXcd Base;
            Base _data;

            // Initialization methods

            GLC() : GLC(0)
            {
                /*! \f{equation*}{ () \rightarrow GLC \f}
                * 
                * Empty initialization function. Enables instantiation like:
                * 
                *     Lielab::domain::GLC x, y, z;
                * 
                */

            }

            GLC(const size_t shape)
            {
                /*! \f{equation*}{ (\mathbb{Z}) \rightarrow GLC \f}
                *
                * Constructor instantiating an \f$GLC\f$ object.
                * 
                * Enables instantiation like:
                * 
                *     Lielab::domain::GLC x(2), y(3), z(4);
                * 
                * @param[in] shape The shape of the data matrix.
                */
                
                this->_data = Base::Identity(shape, shape);
                this->shape = shape;
            }

            template<typename OtherDerived>
            GLC(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& GLC \\ (\mathbb{R}^{n \times 1}) &\rightarrow& GLC \f}
                *
                * Constructor instantiating an \f$GLC\f$ object from either an
                * \f$n \times n\f$ imaginary matrix or \f$n \times 1\f$ imaginary vector.
                *
                * @param[in] other The object to instantiate from as a real matrix.
                */

                if (other.rows() == other.cols())
                {
                    this->_data = other;
                    this->shape = other.rows();
                }
                else
                {
                    throw Errorx("Input data to GLC malformed.");
                }
            }

            template<typename OtherDerived>
            GLC & operator=(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{eqnarray*}{ GLC &:= \mathbb{C}^{n \times n} \\ &:= \mathbb{C}^{n \times 1} \f}
                * 
                * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `GLC`.
                */

                if (other.rows() == other.cols())
                {   
                    shape = other.rows();
                    _data = other;
                }
                return *this;
            }

            static Eigen::MatrixXcd project(const Eigen::MatrixXcd & other)
            {
                /*! \f{equation*}{ (\mathbb{C}^{n \times n}) \rightarrow \mathbb{C}^{n \times n} \in GLC \f}
                *
                * Projects a matrix suitable for data.
                */

                if (other.rows() != other.cols())
                {
                    throw Errorx("Size of the matrix must be square.");
                }

                return other;
            }

            size_t get_dimension() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
                * 
                * Gets the dimension of the group.
                */

                return static_cast<size_t>(2*std::pow(this->shape, 2));
            }

            size_t get_size() const
            {
                /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
                 *
                 * Gets the size of the data representation.
                 */

                return static_cast<size_t>(2*std::pow(this->shape, 2));
            }

            Eigen::MatrixXcd get_matrix() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times n} \f}
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

                return _data;
            }

            GLC inverse() const
            {
                /*! \f{equation*}{ (GLC) \rightarrow GLC \f}
                * 
                * Returns the inverse.
                */

                return this->_data.inverse();
            }

            // Data representation

            Eigen::VectorXd serialize() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times 1} \f}
                * 
                * Returns a serialized representation.
                */
                
                const ptrdiff_t dim = this->get_dimension();
                const Eigen::MatrixXcd A = this->get_matrix();

                Eigen::VectorXd out = Eigen::VectorXd::Zero(this->get_dimension());
                ptrdiff_t kk = 0;
                for (ptrdiff_t ii = 0; ii < this->shape; ii++)
                {
                    for (ptrdiff_t jj = 0; jj < this->shape; jj++)
                    {
                        out(kk) = std::real(A(ii,jj));
                        out(kk+1) = std::imag(A(ii,jj));
                        kk = kk + 2;
                    }
                }

                return out;
            }

            void unserialize(const Eigen::VectorXd &vector)
            {
                /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow () \f}
                * 
                * Sets the GL object from a serialized vector.
                */
                
                const ptrdiff_t dim = vector.size();
                this->shape = static_cast<size_t>(std::sqrt(dim/2));
                this->_data = Eigen::MatrixXcd::Zero(this->shape, this->shape);
                
                ptrdiff_t kk = 0;
                for (ptrdiff_t ii = 0; ii < this->shape; ii++)
                {
                    for (ptrdiff_t jj = 0; jj < this->shape; jj++)
                    {
                        this->_data(ii,jj) = std::complex<double>(vector(kk), vector(kk+1));
                        kk = kk + 2;
                    }
                }
            }

            std::complex<double> operator()(const size_t index1, const size_t index2) const
            {
                /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{C} \f}
                *
                * Gets a value in the square matrix representation.
                */

                return _data(index1, index2);
            }

            std::complex<double> & operator()(const size_t index1, const size_t index2)
            {
                /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{C} \f}
                *
                * Gets a value in the square matrix representation.
                */

                return _data(index1, index2);
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

            GLC operator*(const GLC & other) const
            {
                /*! \f{equation*}{ (GLC, GLC) \rightarrow GLC \f}
                *
                * Group product.
                */

                assert(this->shape == other.shape);
                Base out = _data * other._data;
                return out;
            }

            GLC & operator*=(const GLC & other)
            {
                /*! \f{equation*}{ (GLC, GLC) \rightarrow GLC \f}
                *
                * In place group product.
                */

                assert(this->shape == other.shape);
                this->_data *= other._data;
                return *this;
            }

            friend std::ostream & operator<<(std::ostream& os, const GLC & other);
        };
        
        std::ostream & operator<<(std::ostream& os, const GLC & other)
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
