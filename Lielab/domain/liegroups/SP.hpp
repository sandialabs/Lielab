#ifndef _LIELAB_DOMAIN_SP_HPP
#define _LIELAB_DOMAIN_SP_HPP

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab
{
    namespace domain
    {
        class SP
        {
            /*!
            * The SP class.
            */
            public:
            static constexpr bool abelian = false;
            size_t shape = 0;

            Eigen::MatrixXd _data;

            // Initialization methods

            SP() : SP(0)
            {
                /*! \f{equation*}{ () \rightarrow SP \f}
                * 
                * Empty initialization function. Enables instantiation like:
                * 
                *     Lielab::domain::SP x, y, z;
                * 
                */

            }

            SP(const size_t shape)
            {
                /*! \f{equation*}{ (\mathbb{Z}) \rightarrow SP \f}
                *
                * Constructor instantiating an \f$SP\f$ object.
                * 
                * Enables instantiation like:
                * 
                *     Lielab::domain::SP x(2), y(4), z(6);
                * 
                * @param[in] shape The shape of the data matrix.
                */

                if (shape % 2 != 0)
                {
                    throw Errorx("Shape of sp must be even dimensional.");
                }

                this->_data = Eigen::MatrixXd::Identity(shape, shape);
                this->shape = shape;
            }

            template<typename OtherDerived>
            SP(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{equation}{(\mathbb{R}^{n \times n}) \rightarrow SP \f}
                *
                * Constructor instantiating an \f$SP\f$ object from an
                * \f$n \times n\f$ real matrix.
                *
                * @param[in] other The object to instantiate from as a real matrix.
                */

                if (other.rows() != other.cols())
                {
                    throw Errorx("Size of the matrix must be square.");
                }

                if (other.rows() % 2 != 0)
                {
                    throw Errorx("Input matrix must be even dimensional.");
                }

                this->_data = Eigen::MatrixXd(other);
                this->shape = other.rows();
            }

            template<typename OtherDerived>
            SP & operator=(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{equation}{ SP := \mathbb{R}^{n \times n} \f}
                * 
                * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `SP`.
                */

                if (other.rows() != other.cols())
                {
                    throw Errorx("Size of the matrix must be square.");
                }

                if (other.rows() % 2 != 0)
                {
                    throw Errorx("Input matrix must be even dimensional.");
                }
                
                this->_data = Eigen::MatrixXd(other);
                this->shape = this->_data.rows();
                return *this;
            }

            // TODO: project

            size_t get_dimension() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
                * 
                * Gets the dimension of the group.
                */

                return this->shape * (this->shape + 1) / 2;
            }

            size_t get_size() const
            {
                /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
                 *
                 * Gets the size of the data representation.
                 */

                return static_cast<size_t>(std::pow(this->shape, 2));
            }

            Eigen::MatrixXd get_ados_representation() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times n} \f}
                * 
                * Returns a matrix representation.
                */

                return this->_data;
            }

            SP inverse() const
            {
                /*! \f{equation*}{ (SP) \rightarrow SP \f}
                * 
                * Returns the inverse.
                */

                return this->_data.inverse();
            }

            Eigen::VectorXd serialize() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
                * 
                * Returns a serialized representation.
                */

                Eigen::VectorXd out(this->shape * this->shape);

                for (size_t jj = 0; jj < this->shape; jj++)
                {
                    for (size_t ii = 0; ii < this->shape; ii++)
                    {
                        out(jj*this->shape + ii) = this->_data(ii, jj);
                    }
                }

                return out;
            }

            void unserialize(const Eigen::VectorXd &vec)
            {
                /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow () \f}
                * 
                * Sets the SP object from a serialized vector.
                */

                this->_data = vec(Eigen::seqN(0, this->shape*this->shape)).reshaped(this->shape, this->shape);
            }

            SP adjoint() const
            {
                // TODO: Delete this?
                return _data.adjoint();
            }

            size_t rows() const
            {
                // TODO: Delete this method.
                return this->_data.rows();
            }

            size_t cols() const
            {
                // TODO: Delete this method.
                return this->_data.cols();
            }

            // TODO: Unserialize

            double operator()(const size_t index1, const size_t index2) const
            {
                /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
                *
                * Gets a value in the square matrix representation.
                */
                
                return this->_data(index1, index2);
            }

            double & operator()(const size_t index1, const size_t index2)
            {
                /*! \f{equation*}{ SP(\mathbb{Z}, \mathbb{Z}) := \mathbb{R} \f}
                *
                * Assignment of a value in the square matrix representation.
                */

                return this->_data(index1, index2);
            }

            SP operator*(const SP & other) const
            {
                /*! \f{equation*}{ (SP, SP) \rightarrow SP \f}
                *
                * Group product.
                */

                assert(this->shape == other.shape);
                return this->_data * other._data;
            }

            SP & operator*=(const SP & other)
            {
                /*! \f{equation*}{ (SP, SP) \rightarrow SP \f}
                *
                * In place group product.
                */

                assert(this->shape == other.shape);
                this->_data *= other._data;
                return *this;
            }

            friend std::ostream & operator<<(std::ostream& os, const SP & other);
        };

        std::ostream & operator<<(std::ostream& os, const SP & other)
        {
            /*!
            * Overloads the "<<" stream insertion operator.
            */

            os << static_cast<const Eigen::MatrixXd>(other._data);
            return os;
        }
    }
}

#endif
