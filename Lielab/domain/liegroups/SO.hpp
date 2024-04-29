#ifndef _LIELAB_DOMAIN_SO_HPP
#define _LIELAB_DOMAIN_SO_HPP

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab
{
    namespace domain
    {
        class SO
        {
            /*!
            * The SO class.
            */
            public:
            static constexpr bool abelian = false;
            size_t shape = 0;
            
            Eigen::MatrixXd _data;

            SO() : SO(0)
            {
                /*! \f{equation*}{ () \rightarrow SO \f}
                * 
                * Empty initialization function. Enables instantiation like:
                * 
                *     Lielab::domain::SO x, y, z;
                * 
                */

            }

            SO(const size_t shape)
            {
                /*! \f{equation*}{ (\mathbb{Z}) \rightarrow SO \f}
                *
                * Constructor instantiating an \f$SO\f$ object.
                * 
                * Enables instantiation like:
                * 
                *     Lielab::domain::SO x(2), y(3), z(4);
                * 
                * @param[in] shape The shape of the data matrix.
                */

                this->_data = Eigen::MatrixXd::Identity(shape, shape);
                this->shape = shape;
            }

            template<typename OtherDerived>
            SO(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{equation}{ (\mathbb{R}^{n \times n}) \rightarrow SO \f}
                *
                * Constructor instantiating an \f$SO\f$ object from an
                * \f$n \times n\f$ real matrix.
                *
                * @param[in] other The object to instantiate from as a real matrix.
                */

                if (other.rows() != other.cols())
                {
                    throw Errorx("Size of the matrix must be square.");
                }

                this->_data = Eigen::MatrixXd(other);
                this->shape = this->_data.rows();
            }

            template<typename OtherDerived>
            SO & operator=(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{equation}{ SO := \mathbb{R}^{n \times n} \f}
                * 
                * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `SO`.
                */

                if (other.rows() != other.cols())
                {
                    throw Errorx("Size of the matrix must be square.");
                }
                
                this->_data = Eigen::MatrixXd(other);
                this->shape = this->_data.rows();
                return *this;
            }

            static Eigen::MatrixXd project(const Eigen::MatrixXd & other)
            {
                /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{SO} \f}
                *
                * Projects a matrix suitable for data.
                */

                if (other.rows() != other.cols())
                {
                    throw Errorx("Size of the matrix must be square.");
                }

                Eigen::HouseholderQR<Eigen::MatrixXd> qr(other);
                Eigen::MatrixXd Q = qr.householderQ();
                Eigen::MatrixXd R = qr.matrixQR();
                Eigen::MatrixXd P = Eigen::MatrixXd::Zero(other.rows(), other.cols());

                for (int ii = 0; ii < other.rows(); ii++)
                {
                    if (R(ii,ii) < 0)
                    {
                        P(ii, ii) = -1.0;
                    }
                    else
                    {
                        P(ii, ii) = 1.0;
                    }
                }

                return Q*P;
            }

            size_t get_dimension() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
                * 
                * Gets the dimension of the group.
                */

                return this->shape * (this->shape - 1) / 2;
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

                return _data;
            }

            SO inverse() const
            {
                /*! \f{equation*}{ (SO) \rightarrow SO \f}
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
                * Sets the GL object from a serialized vector.
                */

                this->_data = vec(Eigen::seqN(0, this->shape*this->shape)).reshaped(this->shape, this->shape);
            }

            double operator()(const size_t index1, const size_t index2) const
            {
                /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
                *
                * Gets a value in the square matrix representation.
                */

                return _data(index1, index2);
            }

            double & operator()(const size_t index1, const size_t index2)
            {
                /*! \f{equation*}{ SO(\mathbb{Z}, \mathbb{Z}) := \mathbb{R} \f}
                *
                * Assignment of a value in the square matrix representation.
                */

                return _data(index1, index2);
            }

            size_t rows() const
            {
                // TODO: Remove this
                return this->_data.rows();
            }

            size_t cols() const
            {
                // TODO: Remove this
                return this->_data.cols();
            }

            SO operator*(const SO & other) const
            {
                /*! \f{equation*}{ (SO, SO) \rightarrow SO \f}
                *
                * Group product.
                */

                assert(this->shape == other.shape);
                return this->_data * other._data;
            }

            SO & operator*=(const SO & other)
            {
                /*! \f{equation*}{ (SO, SO) \rightarrow SO \f}
                *
                * In place group product.
                */

                assert(this->shape == other.shape);
                this->_data *= other._data;
                return *this;
            }

            friend std::ostream & operator<<(std::ostream & os, const SO & other);
        };

        std::ostream & operator<<(std::ostream& os, const SO & other)
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
