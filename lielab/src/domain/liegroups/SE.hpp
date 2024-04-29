#ifndef _LIELAB_DOMAIN_SE_HPP
#define _LIELAB_DOMAIN_SE_HPP

#include <Eigen/Dense>

namespace lielab
{
    namespace domain
    {
        class SE
        {
            /*!
            * The SE class.
            */
            public:
            static constexpr bool abelian = false;
            size_t shape = 0;
            
            Eigen::MatrixXd _data;

            SE() : SE(0)
            {
                /*! \f{equation*}{ () \rightarrow SO \f}
                * 
                * Empty initialization function. Enables instantiation like:
                * 
                *     lielab::domain::SE x, y, z;
                * 
                */

            }

            SE(const size_t shape)
            {
                /*! \f{equation*}{ (\mathbb{Z}) \rightarrow SE \f}
                *
                * Constructor instantiating an \f$SO\f$ object.
                * 
                * Enables instantiation like:
                * 
                *     lielab::domain::SE x(2), y(3), z(4);
                * 
                * @param[in] shape The shape of the data matrix.
                */

                this->_data = Eigen::MatrixXd::Identity(shape, shape);
                this->shape = shape;
            }

            template<typename OtherDerived>
            SE(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{equation}{ (\mathbb{R}^{n \times n}) \rightarrow SE \f}
                *
                * Constructor instantiating an \f$SE\f$ object from an
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
            SE & operator=(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{equation}{ SE := \mathbb{R}^{n \times n} \f}
                * 
                * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `SE`.
                */

                if (other.rows() != other.cols())
                {
                    throw Errorx("Size of the matrix must be square.");
                }
                
                this->_data = Eigen::MatrixXd(other);
                this->shape = this->_data.rows();
                return *this;
            }

            // static Eigen::MatrixXd project(const Eigen::MatrixXd & other)
            // {
            //     /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{SE} \f}
            //     *
            //     * Projects a matrix suitable for data.
            //     */

            //     if (other.rows() != other.cols())
            //     {
            //         throw Errorx("Size of the matrix must be square.");
            //     }

            //     Eigen::HouseholderQR<Eigen::MatrixXd> qr(other);
            //     Eigen::MatrixXd Q = qr.householderQ();
            //     Eigen::MatrixXd R = qr.matrixQR();
            //     Eigen::MatrixXd P = Eigen::MatrixXd::Zero(other.rows(), other.cols());

            //     for (int ii = 0; ii < other.rows(); ii++)
            //     {
            //         if (R(ii,ii) < 0)
            //         {
            //             P(ii, ii) = -1.0;
            //         }
            //         else
            //         {
            //             P(ii, ii) = 1.0;
            //         }
            //     }

            //     return Q*P;
            // }

            size_t get_dimension() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
                * 
                * Gets the dimension of the group.
                */

                return (this->shape - 1) * (this->shape - 2) / 2 + this->shape - 1;
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

            SE inverse() const
            {
                /*! \f{equation*}{ (SE) \rightarrow SE \f}
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

            SE operator*(const SE & other) const
            {
                /*! \f{equation*}{ (SE, SE) \rightarrow SE \f}
                *
                * Group product.
                */

                assert(this->shape == other.shape);
                return this->_data * other._data;
            }

            SE & operator*=(const SE & other)
            {
                /*! \f{equation*}{ (SO, SO) \rightarrow SO \f}
                *
                * In place group product.
                */

                assert(this->shape == other.shape);
                this->_data *= other._data;
                return *this;
            }

            friend std::ostream & operator<<(std::ostream & os, const SE & other);
        };

        std::ostream & operator<<(std::ostream& os, const SE & other)
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
