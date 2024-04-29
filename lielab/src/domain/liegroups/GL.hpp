#ifndef _LIELAB_DOMAIN_GL_HPP
#define _LIELAB_DOMAIN_GL_HPP

namespace lielab
{
    namespace domain
    {
        class GL
        {
            /*!
            * The GL class.
            */
            public:
            static constexpr bool abelian = false;
            size_t shape = 0;

            typedef Eigen::MatrixXd Base;
            Base _data;

            // Initialization methods

            GL() : GL(0)
            {
                /*! \f{equation*}{ () \rightarrow GL \f}
                * 
                * Empty initialization function. Enables instantiation like:
                * 
                *     lielab::domain::GL x, y, z;
                * 
                */

            }

            GL(const size_t shape)
            {
                /*! \f{equation*}{ (\mathbb{Z}) \rightarrow GL \f}
                *
                * Constructor instantiating an \f$GL\f$ object.
                * 
                * Enables instantiation like:
                * 
                *     lielab::domain::GL x(2), y(3), z(4);
                * 
                * @param[in] shape The shape of the data matrix.
                */
                
                this->_data = Base::Identity(shape, shape);
                this->shape = shape;
            }

            template<typename OtherDerived>
            GL(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& GL \\ (\mathbb{R}^{n \times 1}) &\rightarrow& GL \f}
                *
                * Constructor instantiating an \f$GL\f$ object from either an
                * \f$n \times n\f$ real matrix or \f$n \times 1\f$ real vector.
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
                    throw Errorx("Input data to GL malformed.");
                }
            }

            template<typename OtherDerived>
            GL & operator=(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{eqnarray*}{ GL &:= \mathbb{R}^{n \times n} \\ &:= \mathbb{R}^{n \times 1} \f}
                * 
                * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `GL`.
                */

                if (other.rows() == other.cols())
                {   
                    shape = other.rows();
                    _data = other;
                }
                return *this;
            }

            static Eigen::MatrixXd project(const Eigen::MatrixXd & other)
            {
                /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in GL \f}
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

                return static_cast<size_t>(std::pow(this->shape, 2));
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

            GL inverse() const
            {
                /*! \f{equation*}{ (GL) \rightarrow GL \f}
                * 
                * Returns the inverse.
                */

                return this->_data.inverse();
            }

            // Data representation

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
                /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
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

            GL operator*(const GL & other) const
            {
                /*! \f{equation*}{ (GL, GL) \rightarrow GL \f}
                *
                * Group product.
                */

                assert(this->shape == other.shape);
                Base out = _data * other._data;
                return out;
            }

            GL & operator*=(const GL & other)
            {
                /*! \f{equation*}{ (GL, GL) \rightarrow GL \f}
                *
                * In place group product.
                */

                assert(this->shape == other.shape);
                this->_data *= other._data;
                return *this;
            }

            friend std::ostream & operator<<(std::ostream& os, const GL & other);
        };
        
        std::ostream & operator<<(std::ostream& os, const GL & other)
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
