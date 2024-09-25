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

    /*
     * Additional static initializers. Not a part of the core Lie group, but are convenient.
     */
    template <typename T>
    static SO from_eulerangles_body123(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_2 \cos \theta_3 & -\cos \theta_2 \sin \theta_3 & \sin \theta_2 \\
        * \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \sin \theta_2 \sin \theta_3 & - \sin \theta_1 \cos \theta_2 \\
        * \sin \theta_1 \sin \theta_3 - \cos \theta_1 \sin \theta_2 \cos \theta_3 & \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_2
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Body-123 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        *
        * TODO: Citation needed.
        *
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta2*ctheta3;
        dcm(0,1) = -ctheta2*stheta3;
        dcm(0,2) = stheta2;
        dcm(1,0) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
        dcm(1,1) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;
        dcm(1,2) = -stheta1*ctheta2;
        dcm(2,0) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
        dcm(2,1) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
        dcm(2,2) = ctheta1*ctheta2;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_body231(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_1 \cos \theta_2 & \sin \theta_1 \sin \theta_3 - \cos \theta_1 \sin \theta_2 \cos \theta_3 & \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 \\
        * \sin \theta_2 & \cos \theta_2 \cos \theta_3 & - \cos \theta_2 \sin \theta_3 \\
        * - \sin \theta_1 \cos \theta_2 & \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \sin \theta_2 \sin \theta_3
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Body-231 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta1*ctheta2;
        dcm(0,1) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
        dcm(0,2) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
        dcm(1,0) = stheta2;
        dcm(1,1) = ctheta2*ctheta3;
        dcm(1,2) = -ctheta2*stheta3;
        dcm(2,0) = -stheta1*ctheta2;
        dcm(2,1) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
        dcm(2,2) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_body312(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_1 \cos \theta_3 - \sin \theta_1 \sin \theta_2 \sin \theta_3 & - \sin \theta_1 \cos \theta_2 & \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \sin \theta_3 \\
        * \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_2 & \sin \theta_1 \sin \theta_3 - \cos \theta_1 \sin \theta_2 \cos \theta_3 \\
        * - \cos \theta_2 \sin \theta_3 & \sin \theta_2 & \cos \theta_2 \cos \theta_3
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Body-312 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;
        dcm(0,1) = -stheta1*ctheta2;
        dcm(0,2) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
        dcm(1,0) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
        dcm(1,1) = ctheta1*ctheta2;
        dcm(1,2) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
        dcm(2,0) = -ctheta2*stheta3;
        dcm(2,1) = stheta2;
        dcm(2,2) = ctheta2*ctheta3;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_body132(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_2 \cos \theta_3 & - \sin \theta_2 & \cos \theta_2 \sin \theta_3 \\
        * \sin \theta_1 \sin \theta_3 + \cos \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_1 \cos \theta_2 & - \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 \\
        * - \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 & \sin \theta_1 \cos \theta_2 & \cos \theta_1 \cos \theta_3 + \sin \theta_1 \sin \theta_2 \sin \theta_3
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Body-132 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta2*ctheta3;
        dcm(0,1) = -stheta2;
        dcm(0,2) = ctheta2*stheta3;
        dcm(1,0) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
        dcm(1,1) = ctheta1*ctheta2;
        dcm(1,2) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
        dcm(2,0) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
        dcm(2,1) = stheta1*ctheta2;
        dcm(2,2) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_body213(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_1 \cos \theta_3 + \sin \theta_1 \sin \theta_2 \sin \theta_3 & - \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 & \sin \theta_1 \cos \theta_2 \\
        * \cos \theta_2 \sin \theta_3 & \cos \theta_2 \cos \theta_3 & - \sin \theta_2 \\
        * - \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 & \sin \theta_1 \sin \theta_3 + \cos \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_1 \cos \theta_2
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Body-213 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;
        dcm(0,1) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
        dcm(0,2) = stheta1*ctheta2;
        dcm(1,0) = ctheta2*stheta3;
        dcm(1,1) = ctheta2*ctheta3;
        dcm(1,2) = -stheta2;
        dcm(2,0) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
        dcm(2,1) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
        dcm(2,2) = ctheta1*ctheta2;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_body321(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_1 \cos \theta_2 & - \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 & \sin \theta_1 \sin \theta_3 + \cos \theta_1 \sin \theta_2 \cos \theta_3 \\
        * \sin \theta_1 \cos \theta_2 & \cos \theta_1 \cos \theta_3 + \sin \theta_1 \sin \theta_2 \sin \theta_3 & - \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 \\
        * - \sin \theta_2 & \cos \theta_2 \sin \theta_3 & \cos \theta_2 \cos \theta_3
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Body-321 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta1*ctheta2;
        dcm(0,1) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
        dcm(0,2) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
        dcm(1,0) = stheta1*ctheta2;
        dcm(1,1) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;
        dcm(1,2) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
        dcm(2,0) = -stheta2;
        dcm(2,1) = ctheta2*stheta3;
        dcm(2,2) = ctheta2*ctheta3;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_body121(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_2 & \sin \theta_2 \sin \theta_3 & \sin \theta_2 \cos \theta_3 \\
        * \sin \theta_1 \sin \theta_2 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3 & - \cos \theta_1 \sin \theta_3 - \sin \theta_1 \cos \theta_2 \cos \theta_3 \\
        * - \cos \theta_1 \sin \theta_2 & \sin \theta_1 \cos \theta_3 + \cos \theta_1 \cos \theta_2 \sin \theta_3 & - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Body-121 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta2;
        dcm(0,1) = stheta2*stheta3;
        dcm(0,2) = stheta2*ctheta3;
        dcm(1,0) = stheta1*stheta2;
        dcm(1,1) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
        dcm(1,2) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
        dcm(2,0) = -ctheta1*stheta2;
        dcm(2,1) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
        dcm(2,2) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_body131(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_2 & - \sin \theta_2 \cos \theta_3 & \sin \theta_2 \sin \theta_3 \\
        * \cos \theta_1 \sin \theta_2 & - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3 & - \sin \theta_1 \cos \theta_3 - \cos \theta_1 \cos \theta_2 \sin \theta_3 \\
        * \sin \theta_1 \sin \theta_2 & \cos \theta_1 \sin \theta_3 + \sin \theta_1 \cos \theta_2 \cos \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Body-131 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta2;
        dcm(0,1) = -stheta2*ctheta3;
        dcm(0,2) = stheta2*stheta3;
        dcm(1,0) = ctheta1*stheta2;
        dcm(1,1) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
        dcm(1,2) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
        dcm(2,0) = stheta1*stheta2;
        dcm(2,1) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
        dcm(2,2) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_body212(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3 & \sin \theta_1 \sin \theta_2 & \cos \theta_1 \sin \theta_3 + \sin \theta_1 \cos \theta_2 \cos \theta_3 \\
        * \sin \theta_2 \sin \theta_3 & \cos \theta_2 & - \sin \theta_2 \cos \theta_3 \\
        * - \sin \theta_1 \cos \theta_3 - \cos \theta_1 \cos \theta_2 \sin \theta_3 & \cos \theta_1 \sin \theta_2 & - \sin \theta_1 \cos \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Body-212 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
        dcm(0,1) = stheta1*stheta2;
        dcm(0,2) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
        dcm(1,0) = stheta2*stheta3;
        dcm(1,1) = ctheta2;
        dcm(1,2) = -stheta2*ctheta3;
        dcm(2,0) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
        dcm(2,1) = ctheta1*stheta2;
        dcm(2,2) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_body232(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3 & - \cos \theta_1 \sin \theta_2 & \sin \theta_1 \cos \theta_3 + \cos \theta_1 \cos \theta_2 \sin \theta_3 \\
        * \sin \theta_2 \cos \theta_3 & \cos \theta_2 & \sin \theta_2 \sin \theta_3 \\
        * - \cos \theta_1 \sin \theta_3 - \sin \theta_1 \cos \theta_2 \cos \theta_3 & \sin \theta_1 \sin \theta_2 & - \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Body-232 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
        dcm(0,1) = -ctheta1*stheta2;
        dcm(0,2) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
        dcm(1,0) = stheta2*ctheta3;
        dcm(1,1) = ctheta2;
        dcm(1,2) = stheta2*stheta3;
        dcm(2,0) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
        dcm(2,1) = stheta1*stheta2;
        dcm(2,2) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_body313(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3 & - \cos \theta_1 \sin \theta_3 - \sin \theta_1 \cos \theta_2 \cos \theta_3 & \sin \theta_1 \sin \theta_2 \\
        * \sin \theta_1 \cos \theta_3 + \cos \theta_1 \cos \theta_2 \sin \theta_3 & - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3 - \cos \theta_1 \sin \theta_2 \\
        * \sin \theta_2 \sin \theta_3 & \sin \theta_2 \cos \theta_3 & \cos \theta_2
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Body-313 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
        dcm(0,1) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
        dcm(0,2) = stheta1*stheta2;
        dcm(1,0) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
        dcm(1,1) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
        dcm(1,2) = -ctheta1*stheta2;
        dcm(2,0) = stheta2*stheta3;
        dcm(2,1) = stheta2*ctheta3;
        dcm(2,2) = ctheta2;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_body323(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3 & - \sin \theta_1 \cos \theta_3 - \cos \theta_1 \cos \theta_2 \sin \theta_3 & \cos \theta_1 \sin \theta_2 \\
        * \cos \theta_1 \sin \theta_3 + \sin \theta_1 \cos \theta_2 \cos \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3 & \sin \theta_1 \sin \theta_2 \\
        * - \sin \theta_2 \cos \theta_3 & \sin \theta_2 \sin \theta_3 & \cos \theta_2
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Body-323 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
        dcm(0,1) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
        dcm(0,2) = ctheta1*stheta2;
        dcm(1,0) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
        dcm(1,1) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
        dcm(1,2) = stheta1*stheta2;
        dcm(2,0) = -stheta2*ctheta3;
        dcm(2,1) = stheta2*stheta3;
        dcm(2,2) = ctheta2;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_space123(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_2 \cos \theta_3 & - \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 & \sin \theta_1 \sin \theta_3 + \cos \theta_1 \sin \theta_2 \cos \theta_3 \\
        * \cos \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_3 + \sin \theta_1 \sin \theta_2 \sin \theta_3 & - \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 \\
        * - \sin \theta_2 & \sin \theta_1 \cos \theta_2 & \cos \theta_1 \cos \theta_2
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Space-123 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta2*ctheta3;
        dcm(0,1) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
        dcm(0,2) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
        dcm(1,0) = ctheta2*stheta3;
        dcm(1,1) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;
        dcm(1,2) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
        dcm(2,0) = -stheta2;
        dcm(2,1) = stheta1*ctheta2;
        dcm(2,2) = ctheta1*ctheta2;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_space231(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_1 \cos \theta 2 & - \sin \theta_2 & \sin \theta_1 \cos \theta_2 \\
        * \sin \theta_1 \sin \theta_3 + \cos \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_2 \cos \theta_3 & - \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 \\
        * - \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 & \cos \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_3 + \sin \theta_1 \sin \theta_2 \sin \theta_3
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Space-231 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta1*ctheta2;
        dcm(0,1) = -stheta2;
        dcm(0,2) = stheta1*ctheta2;
        dcm(1,0) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
        dcm(1,1) = ctheta2*ctheta3;
        dcm(1,2) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
        dcm(2,0) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
        dcm(2,1) = ctheta2*stheta3;
        dcm(2,2) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_space312(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_1 \cos \theta_3 + \sin \theta_1 \sin \theta_2 \sin \theta_3 & - \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 & \cos \theta_2 \sin \theta_3 \\
        * \sin \theta_1 \cos \theta_2 & \cos \theta_1 \cos \theta_2 & - \sin \theta_2 \\
        * - \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 & \sin \theta_1 \sin \theta_3 + \cos \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_2 \cos \theta_3
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Space-312 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;
        dcm(0,1) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
        dcm(0,2) = ctheta2*stheta3;
        dcm(1,0) = stheta1*ctheta2;
        dcm(1,1) = ctheta1*ctheta2;
        dcm(1,2) = -stheta2;
        dcm(2,0) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
        dcm(2,1) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
        dcm(2,2) = ctheta2*ctheta3;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_space132(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_2 \cos \theta_3 & \sin \theta_1 \sin \theta_3 - \cos \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 \\
        * \sin \theta_2 & \cos \theta_1 \cos \theta_2 & \sin \theta_1 \cos \theta_2 \\
        * - \cos \theta_2 \sin \theta_3 & \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \sin \theta_2 \sin \theta_3
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Space-132 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta2*ctheta3;
        dcm(0,1) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
        dcm(0,2) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
        dcm(1,0) = stheta2;
        dcm(1,1) = ctheta1*ctheta2;
        dcm(1,2) = -stheta1*ctheta2;
        dcm(2,0) = -ctheta2*stheta3;
        dcm(2,1) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
        dcm(2,2) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_space213(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_1 \cos \theta_3 - \sin \theta_1 \sin \theta_2 \sin \theta_3 & - \cos \theta_2 \sin \theta_3 & \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 \\
        * \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_2 \cos \theta_3 & \sin \theta_1 \sin \theta_3 - \cos \theta_1 \sin \theta_2 \cos \theta_3 \\
        * - \sin \theta_1 \cos \theta_2 & \sin \theta_2 & \cos \theta_1 \cos \theta_2
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Space-213 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;
        dcm(0,1) = -ctheta2*stheta3;
        dcm(0,2) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
        dcm(1,0) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
        dcm(1,1) = ctheta2*ctheta3;
        dcm(1,2) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
        dcm(2,0) = -stheta1*ctheta2;
        dcm(2,1) = stheta2;
        dcm(2,2) = ctheta1*ctheta2;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_space321(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_1 \cos \theta_2 & - \sin \theta_1 \cos \theta_2 & \sin \theta_2 \\
        * \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \sin \theta_2 \sin \theta_3 & - \cos \theta_2 \sin \theta_3 \\
        * \sin \theta_1 \sin \theta_3 - \cos \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_2 \cos \theta_3
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Space-321 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta1*ctheta2;
        dcm(0,1) = -stheta1*ctheta2;
        dcm(0,2) = stheta2;
        dcm(1,0) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
        dcm(1,1) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;
        dcm(1,2) = -ctheta2*stheta3;
        dcm(2,0) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
        dcm(2,1) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
        dcm(2,2) = ctheta2*ctheta3;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_space121(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_2 & \sin \theta_1 \sin \theta_2 & \cos \theta_1 \sin \theta_2 \\
        * \sin \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3 & - \sin \theta_1 \cos \theta_3 - \cos \theta_1 \cos \theta_2 \sin \theta_3 \\
        * - \sin \theta_2 \cos \theta_3 & \cos \theta_1 \sin \theta_3 + \sin \theta_1 \cos \theta_2 \cos \theta_3 & - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Space-121 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta2;
        dcm(0,1) = stheta1*stheta2;
        dcm(0,2) = ctheta1*stheta2;
        dcm(1,0) = stheta2*stheta3;
        dcm(1,1) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
        dcm(1,2) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
        dcm(2,0) = -stheta2*ctheta3;
        dcm(2,1) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
        dcm(2,2) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_space131(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_2 & - \cos \theta_1 \sin \theta_2 & \sin \theta_1 \sin \theta_2 \\
        * \sin \theta_2 \cos \theta_3 & - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3 & - \cos \theta_1 \sin \theta_3 - \sin \theta_1 \cos \theta_2 \cos \theta_3 \\
        * \sin \theta_2 \sin \theta_3 & \sin \theta_1 \cos \theta_3 + \cos \theta_1 \cos \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Space-131 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta2;
        dcm(0,1) = -ctheta1*stheta2;
        dcm(0,2) = stheta1*stheta2;
        dcm(1,0) = stheta2*ctheta3;
        dcm(1,1) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
        dcm(1,2) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
        dcm(2,0) = stheta2*stheta3;
        dcm(2,1) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
        dcm(2,2) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_space212(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3 & \sin \theta_2 \sin \theta_3 & \sin \theta_1 \cos \theta_3 + \cos \theta_1 \cos \theta_2 \sin \theta_3 \\
        * \sin \theta_1 \sin \theta_2 & \cos \theta_2 & - \cos \theta_1 \sin \theta_2 \\
        * - \cos \theta_1 \sin \theta_3 - \sin \theta_1 \cos \theta_2 \cos \theta_3 & \sin \theta_2 \cos \theta_3 & - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Space-212 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
        dcm(0,1) = stheta2*stheta3;
        dcm(0,2) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
        dcm(1,0) = stheta1*stheta2;
        dcm(1,1) = ctheta2;
        dcm(1,2) = -ctheta1*stheta2;
        dcm(2,0) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
        dcm(2,1) = stheta2*ctheta3;
        dcm(2,2) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_space232(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3 & - \sin \theta_2 \cos \theta_3 & \cos \theta_1 \sin \theta_3 + \sin \theta_1 \cos \theta_2 \cos \theta_3 \\
        * \cos \theta_1 \sin \theta_2 & \cos \theta_2 & \sin \theta_1 \sin \theta_2 \\
        * - \sin \theta_1 \cos \theta_3 - \cos \theta_1 \cos \theta_2 \sin \theta_3 & \sin \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Space-232 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
        dcm(0,1) = -stheta2*ctheta3;
        dcm(0,2) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
        dcm(1,0) = ctheta1*stheta2;
        dcm(1,1) = ctheta2;
        dcm(1,2) = stheta1*stheta2;
        dcm(2,0) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
        dcm(2,1) = stheta2*stheta3;
        dcm(2,2) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_space313(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3 & - \sin \theta_1 \cos \theta_3 - \cos \theta_1 \cos \theta_2 \sin \theta_3 & \sin \theta_2 \sin \theta_3 \\
        * \cos \theta_1 \sin \theta_3 + \sin \theta_1 \cos \theta_2 \cos \theta_3 & - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3 & - \sin \theta_2 \cos \theta_3 \\
        * \sin \theta_1 \sin \theta_2 & \cos \theta_1 \sin \theta_2 & \cos \theta_2
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Space-313 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
        dcm(0,1) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
        dcm(0,2) = stheta2*stheta3;
        dcm(1,0) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
        dcm(1,1) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
        dcm(1,2) = -stheta2*ctheta3;
        dcm(2,0) = stheta1*stheta2;
        dcm(2,1) = ctheta1*stheta2;
        dcm(2,2) = ctheta2;

        return dcm;
    }

    template <typename T>
    static SO from_eulerangles_space323(const T theta1, const T theta2, const T theta3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \begin{bmatrix}
        * - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3 & - \cos \theta_1 \sin \theta_3 - \sin \theta_1 \cos \theta_2 \cos \theta_3 & \sin \theta_2 \cos \theta_3 \\
        * \sin \theta_1 \cos \theta_3 + \cos \theta_1 \cos \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3 & \sin \theta_2 \sin \theta_3 \\
        * - \cos \theta_1 \sin \theta_2 & \sin \theta_1 \sin \theta_2 & \cos \theta_2
        * \end{bmatrix}\f}
        *
        * Transforms an Euler Angle Space-323 rotation sequence into a direction co-sine matrix.
        * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
        * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
        * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
        * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
        */

        SO dcm(3);

        const T stheta1 = std::sin(theta1);
        const T ctheta1 = std::cos(theta1);
        const T stheta2 = std::sin(theta2);
        const T ctheta2 = std::cos(theta2);
        const T stheta3 = std::sin(theta3);
        const T ctheta3 = std::cos(theta3);

        dcm(0,0) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
        dcm(0,1) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
        dcm(0,2) = stheta2*ctheta3;
        dcm(1,0) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
        dcm(1,1) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
        dcm(1,2) = stheta2*stheta3;
        dcm(2,0) = -ctheta1*stheta2;
        dcm(2,1) = stheta1*stheta2;
        dcm(2,2) = ctheta2;

        return dcm;
    }

    template <typename T>
    static SO from_quaternion(const T e0, const T e1, const T e2, const T e3)
    {
        /*! \f{equation*}{ (\mathbb{R}^4) \rightarrow SO \f}
         *
         * Constructor instantiating a Quaternion as an \f$SO\f$ object.
         * 
         * Enables instatiation like:
         * 
         *     Lielab::domain::SO DCM0 = Lielab::domain::SO::from_quaternion(1.0, 0.0, 0.0, 0.0);
         *
         *
         * \f{equation*}{dcm = \begin{bmatrix}
         * e_0^2 + e_1^2 - e_2^2 - e_3^2 & 2(e1 e2 - e3 e0) & 2(e1 e3 + e2 e0) \\
         * 2(e1 e2 + e3 e0) & e_0^2 - e_1^2 + e_2^2 - e_3^2 & 2(e2 e3 - e1 e0) \\
         * 2(e1 e3 - e2 e0) & 2(e2 e3 + e1 e0) & e_0^2 - e_1^2 - e_2^2 + e_3^2
         * \end{bmatrix}\f}
         * 
         * @param[out] dcm An SO object representing the input quaternion.
         *
         * TODO: Citation needed
         *
         */


        SO dcm(3);

        dcm(0,0) = e0*e0 + e1*e1 - e2*e2 - e3*e3;
        dcm(0,1) = 2.0*e1*e2 - 2.0*e3*e0;
        dcm(0,2) = 2.0*e1*e3 + 2.0*e2*e0;
        dcm(1,0) = 2.0*e1*e2 + 2.0*e3*e0;
        dcm(1,1) = e0*e0 - e1*e1 + e2*e2 - e3*e3;
        dcm(1,2) = 2.0*e2*e3 - 2.0*e1*e0;
        dcm(2,0) = 2.0*e1*e3 - 2.0*e2*e0;
        dcm(2,1) = 2.0*e2*e3 + 2.0*e1*e0;
        dcm(2,2) = e0*e0 - e1*e1 - e2*e2 + e3*e3;

        return dcm;
    }

    template <typename T>
    static SO from_rodriguesvector(const T g1, const T g2, const T g3)
    {
        /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
        *
        * \f{equation*}{dcm = \frac{1}{1 + g_1^2 + g_2^2 + g_3^2}\begin{bmatrix}
        * 1 + g_1^2 - g_2^2 - g_3^2 & 2(g_1 g_2 + g_3) & 2(g_1 g_3 - g_2) \\
        * 2(g_1 g_2 - g_3) & 1 - g_1^2 + g_2^2 - g_3^2 & 2(g_2 g_3 + g_1) \\
        * 2(g_1 g_3 + g_2) & 2(g_2 g_3 - g_1) & 1 - g_1^2 - g_2^2 + g_3^2
        * \end{bmatrix}\f}
        * 
        * Transforms a Rodrigues (Gibbs) vector into a direction co-sine matrix.
        * @param[in] g1 The 1st element in Rodrigues representation.
        * @param[in] g2 The 2nd element in Rodrigues representation.
        * @param[in] g3 The 3rd element in Rodrigues representation.
        * @param[out] dcm The direction co-sine matrix.
        * 
        * TODO: Citation needed
        * TODO: Check conventions. This might be transposed from what we use elsewhere.
        * 
        */

        SO dcm(3);

        const T g12 = std::pow(g1, T(2.0));
        const T g22 = std::pow(g2, T(2.0));
        const T g32 = std::pow(g3, T(2.0));
        const T mul = 1.0/(1.0 + g12 + g22 + g32);

        dcm(0,0) = mul*(1.0 + g12 - g22 - g32);
        dcm(0,1) = mul*(2.0*(g1*g2 + g3));
        dcm(0,2) = mul*(2.0*(g1*g3 - g2));
        dcm(1,0) = mul*(2.0*(g1*g2 - g3));
        dcm(1,1) = mul*(1.0 - g12 + g22 - g32);
        dcm(1,2) = mul*(2.0*(g2*g3 + g1));
        dcm(2,0) = mul*(2.0*(g1*g3 + g2));
        dcm(2,1) = mul*(2.0*(g2*g3 - g1));
        dcm(2,2) = mul*(1.0 - g12 - g22 + g32);

        return dcm;
    }
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
