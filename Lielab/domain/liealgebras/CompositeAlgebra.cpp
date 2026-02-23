#include "CompositeAlgebra.hpp"

#include "Lielab/domain/liealgebras.hpp"

#include "Lielab/testing.hpp"

#include <Eigen/Core>

#include <initializer_list>
#include <numeric>
#include <variant>

namespace Lielab::domain
{

CompositeAlgebra::CompositeAlgebra()
{
    /*! \f{equation*}{() \rightarrow \mathfrak{CompositeAlgebra} \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::CompositeAlgebra x, y, z;
    * 
    */
}

CompositeAlgebra::CompositeAlgebra(const CompositeAlgebra::matrix_t& matrix)
{
    this->point.push_back(glc(matrix));
}

CompositeAlgebra CompositeAlgebra::basis(const int index, const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{CompositeAlgebra} \f}
    *
    * Returns the i'th basis element of the CompositeAlgebra algebra.
    * 
    * @param[in] index The basis vector.
    * @param[in] shape The shape of the algebra.
    * @param[out] out The CompositeAlgebra element.
    */

    return CompositeAlgebra({glc::basis(index, shape)});
}

CompositeAlgebra CompositeAlgebra::zero(const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{CompositeAlgebra} \f}
    *
    * Returns a zero algebra with given shape.
    * 
    * @param[in] shape The shape of the algebra.
    * @param[out] out The CompositeAlgebra element. 
    */

    return CompositeAlgebra({glc::zero(shape)});
}

CompositeAlgebra CompositeAlgebra::project(const CompositeAlgebra::matrix_t& matrix)
{
    return CompositeAlgebra({glc::project(matrix)});
}

CompositeAlgebra::CompositeAlgebra(const int n)
{
    /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{CompositeAlgebra} \f}
    *
    * Constructor instantiating a \f$\mathfrak{CompositeAlgebra}\f$ as multiple empty glc objects.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::CompositeAlgebra x(3), y(4), z(5);
    * 
    * @param[in] n
    */

    this->point.reserve(n);
    for (int ii = 0; ii < n; ii++)
    {
        this->point.push_back(glc(0));
    }
}



CompositeAlgebra::CompositeAlgebra(std::initializer_list<CompositeAlgebra::TYPES> others)
{
    /*!
    * Initializer list constructor for CompositeAlgebra
    *
    * Enables construction like:
    *     Lielab::domain::rn R(3);
    *     Lielab::domain::so O(3);
    *     Lielab::domain::CompositeAlgebra M{R, O};
    */

    this->point = std::vector<CompositeAlgebra::TYPES>{std::move(others)};
}

CompositeAlgebra::CompositeAlgebra(const std::vector<CompositeAlgebra::TYPES>& others)
{
    /*!
    * Vector list constructor for CompositeAlgebra
    *
    * Not needed for C++, but enables construction in Python like:
    *     R = lielab.domain.rn(3)
    *     O = lielab.domain.so(3)
    *     M = lielab.domain.CompositeAlgebra([R, O])
    */

    this->point = others;
}

std::string CompositeAlgebra::to_string() const
{
    std::string out = "";
    const int sz = static_cast<int>(this->point.size());

    int ii = 0;
    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            out += _element.to_string();
        }, element);

        if (ii < sz - 1) out += " ⊕ ";
        ii++;
    }

    return out;
}

int CompositeAlgebra::get_dimension() const
{
    const std::vector<int> dims = this->get_dimensions();
    return std::accumulate(dims.begin(), dims.end(), 0);
}

int CompositeAlgebra::get_size() const
{
    const std::vector<int> sizes = this->get_sizes();
    return std::accumulate(sizes.begin(), sizes.end(), 0);
}

bool CompositeAlgebra::is_abelian() const
{
    bool abelian = true;
    for (auto& element : this->point)
    {
        std::visit([&](auto& _element)
        {
            if (!_element.is_abelian()) abelian = false;
        }, element);
    }
    return abelian;
}

int CompositeAlgebra::get_shape() const
{
    const std::vector<int> shapes = this->get_shapes();
    return std::accumulate(shapes.begin(), shapes.end(), 0);
}

std::vector<CompositeAlgebra::TYPES>::iterator CompositeAlgebra::begin()
{
    return this->point.begin();
}

std::vector<CompositeAlgebra::TYPES>::iterator CompositeAlgebra::end()
{
    return this->point.end();
}

std::vector<CompositeAlgebra::TYPES>::const_iterator CompositeAlgebra::begin() const
{
    return this->point.begin();
}

std::vector<CompositeAlgebra::TYPES>::const_iterator CompositeAlgebra::end() const
{
    return this->point.end();
}

std::vector<int> CompositeAlgebra::get_dimensions() const
{
    std::vector<int> dimensions(this->point.size());

    int ii = 0;
    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            dimensions[ii] = _element.get_dimension();
        }, element);
        ii++;
    }

    return dimensions;
}

std::vector<int> CompositeAlgebra::get_sizes() const
{
    std::vector<int> sizes(this->point.size());

    int ii = 0;
    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            sizes[ii] = _element.get_size();
        }, element);
        ii++;
    }

    return sizes;
}

std::vector<int> CompositeAlgebra::get_shapes() const
{
    std::vector<int> shapes(this->point.size());

    int ii = 0;
    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            shapes[ii] = _element.get_shape();
        }, element);
        ii++;
    }

    return shapes;
}

CompositeAlgebra::point_t CompositeAlgebra::get_point() const
{
    return this->point;
}

Eigen::VectorXd CompositeAlgebra::serialize() const
{
    return this->get_vector();
}

void CompositeAlgebra::unserialize(const Eigen::VectorXd& serialized)
{
    this->set_vector(serialized);
}

void CompositeAlgebra::unserialize(std::initializer_list<double> serialized)
{
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
}

CompositeAlgebra::matrix_t CompositeAlgebra::get_matrix() const
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
    
    const int shape = this->get_shape();
    Eigen::MatrixXcd out = Eigen::MatrixXcd::Zero(shape, shape);
    
    int start = 0;
    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            const int _shape = _element.get_shape();
            out(Eigen::seqN(start, _shape), Eigen::seqN(start, _shape)) = _element.get_matrix();
            start += _shape;
        }, element);
    }

    return out;
}

Eigen::VectorXd CompositeAlgebra::get_vector() const
{
    return Lielab::utils::concatenate(this->get_vectors());
}

void CompositeAlgebra::set_vector(const Eigen::VectorXd& vector)
{
    // TODO: Error check for lengths here?

    int start = 0;
    for (auto& element : this->point)
    {
        std::visit([&](auto& _element)
        {
            const int dim = _element.get_dimension();
            _element.set_vector(vector(Eigen::seqN(start, dim)));
            start += dim;
        }, element);
    }
}

void CompositeAlgebra::set_vector(std::initializer_list<double> vector)
{
    /*!
    *
    * @param[in] vector
    */
   
    this->set_vector(Eigen::VectorXd{std::move(vector)});
}

double CompositeAlgebra::operator()(const int index) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the vector representation.
    */
    
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const int length = static_cast<int>(this->point.size());

    if (length == 0) return nan;

    const int dimension = this->get_dimension();

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? dimension + index : index;

    // Error check for out of bounds
    if (_index < 0) return nan;
    if (_index >= dimension) return nan;
    
    const std::vector<int> dimensions = this->get_dimensions();

    int base_index = 0;
    for (int ii = 0; ii < length; ii++)
    {
        if ((_index - base_index) < dimensions[ii])
        {
            return std::visit([&](const auto& _element)
            {
                return _element(_index - base_index);
            }, this->point[ii]);
        }

        base_index += dimensions[ii];
    }

    // This should never be returned.
    return nan;
}

CompositeAlgebra::field_t CompositeAlgebra::operator()(const int index1, const int index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{C} \f}
    *
    * Gets a value in the square matrix representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const int shape = this->get_shape();
    if (shape == 0) return std::complex<double>(nan, nan);

    // If input index is negative, index from the back of the array
    const int _index1 = (index1 < 0) ? shape + index1 : index1;
    const int _index2 = (index2 < 0) ? shape + index2 : index2;

    // Error check for out of bounds
    if (_index1 < 0) return std::complex<double>(nan, nan);
    if (_index1 >= shape) return std::complex<double>(nan, nan);
    if (_index2 < 0) return std::complex<double>(nan, nan);
    if (_index2 >= shape) return std::complex<double>(nan, nan);
    
    const std::vector<int> shapes = this->get_shapes();

    int relind = 0;
    for (int ii = 0; ii < static_cast<int>(shapes.size()); ii++)
    {
        // Sparse components
        if ((_index1 - relind) < 0 || (_index2 - relind) < 0)
        {
            return std::complex<double>(0.0, 0.0);
        }

        if ((_index1 - relind) < shapes[ii] && (_index2 - relind) < shapes[ii])
        {
            return std::visit([&](const auto& _element)
            {
                return static_cast<std::complex<double>>(_element(_index1 - relind, _index2 - relind));
            }, this->point[ii]);
        }

        relind += shapes[ii];
    }

    // This should never be called.
    return std::complex<double>(nan, nan);
}

std::vector<Eigen::VectorXd> CompositeAlgebra::get_vectors() const
{
    std::vector<Eigen::VectorXd> vectors = std::vector<Eigen::VectorXd>(this->point.size());

    int ii = 0;
    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            vectors[ii] = _element.get_vector();
        }, element);
        ii++;
    }

    return vectors;
}

const CompositeAlgebra::dataproxy CompositeAlgebra::operator[](const int index) const
{
    const int len = static_cast<int>(this->point.size());

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? len + index : index;

    // Error check for out of bounds
    lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for CompositeAlgebra of length " + std::to_string(len));

    return CompositeAlgebra::dataproxy{const_cast<CompositeAlgebra::TYPES&>(this->point[_index])};
}

CompositeAlgebra::dataproxy CompositeAlgebra::operator[](const int index)
{
    const int len = static_cast<int>(this->point.size());

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? len + index : index;

    // Error check for out of bounds
    lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for CompositeAlgebra of length " + std::to_string(len));

    return CompositeAlgebra::dataproxy{this->point[_index]};
}

CompositeAlgebra CompositeAlgebra::operator+(const CompositeAlgebra& other) const
{
    using Lielab::testing::check_topology;

    lielab_assert(check_topology(*this, other), "Unable to add topologically inconsistent algebras: (" + this->to_string() + ") !≅ (" + other.to_string() + ").");

    // const int length = this->point.size();

    CompositeAlgebra out;
    int index = 0;
    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            using other_t = std::decay_t<decltype(_element)>;
            out.point.emplace_back(_element + std::get<other_t>(other.point[index]));
        }, element);
        index++;
    }

    return out;
}

CompositeAlgebra& CompositeAlgebra::operator+=(const CompositeAlgebra& other)
{
    using Lielab::testing::check_topology;

    lielab_assert(check_topology(*this, other), "Unable to add topologically inconsistent algebras: (" + this->to_string() + ") !≅ (" + other.to_string() + ").");

    int index = 0;
    for (auto& element : this->point)
    {
        std::visit([&](auto& _element)
        {
            using other_t = std::decay_t<decltype(_element)>;
            _element += std::get<other_t>(other.point[index]);
        }, element);
        index++;
    }

    return *this;
}

CompositeAlgebra CompositeAlgebra::operator-(const CompositeAlgebra& other) const
{
    using Lielab::testing::check_topology;

    lielab_assert(check_topology(*this, other), "Unable to subtract topologically inconsistent algebras: (" + this->to_string() + ") !≅ (" + other.to_string() + ").");

    CompositeAlgebra out;
    int index = 0;
    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            using other_t = std::decay_t<decltype(_element)>;
            out.point.push_back(_element - std::get<other_t>(other.point[index]));
        }, element);
        index++;
    }

    return out;
}

CompositeAlgebra& CompositeAlgebra::operator-=(const CompositeAlgebra& other)
{
    using Lielab::testing::check_topology;

    lielab_assert(check_topology(*this, other), "Unable to subtract topologically inconsistent algebras: (" + this->to_string() + ") !≅ (" + other.to_string() + ").");

    int index = 0;
    for (auto& element : this->point)
    {
        std::visit([&](auto& _element)
        {
            using other_t = std::decay_t<decltype(_element)>;
            _element -= std::get<other_t>(other.point[index]);
        }, element);
        index++;
    }

    return *this;
}

CompositeAlgebra CompositeAlgebra::operator-() const
{
    CompositeAlgebra out;

    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            out.point.push_back(-_element);
        }, element);
    }

    return out;
}

CompositeAlgebra CompositeAlgebra::operator*(const double other) const
{
    CompositeAlgebra out;

    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            out.point.push_back(_element*other);
        }, element);
    }

    return out;
}

CompositeAlgebra operator*(const double other, const CompositeAlgebra& rhs)
{
    return rhs*other;
}

CompositeAlgebra& CompositeAlgebra::operator*=(const double other)
{
    for (auto& element : this->point)
    {
        std::visit([&](auto& _element)
        {
            _element *= other;
        }, element);
    }

    return *this;
}

CompositeAlgebra CompositeAlgebra::operator/(const double other) const
{
    CompositeAlgebra out;

    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            out.point.push_back(_element / other);
        }, element);
    }

    return out;
}

CompositeAlgebra& CompositeAlgebra::operator/=(const double other)
{
    for (auto& element : this->point)
    {
        std::visit([&](auto& _element)
        {
            _element /= other;
        }, element);
    }

    return *this;
}

}
