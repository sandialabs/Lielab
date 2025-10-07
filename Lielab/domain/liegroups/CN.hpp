#ifndef LIELAB_DOMAIN_CN_HPP
#define LIELAB_DOMAIN_CN_HPP

#include "LieGroup.hpp"
#include "GLC.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

class CN : public GLC
{
    /*!
    * The CN class.
    */
    public:
    static constexpr bool abelian = true;

    using data_t = Eigen::VectorXcd;
    data_t data;

    std::string to_string() const;
    // Initialization methods

    CN();
    CN(const size_t n);
    static CN from_shape(const size_t shape);

    template<typename OtherDerived>
    CN(const Eigen::MatrixBase<OtherDerived> & other);

    template<typename OtherDerived>
    CN & operator=(const Eigen::MatrixBase<OtherDerived> & other);

    size_t get_dimension() const;
    size_t get_shape() const;
    size_t get_size() const;

    Eigen::MatrixXcd get_matrix() const;

    CN inverse() const;

    // Data representation

    Eigen::VectorXd serialize() const override;

    void unserialize(const Eigen::VectorXd& vector) override;
    void unserialize(const std::initializer_list<double> vector); // override;

    double operator()(const ptrdiff_t index) const;
    // std::complex<double> & operator()(const size_t index);
    std::complex<double> operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;

    CN operator*(const CN & other) const;
    CN & operator*=(const CN & other);

    std::complex<double> operator[](const ptrdiff_t index) const;

    static CN from_vector(const Eigen::VectorXd& other);
    static CN from_vector(std::initializer_list<double> other);
    static CN from_complex_vector(const Eigen::VectorXcd& other);
    static CN from_complex_vector(const std::initializer_list<std::complex<double>> other);
    Eigen::VectorXcd to_complex_vector() const;

    static Eigen::MatrixXcd project(const Eigen::MatrixXcd & other);

    friend std::ostream & operator<<(std::ostream& os, const CN & other);
};

}

#include "CN.tpp"

#endif
