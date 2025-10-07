#ifndef LIELAB_DOMAIN_GLC_HPP
#define LIELAB_DOMAIN_GLC_HPP

#include "LieGroup.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

class GLC : public LieGroup<std::complex<double>>
{
    /*!
    * The GLC class.
    */
    public:
    static constexpr bool abelian = false;
    size_t _shape = 0;

    std::string to_string() const;

    // Initialization methods

    GLC();

    GLC(const size_t shape);
    static GLC from_shape(const size_t shape);

    template<typename OtherDerived>
    GLC(const Eigen::MatrixBase<OtherDerived> & other);

    template<typename OtherDerived>
    GLC & operator=(const Eigen::MatrixBase<OtherDerived> & other);

    static Eigen::MatrixXcd project(const Eigen::MatrixXcd & other);

    size_t get_dimension() const;
    size_t get_shape() const;

    size_t get_size() const;

    Eigen::MatrixXcd get_matrix() const;

    GLC inverse() const;

    // Data representation

    Eigen::VectorXd serialize() const override;

    void unserialize(const Eigen::VectorXd &vector) override;
    void unserialize(std::initializer_list<double> vector);

    std::complex<double> operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;

    GLC operator*(const GLC & other) const;

    GLC & operator*=(const GLC & other);

    friend std::ostream & operator<<(std::ostream& os, const GLC & other);
};

}

#include "GLC.tpp"

#endif
