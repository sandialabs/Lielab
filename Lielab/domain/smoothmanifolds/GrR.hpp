#ifndef LIELAB_DOMAIN_GRR_HPP
#define LIELAB_DOMAIN_GRR_HPP

#include <Eigen/Core>

namespace Lielab::domain
{

class GrR
{
    public:

    Eigen::MatrixXd data;

    std::string to_string() const;

    GrR();
    GrR(const size_t k, const size_t n);
    
    template<typename OtherDerived>
    GrR(const Eigen::MatrixBase<OtherDerived>& other);

    template<typename OtherDerived>
    GrR& operator=(const Eigen::MatrixBase<OtherDerived>& other);

    Eigen::VectorXd serialize() const;

    size_t get_dimension() const;
    size_t get_size() const;

    double operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;

    void unserialize(const Eigen::VectorXd& vec);
    void unserialize(std::initializer_list<double> vec);

    Eigen::MatrixXd get_matrix() const;
    Eigen::VectorXd project_onto(const Eigen::VectorXd& vec) const;

};

}

#include "GrR.tpp"

#endif
