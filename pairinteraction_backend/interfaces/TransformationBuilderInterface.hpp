#pragma once

#include "utils/eigen_assertion.hpp"
#include "utils/traits.hpp"

#include <Eigen/SparseCore>
#include <array>
#include <complex>
#include <vector>

enum class TransformationType : unsigned char;

template <typename Scalar>
struct Transformation {
    Transformation() = default;
    Transformation(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                   const std::vector<TransformationType> &transformation_type);
    Transformation(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix);
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> matrix;
    std::vector<TransformationType> transformation_type;
};

struct Sorting {
    Sorting() = default;
    Sorting(const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> &matrix,
            const std::vector<TransformationType> &transformation_type);
    Sorting(const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> &matrix);
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> matrix;
    std::vector<TransformationType> transformation_type;
};

struct Blocks {
    Blocks() = default;
    Blocks(const std::vector<int> &start,
           const std::vector<TransformationType> &transformation_type);
    Blocks(const std::vector<int> &start);
    std::vector<int> start;
    std::vector<TransformationType> transformation_type;
};

template <typename Scalar>
class TransformationBuilderInterface {
public:
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using real_t = typename traits::NumTraits<Scalar>::real_t;

    virtual const Transformation<Scalar> &get_transformation() const = 0;
    virtual Transformation<Scalar> get_rotator(real_t alpha, real_t beta, real_t gamma) const = 0;
    virtual Sorting get_sorter(TransformationType label) const = 0;
    virtual Blocks get_blocks(TransformationType label) const = 0;

    Transformation<Scalar> get_rotator(std::array<real_t, 3> to_z_axis,
                                       std::array<real_t, 3> to_y_axis) const;
};

extern template class Transformation<float>;
extern template class Transformation<double>;
extern template class Transformation<std::complex<float>>;
extern template class Transformation<std::complex<double>>;

extern template class TransformationBuilderInterface<float>;
extern template class TransformationBuilderInterface<double>;
extern template class TransformationBuilderInterface<std::complex<float>>;
extern template class TransformationBuilderInterface<std::complex<double>>;
