#include "pairinteraction/interfaces/TransformationBuilderInterface.hpp"

#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/utils/euler.hpp"

namespace pairinteraction {

IndicesOfBlock::IndicesOfBlock(size_t start, size_t end) : start(start), end(end) {}

size_t IndicesOfBlock::size() const { return end - start; }

IndicesOfBlocks::IndicesOfBlocks(std::initializer_list<size_t> boundaries)
    : boundaries(boundaries) {}

void IndicesOfBlocks::add(size_t boundary) { boundaries.insert(boundary); }

std::vector<IndicesOfBlock> IndicesOfBlocks::get() const {
    std::vector<IndicesOfBlock> blocks;
    if (boundaries.empty()) {
        return blocks;
    }

    auto it = boundaries.begin();
    size_t start = *it++;

    while (it != boundaries.end()) {
        blocks.emplace_back(start, *it);
        start = *it++;
    }

    return blocks;
}

size_t IndicesOfBlocks::size() const { return boundaries.size(); }

template <typename Scalar>
Transformation<Scalar>::Transformation(Eigen::SparseMatrix<Scalar, Eigen::RowMajor> matrix,
                                       std::vector<TransformationType> transformation_type)
    : matrix(std::move(matrix)), transformation_type(std::move(transformation_type)) {}

template <typename Scalar>
Transformation<Scalar>::Transformation(Eigen::SparseMatrix<Scalar, Eigen::RowMajor> matrix)
    : matrix(std::move(matrix)), transformation_type({TransformationType::ARBITRARY}) {}

Sorting::Sorting(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> matrix,
                 std::vector<TransformationType> transformation_type)
    : matrix(std::move(matrix)), transformation_type(std::move(transformation_type)) {}

Sorting::Sorting(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> matrix)
    : matrix(std::move(matrix)), transformation_type({TransformationType::ARBITRARY}) {}

template <typename Scalar>
Transformation<Scalar>
TransformationBuilderInterface<Scalar>::get_rotator(const std::array<real_t, 3> &to_z_axis,
                                                    const std::array<real_t, 3> &to_y_axis) const {
    auto euler_zyz_angles = euler::get_euler_angles(to_z_axis, to_y_axis);
    return this->get_rotator(euler_zyz_angles[0], euler_zyz_angles[1], euler_zyz_angles[2]);
}

// Explicit instantiations
template struct Transformation<float>;
template struct Transformation<double>;
template struct Transformation<std::complex<float>>;
template struct Transformation<std::complex<double>>;

template class TransformationBuilderInterface<float>;
template class TransformationBuilderInterface<double>;
template class TransformationBuilderInterface<std::complex<float>>;
template class TransformationBuilderInterface<std::complex<double>>;
} // namespace pairinteraction
