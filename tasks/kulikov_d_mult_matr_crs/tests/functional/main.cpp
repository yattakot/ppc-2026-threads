#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "kulikov_d_mult_matr_crs/common/include/common.hpp"
#include "kulikov_d_mult_matr_crs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kulikov_d_mult_matr_crs {

namespace {
CRSMatrix CreateCRSMatrix(int rows, int cols, int nnz, const std::vector<double> &values,
                          const std::vector<int> &col_ind, const std::vector<int> &row_ptr) {
  CRSMatrix matrix;
  matrix.rows = rows;
  matrix.cols = cols;
  matrix.nnz = nnz;
  matrix.values = values;
  matrix.col_ind = col_ind;
  matrix.row_ptr = row_ptr;
  return matrix;
}
}  // namespace

class KulikovDMultMatrCrsRunFuncTestsSEQ : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    const int test_case = std::get<0>(params);

    switch (test_case) {
      case 1: {
        const CRSMatrix a = CreateCRSMatrix(3, 3, 5, {2.5, 1.0, 3.0, 1.5, 2.0}, {0, 2, 1, 0, 1}, {0, 2, 3, 5});
        const CRSMatrix b = CreateCRSMatrix(3, 3, 5, {1.0, 2.0, 4.0, 3.0, 1.0}, {0, 2, 1, 0, 2}, {0, 2, 3, 5});

        expected_output_ = CreateCRSMatrix(3, 3, 6, {5.5, 6.0, 12.0, 1.5, 8.0, 3.0}, {0, 2, 1, 0, 1, 2}, {0, 2, 3, 6});
        input_data_ = std::make_tuple(a, b);
        break;
      }

      case 2: {
        const CRSMatrix a = CreateCRSMatrix(2, 4, 4, {1.0, 2.0, 3.0, 4.0}, {0, 2, 1, 3}, {0, 2, 4});

        const CRSMatrix b = CreateCRSMatrix(4, 3, 4, {1.0, 2.0, 3.0, 4.0}, {1, 0, 2, 1}, {0, 1, 2, 3, 4});

        expected_output_ = CreateCRSMatrix(2, 3, 4, {1.0, 6.0, 6.0, 16.0}, {1, 2, 0, 1}, {0, 2, 4});
        input_data_ = std::make_tuple(a, b);
        break;
      }

      case 3: {
        const CRSMatrix a = CreateCRSMatrix(2, 2, 1, {5.0}, {0}, {0, 1, 1});

        const CRSMatrix b = CreateCRSMatrix(2, 2, 1, {7.0}, {1}, {0, 0, 1});

        expected_output_ = CreateCRSMatrix(2, 2, 0, {}, {}, {0, 0, 0});
        input_data_ = std::make_tuple(a, b);
        break;
      }

      case 4: {
        const CRSMatrix a = CreateCRSMatrix(3, 3, 5, {0.5, 1.5, 2.0, 3.0, 4.0}, {0, 2, 0, 1, 2}, {0, 2, 3, 5});
        const CRSMatrix b = CreateCRSMatrix(3, 3, 3, {1.0, 1.0, 1.0}, {0, 1, 2}, {0, 1, 2, 3});

        expected_output_ = a;
        input_data_ = std::make_tuple(a, b);
        break;
      }

      case 5: {
        const CRSMatrix a = CreateCRSMatrix(4, 4, 4, {10.0, 20.0, 30.0, 40.0}, {0, 1, 2, 3}, {0, 1, 2, 3, 4});
        const CRSMatrix b = CreateCRSMatrix(4, 4, 4, {0.1, 0.2, 0.3, 0.4}, {0, 1, 2, 3}, {0, 1, 2, 3, 4});

        expected_output_ = CreateCRSMatrix(4, 4, 4, {1.0, 4.0, 9.0, 16.0}, {0, 1, 2, 3}, {0, 1, 2, 3, 4});
        input_data_ = std::make_tuple(a, b);
        break;
      }

      default:
        break;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.rows != expected_output_.rows || output_data.cols != expected_output_.cols ||
        output_data.nnz != expected_output_.nnz) {
      return false;
    }

    if (output_data.values.size() != expected_output_.values.size() ||
        output_data.col_ind.size() != expected_output_.col_ind.size() ||
        output_data.row_ptr.size() != expected_output_.row_ptr.size()) {
      return false;
    }

    for (size_t i = 0; i < output_data.row_ptr.size(); ++i) {
      if (output_data.row_ptr[i] != expected_output_.row_ptr[i]) {
        return false;
      }
    }

    const double epsilon = 1e-10;
    for (size_t i = 0; i < output_data.values.size(); ++i) {
      if (std::abs(output_data.values[i] - expected_output_.values[i]) > epsilon) {
        return false;
      }
      if (output_data.col_ind[i] != expected_output_.col_ind[i]) {
        return false;
      }
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  CRSMatrix expected_output_;
};

namespace {

TEST_P(KulikovDMultMatrCrsRunFuncTestsSEQ, SparseMatrixMultiplicationCRS) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParam = {
    std::make_tuple(1, "simple_3x3_"), std::make_tuple(2, "rectangular_2x4_4x3"), std::make_tuple(3, "zero_result_"),
    std::make_tuple(4, "identity_mult_"), std::make_tuple(5, "diagonal_scaled")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<KulikovDMultMatrCrsSEQ, InType>(kTestParam, PPC_SETTINGS_kulikov_d_mult_matr_crs);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = KulikovDMultMatrCrsRunFuncTestsSEQ::PrintFuncTestName<KulikovDMultMatrCrsRunFuncTestsSEQ>;

INSTANTIATE_TEST_SUITE_P(SparseMatrixMultCRS, KulikovDMultMatrCrsRunFuncTestsSEQ, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kulikov_d_mult_matr_crs
