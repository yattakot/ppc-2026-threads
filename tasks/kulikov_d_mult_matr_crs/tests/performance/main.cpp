#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <tuple>
#include <vector>

#include "kulikov_d_mult_matr_crs/common/include/common.hpp"
#include "kulikov_d_mult_matr_crs/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kulikov_d_mult_matr_crs {

class KulikovDMultMatrCrsRunPerfTestSEQ : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kRowsA_ = 100;        // строки матрицы A
  const int kColsA_ = 100;        // столбцы A (= строки B)
  const int kColsB_ = 100;        // столбцы матрицы B
  const double kDensity_ = 0.05;  // плотность ненулевых элементов (5%)
  const int kRandomSeed_ = 42;    // фиксированный seed для воспроизводимости

  InType input_data_{};

  CRSMatrix GenerateSparseCRS(int rows, int cols, double density, int seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> val_dist(-10.0, 10.0);
    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

    CRSMatrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.row_ptr.push_back(0);

    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        if (prob_dist(gen) < density) {
          mat.values.push_back(val_dist(gen));
          mat.col_ind.push_back(j);
        }
      }
      mat.row_ptr.push_back(static_cast<int>(mat.values.size()));
    }

    mat.nnz = static_cast<int>(mat.values.size());
    return mat;
  }

  void SetUp() override {
    CRSMatrix matA = GenerateSparseCRS(kRowsA_, kColsA_, kDensity_, kRandomSeed_);
    CRSMatrix matB = GenerateSparseCRS(kColsA_, kColsB_, kDensity_, kRandomSeed_ + 1);

    input_data_ = std::make_tuple(matA, matB);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const auto &[matA, matB] = input_data_;

    if (output_data.rows != matA.rows || output_data.cols != matB.cols) {
      return false;
    }

    if (output_data.row_ptr.size() != static_cast<size_t>(output_data.rows + 1)) {
      return false;
    }

    if (output_data.nnz != static_cast<int>(output_data.values.size()) ||
        output_data.nnz != static_cast<int>(output_data.col_ind.size())) {
      return false;
    }

    if (output_data.row_ptr[0] != 0 || output_data.row_ptr.back() != output_data.nnz) {
      return false;
    }

    for (int i = 1; i <= output_data.rows; ++i) {
      if (output_data.row_ptr[i] < output_data.row_ptr[i - 1]) {
        return false;
      }
    }

    for (int col : output_data.col_ind) {
      if (col < 0 || col >= output_data.cols) {
        return false;
      }
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KulikovDMultMatrCrsRunPerfTestSEQ, RunPerfSEQ) {
  ExecuteTest(GetParam());
}

namespace {

const auto kSeqPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KulikovDMultMatrCrsSEQ>(PPC_SETTINGS_kulikov_d_mult_matr_crs);

const auto kGtestValues = ppc::util::TupleToGTestValues(kSeqPerfTasks);

const auto kPerfTestName = KulikovDMultMatrCrsRunPerfTestSEQ::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(SEQ_PerfTests, KulikovDMultMatrCrsRunPerfTestSEQ, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kulikov_d_mult_matr_crs
