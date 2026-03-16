#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "kulikov_d_mult_matr_crs/common/include/common.hpp"
#include "kulikov_d_mult_matr_crs/common/include/test_reader.hpp"
#include "kulikov_d_mult_matr_crs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kulikov_d_mult_matr_crs {

class KulikovDMultMatrCrsFuncTest
    : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return test_param;
  }

 protected:
  void SetUp() override {
    const auto &filename = ppc::util::GetAbsoluteTaskPath(
        std::string(PPC_ID_kulikov_d_mult_matr_crs),
        std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam()) + ".txt");

    const auto &[a, b, c] = ReadTestFromFile(filename);
    input_data_ = std::make_tuple(a, b);
    output_data_ = c;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return Equal(output_data_, output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType output_data_;
};

namespace {

TEST_P(KulikovDMultMatrCrsFuncTest, MatMulCRS) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 6> kTestParam = {
};

const auto kTestTasksList =
    std::tuple_cat(
        ppc::util::AddFuncTask<KulikovDMultMatrCrsSEQ, InType>(
            kTestParam, PPC_SETTINGS_kulikov_d_mult_matr_crs));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    KulikovDMultMatrCrsFuncTest::PrintFuncTestName<KulikovDMultMatrCrsFuncTest>;

INSTANTIATE_TEST_SUITE_P(
    CRSMatrixMul,
    KulikovDMultMatrCrsFuncTest,
    kGtestValues,
    kPerfTestName);

}  // namespace
}  // namespace kulikov_d_mult_matr_crs
