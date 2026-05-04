#pragma once

#include "kulikov_d_mult_matr_crs/common/include/common.hpp"

namespace kulikov_d_mult_matr_crs {

class KulikovDMultMatrCrsSEQ : public BaseTask {
 public:
  explicit KulikovDMultMatrCrsSEQ(const InType &in);

  static std::string GetStaticTypeOfTask() {
    return "kulikov_d_mult_matr_crs_seq";
  }

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static bool IsCRSValid(const CRSMatrix &mat);
  static CRSMatrix MultiplyCRS(const CRSMatrix &a, const CRSMatrix &b);
};

}  // namespace kulikov_d_mult_matr_crs
