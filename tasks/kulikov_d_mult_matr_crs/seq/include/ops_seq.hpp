#pragma once

#include "example_threads/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kulikov_d_mult_matr_crs {

class KulikovDMultMatrCrsSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KulikovDMultMatrCrsSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace kulikov_d_mult_matr_crs
