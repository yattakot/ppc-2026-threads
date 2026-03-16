#include "kulikov_d_mult_matr_crs/seq/include/ops_seq.hpp"

#include <unordered_map>
#include <algorithm>
#include <cmath>

#include "kulikov_d_mult_matr_crs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace kulikov_d_mult_matr_crs {

KulikovDMultMatrCrsSEQ::KulikovDMultMatrCrsSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CRSMatrix{0, 0, 0, {}, {}, {}};
}

bool KulikovDMultMatrCrsSEQ::ValidationImpl() {
  const auto& [matA, matB] = GetInput();

  if (matA.nnz != static_cast<int>(matA.values.size()) ||
      matB.nnz != static_cast<int>(matB.values.size())) {
    return false;
  }

  if (matA.cols != matB.rows) {
    return false;
  }

  if (!IsCRSValid(matA) || !IsCRSValid(matB)) {
    return false;
  }

  return true;
}

bool KulikovDMultMatrCrsSEQ::PreProcessingImpl() {
  return true;
}

bool KulikovDMultMatrCrsSEQ::RunImpl() {
  const auto& [matA, matB] = GetInput();

  try {
    GetOutput() = MultiplyCRS(matA, matB);
  } catch (...) {
    return false;
  }

  return !GetOutput().values.empty() || GetOutput().nnz == 0;
}

bool KulikovDMultMatrCrsSEQ::PostProcessingImpl() {
  const auto& result = GetOutput();
  return result.row_ptr.size() == static_cast<size_t>(result.rows + 1);
}

bool KulikovDMultMatrCrsSEQ::IsCRSValid(const CRSMatrix& mat) {
  if (mat.rows < 0 || mat.cols < 0 || mat.nnz < 0) return false;
  if (mat.row_ptr.size() != static_cast<size_t>(mat.rows + 1)) return false;
  if (mat.values.size() != static_cast<size_t>(mat.nnz)) return false;
  if (mat.col_ind.size() != static_cast<size_t>(mat.nnz)) return false;
  if (mat.row_ptr[0] != 0) return false;
  if (mat.row_ptr.back() != mat.nnz) return false;

  for (int i = 1; i <= mat.rows; ++i) {
    if (mat.row_ptr[i] < mat.row_ptr[i-1]) return false;
  }

  for (int col : mat.col_ind) {
    if (col < 0 || col >= mat.cols) return false;
  }

  return true;
}

CRSMatrix KulikovDMultMatrCrsSEQ::MultiplyCRS(const CRSMatrix& A,
                                               const CRSMatrix& B) {
  CRSMatrix C;
  C.rows = A.rows;
  C.cols = B.cols;
  C.row_ptr.push_back(0);

  for (int i = 0; i < A.rows; ++i) {
    std::unordered_map<int, double> row_accum;

    for (int idx_a = A.row_ptr[i]; idx_a < A.row_ptr[i + 1]; ++idx_a) {
      int col_a = A.col_ind[idx_a];
      double val_a = A.values[idx_a];

      for (int idx_b = B.row_ptr[col_a]; idx_b < B.row_ptr[col_a + 1]; ++idx_b) {
        int col_b = B.col_ind[idx_b];
        double val_b = B.values[idx_b];

        row_accum[col_b] += val_a * val_b;
      }
    }

    std::vector<std::pair<int, double>> sorted_elems;
    for (const auto& [col, val] : row_accum) {
      if (std::abs(val) > 1e-15) {
        sorted_elems.emplace_back(col, val);
      }
    }

    std::sort(sorted_elems.begin(), sorted_elems.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    for (const auto& [col, val] : sorted_elems) {
      C.col_ind.push_back(col);
      C.values.push_back(val);
    }

    C.row_ptr.push_back(static_cast<int>(C.values.size()));
  }

  C.nnz = static_cast<int>(C.values.size());

  return C;
}

}  // namespace kulikov_d_mult_matr_crs
