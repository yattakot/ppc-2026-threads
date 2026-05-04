#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace kulikov_d_mult_matr_crs {

struct CRSMatrix {
  int rows{0};                 // количество строк
  int cols{0};                 // количество столбцов
  int nnz{0};                  // количество ненулевых элементов
  std::vector<double> values;  // ненулевые значения
  std::vector<int> col_ind;    // индексы столбцов
  std::vector<int> row_ptr;    // указатели на строки
};

using InType = std::tuple<CRSMatrix, CRSMatrix>;  // две входные матрицы
using OutType = CRSMatrix;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kulikov_d_mult_matr_crs
