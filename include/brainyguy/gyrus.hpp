#pragma once
#ifndef BRAINYGUY_GYRUS_H
#define BRAINYGUY_GYRUS_H 1

/*
 *   Copyright 2022 Carlos Reyes
 *
 *   Licensed under the Apache License, Version 2.0 (the "License");
 *   you may not use this file except in compliance with the License.
 *   You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *   See the License for the specific language governing permissions and
 *   limitations under the License.
 */

#include <array>
#include <memory>
#include <ostream>
#include <tuple>
#include <cassert>
#include <cfloat>

// -------------------------------------------------------------------
constexpr bool approx_eq(const double a, const double b)
{
  return (std::abs(a - b) < (DBL_EPSILON * std::abs(a + b)));
}

// -------------------------------------------------------------------
constexpr bool is_power_of_two(const size_t n)
{
  return n && (!(n & (n-1)));
}

// -------------------------------------------------------------------
template<typename Scalar, size_t rows, size_t cols>
class Matrix
{
public:
  Matrix() : _elements{std::make_unique<std::array<Scalar, rows*cols>>()} {
    static_assert(is_power_of_two(rows));
    static_assert(is_power_of_two(cols));
  }

  Scalar get(const size_t row, const size_t col) const {
    return (*_elements)[row+col*cols];
  }

  Scalar operator()(const size_t row, const size_t col) const {
    return (*_elements)[row+col*cols];
  }

  Scalar& operator()(const size_t row, const size_t col) {
    return (*_elements)[row+col*cols];
  }

  static Matrix<Scalar, rows, cols> identity() {
    static_assert(rows == cols);
    Matrix<Scalar, rows, cols> result;
    for (size_t diag = 0; diag < rows; ++diag) {
      result(diag, diag) = 1.0;
    }
    return result;
  }

  static Matrix<Scalar, rows, cols>
  zeros() {
    Matrix<Scalar, rows, cols> result;
    return result;
  }

  static Matrix<Scalar, rows, cols>
  ones() {
    Matrix<Scalar, rows, cols> result;
    for (size_t row = 0; row < rows; ++row) {
      for (size_t col = 0; col < cols; ++col) {
	result(row, col) = 1.0;
      }
    }
    return result;
  }

  template<size_t r, size_t c>
  Matrix<Scalar, r, c>
  block(const size_t start_row, const size_t start_col) const {
    static_assert(r <= rows);
    static_assert(c <= cols);
    assert(start_row >= 0 && start_row < rows);
    assert(start_col >= 0 && start_col < cols);
    assert((start_row+r) <= rows);
    assert((start_col+c) <= cols);

    Matrix<Scalar, r, c> result;
    for (size_t row = 0; row < r; ++row) {
      for (size_t col = 0; col < c; ++col) {
	result(row, col) = get(start_row+row, start_col+col);
      }
    }
    return result;
  }

private:
  std::unique_ptr<std::array<Scalar, rows*cols>> _elements;
};

// -------------------------------------------------------------------
template<typename Scalar, size_t rows, size_t cols>
std::ostream&
operator<<(std::ostream& os, const Matrix<Scalar, rows, cols>& matrix)
{
  os << '[';
  for (size_t row = 0; row < rows; ++row) {
    os << '[';
    for (size_t col = 0; col < cols; ++col) {
      if (col != 0) {
	os << ',';
      }
      os << matrix(row, col);
    }
    os << ']';
  }

  os << ']' << std::endl;
  return os;
}

// -------------------------------------------------------------------
template<typename Scalar>
Scalar
determinant(const Matrix<Scalar, 1, 1>& matrix)
{
  return matrix(0, 0);
}

// -------------------------------------------------------------------
template<typename Scalar, size_t rows, size_t cols>
std::tuple<Matrix<Scalar, rows/2, cols/2>, Matrix<Scalar, rows/2, cols/2>,
  Matrix<Scalar, rows/2, cols/2>, Matrix<Scalar, rows/2, cols/2>>
four_way_split(const Matrix<Scalar, rows, cols>& matrix)
{
  static_assert(rows == cols);
  constexpr size_t r = rows/2;
  constexpr size_t c = cols/2;
  const Matrix<Scalar, r, c> A = matrix.template block<r, c>(0,      0);
  const Matrix<Scalar, r, c> B = matrix.template block<r, c>(0,      cols/2);
  const Matrix<Scalar, r, c> C = matrix.template block<r, c>(rows/2, 0);
  const Matrix<Scalar, r, c> D = matrix.template block<r, c>(rows/2, cols/2);
  return std::make_tuple(A, B, C, D);
}

// -------------------------------------------------------------------
template<typename Scalar, size_t rows, size_t cols>
Matrix<Scalar, rows, cols>
four_way_join(const std::tuple<Matrix<Scalar, rows/2, cols/2>, Matrix<Scalar, rows/2, cols/2>,
	      Matrix<Scalar, rows/2, cols/2>, Matrix<Scalar, rows/2, cols/2>>& matrices)
{
  static_assert(rows == cols);
  constexpr size_t r = rows/2;
  constexpr size_t c = cols/2;
  Matrix<Scalar, rows, cols> result;

  for (size_t row = 0; row < r; ++row) {
    for (size_t col = 0; col < c; ++col) {
      result(  row,   col) = std::get<0>(matrices)(row, col);
      result(  row, c+col) = std::get<1>(matrices)(row, col);
      result(r+row,   col) = std::get<2>(matrices)(row, col);
      result(r+row, c+col) = std::get<3>(matrices)(row, col);
    }
  }

  return result;
}

// -------------------------------------------------------------------
template<typename Scalar>
Scalar
determinant(const Matrix<Scalar, 2, 2>& matrix)
{
  return matrix(0, 0)*matrix(1, 1) - matrix(0, 1)*matrix(1, 0);
}

// -------------------------------------------------------------------
template<typename Scalar>
Matrix<Scalar, 1, 1>
inverse(const Matrix<Scalar, 1, 1>& matrix)
{
  Matrix<Scalar, 1, 1> result;
  result(0,0) = matrix(0, 0);
  return result;
}

// -------------------------------------------------------------------
template<typename Scalar>
Matrix<Scalar, 2, 2>
inverse(const Matrix<Scalar, 2, 2>& matrix)
{
  Matrix<Scalar, 2, 2> result;
  const Scalar inv_det = 1.0 / determinant(matrix);
  result(0, 0) = inv_det *  matrix(1, 1);
  result(0, 1) = inv_det * -matrix(0, 1);
  result(1, 0) = inv_det * -matrix(1, 0);
  result(1, 1) = inv_det *  matrix(0, 0);
  return result;
}

// -------------------------------------------------------------------
template<typename Scalar, size_t rows, size_t cols>
Matrix<Scalar, rows, cols>
inverse(const Matrix<Scalar, rows, cols>& matrix)
{
  auto [A, B, C, D] = four_way_split(matrix);



  Matrix<Scalar, rows, cols> result;
  return result;
}

// -------------------------------------------------------------------
template<typename Scalar, size_t rows1, size_t cols1, size_t rows2, size_t cols2>
Matrix<Scalar, rows1, cols2>
operator*(const Matrix<Scalar, rows1, cols1>& matrix1,
	  const Matrix<Scalar, rows2, cols2>& matrix2) {
  static_assert(cols1 == rows2);
  Matrix<Scalar, rows1, cols2> result;

  for (size_t row = 0; row < rows1; ++row) {
    for (size_t col = 0; col < cols2; ++col) {
      result(row, col) = matrix1(row, 0) * matrix2(0, col);
      for (size_t cell = 1; cell < cols1; ++cell) {
	result(row, col) += matrix1(row, cell) * matrix2(cell, col);
      }
    }
  }
}

#endif // BRAINYGUY_GYRUS_H
