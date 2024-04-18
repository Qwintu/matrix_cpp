#include <gtest/gtest.h>
// #include <iostream>
#include "s21_matrix_oop.h"

class MatrixTest : public ::testing::Test {
  void SetUp() {
    A = S21Matrix(3, 3);
    A.fill_matrix();
    B = S21Matrix(3, 3);
    B.fill_matrix();
    C = S21Matrix(4, 4);
    C.fill_matrix();
  }
  void TearDown() {}

 protected:
  S21Matrix A, B, C, result;
};

TEST_F(MatrixTest, EqMatrix) {
  EXPECT_TRUE(A == A);
  B.change_val(0, 0, 1);
  EXPECT_FALSE(A == B);
  EXPECT_FALSE(A == C);
};

TEST_F(MatrixTest, Constr) {
  S21Matrix t;
  ASSERT_EQ(t.GetCols(), 0);
  ASSERT_EQ(t.GetRows(), 0);

  S21Matrix D(3, 3);
  result = S21Matrix(3, 3);
  double mass[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  result.fill_matrix_by_massive(3, 3, mass);
  EXPECT_TRUE(D == result);

  S21Matrix A2 = A;
  EXPECT_TRUE(A2 == A);

  S21Matrix A3 = std::move(A2);
  EXPECT_TRUE(A == A3);
  ASSERT_EQ(A2.GetCols(), 0);
  ASSERT_EQ(A2.GetRows(), 0);
};

TEST_F(MatrixTest, sum) {
  result = S21Matrix(3, 3);
  double mass[9] = {0, 1, 2, 1, 2, 3, 2, 3, 4};
  result.fill_matrix_by_massive(3, 3, mass);
  S21Matrix tmp0(3, 3);
  S21Matrix tmp1(3, 3);
  S21Matrix tmp2(3, 3);
  tmp0.SumMatrix(A);
  tmp1 = tmp1 + A;
  tmp2 += A;
  EXPECT_TRUE(tmp0 == result && tmp1 == result && tmp2 == result);
  EXPECT_ANY_THROW(C.SumMatrix(A));
};

TEST_F(MatrixTest, sub) {
  result = S21Matrix(3, 3);
  S21Matrix tmp0(3, 3);
  S21Matrix tmp1(3, 3);
  S21Matrix tmp2(3, 3);
  tmp0 = B;
  tmp1 = B;
  tmp2 = B;
  tmp0.SubMatrix(A);
  tmp1 = tmp1 - A;
  tmp2 -= A;
  EXPECT_TRUE(tmp0 == result && tmp1 == result && tmp2 == result);
  EXPECT_ANY_THROW(C.SubMatrix(A));
};

TEST_F(MatrixTest, MulNum) {
  result = S21Matrix(3, 3);
  double mass[9] = {0, 10, 20, 10, 20, 30, 20, 30, 40};
  result.fill_matrix_by_massive(3, 3, mass);
  A.MulNumber(10);
  EXPECT_TRUE(A == result);
};

TEST_F(MatrixTest, MulMatr) {
  result = S21Matrix(3, 4);
  double mass[12] = {14, 20, 26, 32, 20, 30, 40, 50, 26, 40, 54, 68};
  result.fill_matrix_by_massive(3, 4, mass);
  S21Matrix tmp0(3, 4);
  tmp0.fill_matrix();
  S21Matrix tmp1(3, 4);
  tmp1.fill_matrix();
  S21Matrix tmp2(3, 4);
  tmp2.fill_matrix();

  tmp0.MulMatrix(C);
  tmp1 = tmp1 * C;
  tmp2 *= C;
  EXPECT_TRUE(tmp0 == result && tmp1 == result && tmp2 == result);
  EXPECT_ANY_THROW(C.MulMatrix(A));
};

TEST_F(MatrixTest, MulMatrNum) {
  result = S21Matrix(3, 3);
  double mass[9] = {0, 10, 20, 10, 20, 30, 20, 30, 40};
  result.fill_matrix_by_massive(3, 3, mass);
  S21Matrix tmp0(3, 3);
  tmp0 = A * 10;
  A *= 10;
  EXPECT_TRUE(tmp0 == result && A == result);
};

TEST_F(MatrixTest, Transpose) {
  S21Matrix test(3, 4);
  test.fill_matrix();
  S21Matrix tmp(4, 3);
  tmp = test.Transpose();
  result = S21Matrix(4, 3);
  double mass[12] = {0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5};
  result.fill_matrix_by_massive(4, 3, mass);
  EXPECT_TRUE(tmp == result);
};

TEST_F(MatrixTest, RowsColls) {
  result = S21Matrix(4, 4);
  A.SetRows(4);
  A.SetCols(4);
  double mass[16] = {0, 1, 2, 0, 1, 2, 3, 0, 2, 3, 4, 0, 0, 0, 0, 0};
  result.fill_matrix_by_massive(4, 4, mass);
  EXPECT_TRUE(A == result);

  S21Matrix result_2(3, 3);
  A.SetRows(3);
  A.SetCols(3);
  double mass_1[9] = {0, 1, 2, 1, 2, 3, 2, 3, 4};
  result_2.fill_matrix_by_massive(3, 3, mass_1);
  EXPECT_TRUE(A == result_2);

  EXPECT_TRUE(A(1, 1) == 2);
  EXPECT_ANY_THROW(A(-1, 1));
  EXPECT_ANY_THROW(A(1, 10));
  EXPECT_TRUE(A.GetRows() == 3);
  EXPECT_TRUE(A.GetCols() == 3);
};

TEST_F(MatrixTest, det) {
  double det = 0;
  C.change_val(0, 0, 1);
  C.change_val(1, 2, 13);
  det = C.Determinant();
  EXPECT_TRUE(det == 20);

  S21Matrix tmp0(3, 4);
  EXPECT_ANY_THROW(tmp0.Determinant());
};

TEST_F(MatrixTest, CalcCompl) {
  S21Matrix tmp0(3, 3);
  double mass_1[9] = {1, 2, 3, 0, 4, 2, 5, 2, 1};
  A.fill_matrix_by_massive(3, 3, mass_1);
  result = S21Matrix(3, 3);
  double mass_2[9] = {0, 10, -20, 4, -14, 8, -8, -2, 4};
  result.fill_matrix_by_massive(3, 3, mass_2);
  tmp0 = A.CalcComplements();
  EXPECT_TRUE(tmp0 == result);

  S21Matrix tmp1(3, 4);
  EXPECT_ANY_THROW(tmp1.CalcComplements());
};

TEST_F(MatrixTest, Inverse) {
  S21Matrix tmp0(3, 3);
  double mass_1[9] = {3, 2, 5, 4, 5, 6, 7, 7, 7};
  A.fill_matrix_by_massive(3, 3, mass_1);
  result = S21Matrix(3, 3);
  double mass_2[9] = {0.25,      -0.75, 0.464286, -0.5, 0.5,
                      -0.071429, 0.25,  0.25,     -0.25};
  result.fill_matrix_by_massive(3, 3, mass_2);
  tmp0 = A.InverseMatrix();
  EXPECT_TRUE(tmp0 == result);

  S21Matrix tmp1(4, 4);
  tmp1.change_val(0, 0, 1);
  EXPECT_ANY_THROW(tmp1.InverseMatrix());
};

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}