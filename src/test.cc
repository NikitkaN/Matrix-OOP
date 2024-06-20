#include <gtest/gtest.h>

#include "S21_Matrix_oop.h"

TEST(test1, EqMatrix) {
  School21::S21Matrix m1(2, 2);
  School21::S21Matrix m2(2, 2);
  School21::S21Matrix m3(4, 4);
  m1.fillMatrix();
  m2.fillMatrix();
  EXPECT_TRUE(m1.EqMatrix(m2));
  EXPECT_FALSE(m1.EqMatrix(m3));

  m1.setNum(3.0000007, 0, 0);
  m2.setNum(3.0000007, 0, 0);
  EXPECT_TRUE(m1.EqMatrix(m2));

  m2.setNum(3.0000001, 0, 0);
  EXPECT_FALSE(m1.EqMatrix(m2));
}

TEST(test2, SumMatrix) {
  School21::S21Matrix m1(6, 6);
  School21::S21Matrix m2(6, 6);
  m1.fillMatrix();
  m2.fillMatrix();
  m1.SumMatrix(m2);
  m2.MulNumber(2);
  EXPECT_TRUE(m1 == m2);
}

TEST(test3, SubMatrix) {
  School21::S21Matrix m1(4, 4);
  School21::S21Matrix m2(4, 4);
  m1.fillMatrix();
  m2.fillMatrix();
  m1.SubMatrix(m2);
  School21::S21Matrix m3(4, 4);
  EXPECT_TRUE(m1 == m3);
}

TEST(test4, MulNumber) {
  School21::S21Matrix m1(2, 4);
  m1.fillMatrix();
  m1.MulNumber(2);
  double tmp[2][4] = {{0, 2, 4, 6}, {8, 10, 12, 14}};
  for (int i = 0; i < m1.getRows(); i++) {
    for (int j = 0; j < m1.getCols(); j++) {
      EXPECT_DOUBLE_EQ(m1(i, j), tmp[i][j]);
    }
  }
}

TEST(test5, MulMatrix) {
  School21::S21Matrix m1(3, 3);
  School21::S21Matrix m2(3, 3);
  m1.fillMatrix();
  m2.fillMatrix();
  m1.MulMatrix(m2);
  School21::S21Matrix m3(3, 3);
  m3.setNum(15, 0, 0);
  m3.setNum(18, 0, 1);
  m3.setNum(21, 0, 2);
  m3.setNum(42, 1, 0);
  m3.setNum(54, 1, 1);
  m3.setNum(66, 1, 2);
  m3.setNum(69, 2, 0);
  m3.setNum(90, 2, 1);
  m3.setNum(111, 2, 2);

  for (int i = 0; i < m3.getRows(); i++) {
    for (int j = 0; j < m1.getCols(); j++) {
      EXPECT_DOUBLE_EQ(m1(i, j), m3(i, j));
    }
  }
}

TEST(test6, Transpose) {
  School21::S21Matrix m1(3, 3);
  m1.fillMatrix();
  School21::S21Matrix m2 = m1.Transpose();
  double tmp[3][3] = {{0, 3, 6}, {1, 4, 7}, {2, 5, 8}};
  for (int i = 0; i < m1.getRows(); i++) {
    for (int j = 0; j < m1.getCols(); j++) {
      EXPECT_DOUBLE_EQ(m2(i, j), tmp[i][j]);
    }
  }
}

TEST(test7, Determinat) {
  School21::S21Matrix m1(3, 3), m2(3, 3);
  m2.fillMatrix();
  double tmp = m2.Determinant();
  EXPECT_DOUBLE_EQ(tmp, 0);

  m1.setNum(98, 0, 0);
  m1.setNum(52, 0, 1);
  m1.setNum(50, 0, 2);
  m1.setNum(30, 1, 0);
  m1.setNum(18, 1, 1);
  m1.setNum(34, 1, 2);
  m1.setNum(1, 2, 0);
  m1.setNum(2, 2, 1);
  m1.setNum(3, 2, 2);
  tmp = m1.Determinant();
  EXPECT_DOUBLE_EQ(tmp, -2184);
}

TEST(test8, CalcComplements) {
  School21::S21Matrix m1(3, 3);
  m1.fillMatrix();
  School21::S21Matrix m2 = m1.CalcComplements();
  double tmp[3][3] = {{-3, 6, -3}, {6, -12, 6}, {-3, 6, -3}};
  for (int i = 0; i < m2.getRows(); i++) {
    for (int j = 0; j < m2.getCols(); j++) {
      EXPECT_DOUBLE_EQ(m2(i, j), tmp[i][j]);
    }
  }
}

TEST(test9, InverseMatrix) {
  School21::S21Matrix m1(3, 3);
  m1.setNum(68, 0, 0);
  m1.setNum(52, 0, 1);
  m1.setNum(40, 0, 2);
  m1.setNum(23, 1, 0);
  m1.setNum(18, 1, 1);
  m1.setNum(14, 1, 2);
  m1.setNum(34, 2, 0);
  m1.setNum(2, 2, 1);
  m1.setNum(8, 2, 2);
  School21::S21Matrix m2 = m1.InverseMatrix();
  School21::S21Matrix m3(3, 3);
  m3.setNum(0.006410, 0, 0);
  m3.setNum(0.025641, 0, 1);
  m3.setNum(-0.397436, 0, 2);
  m3.setNum(0.025641, 1, 0);
  m3.setNum(-0.111722, 1, 1);
  m3.setNum(0.838828, 1, 2);
  m3.setNum(-0.019231, 2, 0);
  m3.setNum(0.065934, 2, 1);
  m3.setNum(-0.093407, 2, 2);
  for (int i = 0; i < m2.getRows(); i++) {
    for (int j = 0; j < m2.getCols(); j++) {
      double val_1 = m2(i, j);
      double val_2 = m3(i, j);
      if (val_1 == val_2) {
        EXPECT_TRUE(1);
      } else {
        EXPECT_FALSE(false);
      }
    }
  }
}

TEST(test10, Sum_operator) {
  School21::S21Matrix m1(4, 4);
  School21::S21Matrix m2(4, 4);
  School21::S21Matrix m3;
  m1.fillMatrix();
  m2.fillMatrix();
  m3 = m1 + m2;
  m2 *= 2;
  EXPECT_TRUE(m3 == m2);
}

TEST(test11, Sub_operator) {
  School21::S21Matrix m1(5, 5);
  School21::S21Matrix m2(5, 5);
  School21::S21Matrix m3;
  m1.fillMatrix();
  m2.fillMatrix();
  m3 = m1 - m2;
  m2.zeroMatrix();
  EXPECT_TRUE(m3 == m2);
}

TEST(test12, operators1) {
  School21::S21Matrix m1(4, 4);
  School21::S21Matrix m2(4, 4);
  School21::S21Matrix m3(4, 4);
  m1.fillMatrix();
  m2.fillMatrix();
  m3.fillMatrix();
  m1 += m2;
  m3 *= 2;
  EXPECT_TRUE(m1 == m3);
  m3 -= m1;
  m1.zeroMatrix();
  EXPECT_TRUE(m1 == m3);
}

TEST(test_13, operators2) {
  School21::S21Matrix m1(3, 3);
  School21::S21Matrix m2(3, 3);
  School21::S21Matrix m3(3, 3);
  m1.fillMatrix();
  m2.fillMatrix();
  School21::S21Matrix m4 = m1 * m2;
  m3.setNum(15, 0, 0);
  m3.setNum(18, 0, 1);
  m3.setNum(21, 0, 2);
  m3.setNum(42, 1, 0);
  m3.setNum(54, 1, 1);
  m3.setNum(66, 1, 2);
  m3.setNum(69, 2, 0);
  m3.setNum(90, 2, 1);
  m3.setNum(111, 2, 2);
  EXPECT_TRUE(m3 == m4);

  m1 *= m2;
  EXPECT_TRUE(m1 == m4);

  School21::S21Matrix m5(m2);
  m5 = m5 * 2;
  m2 += m2;
  EXPECT_TRUE(m2 == m5);
}

TEST(test14, Move) {
  School21::S21Matrix m1(4, 4);
  m1.fillMatrix();
  School21::S21Matrix m2 = std::move(m1);
  School21::S21Matrix m3(m2);
  EXPECT_TRUE(m2 == m3);
}

TEST(test15, Set) {
  double val1[2][3] = {{0, 1, 2}, {3, 4, 5}};
  School21::S21Matrix m1(3, 3);
  m1.fillMatrix();
  m1.setRows(2);
  for (int i = 0; i < m1.getRows(); i++) {
    for (int j = 0; j < m1.getCols(); j++) {
      EXPECT_EQ(m1(i, j), val1[i][j]);
    }
  }

  double val2[4][3] = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}, {0, 0, 0}};
  School21::S21Matrix m2(3, 3);
  m2.fillMatrix();
  m2.setRows(4);
  for (int i = 0; i < m2.getRows(); i++) {
    for (int j = 0; j < m2.getCols(); j++) {
      EXPECT_EQ(m2(i, j), val2[i][j]);
    }
  }

  double val3[2][1] = {{0}, {2}};
  School21::S21Matrix m3(2, 2);
  m3.fillMatrix();
  m3.setCols(1);
  for (int i = 0; i < m3.getRows(); i++) {
    for (int j = 0; j < m3.getCols(); j++) {
      EXPECT_EQ(m3(i, j), val3[i][j]);
    }
  }

  double val4[2][3] = {{0, 1, 0}, {2, 3, 0}};
  School21::S21Matrix m4(2, 2);
  m4.fillMatrix();
  m4.setCols(3);
  for (int i = 0; i < m4.getRows(); i++) {
    for (int j = 0; j < m4.getCols(); j++) {
      EXPECT_EQ(m4(i, j), val4[i][j]);
    }
  }

  School21::S21Matrix m5(3, 3);
  m5.fillMatrix();
  double **test_body1 = m5.getMatrix_();
  EXPECT_EQ(test_body1[1][1], 4);

  m5.setNum(0, 1, 1);
  double **test_body2 = m5.getMatrix_();
  EXPECT_EQ(test_body2[1][1], 0);
}

int main(int argc, char *argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}