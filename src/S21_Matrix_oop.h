#ifndef SRC_S21_MATRIX_OOP_H_
#define SRC_S21_MATRIX_OOP_H_

#include <cmath>
#include <cstring>
#include <iostream>
#include <typeinfo>
using namespace std;

namespace School21 {
class S21Matrix {
 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);
  ~S21Matrix();

  bool EqMatrix(const S21Matrix& other);
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  double Determinant();
  S21Matrix CalcComplements();
  S21Matrix InverseMatrix();

  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const double num);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix& operator*=(const double num);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  double& operator()(int row, int col);
  bool operator==(const S21Matrix& other);

  int getRows() const;
  int getCols() const;
  double** getMatrix_() const;

  void setRows(int rows);
  void setCols(int cols);
  void setNum(double num, int rows, int cols);

  void fillMatrix();
  void zeroMatrix();

 private:
  int rows_, cols_;
  double** matrix_;

  void createMatrix();
  void removeMatrix();
  void copyMatrix(const S21Matrix& other, int i);
  void copyObj(const S21Matrix& other);
  void Cutter(int i_val, int j_val, S21Matrix& minor) const;
  double MFinder() const;
  void CheckMatrix(const S21Matrix& other) const;
  void CheckMulMatrix(const S21Matrix& other) const;
  void CheckDeter() const;
  void CheckZero() const;
};
}  // namespace School21

#endif  // SRC_S21_MATRIX_OOP_H_