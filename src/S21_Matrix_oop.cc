#include "S21_Matrix_oop.h"

School21::S21Matrix::S21Matrix() {
  rows_ = 0;
  cols_ = 0;
  this->createMatrix();
}

School21::S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  this->createMatrix();
}

School21::S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  this->createMatrix();
  this->copyObj(other);
}

School21::S21Matrix::S21Matrix(S21Matrix&& other) {
  if (this != &other) {
    this->matrix_ = other.matrix_;
    this->rows_ = other.rows_;
    this->cols_ = other.cols_;
    other.matrix_ = NULL;
    other.rows_ = 0;
    other.cols_ = 0;
  }
}

School21::S21Matrix::~S21Matrix() {
  this->removeMatrix();
  this->rows_ = 0;
  this->cols_ = 0;
}

void School21::S21Matrix::createMatrix() {
  this->matrix_ = new double*[this->rows_];
  for (int i = 0; i < this->rows_; i++) {
    this->matrix_[i] = new double[this->cols_]();
  }
}

void School21::S21Matrix::removeMatrix() {
  if (this->matrix_) {
    for (int i = 0; i < this->rows_; i++) {
      delete this->matrix_[i];
    }
    delete[] this->matrix_;
  }
}

void School21::S21Matrix::copyMatrix(const S21Matrix& other, int i) {
  if (i) {
    for (int i = 0; i < other.rows_; i++) {
      for (int j = 0; j < other.cols_; j++) {
        this->matrix_[i][j] = other.matrix_[i][j];
      }
    }
  } else {
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++) {
        this->matrix_[i][j] = other.matrix_[i][j];
      }
    }
  }
}

void School21::S21Matrix::copyObj(const S21Matrix& other) {
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

int School21::S21Matrix::getRows() const { return rows_; }

int School21::S21Matrix::getCols() const { return cols_; }

double** School21::S21Matrix::getMatrix_() const { return matrix_; }

void School21::S21Matrix::setRows(int rows) {
  int i = 0;
  if (rows > this->rows_) {
    i = 1;
  }
  S21Matrix A = S21Matrix(*this);
  this->removeMatrix();
  this->rows_ = rows;
  this->createMatrix();
  this->copyMatrix(A, i);
}

void School21::S21Matrix::setCols(int cols) {
  int i = 0;
  if (cols > this->cols_) {
    i = 1;
  }
  S21Matrix A = S21Matrix(*this);
  this->removeMatrix();
  this->cols_ = cols;
  this->createMatrix();
  this->copyMatrix(A, i);
}

void School21::S21Matrix::setNum(double num, int rows, int cols) {
  if (this->rows_ >= rows || this->cols_ >= cols) {
    matrix_[rows][cols] = num;
  }
}

bool School21::S21Matrix::EqMatrix(const S21Matrix& other) {
  bool result = true;

  if (this->rows_ != other.rows_ || this->cols_ != other.cols_) {
    result = false;
  } else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (fabs(matrix_[i][j] - other.matrix_[i][j]) >
            1e-7 * fabs(other.matrix_[i][j])) {
          result = false;
          break;
        }
        if (!result) break;
      }
    }
  }
  return result;
}

void School21::S21Matrix::SumMatrix(const S21Matrix& other) {
  CheckMatrix(other);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void School21::S21Matrix::SubMatrix(const S21Matrix& other) {
  CheckMatrix(other);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void School21::S21Matrix::MulNumber(const double num) {
  for (auto i = 0; i < rows_; i++) {
    for (auto j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void School21::S21Matrix::MulMatrix(const S21Matrix& other) {
  CheckMulMatrix(other);
  S21Matrix res = S21Matrix(this->rows_, other.cols_);
  for (int i = 0; i < res.rows_; i++) {
    for (int j = 0; j < res.cols_; j++) {
      for (int n = 0; n < this->cols_; n++) {
        res.matrix_[i][j] += this->matrix_[i][n] * other.matrix_[n][j];
      }
    }
  }
  *this = res;
}

School21::S21Matrix School21::S21Matrix::Transpose() {
  S21Matrix trans = S21Matrix(cols_, rows_);
  for (int i = 0; i < trans.rows_; i++) {
    for (int j = 0; j < trans.cols_; j++) {
      trans.matrix_[i][j] = matrix_[j][i];
    }
  }

  return trans;
}

double School21::S21Matrix::Determinant() {
  CheckDeter();
  double res = 0;

  if (rows_ == 2 && cols_ == 2) {
    res = MFinder();
  } else {
    for (int j = 0; j < cols_; j++) {
      S21Matrix minor_m(rows_ - 1, cols_ - 1);
      Cutter(0, j, minor_m);
      int m = j % 2 == 0 ? 1 : -1;
      if (minor_m.rows_ == 2) {
        res += m * matrix_[0][j] * minor_m.MFinder();
      } else {
        res += m * matrix_[0][j] * minor_m.Determinant();
      }
    }
  }
  return res;
}

double School21::S21Matrix::MFinder() const {
  return matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
}

void School21::S21Matrix::Cutter(int i_val, int j_val, S21Matrix& minor) const {
  int n = 0, m = 0;
  for (int i = 0; i < rows_; i++) {
    if (i_val != i) {
      for (int j = 0; j < cols_; j++) {
        if (j_val != j) {
          minor.matrix_[n][m] = matrix_[i][j];
          m++;
        }
      }
      m = 0;
      n++;
    }
  }
}

School21::S21Matrix School21::S21Matrix::CalcComplements() {
  CheckDeter();
  S21Matrix res = S21Matrix(rows_, cols_);
  for (int i = 0; i < res.rows_; i++) {
    for (int j = 0; j < res.cols_; j++) {
      S21Matrix minor_m(rows_ - 1, cols_ - 1);
      Cutter(i, j, minor_m);
      if (minor_m.rows_ >= 2 && minor_m.cols_ >= 2) {
        res.matrix_[i][j] = minor_m.Determinant() * std::pow(-1, i + j);
      } else if (minor_m.rows_ == 1 && minor_m.cols_ == 1) {
        res.matrix_[i][j] = minor_m.matrix_[0][0] * std::pow(-1, i + j);
      }
    }
  }

  return res;
}

School21::S21Matrix School21::S21Matrix::InverseMatrix() {
  double res = Determinant();
  if (res == 0) {
    throw std::out_of_range("Determinant == 0");
  } else {
    CalcComplements();
    Transpose();
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] = matrix_[i][j] / res;
      }
    }
    return *this;
  }
}

School21::S21Matrix School21::S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix res(*this);
  res.SumMatrix(other);
  return res;
}

School21::S21Matrix& School21::S21Matrix::operator+=(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

School21::S21Matrix School21::S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix res(*this);
  res.SubMatrix(other);
  return res;
}

School21::S21Matrix& School21::S21Matrix::operator-=(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

School21::S21Matrix& School21::S21Matrix::operator=(const S21Matrix& other) {
  if (this->EqMatrix(other)) {
    return *this;
  }
  this->removeMatrix();
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  this->createMatrix();
  this->copyObj(other);
  return *this;
}

School21::S21Matrix School21::S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix res(*this);
  res.MulMatrix(other);
  return res;
}

School21::S21Matrix& School21::S21Matrix::operator*=(const S21Matrix& other) {
  this->MulMatrix(other);
  return *this;
}

School21::S21Matrix School21::S21Matrix::operator*(const double num) {
  S21Matrix res(*this);
  res.MulNumber(num);
  return res;
}

School21::S21Matrix& School21::S21Matrix::operator*=(const double num) {
  this->MulNumber(num);
  return *this;
}

bool School21::S21Matrix::operator==(const S21Matrix& other) {
  bool res = this->EqMatrix(other);
  return res;
}

double& School21::S21Matrix::operator()(int row, int col) {
  if (row < 0 || col < 0 || row > rows_ || col > cols_) {
    throw out_of_range("Index is outside the matrix");
  }
  return matrix_[row][col];
}

void School21::S21Matrix::CheckZero() const {
  if (rows_ <= 0 || cols_ <= 0) {
    throw runtime_error("Rows or cols <= 0");
  }
}

void School21::S21Matrix::CheckMatrix(const S21Matrix& other) const {
  CheckZero();
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw runtime_error("Wrong Matrix");
  }
}

void School21::S21Matrix::CheckMulMatrix(const S21Matrix& other) const {
  CheckZero();
  other.CheckZero();
  if (cols_ != other.rows_) {
    throw runtime_error("Rows of the second matrix");
  }
}

void School21::S21Matrix::CheckDeter() const {
  CheckZero();
  if (cols_ != rows_) {
    throw runtime_error("Need a square matrix");
  }
}

void School21::S21Matrix::fillMatrix() {
  int num = 0;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = num++;
    }
  }
}

void School21::S21Matrix::zeroMatrix() {
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] = 0;
    }
  }
}