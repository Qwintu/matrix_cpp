#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : rows_(0), cols_(0) { matrix_ = nullptr; }

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows < 1 || cols < 1) {
    throw std::out_of_range("Wrong size of matrix!");
  }
  matrix_ = new double* [rows_] {};
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]{};
  };
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : S21Matrix(other.rows_, other.cols_) {
  for (int i{}; i < rows_; i++) {
    for (int j{}; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix&& other) : S21Matrix() {
  rows_ = other.rows_, cols_ = other.cols_;
  matrix_ = other.matrix_;
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

// S21Matrix::S21Matrix(S21Matrix &&other) noexcept
//     : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
//   other.rows_ = 0;
//   other.cols_ = 0;
//   other.matrix_ = nullptr;
//     std::cout << "constr move" << std::endl;
// }

// S21Matrix::S21Matrix(S21Matrix&& other) : S21Matrix() {
//   *this = std::move(other);
//   other.rows_ = 0;
//   other.cols_ = 0;
//   other.matrix_ = nullptr;
// }

S21Matrix::~S21Matrix() { DelMatrix(); };

bool S21Matrix::EqMatrix(const S21Matrix& other) {
  double eps = 0.000001;
  bool error_code = true;

  if ((rows_ != other.rows_) || (cols_ != other.cols_)) return false;
  if (this == &other) return true;

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (fabs(matrix_[i][j] - other.matrix_[i][j]) >= eps) {
        error_code = false;
        break;
      }
    }
    if (error_code == false) break;
  }
  return error_code;
}

void S21Matrix::MulNumber(const double num) {
  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] = matrix_[i][j] * num;
      }
    }
  }
}

void S21Matrix::SumMatrix(
    const S21Matrix& other) {  // различная размерность матриц
  if (matrix_ == nullptr || other.matrix_ == nullptr)
    throw std::out_of_range("Invalid matrix");
  if ((rows_ != other.rows_) || (other.cols_ != cols_))
    throw std::out_of_range("Incorrect matrix size");

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] + other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(
    const S21Matrix& other) {  // различная размерность матриц
  if (matrix_ == nullptr || other.matrix_ == nullptr)
    throw std::out_of_range("Invalid matrix");
  if ((rows_ != other.rows_) || (other.cols_ != cols_))
    throw std::out_of_range("Incorrect matrix size");

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] - other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulMatrix(
    const S21Matrix&
        other) {  // Умножает текущую матрицу на вторую / число столбцов первой
                  // матрицы не равно числу строк второй матрицы
  if (matrix_ == nullptr || other.matrix_ == nullptr)
    throw std::out_of_range("Invalid matrix");
  if ((cols_ != other.rows_)) throw std::out_of_range("Incorrect matrix size");

  S21Matrix temp(rows_, other.cols_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      double value_of_calc = 0;
      for (int n = 0; n < cols_; n++) {
        value_of_calc += matrix_[i][n] * other.matrix_[n][j];
      }
      temp.matrix_[i][j] = value_of_calc;
    }
  }
  *this = temp;
}

S21Matrix S21Matrix::Transpose() {
  if (matrix_ == nullptr) throw std::out_of_range("Invalid matrix");

  S21Matrix tmp(cols_, rows_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      tmp.matrix_[j][i] = matrix_[i][j];
    }
  }
  return tmp;
}

double
S21Matrix::Determinant() {  // Вычисляет и возвращает определитель текущей
                            // матрицы / матрица не является квадратной
  if (matrix_ == nullptr || cols_ != rows_)
    throw std::out_of_range("Invalid matrix");
  double det = 0;
  if (rows_ <= 2 && rows_ > 0) {
    det = determinant_2_1();
  } else if (rows_ > 2) {
    for (int j = 0; j < cols_; j++) {
      double minor_det = 0;
      S21Matrix minor(rows_ - 1, cols_ - 1);
      minor = create_minor(0, j);
      minor_det = minor.Determinant();
      det += pow((-1), (0 + j)) * matrix_[0][j] * minor_det;
    }
  }

  return det;
}

S21Matrix S21Matrix::CalcComplements() {
  if (matrix_ == nullptr || cols_ != rows_)
    throw std::out_of_range("Invalid matrix");

  S21Matrix ComplMatrix(rows_, cols_);

  if (rows_ == 1) {
    ComplMatrix.matrix_[0][0] = matrix_[0][0];
  } else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        double minor_det = 0;
        S21Matrix minor(rows_ - 1, cols_ - 1);
        minor = create_minor(i, j);
        minor_det = minor.Determinant();
        minor_det *= pow((-1), (i + j));
        ComplMatrix.matrix_[i][j] = minor_det;
      }
    }
  }
  return ComplMatrix;
}

S21Matrix S21Matrix::InverseMatrix() {
  double det = this->Determinant();
  if (fabs(det - 0.0) < 1e-6)
    throw std::out_of_range("The determinant shouldn't be equals zero");

  S21Matrix InvMatr(rows_, cols_);
  InvMatr = CalcComplements();
  InvMatr = InvMatr.Transpose();
  InvMatr.MulNumber(1.0 / det);
  return InvMatr;
}

double S21Matrix::determinant_2_1() {
  double det = 0;

  if (rows_ == 1) {
    det = matrix_[0][0];
  } else if (rows_ == 2) {
    // определитель матрицы 2*2
    det = matrix_[0][0] * matrix_[1][1] - matrix_[1][0] * matrix_[0][1];
  }
  return det;
}

S21Matrix S21Matrix::create_minor(int i, int j) {
  S21Matrix minor(rows_ - 1, cols_ - 1);
  int minor_row = 0;
  int minor_columns = 0;

  for (int k = 0; k < rows_; k++) {
    if (k == i) continue;
    minor_columns = 0;
    for (int m = 0; m < cols_; m++) {
      if (m != j) {
        minor.matrix_[minor_row][minor_columns] = matrix_[k][m];
        minor_columns++;
      }
    }
    minor_row++;
  }
  return minor;
}

// _____overloaded oper_____________________

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this != &other) {
    SetRows(other.rows_);
    SetCols(other.cols_);
    rows_ = other.rows_;
    cols_ = other.cols_;
    if (other.matrix_ != nullptr) {
      for (int i = 0; i < other.rows_; i++) {
        for (int j = 0; j < other.cols_; j++) {
          matrix_[i][j] = other.matrix_[i][j];
        }
      }
    } else {
      matrix_ = nullptr;
    }
  }
  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) {
  if (this != &other) {
    DelMatrix();
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = other.matrix_;
    other.rows_ = 0;
    other.cols_ = 0;
    other.matrix_ = nullptr;
  }
  return *this;
}

S21Matrix S21Matrix::operator+(
    const S21Matrix& other) {  //    +   Сложение двух матриц  /  различная
                               //    размерность матриц
  S21Matrix tmp = *this;
  tmp.SumMatrix(other);
  return tmp;
}

S21Matrix S21Matrix::operator-(
    const S21Matrix& other) {  //    -   Вычитание одной матрицы из другой  /
                               //    различная размерность матриц
  S21Matrix tmp = *this;
  tmp.SubMatrix(other);
  return tmp;
}

S21Matrix S21Matrix::operator*(
    const S21Matrix& other) {  //    *   Умножение матриц и умножение матрицы на
                               //    число /  число столбцов первой матрицы не
                               //    равно числу строк второй матрицы
  S21Matrix tmp = *this;
  tmp.MulMatrix(other);
  return tmp;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix tmp = *this;
  tmp.MulNumber(num);
  return tmp;
}

bool S21Matrix::operator==(
    const S21Matrix&
        other) {  //    ==  Проверка на равенство матриц (EqMatrix) / -
  return EqMatrix(other);
}

S21Matrix& S21Matrix::operator+=(
    const S21Matrix& other) {  //    +=  Присвоение сложения (SumMatrix) /
                               //    различная размерность матриц
  this->SumMatrix(other);
  return *this;
}

void S21Matrix::operator-=(
    const S21Matrix& other) {  //    -=  Присвоение разности (SubMatrix) /
                               //    различная размерность матриц
  SubMatrix(other);
}

void S21Matrix::operator*=(
    const S21Matrix&
        other) {  //    *=  Присвоение умножения (MulMatrix/MulNumber) / число
                  //    столбцов первой матрицы не равно числу строк второй
                  //    матрицы
  MulMatrix(other);
}

void S21Matrix::operator*=(const double num) { MulNumber(num); }

double& S21Matrix::operator()(int row,
                              int col) {  // индекс за пределами матрицы
  if (matrix_ == nullptr) throw std::out_of_range("Invalid matrix");
  if (row < 1 || col < 1 || row > rows_ || col > cols_) {
    throw std::out_of_range("Wrong matrix index");
  }
  return matrix_[row][col];
}

// _____add func_____________________

int S21Matrix::GetRows() { return rows_; }
int S21Matrix::GetCols() { return cols_; }

void S21Matrix::SetRows(int row) {
  if (row <= 0) throw std::length_error("Wrong rows quantity");
  S21Matrix temp(row, cols_);
  for (int i = 0; i < std::min(rows_, row); i++) {
    for (int j = 0; j < cols_; j++) {
      temp.matrix_[i][j] = matrix_[i][j];
    }
  }
  DelMatrix();
  *this = std::move(temp);
}

void S21Matrix::SetCols(int cols) {
  if (cols <= 0) throw std::length_error("Wrong cols quantity");
  S21Matrix temp(rows_, cols);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < std::min(cols_, cols); j++) {
      temp.matrix_[i][j] = matrix_[i][j];
    }
  }
  DelMatrix();
  *this = std::move(temp);
}

void S21Matrix::DelMatrix() {
  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    };
    delete[] matrix_;
    matrix_ = nullptr;
  }
}

// _____serv func_____________________

// void S21Matrix::print() {
//     std::cout << "Matr " << "r-" << rows_ << " c-"<< cols_ << std::endl;
//     for(int i{}; i < rows_; i++){
//         for(int j{}; j < cols_; j++){
//             std::cout << matrix_[i][j] << "\t";
//         }
//         std::cout << std::endl;
//     }
// }

void S21Matrix::fill_matrix() {
  for (int i{}; i < rows_; i++) {
    for (int j{}; j < cols_; j++) {
      matrix_[i][j] = (i + j);
    }
  }
}

void S21Matrix::fill_matrix_by_massive(int rows, int columns, double* values) {
  int mass_ind = 0;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      matrix_[i][j] = values[mass_ind];
      mass_ind++;
      ;
    }
  }
}

void S21Matrix::change_val(int i, int j, double v) { matrix_[i][j] = v; }
