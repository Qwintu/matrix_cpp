#ifndef S21_MATRIX_OOP_H
#define S21_MATRIX_OOP_H

#include <math.h>

#include <algorithm>
#include <stdexcept>

class S21Matrix {
 private:
  int rows_, cols_;  // Rows and columns
  double** matrix_;  // Pointer to the memory where the matrix is allocated

 protected:
  double determinant_2_1();
  S21Matrix create_minor(int i, int j);

 public:
  S21Matrix();  // Базовый конструктор, инициализирующий матрицу некоторой
                // заранее заданной размерностью
  S21Matrix(int rows, int cols);  // Параметризированный конструктор с
                                  // количеством строк и столбцов
  S21Matrix(const S21Matrix& other);  // Конструктор копирования
  S21Matrix(S21Matrix&& other);  // Конструктор переноса

  ~S21Matrix();  // Деструктор

  bool EqMatrix(
      const S21Matrix& other);  // Проверяет матрицы на равенство между собой
  void SumMatrix(
      const S21Matrix& other);  // Прибавляет вторую матрицы к текущей /
                                // различная размерность матриц
  void SubMatrix(const S21Matrix& other);  // Вычитает из текущей матрицы другую
                                           // / различная размерность матриц
  void MulNumber(const double num);  // Умножает текущую матрицу на число
  void MulMatrix(const S21Matrix& other);  // Умножает текущую матрицу на вторую
                                           // / число столбцов первой матрицы не
                                           // равно числу строк второй матрицы
  S21Matrix Transpose();  // Создает новую транспонированную матрицу из текущей
                          // и возвращает ее
  S21Matrix CalcComplements();  // Вычисляет матрицу алгебраических дополнений
                                // текущей матрицы и возвращает ее / матрица не
                                // является квадратной
  double Determinant();  // Вычисляет и возвращает определитель текущей матрицы
                         // / матрица не является квадратной
  S21Matrix InverseMatrix();  // Вычисляет и возвращает обратную матрицу /
                              // определитель матрицы равен 0

  // OVERLOADED operators

  S21Matrix& operator=(
      const S21Matrix&
          other);  //    =   Присвоение матрице значений другой матрицы / -
  S21Matrix& operator=(S21Matrix&& other);
  S21Matrix operator+(
      const S21Matrix& other);  //    +   Сложение двух матриц  /  различная
                                //    размерность матриц
  S21Matrix operator-(
      const S21Matrix& other);  //    -   Вычитание одной матрицы из другой  /
                                //    различная размерность матриц
  S21Matrix operator*(
      const S21Matrix& other);  //    *   Умножение матриц и умножение матрицы
                                //    на число /  число столбцов первой матрицы
                                //    не равно числу строк второй матрицы
  S21Matrix operator*(const double num);
  bool operator==(const S21Matrix& other);  //    ==  Проверка на равенство
                                            //    матриц (EqMatrix) / -
  S21Matrix& operator+=(
      const S21Matrix& other);  //    +=  Присвоение сложения (SumMatrix) /
                                //    различная размерность матриц
  void operator-=(
      const S21Matrix& other);  //    -=  Присвоение разности (SubMatrix) /
                                //    различная размерность матриц
  void operator*=(
      const S21Matrix&
          other);  //    *=  Присвоение умножения (MulMatrix/MulNumber) / число
                   //    столбцов первой матрицы не равно числу строк второй
                   //    матрицы
  void operator*=(const double num);
  double& operator()(
      int row, int col);  //    (int i, int j)  Индексация по элементам матрицы
                          //    (строка, колонка) /  индекс за пределами матрицы

  // Additional func

  int GetRows();
  int GetCols();
  void SetRows(int row);
  void SetCols(int cols);
  void DelMatrix();

  void print();
  void fill_matrix();
  void fill_matrix_by_massive(int rows, int columns, double* values);
  void change_val(int i, int j, double v);
};

#endif  // S21_MATRIX_OOP_H