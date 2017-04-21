// Class for Matrix calculations of arbitrary size, used in orcaplugin
// Header file

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>

class Matrix {

public:
    double **value;
    int columns, rows;

    Matrix();
    Matrix(std::string input);
    // template <typename T>
    Matrix(std::vector<std::vector<float>> input);
    ~Matrix();
    double** parseInput(std::string input);

    static double stringToDouble(std::string& s);
    static void printMatrix(Matrix *matrix);
    void printMatrix();
    static std::vector<double>* rowAtIndex(Matrix *input, unsigned int index);
    static std::vector<double>* columnAtIndex(Matrix *input, unsigned int index);
    static double dotProduct(std::vector<double> *firstVector, std::vector<double> *secondVector);
    static std::vector<double>* crossProduct(std::vector<double> *firstVector, std::vector<double> *secondVector);
    static double getAngle(std::vector<double> *firstVector,std::vector<double> *secondVector);
    static double getNorm(std::vector<double> *firstVector);
    static std::vector<double>* makeVector(std::string input);
    static Matrix* copyMatrix(Matrix *input);
    static Matrix* setRow(Matrix *input, std::vector<double> *rowValues, int row);
    static Matrix* setColumn(Matrix *input, std::vector<double> *columnValues, int column);
    static Matrix* switchRows(Matrix *input, int firstRow, int secondRow);
    static Matrix* switchColumns(Matrix *input, int firstColumn, int secondColumn);
    static Matrix* multiplyRow(Matrix *input, int row, double factor);
    static Matrix* multiplyColumn(Matrix *input, int column, double factor);
    static Matrix* addRowToRow(Matrix *input, int addRow, int modifiedRow, double factor);
    static Matrix* addColumnToColumn(Matrix *input, int addColumn, int modifiedColumn, double factor);
    static Matrix* multiply(Matrix *firstMatrix, Matrix *secondMatrix);
    static Matrix* zeilenStufenForm(Matrix *input);
    static Matrix* triagonal(Matrix *input);
    static double determinante(Matrix *input);
    static int rang(Matrix *input);
    static void printVector(std::vector<double> *vec);
};
