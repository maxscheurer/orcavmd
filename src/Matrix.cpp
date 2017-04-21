// Class for Matrix calculations of arbitrary size, used in orcaplugin
// Matrix.cpp
//

#include "Matrix.h"

using namespace std;

double Matrix::stringToDouble(std::string& s)
{
    istringstream i(s);
    double x;
    if (!(i >> x))
        return 0;
    return x;
}

double** Matrix::parseInput(std::string input)
{
    rows = 0;
    columns = 0;
    vector<double> tempValues;
    string tempString;
    bool firstBracket = false;
    bool closedBracket = false;
    char previousChar = '\0';

    for(char &c : input) {
        if (firstBracket && c == '{') {
            rows++;
        } else if(c == '{') {
            firstBracket = true;
        } else if(c == ',') {
            tempValues.push_back(stringToDouble(tempString));
            tempString.clear();
            if (closedBracket == false) {
                columns++;
            }
        } else if(c == '}' && closedBracket == false) {
            closedBracket = true;
            columns++;
        } else if (c != ',' && c != '{' && c != '}') {
            tempString.append(&c);
        } else if (c == '}' && previousChar == '}') {
            tempValues.push_back(stringToDouble(tempString));
        }
        previousChar = c;
    }

    double **matrix = 0;
    matrix = new double *[rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double [columns];
    }


    int n = 0;
    for (vector<double>::iterator i = tempValues.begin(); i < tempValues.end(); i++) {
        div_t mod = div(n, columns);
        matrix[mod.quot][mod.rem] = *i;
        n++;
    }


    return matrix;
}

Matrix::Matrix(std::string input)
{
    value = parseInput(input);
}

Matrix::Matrix()
{
    value = nullptr;
    rows = 0;
    columns = 0;
}

// template <typename T>
Matrix::Matrix(std::vector<std::vector<float>> input) {
    rows = input.size();
    columns = input[0].size();
    value = new double *[rows];
    for (size_t i = 0; i < rows; i++) {
        value[i] = new double [columns];
        if (columns != input[i].size()) {
            cout << "Matrix rows have different number of columns." << endl;
            return;
        }
    }
    int row = 0, col = 0;
    for (auto r : input) {
        for (auto v : r) {
            value[row][col] = v;
            col++;
        }
        row++;
        col = 0;
    }
}

Matrix::~Matrix()
{
    for (size_t i = 0; i < rows; i++) {
        delete[] value[i];
    }
    delete[] value;
}

std::vector<std::vector<float>> Matrix::toVector() {
    std::vector<std::vector<float>> result;
    for (size_t i = 0; i < rows; i++) {
        std::vector<float> tmpRow;
        for (size_t j = 0; j < columns; j++) {
            tmpRow.push_back(value[i][j]);
        }
        result.push_back(tmpRow);
    }
    return result;
}

Matrix* Matrix::copyMatrix(Matrix *input)
{
    double **matrix = 0;
    matrix = new double *[input->rows];
    for (int i = 0; i < input->rows; i++) {
        matrix[i] = new double [input->columns];
    }
    Matrix *result = new Matrix();
    result->value = matrix;
    result->rows = input->rows;
    result->columns = input->columns;

    for (int y = 0; y < input->rows; y++) {
        for (int x = 0; x < input->columns; x++) {
            result->value[y][x] = input->value[y][x];
        }
    }
    return result;
}

std::vector<double>* Matrix::rowAtIndex(Matrix *input, unsigned int index)
{
    vector<double> *result = new vector<double>;
    for (int i = 0; i < input->columns; i++) {
        result->push_back(input->value[index][i]);
    }
    return result;
}

std::vector<double>* Matrix::columnAtIndex(Matrix *input, unsigned int index)
{
    vector<double> *result = new vector<double>;
    for (int i = 0; i < input->rows; i++) {
        result->push_back(input->value[i][index]);
    }
    return result;
}

std::vector<double>* Matrix::makeVector(std::string input) {
    vector<double>* tempValues = new vector<double>;
    string tempString;
    bool firstIteration = true;
    char previousChar = '\0';
    for(char &c : input) {
        if((c == ',' && firstIteration) || (c == ',' && previousChar == ',')) {
            cout << "Syntax error.";
            break;
        }
        else if(c == ',' && !firstIteration && previousChar != ',') {
            tempValues->push_back(stringToDouble(tempString));
            tempString.clear();
        }
        else if(c!= ',') {
            tempString.append(&c);
        }
        previousChar = c;
        firstIteration = false;
    }
    tempValues->push_back(stringToDouble(tempString));
    tempString.clear();
    return tempValues;
}

double Matrix::dotProduct(std::vector<double> *firstVector, std::vector<double> *secondVector)
{
    double result = 0.0;
    if (firstVector == nullptr || secondVector == nullptr){
        cout << "Nullpointer Exception. \n";
    }
    else if (firstVector->size() == secondVector->size()) {
        for (int i = 0; i < firstVector->size(); i++) {
            result += firstVector->at(i) * secondVector->at(i);
        }
    }
    return result;
}

double Matrix::getNorm(std::vector<double> *firstVector){
    double result = 0.0;
    if (firstVector != nullptr) {
        double temp = 0.0;
        for (int i = 0; i < firstVector->size(); i++) {
            temp += pow(firstVector->at(i), 2);
        }
        result = sqrt(temp);
    }
    else {
        cout << "Nullpointer Exception.\n";
    }
    return result;
}

double Matrix::getAngle(std::vector<double> *firstVector,std::vector<double> *secondVector) {
    double result = 0.0;
    if (firstVector != nullptr && secondVector != nullptr) {
        result = acos((dotProduct(firstVector, secondVector))/(getNorm(firstVector)*getNorm(secondVector)));
    }
    else {
       cout << "Nullpointer Exception.\n";
    }
    return result;
}

std::vector<double>* Matrix::crossProduct(std::vector<double> *firstVector, std::vector<double> *secondVector){
    vector<double> *result = new vector<double>;
    if (firstVector->size() == 3 && secondVector->size() == 3) {
        double a,b,c;
        a = firstVector->at(1)*secondVector->at(2)-firstVector->at(2)*secondVector->at(1); // 23 - 32  31-13 12-21
        b = firstVector->at(2)*secondVector->at(0)-firstVector->at(0)*secondVector->at(2);
        c = firstVector->at(0)*secondVector->at(1)-firstVector->at(1)*secondVector->at(0);
        result->push_back(a);
        result->push_back(b);
        result->push_back(c);
    }
    else {
        cout << "Cross product not possible \n";
        return nullptr;
    }
    return result;
}

Matrix* Matrix::setRow(Matrix *input, std::vector<double> *rowValues, int row)
{
    Matrix* result = copyMatrix(input);
    for (int i = 0; i < rowValues->size(); i++) {
        result->value[row][i] = rowValues->at(i);
    }
    return result;
}

Matrix* Matrix::setColumn(Matrix *input, std::vector<double> *columnValues, int column)
{
    Matrix* result = copyMatrix(input);
    for (int i = 0; i < columnValues->size(); i++) {
        result->value[i][column] = columnValues->at(i);
    }
    return result;
}

Matrix* Matrix::switchRows(Matrix *input, int firstRow, int secondRow)
{
    Matrix *result = copyMatrix(input);
    vector<double> *tempRow = rowAtIndex(result, firstRow);
    result = setRow(result, rowAtIndex(result, secondRow), firstRow);
    result = setRow(result, tempRow, secondRow);
    return result;
}

Matrix* Matrix::switchColumns(Matrix *input, int firstColumn, int secondColumn)
{
    Matrix *result = copyMatrix(input);
    vector<double> *tempColumn = columnAtIndex(result, firstColumn);
    result = setColumn(result, columnAtIndex(result, secondColumn), firstColumn);
    result = setColumn(result, tempColumn, secondColumn);
    return result;
}

Matrix* Matrix::multiplyRow(Matrix *input, int row, double factor)
{
    Matrix *result = copyMatrix(input);
    vector<double> *tempRow = rowAtIndex(result, row);
    for (int i = 0; i < result->columns; i++) {
        tempRow->at(i) *= factor;
    }
    result = setRow(result, tempRow, row);
    return result;
}

Matrix* Matrix::multiplyColumn(Matrix *input, int column, double factor)
{
    Matrix *result = copyMatrix(input);
    vector<double> *tempColumn = columnAtIndex(result, column);
    for (int i = 0; i < result->rows; i++) {
        tempColumn->at(i) *= factor;
    }
    result = setColumn(result, tempColumn, column);
    return result;
}

Matrix* Matrix::addRowToRow(Matrix *input, int addRow, int modifiedRow, double factor)
{
    Matrix* result = copyMatrix(input);
    vector<double> *tempRow = rowAtIndex(result, addRow);
    vector<double> *modRow = rowAtIndex(result, modifiedRow);
    for (int i = 0; i < result->columns; i++) {
        tempRow->at(i) *= factor;
        modRow->at(i) += tempRow->at(i);
    }
    result = setRow(result, modRow, modifiedRow);
    return result;
}

Matrix* Matrix::addColumnToColumn(Matrix *input, int addColumn, int modifiedColumn, double factor)
{
    Matrix* result = copyMatrix(input);
    vector<double> *tempColumn = columnAtIndex(result, addColumn);
    vector<double> *modColumn = columnAtIndex(result, modifiedColumn);
    for (int i = 0; i < result->rows; i++) {
        tempColumn->at(i) *= factor;
        modColumn->at(i) += tempColumn->at(i);
    }
    result = setColumn(result, modColumn, modifiedColumn);
    return result;
}

Matrix* Matrix::zeilenStufenForm(Matrix *input)
{
    Matrix* result = copyMatrix(input);

    for (int x = 0; x < result->columns; x++) {
        vector<double> *currentColumn = columnAtIndex(result, x);
        int startIndex = -1;
        double startValue = 0.0;
        for (int y = x; y < result->rows; y++) {
            if (currentColumn->at(y) != 0.0) {
                startIndex = y;
                startValue = currentColumn->at(y);
                break;
            }
        }


        if (startIndex >= x) {
            if (startIndex > x) {
                result = switchRows(result, x, startIndex);
                currentColumn = columnAtIndex(result, x);
            }
            result = multiplyRow(result, x, 1/startValue);

            for (int y = x + 1; y < result->rows; y++) {
                if (currentColumn->at(y) != 0.0) {
                    result = addRowToRow(result, x, y, -currentColumn->at(y));
                }
            }
        }
    }

//    printMatrix(result);
    int currentY = 0;
    for (int x = 0; x < result->columns; x++) {
        if (x < result->rows) {
            vector<double> *currentColumn = columnAtIndex(result, x);
            if (currentColumn->at(currentY) == 1) {
                for (int y = currentY - 1; y >= 0; y--) {
                    result = addRowToRow(result, currentY, y, -currentColumn->at(y));
                }
                currentY++;
            }
        } else {
//            multiplyColumn(input, x, 0);
        }
    }

    return result;
}

Matrix* Matrix::triagonal(Matrix *input)
{
    Matrix *result = copyMatrix(input);
    if (result->rows == result->columns) {
        for (int x = 0; x < result->columns; x++) {
            if (result->value[x][x] == 0.0) {
                for (int i = x + 1; i < result->rows; i++) {
                    if (result->value[i][x] != 0.0) {
                        result = switchRows(result, x, i);
                        result = multiplyRow(result, x, -1.0);
                        break;
                    }
                }
            }

            for (int y = x + 1; y < result->rows; y++) {
                if (result->value[y][x] != 0.0) {
                    result = addRowToRow(result, x, y, -(result->value[y][x]/result->value[x][x]));
                }
            }
        }
        return result;
    } else {
        cout << "Matrix ist nicht quadratisch" << endl;
        return nullptr;
    }
}

double Matrix::determinante(Matrix *input)
{
    Matrix *temp = triagonal(input);
    double result = 1.0;
    for (int i = 0; i < temp->rows; i++) {
        result *= temp->value[i][i];
    }
    return result;
}

int Matrix::rang(Matrix *input)  //Funktioniert noch nicht zuverlässig.
{
    Matrix *temp = zeilenStufenForm(input);

    int currentY = 0;
    for (int x = 0; x < temp->columns; x++) {
        if (currentY < temp->rows) {
            vector<double> *currentColumn = columnAtIndex(temp, x);
            if (currentColumn->at(currentY) == 1) {
                currentY++;
            }
        }
    }

    return currentY;
}

Matrix* Matrix::multiply(Matrix *firstMatrix, Matrix *secondMatrix)
{
    if (firstMatrix->columns == secondMatrix->rows && firstMatrix->rows == secondMatrix->columns) {
        double **matrix = 0;
        matrix = new double *[firstMatrix->rows];
        for (int i = 0; i < firstMatrix->rows; i++) {
            matrix[i] = new double [secondMatrix->columns];
        }

        for (int y = 0; y < firstMatrix->rows; y++) {
            for (int x = 0; x < secondMatrix->columns; x++) {
                auto r = rowAtIndex(firstMatrix, y);
                auto c = columnAtIndex(secondMatrix, x);
                matrix[y][x] = dotProduct(r, c);
                delete r;
                delete c;
            }
        }
        Matrix *result = new Matrix();
        result->rows = firstMatrix->rows;
        result->columns = secondMatrix->columns;
        result->value = matrix;
        return result;
    } else {
        cout << "Multiplication not possible." << endl;
        return nullptr;
    }
}

void Matrix::printMatrix(Matrix *matrix)
{
    if (matrix != nullptr) {
        for (int y = 0; y < matrix->rows; y++) {
            for (int x = 0; x < matrix->columns; x++) {
                if (fabs(matrix->value[y][x]) == 0) {
                    matrix->value[y][x] = fabs(matrix->value[y][x]);
                }
                cout << matrix->value[y][x] << " ";
            }
            cout << endl;
        }
        cout << endl;
    } else {
        cout << "Matrix empty" << endl << endl;
    }
}

void Matrix::printMatrix()
{
    if (value != nullptr) {
        for (int y = 0; y < rows; y++) {
            for (int x = 0; x < columns; x++) {
                if (fabs(value[y][x]) == 0) {
                    value[y][x] = fabs(value[y][x]);
                }
                cout << value[y][x] << " ";
            }
            cout << endl;
        }
        cout << endl;
    } else {
        cout << "Matrix empty" << endl << endl;
    }
}

void Matrix::printVector(vector<double> *vec)
{
    if (vec != nullptr) {
        for (int i = 0; i<vec->size(); i++) {
            cout << vec->at(i) << "\n";
        }
    } else {
        cout << "Vector empty" << endl << endl;
    }
    cout << "\n";
}
