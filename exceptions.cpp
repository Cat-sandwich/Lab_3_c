#include "exceptions.h"
#include "matrix.h"
#include <iostream>
using namespace std;


Exception::Exception() {}

Invalid_Index::Invalid_Index() {}

void Invalid_Index::print() const {
	cout << "Произошла ошибка: индекс не найден" << endl;
}


Empty::Empty() {}

void Empty::print() const
{
	cout << "Произошла ошибка: в матрице нет элементов" << endl;
}

Different_Dimensions::Different_Dimensions() {}

void Different_Dimensions::print() const
{
	cout << "Произошла ошибка: размерности матриц не совпадают" << endl;
}

Divizion_By_Zero::Divizion_By_Zero() {}

void Divizion_By_Zero::print() const
{
	cout << "Произошла ошибка: деление матрицы на ноль" << endl;
}

Zero_Determinant::Zero_Determinant() {}

void Zero_Determinant::print() const
{
	cout << "Произошла ошибка: деление матрицы на нулевой определитель" << endl;
}

Dimensions_Incorrect::Dimensions_Incorrect() {}

void Dimensions_Incorrect::print() const
{
	cout << "Произошла ошибка: матрица должна быть квадратной" << endl;
}