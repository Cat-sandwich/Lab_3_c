#pragma once
#include "matrix.h"
#include "exceptions.h"
#include <iostream>
#include <random>
#include <complex>
using namespace std;

template <class T>
Matrix<T>::Matrix()
{
	m = 0;
	n = 0;

}

template <class T>
Matrix<T>::Matrix(int m, int n)
{
	this->m = m;
	this->n = n;

}

template <class T>
Matrix<T>::Matrix(int m, int n, const T& value)
{

	Set_m(m);
	Set_n(n);

	vector<vector<T>> new_data(m, vector<T>(n, value));

	*this = new_data;
}

template <class T>
Matrix<T>::Matrix(const Matrix<T>& Matrix)
{
	m = Matrix.m;
	n = Matrix.n;


	auto it = Matrix.cbegin();

	auto jt = Matrix.cbegin(vector<T>);

	for (it; it != Matrix.cend(); it++)
		for (jt; jt != Matrix.cend(vector<T>); jt++)
			data.push_back(*it(vector<T>(*jt)));


}

template <class T>
void Matrix<T>::Set_m(int m)
{
	if (m > 0)
		this->m = m;
}

template <class T>
void Matrix<T>::Set_n(int n)
{
	if (n > 0)
		this->n = n;
}

template <class T>
int Matrix<T>::Get_m()
{
	return m;
}

template <class T>
int Matrix<T>::Get_n()
{
	return n;
}

template <class T>
void Matrix<T>::Set_Data(const T& value)
{
	vector<vector<T>> new_data(m, vector<T>(n, value));

	*this = new_data;
}

template <class T>
T Matrix<T>::Get_Data(int i, int j) const //todo+
{
	if ((m <= i) && (n <= j)) throw Invalid_Index();
	auto it = data.cbegin();
	auto jt = data.cbegin(<vector<T>>);
	it += i;
	jt += j;
	return *it(*jt);

}




template <class T>
Matrix<T>& Matrix<T>::operator = (const Matrix<T>& M) //ToDO+
{
	if ((m == M.m) && (n == M.n))
	{
		data = M;
		return *this;
	}
	if ((m < M.m) || (n < M.n))
	{
		data.clear();

	}
	m = M.m;
	n = M.n;

	auto it = Matrix.cbegin();

	auto jt = Matrix.begin(vector<T>);


	for (it; it != Matrix.cend(); it++)
		for (jt; jt != Matrix.cend(vector<T>); jt++)
			data.push_back(*it(vector<T>(*jt)));
	return *this;
}

template <class T>
Matrix<T>::~Matrix()
{
	data.clear();
}

template <class T>
T& Matrix<T>::operator ()(int m, int n) const
{

	return data.at(m, (vector<T>(n)));
}

template <class T>
Matrix<T>& Matrix<T>::operator () (int m, int n, const T& value)
{
	if ((m > this->m) || (n > this->n)) throw Invalid_Index();


	data.at(m, (vector<T> n)) = value;
	return *this;

}

template <class T>
Matrix<T> Matrix<T>::operator + (const Matrix<T>& New_Matrix) {
	if (this->m != New_Matrix.m || this->n != New_Matrix.n) throw Different_Dimensions();
	Matrix<T> res = data;
	auto it = res.cbegin();
	auto jt = res.cbegin(<vector<T>>);

	auto iter = New_Matrix.cbegin();
	auto jter = (*iter).cbegin();
	for (it; it != res.cend(); it++) {
		for (jt; jt != res.cend(vector<T>); jt++)
		{
			(*jt) = (*jt) + (*jter);

			jter++;
		}
		iter++;
	}
	return res;
}

template <class T>
Matrix<T> Matrix<T>::operator - (const Matrix<T>& New_Matrix) {
	if (this->m != New_Matrix.m || this->n != New_Matrix.n) throw Different_Dimensions();

	Matrix<T> res = data;
	auto it = res.cbegin();
	auto jt = res.cbegin(<vector<T>>);

	auto iter = New_Matrix.cbegin();
	auto jter = (*iter).cbegin();
	for (it; it != res.cend(); it++) {
		for (jt; jt != res.cend(vector<T>); jt++)
		{
			(*jt) = (*jt) - (*jter);

			jter++;
		}
		iter++;
	}
	return res;
}

template <class T>
Matrix<T> Matrix<T>::operator * (const Matrix<T>& New_Matrix)
{
	if (n != New_Matrix.m) throw Different_Dimensions();

	Matrix<T> res(m, New_Matrix.n, T(0));

	for (auto it = res.cbegin(); it != res.cend(); ++it)
		for (auto jt = (*it).cbegin(); jt != (*jt).cend(); ++jt)
		{
			for (int k = 0; k < m; ++k)
			{
				(*jt) += (*it)[k] ;
			}
		}


	for (int i = 0; i < res.m; ++i)
		for (int j = 0; j < res.n; ++j)
		{

			for (int k = 0; k < m; ++k)

				res.data[i][j] += data[i][k] * New_Matrix.data[k][j];
		}
	return res;
}

template <class T>
Matrix<T> Matrix<T>::operator * (const T& scalar)
{
	Matrix<T> res = data;

	for (auto it = res.cbegin(); it != res.cend(); it++)
		for (auto jt = (*it).cbegin(); jt != (*jt).cend(); jt++)
			(*jt) = (*jt) * scalar;


	return res;
}

template <class T>
Matrix<T> Matrix<T>::operator / (const T& scalar)
{
	if (scalar == T(0)) throw Divizion_By_Zero();
	Matrix<T> res = data;

	for (auto it = res.cbegin(); it != res.cend(); it++)
		for (auto jt = (*it).cbegin(); jt != (*jt).cend(); jt++)
			(*jt) = (*jt) / scalar;

	return res;
}

template <class T>
T Matrix<T>::Сalculating_trace_matrix()
{
	if (n != m) throw Dimensions_Incorrect();
	T trace = 0;
	for (auto it = data.cbegin(); it != data.cend(); it++)
		for (auto jt = (*it).cbegin(); jt != (*jt).cend(); jt++)
		{
			if (it == jt)
				trace += (*jt);
		}

	return trace;
}

template <class T>
Matrix<T> Matrix<T>::Transpose()
{
	Matrix<T> Transposed(m, n);
	for (int i = 0; i < m; ++i)
	{
		for (int j = i; j < n; ++j)
		{
			Transposed.data[i][j] = data[j][i];
			Transposed.data[j][i] = data[i][j];
		}
		cout << '\n';
	}

	return Transposed;
}

template <class T>
void Matrix<T>::Random()
{
	srand(time(0));
	for (auto it = data.cbegin(); it != data.cend(); it++)
		for (auto jt = (*it).cbegin(); jt != (*jt).cend(); jt++)
			(*jt) = T((1 + rand() % 100) / 10.0);

}

template <class T>
Matrix<T> Matrix<T>::Pre_Minor(int row, int col) const
{
	if (n != m) throw Different_Dimensions();
	Matrix<T> New_Matrix(m - 1, n - 1);
	int in = 0, jn = 0;

	for (int i = 0; i < m; i++)
	{
		if (i != row)
		{
			jn = 0;
			for (int j = 0; j < m; j++) {
				if (j != col) {
					New_Matrix.data[in][jn] = data[i][j];
					jn++;
				}
			}
			in++;
		}
	}

	return New_Matrix;
}

template <class T>
T Matrix<T>::NDeterminant(int size)
{
	if (n != m) throw Different_Dimensions();
	Matrix<T> TmpMatrix(m, m);
	int new_size = size - 1;
	T d = 0;
	int k = 1; //(-1) в степени i
	if (size < 1) cout << "Определитель вычислить невозможно!";
	if (m == 1)
	{
		d = data[0][0];
		return(d);
	}
	if (m == 2)
	{
		d = data[0][0] * data[1][1] - data[1][0] * data[0][1];
		return(d);
	}
	if (m > 2)
	{


		for (int i = 0; i < size; i++)
		{

			TmpMatrix = (*this).Pre_Minor(i, 0);

			d += T(k) * data[i][0] * TmpMatrix.NDeterminant(new_size);
			k = -k;
		}

	}

	return(d);
}


template <class T>
Matrix<T> Matrix<T>::Search_Matrix_X(const Matrix<T>& Vector)
{
	if (this->n != Vector.m) throw Different_Dimensions();

	T det = (*this).NDeterminant(m);
	if (det == T(0)) throw Zero_Determinant();

	cout << det << endl;
	Matrix<T> Minors(n, m);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			Minors.data[i][j] = (T(pow(-1, (i + j)))) * Pre_Minor(i, j).NDeterminant(m - 1);

		}
	}
	cout << Minors << endl;

	Matrix<T> Minors_Transpose = Minors.Transpose();
	cout << Minors_Transpose << endl;
	Matrix<T> Ans = (T(1) / det) * Minors_Transpose * Vector;


	return Ans;
}

//template <>
//void Matrix<complex<double>>::Random()
//{
//	srand(time(0));
//	for (int i = 0; i < m; ++i)
//	{
//		for (int j = 0; j < n; ++j)
//		{
//			complex<double> value(((1 + rand() % 100) / 10.0), ((1 + rand() % 100) / 10.0));
//			data[i][j] = value;
//		}
//	}
//}

//template <>
//void Matrix<complex<float>>::Random()
//{
//
//	srand(time(0));
//
//	for (auto it = data.cbegin(); it != data.cend(); it++)
//		for (auto jt = (*it).cbegin(); jt != (*jt).cend(); jt++)
//		{
//			complex<float> value(((1 + rand() % 100) / 10.0), ((1 + rand() % 100) / 10.0));
//			(*jt) = value;
//		}
//
//}

template class Matrix<int>;
template class Matrix<float>;
template class Matrix<double>;
template class Matrix<complex<float>>;
template class Matrix<complex<double>>;