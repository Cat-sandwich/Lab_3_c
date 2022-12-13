#include <iostream>
#include <math.h>
#include <complex>
#include <iomanip>
#include <vector>
using namespace std;

template <class T>
class Matrix
{
private:
	vector<vector<T>> data;
	int m, n;

	auto begin() { return data.begin(); }
	auto end() { return data.end(); }

public:
	Matrix() = default;

	Matrix(int m, int n);

	Matrix(int m, int n, const T& value);

	Matrix(const Matrix& Matrix) = default;

	void Set_m(int m = 1);

	void Set_n(int n = 1);

	int Get_m();

	int Get_n();

	void Set_Data(const T& value); 

	T Get_Data(int i, int j) const;

	auto cbegin() const { return data.cbegin(); }

	auto cend() const { return data.cend(); }

	~Matrix() = default;

	T operator () (int m, int n) const;

	Matrix& operator = (const Matrix& M) = default;

	Matrix& operator () (int m, int n, const T& value);

	Matrix operator + (const Matrix& New_Matrix);

	Matrix operator - (const Matrix& New_Matrix);

	Matrix operator * (const Matrix& New_Matrix);

	Matrix operator * (const T& scalar);


	friend Matrix operator * (const T& scalar, Matrix& Matrix)
	{
		return Matrix * scalar;
	}

	Matrix operator / (const T& scalar);

	T Ñalculating_trace_matrix();

	Matrix  Transpose();

	void Random();

	Matrix Pre_Minor(int row, int col) const;

	T NDeterminant(int size);

	Matrix Search_Matrix_X(const Matrix& Vector);

	friend ostream& operator << (ostream& os, const Matrix& New_Matrix)
	{
		for (auto it = New_Matrix.cbegin(); it != New_Matrix.cend(); it++) 
		{
			for (auto iter : (*it)) 
			{
				os << setw(10);
				os<<  iter;
			}
			cout << endl;
		}
		return os;
	}

};

