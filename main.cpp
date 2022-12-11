#include <stdio.h>
#include "matrix.h"
#include "exceptions.h"
#include <iostream>
#include <stdio.h>
#include <conio.h>
#include <Windows.h>
#include <complex>




int Check_Int()
{
	int number = 0;
	while (number <= 0)
	{
		while (!(cin >> number) || (cin.peek() != '\n'))
		{
			cin.clear();
			while (cin.get() != '\n');
			cout << "������� ���������� ��������...\n";
		}
		if (number <= 0) cout << "������� ���������� ��������...\n";

	}

	return number;
}
template<class T>
T Check_Double()
{
	T number;

	while (!(cin >> number) || (cin.peek() != '\n'))
	{
		cin.clear();
		while (cin.get() != '\n');
		cout << "������� ���������� ��������...\n";
	}


	return number;
}

template <class T>
void �hange_Matrix(Matrix<T>** Many_Matrix, const int& current)
{
	cout << "������� ������� ��������, ������� ����� ��������:" << endl;
	cout << "������� ����� ������: ";
	int i = Check_Int() - 1;
	cout << "������� ����� �������: ";
	int j = Check_Int() - 1;

	try
	{
		((*Many_Matrix)[current]).Get_Data(i, j);
		cout << "������� ��������, �� ������� �� ������ �������� �����:";
		T value = Check_Double<T>();
		((*Many_Matrix)[current])(i, j, value);
	}
	catch (Exception& error)
	{
		error.print();

	}
}

template <class T>
Matrix<T> Random_Matrix(int m, int n)
{
	Matrix<T> New_matrix(m, n);
	New_matrix.Random();
	return  New_matrix;
}

template <class T>
void Add_Matrix(int* size, Matrix<T>** Many_Matrix, Matrix<T> New_matrix)
{
	*size += 1;
	Matrix<T>* tmp = new Matrix<T>[*size];
	if (*size - 1 != 0)
	{
		for (int i = 0; i < (*size) - 1; i++)
			tmp[i] = (*Many_Matrix)[i];
		delete[](*Many_Matrix);
	}
	tmp[*size - 1] = New_matrix;
	*Many_Matrix = tmp;

}

template <class T>
void Print_Matrix(Matrix<T>* Many_Matrix, int current, int size)
{

	if (Many_Matrix == NULL) cout << "������ ���(\n\n";
	else cout << "������� �" << (current + 1) << "\n" << (Many_Matrix)[current];
}

int get_key()
{
	int key = _getch();
	if ((key == 0) || (key == 224)) key = _getch();
	return key;
}

template <class T>
void menu()
{
	int key = 0;
	bool menu1 = true;
	Matrix<T>* Many_Matrix = NULL;
	Matrix<T> New_matrix;
	Matrix<T>* Vector = NULL;
	int current = 0, size = 0;
	T value = 0;
	while (menu1)
	{

		system("cls");
		cout << "\t��� �������\n" << endl;
		Print_Matrix<T>(Many_Matrix, current, size);

		cout << "1 - ������ ������� ��������\n2 - ������� ��� ������� \n3 - ������� �� ����� ������� ������" << endl;
		cout << "4 - �������� ���� ������� �� ������\n5 - ��������� �� ������\n6 - ��������� �� ������" << endl;
		cout << "7 - ��������� ���� �������\n8 - ��������� �������" << endl;
		cout << "9 - �������� ���� �������� � ������� �������\n0 - ��������� ������" << endl;
		cout << "\"+\" - ������ ������� � ������������ ������" << endl;
		cout << "-> ������\n-< �����\n" << endl;

		T Scalar = 0;

		key = get_key();
		int m = 0, n = 0;
		switch (key)
		{
		case(49):

			cout << "\t������� ����������� �������: \n";
			cout << "\n������� ���������� ��������: ";
			m = Check_Int();

			cout << "������� ���������� �����: ";
			n = Check_Int();
			printf("������ �������\n ");
			system("pause");

			New_matrix = Random_Matrix<T>(m, n);
			Add_Matrix<T>(&size, &Many_Matrix, New_matrix);
			current = size - 1;

			break;
		case(50):
			if ((size + 1) < 2)
			{
				cout << "\n� ��� ������ ���� ���� �� ��� �������...\n" << endl;
				system("pause");
				break;
			}
			cout << "������� ����� �������, ������� �� ������ ������� � �������:\n" << endl;
			m = Check_Int();
			while (m > size)
			{
				cout << "����� ������� ���..." << endl;
				cout << "���������� ��� ���..." << endl;
				m = Check_Int();

			}
			try
			{
				New_matrix = Many_Matrix[current] + Many_Matrix[m - 1];
				Add_Matrix<T>(&size, &Many_Matrix, New_matrix);
				current = size - 1;
			}
			catch (Exception& error)
			{
				error.print();
				system("pause");
			}

			break;
		case(51):
			if ((size + 1) < 2)
			{
				cout << "\n� ��� ������ ���� ���� �� ��� �������...\n" << endl;
				system("pause");
				break;
			}
			cout << "������� ����� �������, ������� �� ������ �������:\n" << endl;
			m = Check_Int();

			while (m > size)
			{
				cout << "����� ������� ���..." << endl;
				cout << "���������� ��� ���..." << endl;
				m = Check_Int();
			}
			try
			{
				New_matrix = Many_Matrix[current] - Many_Matrix[m - 1];
				Add_Matrix<T>(&size, &Many_Matrix, New_matrix);
				current = size - 1;
			}
			catch (Exception& error)
			{
				error.print();
				system("pause");
			}

			break;
		case(52):
			if ((size + 1) < 2)
			{
				cout << "\n� ��� ������ ���� ���� �� ��� �������...\n" << endl;
				system("pause");
				break;
			}
			do
			{
				cout << "������� ����� �������, ������� �� ������ �������� �� �������:\n" << endl;
				m = Check_Int();
				system("pause");
			} while (m > size);
			try
			{
				New_matrix = Many_Matrix[current] * Many_Matrix[m - 1];
				Add_Matrix<T>(&size, &Many_Matrix, New_matrix);
				current = size - 1;
			}
			catch (Exception& error)
			{
				error.print();
				system("pause");
			}

			break;

		case(53):
			if (size == 0)
			{
				cout << "\n� ��� ������ ���� ���� �� ���� �������...\n" << endl;
				system("pause");
				break;
			}
			cout << "������� �����: ";
			Scalar = Check_Double<T>();
			New_matrix = Many_Matrix[current] * Scalar;
			Add_Matrix<T>(&size, &Many_Matrix, New_matrix);
			system("pause");
			current = size - 1;
			break;
		case(54):
			if (size == 0)
			{
				cout << "\n� ��� ������ ���� ���� �� ���� �������...\n" << endl;
				system("pause");
				break;
			}
			cout << "������� �����:" << endl;
			Scalar = Check_Double<T>();
			try
			{
				New_matrix = Many_Matrix[current] / Scalar;

				Add_Matrix<T>(&size, &Many_Matrix, New_matrix);
				current = size - 1;
			}
			catch (Exception& Error)
			{
				Error.print();
				system("pause");
			}


			break;

		case(55):
			if (size == 0)
			{
				cout << "\n� ��� ������ ���� ���� �� ���� �������...\n" << endl;
				system("pause");
				break;
			}
			try
			{
				cout << "���� ������� = " << Many_Matrix[current].�alculating_trace_matrix() << endl;
				system("pause");
			}
			catch (Exception& Error)
			{
				Error.print();
				system("pause");
			}
			break;

		case(56):

			if (size == 0)
			{
				cout << "\n� ��� ������ ���� ���� �� ���� �������...\n" << endl;
				system("pause");
				break;
			}
			if ((Many_Matrix[current].Get_m() == Many_Matrix[current].Get_n()) && (Many_Matrix[current].Get_m() != 3))
			{
				cout << "\n�������� ������� � �������� ��� �� ���..." << endl;
				system("pause");
				break;
			}

			cout << "������� ������ ������������ 3 �� 1...\n" << endl;
			Vector = new Matrix<T>(3, 1);
			n = 0;
			for (int m = 0; m < 3; m++)
			{
				cout << "������� B[" << m + 1 << "][" << n + 1 << "]: ";
				value = Check_Double<T>();
				(*Vector)(m, n, value);


			}
			cout << Vector;


			Add_Matrix<T>(&size, &Many_Matrix, *Vector);
			try
			{
				New_matrix = Many_Matrix[current].Search_Matrix_X(*Vector);
				Add_Matrix<T>(&size, &Many_Matrix, New_matrix);
			}
			catch (Exception& Error)
			{
				Error.print();
				system("pause");
			}

			system("pause");

			break;
		case(57):
			if (size == 0)
			{
				cout << "\n� ��� ������ ���� ���� �� ���� �������...\n" << endl;
				system("pause");
				break;
			}

			�hange_Matrix<T>(&Many_Matrix, current);
			system("pause");

			break;
		case(61):
		{
			cout << "\t������� ����������� �������: \n";
			cout << "\n������� ���������� ��������: ";
			m = Check_Int();

			cout << "������� ���������� �����: ";
			n = Check_Int();

			cout << "������� �����, ������� ����� ��������� ��� �������:";
			value = Check_Double<T>();
			New_matrix = Matrix<T>(m, n, value);
			Add_Matrix<T>(&size, &Many_Matrix, New_matrix);
			current = size - 1;
			cout << "������ �������\n ";
			system("pause");
		}
		break;
		case(48):
			system("cls");


			menu1 = false;
			break;
		case 75:
			if (current > 0) current--;
			break;
		case 77:
			if (current < size - 1) current++;
			break;
		default:
			break;
		}
	}
}

int main()
{
	setlocale(LC_ALL, "RUS");
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	bool exit = true;
	while (exit) {
		cout << "������������! ��� ������������ ��������� \"����� ������\"\n" << endl;
		cout << "�������:" << endl;
		cout << "1 - ����� �������� � ����� int" << endl;
		cout << "2 - ����� �������� � ����� double" << endl;
		cout << "3 - ����� �������� � ����� float" << endl;
		cout << "4 - ����� �������� � ����� complex-double" << endl;
		cout << "5 - ����� �������� � ����� complex-float" << endl;
		cout << "0 - ����� ��������� ������" << endl;

		int key = get_key();

		switch (key)
		{
		case 49:
			menu<int>();
			break;
		case 50:
			menu<double>();
			break;
		case 51:
			menu<float>();
			break;
		case 52:
			menu<complex<double>>();
			break;
		case 53:
			menu<complex<float>>();
			break;
		case 48:
			exit = false;
			break;
		default:
			break;
		}

	}

	return 0;
}