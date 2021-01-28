#pragma once
#ifndef _SIMPLEMATRIX_H
#define _SIMPLEMATRIX_H
#include "Vector.h"
#include <iomanip>
template <class T>
class SimpleMatrix {
protected:
	T** matrix;
	const int rows;
	const int cols;
	SimpleMatrix(const int& rows, const int& cols):rows(rows),cols(cols)
	{
	   matrix = new T * [rows];
		for (register int i = 0; i < rows; i++)
			matrix[i] = new T[cols];
		for (register int i = 0; i < this->rows; i++)
		{
			for (register int j = 0; j < this->cols; j++)
			{
				if (i == j)
					matrix[i][j] = 1;
				else
					matrix[i][j] = 0;
			}
		}
	};
	SimpleMatrix(const SimpleMatrix<T>& mat) :rows(mat.rows), cols(mat.cols)
	{
		matrix = new T * [rows];
		for (register int i = 0; i < rows; i++)
			matrix[i] = new T[cols];
		for (register int i = 0; i < this->rows; i++)
		{
			for (register int j = 0; j < this->cols; j++)
			{
				matrix[i][j] = mat(i, j);
			}
		}
	};
	SimpleMatrix(const T* const* const arr, const int& rowSize, const int& colsSize):rows(rowSize),cols(colsSize)
	{
		this->matrix = new T * [rows];
		for (register int i = 0; i < rows; i++)
			matrix[i] = new T[cols];
		for (register int i = 0; i < this->rows; i++)
		{
			for (register int j = 0; j < this->cols; j++)
			{
				this->matrix[i][j] = arr[i][j];
			}
		}
	};
	
	T getElement(const int& i, const int& j) const
	{
		int N = i;
		int M = j;
		while (N < 0 || N > this->rows)
		{
			cout << "Wrong row number!Try new one" << endl;
			cin >> N;
		}
		while (M < 0 || M > this->cols)
		{
			cout << "Wrong col number!Try new one" << endl;
			cin >> M;
		}
		return this->matrix[N][M];
	};
	SimpleMatrix<T> setElement(const int& i, const int& j, const T& elem)
	{

		int N = i;
		int M = j;
		while (N <0 || N > this->rows)
		{
			cout << "Wrong row number!Try new one" << endl;
			cin >> N;
		}
		while (M < 0 || M > this->cols)
		{
			cout << "Wrong col number!Try new one" << endl;
			cin >> M;
		}
		this->matrix[N][M] = elem;
		return *this;
	};
	const int getRowSize() const { return rows; };
	const int getColSize() const { return cols; };
	void getRow(const int& i, Vector<T>& arr) const
	{

		int N = i;
		while (N < 0 || N > this->rows)
		{
			cout << "Wrong row number!Try new one" << endl;
			cin >> N;
		}
		for (register int j = 0; j < this->cols; j++)
			arr[j] = this->matrix[N][j];
	};
	void getCol(const int& j, Vector<T>& arr) const
	{

		int M = j;
		while (M < 0 || M > this->cols)
		{
			cout << "Wrong col number!Try new one" << endl;
			cin >> M;
		}
		for (register int i = 0; i < this->rows; i++)
			arr[i] = this->matrix[i][M];
	};
	SimpleMatrix<T>& setRow(const int& i, const Vector<T>& arr)
	{

		int N = i;
		while (N < 0 || N > this->rows)
		{
			cout << "Wrong row number!Try new one" << endl;
			cin >> N;
		}
		for (register int j = 0; j < this->cols; j++)
			this->matrix[N][j] = arr[j];
		return *this;
	};
	SimpleMatrix<T>& setCol(const int& j, const Vector<T>& arr)
	{
		int M = j;
		while (M < 0 || M > this->cols)
		{
			cout << "Wrong col number!Try new one" << endl;
			cin >> M;
		}
		for (register int i = 0; i < this->cols; i++)
			this->matrix[i][M] = arr[i];
		return *this;
	};

	T operator()(const int& i, const int& j) const
	{
		int N = i;
		int M = j;
		while (N < 0 || N > this->rows)
		{
			cout << "Wrong row number!Try new one" << endl;
			cin >> N;
		}
		while (M < 0 || M > this->cols)
		{
			cout << "Wrong col number!Try new one" << endl;
			cin >> M;
		}
		return matrix[N][M];
	};
	T& operator()(const int& i, const int& j)
	{
		int N = i;
		int M = j;
		while (N < 0 || N > this->rows)
		{
			cout << "Wrong row number!Try new one" << endl;
			cin >> N;
		}
		while (M < 0 || M > this->cols)
		{
			cout << "Wrong col number!Try new one" << endl;
			cin >> M;
		}
		return matrix[N][M];
	};
	SimpleMatrix<T> transpMatrix() const
	{
		SimpleMatrix<T> helpMatrix(this->rows,this->cols);
		for (register int i = 0; i < this->cols; i++)
		{
			for (register int j = 0; j < this->rows; j++)
				helpMatrix(i, j) = this->matrix[j][i];
		}
		return helpMatrix;

	};
	T findNormL1() const
	{
		Vector<T> maxVector(this->cols);
		T max = 0;
		for (register int i = 0; i < this->rows; i++)
		{
			max = this->matrix[i][0];
			for (register int j = 0; j < this->cols; j++)
			{
				if (this->matrix[i][j] > max)
					max = this->matrix[i][j];
			}
			maxVector[i] = max;
		}
		max = maxVector[0];
		for (register int i = 0; i < maxVector.getSize(); i++)
		{
			if (maxVector[i] > max)
				max = maxVector[i];
		}
		return max;
	};
	float findNormL2() const
	{
		float sum = 0;
		for (register int i = 0; i < this->rows; i++)
		{
			for (register int j = 0; j < this->cols; j++)
				sum += pow(this->matrix[i][j], 2);
		}
		return sqrt(sum);
	};

	template<class T1>
	SimpleMatrix<T>& operator+=(const SimpleMatrix<T1>& mat)
	{
		if (this->cols != mat.cols || this->rows != mat.rows)
		{
			cout << "Can not add matrixes with different dimensions" << endl;
			return *this;
		}
		for (register int i = 0; i < this->rows; i++)
		{
			for (register int j = 0; j < this->cols; j++)
				this->matrix[i][j] += mat(i, j);
		}
		return *this;
	};
	template<class T1>
	SimpleMatrix<T>& operator -=(const SimpleMatrix<T1>& mat)
	{
		if (this->cols != mat.cols || this->rows != mat.rows)
		{
			cout << "Can not add matrixes wuth different dimension" << endl;
			return *this;
		}
		for (register int i = 0; i < this->rows; i++)
		{
			for (register int j = 0; j < this->cols; j++)
				this->matrix[i][j] -= mat(i, j);
		}
		return *this;
	};
	template<class T1>
	SimpleMatrix<T> operator *=(const SimpleMatrix<T1>& mat)
	{
		if (this->cols != mat.cols || this->rows != mat.rows)
		{
			cout << "Can not add matrixes wuth different dimension" << endl;
			return *this;
		}
		SimpleMatrix<T> helpMatrix(*this);
		for (register int i = 0; i < this->rows; i++)
		{
			for (register int j = 0; j < mat.cols; j++)
			{
				helpMatrix(i, j) = 0;
				for (register int k = 0; k < helpMatrix.cols; k++)
					helpMatrix(i, j) += this->matrix[i][k] * mat(k, j);
			}
			*this = helpMatrix;
		}
		return *this;
	};

	template<class T1>
	SimpleMatrix<T> operator /=(const T1& val)
	{
		for (register int i = 0; i < this->rows; i++)
		{
			for (register int j = 0; j < this->cols; j++)
				this->matrix[i][j] /= val;
		}
		return *this;
	};
	SimpleMatrix<T> Pow(const int& i)
	{
		int j = i;
		SimpleMatrix<T> helpMatrix(*this);
		while (j > 1)
		{
			(*this) *= helpMatrix;
			j--;
		}
		return *this;
	};
	template<class T1>
	Vector<T> vectorMultiplication(const Vector<T1>& vec)
	{
		Vector<T> helpVector(this->cols);
		for (register int i = 0; i < helpVector.getSize(); i++)
		{
			for (register int j = 0; j < this->cols; j++)
			{
				helpVector[i] += this->matrix[i][j] * vec[j];
			}
		}
		return helpVector;
	};
	SimpleMatrix<T>& operator=(const SimpleMatrix<T>& matr)
	{
		if (this->cols != matr.cols || this->rows != matr.rows)
		{
			cout << "Can not make equal matrixes wuth different dimension" << endl;
			return *this;
		}
		if (this == &matr)
			return *this;
		for (register int i = 0; i < this->rows; i++)
		{
			for (register int j = 0; j < this->cols; j++)
				this->matrix[i][j] = matr(i, j);
		}
		return*this;
	};
	SimpleMatrix<T>& operator=(SimpleMatrix<T>&& matr)
	{
		if (this->cols != matr.cols || this->rows != matr.rows)
		{
			cout << "Can not add matrixes wuth different dimension" << endl;
			return *this;
		}
		if (this == &matr)
			return *this;
	    this->matrix = matr.matrix;
		matr.matrix = nullptr;
		return *this;
	};
	void makeZeroMatrix()
	{
		for (register int i = 0; i < this->rows; i++)
		{
			for (register int j = 0; j < this->cols; j++)
				this->matrix[i][j] = 0;
		}
	};
public:

	template<class T1>
	friend SimpleMatrix<T> operator *(const SimpleMatrix<T>& mat1, const SimpleMatrix<T1>& mat)
	{
		
		SimpleMatrix<T> helpMatrix(mat1.rows,mat.cols);
		if (this->cols != mat.rows)
		{
			cout << "Can not multiplicate these matrixes!" << endl;
			return helpMatrix;
		}
		for (register int i = 0; i < mat1.rows; i++)
		{
			for (register int j = 0; j < mat.cols; j++)
			{
				helpMatrix(i, j) = 0;
				for (register int k = 0; k < mat1.cols; k++)
					helpMatrix(i, j) += mat1(i, k) * mat(k, j);
			}
		}
		return helpMatrix;
	};
	
	template<class T1>
	friend SimpleMatrix<T> operator*(const T& val, const SimpleMatrix<T1>& mat)
	{
		SimpleMatrix<T> helpMatrix(mat.rows,mat.cols);
		for (register int i = 0; i < mat.rows; i++)
		{
			for (register int j = 0; j < mat.cols; j++)
				helpMatrix(i, j) = val * mat(i, j);
		}
		return helpMatrix;
	};
	template <class T1>
	friend SimpleMatrix<T> operator + (const SimpleMatrix<T>& mat1, const SimpleMatrix<T1>& mat)
	{
		if (mat1.rows != mat.rows || mat1.cols != mat.cols)
		{
			cout << "Can not add matrixes with different dimensions" << endl;
			return mat1;
		}
		SimpleMatrix<T> helpMatrix(mat1.rows,mat1.cols);
		for (register int i = 0; i < mat1.rows; i++)
		{
			for (register int j = 0; j < mat1.cols; j++)
				helpMatrix(i, j) = mat1(i, j) + mat(i, j);
		}
		return helpMatrix;
	};
	template <class T1>
	friend SimpleMatrix<T> operator - (const SimpleMatrix<T>& mat1, const SimpleMatrix<T1>& mat)
	{
		if (mat1.rows != mat.rows || mat1.cols != mat.cols)
		{
			cout << "Can not add matrixes with different dimensions" << endl;
			return mat1;
		}
		SimpleMatrix<T> helpMatrix(mat1.rows, mat1.cols);
		for (register int i = 0; i < mat1.rows; i++)
		{
			for (register int j = 0; j < mat1.cols; j++)
				helpMatrix(i, j) = mat1(i, j) - mat(i, j);
		}
		return helpMatrix;
	};
	template <class T1>
	friend SimpleMatrix<T> operator / (const SimpleMatrix<T>& mat1, const T1& val)
	{
		SimpleMatrix<T> helpMatrix(mat1.rows, mat1.cols);
		for (register int i = 0; i < mat1.rows; i++)
		{
			for (register int j = 0; j < mat1.cols; j++)
				helpMatrix(i, j) = mat1(i, j) / val;
		}
		return helpMatrix;
	};

	friend ostream& operator<<(ostream& out, const SimpleMatrix<T>& mat)
	{
		for (register int i = 0; i < mat.rows; i++)
		{
			for (register int j = 0; j < mat.cols; j++)
			{
				if (fabs(mat(i, j)) > 0.0001)
					out << setw(11) << mat.getElement(i, j);
				else out << setw(11) << 0;
			}
			out << endl;
		}
		return out;
	};
	virtual ~SimpleMatrix() {
		for (register int i = 0; i < this->rows; i++)
			delete[] matrix[i];
		delete[] matrix;
		matrix = nullptr;
	};
};
#endif
