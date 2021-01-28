#pragma once
#ifndef _VECTOR_H
#define _VECTOR_H
#include <cmath>
#include <iostream>
#include <new>
#include <vector>
#include <iterator>
#include <algorithm>
using namespace std;
const float eps = 0.000001;
template <class T>
class Vector
{
private:
	T* vector;
	const int size;
public:
	Vector(const int& n):size(n) {
		
		vector = new T[size];
		for (register int i = 0; i < this->size; i++)
			this->vector[i] = 0;
	};
	Vector(const Vector<T>& vec) :size(vec.size)
	{
		vector = new T[size];
		for (register int i = 0; i < size; i++)
			this->vector[i] = vec[i];
	};
	Vector(const T* const arr, const int arrSize):size(arrSize)
	{
		vector = new T[size];
		for (register int i = 0; i < size; i++)
			this->vector[i] = arr[i];
	};
	~Vector() {
		delete[] vector; vector = nullptr;
	};
	T& getPoint(const int& i) const {
		int j = i;
		while (j < 0 || j >= this->size)
		{
			cout << "Wrong index!Try new one" << endl;
			cin >> j;
		}
		return vector[i];
	};
	int getSize() const
	{
		return size;
	}
	Vector<T>& setPoint(const int& i, const T& val) {
		int j = i;
		while (j < 0 || j >= this->size)
		{
			cout << "Wrong index!Try new one" << endl;
			cin >> j;
		}
		vector[i] = val;
		return *this;
	};
	T operator[](const int& i) const
	{
		int j = i;
		while (j >= size && j < 0)
		{
			cout << "Wrong index!Try new one" << endl;
			cin >> j;
		}
		return vector[j];
	};
	T& operator[](const int& i)
	{
		int j = i;
		while (j >= size && j < 0)
		{
			cout << "Wrong index!Try new one" << endl;
			cin >> j;
		}
		return vector[j];
	};

	friend std::ostream& operator<<(std::ostream& out, const Vector<T>& vec)
	{    out << "(";
		for (register int i = 0; i < vec.size; i++)
		{
			if(i!=vec.size-1)
			{ if (fabs(vec[i]) >eps)
				out << vec[i] << ", ";
			else out << 0 << ", ";}
			if (i == vec.size - 1)
			{
				if (fabs(vec[i]) >eps)
					out << vec[i] << "";
				else out << 0 << "";
			}
			
		}
		out << ")";
		return out;
	};
	template <class T1>
	Vector<T>& operator +=(const Vector<T1>& v2)
	{
		if (this->size != v2.size)
		{
			cout << "These vectors have different sizes" << endl;
			return *this;
		}
		for (register int i = 0; i < this->size; i++)
			this->vector[i] += v2[i];
		return *this;
	};
	template <class T1>
	Vector<T>& operator -=(const Vector<T1>& v2)
	{
		if (this->size != v2.size)
		{
			cout << "These vectors have different sizes" << endl;
			return *this;
		}
		for (register int i = 0; i < this->size; i++)
			this->vector[i] -= v2[i];
		return *this;
	};
	template <class T1>
	Vector<T>& operator *=(const T1& val)
	{
		for (register int i = 0; i < this->size; i++)
		{
			this->vector[i] *= val;
		}
		return *this;
	};

	T findNormL1() const
	{
		T max = vector[0];
		for (register int i = 0; i < this->size; i++)
			if (vector[i] > max)
				max = vector[i];
		return max;
	};
	float findNormL2() const
	{
		float sum = 0;
		for (register int i = 0; i < this->size; i++)
			sum += pow(vector[i], 2);
		return sqrt(sum);
	};
	Vector<T>& findOrtogonal() const
	{
		if (this->size == 1)
		{
			cout << "Can not find ortogonal to the vector with size: 1" << endl;
		}
		bool f = true;
		int j = 0;
		for (register int i = 0; i < this->size; i++) {
			if (vector[i] != 0)
			{
				f = false;
				j = i;
			}
		}
		if (f)
		{
			for (register int i = 0; i < this->size; i++)
				vector[i] = 1;
			return *this;
		}

		Vector<T> helpVec(*this);
		for (register int i = 0; i < this->size; i++)
		{
			if (i == j)
			{
				vector[i] = 0;
			}
			else
				vector[i] = 1;
		}

		for (register int i = 0; i < this->size - 1; i++)
		{
			if (i != j)
				vector[this->size - 1] -= (helpVec[i]);
		}

		vector[j] /= helpVec[j];
		return *this;
	};
	template <class T1>
	bool checkCollinear(const Vector<T1>& vec) const
	{
		if (this->size != vec.size)
		{
			cout << "These vectors have different sizes" << endl;
			return false;
		}
		else {
			int j = 0;
			if ((this->vector[j] == 0 && vec[j] != 0) || (this->vector[j] != 0 && vec[j] == 0))
				return false;

			while (this->vector[j] == 0 && vec[j] == 0 && j < this->size)
			{
				j++;
				if ((this->vector[j] == 0 && vec[j] != 0) || (this->vector[j] != 0 && vec[j] == 0))
					return false;
			}
			if (j == size - 1)
			{
				cout << "There are two NULL vectors" << endl;
				return false;
			}
			float alpha = this->vector[j] / vec[j];
			for (register int i = j + 1; i < this->size; i++)
			{
				if ((this->vector[j] == 0 && vec[j] != 0) || (this->vector[j] != 0 && vec[j] == 0))
					return false;
				else if (this->vector[i] != 0 && vec[i] != 0)
				{
					if (this->vector[i] / vec[i] != alpha)
						return false;
				}
			}


			return true;
		}
	};
	template <class T1>
	float scalarMultiplication(const Vector<T1>& vec) const
	{
		if (this->size != vec.size)
		{
			cout << "These vectors have different sizes" << endl; return 0;
		}
		float rez = 0;

		for (register int i = 0; i < this->size; i++)
			rez += this->vector[i] * vec[i];

		return rez;
	};
	template <class T1>
	friend Vector<T> operator+(const Vector<T>& vec1, const Vector<T1>& vec2)
	{
		if (vec1.size != vec2.size)
		{
			cout << "These vectors have different sizes" << endl; return vec1;
		}
		Vector<T> helpVector(vec1.size);
		for (register int i = 0; i < vec1.size; i++)
			helpVector[i] = vec1[i] + vec2[i];
		return helpVector;
	};
	template <class T1>
	friend Vector<T> operator-(const Vector<T>& vec1, const Vector<T1>& vec2)
	{
		if (vec1.size != vec2.size)
		{
			cout << "These vectors have different sizes" << endl; return vec1;
		}
		Vector<T> helpVector(vec1.getSize());
		for (register int i = 0; i < vec1.size; i++)
			helpVector[i] = vec1[i] - vec2[i];
		return helpVector;
	};
	template <class T1>
	friend Vector<T> operator *(const T1& val, const Vector<T>& vec1)
	{
		Vector<T> helpVector(vec1.size);
		for (register int i = 0; i < vec1.size; i++)
		{
			helpVector[i] = val * vec1[i];
		}
		return helpVector;
	};
	template <class T1>
	Vector<T>& operator /=(const T1& val)
	{
		for (register int i = 0; i < this->size; i++)
		{
			this->vector[i] /= val;
		}
		return *this;
	};
	template <class T1>
	friend Vector<float> operator /(const Vector<T>& vec1, const T1& val)
	{
		Vector<float> helpVector(vec1.size);
		for (register int i = 0; i < vec1.size; i++)
		{
			helpVector[i] = vec1[i] / val;
		}
		return helpVector;
	};
	bool checkZeroVector() const
	{
		bool f = true;
		for (register int i = 0; i < this->size; i++)
		{
			if (fabs(vector[i]) > eps)
				f = false;
		}
		return f;
	};
	Vector<T>& makeZeroVec()
	{
		for (register int i = 0; i < this->size; i++)
			vector[i] = 0;
		return *this;
	};
	Vector<T>& operator=(const Vector<T>& vec)
	{
		if (this->size != vec.size)
		{
			cout << "Can not make equal vectors with different dimension!" << endl;
			return *this;
		}
		if (this == &vec)
			return *this;
		for (register int i = 0; i < this->size; i++)
			this->vector[i] = vec[i];
		return *this;
	};
	Vector<T>& operator=(Vector<T>&& vec) noexcept
	{
		if (this->size != vec.size)
		{
			cout << "Can not make equal vectors with different dimension!" << endl;
			return *this;
		}
		if (this == &vec)
			return *this;
		this->vector = vec.vector;
		vec.vector = nullptr;
		return *this;
	};


};


#endif