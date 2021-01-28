#pragma once
#ifndef _SQUAREMATRIX_H
#define _SQUAREMATRIX_H
#include "SimpleMatrix.h"
#include "IMatrix.h"
template<class T>
class SquareMatrix :public SimpleMatrix<T>, public IMatrix {
private:
	int findNoZeroElement(const int& colIndex) const
	{
		for (register int i = colIndex + 1; i < this->rows; i++)
		{
			if (fabs(this->matrix[i][colIndex]) > 0.0001)
				return i;
		}
		return 0;
	};
	void swap(const int& firstRowIndex, const int& secondRowIndex)
	{
		Vector<T>  vec1(this->cols);
		Vector<T> vec2(this->cols);
		this->getRow(firstRowIndex, vec1);
		this->getRow(secondRowIndex, vec2);
		this->setRow(firstRowIndex, vec2);
		this->setRow(secondRowIndex, vec1);
		~flag;
	};
	int findAmountOfZeroRows() const
	{
		int amount = 0;
		bool f = true;
		for (register int i = 0; i < this->rows; i++)
		{
			f = true;
			for (register int j = 0; j < this->cols; j++)
			{
				if (fabs(this->matrix[i][j]) > 0.0001)
					f = false;

			}
			if (f) amount++;
		}
		return amount;
	};
	template<class T1>
	void findTriangleForm(Vector<T1>& freeVec)
	{

		for (register int j = 0; j < this->cols - 1; j++)
		{
			Vector<T> vec1(this->cols);
			Vector<T> vec2(this->cols);
			bool f = true;
			if (fabs(this->matrix[j][j]) < 0.0001)
			{
				int i = this->findNoZeroElement(j);
				if (i)
				{
					this->swap(j, i);
					T1 temp = freeVec[j];
					freeVec[j] = freeVec[i];
					freeVec[i] = temp;
				}
				else f = false;
			}
			if (f)
			{
				for (register int i = j + 1; i < this->rows; i++)
				{
					float mCoef = this->matrix[i][j] / this->matrix[j][j];
					Vector<T> vec1(this->cols);
					Vector<T> vec2(this->cols);
					this->getRow(j, vec1);
					this->getRow(i, vec2);
					vec2 = vec2 - mCoef * vec1;
					this->setRow(i, vec2);
					freeVec[i] = freeVec[i] - mCoef * freeVec[j];

				}

			}
		}

	};
	bool checkDiagonal() const
	{
		for (register int i = 0; i < this->rows; i++)
		{
			for (register int j = 0; j < i; j++)
			{
				if (fabs(this->matrix[i][j]) >= eps)
					return false;
			}
		}
		return true;
	};
	bool QR_methodOwnVal(Vector<float>& ownValVec) const
	{
		SquareMatrix<T> Q_matr(this->cols);
		SquareMatrix<T> R_matr(this->cols);
		SquareMatrix<T> helpMatrix(*this);
		if (!helpMatrix.QR_decomposition(Q_matr, R_matr)) return false;
		R_matr *= Q_matr;
		helpMatrix = R_matr;
		long int checkIterations = 0;
		while (!helpMatrix.checkDiagonal() && checkIterations < 100000)
		{
			if (!helpMatrix.QR_decomposition(Q_matr, R_matr)) return false;
			R_matr *= Q_matr;
			helpMatrix = R_matr;
			checkIterations++;
		}
		if (checkIterations == 1000000)
		{
			cout << "Can not find own numbers with this method!" << endl;
			return false;
		}
		for (register int i = 0; i < helpMatrix.cols; i++)
		{
			ownValVec[i] = helpMatrix(i, i);
		}
		return true;

	};
	bool LU_methodOwnVal(Vector<float>& ownValVec) const
	{
		SquareMatrix<T> L_matr(this->cols);
		SquareMatrix<T> R_matr(this->cols);
		SquareMatrix<T> helpMatrix(*this);
		if (!helpMatrix.LU_decomposition(L_matr, R_matr)) return false;
		R_matr *= L_matr;
		helpMatrix = R_matr;
		long int checkIterations = 0;
		while (!helpMatrix.checkDiagonal() && checkIterations < 100000)
		{
			if (helpMatrix.LU_decomposition(L_matr, R_matr)) {
				R_matr *= L_matr;
				helpMatrix = R_matr;
				checkIterations++;
			}
			else return false;

		}
		if (checkIterations == 10000000)
		{
			cout << "Can not find own nu,bers with this method!" << endl;
			return false;
		}
		for (register int i = 0; i < helpMatrix.cols; i++)
		{
			ownValVec[i] = helpMatrix(i, i);
		}
		return true;
	};
	bool fullSpectralCaseQR(Vector< float>& ownValues, SquareMatrix< float>& ownVectors) const
	{
		if (!this->QR_methodOwnVal(ownValues)) return false;
		Vector<int> freeVec(this->cols);
		SquareMatrix<T> helpMatrix(*this);
		Vector<float>* solutionVec = new Vector<float>(this->cols);
		int ownVectorMatrixIndex = 0;
		for (register int i = 0; i < ownValues.getSize(); i++)
		{
			helpMatrix = *this;
			for (register int k = 0; k < helpMatrix.cols; k++)
			{
				helpMatrix(k, k) -= ownValues[i];
				freeVec[k] = 0;
			}
			helpMatrix.findSlauSolution(freeVec, solutionVec, true);
			if (!solutionVec->checkZeroVector())
				{
					ownVectors.setCol(ownVectorMatrixIndex, *solutionVec);
					ownVectorMatrixIndex++;

				}
			solutionVec->makeZeroVec();
		}
		
		 delete solutionVec;
		 solutionVec = nullptr;
		return true; 
	};
	bool fullSpectralCaseLU(Vector< float>& ownValues, SquareMatrix< float>& ownVectors) const
	{
		if (!this->LU_methodOwnVal(ownValues)) return false;
		Vector<int> freeVec(this->cols);
		SquareMatrix<T> helpMatrix(*this);
		Vector<float>* solutionVec = new Vector<float>(this->cols);
		int ownVectorMatrixIndex = 0;
		for (register int i = 0; i < ownValues.getSize(); i++)
		{
			helpMatrix = *this;
			for (register int k = 0; k < helpMatrix.cols; k++)
			{
				helpMatrix(k, k) -= ownValues[i];
				freeVec[k] = 0;
			}
			helpMatrix.findSlauSolution(freeVec, solutionVec, true);
			if (!solutionVec->checkZeroVector())
				{
					ownVectors.setCol(ownVectorMatrixIndex, *solutionVec);
					ownVectorMatrixIndex++;

				}
			solutionVec->makeZeroVec();

			
		}
		
		 delete solutionVec;
		 solutionVec = nullptr;
		return true;
	};
	bool LU_decomposition(SquareMatrix<T>& Lower, SquareMatrix<T>& Upper) const
	{
		Upper = *this;
		for (register int i = 0; i < Lower.rows; i++)
		{
			for (register int j = i; j < Lower.cols; j++)
			{
				Lower(i, j) = 0;
				Lower(i, i) = 1;
			}
		}
		float mCoef = 0;
		for (register int j = 0; j < Upper.cols - 1; j++)
		{
			if (Upper(j, j) == 0)
			{
				cout << "Can not find LU decomposition of this matrix!" << endl;
				return false;
			}
			for (register int i = j + 1; i < Upper.rows; i++)
			{
				mCoef = Upper(i, j) / Upper(j, j);
				Lower(i, j) = mCoef;
				Vector<T> vec1(this->cols);
				Vector<T> vec2(this->cols);
				Upper.getRow(j, vec1);
				Upper.getRow(i, vec2);
				vec2 = vec2 - mCoef * vec1;
				Upper.setRow(i, vec2);
			}
		}
		return true;
	};
	bool QR_decomposition(SquareMatrix<T>& Q, SquareMatrix<T>& R) const
	{
		float r_coef;
		Vector<T> helpVec1(this->cols);
		Vector<T> helpVec2(this->cols);
		Vector<T> sumVec(this->cols);
		this->getCol(0, helpVec1);
		float norm = helpVec1.findNormL2();
		if (norm < eps)
		{
			cout << "Can not find QR decomposition for this matrix!" << endl;
			return false;
		}
		helpVec1 /= norm;
		Q.setCol(0, helpVec1);
		for (register int i = 0; i < R.rows; i++)
		{
			for (register int j = 0; j < R.cols; j++)
			{
				R(i, j) = 0;
			}
			R(i, i) = 1;
		}
		for (register int i = 1; i < Q.cols; i++)
		{
			sumVec.makeZeroVec();
			for (register int j = 0; j < i; j++)
			{
				this->getCol(i, helpVec1);
				Q.getCol(j, helpVec2);
				r_coef = helpVec1.scalarMultiplication(helpVec2) / helpVec2.scalarMultiplication(helpVec2);
				helpVec2 *= r_coef;
				sumVec += helpVec2;
			}
			helpVec1 -= sumVec;
			norm = helpVec1.findNormL2();
			if (norm < eps)
			{
				cout << "Can not find QR decomposition for this matrix!" << endl;
				return false;
			}
			helpVec1 /= norm;
			Q.setCol(i, helpVec1);
		}
		SquareMatrix<T> helpMatrix(Q);
		helpMatrix.transpMatrix();
		helpMatrix *= (*this);
		R = helpMatrix;
		return true;
	};
	bool findReversedMatrix()
	{
		if (fabs(this->findDeterminant()) < eps)
		{
			cout << "This matrix can not be reversed!" << endl;
			return false;
		}
		SquareMatrix<T> helpMatrix(*this);
		Vector<float>* solutionVec = (Vector<float>*)operator new((this->cols * sizeof(Vector<float>)));
		for (register int i = 0; i < this->cols; i++)
			new(solutionVec + i) Vector<float>(this->cols);
		for (register int i = 0; i < this->rows; i++)
		{
			Vector<double> vec(this->cols);
			vec[i] = 1;
			helpMatrix.findSlauSolution(vec, solutionVec);
			solutionVec[this->cols - i - 1] = solutionVec[0];
			if (i != this->rows - 1)
			{
				for (register int j = 0; j < this->rows; j++)
					solutionVec[0][j] = 0;
			}
			helpMatrix = *this;
		}
		for (register int j = 0; j < this->rows; j++)
		{
			this->setCol(j, solutionVec[this->cols - j - 1]);
		}
		for (register int i = 0; i < this->cols; i++)
			solutionVec[i].~Vector();
		delete(solutionVec);
		return true;
	};
	void findTriangleForm()
	{
		for (register int j = 0; j < this->cols - 1; j++)
		{
			Vector<T> vec1(this->cols);
			Vector<T> vec2(this->cols);
			bool f = true;

			if (fabs(this->matrix[j][j]) < 0.0001)
			{
				int i = this->findNoZeroElement(j);
				if (i)
				{
					this->swap(j, i);

				}
				else f = false;
			}
			if (f)
			{

				for (register int i = j + 1; i < this->rows; i++)
				{
					float mCoef = this->matrix[i][j] / this->matrix[j][j];
					Vector<T> vec1(this->cols);
					Vector<T> vec2(this->cols);
					Vector<T> vec3(this->cols);
					this->getRow(j, vec1);;
					this->getRow(i, vec2);
					vec3 = mCoef * vec1;
					vec2 = vec2 - vec3;
					this->setRow(i, vec2);

				}

			}
		}
	};
	bool findJordanForm()
	{
		Vector<int> rankTable(this->cols + 2);
		Vector<int> jordanBlockTable(this->cols);
		Vector<float> ownValVec(this->cols);
		SquareMatrix<T> helpMatrix1(*this);
		SquareMatrix<T> helpMatrix2(this->cols);
		bool areDif = true;
		if (!this->LU_methodOwnVal(ownValVec))
		{
			if (!this->QR_methodOwnVal(ownValVec))
				return false;
		}
		cout << ownValVec << endl;

		for (register int i = 0; i < this->rows; i++)
		{
			for (register int j = 0; j < this->cols; j++)
				helpMatrix2(i, j) = 0;
		}
		for (register int i = 0; i < this->cols-1; i++)
		{
			for (register int j = i+1; j < this->cols; j++)
			{
				if (fabs(ownValVec[i] - ownValVec[j]) < 0.0001)
				{
					areDif = false;
					break;
				}
			}
		}
		if (areDif)
		{
			for (register int i = 0; i < this->cols; i++)
				helpMatrix2(i, i) = ownValVec[i];
			*this = helpMatrix2;
			return true;
		}
		bool checkEqualOwnVal = false;
		int qIndex = 0;
		int jordanBoxIndex = 0;
		for (register int i = 0; i < ownValVec.getSize(); i++)
		{
			checkEqualOwnVal = false;
			for (register int m = 0; m < jordanBoxIndex; m++)
			{
				if (fabs(helpMatrix2(m, m) - ownValVec[i]) <= 0.0001)
				{
					checkEqualOwnVal = true;
					break;
				}
			}
			jordanBoxIndex = 0;
			if (checkEqualOwnVal) continue;
			rankTable[0] = helpMatrix1.findRang();
			for (register int k = 0; k < helpMatrix2.rows; k++)
				helpMatrix1(k, k) -= ownValVec[i];
			helpMatrix1.findTriangleForm();
			rankTable[1] = helpMatrix1.findRang();
			if (rankTable[1] == rankTable[0])
			{
				cout << "Can not find Jourdan form for this matrix!" << endl;
				return false;
			}
			for (register int j = 2; j < this->cols + 2; j++)
			{
				helpMatrix1 *= helpMatrix1;
				rankTable[j] = helpMatrix1.findRang();
			}
			for (register int i = 0; i < this->cols - 1; i++)
			{
				if (rankTable[i + 1] > rankTable[i])
				{
					cout << "Can not find Jourdan form for this matrix!" << endl;
					return false;
				}
			}

			for (register int j = 0; j < this->cols; j++)
			{
				jordanBlockTable[j] = rankTable[j] - 2 * rankTable[j + 1] + rankTable[j + 2];
				jordanBoxIndex += jordanBlockTable[j] * (j + 1);
			}
			for (register int j = qIndex; j < qIndex + jordanBoxIndex; j++)
				helpMatrix2(j, j) = ownValVec[i];
			for (register int j = 0; j < this->cols; j++)
			{
				if (jordanBlockTable[j])
				{
					for (register int v = 1; v <= jordanBlockTable[j]; v++)
					{
						for (register int w = 0; w < j; w++)
						{
							helpMatrix2(qIndex + v, qIndex + v + 1) = 1;
						}
					}
					qIndex += jordanBlockTable[j] * (j + 1);

				}
			}
			helpMatrix1 = *this;
		}
		*this = helpMatrix2;
		return true;
	};
	void transpMatrix()
	{
		SquareMatrix<T> helpmatrix(this->cols);
		for (register int i = 0; i < this->rows; i++)
		{
			for (register int j = 0; j < this->cols; j++)
				helpmatrix(i, j) = this->matrix[j][i];
		}
		(*this) = helpmatrix;
	};
	int  findRang() const 
	{
		SquareMatrix<T> helpMatrix(*this);
		helpMatrix.findTriangleForm();
		return helpMatrix.rows - helpMatrix.findAmountOfZeroRows();
	};
	float findDeterminant() const
	{
		SquareMatrix<T> helpMatrix(*this);
		helpMatrix.findTriangleForm();
		float determinant = 1;
		for (register int i = 0; i < this->cols; i++)
		{
			determinant *= helpMatrix(i, i);
		}
		if (!flag)
			determinant *= -1;
		return determinant ;

	};
	

public:
	SquareMatrix( const int& n) :SimpleMatrix<T>(n,n) {};
	SquareMatrix(const SquareMatrix<T>& mat) :SimpleMatrix<T>(mat) {};
	SquareMatrix(const T* const* const arr, const int& size) :SimpleMatrix<T>(arr, size, size) {};
	~SquareMatrix() override { };
	bool flag = true;

	void findMatrixRang() const override
	{
		cout << "Matrix Rang is " << this->findRang() << endl;
	}
	template<class T1>
	void findSlauSolution(Vector<T1>& freeVec, Vector<float>* solutionVec, bool forOwnVec = false) const
	{
		SquareMatrix<T> helpMatrix(*this);
		helpMatrix.findTriangleForm(freeVec);
		if (!forOwnVec)
		{
			int amountOfZeroRows = helpMatrix.findAmountOfZeroRows();
			if (amountOfZeroRows == 0)
			{

				solutionVec[0][this->cols - 1] = freeVec[this->cols - 1] / helpMatrix(helpMatrix.rows - 1, helpMatrix.rows - 1);
				for (register int i = helpMatrix.rows - 2; i >= 0; i--)
				{
					for (register int j = helpMatrix.rows - 1; j > i; j--)
					{

						solutionVec[0][i] -= solutionVec[0][j] * helpMatrix(i, j);
					}
					if (helpMatrix(i, i) == 0)
						solutionVec[0][i] = 1;
					else
						solutionVec[0][i] = (freeVec[i] + solutionVec[0][i]) / helpMatrix(i, i);
				}
			}
			else if (amountOfZeroRows > 0)
			{
				
				bool f = true;
				for (register int i = this->cols - amountOfZeroRows; i < this->cols; i++)
				{
					if (freeVec[i] != 0)
					{
						f = false;
						break;
					}
				}
				if (!f)
				{
					cout << "System has no solutions!" << endl;
					return;
				}
				else
				{
					cout << "System has a lot of solutions!Now we will find foundamental system of solutions! There are " << amountOfZeroRows+1 << " solution vectors" << endl;
					int k = amountOfZeroRows;
					if (amountOfZeroRows == this->rows)
					{
						
							cout << "Each combination of numbers is solution!" << endl;
							for (register int i = 0; i < this->cols; i++)
								solutionVec[i][i] = 1;

							return;
						
						
					}
					int c = 1;
					while (c <= k)
					{

						solutionVec[c - 1][this->cols - c] = 1;

						for (register int i = this->cols - k - 1; i >= 0; i--)
						{
							for (register int j = this->cols - 1; j > i; j--)
							{

								solutionVec[c - 1][i] -= solutionVec[c - 1][j] * helpMatrix(i, j);
							}
							if (helpMatrix(i, i) == 0)
								solutionVec[c - 1][i] = 1;
							else
								solutionVec[c - 1][i] = (freeVec[i] + solutionVec[c - 1][i]) / helpMatrix(i, i);
						}
						c++;

					}
					for (register int i = 1; i <= k; i++)
						solutionVec[k][this->cols-i] = 0;
					for (register int i = this->cols - k - 1; i >= 0; i--)
					{
						for (register int j = this->cols - 1; j > i; j--)
						{

							solutionVec[c - 1][i] -= solutionVec[c - 1][j] * helpMatrix(i, j);
						}
						if (helpMatrix(i, i) == 0)
							solutionVec[c - 1][i] = 1;
						else
							solutionVec[c - 1][i] = (freeVec[i] + solutionVec[c - 1][i]) / helpMatrix(i, i);
					}

				}
			}
		}
		else if (forOwnVec)
		{
			solutionVec[0][this->cols - 1] = 1;
			for (register int i = helpMatrix.rows - 2; i >= 0; i--)
			{
				for (register int j = helpMatrix.rows - 1; j > i; j--)
				{

					solutionVec[0][i] -= solutionVec[0][j] * helpMatrix(i, j);
				}
				if (helpMatrix(i, i) == 0)
					solutionVec[0][i] = 1;
				else
					solutionVec[0][i] = (freeVec[i] + solutionVec[0][i]) / helpMatrix(i, i);
			}
		}

	};
	void  findMatrixDeterminant() const override
	{
		
		cout<<"Determinant is "<< this->findDeterminant()<<endl;

	};

	void findConditionNumber() const override
	{
		SquareMatrix<T> helpMatrix(*this);
		helpMatrix.findReversedMatrix();
		cout<<"Condition number is "<< this->findNormL2() * helpMatrix.findNormL2()<<endl;
	};
	
	void findOwnValues() const override
	{
		Vector<T> ownValVec1(this->cols);
		Vector<T> ownValVec2(this->cols);
		if (this->LU_methodOwnVal(ownValVec1))
		{
			cout << "Own values found with LU method:\n";
			cout << ownValVec1 << endl;
		}
		if (this->QR_methodOwnVal(ownValVec2))
		{
			cout << "Own values found with QR method:\n";
			cout << ownValVec2 << endl;
		}
	};
	void fullSpectralCase() const override
	{
		Vector<float> ownValVec(this->cols);
		Vector<float> helpVector(this->cols);
		SquareMatrix<float> ownVectors(this->cols);
		if (this->fullSpectralCaseLU(ownValVec, ownVectors))
		{
			cout << "FULL SPECTRAL CASE LU IS:\n";
			cout << "Own values are:" << ownValVec << endl;
			for (register int i = 0; i < this->cols; i++)
			{
				ownVectors.getCol(i, helpVector);
				cout << "For own value " << ownValVec[i] << " own vector is " << helpVector << endl;
			}
		}
		ownValVec.makeZeroVec();
		ownVectors.makeZeroMatrix();
		if (this->fullSpectralCaseQR(ownValVec, ownVectors))
		{
			cout << "FULL SPECTRAL CASE QR IS:\n";
			cout << "Own values are:" << ownValVec << endl;
			for (register int i = 0; i < this->cols; i++)
			{
				ownVectors.getCol(i, helpVector);
				cout << "For own value " << ownValVec[i] << " own vector is " << helpVector << endl;
			}
		}
	};
	void findLU() const override
	{
		SquareMatrix<T> helpMatrix(*this);
		SquareMatrix<T> matr1(this->cols);
		SquareMatrix<T> matr2(this->cols);
		if (helpMatrix.LU_decomposition(matr1, matr2))
			cout << endl << "For matrix " << endl << *this << endl << "Lower=" << endl << matr1 << endl << "Upper=" << endl << matr2 << endl;
		else return;
	};
	void findQR() const override
	{
		SquareMatrix<T> helpMatrix(*this);
		SquareMatrix<T> matr1(this->cols);
		SquareMatrix<T> matr2(this->cols);
		if (helpMatrix.QR_decomposition(matr1, matr2))
			cout << endl << "For matrix " << endl << *this << endl << "Q=" << endl << matr1 << endl << "R=" << endl << matr2 << endl;
		else return;
	};
	void findReversed() const override
	{
		SquareMatrix<T> helpMatrix(*this);
		if (helpMatrix.findReversedMatrix())
			cout << endl << "Reversed for matrix " << endl << *this << endl << " is matrix " << endl << helpMatrix << endl;
		else return;
	};
	void transpMatrixForm() const override
	{
		SquareMatrix<T> helpMatrix(*this);
		helpMatrix.transpMatrix();
		cout << endl << "Transp matrix for " << endl << *this << endl << " is " << endl << helpMatrix << endl;

	};
	void findJordan() const
	{
		SquareMatrix<T> helpMatrix(*this);
		if (helpMatrix.findJordanForm())
			cout << "Jordan form for matrix " << endl << *this << endl << " is " << endl << helpMatrix << endl;
		else return;
	};
	void findTriangle() const
	{
		SquareMatrix<T> helpMatrix(*this);
		helpMatrix.findTriangleForm();
		cout << endl << "Triangle form for matrix " << endl << *this << endl << " is " << endl << helpMatrix << endl;
	};
	void findMatrixNormL1() const override {
		cout << "L1 norm of this matrix is " << this->findNormL1() << endl;
	};
	void findMatrixNormL2() const  override{
		cout << "L2 norm of this matrix is " << this->findNormL2() << endl;
	};
	using SimpleMatrix<T>:: operator *=;
	using SimpleMatrix<T>:: vectorMultiplication;
	using SimpleMatrix<T>::Pow;
	using SimpleMatrix<T>::operator ();
	using SimpleMatrix<T>::getColSize;
};
#endif
