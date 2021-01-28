#pragma once
#ifndef _IMATRIX_H
#define _IMATRIX_H

class IMatrix {
public:
	virtual void  findMatrixDeterminant() const  = 0;
	virtual void findMatrixRang() const=0;
	virtual void transpMatrixForm() const  = 0;
	virtual void findConditionNumber() const = 0;
	virtual void findReversed() const = 0;
	virtual void findTriangle() const = 0;
	virtual void findOwnValues() const=0;
	virtual void fullSpectralCase() const = 0;
	virtual void findLU() const = 0;
	virtual void findQR() const = 0;
	virtual void findJordan() const = 0;
	virtual void findMatrixNormL1() const = 0;
	virtual void findMatrixNormL2() const = 0;
};
#endif