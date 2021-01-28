#pragma once
#ifndef _CONSOLE_INTERFACE_
#define _CONSOLE_INTERFACE_
#include "SquareMatrix.h"
#include <fstream>
class ConsoleInterface {
public:
	template <class T>
	static void checkAccurancy(T& val, const int& lowVal = -30000, const int& highVal = 30000)
	{
		while (cin.fail() || val < lowVal || val>highVal)
		{
			cout << "Wrong choise number!Try new one" << endl;
			cin.clear();
			cin.ignore(32767, '\n');
			cin >> val;
		}
	};
	static void randomMatrix(SquareMatrix<float>& mat)
	{
		srand((unsigned)time(NULL));
		for (register int i = 0; i < mat.getColSize(); i++)
		{
			for (register int j = 0; j < mat.getColSize(); j++)
				mat(i,j) = rand() % 100 - 50;
		}
	};
	static void fillMatrix(SquareMatrix<float>& mat)
	{
		for (register int i = 0; i < mat.getColSize(); i++)
		{
			for (register int j = 0; j < mat.getColSize(); j++)
			{
				cout << "matrix[" << i << "][" << j << "]=";
				cin >> mat(i,j);
				checkAccurancy(mat(i,j));
			}

		}
	};
	static void fillFreeVec(Vector<float>& freeVec)
	{
		for (register int i = 0; i < freeVec.getSize(); i++)
		{
			cout << "FreeVec [" << i << "] =";
			cin >> freeVec[i];
			checkAccurancy(freeVec[i]);

		}
	};
	static void printChoiseTable()
	{
		cout <<
			"Now you will be able to manipulate with your matrix." << endl << " Press 1 to find determinant;" <<
			endl << "Press 2 to find rang;" << endl << "Press 3 to find transp matrix;" << endl <<
			"Press 4 to find condition number;" << endl << "Press 5 to find reversed matrix;" <<
			endl << "Press 6 to find triangle form;" <<
			endl << "Press 7 to find own values;" <<
			endl << "Press 8 to find full spectral case solution;" <<
			endl << "Press 9 to find LU decomposition;" <<
			endl << "Press 10 to find QR decomposition;" <<
			endl << "Press 11 to find Jordan form;" <<
			endl << "Press 12 to find SLAU solution with this matrix;" <<
			endl << "Press 13 to find L1 norm of  this matrix;" <<
			endl << "Press 14 to find L2 norm of  this matrix;" <<
			endl << "Press 15 to find matrix on matrix multiplication;" <<
			endl << "Press 16 to find matrix on vector multiplication;" <<
			endl << "Press 17 to find pow for matrix;" <<
			endl << "Press 0 to exit" << endl
			<< "Press -1 to change matrix" << endl;
	};
	static int readFromFile(ifstream& in)
	{
		int count = 0;
		int temp;

		while (!in.eof())
		{
			in >> temp;
			count++;
		}
		in.seekg(0, ios::beg);
		in.clear();
		int count_space = 0;
		char symbol;
		while (!in.eof())
		{
			
			in.get(symbol);
			if (symbol == ' ') count_space++;
			if (symbol == '\n') break;
		}
		in.seekg(0, ios::beg);
		in.clear();
		return count / (count_space + 1);
	}
	static void realizationProg()
	{
		int n;
		SquareMatrix<float>* mainObj = static_cast<SquareMatrix<float>*>(operator new[](sizeof(SquareMatrix<float>)));
		cout << "Choose 1 to get matrix from file an 0 to generate matrix in programm: ";
		int isFromFile;
		cin >> isFromFile;
		checkAccurancy(isFromFile, 0, 1);
		if (isFromFile)
		{
			ifstream in("D:\\Проекты 2 курс\\курсач\\OOPkurs\\OOP.txt");
			if (!in.is_open())
			{
				cout << "Can not open the file!\n";
			}
			n = readFromFile(in);
			SquareMatrix<float> mainObj1(n);

			for (register int i = 0; i < n; i++)
			{
				for (register int j = 0; j < n; j++)
				{
					in >> mainObj1(i, j);
				}
			}
			new(mainObj) SquareMatrix<float>(mainObj1);
		}
		else {
			cout << "Please choose dimension of your matrix!" << endl;
			cin >> n;
			checkAccurancy(n, 2, 10);
			system("cls");
			cout << "Dimension is " << n << endl;
			cout << "Press 0 if you want to fill matrix automatically or 1 if you want to fill matrix by your own " << endl;
			int waytoCreate;
			cin >> waytoCreate;
			checkAccurancy(waytoCreate, 0, 1);
			system("cls");
			SquareMatrix<float> mainObj1(n);

			if (waytoCreate == 1) {
				cout << "Filling the matrix" << endl; fillMatrix(mainObj1);
			}
			else if (waytoCreate == 0) {
				cout << "Creating random matrix!" << endl; randomMatrix(mainObj1);
			}
			new(mainObj) SquareMatrix<float>(mainObj1);
		
		}
		system("cls");
		cout << "Your matrix is:" << endl << *mainObj;
		cout << endl;
		printChoiseTable();
		int choise;
		cin >> choise;
		checkAccurancy(choise, -1, 17);

		system("cls");
		Vector<float>* solVec = (Vector<float>*)operator new((n * sizeof(Vector<float>)));
		for (register int i = 0; i < n; i++)
			new(solVec + i) Vector<float>(n);
		Vector<float> freeVector(n);
		Vector<float> multVector(n);
		Vector<float> resultVector(n);
		SquareMatrix<float> helpObj(n);
		SquareMatrix<float> multObj(n);
		while (choise)
		{
			switch (choise)
			{
			case 1: cout << "Your matrix is " << endl << *mainObj << endl; 
				cout << endl << "........................................................................\n"; 
				mainObj->findMatrixDeterminant();
				cout << endl << "........................................................................\n"; 
				break;
			case 2:cout << "Your matrix is " << endl << *mainObj << endl;
				cout << endl << "........................................................................\n"; 
				mainObj->findMatrixRang();
				cout << endl << "........................................................................\n"; 
				break;
			case 3: cout << "Your matrix is " << endl << *mainObj << endl;
				cout << endl << "........................................................................\n"; 
				mainObj->transpMatrixForm();
				cout << endl << "........................................................................\n"; 
				break;
			case 4: cout << "Your matrix is " << endl << *mainObj << endl;
				cout << endl << "........................................................................\n"; 
				mainObj->findConditionNumber();
				cout << endl << "........................................................................\n"; 
				break;
			case 5: cout << "Your matrix is " << endl << *mainObj << endl;
				cout << endl << "........................................................................\n"; 
				mainObj->findReversed();
				cout << endl << "........................................................................\n"; 
				break;
			case 6: cout << "Your matrix is " << endl << *mainObj << endl;
				cout << endl << "........................................................................\n"; 
				mainObj->findTriangle();
				cout << endl << "........................................................................\n"; 
				break;
			case 7: cout << "Your matrix is " << endl << *mainObj << endl;
				cout << endl << "........................................................................\n"; 
				mainObj->findOwnValues();
				cout << endl << "........................................................................\n"; 
				break;
			case 8: cout << "Your matrix is " << endl << *mainObj << endl;
				cout << endl << "........................................................................\n"; 
				mainObj->fullSpectralCase();
				cout << endl << "........................................................................\n"; 
				break;
			case 9: cout << "Your matrix is " << endl << *mainObj << endl;
				cout << endl << "........................................................................\n"; 
				mainObj->findLU();
				cout << endl << "........................................................................\n"; 
				break;
			case 10: cout << "Your matrix is " << endl << *mainObj << endl;
				cout << endl << "........................................................................\n"; 
				mainObj->findQR();
				cout << endl << "........................................................................\n"; 
				break;
			case 11:cout << "Your matrix is " << endl << *mainObj << endl;
				cout << endl << "........................................................................\n"; 
				mainObj->findJordan();
				cout << endl << "........................................................................\n"; 
				break;
			case 12: cout << "Your matrix is " << endl << *mainObj << endl;
				fillFreeVec(freeVector);
				mainObj->findSlauSolution(freeVector, solVec);
				cout << endl << "........................................................................\n";
				for (register int i = 0; i < n; i++)
				{
					if (!solVec[i].checkZeroVector())
					{
						cout << "Solution is: " << solVec[i] << endl;
					}
				}
				cout << endl << "........................................................................\n";
				break;
			case 13: cout << "Your matrix is " << endl << *mainObj << endl;
				cout << endl << "........................................................................\n"; 
				mainObj->findMatrixNormL1();
				cout << endl << "........................................................................\n"; 
				break;
			case 14: cout << "Your matrix is " << endl << *mainObj << endl;
				cout << endl << "........................................................................\n"; 
				mainObj->findMatrixNormL2();
				cout << endl << "........................................................................\n"; 
				break;
			case 15: cout << "Your matrix is " << endl << *mainObj << endl;
				helpObj = *mainObj;
				cout << endl << "........................................................................\n";
				cout << "Please enter multiplication matrix!" << endl;
				fillMatrix(multObj);
				helpObj *= multObj;
				cout << endl << "........................................................................\n";
				cout << "Result matrix is \n" << helpObj << endl;
				cout << endl << "........................................................................\n";
				break;
			case 16: cout << "Your matrix is " << endl << *mainObj << endl;
				fillFreeVec(multVector);
				resultVector = mainObj->vectorMultiplication(multVector);
				cout << endl << "........................................................................\n";
				cout << "Result is " << resultVector << endl;
				resultVector.makeZeroVec();
				cout << endl << "........................................................................\n";
				break;
			case 17: cout << "Your matrix is " << endl << *mainObj << endl;
				cout << "Enter pow degree:" << endl;
				int degree;
				cin >> degree;
				checkAccurancy(degree, 0, 10);
				cout << endl << "........................................................................\n";
				if (degree == 0) { cout << SquareMatrix<float>(n) << endl; break; }
				else if (degree > 0) { helpObj = *mainObj; helpObj.Pow(n); cout << "Result is \n" << helpObj;
				cout << endl << "........................................................................\n"; 
				break; }
			case -1: cout << "Your matrix is " << endl << *mainObj << endl;
				cout << endl << "........................................................................\n"; 
				cout << "Choose row, col for the element you want to change:" << endl;
				int i, j; float val;
				cout << "i= ";
				cin >> i;
				checkAccurancy(i, 0, n-1);
				cout << endl;
				cout << "j= ";
				cin >> j;
				checkAccurancy(j, 0, n-1);
				cout << endl;
				cout << "Choose exchanging value:" << endl;
				cout << "val= ";
				cin >> val;
				checkAccurancy(val);
				changeElement(*mainObj, i, j, val);
				cout << endl << "........................................................................\n";
				break;
			case 0: system("cls"); return ;

			}
			printChoiseTable();
			cin >> choise;
			checkAccurancy(choise, -1, 17);
			system("cls");
		}

		cout << "Goodbye!" << endl;
		for (register int i = 0; i < n; i++)
			solVec[i].~Vector();
		delete(solVec);
		mainObj->~SquareMatrix();
	};
	static void changeElement(SquareMatrix<float>& obj, const int& i, const int& j, const float& val)
	{
		obj(i, j) = val;
	};
	
};
#endif
