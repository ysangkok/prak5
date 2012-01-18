/*
 * Matrix.h
 *
 * A template matrix class.
 *
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <cstdlib>
#include <cstdio>
#include <iostream>
using namespace std;

template<class T> class Matrix {
public:
	/**
	 * Initializes the matrix with zeros.
	 */
	Matrix(int xWidth, int yWidth);
	//copy constructor
	Matrix(const Matrix<T>& other);
	//copy elements into matrix
	void copyElements(const T* source);
	//assignment operator
	Matrix<T>& operator=(const Matrix<T>& other);
	//destructor
	virtual ~Matrix();
	//returns the element (x,y)
	T at(int x, int y);
	//sets the element (x,y) to the given value
	void set(int x, int y, const T& value);
	//copies the values from the lower into the upper diagonal matrix
	void reflectOnMainDiagonal();
	//prints the matrix (e.g. to the console or to a file)
	void print(ostream& os, bool lowerTriangle);
	void fastprint(ostream& os, bool lowerTriangle);
	
	int size();
private:
	int xWid, yWid;
	T* dat;
};

template<class T> Matrix<T>::Matrix(int xWidth, int yWidth) :
	xWid(xWidth), yWid(yWidth){
	dat = new T[xWid * yWid];
	fill(dat, dat + xWid * yWid, 0);
}

template<class T> Matrix<T>::Matrix(const Matrix<T>& other) :
	xWid(other.xWid), yWid(other.yWid){
	dat = new T[xWid * yWid];
	copyElements(other.dat);
}

template<class T> void Matrix<T>::copyElements(const T* source) {
	copy(source, source + xWid * yWid, dat);
}

template<class T> Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other) {
	if (this != &other) {
		xWid = other.xWid;
		yWid = other.yWid;
		T* new_dat = new T[xWid * yWid];
		copy(other.dat, other.dat + xWid * yWid, new_dat);
		delete[] dat;
		dat = new_dat;
	}
	return *this;
}

template<class T> Matrix<T>::~Matrix() {
	delete[] dat;
}

template<class T> inline T Matrix<T>::at(int x, int y) {
	return dat[y * xWid + x];
}

template<class T> void Matrix<T>::set(int x, int y, const T& value) {
	//printf("%d %d %d\n", x, y, xWid);
	dat[y * xWid + x] = value;
}

template<class T> void Matrix<T>::reflectOnMainDiagonal() {
	for (int y = 0; y < yWid; y++)
		for (int x = 0; x <= y; x++)
			set(y, x, at(x, y));
}

template<class T> void Matrix<T>::print(ostream& os, bool lowerTriangle) {
	for (int y = 0; y < yWid; y++)
		for (int x = 0; lowerTriangle ? x <= y : x < xWid; x++)
			os << "x=" << x << " y=" << y << ": " << at(x, y) << "\n";
}

template<class T> void Matrix<T>::fastprint(ostream& os, bool lowerTriangle) {
	if (!lowerTriangle) abort();

	os.write((char*) dat, (xWid*yWid)*(sizeof(float)/sizeof(char)));
}

template<class T> int Matrix<T>::size()
{
	return xWid * yWid;
}

#endif /* MATRIX_H_ */
