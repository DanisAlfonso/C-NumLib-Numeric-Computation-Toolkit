#pragma once

// all the system #include's we'll ever need
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>

namespace nm {
// macro-like inline functions

    template<typename T>
    inline T SQR(const T a) {return a*a;}

    template<typename T>
    inline const T &MAX(const T &a, const T &b)
    {return b > a ? (b) : (a);}

    inline float MAX(const double &a, const float &b)
    {return b > a ? (b) : float(a);}

    inline float MAX(const float &a, const double &b)
    {return b > a ? float(b) : (a);}

    template<typename T>
    inline const T &MIN(const T &a, const T &b)
    {return b < a ? (b) : (a);}

    inline float MIN(const double &a, const float &b)
    {return b < a ? (b) : float(a);}

    inline float MIN(const float &a, const double &b)
    {return b < a ? float(b) : (a);}

    template<typename T>
    inline T SIGN(const T &a, const T &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

    inline float SIGN(const float &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

    inline float SIGN(const double &a, const float &b)
	{return (float) (b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}

    template<typename T>
    inline void SWAP(T &a, T &b)
	{T dum = a; a = b; b = dum;}

// exception handling

#ifndef _USEERRORCLASS_
#define throw(message) \
{printf("ERROR: %s\n     in file %s at line %d\n", message,__FILE__,__LINE__); throw(1);}
#else
struct error {
	char *message;
	char *file;
	int line;
	error(char *m, char *f, int l) : message(m), file(f), line(l) {}
};
#define throw(message) throw(error(message,__FILE__,__LINE__));
void catch(error err) {
	printf("ERROR: %s\n     in file %s at line %d\n",
		err.message, err.file, err.line);
	exit(1);
}
#endif

// Vector and Matrix Classes

#ifdef _USESTDVECTOR_
#define vector Vector
#else

    template <class T>
    class Vector {
    private:
        int nn;	// size of array. upper index is nn-1
        T *v;
    public:
        Vector();
        explicit Vector(int n);		// Zero-based array
        Vector(int n, const T &a);	//initialize to constant value
        Vector(int n, const T *a);	// Initialize to array
        Vector(const Vector &rhs);	// Copy constructor
        Vector & operator=(const Vector &rhs);	//assignment
        typedef T value_type; // make T available externally
        inline T & operator[](const int i);	//i'th element
        inline const T & operator[](const int i) const;
        inline int size() const;
        void resize(int newn); // resize (contents not preserved)
        void assign(int newn, const T &a); // resize and assign a constant value
        ~Vector();
    };

    // nm Vector definitions

    template <typename T>
    Vector<T>::Vector() : nn(0), v(nullptr) {}

    template <typename T>
    Vector<T>::Vector(int n) : nn(n), v(n > 0 ? new T[n] : nullptr) {}

    template <typename T>
    Vector<T>::Vector(int n, const T& a) : nn(n), v(n > 0 ? new T[n] : nullptr)
    {
        for(int i = 0; i < n; i++) v[i] = a;
    }

    template <typename T>
    Vector<T>::Vector(int n, const T *a) : nn(n), v(n > 0 ? new T[n] : nullptr)
    {
        for(int i = 0; i < n; i++) v[i] = *a++;
    }

    template <typename T>
    Vector<T>::Vector(const Vector<T> &rhs) : nn(rhs.nn), v(nn > 0 ? new T[nn] : nullptr)
    {
        for(int i = 0; i < nn; i++) v[i] = rhs[i];
    }

    template <typename T>
    Vector<T> & Vector<T>::operator=(const Vector<T> &rhs)
    // postcondition: normal assignment via copying has been performed;
    //		if Vector and rhs were different sizes, Vector
    //		has been resized to match the size of rhs
    {
        if (this != &rhs)
        {
            if (nn != rhs.nn) {
                if (v != nullptr) delete [] (v);
                nn = rhs.nn;
                v = nn > 0 ? new T[nn] : nullptr;
            }
            for (int i = 0; i < nn; i++)
                v[i] = rhs[i];
        }
        return *this;
    }

    template <typename T>
    inline T & Vector<T>::operator[](const int i)	//subscripting
    {
#ifdef _CHECKBOUNDS_
        if (i < 0 or i >= nn) {
            throw("Vector subscript out of bounds");
        }
#endif
        return v[i];
    }

    template <typename T>
    inline const T & Vector<T>::operator[](const int i) const	//subscripting
    {
#ifdef _CHECKBOUNDS_
        if (i < 0 or i >= nn) {
            throw("Vector subscript out of bounds");
        }
#endif
        return v[i];
    }

    template <typename T>
    inline int Vector<T>::size() const
    {
        return nn;
    }

    template <typename T>
    void Vector<T>::resize(int newn)
    {
        if (newn != nn) {
            if (v != nullptr) delete[] (v);
            nn = newn;
            v = nn > 0 ? new T[nn] : nullptr;
        }
    }

    template <typename T>
    void Vector<T>::assign(int newn, const T& a)
    {
        if (newn != nn) {
            if (v != nullptr) delete[] (v);
            nn = newn;
            v = nn > 0 ? new T[nn] : nullptr;
        }
        for (int i = 0; i < nn; i++) v[i] = a;
    }

    template <typename T>
    Vector<T>::~Vector()
    {
        if (v != nullptr) delete[] (v);
    }

// end of nm Vector definitions

#endif //ifdef _USESTDVECTOR_

template <typename T>
class Matrix {
private:
	int nn;
	int mm;
	T **v;
public:
	Matrix();
	Matrix(int n, int m);			             // Zero-based array
	Matrix(int n, int m, const T &a);	         //Initialize to constant
	Matrix(int n, int m, const T *a);       	 // Initialize to array
	Matrix(const Matrix &rhs);            	     // Copy constructor
	Matrix & operator=(const Matrix &rhs);    	 //assignment
	typedef T value_type;                        // make T available externally
	inline T* operator[](const int i);           //subscripting: pointer to row i
	inline const T* operator[](const int i) const;
	inline int nrows() const;
	inline int ncols() const;
	void resize(int newn, int newm);             // resize (contents not preserved)
	void assign(int newn, int newm, const T &a); // resize and assign a constant value
	~Matrix();
};

template <typename T>
Matrix<T>::Matrix() : nn(0), mm(0), v(nullptr) {}

template <typename T>
Matrix<T>::Matrix(int n, int m) : nn(n), mm(m), v(n > 0 ? new T*[n] : nullptr)
{
	int i, nel = m*n;
	if (v) v[0] = nel>0 ? new T[nel] : nullptr;
	for (i = 1;i < n;i++) v[i] = v[i-1] + m;
}

template <typename T>
Matrix<T>::Matrix(int n, int m, const T &a) : nn(n), mm(m), v(n>0 ? new T*[n] : nullptr)
{
	int i, j, nel = m*n;
	if (v) v[0] = nel>0 ? new T[nel] : nullptr;
	for (i = 1; i < n; i++) v[i] = v[i-1] + m;
	for (i = 0; i < n; i++) for (j = 0; j < m; j++) v[i][j] = a;
}

template <typename T>
Matrix<T>::Matrix(int n, int m, const T *a) : nn(n), mm(m), v(n>0 ? new T*[n] : nullptr)
{
	int i, j, nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : nullptr;
	for (i = 1; i < n; i++) v[i] = v[i-1] + m;
	for (i = 0; i < n; i++) for (j = 0; j < m; j++) v[i][j] = *a++;
}

template <typename T>
Matrix<T>::Matrix(const Matrix &rhs) : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? new T*[nn] : nullptr)
{
	int i, j, nel = mm*nn;
	if (v) v[0] = nel>0 ? new T[nel] : nullptr;
	for (i = 1; i < nn; i++) v[i] = v[i-1] + mm;
	for (i = 0; i < nn; i++) for (j = 0; j < mm; j++) v[i][j] = rhs[i][j];
}

template <typename T>
Matrix<T> & Matrix<T>::operator=(const Matrix<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	if (this != &rhs) {
		int i, j, nel;
		if (nn != rhs.nn or mm != rhs.mm) {
			if (v != nullptr) {
				delete[] (v[0]);
				delete[] (v);
			}
			nn = rhs.nn;
			mm = rhs.mm;
			v = nn > 0 ? new T*[nn] : nullptr;
			nel = mm*nn;
			if (v) v[0] = nel>0 ? new T[nel] : nullptr;
			for (i = 1; i < nn; i++) v[i] = v[i-1] + mm;
		}
		for (i = 0; i < nn; i++) for (j = 0; j < mm; j++) v[i][j] = rhs[i][j];
	}
	return *this;
}

template <typename T>
inline T* Matrix<T>::operator[](const int i)	//subscripting: pointer to row i
{
#ifdef _CHECKBOUNDS_
    if (i < 0 or i>=nn) {
        throw("Matrix subscript out of bounds");
}
#endif
	return v[i];
}

template <typename T>
inline const T* Matrix<T>::operator[](const int i) const
{
#ifdef _CHECKBOUNDS_
    if (i < 0 or i >= nn) {
        throw("Matrix subscript out of bounds");
}
#endif
	return v[i];
}

template <typename T>
inline int Matrix<T>::nrows() const
{
	return nn;
}

template <typename T>
inline int Matrix<T>::ncols() const
{
	return mm;
}

template <typename T>
void Matrix<T>::resize(int newn, int newm)
{
	int i, nel;
	if (newn != nn or newm != mm) {
		if (v != nullptr) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn>0 ? new T*[nn] : nullptr;
		nel = mm*nn;
		if (v) v[0] = nel>0 ? new T[nel] : nullptr;
		for (i = 1; i < nn; i++) v[i] = v[i-1] + mm;
	}
}

template <typename T>
void Matrix<T>::assign(int newn, int newm, const T& a)
{
	int i, j, nel;
	if (newn != nn or newm != mm) {
		if (v != nullptr) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn>0 ? new T*[nn] : nullptr;
		nel = mm*nn;
		if (v) v[0] = nel>0 ? new T[nel] : nullptr;
		for (i = 1; i < nn; i++) v[i] = v[i-1] + mm;
	}
	for (i = 0; i < nn; i++) for (j = 0; j < mm; j++) v[i][j] = a;
}

template <typename T>
Matrix<T>::~Matrix()
{
	if (v != nullptr) {
		delete[] (v[0]);
		delete[] (v);
	}
}

template <typename T>
class Matrix3D {
private:
	int nn;
	int mm;
	int kk;
	T ***v;
public:
	Matrix3D();
	Matrix3D(int n, int m, int k);
	inline T** operator[](const int i);	//subscripting: pointer to row i
	inline const T* const * operator[](const int i) const;
	inline int dim1() const;
	inline int dim2() const;
	inline int dim3() const;
	~Matrix3D();
};

template <typename T>
Matrix3D<T>::Matrix3D(): nn(0), mm(0), kk(0), v(nullptr) {}

template <typename T>
Matrix3D<T>::Matrix3D(int n, int m, int k) : nn(n), mm(m), kk(k), v(new T**[n])
{
	int i,j;
	v[0] = new T*[n*m];
	v[0][0] = new T[n*m*k];
	for(j = 1; j < m; j++) v[0][j] = v[0][j-1] + k;
	for(i = 1; i < n; i++) {
		v[i] = v[i-1] + m;
		v[i][0] = v[i-1][0] + m*k;
		for(j = 1; j < m; j++) v[i][j] = v[i][j-1] + k;
	}
}

template <typename T>
inline T** Matrix3D<T>::operator[](const int i) //subscripting: pointer to row i
{
	return v[i];
}

template <typename T>
inline const T* const * Matrix3D<T>::operator[](const int i) const
{
	return v[i];
}

template <typename T>
inline int Matrix3D<T>::dim1() const
{
	return nn;
}

template <typename T>
inline int Matrix3D<T>::dim2() const
{
	return mm;
}

template <typename T>
inline int Matrix3D<T>::dim3() const
{
	return kk;
}

template <typename T>
Matrix3D<T>::~Matrix3D()
{
	if (v != nullptr) {
		delete[] (v[0][0]);
		delete[] (v[0]);
		delete[] (v);
	}
}

// NaN: uncomment one of the following 3 methods of defining a global NaN
// you can test by verifying that (NaN != NaN) is true

static const double NaN = std::numeric_limits<double>::quiet_NaN();

//Uint proto_nan[2]={0xffffffff, 0x7fffffff};
//double NaN = *( double* )proto_nan;

//Doub NaN = sqrt(-1.);
// complex
typedef std::complex<double> Complex;

// Vector types


typedef const Vector<int> VectorIntI;
typedef Vector<int> VectorInt;
typedef Vector<int> VectorIntO;
typedef Vector<int> VectorIntIO;

typedef const Vector<unsigned int> VectorUsignedIntI;
typedef Vector<unsigned int> VectorUnsignedInt;
typedef Vector<unsigned int> VectorUnsignedIntO;
typedef Vector<unsigned int> VectorUnsignedIntIO;

typedef const Vector<long long int> VectorLLongIntI;
typedef Vector<long long int> VectorLongLongInt;
typedef Vector<long long int> VectorLongLongIntO;
typedef Vector<long long int> VectorLongLongIntIO;

typedef const Vector<unsigned long long int> VectorUnsignedLongLongIntI;
typedef Vector<unsigned long long int> VectorUnsignedLongLongInt;
typedef Vector<unsigned long long int> VectorUnsignedLongLongIntO;
typedef Vector<unsigned long long int> VectorUnsignedLongLongIntIO;

typedef const Vector<char> VectorCharI;
typedef Vector<char> VectorChar;
typedef Vector<char> VectorCharO;
typedef Vector<char> VectorCharIO;

typedef const Vector<char*> VectorCharPI;
typedef Vector<char*> VectorCharP;
typedef Vector<char*> VectorCharPO;
typedef Vector<char*> VectorCharPIO;

typedef const Vector<unsigned char> VectorUnsignedCharI;
typedef Vector<unsigned char> VectorUnsignedChar;
typedef Vector<unsigned char> VectorUnsignedCharO;
typedef Vector<unsigned char> VectorUnsignedCharIO;

typedef const Vector<double> VectorDoubleI;
typedef Vector<double> VectorDouble;
typedef Vector<double> VectorDoubleO;
typedef Vector<double> VectorDoubleIO;

typedef const Vector<double*> VectorPI;
typedef Vector<double*> VectorDoubleP;
typedef Vector<double*> VectorDoublePO;
typedef Vector<double*> VectorDoublePIO;

typedef const Vector<Complex > VectorComplexI;
typedef Vector<Complex> VectorComplex;
typedef Vector<Complex> VectorComplexO;
typedef Vector<Complex> VectorComplexIO;

typedef const Vector<bool> VectorBoolI;
typedef Vector<bool> VectorBool;
typedef Vector<bool> VectorBoolO;
typedef Vector<bool> VectorBoolIO;

// matrix types
typedef const Matrix<int> MatrixIntI;
typedef Matrix<int> MatrixInt;
typedef Matrix<int> MatrixIntO;
typedef Matrix<int> MatrixIntIO;

typedef const Matrix<unsigned int> MatrixUIntI;
typedef Matrix<unsigned int> MatrixUnsignedInt;
typedef Matrix<unsigned int> MatrixUnsignedIntO;
typedef Matrix<unsigned int> MatrixUnsignedIntIO;

typedef const Matrix<long long int> MatrixLongLongIntI;
typedef Matrix<long long int> MatrixLongLongInt;
typedef Matrix<long long int> MatrixLongLongIntO;
typedef Matrix<long long int> MatrixLongLongIntIO;

typedef const Matrix<unsigned long long int> MatrixUnsignedLongLongIntI;
typedef Matrix<unsigned long long int> MatrixUnsignedLongLongInt;
typedef Matrix<unsigned long long int> MatrixUnsignedLongLongIntO;
typedef Matrix<unsigned long long int> MatrixUnsignedLongLongIntIO;

typedef const Matrix<char> MatrixCharI;
typedef Matrix<char> MatrixChar;
typedef Matrix<char> MatrixCharO;
typedef Matrix<char> MatrixCharIO;

typedef const Vector<unsigned char> MatrixUnsignedCharI;
typedef Matrix<unsigned char> MatrixUnsignedChar;
typedef Matrix<unsigned char> MatrixUnsignedCharO;
typedef Matrix<unsigned char> MatrixUnsignedCharIO;


typedef const Matrix<double> MatrixDoubleI;
typedef Matrix<double> MatrixDouble;
typedef Matrix<double> MatrixDoubleO;
typedef Matrix<double> MatrixDoubleIO;

typedef const Matrix<bool> MatrixBoolI;
typedef Matrix<bool> MatrixBool;
typedef Matrix<bool> MatrixBoolO;
typedef Matrix<bool> MatrixBoolIO;

// 3D matrix types
typedef const Matrix3D<double> Matrix3DDoubleI;
typedef Matrix3D<double> Matrix3DDouble;
typedef Matrix3D<double> Matrix3DDoubleO;
typedef Matrix3D<double> Matrix3DDoubleIO;

}

