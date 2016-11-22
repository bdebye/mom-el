
#include "mathlib.h"
#include <iostream>

#include <omp.h>

VectorXd linspace(double a, double b, int n) {
	VectorXd vect(n + 1);
	double h = (b - a) / n;
	double pos = a;
	for (int i = 0; i <= n; i++) {
		vect(i) = pos;
		pos += h;
	}
	return vect;
}

MatrixXd transpose(const MatrixXd & matrix) {
	return matrix.transpose();
}

MatrixXd transpose(const VectorXd & vect) {
	int leng = vect.size();
	MatrixXd vvec(leng, 1);
	for (int i = 0; i < leng; i++)
		vvec(i, 0) = vect(i);
	return vvec;
}

MatrixXcd transpose(const VectorXcd &vect) {
	int leng = vect.size();
	MatrixXcd vvec(leng, 1);
	for (int i = 0; i < leng; i++)
		vvec(i, 0) = vect(i);
	return vvec;
}

double norm(const Vector & vect) {
	int n = vect.size();
	double result = 0.0;
	for (int i = 0; i < n; i++)
		result += vect(i) * vect(i);
	return sqrt(result);
}

double trian_area(const Triangle & trian) {

	Vector vect1(trian[0] - trian[1]);
	Vector vect2(trian[2] - trian[1]);
	return sqrt(norm(cross(vect1, vect2))) / 2;
}

Point trian_center(const Triangle & trian) {
	return Point((trian[0] + trian[1] + trian[2]) / 3);
}

Point trian_center(const Point &p1, const Point &p2, const Point &p3) {
	return  Point((p1 + p2 + p3) / 3);
}

double distance(const Point & point1, const Point & point2) {
	return norm(Point(point1 - point2));
}

bool diagonally_dominant(const MatrixXcd &Z) {
	bool flag = true;
	for (int i = 0; i < Z.rows(); i++) {
		double diag = std::abs(Z(i, i));
		double remain = 0;
		for (int j = 0; j < i; j++)
			remain += std::abs(Z(i, j));
		for (int j = i + 1; j < Z.cols(); j++)
			remain += std::abs(Z(i, j));
		flag = flag && (diag > remain);
	}
	return flag;
}


double error_estimate(const VectorXcd &I, const VectorXcd &I_old) {
	VectorXcd diff = I - I_old;
	double md = 0;
	for (int i = 0; i < diff.size(); i++) {
		if (std::abs(diff(i)) > md) {
			md = std::abs(diff(i));
		}
	}
	return md;
}

void jacobi_iteration(VectorXcd &I, MatrixXcd &Z, VectorXcd& V) {
	int N = 1000;
	double error;
	std::cout << "DD: " << diagonally_dominant(Z) << std::endl;
	VectorXcd I_old = I;
	for (int i = 0; i < N; i++) {
		#pragma omp parallel for
		for (int j = 0; j < Z.rows(); j++) {
			Complex temp = 0;
			for (int k = 0; k < j; k++)
				temp += Z(j, k) * I_old(k);
			for (int k = j + 1; k < Z.cols(); k++)
				temp += Z(j, k) * I_old(k);
			I(j) = (V(j) - temp) / Z(j, j);
		}
		error = error_estimate(I, I_old);
		I_old = I;
		printf("Iteration number: %d, error estimated: %f.\n", i, error);
	}
}


#include <mkl.h>

void solve_linear(VectorXcd &I, MatrixXcd &Z, VectorXcd& V) {
	I = Z.lu().solve(V);
	//gauss_elimination(I, Z, V);
	//int N = Z.rows();
	//VectorXi ipiv(N);
	//LAPACKE_zgetrf(LAPACK_ROW_MAJOR, Z.rows(), Z.cols(),
	//	(MKL_Complex16 *)Z.data(), Z.rows(), ipiv.data());
	//LAPACKE_zgetrs(LAPACK_ROW_MAJOR, 'N', Z.cols(), 1, (MKL_Complex16 *)Z.data(),
	//	Z.rows(), ipiv.data(), (MKL_Complex16 *)V.data(), V.size());
	//I = V;
}

void gauss_elimination(VectorXcd &I, MatrixXcd &Z, VectorXcd& V) {
	int m = Z.rows();
	int n = Z.cols() + 1;

	if (Z.rows() != Z.cols())
		std::cout << "Not a square matrix, linear system cannot be solved..." << std::endl;

	//Array<Complex, 2> extend(m, n);

	//extend(Range::all(), Range(0, n - 2)) = Z;
	//extend(Range::all(), Range(n - 1, n - 1)) = transpose(V);

	for (int i = 0; i < m; i++) {
		printf("Elimination count %d: %d\n", i, m);
		#pragma omp parallel for
		for (int j = i + 1; j < m; j++) {
			Complex  scala = Z(j, i) / Z(i, i);
			for (int k = i; k < m; k++) {
				Z(j, k) = Z(j, k) - scala * Z(i, k);
			}

			//Z(Range(j, j), Range(i, m - 1))
			//	= Z(Range(j, j), Range(i, m - 1)) - scala * Z(Range(i, i), Range(i, m - 1));
			V(j) = V(j) - scala * V(i);
		}
	}

	// cout << extend << endl;
	for (int i = m - 1; i >= 0; i--) {
		Complex value = Complex(0);
		double value_real = 0;
		double value_imag = 0;
		//#pragma omp parallel for
		printf("Back-substitute count %d: %d\n", m - i, m);
		#pragma omp parallel for reduction(+:value_real)
		for (int j = i + 1; j < m; j++) {
			Complex accu = I(j) * Z(i, j);
			value_real += accu.real();
		}
		#pragma omp parallel for reduction(+:value_imag)
		for (int j = i + 1; j < m; j++) {
			Complex accu = I(j) * Z(i, j);
			value_imag += accu.imag();
		}
		value = Complex(value_real, value_imag);
		I(i) = (V(i) - value) / Z(i, i);
	}
}

Point edge_center(const Point & point1, const Point & point2) {
	return Point((point1 + point2) / 2);
}

Complex exp(Complex comp) {
	double real = comp.real();
	double image = comp.imag();
	return Complex(std::exp(real) * cos(image), std::exp(real) * sin(image));
}

Vector_c exp(Vector_c comp) {
	Vector_c result;
	for (int i = 0; i < 3; i++)
		result(i) = exp(comp(i));
	return result;
}

Complex operator/(const Complex & c1, const Complex & c2) {
	Complex divid = c1 * Complex(c2.real(), -c2.imag());
	return divid / (square(c2.real()) + square(c2.imag()));
}

Complex transv(const Vector_c & v1, const Vector & v2) {
	return (v1(0) * v2(0) + v1(1) * v2(1) + v1(2) * v2(2));
}

Complex transv(const Vector_c & v1, const Vector_c & v2) {
	return (v1(0) * v2(0) + v1(1) * v2(1) + v1(2) * v2(2));
}

Vector_c expand(const Vector &v) {
	Vector_c result;
	for (int i = 0; i < 3; i++) {
		result(i) = Complex(v(i), 0);
	}
	return result;
}

Vector_c conj(const Vector_c &vect_c) {
	Vector_c cj = vect_c;
	cj(0) = std::conj(cj(0));
	cj(1) = std::conj(cj(1));
	cj(2) = std::conj(cj(2));
	return cj;
}

Vector real(const Vector_c &vect_c) {
	Vector rl;
	rl(0) = std::real(vect_c(0));
	rl(1) = std::real(vect_c(1));
	rl(2) = std::real(vect_c(2));
	return rl;
}
