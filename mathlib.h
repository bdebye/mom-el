

#ifndef MATHLIB_H_
#define MATHLIB_H_

#include "type.h"
using namespace Eigen;

VectorXd linspace(double a, double b, int n);

MatrixXd transpose(const MatrixXd &matrix);

MatrixXd transpose(const VectorXd &vect);
MatrixXcd transpose(const VectorXcd & vect);

double norm(const Vector & vect);

double trian_area(const Triangle & trian);

Point trian_center(const Triangle & trian);
Point trian_center(const Point &p1, const Point &p2, const Point &p3);

double distance(const Point & point1, const Point & point2);

void gauss_elimination(VectorXcd &I, MatrixXcd &Z, VectorXcd &V);

Point edge_center(const Point & point1,
	const Point & point2);



Vector_c exp(Vector_c comp);

Complex operator/(const Complex & c1, const Complex & c2);

inline double square(double v) {
	return v * v;
}

double abs(Complex v);

Complex transv(const Vector_c & v1, const Vector & v2);
Complex transv(const Vector_c & v1, const Vector_c & v2);

Vector_c expand(const Vector &v);

Vector_c conj(const Vector_c &vect_c);

Vector real(const Vector_c &vect_c);

void solve_linear(VectorXcd &I, MatrixXcd &Z, VectorXcd& V);

inline Vector cross(const Vector &v1, const Vector &v2) {
	return (v1.cross(v2));
}

#endif /* MATHLIB_H_ */
