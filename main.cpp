
#include <iostream>
#include <map>
#include <complex>
#include <omp.h>

using namespace std;


#include "antenna.h"


double f = 35.22e9;
double omega = 2 * pi_ * f;

const Point p0(0, 0, 15.2);
const Point p1(-0.14, 0, 15.2);
const Point p2(0.14, 0, 15.2);
const Point p3(0, -0.14, 15.2);
const Point p4(0, 0.14, 15.2);

void print_point_field_information(const Point p) {
	printf("Coordinate: (%f, %f, %f)\n", p(0), p(1), p(2));
	printf("Ex: (%9e, %9e)\n", E_strength(p)(0).real(), E_strength(p)(0).imag());
	printf("Ey: (%9e, %9e)\n", E_strength(p)(1).real(), E_strength(p)(1).imag());
	printf("Ez: (%9e, %9e)\n", E_strength(p)(2).real(), E_strength(p)(2).imag());
	printf("Hx: (%9e, %9e)\n", H_strength(p)(0).real(), H_strength(p)(0).imag());
	printf("Hy: (%9e, %9e)\n", H_strength(p)(1).real(), H_strength(p)(1).imag());
	printf("Hz: (%9e, %9e)\n", H_strength(p)(2).real(), H_strength(p)(2).imag());
	printf("\n");
}


int main()
{
	//freopen("out.txt", "w", stdout);
	load_gmsh_file("conical_array.msh");

	set_radiate_model(omega);
	//print_feed_information();
	//print_shared_edge();
	solve_moment_equation();
	//cout << near_field_radius() << endl;


	//load_current_distribution();

	directional_radiate_pattern_XZ();
	//directional_radiate_pattern_XZ();

	/*print_point_field_information(p0);
	print_point_field_information(p1);
	print_point_field_information(p2);
	print_point_field_information(p3);
	print_point_field_information(p4);*/

	//cout << norm(P_strength(p1)) << endl;
	//cout << norm(P_strength(p2)) << endl;
	//cout << norm(P_strength(p3)) << endl;
	//cout << norm(P_strength(p4)) << endl;

	//cout << E_strength(p1).real().norm() << "+j*" << E_strength(p1).imag().norm() << endl;
	//cout << E_strength(p2).real().norm() << "+j*" << E_strength(p2).imag().norm() << endl;
	//cout << E_strength(p3).real().norm() << "+j*" << E_strength(p3).imag().norm() << endl;
	//cout << E_strength(p4).real().norm() << "+j*" << E_strength(p4).imag().norm() << endl;

	special_output();

	//save_current_distribution();
	system("pause");
}
