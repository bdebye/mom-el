
#include "antenna.h"
#include "type.h"
#include "mathlib.h"
#include "gmsh.h"
#include "rwg_edge.h"
#include "clock.h"
#include "graphic.h"

#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <complex>
#include <algorithm>
#include <numeric>

#include <omp.h>

using namespace std;

VectorXcd Im;

Matrix<Point, Dynamic, Dynamic> integral_center;

MatrixXcd impedance;

VectorXcd Vm;

Vector_c E_inc;

vector<Feed> feed;

vector<RWGedge> edge_element;

vector<vector<int>> trianP_toRWG;

vector<vector<int>> trianM_toRWG;

int trian_num;

int point_num;

int rwg_num;

double omega_;

Clock local_clock;

const RWGedge& get_RWGedge(int n) {
	return edge_element.at(n);
}

int get_RWGedge_num() {
	return edge_element.size();
}


bool same_edge_index(Index2 idx1, Index2 idx2) {
	return (idx1(0) == idx2(0) && idx1(1) == idx2(1)) ||
		(idx1(0) == idx2(1) && idx1(1) == idx2(0));
}

bool share_edge(int m, int n, Index2 &edge, int &vertex1, int &vertex2) {

	Index3 trian1 = get_triangle_index(m);
	Index3 trian2 = get_triangle_index(n);

	int count = 0;
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++) {
		if (trian1(i) == trian2(j)) {
			edge(count) = trian1(i);
			count++;
		}
	}

	if (count < 2)
		return false;
	//if(edge_feed_considered(m, n, edge))
	//	return false;

	for (int i = 0; i < 3; i++) {
		if (trian1(i) != edge(0) && trian1(i) != edge(1))
			vertex1 = trian1(i);
		if (trian2(i) != edge(0) && trian2(i) != edge(1))
			vertex2 = trian2(i);
	}
	//cout << "--------------------------------" << endl;
	//cout << m << " " << n << endl;
	//cout << trian1 << " " << trian2 << endl;
	//cout << edge << endl;
	return true;
}



void search_RWGedge() {
	trian_num = mesh.triang_num;
	trianP_toRWG.resize(trian_num + 1);
	trianM_toRWG.resize(trian_num + 1);
	std::cout << "generating the RWG edge elements..." << std::endl;
	for (int i = 1; i <= trian_num; i++)
	for (int j = i + 1; j <= trian_num; j++) {
		EdgeIndex edge;
		int vertex1, vertex2;
		if (share_edge(i, j, edge, vertex1, vertex2)) {
			edge_element.push_back(RWGedge(edge_element.size(), edge, vertex1, vertex2, i, j));
			trianP_toRWG.at(i).push_back(edge_element.size() - 1);
			trianM_toRWG.at(j).push_back(edge_element.size() - 1);
			//if(find_feed_edge(i, j) >= 0) {
			//	feed[find_feed_edge(i, j)].rwg_index = edge_element.size() - 1;
			//}
		}
	}

	rwg_num = edge_element.size();
	cout << "finished gerenation, " << get_RWGedge_num()
		<< " edge elements found..." << endl
		<< endl;
}

void generate_impedance_matrix(double omega) {
	//预先计算一些常量
	double c_ = 1 / sqrt(mu_ * epsilon_);
	double k = omega / c_;
	double Constant1 = mu_ / (4 * pi_);
	double Factor = 1.0 / 9;

	Complex K = k * j;
	Complex Constant2 = 1 / (j * 4.0 * pi_ * omega * epsilon_);
	Complex FactorA = Factor * (j * omega / 4) * Constant1;
	Complex FactorFi = Factor * Constant2;

	omega_ = omega;
	int dim = get_RWGedge_num();
	impedance.resize(dim, dim);

	//cout << Factor << endl;
	//cout << (j * omega / 4) << endl;
	//cout << Constant1 << endl;
	//cout << Constant2 << endl;
	//cout << k << endl;
	//cout << FactorA << endl;
	//cout << "-------------------------------------------------------------" << endl;

	//初始化阻抗矩阵
	for (int i = 0; i < rwg_num; i++)
	for (int j = 0; j < rwg_num; j++) {
		impedance(i, j) = Complex(0, 0);
	}
	printf("generating impedance matrix...");
	#pragma omp parallel for
	for (int p = 1; p <= trian_num; p++) {
		vector<int> Plus = trianP_toRWG.at(p);
		vector<int> Minus = trianM_toRWG.at(p);
		Point t_center = get_triangle_center(p);
		/*if (p) {
			printf("generating impedance matrix, finishing %d%%\n",
				(int)((double)(p) / trian_num * 100));
		}*/
		for (int n = 0; n < rwg_num; n++) {
			const RWGedge &edge_n = get_RWGedge(n);
			std::array<Complex, 9> gP;
			std::array<Complex, 9> gM;
			Complex ZF;
			Complex Fi;
			int trian_p = edge_n.trian_plus_index();
			int trian_m = edge_n.trian_minus_index();

			for (int h = 0; h < 9; h++) {
				double R = norm(Vector(integral_center(trian_p, h) - t_center));
				gP[h] = exp(-K * R) / R;
				R = norm(Vector(integral_center(trian_m, h) - t_center));
				gM[h] = exp(-K * R) / R;
			}
			Fi = std::accumulate(gP.begin(), gP.end(), Complex(0, 0)) - std::accumulate(gM.begin(), gM.end(), Complex(0, 0));
			ZF = FactorFi * edge_n.share_edge_length() * Fi;

			for (auto iter = Plus.begin(); iter != Plus.end(); ++iter) {
				int m = *iter;
				const RWGedge &edge_m = get_RWGedge(m);
				int trian_index = edge_m.trian_plus_index();
				Complex A(0, 0);
				for (int h = 0; h < 9; h++) {
					A += gP[h] * edge_n.rho_c_plus().dot(edge_m.rho_plus(integral_center(trian_index, h))) +
						gM[h] * edge_n.rho_c_minus().dot(edge_m.rho_plus(integral_center(trian_index, h)));
				}
				Complex Zl = FactorA * A * edge_n.share_edge_length();

				//cout << Zl << endl;
				impedance(m, n) += edge_m.share_edge_length() * (Zl + ZF);
			}

			for (auto iter = Minus.begin(); iter != Minus.end(); ++iter) {
				int m = *iter;
				const RWGedge &edge_m = get_RWGedge(m);
				int trian_index = edge_m.trian_minus_index();
				Complex A(0, 0);
				for (int h = 0; h < 9; h++) {
					A += gP[h] * edge_n.rho_c_plus().dot(edge_m.rho_minus(integral_center(trian_index, h))) +
						gM[h] * edge_n.rho_c_minus().dot(edge_m.rho_minus(integral_center(trian_index, h)));
				}
				Complex Zl = FactorA * A * edge_n.share_edge_length();
				impedance(m, n) += edge_m.share_edge_length() * (Zl - ZF);
			}
		}
	}
	cout << endl;
	cout << "finished generation, time elapsed: " << local_clock.time_elapse_format() << endl;
	cout << endl;
}

const MatrixXcd & get_impedance_matrix() {
	return impedance;
}

Complex scatter_Vm_value(int m) {
	RWGedge edge_m = get_RWGedge(m);
	Vector rmc_p = edge_m.rho_c_plus();
	Vector rmc_m = edge_m.rho_c_minus();
	double lm = edge_m.share_edge_length();
	return lm * (E_inc.dot(rmc_p) / 2 + E_inc.dot(rmc_m) / 2);
}

void generate_integral_pieces() {
	Triangle trian;
	std::array<Point, 10> pseq;
	Matrix<Point, Dynamic, Dynamic> &cseq = integral_center;
	cseq.resize(trian_num + 1, 9);
	for (int i = 1; i <= trian_num; i++) {
		trian = get_triangle(i);
		pseq[0] = trian[0];
		pseq[1] = trian[0] * 2 / 3 + trian[1] / 3;
		pseq[2] = trian[0] * 2 / 3 + trian[2] / 3;
		pseq[3] = trian[0] / 3 + trian[1] * 2 / 3;
		pseq[5] = trian[0] / 3 + trian[2] * 2 / 3;
		pseq[4] = pseq[3] / 2 + pseq[5] / 2;
		pseq[6] = trian[1];
		pseq[7] = trian[1] * 2 / 3 + trian[2] / 3;
		pseq[8] = trian[1] / 3 + trian[2] * 2 / 3;
		pseq[9] = trian[2];

		cseq(i, 0) = trian_center(pseq[0], pseq[1], pseq[2]);
		cseq(i, 1) = trian_center(pseq[1], pseq[3], pseq[4]);
		cseq(i, 2) = trian_center(pseq[1], pseq[2], pseq[4]);
		cseq(i, 3) = trian_center(pseq[2], pseq[4], pseq[5]);
		cseq(i, 4) = trian_center(pseq[3], pseq[6], pseq[7]);
		cseq(i, 5) = trian_center(pseq[3], pseq[4], pseq[7]);
		cseq(i, 6) = trian_center(pseq[4], pseq[7], pseq[8]);
		cseq(i, 7) = trian_center(pseq[4], pseq[5], pseq[8]);
		cseq(i, 8) = trian_center(pseq[5], pseq[8], pseq[9]);
	}
}


void generate_edge_elements() {
	search_RWGedge();
	generate_integral_pieces();
}

/*
void set_scatter_model(Vector_c Ei, double omega) {
feed.clear();
omega_ = omega;
cout << "using omega: " << omega_ << endl;
generate_edge_elements();
generate_impedance_matrix(omega);
E_inc = Ei;
Vm.resize(rwg_num);
for(int i = 0; i < rwg_num; i++) {
Vm(i) = scatter_Vm_value(i);
}
gauss_elimination(Im, impedance, Vm);
}
*/

void solve_moment_equation() {
	generate_impedance_matrix(omega_);
	//save_impedance_matrix();
	Im.resize(rwg_num);
	cout << "Solving the linear system..." << endl;
	//gauss_elimination(Im, impedance, Vm);
	Clock LU_clock;
	solve_linear(Im, impedance, Vm);
	cout << "Finished solving the linear system, time elapsed: " << local_clock.time_elapse_format() << endl;
	cout << "Time used solving linear system: " << LU_clock.time_elapse_format() << endl;
}

int rwg_index_feed_2() {
	int seg_index;
	for (int i = 0; i < mesh.segment_num; i++) {
		Point center = (get_point(mesh.segment[i](0)) + get_point(mesh.segment[i](1))) / 2;
		if (center(1) > 0) {
			seg_index = i;
		}
	}
	//cout << seg_index << endl;
	for (int i = 0; i < get_RWGedge_num(); i++) {
		Index2 edge_share = get_RWGedge(i).edge_share();
		if (same_edge_index(mesh.segment[seg_index], edge_share)) {
			return i;
		}
	}
	return -1;
}

int rwg_index_feed_3() {
	int seg_index;
	for (int i = 0; i < mesh.segment_num; i++) {
		Point center = (get_point(mesh.segment[i](0)) + get_point(mesh.segment[i](1))) / 2;
		//cout << center << endl << endl;
		if (center(1) < 0 && center(0) > 0) {
			seg_index = i;

		}
	}
	//cout << seg_index << endl;
	for (int i = 0; i < get_RWGedge_num(); i++) {
		Index2 edge_share = get_RWGedge(i).edge_share();
		if (same_edge_index(mesh.segment[seg_index], edge_share)) {
			return i;
		}
	}
	return -1;
}

int rwg_index_feed_1() {
	int seg_index;
	for (int i = 0; i < mesh.segment_num; i++) {
		Point center = (get_point(mesh.segment[i](0)) + get_point(mesh.segment[i](1))) / 2;
		if (center(1) < 0 && center(0) < 0) {
			seg_index = i;
		}
	}
	//cout << seg_index << endl;
	for (int i = 0; i < get_RWGedge_num(); i++) {
		Index2 edge_share = get_RWGedge(i).edge_share();
		if (same_edge_index(mesh.segment[seg_index], edge_share)) {
			return i;
		}
	}
	return -1;
}

void unique_feed(Complex A) {
	int feed_index;
	for (int i = 0; i < get_RWGedge_num(); i++) {
		Index2 edge_share = get_RWGedge(i).edge_share();
		if (same_edge_index(mesh.segment[0], edge_share)) {
			feed_index = i;
		}
	}
	feed.push_back(Feed(feed_index, A));
}

void set_feed_1(Complex a) {
	feed.push_back(Feed(rwg_index_feed_1(), a));
}

void set_feed_2(Complex a) {
	feed.push_back(Feed(rwg_index_feed_2(), a));
}

void set_feed_3(Complex a) {
	feed.push_back(Feed(rwg_index_feed_3(), a));
}

Complex feed_Aphase(double A, double phase) {
	return Complex(A * cos(phase), A * sin(phase));
}

void set_feed_strategy() {
	feed.clear();
	//set_feed_1(feed_Aphase(1, 0));
	//set_feed_2(feed_Aphase(1, 0));
	//set_feed_3(feed_Aphase(1, 0));
	unique_feed(Complex(1, 0));
}

void set_radiate_model(double omega) {
	//E_inc = Complex(0, 0);
	omega_ = omega;
	cout << "using omega: " << omega_ << endl;
	generate_edge_elements();
	Vm.resize(rwg_num);
	for (int i = 0; i < rwg_num; i++) {
		Vm(i) = 0;
	}
	set_feed_strategy();
	for (vector<Feed>::iterator iter = feed.begin(); iter != feed.end(); ++iter)
		Vm(iter->rwg_index) = iter->V * get_RWGedge(iter->rwg_index).share_edge_length();

}

void print_shared_edge() {
	for (unsigned i = 0; i < edge_element.size(); i++) {
		Index2 node_index = edge_element.at(i).edge_share();
		if (get_point(node_index(0))(1) == 0.01 / 100.0) {
			cout << i << ": " << get_point(node_index(0)) << "----" << get_point(node_index(1)) << endl;
		}
	}
}

Complex get_input_impedance(int feed_edge) {
	bool is_feed = false;
	if (feed.empty()) {
		cout << "Sorry, not working in radiate mode... no feed edge specified..." << endl;
		return Complex(0, 0);
	}
	for (auto iter = feed.begin(); iter != feed.end(); ++iter) {
		if (iter->rwg_index == feed_edge)
			is_feed = true;
	}
	if (!is_feed) {
		cout << "Soory, the RWG edge specified is not an feeding port..." << endl;
	}
	double ln = get_RWGedge(feed_edge).share_edge_length();
	return Vm(feed_edge) / (pow(ln, 2) * Im(feed_edge));
}

Complex get_input_power(int feed_edge) {
	bool is_feed = false;
	if (feed.empty()) {
		cout << "Sorry, not working in radiate mode... no feed edge specified..." << endl;
		return Complex(0, 0);
	}
	for (auto iter = feed.begin(); iter != feed.end(); ++iter) {
		if (iter->rwg_index == feed_edge)
			is_feed = true;
	}
	if (!is_feed) {
		cout << "Soory, the RWG edge specified is not an feeding port..." << endl;
	}
	Complex V = 0.0;
	for (auto iter = feed.begin(); iter != feed.end(); ++iter) {
		if (iter->rwg_index == feed_edge)
			V = iter->V;
	}
	return pow(V, 2) / get_input_impedance(feed_edge) / 2;
}

double near_field_radius() {
	double lambda = c_ / omega_ * 2 * pi_;
	return 2 * pow(geometry_size(), 2) / lambda;
}

Vector_c RWG_H_strength(const RWGedge &edge, const Point &r) {
	Vector_c result;
	Vector r_ = r - (edge.minus_center() + edge.plus_center()) / 2;
	Vector_c m = edge.share_edge_length() * Im(edge.index) * (edge.minus_center() - edge.plus_center());
	double c_ = 1 / sqrt(mu_ * epsilon_);
	double rd = norm(r_);
	double k = omega_ / c_;
	Complex C = 1 / pow(rd, 2) * (Complex(1, 0) + 1 / (j * k * rd));
	result = (j * k) / (4 * pi_) * m.cross(expand(r_)) * C * exp(-j * k * rd);

	return result;
}

Vector_c RWG_E_strength(const RWGedge &edge, const Point &r) {
	Vector_c result;
	Vector r_ = r - (edge.minus_center() + edge.plus_center()) / 2;
	double c_ = 1 / sqrt(mu_ * epsilon_);
	double rd = norm(r_);
	double k = omega_ / c_;
	Complex C = 1 / pow(rd, 2) * (Complex(1, 0) + 1 / (j * k * rd));
	Vector_c m = edge.share_edge_length() * Im(edge.index) * (edge.minus_center() - edge.plus_center());
	Vector_c M = r_.dot(m) * r_ / (rd * rd);
	result = eta_ / (4.0 * pi_) * ((M - m) * (j * k / rd + C) + 2.0 * C * M) * exp(-j * k * rd);
	return result;
}

void save_impedance_matrix() {
	const ArrayXXcd &data = get_impedance_matrix();
	int m = data.rows();
	int n = data.cols();
	ofstream fout("Z.txt");
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			fout << std::abs(data(i, j));
			if (j != n - 1)
				fout << ",\t";
		}
		fout << endl;
	}
	fout.close();
}

Vector_c E_strength(const Point &r) {
	Vector_c E(Complex(0, 0), Complex(0, 0), Complex(0, 0));
	//#pragma omp parallel for
	for (int i = 0; i < rwg_num; i++) {
		//#pragma omp critical
		E = E + RWG_E_strength(get_RWGedge(i), r);
	}
	return E;
}

Vector_c H_strength(const Point &r) {
	Vector_c H(Complex(0, 0), Complex(0, 0), Complex(0, 0));
	//#pragma omp parallel for
	for (int i = 0; i < rwg_num; i++) {
		//#pragma omp critical
		H = H + RWG_H_strength(get_RWGedge(i), r);
	}
	return H;
}

Vector W_strength(const Point &r) {
	Vector_c temp = E_strength(r).cross(conj(H_strength(r))) / 2.0;
	return temp.real();
}

Vector_c Poynting(const Point &r) {
	return E_strength(r).cross(conj(H_strength(r))) / 2.0;
}

//Vector P_strength(const Point &r) {
//	return E_strength(r).real().cross(H_strength(r).real());
//}

double phase(const Point &r) {
	Vector_c E = H_strength(r);
	return acos(E.real().norm() / E.norm());
}

double U_density(const Point &r) {
	return norm(W_strength(r)) * pow(norm(r), 2);
}

void directional_radiate_pattern_XZ() {
	int feed_edge = feed[0].rwg_index;
	double R = 100;
	cout << R << endl;
	VectorXd X(1001), Y(1001);
	X = linspace(-0.5 * pi_, 1.5 * pi_, 1000);
	for (int i = 0; i < 1001; i++) {
		Y(i) = U_density(Point(R * cos(X(i)), 0, R * sin(X(i))));
	}

	double scalar = get_input_power(feed_edge).real() / (4 * pi_);
	for (int i = 0; i < 1001; i++) {
		Y(i) = 10 * log10(Y(i) / scalar);
	}
	polar(X, Y);
}
//
void directional_radiate_pattern_XY() {
	int feed_edge = feed[0].rwg_index;
	double R = near_field_radius() * 100;
	VectorXd X(1001), Y(1001);
	X = linspace(0, 2 * pi_, 1000);
	for (int i = 0; i < 1001; i++) {
		Y(i) = U_density(Point(R * cos(X(i)), R * sin(X(i)), 0));
	}

	double scalar = get_input_power(feed_edge).real() / (4 * pi_);
	for (int i = 0; i < 1001; i++) {
		Y(i) = 10 * log10(Y(i) / scalar);
	}
	polar(X, Y);
}

void print_feed_information() {
	for (auto iter = feed.begin(); iter != feed.end(); ++iter) {
		cout << iter->rwg_index << endl;
	}
}

/*
void special_output() {
int nx = 100, ny = 100;
double dx = 0.05 / nx, dy = 0.05 / ny;
double wave_len = c_ / omega_ * 2 * pi_;
double z_coord = 8.3 / 100 + wave_len / 2;
Array<double, 2> field(nx, ny);
cout << "Save field distribution..." << endl;
for(int i = 0; i < nx; i++)
for(int j = 0; j < ny; j++) {
double x = -2.5 / 100 + i * dx;
double y = -2.5 / 100 + j * dy;
double z = z_coord;
double Ex = E_strength(Point(x, y, z))(0).real();
double Ey = E_strength(Point(x, y, z))(1).real();
field(i, j) =  Ex * Ex + Ey * Ey;// norm(real(E_strength(Point(x, y, z))));
}
save_csv(field, "field.txt");
}
*/
//
void special_output() {
	int nx = 100, nz = 300;
	double dx = 0.05 / nx, dz = 0.15 / nz;
	double wave_len = c_ / omega_ * 2 * pi_;
	double y_coord = 0;
	MatrixXd field(nx, nz);
	cout << "Save field distribution..." << endl;
	#pragma omp parallel for
	for(int i = 0; i < nx; i++)
	for(int j = 0; j < nz; j++) {
		double x = 0;
		double y = -0.025 + i * dx;
		double z = -0.15 + j * dz;
		double Ex = E_strength(Point(x, y, z))(0).real();
		double Ey = E_strength(Point(x, y, z))(1).real();
		field(i, j) =  Ex * Ex + Ey * Ey;// norm(real(E_strength(Point(x, y, z))));
	}
	save_csv(field, "field.txt");
}



//void special_output() {
//	int nx = 100, ny = 100;
//	double dx = 2.0 / nx, dy = 2.0 / ny;
//	double wave_len = c_ / omega_ * 2 * pi_;
//	MatrixXd field(nx, ny);
//	cout << "Save field distribution..." << endl;
//	#pragma omp parallel for
//	for (int i = 0; i < nx; i++)
//	for (int j = 0; j < ny; j++) {
//		double x = -1.0 + i * dx;
//		double y = -1.0 + j * dy;
//		double z = 15.2;
//		//Vector_c E = E_strength(Point(x, y, z));
//		//Vector_c H = H_strength(Point(x, y, z));
//		//Complex Ex = E(0);
//		//Complex Ey = E(1);
//		field(i, j) = W_strength(Point(x, y, z))(2);
//	}
//	save_csv(field, "field.txt");
//}

//void special_output() {
//	VectorXd Z(30001), A(30001);
//	Z = linspace(0, 15.20, 30000);
//	cout << "Save the field along Z axis..." << endl;
//#pragma omp parallel for
//	for (int i = 0; i < 30001; i++) {
//		double Ex = E_strength(Point(0, 0, Z(i)))(0).real();
//		double Ey = E_strength(Point(0, 0, Z(i)))(1).real();
//		A(i) = Ex;
//	}
//	save(Z, A, "field.txt");
//}

void save_current_distribution() {
	save(Im.real(), Im.imag(), "I.txt");
}

void print_Im() {
	for (int i = 0; i < Im.size(); i++) {
		printf("(%f, %f)\n", Im(i).real(), Im(i).imag());
	}
}

void load_current_distribution() {
	ifstream fin("I.txt");
	std::cout << "Load the current distribution instead of computing it..."
		<< std::endl;
	int N = edge_element.size();
	double real;
	double imag;
	Im.resize(N);
	for (int i = 0; i < N; i++) {
		fin >> real;
		fin >> imag;
		Im(i) = Complex(real, imag);
	}
	// print_Im();
}
