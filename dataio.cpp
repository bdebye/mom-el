
#include "dataio.h"

#include <fstream>
#include <iostream>

using namespace std;

bool detect_file(const string filename, int& m, int& n) {
	ifstream fin(filename);
	string line;
	int count = 1, wcount = 0;
	double value;

	if (fin.bad()) {
		cout << "cannot open file "
			<< filename << " ..." << endl;
		return false;
	}

	getline(fin, line);
	wcount = 0;
	stringstream lins(line);
	while (lins) {
		lins >> value;
		wcount++;
	}
	n = wcount - 1;
	while (getline(fin, line)) {
		if (line.length() > 0) {
			count++;
		}
	}
	m = count;
	fin.close();
	return true;
}


MatrixXd load(const string filename) {
	int m, n;
	detect_file(filename, m, n);
	ifstream fin(filename);
	MatrixXd data(m, n);
	for (int i = 0; i < m; i++)
	for (int j = 0; j < n; j++) {
		fin >> data(i, j);
	}
	fin.close();
	return data;
}


void save_csv(const MatrixXd & data, string filename) {
	int m = data.rows();
	int n = data.cols();
	ofstream fout(filename);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			fout << data(i, j);
			if (j != n - 1)
				fout << ",";
		}
		fout << endl;
	}
	fout.close();
}

void save_csv(const VectorXd &vect_base,
	const VectorXd vect_value,
	const string filename = "temp") {

	int leng = vect_base.size();
	if (leng != vect_value.size()) {
		cout << "Error, lengths of two vectors are not same..." << endl;
		return;
	}
	ofstream fout(filename);
	for (int i = 0; i < leng; i++) {
		fout << vect_base(i) << ",\t"
			<< vect_value(i)
			<< endl;
	}
	fout.close();
}

void save(const VectorXd &vect_base,
	const VectorXd &vect_value,
	const string filename = "temp") {

	int leng = vect_base.size();
	if (leng != vect_value.size()) {
		cout << "Error, lengths of two vectors are not same..." << endl;
		return;
	}
	FILE *fp = fopen(filename.c_str(), "w");
	for (int i = 0; i < leng; i++) {
		fprintf(fp, "%12e\t %12e\n", vect_base(i), vect_value(i));
	}
	fclose(fp);
}

void save(const VectorXd & data, const string filename) {
	int n = data.size();
	ofstream fout(filename);
	for (int i = 0; i < n; i++) {
		fout << data(i) << endl;
	}
	fout.close();
}

const string POS_HEADER = "View \"field\" {";
const string POS_FOOTER = "};";

const char SCALAR_POINT_FORMAT[] = "SP(%f, %f, %f){%e};";
const char VECTOR_POINT_FORMAT[] = "VP(%f, %f, %f){%e, %e, %e};";
char BUFFER[500];

void pos_scalar_field(string filename, const vector<scalar_field> &dist) {
	ofstream ofs(filename);
	ofs << POS_HEADER << endl;
	for (auto iter = dist.begin(); iter != dist.end(); ++iter) {

#ifdef _WIN32
	sprintf_s(BUFFER, SCALAR_POINT_FORMAT, iter->posi(0), iter->posi(1), iter->posi(2),
			iter->value);
#else
	sprintf(BUFFER, SCALAR_POINT_FORMAT, iter->posi(0), iter->posi(1), iter->posi(2),
			iter->value);
#endif
		ofs << BUFFER << endl;
	}
	ofs << POS_FOOTER << endl;
}

void pos_vector_field(string filename, const vector<vector_field> &dist) {
	ofstream ofs(filename);
	ofs << POS_HEADER << endl;
	for (auto iter = dist.begin(); iter != dist.end(); ++iter) {
#ifdef _WIN32
		sprintf_s(BUFFER, VECTOR_POINT_FORMAT, iter->posi(0), iter->posi(1), iter->posi(2),
			iter->value(0), iter->value(1), iter->value(2));
#else
		sprintf(BUFFER, VECTOR_POINT_FORMAT, iter->posi(0), iter->posi(1), iter->posi(2),
			iter->value(0), iter->value(1), iter->value(2));
#endif
		ofs << BUFFER << endl;
	}
	ofs << POS_FOOTER << endl;
}
