/*
* dataio.h
*
*  Created on: Sep 28, 2015
*      Author: lambda
*/

#ifndef DATAIO_H_
#define DATAIO_H_


#include "type.h"
#include <vector>
#include <string>


bool detect_file(const std::string filename, int& m, int& n);


MatrixXd load(const std::string filename);

//template <typename T>
//void load(const string filename, Array<T, 2> & matrix);

//template <typename T>
//void save(const Array<T, 2> & data, string filename);


void save_csv(const MatrixXd & data, std::string filename);

void save_csv(const VectorXd &vect_base,
	const VectorXd &vect_value,
	const std::string filename);

void save(const VectorXd &vect_base,
	const VectorXd &vect_value,
	const std::string filename);

void save(const VectorXd& data, const std::string filename);


struct scalar_field {
	Vector posi;
	double value;
};

struct vector_field {
	Vector posi;
	Vector value;
};

void pos_scalar_field(std::string filename, const std::vector<scalar_field> &dist);
void pos_vector_field(std::string filename, const std::vector<vector_field> &dist);

#endif /* DATAIO_H_ */
