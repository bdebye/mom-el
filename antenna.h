

#ifndef ANTENNA_H_
#define ANTENNA_H_

#include "type.h"
#include "rwg_edge.h"
#include "gmsh.h"
#include "dataio.h"
#include "mathlib.h"
#include "graphic.h"
#include "feed.h"

#include <map>

extern double omega;

void generate_edge_elements();

void set_radiate_model(double omega);

Complex get_input_impedance(int feed_edge);

Complex get_input_power(int feed_edge);

const RWGedge &get_RWGedge(int n);

int get_RWGedge_num();

const MatrixXcd &get_impedance_matrix();

double get_omega();

double near_field_radius();

double U_density(const Point &r);

Vector_c E_strength(const Point &r);

Vector_c H_strength(const Point &r);

Vector_c Poynting(const Point &r);

//Vector P_strength(const Point &r);

double phase(const Point &r);

void save_impedance_matrix();

bool same_edge_index(Index2 idx1, Index2 idx2);

void print_shared_edge();

void directional_radiate_pattern_XZ();

void solve_moment_equation();

void directional_radiate_pattern_XY();

void print_feed_information();

void special_output();

void save_current_distribution();

void load_current_distribution();

#endif /* ANTENNA_H_ */
