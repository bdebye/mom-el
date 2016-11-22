
#include "rwg_edge.h"
#include "mathlib.h"
#include "gmsh.h"


RWGedge::RWGedge(int index,
	const EdgeIndex & edge,
	int vertex1,
	int vertex2,
	int trian_index1,
	int trian_index2
	) :
	index(index),
	eg_share(edge),
	ve_plus(vertex1),
	ve_minus(vertex2),
	t_plus_index(trian_index1),
	t_minus_index(trian_index2) {

}

Point RWGedge::vertex_plus() const {
	return mesh.node[ve_plus];
}

int RWGedge::vertex_plus_index() const {
	return ve_plus;
}

Point RWGedge::vertex_minus() const {
	return mesh.node[ve_minus];
}

int RWGedge::vertex_minus_index() const {
	return ve_minus;
}

double RWGedge::plus_area() const {
	return triangle_area(t_plus_index);
}

double RWGedge::minus_area() const {
	return triangle_area(t_minus_index);
}

Point RWGedge::plus_center() const {
	return trian_center(get_triangle(t_plus_index));
}

Point RWGedge::minus_center() const {
	return trian_center(get_triangle(t_minus_index));
}

EdgeIndex RWGedge::edge_share() const {
	return eg_share;
}

int RWGedge::trian_plus_index() const {
	return t_plus_index;
}

int RWGedge::trian_minus_index() const {
	return t_minus_index;
}

double RWGedge::share_edge_length() const {
	return distance(get_point(eg_share(0)), get_point(eg_share(1)));
}

const Vector RWGedge::rho_c_plus() const {
	return Vector(plus_center() - vertex_plus());
}

const Vector RWGedge::rho_c_minus() const {
	return Vector(vertex_minus() - minus_center());
}

const Vector RWGedge::rho_plus(const Point & point) const {
	return Vector(point - vertex_plus());
}

const Vector RWGedge::rho_minus(const Point & point) const {
	return Vector(vertex_minus() - point);
}


Triangle RWGedge::trian_plus() const {
	return get_triangle(this->trian_plus_index());
}

Triangle RWGedge::trian_minus() const {
	return get_triangle(this->trian_minus_index());
}
