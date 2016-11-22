

#ifndef TYPE_H_
#define TYPE_H_


#include <complex>

#define EIGEN_USE_MKL_ALL

#include <Eigen/Dense>
#include <array>
using namespace Eigen;

typedef Eigen::Vector3d Point;
typedef Eigen::Vector3d Vector;
typedef Eigen::Vector3cd Vector_c;
typedef Eigen::Vector3i TrianIndex;
typedef Eigen::Vector2i EdgeIndex;
typedef std::array<Point, 3> Triangle;

typedef std::complex<double> Complex;

typedef Eigen::Vector3i Index3;
typedef Eigen::Vector2i Index2;

//typedef mxArray *matArray;

extern const double epsilon_;
extern const double mu_;
extern const double c_;	// speed of light
extern const double eta_;  // free-space impedance
//extern const double epsilon_R ;
extern const double pi_;

extern const Complex j;

extern const std::string GRAPH_DATA_FILE;



#endif /* TYPE_H_ */
