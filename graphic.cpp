
#include "graphic.h"
#include "dataio.h"

void polar(const VectorXd &X, const VectorXd &Y) {
	save(X, Y, "graph");
	system("python pplot.py");
}

void plot(const VectorXd &X, const VectorXd &Y) {
	save(X, Y, "graph");
	system("python plot.py");
}
