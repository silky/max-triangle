#include <string>
#include <fstream>
#include "max-triangle.hpp"

bool test_max() {
    std::vector <std::vector <mpq_class> > polygons;

    unsigned int i;
    mpq_class ax,ay,bx,by,cx,cy;
    mpq_class z1,z2;

    std::string line, entry;
    std::ifstream infile("test-polygons");

    unsigned int ret[3];
    InscribedTriangle::state_flag status;

    while (std::getline(infile, line)) {
	std::stringstream line_stream(line);
	polygons.push_back(std::vector <mpq_class>());
	while (line_stream >> entry) {
	    if (!entry.empty()) polygons.back().push_back(mpq_class(entry));
	}
    }

    for (i = 0; i < polygons.size(); i++) {
	InscribedTriangle::maximum_triangle(polygons[i], ret, &status);
	if (status != InscribedTriangle::status_ok) {
	    std::cout << i << " ";
	    print_status(status);
	    return false;
	}
	ax = polygons[i][2*ret[0]+0];
	ay = polygons[i][2*ret[0]+1];
	bx = polygons[i][2*ret[1]+0];
	by = polygons[i][2*ret[1]+1];
	cx = polygons[i][2*ret[2]+0];
	cy = polygons[i][2*ret[2]+1];
	z1 = (bx-ax)*(cy-ay) - (cx-ax)*(by-ay);
	InscribedTriangle::brute_force_maximum_triangle(polygons[i], ret, &status);
    	ax = polygons[i][2*ret[0]+0];
	ay = polygons[i][2*ret[0]+1];
	bx = polygons[i][2*ret[1]+0];
	by = polygons[i][2*ret[1]+1];
	cx = polygons[i][2*ret[2]+0];
	cy = polygons[i][2*ret[2]+1];
	z2 = (bx-ax)*(cy-ay) - (cx-ax)*(by-ay);
	std::cout << z1 << " = " << z2 << std::endl;
	if (z1 != z2) return false;
    }
    return true;
}

int main(int argc, char *argv[]){
    if (test_max()) std::cout << "maximum triangle finding unit test PASSED" << std::endl;
    else std::cout << "maximum triangle finding unit test FAILED" << std::endl;
    return 0;
}

