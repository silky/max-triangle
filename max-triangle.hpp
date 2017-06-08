#include <gmpxx.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>

namespace InscribedTriangle {
    enum state_flag {status_ok, status_not_convex, status_no_interior, status_runtime_error, status_max_not_at_vertex, status_maxiter_exceeded};
    void print_status(state_flag status);
    bool is_convex(std::vector<mpq_class> polygon);
    void maximum_triangle(std::vector<mpq_class> polygon, unsigned int *ret, state_flag *status);
    void brute_force_maximum_triangle(std::vector<mpq_class> polygon, unsigned int *ret, state_flag *status);
}
