#include <gmpxx.h>
#include <vector>

namespace InscribedTriangle {
    enum state_flag {status_ok, status_not_convex, status_no_interior, status_runtime_error, status_max_not_at_vertex, status_maxiter_exceeded};
    bool is_convex(std::vector<mpq_class> polygon);
    void anchored_triangle(std::vector<mpq_class> polygon, mpq_class nx, mpq_class ny, unsigned int *reti, mpq_class *rett, state_flag *status);
    void maximum_triangle(std::vector<mpq_class> polygon, unsigned int *ret, state_flag *status);
    void brute_force_maximum_triangle(std::vector<mpq_class> polygon, unsigned int *ret, state_flag *status);
}
