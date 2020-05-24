#include <pybind11/pybind11.h>
#include "aligner.hpp"

namespace py = pybind11;

PYBIND11_MODULE(npsva, m) {
  m.doc() = "NPSV alignment tools";
  m.def("add", &add);
}