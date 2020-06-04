#include <pybind11/pybind11.h>
#include "aligner.hpp"

namespace py = pybind11;

PYBIND11_MODULE(npsva, m) {
  m.doc() = "NPSV alignment tools";
  py::class_<npsv::AlleleReference>(m, "AlleleReference")
      .def(py::init<const std::string&, double, double>())
      .def("count_alignments", &npsv::AlleleReference::CountAlignments);
}