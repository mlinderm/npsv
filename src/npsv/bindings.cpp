#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "aligner.hpp"

namespace py = pybind11;

PYBIND11_MODULE(npsva, m) {
  m.doc() = "NPSV alignment tools";
  py::class_<npsv::Realigner>(m, "Realigner")
      .def(py::init<const std::string&, double, double, const npsv::InsertSizeDistribution::density_type&>())
      .def("count_alignments", &npsv::Realigner::CountAlignments);
  m.def("alternate_variant_locations", &npsv::AlternateVariantLocations, "Find alternate possible alignments for allele");
  m.def("test_score_alignment", &npsv::test::TestScoreAlignment, "Test interface for scoring alignment");
}