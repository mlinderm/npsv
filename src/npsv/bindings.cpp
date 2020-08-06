#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "aligner.hpp"
#include "simulation.hpp"

namespace py = pybind11;

PYBIND11_MODULE(npsva, m) {
  m.doc() = "NPSV alignment tools";
  py::class_<npsv::Realigner>(m, "Realigner")
      .def(py::init<const std::string&, double, double, const npsv::InsertSizeDistribution::density_type&>())
      .def("count_alignments", &npsv::Realigner::CountAlignments);
  m.def("filter_reads_gc", &npsv::FilterReadsGC, "Filter reads based on GC normalized coverage");
  m.def("filter_reads_gnomad", &npsv::FilterReadsGnomAD, "Filter reads based on gnomAD coverage profile");
  m.def("test_score_alignment", &npsv::test::TestScoreAlignment, "Test interface for scoring alignment");
  m.def("test_alignment_overlap", &npsv::test::TestAlignmentOverlap, "Test interface for testing overlaps");
}