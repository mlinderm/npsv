#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "realigner.hpp"
#include "simulation.hpp"

namespace py = pybind11;

PYBIND11_MODULE(npsva, m) {
  m.doc() = "NPSV alignment tools";
   
  py::class_<npsv::RealignedFragments>(m, "RealignedFragments")
      .def(py::init<const std::string&, double, double, const npsv::InsertSizeDistribution::density_type&, const std::string&>())
      .def("size", &npsv::RealignedFragments::size)
      .def("gather_reads", &npsv::RealignedFragments::GatherReads, py::arg("region"), py::arg("max_reads")=std::numeric_limits<int>::max())
      .def("count_pipeline_straddlers", &npsv::RealignedFragments::CountPipelineStraddlers)
      .def("count_realigned_reads", &npsv::RealignedFragments::CountRealignedReads);

  m.def("filter_reads_gc", &npsv::FilterReadsGC, "Filter reads based on GC normalized coverage");
  m.def("filter_reads_gnomad", &npsv::FilterReadsGnomAD, "Filter reads based on gnomAD coverage profile");
  
  m.def("test_score_alignment", &npsv::test::TestScoreAlignment, "Test interface for scoring alignment");
  m.def("test_alignment_overlap", &npsv::test::TestAlignmentOverlap, "Test interface for testing overlaps");
  m.def("test_straddle", &npsv::test::TestStraddle, "Test interface for testing if read pair straddles regions");
}