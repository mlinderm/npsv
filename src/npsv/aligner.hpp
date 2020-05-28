#pragma once
#include <pybind11/pybind11.h>
#include "SeqLib/BWAWrapper.h"

namespace py = pybind11;
namespace sl = SeqLib;

int add(int i, int j);

namespace npsv {

enum GenomicRegionOverlap {
  NoOverlap = 0,
  PartialOverlap = 1,
  ContainsArg = 2,
  ContainedInArg = 3
};

enum SAMFlags { SecondaryAlignment = 256, SupplementaryAlignment = 2048 };

class AlleleReference {
 public:
  AlleleReference(const std::string& fasta_path);

  py::dict CountAlignments(const std::string& bam_path,
                           const std::string& rl_region,
                           const std::string& al_region,
                           py::kwargs kwargs) const;

 private:
  sl::BWAWrapper bwa_;
};

}  // namespace npsv