#pragma once
#include <pybind11/pybind11.h>
#include <map>
#include "SeqLib/UnalignedSequence.h"

namespace py = pybind11;
namespace sl = SeqLib;

namespace npsv {

typedef std::map<int, float> GCNormalizedCoverage;

void FilterReads(const std::string& fasta_path, const std::string& sam_path, const std::string& fastq_path, const GCNormalizedCoverage& gc_map);


}