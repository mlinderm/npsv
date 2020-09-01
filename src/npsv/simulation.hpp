#pragma once
#include <pybind11/pybind11.h>
#include <map>
#include "SeqLib/UnalignedSequence.h"

namespace py = pybind11;
namespace sl = SeqLib;

namespace npsv {

void FilterReadsGC(const std::string& fasta_path, const std::string& sam_path, const std::string& fastq_path, const std::vector<float>& gc_covg);

void FilterReadsGnomAD(const std::string& covg_path, const std::string& sam_path, const std::string& fastq_path, float max_covg);

}