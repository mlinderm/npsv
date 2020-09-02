#include "simulation.hpp"

#include <fstream>
#include <random>
#include <stdexcept>
#include <vector>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/FastqReader.h"
#include "SeqLib/SeqLibUtils.h"
#include "SeqLib/UnalignedSequence.h"

namespace {
void assert_throw(const bool cond, const std::string& text,
                  const std::string& file, const int line) {
  if (!cond) {
    throw std::runtime_error(text + ". In file: " + file +
                             " on line: " + std::to_string(line));
  }
}

#define pyassert(cond, text) assert_throw(cond, text, __FILE__, __LINE__)
}  // namespace

namespace npsv {

void WriteFastQRead(std::ofstream& writer, const sl::BamRecord& read, int num) {
  writer << "@" << read.Qname() << "/" << num << std::endl;
  if (read.ReverseFlag()) {
    std::string sequence(read.Sequence()), qualities(read.Qualities());
    sl::rcomplement(sequence);
    std::reverse(qualities.begin(), qualities.end());
    writer << sequence << std::endl
           << "+" << std::endl
           << qualities << std::endl;
  } else {
    writer << read.Sequence() << std::endl
           << "+" << std::endl
           << read.Qualities() << std::endl;
  }
}

void WriteFastQ(std::ofstream& writer, const sl::BamRecord& read1,
                const sl::BamRecord& read2) {
  if (read1.FirstFlag()) {
    WriteFastQRead(writer, read1, 1);
    WriteFastQRead(writer, read2, 2);
  } else {
    WriteFastQRead(writer, read2, 1);
    WriteFastQRead(writer, read1, 2);
  }
}

constexpr bool IsGC(char base) {
  return base == 'G' || base == 'C' || base == 'g' || base == 'c';
}

void FilterReadsGC(const std::string& fasta_path, const std::string& sam_path,
                   const std::string& fastq_path,
                   const GCNormalizedCoverage& gc_map) {
  // Open the input SAM file and any output files
  sl::BamReader reader;
  reader.Open(sam_path);

  std::ofstream writer(fastq_path);

  // Load alleles from a FASTA file
  std::vector<std::vector<int> > contigs;
  {
    const auto& header = reader.Header();
    contigs.resize(header.NumSequences());

    sl::FastqReader contig_reader(fasta_path);
    sl::UnalignedSequence next_sequence;
    while (contig_reader.GetNextSequence(next_sequence)) {
      auto& sequence = next_sequence.Seq;
      auto& cuml_gc_count = contigs[header.Name2ID(next_sequence.Name)];

      // Precompute the GC count in each sequence forming vector of cumulative
      // gc counts that can be used to calculate GC fraction for fragments
      cuml_gc_count.resize(sequence.size() + 1);
      cuml_gc_count[0] = 0;
      for (int i = 0; i < sequence.size(); i++) {
        if (IsGC(sequence[i]))
          cuml_gc_count[i + 1] = cuml_gc_count[i] + 1;
        else
          cuml_gc_count[i + 1] = cuml_gc_count[i];
      }
    }
  }

  // Setup random number generator
  std::default_random_engine engine;
  std::uniform_real_distribution<> dist(0.0, 1.0);

  sl::BamRecord read1, read2;
  while (true) {
    if (!reader.GetNextRecord(read1)) break;
    pyassert(reader.GetNextRecord(read2), "Missing second read in pair");

    // Compute GC fraction for insert
    auto& cuml_gc_count = contigs[read1.ChrID()];
    
    // Due to insertions it is possible for the fragment to extend beyond the legnth of the sequence, so
    // clamp at the sequence length
    auto start = std::min(read1.Position(), read2.Position());
    pyassert(start < cuml_gc_count.size(), "Fragment start coordinate outside of sequence");
    
    auto length = std::abs(read1.InsertSize());
    if (start + length >= cuml_gc_count.size())
      length = cuml_gc_count.size() - start - 1;

    int gc = cuml_gc_count[start + length] - cuml_gc_count[start];
    int gc_fraction = std::lround(static_cast<float>(gc * 100) / length);
    pyassert(gc_fraction >= 0 && gc_fraction <= 100, "GC fraction outside of expected range");

    // Downsample reads based on GC normalized coverage
    auto iter = gc_map.find(gc_fraction);
    pyassert(iter != gc_map.end(), "No data for GC fraction");
    
    float gc_norm_covg = iter->second;
    if (dist(engine) < gc_norm_covg)
      WriteFastQ(writer, read1, read2);
  }
}

void FilterReadsGnomAD(const std::string& covg_path,
                       const std::string& sam_path,
                       const std::string& fastq_path, float max_covg) {
  // Open the input SAM file and any output files
  sl::BamReader reader;
  reader.Open(sam_path);

  std::ofstream writer(fastq_path);

  // Load coverage profile
  std::vector<std::vector<int> > contigs;
  {
    const auto& header = reader.Header();
    contigs.resize(header.NumSequences());

    std::ifstream coverage_reader(covg_path);
    while (coverage_reader) {
      std::string line;
      std::getline(coverage_reader, line);
      if (line.empty())
        break;
      auto sep = line.find_first_of('\t');  // Find seperator between name and coverage
      pyassert(sep != std::string::npos, "Couldn't find coverage profile separator");

      auto& cuml_coverage = contigs[header.Name2ID(line.substr(0, sep))];
      sep += 1; // Point to start of the coverage data
      cuml_coverage.resize(line.size() - sep + 1);
      cuml_coverage[0] = 0;
      for (int i = 0; i < (line.size() - sep); i++) {
        cuml_coverage[i+1] = cuml_coverage[i] + static_cast<int>(line[i+sep]) - 33;
      }
    }
  }

  // Setup random number generator
  std::default_random_engine engine;
  std::uniform_real_distribution<> dist(0.0, 1.0);

  sl::BamRecord read1, read2;
  while (true) {
    if (!reader.GetNextRecord(read1)) break;
    pyassert(reader.GetNextRecord(read2), "Missing second read in pair");

    // Compute mean coverage for reads
    auto& cuml_coverage = contigs[read1.ChrID()];
    int total_coverage = (cuml_coverage[read1.PositionEnd()] - cuml_coverage[read1.Position()]) + (cuml_coverage[read2.PositionEnd()] - cuml_coverage[read2.Position()]);
    int total_length = (read1.PositionEnd() - read1.Position()) + (read2.PositionEnd() - read2.Position());
    std::cerr << "Length: " << total_length << std::endl;
    pyassert(false, "Bail");
    float gnomad_norm_covg = static_cast<float>(total_coverage) / total_length / max_covg;

    if (dist(engine) < gnomad_norm_covg)
      WriteFastQ(writer, read1, read2);
  }
}

}  // namespace npsv