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

// ref-1045\t83\tref\t4255\t99\t148=\t=\t4006\t-397\tCCAGCGTCCCAGTGTGGGTGGCGTTTGCCTGGGCTGGGTACTGAGGCCGAGGTCCCCGCCACATCGTGGGCTCTGGGGTTAGGGCTGGGGAGGACAGCCTTGCCCCCGAGTGCGCTGACTGTCTTGGCCGTCTAGGGGGCATGTGGCC\tCGGCCG8C=GGGGGGGGGCGGGGCGGCGGGG1GGCGGGGG=CGCGGCGCCGGGCCGGCGGGGGGCG=GGGGGCGG(GCCJGGCGGJJJGJCJGJGGGGCGJ=CCGJJJJJJJGGJJGCJJGGJGJCGCJ=GJJJJGGGGGGGGG=C=C\nref-1045\t163\tref\t4006\t99\t148=\t=\t4255\t397\tTTGTAGTGGGTGCACACGCGTGCACTGGGACCCCACACAGCAATACGAGTCCAACTTAATAAACACATTTCTGGGGTTCCTCAGGCTGAGCATCTCTCTCTGGCATGTGGGGCAGCTGCGCAACTTTGGGTCCTGTCTGGGGTCCAGG\tCCCGGGGCGGGGGJGGJJJJGJ1JJJJJGJC==GJJCGJJGJGCGJGJGGJJJJGGGGGJJGJJGGGGGGCGGGJGGGGCGC8GGG=C8CCGG=GGGGC=GCJGCGGGCCG=GGGCGGGGCGGGGGCCCGGGGGG=GCGGCGGGGG=C\n

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
  WriteFastQRead(writer, read1, 1);
  WriteFastQRead(writer, read2, 2);
}

void FilterReads(const std::string& fasta_path, const std::string& sam_path,
                 const std::string& fastq_path,
                 const GCNormalizedCoverage& gc_map) {
  // Open the input SAM file and any output files
  sl::BamReader reader;
  reader.Open(sam_path);

  std::ofstream writer(fastq_path);

  // Load alleles from a FASTA file
  std::vector<std::string> contigs;
  {
    const auto& header = reader.Header();
    contigs.resize(header.NumSequences());

    sl::FastqReader contig_reader(fasta_path);
    sl::UnalignedSequence next_sequence;
    while (contig_reader.GetNextSequence(next_sequence)) {
      contigs[header.Name2ID(next_sequence.Name)] = next_sequence.Seq;
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
    auto start = std::min(read1.Position(), read2.Position());
    auto length = std::abs(read1.InsertSize());

    auto& sequence = contigs[read1.ChrID()];
    int gc = 0;
    for (int i = start, e = start + length; i < e; i++) {
      if (sequence[i] == 'G' || sequence[i] == 'C' || sequence[i] == 'g' ||
          sequence[i] == 'c')
        gc += 1;
    }
    int gc_fraction = std::lround(static_cast<float>(gc * 100) / length);

    // Downsample reads based on GC normalized coverage
    if (gc_map.find(gc_fraction) == gc_map.end())
      std::cerr << "GC fraction: " << gc_fraction << std::endl;
    
    float gc_norm_covg = gc_map.at(gc_fraction);
    if (dist(engine) < gc_norm_covg)
      WriteFastQ(writer, read1.FirstFlag() ? read1 : read2,
                 read1.FirstFlag() ? read2 : read1);
  }
}

}  // namespace npsv