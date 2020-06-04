#include "aligner.hpp"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/FastqReader.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/RefGenome.h"

namespace npsv {
namespace {
bool EndsWith(const std::string& str, const std::string& suffix) {
  // https://stackoverflow.com/questions/874134/find-out-if-string-ends-with-another-string-in-c
  return str.size() >= suffix.size() &&
         str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

}  // namespace

AlleleReference::AlleleReference(const std::string& fasta_path) {
  // Load alleles as a FASTA
  sl::FastqReader contigs(fasta_path);
  sl::UnalignedSequenceVector sequences;
  {  // Load contigs
    sl::UnalignedSequence next_sequence;
    while (contigs.GetNextSequence(next_sequence)) {
      sequences.push_back(next_sequence);
    }
  }

  // Construct index and tag alternate contigs
  bwa_.ConstructIndex(sequences);
  {
    bwaidx_t* index = bwa_.GetIndex();
    for (int i = 0; i < index->bns->n_seqs; i++) {
      if (EndsWith(index->bns->anns[i].name, "_alt")) {
        index->bns->anns[i].is_alt = 1;
      }
    }
  }
}

py::dict AlleleReference::CountAlignments(const std::string& bam_path,
                                      const std::string& rl_region,
                                      const std::string& al_region,
                                      py::kwargs kwargs) const {
  // Construct query intervals
  sl::GenomicRegion rl(rl_region, bwa_.HeaderFromIndex());  // Reference "left"
  sl::GenomicRegion al(al_region, bwa_.HeaderFromIndex());  // Alternate "left"

  // The "right" breakpoints are optional
  sl::GenomicRegion rr = sl::GenomicRegion();  // Alternate "right"
  sl::GenomicRegion ar = sl::GenomicRegion();  // Alternate "right"
  if (kwargs && kwargs.contains("rr_region")) {
    rr = sl::GenomicRegion(py::cast<std::string>(kwargs["rr_region"]),
                           bwa_.HeaderFromIndex());
  }
  if (kwargs && kwargs.contains("ar_region")) {
    ar = sl::GenomicRegion(py::cast<std::string>(kwargs["rr_region"]),
                           bwa_.HeaderFromIndex());
  }

  // Open the input BAM/SAM/CRAM
  sl::BamReader reader;
  reader.Open(bam_path);

  // Subset the reader by the specified regions
  if (kwargs && kwargs.contains("region")) {
    // Restrict input BAM to specific regions
    sl::GRC reader_regions;
    for (auto& region : py::cast<py::list>(kwargs["region"])) {
      reader_regions.add(
          sl::GenomicRegion(py::cast<std::string>(region), reader.Header()));
    }
    reader.SetMultipleRegions(reader_regions);
  }

  int rl_reads = 0, al_reads = 0, rr_reads = 0, ar_reads = 0;

  sl::BamRecord read;
  while (reader.GetNextRecord(read)) {
    int32_t rl_score = 0, rr_score = 0, al_score = 0, ar_score = 0;

    // hardclip=false, secondary score cutoff=0.9, max secondary alignments=10
    sl::BamRecordVector results;
    bwa_.AlignSequence(read.Sequence(), read.Qname(), results, false, 0.9, 10);
    for (auto& result : results) {
      if (result.AlignmentFlag() &
          (SAMFlags::SecondaryAlignment | SAMFlags::SupplementaryAlignment)) {
        continue;
      }

      auto read_region = result.AsGenomicRegion();
      if (read_region.GetOverlap(rl) == GenomicRegionOverlap::ContainsArg) {
        result.GetIntTag("AS", rl_score);
      }
      if (read_region.GetOverlap(al) == GenomicRegionOverlap::ContainsArg) {
        result.GetIntTag("AS", al_score);
      }
      if (!rr.IsEmpty() &&
          read_region.GetOverlap(rr) == GenomicRegionOverlap::ContainsArg) {
        result.GetIntTag("AS", rr_score);
      }
      if (!ar.IsEmpty() &&
          read_region.GetOverlap(ar) == GenomicRegionOverlap::ContainsArg) {
        result.GetIntTag("AS", ar_score);
      }
    }

    // Record counts of reads spanning the breakpoints
    if (al_score > 0 && al_score > rl_score) {
      al_reads += 1;
    }
    if (rl_score > 0 && rl_score > al_score) {
      rl_reads += 1;
    }
    if (ar_score > 0 && ar_score > rr_score) {
      ar_reads += 1;
    }
    if (rr_score > 0 && rr_score > ar_score) {
      rr_reads += 1;
    }
  }

  auto results = py::dict();
  results["rl_reads"] = rl_reads;
  results["rr_reads"] = rr_reads;
  results["al_reads"] = al_reads;
  results["ar_reads"] = ar_reads;
  return results;
}

}  // namespace npsv