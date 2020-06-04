#include "aligner.hpp"

#include <iostream>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/FastqReader.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/RefGenome.h"

void assert_throw(const bool cond, const std::string& text,
                  const std::string& file, const int line) {
  if (!cond) {
    throw std::runtime_error(text + ". In file: " + file +
                             " on line: " + std::to_string(line));
  }
}

#define pyassert(cond, text) assert_throw(cond, text, __FILE__, __LINE__)

namespace npsv {
namespace {
bool EndsWith(const std::string& str, const std::string& suffix) {
  // https://stackoverflow.com/questions/874134/find-out-if-string-ends-with-another-string-in-c
  return str.size() >= suffix.size() &&
         str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

}  // namespace

AlignmentPair::AlignmentPair(const sl::BamRecord& first,
                             const sl::BamRecord& second,
                             const InsertSizeDistribution& insert_size_dist)
    : first_(&first), second_(&second), score_(0.) {
  if (second_->Position() <= first_->Position()) std::swap(first_, second_);

  // Scoring algorithm adapted from svviz2:
  // https://github.com/nspies/svviz2/blob/44f7bfc75bf84c1db4563d9fd30bf20967d1c825/src/svviz2/io/readstatistics.py
  score_ = first_->MapQuality() + second_->MapQuality();
  if (!Concordant()) {
    score_ -= 10;
    return;
  }
  auto insert_size_prob = insert_size_dist(InsertSize());
  if (insert_size_prob == 0.) {
    score_ -= 10;
    return;
  }

  score_ += log10(insert_size_prob);
}

bool AlignmentPair::Concordant() const {
  if (!Valid()) return false;
  if (first_->ChrID() != second_->ChrID()) return false;

  // TODO: Check orientation

  return true;
}

int32_t AlignmentPair::InsertSize() const {
  pyassert(first_->Position() <= second_->Position(),
           "First and second in sorted order");
  return second_->Position() + second_->Length() - first_->Position();
}

std::ostream& operator<<(std::ostream& os, const AlignmentPair& pair) {
  return (os << *(pair.first_) << *(pair.second_) << pair.score_ << std::endl);
}

void Fragment::SetSecond(const sl::BamRecord& read) {
  pyassert(!first_.isEmpty(),
           "Can't add second read to fragment if first not defined");
  pyassert(first_.PairedFlag() && read.PairedFlag(),
           "To create fragment reads must be paired");
  second_ = read;
}

void Fragment::SetBestPair(const InsertSizeDistribution& insert_size_dist) {
  std::priority_queue<AlignmentPair> alignment_pairs;
  for (auto& align1 : first_alignments_) {
    for (auto& align2 : second_alignments_) {
      alignment_pairs.emplace(align1, align2, insert_size_dist);
    }
  }
  if (!alignment_pairs.empty()) {
    best_pair_ = alignment_pairs.top();
  }
}

bool AlignmentPair::ReadsContain(const sl::GenomicRegion& region) const {
  if (!Valid()) return false;

  // TODO: Add minimum overlap?
  auto first_region = first_->AsGenomicRegion();
  if (first_region.GetOverlap(region) == GenomicRegionOverlap::ContainsArg)
    return true;
  auto second_region = second_->AsGenomicRegion();
  if (second_region.GetOverlap(region) == GenomicRegionOverlap::ContainsArg)
    return true;

  return false;
}

void AlleleAlignments::Initialize(const sl::UnalignedSequence& sequence) {
  pyassert(!IsInitialized(), "BWA should not previously have been initialized");
  bwa_.ConstructIndex({sequence});
}

void AlleleAlignments::Align(const sl::BamRecord& read,
                             const InsertSizeDistribution& insert_size_dist) {
  auto emplace_result = fragments_.emplace(std::piecewise_construct,
                                           std::forward_as_tuple(read.Qname()),
                                           std::forward_as_tuple(read));
  auto& fragment = emplace_result.first->second;
  if (emplace_result.second) {
    // A new fragment
    bwa_.AlignSequence(read.Sequence(), read.Qname(),
                       fragment.first_alignments_, false, 0.9, 10);
  } else {
    // An existing fragment
    fragment.SetSecond(read);
    bwa_.AlignSequence(read.Sequence(), read.Qname(),
                       fragment.second_alignments_, false, 0.9, 10);
    fragment.SetBestPair(insert_size_dist);
  }
}

AlleleReference::AlleleReference(const std::string& fasta_path,
                                 double insert_size_mean,
                                 double insert_size_std)
    : insert_size_dist_(insert_size_mean, insert_size_std) {
  pyassert(!ref_.IsInitialized() && !alt_.IsInitialized(),
           "BWA wrapper should start out uninitialized");

  // Load alleles as a FASTA
  sl::FastqReader contigs(fasta_path);
  {
    sl::UnalignedSequence next_sequence;
    while (contigs.GetNextSequence(next_sequence)) {
      if (EndsWith(next_sequence.Name, "_alt")) {
        alt_.Initialize(next_sequence);
      } else {
        ref_.Initialize(next_sequence);
      }
    }
  }
}

void AlleleReference::Clear() { 
  ref_.Clear();
  alt_.Clear();
}

void AlleleReference::Realign(const sl::BamRecord& read) {
  ref_.Align(read, insert_size_dist_);
  alt_.Align(read, insert_size_dist_);
}

py::dict AlleleReference::CountAlignments(const std::string& bam_path,
                                          const std::string& rl_region,
                                          const std::string& al_region,
                                          py::kwargs kwargs) {
  // Clear any previous aligned fragments
  Clear();

  // Construct query intervals
  sl::GenomicRegion rl(rl_region, RefHeader());  // Reference "left"
  sl::GenomicRegion al(al_region, AltHeader());  // Alternate "left"

  // The "right" breakpoints are optional
  sl::GenomicRegion rr = sl::GenomicRegion();  // Reference "right"
  sl::GenomicRegion ar = sl::GenomicRegion();  // Alternate "right"
  if (kwargs && kwargs.contains("rr_region")) {
    rr = sl::GenomicRegion(py::cast<std::string>(kwargs["rr_region"]),
                           RefHeader());
  }
  if (kwargs && kwargs.contains("ar_region")) {
    ar = sl::GenomicRegion(py::cast<std::string>(kwargs["rr_region"]),
                           AltHeader());
  }

  int rl_reads = 0, al_reads = 0, rr_reads = 0, ar_reads = 0;

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

  sl::BamRecord read;
  while (reader.GetNextRecord(read)) {
    if (read.DuplicateFlag() || read.SecondaryFlag() || read.SupplementaryFlag())
      continue; // Skip duplicate or secondary/supplemental alignments
    
    Realign(read);
  }
  pyassert(ref_.Size() == alt_.Size(),
           "All reads should be aligned to both alleles");

  auto ref_iter = ref_.begin(), ref_end = ref_.end();
  auto alt_iter = alt_.begin(), alt_end = alt_.end();
  while (ref_iter != ref_end && alt_iter != alt_end) {
    pyassert(ref_iter->first == alt_iter->first,
             "Considering alignments for same underlying fragment");
    auto& ref_fragment = ref_iter->second;
    auto& alt_fragment = alt_iter->second;

    Fragment::score_type rl_score = 0, rr_score = 0, al_score = 0, ar_score = 0;

    if (ref_fragment.BestAlignmentContains(rl)) {
      rl_score = ref_fragment.BestAlignmentScore();
    }
    if (alt_fragment.BestAlignmentContains(al)) {
      al_score = alt_fragment.BestAlignmentScore();
    }
    if (!rr.IsEmpty() && ref_fragment.BestAlignmentContains(rr)) {
      rr_score = ref_fragment.BestAlignmentScore();
    } else if (rr.IsEmpty()) {
      rr_score = rl_score;
    }
    if (!ar.IsEmpty() && alt_fragment.BestAlignmentContains(ar)) {
      ar_score = alt_fragment.BestAlignmentScore();
    } else if (ar.IsEmpty()) {
      ar_score = al_score;
    }

    // Record counts of reads spanning the breakpoints (using svviz2 criteria
    // for assigning reads)
    if (al_score > 0 && al_score - rl_score > 1) {
      al_reads += 1;
    }
    if (rl_score > 0 && rl_score - al_score > 1) {
      rl_reads += 1;
    }
    // TODO: Make this work for different variant types, not just deletions
    if (!ar.IsEmpty() && ar_score > 0 && ar_score > rr_score) {
      ar_reads += 1;
    }
    if (!rr.IsEmpty() && rr_score > 0 && rr_score - ar_score > 1) {
      rr_reads += 1;
    }

    ref_iter++;
    alt_iter++;
  }

  auto results = py::dict();
  results["rl_reads"] = rl_reads;
  results["rr_reads"] = rr_reads;
  results["al_reads"] = al_reads;
  results["ar_reads"] = ar_reads;
  return results;
}

// py::dict AlleleReference::CountAlignments(const std::string& bam_path,
//                                       const std::string& rl_region,
//                                       const std::string& al_region,
//                                       py::kwargs kwargs) const {
//   // Construct query intervals
//   sl::GenomicRegion rl(rl_region, bwa_.HeaderFromIndex());  // Reference
//   "left" sl::GenomicRegion al(al_region, bwa_.HeaderFromIndex());  //
//   Alternate "left"

//   // The "right" breakpoints are optional
//   sl::GenomicRegion rr = sl::GenomicRegion();  // Alternate "right"
//   sl::GenomicRegion ar = sl::GenomicRegion();  // Alternate "right"
//   if (kwargs && kwargs.contains("rr_region")) {
//     rr = sl::GenomicRegion(py::cast<std::string>(kwargs["rr_region"]),
//                            bwa_.HeaderFromIndex());
//   }
//   if (kwargs && kwargs.contains("ar_region")) {
//     ar = sl::GenomicRegion(py::cast<std::string>(kwargs["rr_region"]),
//                            bwa_.HeaderFromIndex());
//   }

//   // Open the input BAM/SAM/CRAM
//   sl::BamReader reader;
//   reader.Open(bam_path);

//   // Subset the reader by the specified regions
//   if (kwargs && kwargs.contains("regions")) {
//     // Restrict input BAM to specific regions
//     sl::GRC reader_regions;
//     for (auto& region : py::cast<py::list>(kwargs["regions"])) {
//       reader_regions.add(
//           sl::GenomicRegion(py::cast<std::string>(region), reader.Header()));
//     }
//     reader.SetMultipleRegions(reader_regions);
//   }

//   int rl_reads = 0, al_reads = 0, rr_reads = 0, ar_reads = 0;

//   sl::BamWriter writer(sl::SAM);
//   writer.SetHeader(bwa_.HeaderFromIndex());
//   writer.Open("/dev/stdout");
//   writer.WriteHeader();

//   sl::BamRecord read;
//   while (reader.GetNextRecord(read)) {
//     //std::cerr << read << std::endl;
//     int32_t rl_score = 0, rr_score = 0, al_score = 0, ar_score = 0;

//     // hardclip=false, secondary score cutoff=0.9, max secondary
//     alignments=10 sl::BamRecordVector results;
//     bwa_.AlignSequence(read.Sequence(), read.Qname(), results, false, 0, 10);
//     for (auto& result : results) {
//       if (read.Qname() == "alt-282") {
//         writer.WriteRecord(result);
//         std::cout << result.AlignmentFlag() << std::endl;
//       }

//       if (result.AlignmentFlag() &
//           (SAMFlags::SecondaryAlignment | SAMFlags::SupplementaryAlignment))
//           {
//         continue;
//       }

//       auto read_region = result.AsGenomicRegion();
//       if (read_region.GetOverlap(rl) == GenomicRegionOverlap::ContainsArg) {
//         result.GetIntTag("AS", rl_score);
//         //std::cerr << result << std::endl;
//       }
//       if (read_region.GetOverlap(al) == GenomicRegionOverlap::ContainsArg) {
//         result.GetIntTag("AS", al_score);
//         //std::cerr << result << std::endl;
//       }
//       if (!rr.IsEmpty() &&
//           read_region.GetOverlap(rr) == GenomicRegionOverlap::ContainsArg) {
//         result.GetIntTag("AS", rr_score);
//       }
//       if (!ar.IsEmpty() &&
//           read_region.GetOverlap(ar) == GenomicRegionOverlap::ContainsArg) {
//         result.GetIntTag("AS", ar_score);
//       }
//     }

//     // Record counts of reads spanning the breakpoints
//     if (al_score > 0 && al_score > rl_score) {
//       al_reads += 1;
//     }
//     if (rl_score > 0 && rl_score > al_score) {
//       rl_reads += 1;
//     }
//     if (ar_score > 0 && ar_score > rr_score) {
//       ar_reads += 1;
//     }
//     if (rr_score > 0 && rr_score > ar_score) {
//       rr_reads += 1;
//     }
//   }

//   auto results = py::dict();
//   results["rl_reads"] = rl_reads;
//   results["rr_reads"] = rr_reads;
//   results["al_reads"] = al_reads;
//   results["ar_reads"] = ar_reads;
//   return results;
// }

}  // namespace npsv