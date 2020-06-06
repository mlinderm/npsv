#include "aligner.hpp"

#include <iostream>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <limits>

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

double LogSumPow(double acc, double prob) {
  double diff = prob - acc;
  if (diff > 100)
    return prob;
  else if (diff < -100)
    return acc;
  else
    return acc + log10(1 + pow(10., diff));
}

 double PhredToProb(double phred) {
    return pow(10.0, phred / -10.0);
  }

  double PhredToLogProb(char quality, double penalty=0.) {
    return -(static_cast<double>(quality) / 10) + penalty;
  }

  double LogProbToPhredQual(double prob, double max_qual=std::numeric_limits<double>::max()) {
    return std::min(log10(1.-pow(10.0, prob)) * -10.0, max_qual);
  }
}  // namespace

AlignmentPair::AlignmentPair(const sl::BamRecord& first,
                             const sl::BamRecord& second,
                             const InsertSizeDistribution& insert_size_dist)
    : first_(&first), second_(&second), score_(0.) {
  
  if (second_->Position() <= first_->Position()) {
    std::swap(first_, second_);
  }

  // https://github.com/nspies/svviz2/blob/44f7bfc75bf84c1db4563d9fd30bf20967d1c825/src/svviz2/remap/alignment.py#L235

  // Scoring algorithm adapted from svviz2:
  // https://github.com/nspies/svviz2/blob/44f7bfc75bf84c1db4563d9fd30bf20967d1c825/src/svviz2/io/readstatistics.py
  {
    int tag_value = 0.;
    first_->GetIntTag("as", tag_value);
    score_ += tag_value;
  }
  {
    int tag_value = 0.;
    second_->GetIntTag("as", tag_value);
    score_ += tag_value;
  }
  
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
  // At this point score_ is log10(P(data|alignments))
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
  for (auto& align1 : first_alignments_) {
    for (auto& align2 : second_alignments_) {
      AlignmentPair pair(align1, align2, insert_size_dist);
      if (!HasBestPair() || pair.Score() > BestAlignmentScore()) {
        best_pair_ = pair;
      }
      total_log_prob_ = LogSumPow(total_log_prob_, pair.Score());
    }
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
  sequence_ = sequence;
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
    ScoreAlignments(fragment.first_, fragment.first_alignments_);
  } else {
    // An existing fragment
    fragment.SetSecond(read);
    
    bwa_.AlignSequence(read.Sequence(), read.Qname(),
                       fragment.second_alignments_, false, 0.9, 10);
    ScoreAlignments(fragment.second_, fragment.second_alignments_);
    
    fragment.SetBestPair(insert_size_dist);   
  }
}

namespace {
  const double kGapOpen = -1.;
  const double kGapExtend = -1.;
}

void AlleleAlignments::ScoreAlignments(const sl::BamRecord& read, sl::BamRecordVector& alignments) {
  const std::string read_sequence(read.Sequence());
  const std::string base_qualities(read.Qualities(0));
  const auto & ref_sequence(sequence_.Seq);

  for (auto & alignment : alignments) {    
    int entry_read_pos = 0;
    int entry_ref_pos = alignment.PositionWithSClips();
    double log_prob = 0; // log10(P(data|alignment))

    sl::Cigar cigar = alignment.GetCigar();
    for (const auto & cigar_entry : cigar) {
      int entry_read_end = entry_read_pos + cigar_entry.Length();
      switch (cigar_entry.Type()) { // MIDNSHPX
        default:
          pyassert(false, "CIGAR entry not implemented");
        case 'S':
          // TODO: Don't penalize shorter soft-clip regions (reduce penalty for < 10 bases)
          for (; entry_read_pos < entry_read_end; entry_read_pos++, entry_ref_pos++) {
            auto quality = static_cast<double>(base_qualities[entry_read_pos]);
            log_prob += PhredToLogProb(base_qualities[entry_read_pos]);
          }
          break;
        case 'M':
          for (; entry_read_pos < entry_read_end; entry_read_pos++, entry_ref_pos++) {
            if (read_sequence[entry_read_pos] == ref_sequence[entry_ref_pos]) {
              auto quality = static_cast<double>(base_qualities[entry_read_pos]);
              log_prob += log10(1 - PhredToProb(quality));
            } else {
              log_prob += PhredToLogProb(base_qualities[entry_read_pos]);
            }
          }
          break;
        case 'I':
          log_prob += PhredToLogProb(base_qualities[entry_read_pos++], kGapOpen);
          for (;entry_read_pos < entry_read_end; entry_read_pos++) {
            log_prob += PhredToLogProb(base_qualities[entry_read_pos], kGapExtend);
          }
          break;
        case 'D':
          log_prob += kGapOpen;
          if (cigar_entry.Length() > 1)
            log_prob += (cigar_entry.Length() - 1) * kGapExtend;
          entry_ref_pos += cigar_entry.Length();
          break;

      }
    }
    alignment.AddIntTag("as", static_cast<int>(log_prob));
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
  ref_.clear();
  alt_.clear();
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
    ar = sl::GenomicRegion(py::cast<std::string>(kwargs["ar_region"]),
                           AltHeader());
  }

  int rl_reads = 0, al_reads = 0, rr_reads = 0, ar_reads = 0;

  // Open the input BAM/SAM/CRAM and any output files
  sl::BamReader reader;
  reader.Open(bam_path);

  sl::BamWriter writer(sl::BAM);
  if (kwargs && kwags.contains("output_bam") && !kwargs["output_bam"].get.is_none()) {
    writer.SetHeader(reader.Header());
    writer.Open(py::cast<std::string>(kwargs["output_bam"]));
    writer.WriteHeader();
  }

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
  pyassert(ref_.size() == alt_.size(),
           "All reads should be aligned to both alleles");

  auto ref_iter = ref_.begin(), ref_end = ref_.end();
  auto alt_iter = alt_.begin(), alt_end = alt_.end();
  while (ref_iter != ref_end && alt_iter != alt_end) {
    pyassert(ref_iter->first == alt_iter->first,
             "Considering alignments for same underlying fragment");
    auto& ref_fragment = ref_iter->second;
    auto& alt_fragment = alt_iter->second;

    
    Fragment::score_type total_log_prob = LogSumPow(ref_fragment.total_log_prob_, alt_fragment.total_log_prob_);
    
    // Compute MAPQ-like quality scores for the reads
    Fragment::score_type rl_score = 0, rr_score = 0, al_score = 0, ar_score = 0;
    if (ref_fragment.BestAlignmentContains(rl)) {
      rl_score = LogProbToPhredQual(ref_fragment.BestAlignmentScore() - total_log_prob, 40);
    }
    if (alt_fragment.BestAlignmentContains(al)) {
      al_score = LogProbToPhredQual(alt_fragment.BestAlignmentScore() - total_log_prob, 40);
    }
    if (!rr.IsEmpty() && ref_fragment.BestAlignmentContains(rr)) {
      rr_score = LogProbToPhredQual(ref_fragment.BestAlignmentScore() - total_log_prob, 40);
    }
    if (!ar.IsEmpty() && alt_fragment.BestAlignmentContains(ar)) {
      ar_score = LogProbToPhredQual(alt_fragment.BestAlignmentScore() - total_log_prob, 40);
    }  
  
     //         sl::BamWriter writer(sl::SAM);
  // writer.SetHeader(AltHeader());
  // writer.Open("/dev/stdout");
  // writer.WriteHeader();
  //         writer.WriteRecord(*(alt_fragment.best_pair_.first_));
  //          writer.WriteRecord(*(alt_fragment.best_pair_.second_));
  //          writer.WriteRecord(*(ref_fragment.best_pair_.first_));
  //       writer.WriteRecord(*(ref_fragment.best_pair_.second_));
  //       writer.Close();

    //if (rl_score > 0 || rr_score > 0 || al_score > 0) {
    // if (ref_iter->first == "HISEQ1:18:H8VC6ADXX:1:1208:7555:15718") {
    //   sl::BamWriter writer(sl::SAM);
    //   writer.SetHeader(AltHeader());
    //   writer.Open("/dev/stdout");
    //   writer.WriteHeader();
    //   writer.WriteRecord(*(alt_fragment.best_pair_.first_));
    //   writer.WriteRecord(*(alt_fragment.best_pair_.second_));
    //   writer.WriteRecord(*(ref_fragment.best_pair_.first_));
    //   writer.WriteRecord(*(ref_fragment.best_pair_.second_));
    //   writer.Close();
    //   std::cout << rl_score << " " << rr_score << " " << al_score << " " << ar_score <<std::endl;
    // }

    if (ar.IsEmpty() && !rr.IsEmpty()) {
      // Deletion
      if ((al_score - rl_score) > 1 && (al_score - rr_score) > 1) {
        al_reads += 1;
      } else {
        if ((rl_score - al_score) > 1) {
          rl_reads += 1;
        } 
        if ((rr_score - al_score) > 1) {
          rr_reads += 1;
        }
      }
    } else {
      throw std::invalid_argument("Realignmet not yet supported for this variant type");
    }

    ref_iter++;
    alt_iter++;
  }

  // Close BAM writer if open
  if (writer.IsOpen()) {
    writer.BuildIndex();
    writer.Close();
  }

  auto results = py::dict();
  results["rl_reads"] = rl_reads;
  results["rr_reads"] = rr_reads;
  results["al_reads"] = al_reads;
  results["ar_reads"] = ar_reads;
  return results;
}

}  // namespace npsv