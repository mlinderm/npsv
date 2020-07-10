#include "aligner.hpp"

#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
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

double LogSumPow(double acc, double prob) {
  double diff = prob - acc;
  if (diff > 100)
    return prob;
  else if (diff < -100)
    return acc;
  else
    return acc + log10(1 + pow(10., diff));
}

double PhredToProb(double phred) { return pow(10.0, phred / -10.0); }

double PhredToLogProb(double quality, double penalty = 0.) {
  return (-quality / 10.) + penalty;
}

double GetDoubleTag(const sl::BamRecord& read, const std::string& tag) {
  uint8_t* p = bam_aux_get(read.raw(), tag.data());
  if (!p) throw std::invalid_argument("Tag does not exist");
  double result = bam_aux2f(p);
  int type = *p++;
  if (type != 'd') throw std::invalid_argument("Tag is not of double type");

  return result;
}
}  // namespace

double InsertSizeDistribution::operator()(int insert_size) const {
  auto entry = density_.find(insert_size);
  if (entry != density_.end()) {
    return entry->second;
  } else {
    // https://stackoverflow.com/a/10848293
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (insert_size - mean_) / std_;
    return inv_sqrt_2pi / std_ * std::exp(-0.5 * a * a);
  }
}

double LogProbToPhredQual(double prob, double max_qual) {
  return std::min(log10(1. - pow(10.0, prob)) * -10.0, max_qual);
}

AlignmentPair::AlignmentPair(const sl::BamRecord& first,
                             const sl::BamRecord& second,
                             const InsertSizeDistribution& insert_size_dist)
    : first_(&first), second_(&second), score_(0.) {
  if (second_->Position() <= first_->Position()) {
    std::swap(first_, second_);
  }

  // Scoring algorithm adapted from svviz2:
  // https://github.com/nspies/svviz2/blob/44f7bfc75bf84c1db4563d9fd30bf20967d1c825/src/svviz2/io/readstatistics.py
  score_ += GetDoubleTag(*first_, "as");
  score_ += GetDoubleTag(*second_, "as");

  if (!Concordant()) {
    score_ -= 10.;
    return;
  }
  auto insert_size_prob = insert_size_dist(InsertSize());
  if (insert_size_prob == 0.) {
    score_ -= 10.;
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
  return second_->PositionWithSClips() + second_->Length() -
         first_->PositionWithSClips();
}

std::ostream& operator<<(std::ostream& os, const AlignmentPair& pair) {
  return (os << *(pair.first_) << *(pair.second_) << pair.score_ << std::endl);
}

Fragment::Fragment(const sl::BamRecord& read)
    : total_log_prob_(std::numeric_limits<score_type>::lowest()) {
  SetRead(read);
}

void Fragment::SetRead(const sl::BamRecord& read) {
  if (read.FirstFlag()) {
    pyassert(first_.isEmpty(),
             "First read in pair already present in fragment");
    first_ = read;
  } else {
    pyassert(second_.isEmpty(),
             "Second read in pair already present in fragment");
    second_ = read;
  }
}

void Fragment::SetBestPair(const InsertSizeDistribution& insert_size_dist) {
  for (auto& align1 : first_alignments_) {
    for (auto& align2 : second_alignments_) {
      AlignmentPair pair(align1, align2, insert_size_dist);
      if (!HasBestPair() || pair.Score() > BestPairScore()) {
        best_pair_ = pair;
      }
      total_log_prob_ = LogSumPow(total_log_prob_, pair.Score());
    }
  }
}

sl::GenomicRegion Fragment::MateQueryRegion() const {
  if (HasFirst() && first_.MateMappedFlag()) {
    return sl::GenomicRegion(first_.MateChrID(), first_.MatePosition(),
                             first_.MatePosition() + 1);
  } else if (HasSecond() && second_.MateMappedFlag()) {
    return sl::GenomicRegion(second_.MateChrID(), second_.MatePosition(),
                             second_.MatePosition() + 1);
  } else
    return sl::GenomicRegion();
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
  if (!emplace_result.second) {
    // An existing fragment
    if ((read.FirstFlag() && fragment.HasFirst()) ||
        (!read.FirstFlag() && fragment.HasSecond())) {
      // We already have both reads in this fragment
      return;
    }
    fragment.SetRead(read);
  }

  if (read.FirstFlag()) {
    pyassert(fragment.HasFirst(), "First read in pair is not defined");
    bwa_.AlignSequence(read.Sequence(), read.Qname(),
                       fragment.first_alignments_, false, 0.9, 10);
    ScoreAlignments(fragment.first_, fragment.first_alignments_);
  } else {
    pyassert(fragment.HasSecond(), "Second read in pair is not defined");
    bwa_.AlignSequence(read.Sequence(), read.Qname(),
                       fragment.second_alignments_, false, 0.9, 10);
    ScoreAlignments(fragment.second_, fragment.second_alignments_);
  }
}

namespace {
// Penalties adapted from svviz2
const double kGapOpen = -1.;
const double kGapExtend = -1.;

void AddDoubleTag(sl::BamRecord& read, const std::string& tag, double val) {
  bam_aux_append(read.raw(), tag.data(), 'd', sizeof(double), (uint8_t*)&val);
}

// svviz2 rescales all base qualities
double RescaleQuality(char quality, double scale = 0.25) {
  return scale * static_cast<double>(quality);
}
}  // namespace

AlleleAlignments::score_type AlleleAlignments::ScoreAlignment(
    const std::string& read_sequence, const std::string& base_qualities,
    const std::string& ref_sequence, const sl::BamRecord& alignment) {
  int entry_read_pos = 0;
  int entry_ref_pos = alignment.PositionWithSClips();
  score_type log_prob = 0;  // log10(P(data|alignment))

  sl::Cigar cigar = alignment.GetCigar();
  for (const auto& cigar_entry : cigar) {
    int entry_read_end = entry_read_pos + cigar_entry.Length();
    switch (cigar_entry.Type()) {  // MIDNSHPX
      default:
        throw std::invalid_argument("CIGAR entry not implemented");
      case 'S':
        // TODO: Don't penalize shorter soft-clip regions (reduce penalty for <
        // 10 bases)
        for (; entry_read_pos < entry_read_end;
             entry_read_pos++, entry_ref_pos++) {
          log_prob +=
              PhredToLogProb(RescaleQuality(base_qualities[entry_read_pos]));
        }
        break;
      case 'M':
        for (; entry_read_pos < entry_read_end;
             entry_read_pos++, entry_ref_pos++) {
          if (read_sequence[entry_read_pos] == ref_sequence[entry_ref_pos]) {
            auto quality = RescaleQuality(base_qualities[entry_read_pos]);
            log_prob += log10(1. - PhredToProb(quality));
          } else {
            log_prob +=
                PhredToLogProb(RescaleQuality(base_qualities[entry_read_pos]));
          }
        }
        break;
      case 'I':
        log_prob += PhredToLogProb(
            RescaleQuality(base_qualities[entry_read_pos++]), kGapOpen);
        for (; entry_read_pos < entry_read_end; entry_read_pos++) {
          log_prob += PhredToLogProb(
              RescaleQuality(base_qualities[entry_read_pos]), kGapExtend);
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

  return log_prob;
}

void AlleleAlignments::ScoreAlignments(const sl::BamRecord& read,
                                       sl::BamRecordVector& alignments) {
  const std::string read_sequence(read.Sequence());
  const std::string base_qualities(read.Qualities(0));
  const std::string& ref_sequence(sequence_.Seq);

  for (auto& alignment : alignments) {
    auto log_prob =
        ScoreAlignment(read_sequence, base_qualities, ref_sequence, alignment);
    AddDoubleTag(alignment, "as", log_prob);
  }
}

Overlap::Overlap(const Fragment& fragment, const sl::GenomicRegion& breakpoint,
                 OverlapAllele allele, score_type total_log_prob,
                 bool count_straddle)
    : fragment_(&fragment),
      breakpoint_(&breakpoint),
      allele_(allele),
      kind_(OverlapKind::None),
      quality_score_(0) {
  // TODO: Handle fragments with just a single read?
  const AlignmentPair& best_pair(fragment_->BestPair());
  if (!best_pair.Valid()) {
    // Use default alignment of "None"
    return;
  }

  // TODO: Add minimum overlap?
  auto first_region = best_pair.first_->AsGenomicRegion();
  auto second_region = best_pair.second_->AsGenomicRegion();
  if (first_region.GetOverlap(*breakpoint_) ==
          GenomicRegionOverlap::ContainsArg ||
      second_region.GetOverlap(*breakpoint_) ==
          GenomicRegionOverlap::ContainsArg) {
    kind_ = OverlapKind::ReadContains;
  } else if (count_straddle && first_region.chr == breakpoint_->chr &&
             second_region.chr == breakpoint_->chr) {
    // Read could straddle...
    sl::GenomicRegion fragment_region(
        first_region.chr, std::min(first_region.pos1, second_region.pos1),
        std::max(first_region.pos2, second_region.pos2));
    if (fragment_region.GetOverlap(*breakpoint_) ==
        GenomicRegionOverlap::ContainsArg) {
      kind_ = OverlapKind::FragmentContains;
    }
  }

  if (kind_ != OverlapKind::None) {
    quality_score_ = LogProbToPhredQual(best_pair.Score() - total_log_prob, 40);
  }
}

std::string Overlap::FormattedBreakpoint() const {
  std::stringstream ss;
  ss << (IsRef() ? "ref:" : "alt:") << breakpoint_->pos1 << "-"
     << breakpoint_->pos2;
  return ss.str();
}

void Overlap::PushTaggedReads(sl::BamRecordVector& reads) const {
  std::string ov_value(FormattedBreakpoint());
  reads.push_back(fragment_->first_);
  reads.back().AddZTag("ov", ov_value);
  AddDoubleTag(reads.back(), "as",
               GetDoubleTag(*(fragment_->BestPair().first_), "as"));
  AddDoubleTag(reads.back(), "fs", quality_score_);
  reads.push_back(fragment_->second_);
  reads.back().AddZTag("ov", ov_value);
  AddDoubleTag(reads.back(), "fs", quality_score_);
  AddDoubleTag(reads.back(), "as",
               GetDoubleTag(*(fragment_->BestPair().second_), "as"));
}

Realigner::Realigner(
    const std::string& fasta_path, double insert_size_mean,
    double insert_size_std,
    const InsertSizeDistribution::density_type& insert_size_density)
    : insert_size_dist_(insert_size_mean, insert_size_std,
                        insert_size_density) {
  pyassert(!ref_.IsInitialized() && !alt_.IsInitialized(),
           "BWA wrapper should start out uninitialized");

  // Load alleles from a FASTA file
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

void Realigner::Clear() {
  ref_.clear();
  alt_.clear();
}

void Realigner::Realign(const sl::BamRecord& read) {
  ref_.Align(read, insert_size_dist_);
  alt_.Align(read, insert_size_dist_);
}

std::tuple<std::map<std::string,int>, std::map<std::string,std::vector<std::string> > > Realigner::CountAlignments(const std::string& bam_path,
                                    const std::string& rl_region,
                                    const std::string& al_region,
                                    py::kwargs kwargs) {
  // Get optional arguments passed from Python
  Fragment::score_type min_score_delta = 1.;
  if (kwargs && kwargs.contains("min_score_delta")) {
    min_score_delta = py::cast<Fragment::score_type>(kwargs["min_score_delta"]);
  }

  bool count_straddle = true;
  if (kwargs && kwargs.contains("count_straddle")) {
    count_straddle = py::cast<bool>(kwargs["count_straddle"]);
  }

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

  int rl_reads = 0, al_reads = 0, rr_reads = 0, ar_reads = 0, amb_reads = 0;
  
  // Track the read names overlapping each breakpoint
  std::map<std::string, std::vector<std::string> > overlap_read_names {
    { "rl", std::vector<std::string>() },
    { "rr", std::vector<std::string>() },
    { "al", std::vector<std::string>() },
    { "ar", std::vector<std::string>() },
  };

  // Open the input BAM/SAM/CRAM and any output files
  sl::BamReader reader;
  reader.Open(bam_path);

  sl::BamWriter overlap_writer(sl::BAM);
  sl::BamRecordVector overlap_reads;
  if (kwargs && kwargs.contains("overlap_bam") &&
      !kwargs["overlap_bam"].is_none()) {
    overlap_writer.SetHeader(reader.Header());
    overlap_writer.Open(py::cast<std::string>(kwargs["overlap_bam"]));
    overlap_writer.WriteHeader();
  }

  // Subset the reader by the specified regions
  sl::GenomicRegion reader_region;
  if (kwargs && kwargs.contains("region")) {
    // Restrict input BAM to specific region
    reader_region = sl::GenomicRegion(py::cast<std::string>(kwargs["region"]),
                                      reader.Header());
    reader.SetRegion(reader_region);
  }

  sl::BamRecord read;
  while (reader.GetNextRecord(read)) {
    if (read.DuplicateFlag() || read.SecondaryFlag() ||
        read.SupplementaryFlag())
      continue;  // Skip duplicate or secondary/supplemental alignments

    // Realign read to ref. and/or alt. alleles
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

    // // If un-paired, try to find and align the mate
    // if (!ref_fragment.IsPaired()) {
    //   pyassert(!alt_fragment.IsPaired(), "Ref and Alt fragments should have
    //   the same pairing status");

    //   auto query_region = ref_fragment.MateQueryRegion();
    //   if (!query_region.IsEmpty() && (reader_region.IsEmpty() ||
    //   reader_region.GetOverlap(query_region) ==
    //   GenomicRegionOverlap::NoOverlap)) {
    //     // Only look for mate in regions we haven't already explored
    //     reader.SetRegion(query_region);
    //     while (reader.GetNextRecord(read)) {
    //       if (read.SecondaryFlag() || read.SupplementaryFlag())
    //         continue;  // Skip secondary/supplemental alignments
    //       if (read.Qname() == ref_iter->first) {
    //         // Is this the mate we are looking for?
    //         if ((ref_fragment.HasFirst() && read.FirstFlag()) ||
    //         (ref_fragment.HasSecond() && !read.FirstFlag())) {
    //           continue;
    //         }
    //         if (read.DuplicateFlag()) {
    //           break;  // Mate is a duplicate
    //         }

    //         // We found the other read in the pair
    //         Realign(read);
    //         pyassert(ref_fragment.IsPaired() && alt_fragment.IsPaired(),
    //         "After realignment both ref and alt should be paired"); break;
    //       }
    //     }
    //   }
    // }

    // Find the best pairing among all of the alignments
    ref_fragment.SetBestPair(insert_size_dist_);
    alt_fragment.SetBestPair(insert_size_dist_);

    // Compute total probability across all alignments for this fragment (to
    // compute posterior)
    Fragment::score_type total_log_prob =
        LogSumPow(ref_fragment.total_log_prob_, alt_fragment.total_log_prob_);

    // Determine overlap and scores for all breakpoints
    Overlap rl_ov(ref_fragment, rl, Overlap::OverlapAllele::RefLeft,
                  total_log_prob, count_straddle);
    Overlap al_ov(alt_fragment, al, Overlap::OverlapAllele::AltLeft,
                  total_log_prob, count_straddle);

    Overlap rr_ov, ar_ov;  // Optional breakpoints
    if (!rr.IsEmpty()) {
      rr_ov = {ref_fragment, rr, Overlap::OverlapAllele::RefRight,
               total_log_prob, count_straddle};
    }
    if (!ar.IsEmpty()) {
      ar_ov = {alt_fragment, ar, Overlap::OverlapAllele::AltRight,
               total_log_prob, count_straddle};
    }

    bool any_overlap = rl_ov.HasOverlap() || al_ov.HasOverlap() ||
                       rr_ov.HasOverlap() || ar_ov.HasOverlap();
    if (any_overlap) {
      Overlap max_ref = std::max(rl_ov, rr_ov);
      Overlap max_alt = std::max(al_ov, ar_ov);

      auto delta = max_alt.QualityScore() - max_ref.QualityScore();
      if (max_alt.HasOverlap() && delta > min_score_delta) {
        if (al_ov.HasOverlap()) {
          al_reads += 1;
          overlap_read_names["al"].emplace_back(alt_iter->first);
        }
        if (ar_ov.HasOverlap()) {
          ar_reads += 1;
          overlap_read_names["ar"].emplace_back(alt_iter->first);
        }

        if (overlap_writer.IsOpen()) max_alt.PushTaggedReads(overlap_reads);
      } else if (max_ref.HasOverlap() && delta < -min_score_delta) {
        if (rl_ov.HasOverlap()) {
          rl_reads += 1;
          overlap_read_names["rl"].emplace_back(ref_iter->first);
        }
        if (rr_ov.HasOverlap()) {
          rr_reads += 1;
          overlap_read_names["rr"].emplace_back(ref_iter->first);
        }

        if (overlap_writer.IsOpen()) max_ref.PushTaggedReads(overlap_reads);
      } else {
        // Everything else is an ambiguous overlap
        // TODO: Write out ambiguous alignments to debugging BAM file, record
        // overlap for ambiguous alignments
        amb_reads += 1;
      }
    }

    ref_iter++;
    alt_iter++;
  }

  // Close overlap BAM writer if open (writing out index file)
  if (overlap_writer.IsOpen()) {
    std::sort(overlap_reads.begin(), overlap_reads.end(),
              sl::BamRecordSort::ByReadPosition());
    for (const auto& read : overlap_reads) overlap_writer.WriteRecord(read);
    overlap_writer.Close();
    overlap_writer.BuildIndex();
  }

  std::map<std::string,int> results;
  results["rl_reads"] = rl_reads;
  results["rr_reads"] = rr_reads;
  results["al_reads"] = al_reads;
  results["ar_reads"] = ar_reads;
  results["amb_reads"] = amb_reads;
  
  return std::make_tuple(results, overlap_read_names);
}

namespace test {

std::vector<AlleleAlignments::score_type> TestScoreAlignment(
    const std::string& ref_sequence, const std::string& alignment_path) {
  // Open the input BAM/SAM/CRAM and any output files
  sl::BamReader reader;
  reader.Open(alignment_path);
  std::vector<AlleleAlignments::score_type> scores;

  sl::BamRecord read;
  while (reader.GetNextRecord(read)) {
    auto log_prob = AlleleAlignments::ScoreAlignment(
        read.Sequence(), read.Qualities(0), ref_sequence, read);
    scores.push_back(log_prob);
  }

  return scores;
}

}  // namespace test

}  // namespace npsv