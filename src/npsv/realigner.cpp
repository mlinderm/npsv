#include "utility.hpp"
#include "realigner.hpp"

#include <cmath>

#include "SeqLib/FastqReader.h"


namespace {
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

double LogProbToPhredQual(double prob, double max_qual) {
  return std::min(log10(1. - pow(10.0, prob)) * -10.0, max_qual);
}

double GetDoubleTag(const sl::BamRecord& read, const std::string& tag) {
  uint8_t* p = bam_aux_get(read.raw(), tag.data());
  if (!p) throw std::invalid_argument("Tag does not exist");
  double result = bam_aux2f(p);
  int type = *p++;
  if (type != 'd') throw std::invalid_argument("Tag is not of double type");

  return result;
}  

int OverlapRegion(int32_t read_start, int32_t read_end,
            const sl::GenomicRegion& region) {
  // Convert 1-indexed inclusive to 0-index half-open
  int region_start = region.pos1 - 1;
  int region_end = region.pos2;

  // This is not as precise as actually counting alignment bases, e.g.
  // https://github.com/pysam-developers/pysam/blob/76f6f5c4e3834fcdb7830755dea713423cfe1c15/pysam/libcalignedsegment.pyx#L2020
  // but simpler and does not require us to have loaded both reads.
  int overlap =
      std::min(read_end, region_end) - std::max(read_start, region_start);
  return overlap >= 0 ? overlap : 0;
}

int ReadOverlapRegion(const sl::BamRecord& read, const sl::GenomicRegion& region) {
  if (read.ChrID() != region.chr) return 0;
  return OverlapRegion(read.Position(), read.PositionEnd(), region);
}

int MateOverlapRegion(const sl::BamRecord& read, const sl::GenomicRegion& region) {
  if (read.MateChrID() != region.chr) return 0;
  return OverlapRegion(read.MatePosition(), read.PositionEndMate(), region);
}

}  // namespace

namespace npsv {
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
  
AlignedFragment::AlignedFragment(const sl::BamRecord& read) {
  if (read.FirstFlag()) {
    first_ = read;
    left_ = &first_;
  } else {
    second_ = read;
    left_ = &second_;
  }
}

bool AlignedFragment::IsProperPair() const {
  return (HasFirst() && first_.ProperPair()) ||
         (HasSecond() && second_.ProperPair());
}

int AlignedFragment::InsertSize() const {
  return left_->FullInsertSize();
}

int32_t AlignedFragment::LeftPos1() const {
  // Convert to 1-indexed
  return left_->Position() + 1;
}

int32_t AlignedFragment::LeftPos2() const {
  // Convert to 1-indexed
  return left_->PositionEnd() + 1;
}

int32_t AlignedFragment::RightPos1() const {
  // Convert to 1-indexed
  return left_->MatePosition() + 1;
}

int32_t AlignedFragment::RightPos2() const {
  // Convert to 1-indexed
  return left_->PositionEndMate() + 1;
}

void AlignedFragment::SetRead(const sl::BamRecord& read) {
  if (read.FirstFlag() && first_.isEmpty()) {
    first_ = read;
  } else if (!read.FirstFlag() && second_.isEmpty()) {
    second_ = read;
  }
  if (HasFirst() && HasSecond()) {
    pyassert(first_.MatePosition() == second_.Position(), "Read and mate positions inconsistent");
    left_ = first_.Position() <  second_.Position() ? &first_ : &second_;
  }
}

bool AlignedFragment::Straddles(const sl::GenomicRegion& left_region,
                                const sl::GenomicRegion& right_region, int min_overlap) const {
  // This is an expansive definition of straddling that only requires minimum overlap of the region
  // (as opposed to the two reads starting in the corresponding regions).
  return ReadOverlapRegion(*left_, left_region) >= min_overlap &&
         MateOverlapRegion(*left_, right_region) >= min_overlap;
}

bool AlignedFragment::StrictStraddles(const sl::GenomicRegion& left_region,
                                      const sl::GenomicRegion& right_region,
                                      int min_overlap) const {
  // A stricter definition of straddling that requires the reads to start in the respective regions
  return LeftPos1() >= left_region.pos1 && LeftPos1() <= (left_region.pos2 - min_overlap + 1) &&
         RightPos1() >= right_region.pos1 && RightPos1() <= (right_region.pos2 - min_overlap + 1);
}

double AlignedFragment::ProbMapQ() const {
  pyassert(HasFirst() || HasSecond(), "Calculating MapQ of emtpy fragment");
  double result = 1.;
  if (HasFirst())
    result *= (1. - std::pow(10., -first_.MapQuality() / 10.));
  if (HasSecond())
    result *= (1. - std::pow(10., -second_.MapQuality() / 10.));
  return result;
}

namespace {
  std::tuple<int, int> CountLeftAndRightClippedBases(const sl::BamRecord& read) {
    sl::Cigar cigar = read.GetCigar();
    // Count clips on and left and right edge of read
    int left_clipped_bases = 0;
    for (int i=0; i<cigar.size(); i++) {
      const auto& field = cigar[i];
      if (field.Type() == 'S' || field.Type() == 'H')
        left_clipped_bases += field.Length();
      else
        break;
    }

    int right_clipped_bases = 0;
    for (int i=cigar.size()-1; i>=0; i--) {
      const auto& field = cigar[i];
      if (field.Type() == 'S' || field.Type() == 'H')
        right_clipped_bases += field.Length();
      else
        break;
    }
    return std::make_tuple(left_clipped_bases, right_clipped_bases);
  }
}


std::tuple<int, int, int, int> AlignedFragment::ClipCount(const sl::GenomicRegion& breakpoint, int min_clip) const {
  int total_reads = 0, left_clipped_reads = 0, right_clipped_reads = 0, both_clipped_reads = 0;
  
  if (HasFirst() && ReadOverlapRegion(first_, breakpoint) > 0) {
    total_reads += 1;
    int left_clipped_bases, right_clipped_bases;
    std::tie(left_clipped_bases, right_clipped_bases) = CountLeftAndRightClippedBases(first_);
    bool left_clipped = left_clipped_bases >= min_clip, right_clipped = right_clipped_bases >= min_clip;
    if (left_clipped && right_clipped)
      both_clipped_reads += 1;
    else if (left_clipped)
      left_clipped_reads += 1;
    else if (right_clipped)
      right_clipped_reads += 1;
  }
  
  if (HasSecond() && ReadOverlapRegion(second_, breakpoint) > 0) {
    total_reads += 1;
    int left_clipped_bases, right_clipped_bases;
    std::tie(left_clipped_bases, right_clipped_bases) = CountLeftAndRightClippedBases(second_);
    bool left_clipped = left_clipped_bases >= min_clip, right_clipped = right_clipped_bases >= min_clip;
    if (left_clipped && right_clipped)
      both_clipped_reads += 1;
    else if (left_clipped)
      left_clipped_reads += 1;
    else if (right_clipped)
      right_clipped_reads += 1;
  }
  
  return std::make_tuple(total_reads, left_clipped_reads, right_clipped_reads, both_clipped_reads);
}

std::ostream& operator<<(std::ostream& os, const AlignedFragment& fragment) {
  return (os << fragment.first_ << std::endl << fragment.second_);
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

double ScoreAlignment(const std::string& read_sequence,
                      const std::string& base_qualities,
                      const std::string& ref_sequence,
                      const sl::BamRecord& alignment) {
  int entry_read_pos = 0;
  int entry_ref_pos = alignment.PositionWithSClips();
  double log_prob = 0;  // log10(P(data|alignment))

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

void ScoreAlignments(const IndexedSequence& index, const sl::BamRecord& read,
                     sl::BamRecordVector& alignments) {
  const std::string read_seq(read.Sequence());
  const std::string base_qualities(read.Qualities(0));
  const std::string& ref_seq(index.Sequence());

  for (auto& alignment : alignments) {
    auto log_prob = ScoreAlignment(read_seq, base_qualities, ref_seq, alignment);
    AddDoubleTag(alignment, "as", log_prob);
  }
}

void RealignRead(const IndexedSequence& index, const sl::BamRecord& read,
                 sl::BamRecordVector& alignments) {
  index.AlignSequence(read, alignments);
  ScoreAlignments(index, read, alignments);
}

}  // namespace

RealignedReadPair::RealignedReadPair(const sl::BamRecord& first)
      : left_(&first), right_(nullptr) {
  score_ = GetDoubleTag(*left_, "as");
}

RealignedReadPair::RealignedReadPair(const sl::BamRecord& first,
                             const sl::BamRecord& second,
                             const InsertSizeDistribution& insert_dist)
    : left_(&first), right_(&second), score_(0.) {
  if (left_->Position() > right_->Position()) {
    std::swap(left_, right_);
  }
  // Scoring algorithm adapted from svviz2:
  // https://github.com/nspies/svviz2/blob/44f7bfc75bf84c1db4563d9fd30bf20967d1c825/src/svviz2/io/readstatistics.py
  score_ += GetDoubleTag(*left_, "as");
  score_ += GetDoubleTag(*right_, "as");

  if (!Concordant()) {
    score_ -= 10.;
    return;
  }
  auto insert_size_prob = insert_dist(InsertSize());
  if (insert_size_prob == 0.) {
    score_ -= 10.;
    return;
  }
  score_ += log10(insert_size_prob);
}

namespace {
  sl::GenomicRegion ReadRegion(const sl::BamRecord& read) {
    // The read region is 0-indexed with a non-inclusive end coordinate, but SeqLib's operations on GenomicRegions
    // assume 1-indexed regions with inclusive start and end coordinates. So we convert those coordinates here.
    return sl::GenomicRegion(read.ChrID(), read.Position() + 1, read.PositionEnd());
  }
}

sl::GenomicRegion RealignedReadPair::LeftReadRegion() const {
  return left_ ? ReadRegion(*left_) : sl::GenomicRegion();
}
sl::GenomicRegion RealignedReadPair::RightReadRegion() const {
  return right_ ? ReadRegion(*right_) : sl::GenomicRegion();
}


int32_t RealignedReadPair::InsertSize() const {
  return right_->PositionWithSClips() + right_->Length() - left_->PositionWithSClips();
}

bool RealignedReadPair::Concordant() const {
  if (left_->ChrID() != right_->ChrID()) return false;

  // TODO: Check orientation

  return true;
}

std::ostream& operator<<(std::ostream& os, const RealignedReadPair& pair) {
  if (pair.left_)
    os << *pair.left_ << std::endl;
  if (pair.right_)
    os << *pair.right_ << std::endl;
  std::cerr << pair.score_ << std::endl;
}

BreakpointOverlap::BreakpointOverlap(const RealignedReadPair& read_pair,
                                 const sl::GenomicRegion& breakpoint,
                                  Allele allele,
                                 score_type total_log_prob, bool count_straddle)
    : breakpoint_(&breakpoint),
      allele_(allele),
      kind_(Overlap::None),
      quality_score_(0) {
  sl::GenomicRegion left_read_region(read_pair.LeftReadRegion());
  sl::GenomicRegion right_read_region(read_pair.RightReadRegion());
  
  // TODO: Incorporate minimum anchor
  if (left_read_region.GetOverlap(*breakpoint_) ==
          GenomicRegionOverlap::ContainsArg ||
      right_read_region.GetOverlap(*breakpoint_) ==
          GenomicRegionOverlap::ContainsArg) {
    kind_ = Overlap::ReadContains;
  } else if (count_straddle && left_read_region.chr == breakpoint_->chr &&
             right_read_region.chr == breakpoint_->chr) {
    // Read could straddle...
    pyassert(left_read_region.pos1 <= right_read_region.pos1,
             "Reads in pair not properly sorted");
    sl::GenomicRegion fragment_region(
        left_read_region.chr, left_read_region.pos1, right_read_region.pos2);
    if (fragment_region.GetOverlap(*breakpoint_) ==
        GenomicRegionOverlap::ContainsArg) {
      kind_ = Overlap::FragmentContains;
    }
  }

  if (kind_ != Overlap::None) {
    quality_score_ = LogProbToPhredQual(read_pair.Score() - total_log_prob, 40);
  }
}

AlleleOverlap::AlleleOverlap(const RealignedReadPair& read_pair,
                             int allele_index,
                             const sl::GenomicRegion& left_breakpoint,
                             const sl::GenomicRegion& right_breakpoint,
                             BreakpointOverlap::Allele allele,
                             score_type total_log_prob,
                             bool count_straddle = true)
    : read_pair_(&read_pair),
      allele_index_(allele_index),
      left_overlap_(read_pair, left_breakpoint,
                    allele & BreakpointOverlap::Allele::Left, total_log_prob,
                    count_straddle),
      right_overlap_(read_pair, right_breakpoint,
                     allele & BreakpointOverlap::Allele::Right, total_log_prob,
                     count_straddle)
{}

RealignedFragment::RealignedFragment(const AlignedFragment& original_alignment,
                                     const IndexedSequence& index,
                                     const InsertSizeDistribution& insert_dist)
    : original_alignment_(original_alignment),
      total_log_prob_(std::numeric_limits<score_type>::lowest()) {
  if (original_alignment_.HasFirst()) {
    RealignRead(index, original_alignment_.FirstRead(), first_alignments_);
  }
  if (original_alignment_.HasSecond()) {
    RealignRead(index, original_alignment_.SecondRead(), second_alignments_);
  }

  // Construct and score possible alignment pairs from these individual
  if (!first_alignments_.empty() && !second_alignments_.empty()) {
    for (auto& align1 : first_alignments_) {
      for (auto& align2 : second_alignments_) {
        read_pairs_.emplace_back(align1, align2, insert_dist);
      }
    }
  } 
  // Previous NPSV only considered actual pairs, but incorporating 
  // singletons slightly reduced accuracy
  // else {
  //   // Handle situation with singleton reads
  //   for (auto& align : first_alignments_)
  //     read_pairs_.emplace_back(align);
  //   for (auto& align : second_alignments_)
  //     read_pairs_.emplace_back(align);
  // }

  // Sort alignments in descending order by score
  std::sort(read_pairs_.begin(), read_pairs_.end(), std::greater<>());

  for (const auto & pair : read_pairs_) {
    total_log_prob_ = LogSumPow(total_log_prob_, pair.Score());
  }
}

IndexedSequence::IndexedSequence(const sl::UnalignedSequence& sequence) {
  Initialize(sequence);
}

void IndexedSequence::Initialize(const sl::UnalignedSequence& sequence) {
  pyassert(!IsInitialized(), "BWA should not previously have been initialized");
  sequence_ = sequence;
  bwa_.ConstructIndex({sequence});
}

void IndexedSequence::AlignSequence(const sl::BamRecord& read, sl::BamRecordVector& alignments) const {
  bwa_.AlignSequence(read.Sequence(), read.Qname(), alignments, false, 0.9, 10);
}

RealignedFragments::RealignedFragments(
    const std::string& fasta_path,
    double insert_size_mean, double insert_size_std,
    const InsertSizeDistribution::density_type& insert_size_density,
    const std::string& bam_path)
    : insert_size_dist_(insert_size_mean, insert_size_std,
                        insert_size_density) {
  
  // Load alleles from a FASTA file
  sl::FastqReader contigs(fasta_path);
  sl::UnalignedSequence next_sequence;

  // We assumed the first sequence is the reference sequence
  pyassert(contigs.GetNextSequence(next_sequence), "Reference sequence not present in the FASTA");
  ref_index_.Initialize(next_sequence);
  
  // The remaining sequences at the alternate sequences
  while (contigs.GetNextSequence(next_sequence)) {
    alt_indexes_.emplace_back(next_sequence);
  }

  // Initialize the BAM reader
  reader_.Open(bam_path);
}

int RealignedFragments::GatherReads(
    const std::string& region,
    int max_reads = std::numeric_limits<int>::max()) {
  // Subset the reader by the specified regions
  sl::GenomicRegion reader_region(region, reader_.Header());
  reader_.SetRegion(reader_region);

  int read_count = 0;
  sl::BamRecord read;
  while (read_count < max_reads && reader_.GetNextRecord(read)) {
    if (read.DuplicateFlag() || read.SecondaryFlag() ||
        read.SupplementaryFlag())
      continue;  // Skip duplicate or secondary/supplemental alignments

    AddRead(read);
    read_count++;
  }

  return read_count;
}

void RealignedFragments::AddRead(const sl::BamRecord& read) {
  auto emplace_result = fragments_.emplace(std::piecewise_construct,
                                           std::forward_as_tuple(read.Qname()),
                                           std::forward_as_tuple(read));

  auto& fragment = emplace_result.first->second;
  if (!emplace_result.second) {
    // An existing fragment
    fragment.SetRead(read);
  }
}

std::tuple<double, double, int> RealignedFragments::FragmentInsertProbabilities(const AlignedFragment& fragment, int alt_size_delta) const {
  int ref_insert_size = fragment.InsertSize();
  double ref_prob = insert_size_dist_(ref_insert_size);
      
  int alt_insert_size = ref_insert_size + alt_size_delta;
  double alt_prob = alt_insert_size > 0 ? insert_size_dist_(alt_insert_size) : 0.;

  return std::make_tuple(ref_prob, alt_prob, ref_insert_size);  
}

double RealignedFragments::FragmentProbConcordance(const AlignedFragment& fragment, int alt_size_delta) const {
  double ref_prob, alt_prob;
  std::tie(ref_prob, alt_prob, std::ignore) = FragmentInsertProbabilities(fragment, alt_size_delta);
  return FragmentProbConcordance(ref_prob, alt_prob);
}

double RealignedFragments::FragmentProbConcordance(double ref_prob, double alt_prob) {
  return CONC_PRIOR * ref_prob / (CONC_PRIOR * ref_prob+ DISC_PRIOR * alt_prob);
}

std::map<std::string, double> RealignedFragments::CountPipelineStraddlers(
    const std::string& left_breakpoint, const std::string& right_breakpoint,
    int flank, int alt_size_delta, double z_threshold, int min_overlap) const {
  // SeqLib GenomicRegions are 1-indexed with inclusive start and end coordinates
  sl::GenomicRegion left_breakpoint_region(left_breakpoint, reader_.Header());
  sl::GenomicRegion right_breakpoint_region(right_breakpoint, reader_.Header());
  
  if (left_breakpoint_region.Width() != 2 || right_breakpoint_region.Width() != 2)
    std::cerr << left_breakpoint_region << " " << right_breakpoint_region << std::endl;

  pyassert(left_breakpoint_region.Width() == 2 && right_breakpoint_region.Width() == 2, "Breakpoints should specify inclusive coordinates on each side of the breakpoint");
  pyassert(flank >= 1, "Flank size must be >= 1");

  sl::GenomicRegion left_region(left_breakpoint_region.chr, left_breakpoint_region.pos1-flank+1, left_breakpoint_region.pos1);
  sl::GenomicRegion right_region(right_breakpoint_region.chr, right_breakpoint_region.pos2, right_breakpoint_region.pos2+flank-1);
  
  // Anchor regions within the event (on left and right side)
  sl::GenomicRegion left_event_region, right_event_region;
  if (left_breakpoint_region.chr == right_breakpoint_region.chr && left_breakpoint_region.pos2 <= right_breakpoint_region.pos1) {
    left_event_region = { left_breakpoint_region.chr, left_breakpoint_region.pos2, std::min(left_breakpoint_region.pos2+flank-1, right_breakpoint_region.pos2) };
    right_event_region = { right_breakpoint_region.chr, std::max(right_breakpoint_region.pos1-flank+1, left_breakpoint_region.pos2), right_breakpoint_region.pos1 };
  }

  int insert_count = 0, insert_upper = 0, insert_lower = 0;
  double ref_count = 0., alt_count = 0., ref_weighted_count = 0., alt_weighted_count = 0., ref_conc_count = 0., alt_conc_count = 0.;

  for (auto& named_fragment : fragments_) {
    auto& fragment = named_fragment.second;
    if (!fragment.IsProperPair()) continue;

    // Read straddles the entire event
    if (fragment.StrictStraddles(left_region, right_region, min_overlap)) {
      insert_count += 1;
      
      int ref_insert_size;
      double ref_prob, alt_prob;
      std::tie(ref_prob, alt_prob, ref_insert_size) = FragmentInsertProbabilities(fragment, alt_size_delta);

      if (ref_prob + alt_prob > 0.) {
        if (alt_prob > ref_prob)
          alt_count += 1;
        else
          ref_count += 1;
        
        double p_alt = alt_prob / (ref_prob + alt_prob);
        ref_weighted_count += 1 - p_alt;
        alt_weighted_count += p_alt;

        // Adapted from SVTyper read pair counting
        double p_conc = FragmentProbConcordance(ref_prob, alt_prob);
        if (p_conc > 0.5) {
          ref_conc_count += fragment.ProbMapQ();
        } else {
          alt_conc_count += fragment.ProbMapQ();
        }
      }

      // Adapted from https://www.sciencedirect.com/science/article/pii/S0092867418316337#sec4
      double zscore = insert_size_dist_.ZScore(ref_insert_size);
      if (zscore < -z_threshold)
        insert_lower += 1;
      else if (zscore > z_threshold)
        insert_upper += 1;
    
      continue;
    }

    if (!left_event_region.IsEmpty() && !right_event_region.IsEmpty()) {
      // Read straddles one of the breakpoints (shouldn't straddle both). This is more restrictive than SVTyper
      // which will count fragments where both reads start outside the event, but one read is anchored in the event.
      bool ref_straddle_left = fragment.StrictStraddles(left_region, left_event_region, min_overlap);
      bool ref_straddle_right = fragment.StrictStraddles(right_event_region, right_region, min_overlap);
      
      // ref_straddle_left = ref_straddle_left && left_breakpoint_region.pos2 <= fragment.RightPos1() && fragment.RightPos1() <= right_breakpoint_region.pos1;
      // ref_straddle_right = ref_straddle_right && left_breakpoint_region.pos2 <= fragment.LeftPos2() && fragment.LeftPos2() <= right_breakpoint_region.pos1;
      
      if (ref_straddle_left || ref_straddle_right) {
        pyassert(ref_straddle_left ^ ref_straddle_right, "Reads straddling both breakpoints should have been filtered out");
        ref_count += 0.5;
        ref_weighted_count += 0.5;  // p_alt is by definition 0
        
        double p_conc = FragmentProbConcordance(fragment, alt_size_delta);
        if (p_conc > 0.5) {
          ref_conc_count += fragment.ProbMapQ();
        }
      }
    }

  }

  std::map<std::string, double> results;
  results["insert_count"] = insert_count;
  results["insert_lower"] = insert_lower;
  results["insert_upper"] = insert_upper;
  results["ref_count"] = ref_count;
  results["alt_count"] = alt_count;
  results["ref_weighted_count"] = ref_weighted_count;
  results["alt_weighted_count"] = alt_weighted_count;
  results["ref_conc_count"] = ref_conc_count;
  results["alt_conc_count"] = alt_conc_count;
  return results;
}

namespace {
  sl::GenomicRegion BreakpointToGenomicRegion(const std::string& region, const sl::BamHeader& header) {
      return (region.empty()) ? sl::GenomicRegion() : sl::GenomicRegion(region, header);
  }
}

std::map<std::string, int> RealignedFragments::CountPipelineClippedReads(const std::string& breakpoint, int min_clip) const {
  sl::GenomicRegion breakpoint_region = BreakpointToGenomicRegion(breakpoint, reader_.Header());
  
  int total_reads = 0, left_clipped_reads = 0, right_clipped_reads = 0, both_clipped_reads = 0;
  for (auto& named_fragment : fragments_) {
    auto& fragment = named_fragment.second;

    // Check if read overlaps breakpoint and count clips
    int total, left_clipped, right_clipped, both_clipped;
    std::tie(total, left_clipped, right_clipped, both_clipped) = fragment.ClipCount(breakpoint_region, min_clip);
    total_reads += total;
    left_clipped_reads += left_clipped;
    right_clipped_reads += right_clipped;
    both_clipped_reads += both_clipped;
  }

  std::map<std::string, int> results;
  results["total"] = total_reads;
  results["left"] = left_clipped_reads;
  results["right"] = right_clipped_reads;
  results["both"] = both_clipped_reads;
  return results;
}

std::tuple<std::map<std::string,int>, std::map<std::string,std::vector<std::string> > > RealignedFragments::CountRealignedReads(const BreakpointList& breakpoint_list, py::kwargs kwargs) {
  pyassert(ref_index_.IsInitialized(), "BWA index for reference allele not initialized");
  pyassert(breakpoint_list.size() == NumAltAlleles(), "Incorrect number breakpoints");

  // Convert the list of breakpoint strings into GenomicRegions
  std::vector<std::array<sl::GenomicRegion, 4> > breakpoints;
  for (int i=0; i<NumAltAlleles(); i++) {
    const auto & allele_breakpoints = breakpoint_list[i];
    breakpoints.push_back({
      BreakpointToGenomicRegion(std::get<0>(allele_breakpoints), RefHeader()),
      BreakpointToGenomicRegion(std::get<1>(allele_breakpoints), RefHeader()), 
      BreakpointToGenomicRegion(std::get<2>(allele_breakpoints), AltHeader(i)),
      BreakpointToGenomicRegion(std::get<3>(allele_breakpoints), AltHeader(i)) 
    });
  }


  RealignedFragment::score_type min_score_delta = 1.;
  if (kwargs && kwargs.contains("min_score_delta")) {
    min_score_delta = py::cast<AlleleOverlap::score_type>(kwargs["min_score_delta"]);
  }

  bool count_straddle = true;
  if (kwargs && kwargs.contains("count_straddle")) {
    count_straddle = py::cast<bool>(kwargs["count_straddle"]);
  }


  int rl_reads = 0, al_reads = 0, rr_reads = 0, ar_reads = 0, amb_reads = 0;
  // Track the read names overlapping each breakpoint
  std::map<std::string, std::vector<std::string> > overlap_read_names {
    { "rl", std::vector<std::string>() },
    { "rr", std::vector<std::string>() },
    { "al", std::vector<std::string>() },
    { "ar", std::vector<std::string>() },
  };


  // Realign all of the fragments, counting the overlaps
  for (const auto & named_fragment : fragments_) {
    const auto & read_name = named_fragment.first;
    const auto & fragment = named_fragment.second;
    
    // Realign the fragment to the reference allele
    RealignedFragment ref_realignment(fragment, ref_index_, insert_size_dist_);
    auto total_log_prob = ref_realignment.TotalLogProb();

    std::vector<RealignedFragment> alt_realignments;
    for (int i=0; i<NumAltAlleles(); i++) {
      // Realign the fragment to this alternate allele
      alt_realignments.emplace_back(fragment, alt_indexes_[i], insert_size_dist_);
      total_log_prob = LogSumPow(total_log_prob, alt_realignments.back().TotalLogProb()); 
    }

    // Find all overlaps
    AlleleOverlap max_ref_overlap, max_alt_overlap;
    for (int i=0; i<NumAltAlleles(); i++) {
      const auto & bp = breakpoints[i];
      if (ref_realignment.HasBestPair())
        max_ref_overlap = std::max(
          max_ref_overlap,
          AlleleOverlap(ref_realignment.BestPair(), i, bp[0], bp[1], BreakpointOverlap::Allele::Ref, total_log_prob, count_straddle)
        );
      if (alt_realignments[i].HasBestPair())
        max_alt_overlap = std::max(
          max_alt_overlap,
          AlleleOverlap(alt_realignments[i].BestPair(), i, bp[2], bp[3], BreakpointOverlap::Allele::Alt, total_log_prob, count_straddle)
        );
    }

    // if (read_name == "alt-130") {
    //   std::cerr << fragment << std::endl;
    //   if (ref_realignment.HasBestPair()) {
    //     std::cerr << "Read realignemnt:" << std::endl;
    //     std::cerr << ref_realignment.BestPair() << std::endl;
    //   }
    //   if (alt_realignments[0].HasBestPair()) {
    //     std::cerr << "Alt realignemnt:" << std::endl;
    //     std::cerr << alt_realignments[0].BestPair() << std::endl;
    //   }
    //   std::cerr << max_alt_overlap.MaxQualityScore() << " " << max_ref_overlap.MaxQualityScore() << std::endl;
    // }

    if (!max_ref_overlap.HasOverlap() && !max_alt_overlap.HasOverlap())
      continue; // No overlaps to count among any of the alleles

    auto delta = max_alt_overlap.MaxQualityScore() - max_ref_overlap.MaxQualityScore();
    if (max_alt_overlap.HasOverlap() && delta > min_score_delta) {
      // Alternate alignment
      if (max_alt_overlap.HasLeftOverlap()) {
        al_reads += 1;
        overlap_read_names["al"].emplace_back(read_name);
      }
      if (max_alt_overlap.HasRightOverlap()) {
        ar_reads += 1;
        overlap_read_names["ar"].emplace_back(read_name);
      }
    } else if (max_ref_overlap.HasOverlap() && delta < -min_score_delta) {
      // Reference alignment
      if (max_ref_overlap.HasLeftOverlap()) {
        rl_reads += 1;
        overlap_read_names["rl"].emplace_back(read_name);
      }
      if (max_ref_overlap.HasRightOverlap()) {
        rr_reads += 1;
        overlap_read_names["rr"].emplace_back(read_name);
      }
    } else {
      // Everything else is an ambiguous overlap
      amb_reads += 1;
    }
  }

  std::map<std::string,int> results;
  results["rl"] = rl_reads;
  results["rr"] = rr_reads;
  results["al"] = al_reads;
  results["ar"] = ar_reads;
  results["amb"] = amb_reads;

  return std::make_tuple(results, overlap_read_names);
}

namespace test {


bool TestAlignmentOverlap(const std::string& sam_path, const std::string& breakpoint, bool count_straddle) {
  InsertSizeDistribution insert_dist(570., 163., InsertSizeDistribution::density_type());
  
  sl::BamReader reader;
  reader.Open(sam_path);

  sl::BamRecord read1, read2;
  pyassert(reader.GetNextRecord(read1), "Missing first read");
  pyassert(reader.GetNextRecord(read2), "Missing second read");
  reader.Close();

  AddDoubleTag(read1, "as", -10.);
  AddDoubleTag(read2, "as", -10.);

  RealignedReadPair read_pair(read1, read2, insert_dist);

  // SeqLib reverses samtools 1-indexing to 0-indexing conversion so the GenomicRegion constructed
  // this way matches the region string exactly and assumes inclusive start and exclusive end. But
  // the AsGenomicRegion method uses the BAM specification of a 0-indexed inclusive start and exclusive end.
  // The GenomicRegion overlap computation assumes 1-indexed inclusive start and end coordinates.

  sl::GenomicRegion breakpoint_region(breakpoint, reader.Header());
  BreakpointOverlap overlap(read_pair, breakpoint_region, BreakpointOverlap::Allele::RefLeft, -10., count_straddle);

  return overlap.HasOverlap();
}

std::vector<double> TestScoreAlignment(const std::string& ref_seq,
                                       const std::string& aln_path) {
  // Open the input BAM/SAM/CRAM
  sl::BamReader reader;
  reader.Open(aln_path);
  std::vector<double> scores;

  sl::BamRecord read;
  while (reader.GetNextRecord(read)) {
    auto log_prob = ScoreAlignment(read.Sequence(), read.Qualities(0), ref_seq, read);
    scores.push_back(log_prob);
  }

  reader.Close();
  return scores;
}

bool TestStraddle(const std::string& sam_path, const std::string& left_region,
                  const std::string& right_region, int min_overlap, bool strict) {
  sl::BamReader reader;
  reader.Open(sam_path);

  sl::BamRecord read1, read2;
  pyassert(reader.GetNextRecord(read1), "Missing first read");
  pyassert(reader.GetNextRecord(read2), "Missing second read");
  reader.Close();

  AlignedFragment fragment(read1);
  fragment.SetRead(read2);

  sl::GenomicRegion left(left_region, reader.Header()), right(right_region, reader.Header());
  return strict ? fragment.StrictStraddles(left, right, min_overlap)
                : fragment.Straddles(left, right, min_overlap);
}

}  // namespace test

}  // namespace npsv