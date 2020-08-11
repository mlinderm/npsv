#include "realigner.hpp"

#include <stdexcept>
#include <cmath>

namespace {
void assert_throw(const bool cond, const std::string& text,
                  const std::string& file, const int line) {
  if (!cond) {
    throw std::runtime_error(text + ". In file: " + file +
                             " on line: " + std::to_string(line));
  }
}

#define pyassert(cond, text) assert_throw(cond, text, __FILE__, __LINE__)

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
RealignedFragment::RealignedFragment(const sl::BamRecord& read) {
  if (read.FirstFlag()) {
    first_ = read;
    left_ = &first_;
  } else {
    second_ = read;
    left_ = &second_;
  }
}

bool RealignedFragment::IsProperPair() const {
  return (HasFirst() && first_.ProperPair()) ||
         (HasSecond() && second_.ProperPair());
}

int RealignedFragment::InsertSize() const {
  return left_->FullInsertSize();
}

void RealignedFragment::SetRead(const sl::BamRecord& read) {
  if (read.FirstFlag()) {
    first_ = read;
  } else {
    second_ = read;
  }
  left_ = first_.InsertSize() > 0 ? &first_ : &second_;
  pyassert(left_->Position() <= left_->MatePosition(), "'left' read not in correct orientation");
}

bool RealignedFragment::Straddles(const sl::GenomicRegion& left_region,
                                  const sl::GenomicRegion& right_region,
                                  int min_overlap) const {
  // This is an expansive definition of straddling, a more strict definition, could require the start
  // coordinates for the reads to be the respective regions
  return ReadOverlapRegion(*left_, left_region) >= min_overlap &&
         MateOverlapRegion(*left_, right_region) >= min_overlap;
}

double RealignedFragment::ProbMapQ() const {
  double result = 1.;
  if (HasFirst())
    result *= (1. - std::pow(10., -first_.MapQuality() / 10.));
  if (HasSecond())
    result *= (1. - std::pow(10., -second_.MapQuality() / 10.));
  return result;
}


std::ostream& operator<<(std::ostream& os, const RealignedFragment& fragment) {
  return (os << fragment.first_ << std::endl << fragment.second_);
}


RealignedFragments::RealignedFragments(
    double insert_size_mean, double insert_size_std,
    const InsertSizeDistribution::density_type& insert_size_density,
    const std::string& bam_path)
    : insert_size_dist_(insert_size_mean, insert_size_std,
                        insert_size_density) {
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

std::tuple<double, double, int> RealignedFragments::FragmentInsertProbabilities(const RealignedFragment& fragment, int alt_size_delta) const {
  int ref_insert_size = fragment.InsertSize();
  double ref_prob = insert_size_dist_(ref_insert_size);
      
  int alt_insert_size = ref_insert_size + alt_size_delta;
  double alt_prob = alt_insert_size > 0 ? insert_size_dist_(alt_insert_size) : 0.;

  return std::make_tuple(ref_prob, alt_prob, ref_insert_size);  
}

double RealignedFragments::FragmentProbConcordance(const RealignedFragment& fragment, int alt_size_delta) const {
  double ref_prob, alt_prob;
  std::tie(ref_prob, alt_prob, std::ignore) = FragmentInsertProbabilities(fragment, alt_size_delta);
  return FragmentProbConcordance(ref_prob, alt_prob);
}

double RealignedFragments::FragmentProbConcordance(double ref_prob, double alt_prob) {
  return CONC_PRIOR * ref_prob / (CONC_PRIOR * ref_prob+ DISC_PRIOR * alt_prob);
}


std::map<std::string,double> RealignedFragments::CountBaselineStraddlers(
    const std::string& event, int flank,
    int alt_size_delta, double z_threshold, int min_overlap) const {
  
  // SeqLib GenomicRegions are 1-indexed with inclusive start and end coordinates
  sl::GenomicRegion event_region(event, reader_.Header());
  
  sl::GenomicRegion left_region(event_region.chr, event_region.pos1-flank, event_region.pos1-1);
  sl::GenomicRegion right_region(event_region.chr, event_region.pos2+1, event_region.pos2+flank);
  
  // Anchor regions within the event (on left and right side)
  sl::GenomicRegion left_event_region(event_region.chr, event_region.pos1, std::min(event_region.pos2, event_region.pos1 + flank - 1));
  sl::GenomicRegion right_event_region(event_region.chr, std::max(event_region.pos1, event_region.pos2 - flank + 1), event_region.pos2);

  int insert_count = 0, insert_upper = 0, insert_lower = 0;
  double ref_weighted_count = 0., alt_weighted_count = 0., ref_conc_count = 0., alt_conc_count = 0.;

  for (auto& named_fragment : fragments_) {
    auto& fragment = named_fragment.second;
    if (!fragment.IsProperPair()) continue;

    // Read straddles the entire event
    if (fragment.Straddles(left_region, right_region, min_overlap)) {
      insert_count += 1;
      
      int ref_insert_size;
      double ref_prob, alt_prob;
      std::tie(ref_prob, alt_prob, ref_insert_size) = FragmentInsertProbabilities(fragment, alt_size_delta);

      if (ref_prob + alt_prob > 0.) {
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

    // Read straddles one of the breakpoints (shouldn't straddle both). This is an expansive
    // definition that doesn't require the "right" read to start in the "right" region. And thus
    // counts reads that span the breakpoint as straddlers.
    bool ref_straddle_left = fragment.Straddles(left_region, left_event_region, min_overlap);
    bool ref_straddle_right = fragment.Straddles(right_event_region, right_region, min_overlap);
    if (ref_straddle_left || ref_straddle_right) {
      pyassert(ref_straddle_left ^ ref_straddle_right, "Reads straddling both breakpoints should have been filtered out");
      ref_weighted_count += 0.5;  // p_alt is by definition 0
      
      double p_conc = FragmentProbConcordance(fragment, alt_size_delta);
      if (p_conc > 0.5) {
        ref_conc_count += fragment.ProbMapQ();
      }
    }

  }

  std::map<std::string, double> results;
  results["insert_count"] = insert_count;
  results["insert_lower"] = insert_lower;
  results["insert_upper"] = insert_upper;
  results["ref_weighted_count"] = ref_weighted_count;
  results["alt_weighted_count"] = alt_weighted_count;
  results["ref_conc_count"] = ref_conc_count;
  results["alt_conc_count"] = alt_conc_count;
  return results;
}

}  // namespace npsv