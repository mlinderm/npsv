#pragma once
#include <pybind11/pybind11.h>
#include <functional>
#include <iosfwd>
#include <map>
#include <string>
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/UnalignedSequence.h"

namespace py = pybind11;
namespace sl = SeqLib;

namespace npsv {

enum GenomicRegionOverlap {
  NoOverlap = 0,
  PartialOverlap = 1,
  ContainsArg = 2,
  ContainedInArg = 3
};

enum SAMFlags { SecondaryAlignment = 256, SupplementaryAlignment = 2048 };

double LogProbToPhredQual(double prob,
                          double max_qual = std::numeric_limits<double>::max());

class InsertSizeDistribution {
 public:
  typedef std::map<int, double> density_type;

  InsertSizeDistribution(double mean, double std, const density_type& density)
      : mean_(mean), std_(std), density_(density) {}

  double operator()(int insert_size) const;

 private:
  double mean_, std_;
  density_type density_;
};

class Fragment;
class AlleleAlignments;
class Overlap;
class Realigner;

class AlignmentPair {
 public:
  typedef double score_type;

  AlignmentPair() : first_(nullptr), second_(nullptr), score_(0) {}
  AlignmentPair(const sl::BamRecord& first, const sl::BamRecord& second_,
                const InsertSizeDistribution&);

  bool Valid() const { return first_ && second_; }
  score_type Score() const { return score_; }
  bool Concordant() const;
  int32_t InsertSize() const;

  bool operator<(const AlignmentPair& other) const {
    return score_ < other.score_;
  }
  friend std::ostream& operator<<(std::ostream&, const AlignmentPair&);

 private:
  const sl::BamRecord* first_;
  const sl::BamRecord* second_;
  score_type score_;

  friend class Fragment;
  friend class Overlap;
  friend class Realigner;
};

class Fragment {
 public:
  typedef AlignmentPair::score_type score_type;

  Fragment(const sl::BamRecord& read);

  bool HasFirst() const { return !first_.isEmpty(); }
  bool HasSecond() const { return !second_.isEmpty(); }
  bool IsPaired() const { return HasFirst() && HasSecond(); }

  bool HasBestPair() const { return best_pair_.Valid(); }
  const AlignmentPair& BestPair() const { return best_pair_; }
  score_type BestPairScore() const { return best_pair_.score_; }

 private:
  void SetRead(const sl::BamRecord& read);
  void SetBestPair(const InsertSizeDistribution&);
  sl::GenomicRegion MateQueryRegion() const;

  sl::BamRecord first_;
  sl::BamRecordVector first_alignments_;
  sl::BamRecord second_;
  sl::BamRecordVector second_alignments_;

  AlignmentPair best_pair_;
  score_type total_log_prob_;  // log10(Pr(read))

  friend class AlleleAlignments;
  friend class Overlap;
  friend class Realigner;
};

class AlleleAlignments {
  typedef std::map<std::string, Fragment> fragments_type;

 public:
  typedef Fragment::score_type score_type;
  typedef fragments_type::size_type size_type;
  typedef fragments_type::iterator iterator;

  bool IsInitialized() const { return !bwa_.IsEmpty(); }
  void Initialize(const sl::UnalignedSequence& sequence);

  sl::BamHeader Header() const { return bwa_.HeaderFromIndex(); }

  void Align(const sl::BamRecord& read,
             const InsertSizeDistribution& insert_size_dist);

  static score_type ScoreAlignment(const std::string& read_sequence,
                                   const std::string& base_qualities,
                                   const std::string& ref_sequence,
                                   const sl::BamRecord& alignment);

  void ScoreAlignments(const sl::BamRecord&, sl::BamRecordVector&);

  // Access realigned fragments
  iterator begin() { return fragments_.begin(); }
  iterator end() { return fragments_.end(); }

  void clear() { fragments_.clear(); }
  size_type size() const { return fragments_.size(); }

 private:
  sl::UnalignedSequence sequence_;
  sl::BWAWrapper bwa_;
  fragments_type fragments_;
};

class Overlap {
 public:
  enum class OverlapAllele {
    None = 0,
    RefLeft = 1,
    RefRight = 2,
    AltLeft = 4,
    AltRight = 8,
    Ref = 3,
    Alt = 12,
    Left = 5,
    Right = 10
  };

  enum class OverlapKind {
    None = 0,
    FragmentContains = 1,
    ReadContains = 3,
  };

  typedef Fragment::score_type score_type;

  Overlap()
      : fragment_(nullptr),
        breakpoint_(nullptr),
        allele_(OverlapAllele::None),
        kind_(OverlapKind::None),
        quality_score_(0) {}
  Overlap(const Fragment& fragment, const sl::GenomicRegion& breakpoint,
          OverlapAllele allele, score_type total_log_prob,
          bool count_straddle = true);

  bool operator()() const { return fragment_ && breakpoint_; }

  bool IsRef() const {
    return static_cast<unsigned int>(allele_) &
           static_cast<unsigned int>(OverlapAllele::Ref);
  }
  bool IsAlt() const {
    return static_cast<unsigned int>(allele_) &
           static_cast<unsigned int>(OverlapAllele::Alt);
  }

  bool IsLeft() const {
    return static_cast<unsigned int>(allele_) &
           static_cast<unsigned int>(OverlapAllele::Left);
  }
  bool IsRight() const {
    return static_cast<unsigned int>(allele_) &
           static_cast<unsigned int>(OverlapAllele::Right);
  }

  bool HasOverlap() const { return kind_ != OverlapKind::None; }

  std::string FormattedBreakpoint() const;
  score_type QualityScore() const { return quality_score_; }

  void PushTaggedReads(sl::BamRecordVector&) const;

  bool operator<(const Overlap& rhs) const {
    return quality_score_ < rhs.quality_score_;
  }
  bool operator>(const Overlap& rhs) const {
    return quality_score_ > rhs.quality_score_;
  }

 private:
  const Fragment* fragment_;
  const sl::GenomicRegion* breakpoint_;
  OverlapAllele allele_;
  OverlapKind kind_;
  Fragment::score_type quality_score_;  // Phred-scale quality score

  friend class Realigner;
};

class Realigner {
 public:
  Realigner(const std::string& fasta_path, double insert_size_mean,
            double insert_size_std,
            const InsertSizeDistribution::density_type& insert_size_density);

  sl::BamHeader RefHeader() const { return ref_.Header(); }
  sl::BamHeader AltHeader() const { return alt_.Header(); }

  void Clear();
  void Realign(const sl::BamRecord& read);

  py::dict CountAlignments(const std::string& bam_path,
                           const std::string& rl_region,
                           const std::string& al_region, py::kwargs kwargs);

 private:
  InsertSizeDistribution insert_size_dist_;
  AlleleAlignments ref_, alt_;
};

namespace test {
std::vector<AlleleAlignments::score_type> TestScoreAlignment(
    const std::string& ref_sequence, const std::string& alignment_path);
}

}  // namespace npsv