#pragma once
#include <pybind11/pybind11.h>
#include <functional>
#include <iosfwd>
#include <string>
#include <unordered_map>
#include <tuple>
#include "SeqLib/BamReader.h"

#include "aligner.hpp"

namespace py = pybind11;
namespace sl = SeqLib;

namespace npsv {

const double CONC_PRIOR = 0.95;
const double DISC_PRIOR = 1 - CONC_PRIOR;

class IndexedSequence;

class AlignedFragment {
 public:
  AlignedFragment(const sl::BamRecord& read);

  bool HasFirst() const { return !first_.isEmpty(); }
  bool HasSecond() const { return !second_.isEmpty(); }
  bool IsProperPair() const;
  int InsertSize() const;

  const sl::BamRecord& FirstRead() const { return first_; }
  const sl::BamRecord& SecondRead() const { return second_; }

  int32_t RightPos1() const;
  int32_t LeftPos2() const;

  void SetRead(const sl::BamRecord& read);
  
  bool Straddles(const sl::GenomicRegion& left_region, const sl::GenomicRegion& right_region, int min_overlap) const;
  double ProbMapQ() const;

  friend std::ostream& operator<<(std::ostream&, const AlignedFragment&);
 private:
  sl::BamRecord first_;
  sl::BamRecord second_;
  sl::BamRecord* left_;
};

class RealignedReadPair {
 public:
  typedef double score_type;

  RealignedReadPair() : left_(nullptr), right_(nullptr), score_(0) {}
  RealignedReadPair(const sl::BamRecord& first);
  RealignedReadPair(const sl::BamRecord& first, const sl::BamRecord& second_,
                    const InsertSizeDistribution&);

  bool IsValid() const{ return left_ || right_; }
  score_type Score() const { return score_; }

  bool operator<(const RealignedReadPair& other) const {
    return score_ < other.score_;
  }
  bool operator>(const RealignedReadPair& other) const {
    return score_ > other.score_;
  }

  sl::GenomicRegion LeftReadRegion() const;
  sl::GenomicRegion RightReadRegion() const;

  friend std::ostream& operator<<(std::ostream&, const RealignedReadPair&);

 private:
  const sl::BamRecord* left_;
  const sl::BamRecord* right_;
  score_type score_;

  bool Concordant() const;
  int32_t InsertSize() const;
};

class BreakpointOverlap {
 public:
  enum class Allele {
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

  enum class Overlap {
    None = 0,
    FragmentContains = 1,
    ReadContains = 3,
  };

  typedef RealignedReadPair::score_type score_type;

  BreakpointOverlap()
      : breakpoint_(nullptr),
        allele_(Allele::None),
        kind_(Overlap::None),
        quality_score_(0) {}
  BreakpointOverlap(const RealignedReadPair& fragment, const sl::GenomicRegion& breakpoint,
          Allele allele, score_type total_log_prob,
          bool count_straddle);

  bool HasOverlap() const { return kind_ != Overlap::None; }

  score_type QualityScore() const { return quality_score_; }

 private:
  const sl::GenomicRegion* breakpoint_;
  Allele allele_;
  Overlap kind_;
  score_type quality_score_;  // Phred-scale quality score
};

inline BreakpointOverlap::Allele operator&(
    const BreakpointOverlap::Allele& lhs,
    const BreakpointOverlap::Allele& rhs) {
  return static_cast<BreakpointOverlap::Allele>(
      static_cast<std::underlying_type_t<BreakpointOverlap::Allele> >(lhs) &
      static_cast<std::underlying_type_t<BreakpointOverlap::Allele> >(rhs));
}

class AlleleOverlap {
 public:
  typedef BreakpointOverlap::score_type score_type;

  AlleleOverlap() : read_pair_(nullptr), allele_index_(-1) {}

  AlleleOverlap(const RealignedReadPair& read_pair,
                int allele_index,
                const sl::GenomicRegion& left_breakpoint,
                const sl::GenomicRegion& right_breakpoint,             
                BreakpointOverlap::Allele allele,
                score_type total_log_prob, bool count_straddle);

  bool HasOverlap() const { return left_overlap_.HasOverlap() || right_overlap_.HasOverlap(); }
  bool HasLeftOverlap() const { return left_overlap_.HasOverlap(); }
  bool HasRightOverlap() const { return right_overlap_.HasOverlap(); }

  score_type MaxQualityScore() const { return std::max(left_overlap_.QualityScore(), right_overlap_.QualityScore()); }

  bool operator<(const AlleleOverlap& rhs) const { 
    return MaxQualityScore() < rhs.MaxQualityScore();
  }

  bool operator>(const AlleleOverlap& rhs) const { 
    return MaxQualityScore() > rhs.MaxQualityScore();
  }

 private:
  const RealignedReadPair* read_pair_;
  int allele_index_;
  BreakpointOverlap left_overlap_;
  BreakpointOverlap right_overlap_;
};

class RealignedFragment {
  typedef std::vector<RealignedReadPair> PairSequence;
 public:
  typedef RealignedReadPair::score_type score_type;

  RealignedFragment(const AlignedFragment&, const IndexedSequence&,
                    const InsertSizeDistribution&);

  PairSequence::size_type NumAlignments() const {
    return read_pairs_.size();
  }

  bool HasBestPair() const { return !read_pairs_.empty(); }
  const RealignedReadPair& BestPair() const { return read_pairs_.front(); }

  score_type TotalLogProb() const { return total_log_prob_; }

 private:
  const AlignedFragment& original_alignment_;
  sl::BamRecordVector first_alignments_;
  sl::BamRecordVector second_alignments_;
  PairSequence read_pairs_;
  score_type total_log_prob_;
};

class IndexedSequence {
  public:
  IndexedSequence() {}
  IndexedSequence(const sl::UnalignedSequence& sequence);

  bool IsInitialized() const  { return !bwa_.IsEmpty(); }
  void Initialize(const sl::UnalignedSequence& sequence);

  const std::string& Sequence() const { return sequence_.Seq; }
  sl::BamHeader Header() const { return bwa_.HeaderFromIndex(); }

  void AlignSequence(const sl::BamRecord& read, sl::BamRecordVector& alignments) const;

  private:
    sl::UnalignedSequence sequence_;
    sl::BWAWrapper bwa_;
};


class RealignedFragments {
  typedef std::unordered_map<std::string, AlignedFragment> FragmentsMap;
  typedef std::vector<IndexedSequence> AltIndexesSequence;
  typedef std::vector<std::tuple<std::string, std::string, std::string, std::string>> BreakpointList;

 public:
  RealignedFragments(
      const std::string& fasta_path,
      double insert_size_mean, double insert_size_std,
      const InsertSizeDistribution::density_type& insert_size_density,
      const std::string& bam_path);

  AltIndexesSequence::size_type NumAltAlleles() const { return alt_indexes_.size(); }

  FragmentsMap::size_type size() const { return fragments_.size(); }

  int GatherReads(const std::string& region, int max_reads);

  std::map<std::string, double> CountPipelineStraddlers(
      const std::string& left_breakpoint, const std::string& right_breakpoint,
      int flank, int alt_size_delta, double z_threshold = 1.5,
      int min_overlap = 1) const;

  std::tuple<std::map<std::string,int>, std::map<std::string,std::vector<std::string> > > CountRealignedReads(const BreakpointList& breakpoints, py::kwargs kwargs);

 private:
  InsertSizeDistribution insert_size_dist_;
  
  IndexedSequence ref_index_;
  std::vector<IndexedSequence> alt_indexes_;
  
  sl::BamReader reader_;
  FragmentsMap fragments_;

  sl::BamHeader RefHeader() const { return ref_index_.Header(); }
  sl::BamHeader AltHeader(int index) const { return alt_indexes_[index].Header(); }

  void AddRead(const sl::BamRecord& read);

  std::tuple<double, double, int> FragmentInsertProbabilities(const AlignedFragment& fragment, int alt_size_delta) const;
  double FragmentProbConcordance(const AlignedFragment& fragment, int alt_size_delta) const;
  static double FragmentProbConcordance(double ref_prob, double alt_prob);
};

namespace test {
  std::vector<double> TestScoreAlignment(const std::string& ref_seq,
                                       const std::string& aln_path);
};


}  // namespace npsv