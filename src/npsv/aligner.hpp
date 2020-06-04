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

class InsertSizeDistribution {
 public:
  InsertSizeDistribution(double mean, double std) : mean_(mean), std_(std) {}

  double operator()(int32_t insert_size) const {
    // https://stackoverflow.com/a/10848293 
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (insert_size - mean_) / std_;
    return inv_sqrt_2pi / std_ * std::exp(-0.5f * a * a);
  }

 private:
  double mean_, std_;
};

class Fragment;
class AlleleAlignments;
class AlleleReference;

class AlignmentPair {
 public:
  typedef double score_type;

  AlignmentPair() : first_(nullptr), second_(nullptr), score_(0) {}
  AlignmentPair(const sl::BamRecord& first, const sl::BamRecord& second_,
                const InsertSizeDistribution&);

  bool Valid() const { return first_ && second_; }
  bool Concordant() const;
  int32_t InsertSize() const;

  bool ReadsContain(const sl::GenomicRegion& region) const;

  bool operator<(const AlignmentPair& other) const {
    return score_ < other.score_;
  }
  friend std::ostream& operator<<(std::ostream&, const AlignmentPair&);

 private:
  const sl::BamRecord* first_;
  const sl::BamRecord* second_;
  score_type score_;

  friend class Fragment;
};

class Fragment {
 public:
  typedef AlignmentPair::score_type score_type;

  Fragment(const sl::BamRecord& read) : first_(read) {}

  void SetSecond(const sl::BamRecord& read);

  bool HasBestPair() const { return best_pair_.Valid(); }
  void SetBestPair(const InsertSizeDistribution&);
  bool BestAlignmentContains(const sl::GenomicRegion& region) const {
    return best_pair_.ReadsContain(region);
  };
  score_type BestAlignmentScore() const { return best_pair_.score_; }

 private:
  sl::BamRecord first_;
  sl::BamRecordVector first_alignments_;
  sl::BamRecord second_;
  sl::BamRecordVector second_alignments_;
  AlignmentPair best_pair_;

  friend class AlleleAlignments;
  friend class AlleleReference;
};

class AlleleAlignments {
  typedef std::map<std::string, Fragment> fragments_type;

 public:
  typedef fragments_type::size_type size_type;
  typedef fragments_type::iterator iterator;

  bool IsInitialized() const { return !bwa_.IsEmpty(); }
  void Initialize(const sl::UnalignedSequence& sequence);

  sl::BamHeader Header() const { return bwa_.HeaderFromIndex(); }

  iterator begin() { return fragments_.begin(); }
  iterator end() { return fragments_.end(); }

  void Clear() { fragments_.clear(); }
  size_type Size() const { return fragments_.size(); }

  void Align(const sl::BamRecord& read,
             const InsertSizeDistribution& insert_size_dist);

 private:
  sl::BWAWrapper bwa_;
  fragments_type fragments_;
};

class AlleleReference {
 public:
  AlleleReference(const std::string& fasta_path, double insert_size_mean,
                  double insert_size_std);

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

}  // namespace npsv