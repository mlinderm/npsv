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

class RealignedFragment {
 public:
  RealignedFragment(const sl::BamRecord& read);

  bool HasFirst() const { return !first_.isEmpty(); }
  bool HasSecond() const { return !second_.isEmpty(); }
  bool IsProperPair() const;
  int InsertSize() const;

  void SetRead(const sl::BamRecord& read);
  
  bool Straddles(const sl::GenomicRegion& left_region, const sl::GenomicRegion& right_region, int min_overlap) const;
  double ProbMapQ() const;

  friend std::ostream& operator<<(std::ostream&, const RealignedFragment&);
 private:
  sl::BamRecord first_;
  sl::BamRecord second_;
  sl::BamRecord* left_;
};

class RealignedFragments {
  typedef std::unordered_map<std::string, RealignedFragment> FragmentsMap;

 public:
  RealignedFragments(
      double insert_size_mean, double insert_size_std,
      const InsertSizeDistribution::density_type& insert_size_density,
      const std::string& bam_path);

  FragmentsMap::size_type size() const { return fragments_.size(); }

  int GatherReads(const std::string& region, int max_reads);

  std::map<std::string,double> CountBaselineStraddlers(const std::string& event, int flank,
    int alt_size_delta, double z_threshold, int min_overlap=1) const;

 private:
  InsertSizeDistribution insert_size_dist_;
  sl::BamReader reader_;
  FragmentsMap fragments_;

  void AddRead(const sl::BamRecord& read);

  std::tuple<double, double, int> FragmentInsertProbabilities(const RealignedFragment& fragment, int alt_size_delta) const;
  double FragmentProbConcordance(const RealignedFragment& fragment, int alt_size_delta) const;
  static double FragmentProbConcordance(double ref_prob, double alt_prob);

};

}  // namespace npsv