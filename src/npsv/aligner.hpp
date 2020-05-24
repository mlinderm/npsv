#pragma once

int add(int i, int j);

namespace npsv {
    
enum GenomicRegionOverlap {
  NoOverlap = 0,
  PartialOverlap = 1,
  ContainsArg = 2,
  ContainedInArg = 3
};

enum SAMFlags { SecondaryAlignment = 256, SupplementaryAlignment = 2048 };

void align();

}  // namespace npsv