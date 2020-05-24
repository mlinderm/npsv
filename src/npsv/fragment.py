class Fragment(object):
    """Single read fragment

    Attributes:
        query_name: Read name
        left: Left read
        right: Right read
    """
    def __init__(self, read):
        assert Fragment.is_primary(read)
        self.query_name = read.query_name
        if not read.is_paired or read.template_length > 0:
            self.left = read
            self.right = None
        else:
            self.left = None
            self.right = read

    @staticmethod
    def is_primary(read):
        return not read.is_supplementary and not read.is_secondary

    def is_properly_paired(self):
        return (
            self.left is not None
            and self.right is not None
            and self.left.is_paired
            and self.left.is_reverse != self.right.is_reverse
        )

    def is_pair_straddle(self, left_id, left_span, right_id, right_span, min_aligned):
        if not self.is_properly_paired():
            return False

        if self.left.reference_id != left_id or self.right.reference_id != right_id:
            return False

        return (
            self.left.get_overlap(*left_span) >= min_aligned
            and self.right.get_overlap(*right_span) >= min_aligned
        )

    def get_potential_split_reads_start_anchor(self, contig_id, pos):
        left_potential = (
            self.left
            and self.left.reference_id == contig_id
            and self.left.reference_start
            < pos
            <= self.left.reference_start + self.left.query_length
        )
        right_potential = (
            self.right
            and self.right.reference_id == contig_id
            and self.right.reference_start
            < pos
            <= self.right.reference_start + self.right.query_length
        )

        if left_potential and right_potential:
            return self.left, self.right
        elif left_potential:
            return (self.left,)
        elif right_potential:
            return (self.right,)
        else:
            return ()

    def get_potential_split_reads_end_anchor(self, contig_id, pos):
        left_potential = (
            self.left
            and self.left.reference_id == contig_id
            and self.left.reference_end - self.left.query_length
            < pos
            <= self.left.reference_end
        )
        right_potential = (
            self.right
            and self.right.reference_id == contig_id
            and self.right.reference_end - self.right.query_length
            < pos
            <= self.right.reference_end
        )

        if left_potential and right_potential:
            return self.left, self.right
        elif left_potential:
            return (self.left,)
        elif right_potential:
            return (self.right,)
        else:
            return ()

    def prob_fragment_length(self, sample, event_length=0):
        fragment_length = self.left.template_length
        if event_length < fragment_length:
            # TODO: Use formula from manta?
            # Should this be the cdf (prob of this insert size of longer?)
            return sample.generic_lib.dens[fragment_length - event_length]
        else:
            return 0.0

    def zscore_fragment_length(self, sample, event_length=0):
        fragment_length = self.left.template_length
        lib = sample.generic_lib
        return (fragment_length - lib.mean) / lib.sd

    def add_read(self, read):
        # TODO: Handle reads that overlap each other
        assert Fragment.is_primary(read)
        assert read.is_paired
        if self.left is not None and self.right is None:
            if read.template_length <= 0:
                self.right = read
            elif read.__hash__() == self.left.__hash__():
                return  # Duplicate read
        elif self.left is None and self.right is not None:
            if read.template_length > 0:
                self.left = read
            elif read.__hash__() == self.right.__hash__():
                return  # Duplicate read
        assert self.left.__hash__() != self.right.__hash__()
        # May not be true if "overlapping" fragment
        # assert(self.left.reference_start <= self.right.reference_start)


class SpanningFragments(object):
    """Read fragments
    
    Attributes:
        fragments: Dictionary of a (query_name, Fragement)
    """
    def __init__(self):
        self.fragments = {}

    def __iter__(self):
        return iter(self.fragments.values())

    def add_read(self, read):
        name = read.query_name
        if name in self.fragments:
            self.fragments[name].add_read(read)
        else:
            self.fragments[name] = Fragment(read)


def gather_reads(fragments: SpanningFragments, bam, contig: str, span, max_reads: int=None):
    """Populate SpanningFragments with reads from region
    
    Args:
        fragments (SpanningFragments): Fragments in BAM file
        bam ([type]): Open PySAM.AlignmentFile
        contig (str): Contig
        span ([type]): Tuple of 0-indexed [pos, end) coordinates
        max_reads (int, optional): Maximum number of reads to extract. Defaults to None.
    
    Returns:
        bool: True if hit max_reads limit
    """
    for i, read in enumerate(bam.fetch(contig, *span)):
        # Only construct Fragments out of mapped, non-duplicate, primary reads
        if read.is_unmapped or read.is_duplicate or not Fragment.is_primary(read):
            continue

        if max_reads is not None and i > max_reads:
            return True

        fragments.add_read(read)

    return False

