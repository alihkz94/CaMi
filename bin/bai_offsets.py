#!/usr/bin/env python3
# =============================================================================
# bai_offsets.py
# -----------------------------------------------------------------------------
# PURPOSE:
#   Parse a BAM index (.bai) and report two byte offsets in the companion BAM:
#     * header_end    = smallest BGZF block-start (coffset) of any alignment
#                       record (everything before it is the BAM header).
#     * unmapped_start = largest BGZF block-start (coffset) among all indexed
#                       alignment chunks = where the placed alignments end and
#                       the fully-unmapped ('*') reads begin.
#   With these we can download ONLY [0:header_end) + [unmapped_start:EOF) of a
#   coordinate-sorted BAM to recover every both-mates-unmapped read, instead of
#   the whole (tens-of-GB) file.
#
#   A BAI virtual offset packs (coffset<<16 | uoffset); coffset is the byte
#   position of the containing BGZF block in the BAM. We only need coffset, so
#   the >>16 is exact at a block boundary (safe to start an aria2c range there).
#
# BAI layout (SAMv1 spec): magic 'BAI\1', n_ref int32, then per ref:
#   n_bin int32; per bin: bin uint32, n_chunk int32, chunks[(u64 beg,u64 end)];
#   n_intv int32; ioffset[n_intv] u64.  Optional trailing n_no_coor u64.
#
# RUN: python bai_offsets.py file.bam.bai   -> prints "header_end unmapped_start"
# =============================================================================
import sys, struct

def parse(path):
    with open(path, "rb") as fh:
        data = fh.read()
    off = 0
    def u(fmt):
        nonlocal off
        n = struct.calcsize(fmt)
        v = struct.unpack_from(fmt, data, off); off += n
        return v
    magic, = u("<4s")
    assert magic == b"BAI\x01", f"not a BAI file: {magic!r}"
    n_ref, = u("<i")
    max_voff = 0         # largest alignment END virtual offset = start of the
                         # fully-unmapped ('*') reads (they follow all mapped data)
    for _ in range(n_ref):
        n_bin, = u("<i")
        for _ in range(n_bin):
            _bin, = u("<I")
            n_chunk, = u("<i")
            for _ in range(n_chunk):
                beg, end = u("<QQ")
                if end > max_voff:
                    max_voff = end
        n_intv, = u("<i")
        for _ in range(n_intv):
            io, = u("<Q")           # linear-index voffsets are chunk-begins; skip
    coffset = max_voff >> 16        # BGZF block byte-start to download from
    uoffset = max_voff & 0xFFFF     # bytes to skip within that block's plaintext
    return max_voff, coffset, uoffset

if __name__ == "__main__":
    v, c, u_ = parse(sys.argv[1])
    # prints: full_virtual_offset  block_byte_offset  within_block_offset
    print(f"{v} {c} {u_}")
