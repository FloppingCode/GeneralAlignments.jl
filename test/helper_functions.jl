function checkAlignedSequenceReadingFrame(aligned_seq,readingFrame)
    # confirms that all triplet indels obey the reference readingFrame
    
    i = 1
    # check for A_ref_alignment
    while i+2 <= length(aligned_seq)
            if aligned_seq[i] == DNA_Gap && aligned_seq[i+1] == DNA_Gap && aligned_seq[i+2] == DNA_Gap 
                if (i-1) % 3 == readingFrame # zeroindex
                    i += 3
                else
                    return false
                end
            else
                i += 1
            end
    end
    return true
end

function checkRef2NoisyAlignmentRespectsCodonBoundaries(A_ref_aligned,B_noisy_aligned,readingFrame)
    """
    Inevitably the readingFrame will be ruined by non-triplet indels. But the alignment will still be respecting the codonBoundaries of A.
    To check that this is the case we have to consider how gaps inserted shift the subsequent codon boundaries.
    """
    gapsInA = 0
    i = 1
    while i+2 <= length(A_ref_aligned)
        
        # check for triplet indel
        if (A_ref_aligned[i] == DNA_Gap && A_ref_aligned[i+1] == DNA_Gap && A_ref_aligned[i+2] == DNA_Gap) ||
            (B_noisy_aligned[i] == DNA_Gap && B_noisy_aligned[i+1] == DNA_Gap && B_noisy_aligned[i+2] == DNA_Gap)
            if (i-1) % 3 == (readingFrame+gapsInA) % 3 # zeroindex
                i += 3
            else
                println("index",i)
                return false
            end
        else
            if A_ref_aligned[i] == DNA_Gap
                gapsInA += 1
            end
            i += 1
        end
        
    end
    return true
end