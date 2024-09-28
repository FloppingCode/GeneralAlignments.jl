using NextGenSeqUtils
using BioSequences

include("MSA_noising.jl")
include("seed_chain_align.jl")
include("needleman_wunsch.jl")
# collect sequnces and ungap sequence 
names, refAndseqs = NextGenSeqUtils.read_fasta_with_names("src/P018_subset.fasta")
numOfSequence = length(refAndseqs)-1
#numOfSequence = 5

nameArray = [String(names[i]) for i in 1:length(names)]
#println(length(nameArray))
ref = LongDNA{4}(refAndseqs[1])
seqs = Array{LongDNA{4}}(undef,numOfSequence)
for iseq in 1:numOfSequence
   seqs[iseq] = LongDNA{4}(refAndseqs[iseq+1])
end

# NOTE that there is at least one single indel which might screw reading frame if mutate
# add some noise into sequences (not reference)
#for i in 1:length(seqs)-1
#    seqs[i] = mutateSequence(seqs[i])
#end

# align pairwise via seeding
function msa_ref_alignment(ref, seqs)
    alignment = Array{LongDNA{4}}(undef,length(seqs)+1)
    alignment[1] = ref
    for seqId in 1:length(seqs)
        # seed pairwise alignment with codon respecting triplet moves
        match_moves = [Move(1,.0), Move(3,.0)]
        # important example
        hor_moves =  [Move(1, 2.4, 1, 0, 1,0, false), Move(3, 2.0, 1,0,3,0,true)]
        vert_moves = [Move(1, 2.4, 1, 0, 1,0, false), Move(3, 2.0, 1,0,3,0,true)]
        aligned_ref, aligned_seq = seed_chain_align(ref, ungap(seqs[seqId]), 0.0, 0.5, match_moves, vert_moves, hor_moves, 0.98*(2.0/3), 21)

        # fix readingframes
        println(seqId)
        fixed_aligned_seq = fix_alignment_readingframe(aligned_ref,aligned_seq)
        alignment[seqId+1] = fixed_aligned_seq
    end
    return alignment
end


# we find all single indels and order them and then do the thing
function fix_alignment_readingframe(aligned_ref,aligned_seq)
    # find all single gaps and get their index. Once we have that we insert N in order unsure of how to fix...
    indelDict = Dict()
    N_codon = LongDNA{4}("AAC")
    NumCodons = (length(aligned_seq) ÷ 3)
    # loop through to find all single indels
    for i in 2:length(aligned_seq)
        #find deletion indicies
        if (i <= length(aligned_seq)-4) && (aligned_seq[i] == DNA_Gap) && !(aligned_seq[i+1] == DNA_Gap || aligned_seq[i-1] == DNA_Gap)
            indelDict[i] = -1
        end
        # find insertion indicies
        if (i <= length(aligned_ref)-4) && (aligned_ref[i] == DNA_Gap) && !(aligned_ref[i+1] == DNA_Gap || aligned_ref[i-1] == DNA_Gap)
            indelDict[i] = 1
        end
    end
    indelIndicies = collect(keys(indelDict))
    indelIndicies = sort(indelIndicies)

    insertAddon = 0
    for x in indelIndicies
        # deletion, we have to identify current codonboundaries in order to replace the codon. 
        # NOTE We assume the readingFrame is 0 mod 3
        if indelDict[x] == -1
            r = (x+insertAddon-1)%3
            startInsertPos = insertAddon+x-r
            aligned_seq = aligned_seq[1:startInsertPos-1] * N_codon * aligned_seq[startInsertPos+3:end]

        # insertion
        else
            aligned_seq = aligned_seq[1:x-1-insertAddon] * aligned_seq[x+1-insertAddon:end]
            insertAddon -= 1
        end
    end

    return aligned_seq
end

msa_alignment = msa_ref_alignment(ref, seqs)
#convert msa_alignment into strings
msa_alignment_str = Array{String}(undef,length(msa_alignment))
for i in 1:length(msa_alignment)
    msa_alignment_str[i] = string(msa_alignment[i])
end

NextGenSeqUtils.write_fasta("src/reconstruction-alignment-no-noise.fasta", msa_alignment_str, names =nameArray)