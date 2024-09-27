using BioSequences
using Distributions

function mutateSequence(seq::LongDNA{4})
    n = length(seq)
    dna = LongDNA{4}("ACGT")
    num_indels = rand(Poisson(3))
    for i in 1:num_indels
        indel_pos = rand(1:n)
        if rand(1:2) == 1
            # insertion
            println("insertion!!!")
            seq = seq[1:indel_pos] * randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 1) * seq[indel_pos + 1 : end]
        else
            #deletion
            println("del")
            deleteat!(seq, indel_pos)
        end
    end
    return seq
end