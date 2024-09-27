using Test
using Random
using BioSequences

# files to test
include("../src/needleman_wunsch.jl")
include("../src/seed_chain_align.jl")
#help functions
include("helper_functions.jl")
include("../test_case_creation_HIV/sequence_generator.jl")

# test helper function 
#@testset "test function checkAlignedSequenceReadingFrame" begin
#    # tests function for checking that the alignment respects the readingFrame
##    A = LongDNA{4}("AT---TGGCTACCAGA")
 #   @test checkAlignedSequenceReadingFrame(A,0) == false
 #   B = LongDNA{4}("AT---TGGCTACCAGACGA")
 #   @test checkAlignedSequenceReadingFrame(B,0) == false
 #   C = LongDNA{4}("ATT---GCTACCAGACGA")
 #   @test checkAlignedSequenceReadingFrame(C,0) == true
 #   B = LongDNA{4}("ATTTTT---GCT-ACCAGACGA")
 #   @test checkAlignedSequenceReadingFrame(B,0) == false
  #  B = LongDNA{4}("ATTTTG---GCTA-CAG-CGA")
  #  @test checkAlignedSequenceReadingFrame(B,0) == true
  #  B = LongDNA{4}("ATTTTG---GCTA-CAG-CGA")
  #  @test checkAlignedSequenceReadingFrame(B,2) == false
  #  B = LongDNA{4}("AATTTTG---GCTA-CAG-CGA")
  #  @test checkAlignedSequenceReadingFrame(B,1) == true
#end;

# test another helper function

# NW non-affine unittest
@testset "needleman-wunsch non-affine" begin
        @testset "test for general alignment output" begin
            for _ in 1:1
                A, B = generate_seq_pair(120, 0.1, 0.2, 0.01, 0.01, 10)
                # default moveset (note symmetric in horizontal and vertical)
                match_moves = [Move(1, 0.0,1,0), Move(3, 0.0,1,0)]
                gap_moves =   [Move(1, 2.2,1,0), Move(3, 2.0,1,0)]
                alignment = nw_align(A, B, .0, 0.5, match_moves, gap_moves, gap_moves)
                @testset "ungapped Alignment Is Equal To Input sequence & Alignments Have Same Length" begin
                    ungapped_aligned_A = ungap(alignment[1])
                    ungapped_aligned_B = ungap(alignment[2])
                    # test that the sequences are unmodified
                    @test ungapped_aligned_A == A
                    @test ungapped_aligned_B == B
                    # additionally assert that alignments have the same length
                    @test length(alignment[1]) == length(alignment[2])
                end;

                @testset "alignment is symmetric in A and B if vertical and horizontal moves are the same" begin
                    alignment_symmetric = nw_align(B, A, .0, 0.5, match_moves, gap_moves, gap_moves)
                    @test alignment_symmetric[2] == alignment[1]
                    @test alignment_symmetric[1] == alignment[2]
                end;

                @testset "alignment Is Reverse Compliment Symmetric" begin
                    # TODO see if it can be made reverse compliment symmetric
                    # Fails mainly due to not being an affine alignment, more global solutions
                    # If possible it would depend on how the Backtracking is done...


                    A_reverse_complement = complement(LongDNA{4}(reverse(string(A))))
                    B_reverse_complement = complement(LongDNA{4}(reverse(string(B))))
                    alignment_reverse_complement = nw_align(A_reverse_complement, B_reverse_complement, .0, 0.5, match_moves, gap_moves, gap_moves)
                    # test if the alignment are the same
                    println("A")
                    println(alignment_reverse_complement[1])
                    println(complement(LongDNA{4}(reverse(string(alignment[1])))))
                    @test alignment_reverse_complement[1] == complement(LongDNA{4}(reverse(string(alignment[1])))) skip=true    
                    println("B")
                    println(alignment_reverse_complement[2])
                    println(complement(LongDNA{4}(reverse(string(alignment[2])))))
                    @test alignment_reverse_complement[2] == complement(LongDNA{4}(reverse(string(alignment[2])))) skip = true 
                end;
            
            end
        end;

        #@testset "reference informed alignment (checking stride and phase work correctly)" begin
            # case where A has a perserved reading Frame
        #    @testset "ref2noisy" begin
                # try different parameters
        #        refseq, noisySeq = generate_seq_pair(120, 0.02, 0.0, 0.005, 0.005, 10)
                # we configure the reading Frame
        #        readingFrame = 0
        #        println(readingFrame)
                # moves for alignment
        #        match_moves =          [Move(1, 0.0), Move(3, 0.0)]
        #        insertions_vertical =  [Move(3, 2.0, 3, 0, 3, 0,false), Move(1, 1.5, 1, 0)]
        #       deletions_horizontal = [Move(3, 2.0, 1, 0, 3, 0,false), Move(1, 1.5, 1, 0)]
                # check that insertions and deletions are correct 
        #        alignment_ref_2_noisy = nw_align(refseq, noisySeq, .0, 0.5, match_moves, insertions_vertical, deletions_horizontal)
        #        println(alignment_ref_2_noisy[1])
        #        println(alignment_ref_2_noisy[2])
        #        println("ref")
        #        println(refseq)
        #        @test checkRef2NoisyAlignmentRespectsCodonBoundaries(alignment_ref_2_noisy[1],alignment_ref_2_noisy[2],readingFrame) == true
        #    end;
        #end;   
end;


# NW-affine unittest
@testset "needleman-wunsch affine" begin
    A, B = generate_seq_pair(120, 0.1, 0.2, 0.01, 0.01, 7)
    # TODO add stride and phase example, begin and end extensions
    match_moves = [Move(1, 0.0,1,0), Move(3, 0.0,1,0)]
    gap_moves = [Move(3, 1.0,1,0,true), Move(1, 2.0, 1,0,true)]
    alignment_aff = nw_align(A, B, .0, 0.5, match_moves, gap_moves, gap_moves, 0.3)

    @testset "ungapped Alignment Is Equal To Input sequence & Alignments Have Same Length" begin
        ungapped_aligned_A = ungap(alignment[1])
        ungapped_aligned_B = ungap(alignment[2])
        # test that the sequences are unmodified
        @test ungapped_aligned_A == A
        @test ungapped_aligned_B == B
        # additionally assert that alignments have the same length
        @test length(alignment_aff[1]) == length(alignment_aff[2])
    end;

    @testset "alignment Is Reverse Compliment Symmetric" begin

        # TODO gain better understanding why the alignment isn't reversecomplement symmetric

        A_reverse_complement = complement(LongDNA{4}(reverse(string(A))))
        B_reverse_complement = complement(LongDNA{4}(reverse(string(B))))
        alignment_reverse_complement = nw_align(A_reverse_complement, B_reverse_complement, .0, 0.5, match_moves, gap_moves, gap_moves, 0.3,true,true)
        # test if the alignment are the same
        println("A")
        println(alignment_reverse_complement[1])
        println(complement(LongDNA{4}(reverse(string(alignment_aff[1])))))
        @test alignment_reverse_complement[1] == complement(LongDNA{4}(reverse(string(alignment[1])))) skip = true
        println("B")
        println(alignment_reverse_complement[2])
        println(complement(LongDNA{4}(reverse(string(alignment_aff[2])))))
        @test alignment_reverse_complement[2] == complement(LongDNA{4}(reverse(string(alignment[2])))) skip = true
    end;

    @testset "alignment is symmetric in A and B if vertical and horizontal moves are the same" begin
        alignment_symmetric = nw_align(B, A, .0, 0.5, match_moves, gap_moves, gap_moves,0.5,true,true)
        println(alignment_aff[1])
        println(alignment_symmetric[2])
        @test alignment_symmetric[2] == alignment_aff[1]
        println(alignment_aff[2])
        println(alignment_symmetric[1])
        @test alignment_symmetric[1] == alignment_aff[2]
    end;

    @testset "ref2noisy" begin 
        # TODO confirm that this test is working correctly
        refseq, noisySeq = generate_seq_pair(120, 0.02, 0.0, 0.005, 0.005, 10)
        insertions_vertical =  [Move(3, 2.0, 3, 0, 3, 0, true), Move(1, 2.5, 1, 0,false)]
        deletions_horizontal = [Move(3, 2.0, 1, 0, 3, 0, true), Move(1, 2.5, 1, 0,false)]
        # set readingFrame
        readingFrame = 0
        # check that insertions and deletions are correct 
        alignment_ref_2_noisy = nw_align(refseq, noisySeq, .0, 0.5, match_moves, insertions_vertical, deletions_horizontal, 0.3)
        println(join([string(i % 3) for i in 1:length(alignment_ref_2_noisy[1])]))
        println(alignment_ref_2_noisy[1])
        println(alignment_ref_2_noisy[2])
        println("ref")
        println(refseq)
        @test checkRef2NoisyAlignmentRespectsCodonBoundaries(alignment_ref_2_noisy[1],alignment_ref_2_noisy[2],readingFrame) == true
    end;
end;

# TODO check stride and phase are working as expected

# TODO check that extensionable/non-extensionable works also how to score good

# seedChainAlign
@testset "seed_chain_alignment" begin
    A, B = generate_seq_pair(120, 0.1, 0.2, 0.01, 0.01, 7)
    # TODO add stride and phase example, begin and end extensions
    match_moves = [Move(1, 0.0,1,0), Move(3, 0.0,1,0)]
    gap_moves = [Move(3, 1.0,1,0), Move(1, 2.0, 1,0)]
    alignment = seed_chain_align(A, B, .0, 0.5, match_moves, gap_moves, gap_moves, 0.3)
    @testset "ungapped Alignment Is Equal To Input" begin
        ungapped_aligned_A = ungap(alignment[1])
        ungapped_aligned_B = ungap(alignment[2])
        @test ungapped_aligned_A == A
        @test ungapped_aligned_B == B
         # additionally assert that alignments have the same length
        @test length(alignment[1]) == length(alignment[2])
    end;
    @testset "alignment is symmetric in A and B if vertical and horizontal moves are the same" begin
        alignment_symmetric = seed_chain_align(B, A, .0, 0.5, match_moves, gap_moves, gap_moves, 0.3)
        println(alignment[1])
        println(alignment_symmetric[2])
        @test alignment_symmetric[2] == alignment[1]
        println(alignment[2])
        println(alignment_symmetric[1])
        @test alignment_symmetric[1] == alignment[2]
    end;

    @testset "alignment Is Reverse Compliment Symmetric" begin

        # TODO gain better understanding why the alignment isn't reversecomplement symmetric

        A_reverse_complement = complement(LongDNA{4}(reverse(string(A))))
        B_reverse_complement = complement(LongDNA{4}(reverse(string(B))))
        alignment_reverse_complement = seed_chain_align(A_reverse_complement, B_reverse_complement, .0, 0.5, match_moves, gap_moves, gap_moves, 0.3)
        # test if the alignment are the same
        println("A")
        println(alignment_reverse_complement[1])
        println(complement(LongDNA{4}(reverse(string(alignment[1])))))
        @test alignment_reverse_complement[1] == complement(LongDNA{4}(reverse(string(alignment[1])))) skip = true
        println("B")
        println(alignment_reverse_complement[2])
        println(complement(LongDNA{4}(reverse(string(alignment[2])))))
        @test alignment_reverse_complement[2] == complement(LongDNA{4}(reverse(string(alignment[2])))) skip = true
    end;

    @testset "ref2Noisy" begin
        refseq, noisySeq = generate_seq_pair(120, 0.02, 0.0, 0.005, 0.005, 10)
        match_moves = [Move(1, 0.0,1,0), Move(3, 0.0,1,0)]
        insertions_vertical =  [Move(3, 2.0, 3, 0, 3, 0, false), Move(1, 1.5, 1, 0,true)] # TODO set correct stride
        deletions_horizontal = [Move(3, 2.0, 1, 0, 3, 0, false), Move(1, 1.5, 1, 0,true)] # TODO set correct stride
        alignment_ref_2_noisy_seed = seed_chain_align(refseq, noisySeq, .0, 0.5, match_moves, insertions_vertical, deletions_horizontal, 0.5)
        println(join([string(i % 3) for i in 1:length(alignment[1])]))
        println(alignment_ref_2_noisy_seed[1])
        println(alignment_ref_2_noisy_seed[2])
        println("ref")
        println(refseq)
        @test checkRef2NoisyAlignmentRespectsCodonBoundaries(alignment_ref_2_noisy_seed[1],alignment_ref_2_noisy_seed[2],0) == true
    end;
end;

# TODO check kmers