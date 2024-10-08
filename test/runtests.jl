using Test
using Random
using BioSequences

# files to test
include("../src/needleman_wunsch.jl")
include("../src/seed_chain_align.jl")
#help functions
include("../test_case_creation_HIV/sequence_generator.jl")

# NW non-affine unittest
@testset "needleman-wunsch non-affine" begin
        @testset "test for general alignment output" begin
            for _ in 1:1
                A, B = generate_seq_pair(120, 0.1, 0.2, 0.01, 0.01, 10)
                # default moveset (note symmetric in horizontal and vertical)
                match_moves = [Move(1, 0.0,1,0), Move(3, 0.0,1,0)]
                gap_moves =   [Move(1, 2.2,1,0), Move(3, 2.0,1,0)]
                alignment = nw_align(A, B, .0, 0.5, match_moves, gap_moves, gap_moves,true)
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
                    alignment_reverse_complement = nw_align(A_reverse_complement, B_reverse_complement, .0, 0.5, match_moves, gap_moves, gap_moves,true)

                    # test if the alignment are the same
                    println("reverse complement")
                    println(alignment_reverse_complement[1])
                    println(alignment_reverse_complement[2])
                    println("normal")
                    println(complement(LongDNA{4}(reverse(string(alignment[1])))))
                    println(complement(LongDNA{4}(reverse(string(alignment[2])))))
                    println("test scores match")
                    println("reversecomplement: ", alignment_reverse_complement[3])
                    println("normal: ", alignment[3])
                    @test alignment_reverse_complement[3] == alignment[3]
                    @test alignment_reverse_complement[1] == complement(LongDNA{4}(reverse(string(alignment[1])))) skip=true    
                    @test alignment_reverse_complement[2] == complement(LongDNA{4}(reverse(string(alignment[2])))) skip=true    
                end;
            
            end
        end;
end;


# NW-affine unittest
@testset "needleman-wunsch affine" begin
    A, B = generate_seq_pair(120, 0.01, 0.02, 0.01, 0.01, 7)
    # TODO add stride and phase example, begin and end extensions
    match_moves = [Move(1, 0.0,1,0), Move(3, 0.0,1,0)]
    gap_moves = [Move(3, 1.0,1,0,true), Move(1, 0.9, 1,0,true)]
    alignment_aff = nw_align(A, B, .0, 0.5, match_moves, gap_moves, gap_moves, 0.3,true,true,true)

    @testset "ungapped Alignment Is Equal To Input sequence & Alignments Have Same Length" begin
        ungapped_aligned_A = ungap(alignment_aff[1])
        ungapped_aligned_B = ungap(alignment_aff[2])
        # test that the sequences are unmodified
        @test ungapped_aligned_A == A
        @test ungapped_aligned_B == B
        # additionally assert that alignments have the same length
        @test length(alignment_aff[1]) == length(alignment_aff[2])
    end;
    @testset "alignment is symmetric in A and B if vertical and horizontal moves are the same" begin
        alignment_symmetric = nw_align(B, A, .0, 0.5, match_moves, gap_moves, gap_moves,0.3,true,true,true)
        println(alignment_aff[1])
        println(alignment_symmetric[2])
        @test alignment_symmetric[2] == alignment_aff[1]
        println(alignment_aff[2])
        println(alignment_symmetric[1])
        @test alignment_symmetric[1] == alignment_aff[2]
    end;
    @testset "alignment Is Reverse Compliment Symmetric" begin
        # TODO gain better understanding why the alignment isn't reversecomplement symmetric
        # I think this might work better on real data. This example might be too synthetic
        A_reverse_complement = complement(LongDNA{4}(reverse(string(A))))
        B_reverse_complement = complement(LongDNA{4}(reverse(string(B))))
        alignment_reverse_complement = nw_align(A_reverse_complement, B_reverse_complement, .0, 0.5, match_moves, gap_moves, gap_moves, 0.3,true,true,true)
        # test if the alignment are the same
        println("reverse complimented")
        println(alignment_reverse_complement[1])
        println(alignment_reverse_complement[2])
        println("normal")
        println(complement(LongDNA{4}(reverse(string(alignment_aff[1])))))
        println(complement(LongDNA{4}(reverse(string(alignment_aff[2])))))
        println("test scores match")
        println("reversecomplement: ", alignment_reverse_complement[3])
        println("normal: ", alignment_aff[3])
        @test alignment_reverse_complement[3] == alignment_aff[3]
        @test alignment_reverse_complement[1] == complement(LongDNA{4}(reverse(string(alignment_aff[1])))) skip=true
        @test alignment_reverse_complement[2] == complement(LongDNA{4}(reverse(string(alignment_aff[2])))) skip=true
    end;
end;

# TODO check stride and phase are working as expected

# TODO check that extensionable/non-extensionable works and is the same as non-affine

# seedChainAlign
@testset "seed_chain_alignment" begin
    A, B = generate_seq_pair(120, 0.05, 0.1, 0.01, 0.01, 7)
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
        alignment_symmetric = seed_chain_align(B, A, .0, 0.5, match_moves, gap_moves, gap_moves, 0.3, 18)
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
end;
