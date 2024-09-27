using BioSequences

mutable struct Move
    step::Int64
    score::Float64
    # refers to readingFrame of top sequence
    horizontal_stride::Int64
    horizontal_phase::Int64
    # refers to readingFrame of bottom sequence
    vertical_stride::Int64
    vertical_phase::Int64
    # allows moves to be extended from 
    extensionAble::Bool

    # moves not considering stride and phase
    Move(step::Int64,score::Float64) = new(step,score,1,0,1,0,false)
    Move(step::Int64,score::Float64,extensionAble::Bool) = new(step,score,1,0,1,0,extensionAble)
    # moves that assume both sequences have the same reading frame # TODO check this thinking is correct
    Move(step::Int64,score::Float64,stride::Int64,phase::Int64) = new(step,score,stride,phase,stride,phase,false)
    Move(step::Int64,score::Float64,stride::Int64,phase::Int64,extensionAble::Bool) = new(step,score,stride,phase,stride,phase,false)
    # normal definition 
    Move(step::Int64,score::Float64,horizontal_stride::Int64,horizontal_phase::Int64,vertical_stride::Int64,vertical_phase::Int64,
        extensionAble::Bool) = new(step,score,horizontal_stride,horizontal_phase,vertical_stride,vertical_phase,extensionAble)
end


# Convert NucleicAcid to integer A -> 1, C -> 2, G -> 3, T -> 4
function toInt(x::NucleicAcid)
    trailing_zeros(reinterpret(UInt8,x)) + 1
end

function simple_match_penalty_matrix(match_score, mismatch_score)
    m = zeros(4, 4)
    for i in 1:4, j in 1:4
        if i == j
            m[i, j] = match_score
        else
            m[i, j] = mismatch_score
        end
    end
    return m
end

# match and mismatch matrix
function nw_align(A::LongDNA{4}, B::LongDNA{4}, match_score::Float64, mismatch_score::Float64, 
                match_moves::Vector{Move}, vgap_moves::Vector{Move}, hgap_moves::Vector{Move}) 

    nw_align(A, B, simple_match_penalty_matrix(match_score, mismatch_score), match_moves, vgap_moves, hgap_moves) 
end

#Needleman Wunsch alignment without affine scoring
function nw_align(A::LongDNA{4}, B::LongDNA{4}, match_score_matrix::Array{Float64, 2},
                match_moves::Vector{Move}, vgap_moves::Vector{Move}, hgap_moves::Vector{Move})

    n, m = length(A), length(B)

    #Margins in the dp_matrix streamlines the code by avoiding boundschecking
    column_offset = maximum(k -> k.step, vcat(match_moves, hgap_moves)) + 1
    row_offset = maximum(k -> k.step, vcat(match_moves, vgap_moves)) + 1

    #The matrix and the sequences are expanded according to the offsets
    column_boundary = n + column_offset
    row_boundary = m + row_offset
    
    #extend sequences to fix offset
    A2 = LongDNA{4}("A")^(column_offset - 1) * A
    B2 = LongDNA{4}("A")^(row_offset - 1) * B

    #Initialize DP matrix
    #The cell at [x + column_offset, y + column_offset] is the score of the best alignment of A[1 : x] with B[1 : y]
    dp_matrix = fill(Inf64, row_boundary,column_boundary)
    # Assign score 0 to starting position
    dp_matrix[row_offset, column_offset] = 0.0
    # itterate through dp_matrix
    for row_index ∈ 1 + row_offset:row_boundary
        for column_index ∈ 1 + column_offset:column_boundary
            # calculate position we are moving to (1 indexed)
            # TODO understand why this works
            top_sequence_pos = column_index-column_offset
            left_sequence_pos = row_index-row_offset

            # find the best diagonal move
            for k ∈ match_moves
                # TODO remove extension of A2 and B2 and use top_seq and left_seq coords
                mismatch_sum = sum(t -> match_score_matrix[toInt(A2[column_index - t]), toInt(B2[row_index - t])], 1 : k.step)
                dp_matrix[row_index, column_index] = min(
                    dp_matrix[row_index, column_index],
                    dp_matrix[row_index-k.step,column_index-k.step]+k.score+mismatch_sum
                )
            end

            for k ∈ vgap_moves
                # the lowest score of vertical gaps is assigned to the current position in the matrix
                if (top_sequence_pos) % k.vertical_stride == k.vertical_phase && #TODO rename
                        (left_sequence_pos - k.step) % k.horizontal_stride == k.horizontal_phase && #TODO rename
                        dp_matrix[row_index-k.step, column_index] + k.score < dp_matrix[row_index, column_index]
                    # update current matrix position with lower score
                    dp_matrix[row_index, column_index] = dp_matrix[row_index - k.step, column_index] + k.score
                end
            end
            
            for k ∈ hgap_moves
                # the lowest score of horizontal gaps is asigned to the current position in the matrix if the score is lower than the current one
                if (top_sequence_pos - k.step) % k.vertical_stride == k.vertical_phase && #TODO rename
                        (left_sequence_pos) % k.horizontal_stride == k.horizontal_phase && #TODO rename
                        dp_matrix[row_index, column_index - k.step] + k.score < dp_matrix[row_index, column_index]
                    # update current matrix position with lower score
                    dp_matrix[row_index, column_index] = dp_matrix[row_index, column_index - k.step] + k.score
                end
            end
        end
    end
    display(dp_matrix[1:end,1:end])
    # Backtracking
    x = column_boundary
    y = row_boundary
    res_A = LongDNA{4}("")
    res_B = LongDNA{4}("")
    while x > column_offset || y > row_offset
        #println("col ",x)
        #println("row ",y)
        top_sequence_pos = x-column_offset
        left_sequence_pos = y-row_offset
        if x == column_offset # first row
            push!(res_A, DNA_Gap)
            push!(res_B, B2[y - 1])
            y -= 1
        elseif y == row_offset # first column
            push!(res_A, A2[x - 1])
            push!(res_B, DNA_Gap)
            x -= 1
        else
            # iterate through diagonal (match) moves
            for k ∈ match_moves
                # calculate total (mis-)match score from move
                s = sum(t -> match_score_matrix[toInt(A2[x-t]), toInt(B2[y-t])], 1 : k.step)
                # check if the move leads to the current cell
                if dp_matrix[y,x] == dp_matrix[y - k.step,x - k.step] + k.score + s
                    
                    # write the resulting sequences
                    for i ∈ 1:k.step
                        push!(res_A, A2[x - i])
                        push!(res_B, B2[y - i])
                    end
                    x -= k.step
                    y -= k.step
                    break
                end
            end

            for k ∈ vgap_moves
                # iterate through possible vertical moves and check if the move leads to the current cell
                if (top_sequence_pos) % k.vertical_stride == k.vertical_phase && 
                        (left_sequence_pos - k.step) % k.horizontal_stride == k.horizontal_phase && 
                        dp_matrix[y, x] == dp_matrix[y- k.step,x] + k.score
                    # perform move
                    for i ∈ 1:k.step
                        push!(res_A, DNA_Gap)
                        push!(res_B, B2[y - i])
                    end
                    # update
                    y -= k.step
                    # perform move
                    break
                end
            end

            for k ∈ hgap_moves
                # iterate through possible horizontal moves and check if move leads to the current cell
                if (top_sequence_pos - k.step) % k.vertical_stride == k.vertical_phase && 
                    (left_sequence_pos) % k.horizontal_stride == k.horizontal_phase && 
                    dp_matrix[y,x] == dp_matrix[y,x-k.step] + k.score
                    #perform move
                    for i ∈ 1:k.step
                        push!(res_A, A2[x - i])
                        push!(res_B, DNA_Gap)
                    end
                    # update
                    x -= k.step
                    break
                end
            end

        end
    end
    println(dp_matrix[row_boundary,column_boundary])
    return reverse(res_A), reverse(res_B)
end

# Needleman Wunsch alignment with affine scoring
function nw_align(A::LongDNA{4}, B::LongDNA{4}, match_score::Float64, mismatch_score::Float64, match_moves::Vector{Move}, 
        vgap_moves::Vector{Move}, hgap_moves::Vector{Move}, extension_score::Float64, 
        edge_extension_begin=false::Bool, edge_extension_end=false::Bool)
 
    nw_align(A, B, simple_match_penalty_matrix(match_score, mismatch_score), match_moves, vgap_moves, hgap_moves, extension_score, edge_extension_begin, edge_extension_end) 
end

function nw_align(A::LongDNA{4}, B::LongDNA{4}, match_score_matrix::Array{Float64, 2}, match_moves::Vector{Move}, vgap_moves::Vector{Move}, hgap_moves::Vector{Move}, extension_score::Float64, edge_extension_begin=false::Bool, edge_extension_end=false::Bool)
    
    n, m = length(A), length(B)

    if (extension_score < 0)
        # Do non-affine NW
        return nw_align(A, B, match_score_matrix, match_moves, vgap_moves, hgap_moves)
    end

    # Offset indeces to avoid bounds-checking
    column_offset = maximum(k -> k.step, vcat(match_moves, hgap_moves)) + 1
    row_offset = maximum(k -> k.step, vcat(match_moves, vgap_moves)) + 1
    column_boundary = n + column_offset
    row_boundary = m + row_offset

    # Length of sequences and matrices are increased according to offset
    A2 = LongDNA{4}("A")^(column_offset - 1) * A
    B2 = LongDNA{4}("A")^(row_offset - 1) * B

    # Initialize DP matrix
    # The cell at [x + column_offset, y + column_offset] is the score of the best alignment of A[1 : x] with B[1 : y]
    dp_matrix = fill(Inf64, row_boundary, column_boundary)

    # Assign score 0 to the empty alignment
    dp_matrix[row_offset, column_offset] = 0.0

    # Affine moves requires two extra DP matrices
    vaffine_matrix = fill(Inf64, row_boundary, column_boundary)
    haffine_matrix = fill(Inf64, row_boundary, column_boundary)

    # end/start extension biase (optional)
    start_extension_biase = 0.0 + extension_score
    end_extension_biase = 0.0 + extension_score
    
    # allow starting in extending
    if edge_extension_begin
        vaffine_matrix[row_offset,column_offset] = start_extension_biase
        haffine_matrix[row_offset,column_offset] = start_extension_biase
    end

    #Main DP -step
    for row_index ∈ 1 + row_offset : row_boundary
        for column_index ∈ 1 + column_offset : column_boundary
            top_sequence_pos = column_index-column_offset
            left_sequence_pos = row_index-row_offset

            # find the best diagonal move
            for k ∈ match_moves
                mismatch_sum = sum(t -> match_score_matrix[toInt(A2[column_index - t]), toInt(B2[row_index - t])], 1 : k.step)
                dp_matrix[row_index, column_index] = min(
                    dp_matrix[row_index, column_index],dp_matrix[row_index-k.step,column_index-k.step]+k.score+mismatch_sum)
            end

            # finds the best vertical move
            for k ∈ vgap_moves
                if (top_sequence_pos) % k.vertical_stride == k.vertical_phase &&
                    (left_sequence_pos-k.step) % k.horizontal_stride == k.horizontal_phase
                    if k.extensionAble # TODO rename
                        vaffine_matrix[row_index, column_index] = min(
                            vaffine_matrix[row_index, column_index],
                            vaffine_matrix[row_index - k.step, column_index] + extension_score * k.step,
                            dp_matrix[row_index - k.step, column_index] + k.score
                        )
                    else
                        dp_matrix[row_index,column_index] = min(
                            dp_matrix[row_index,column_index],
                            dp_matrix[row_index - k.step, column_index] + k.score
                        )
                    end
                end
            end

            # finds the best horizontal move
            for k ∈ hgap_moves
                if (top_sequence_pos-k.step) % k.vertical_stride == k.vertical_phase && 
                    (left_sequence_pos) % k.horizontal_stride == k.horizontal_phase
                    if k.extensionAble # TODO rename
                        haffine_matrix[row_index, column_index] = min(
                            haffine_matrix[row_index, column_index],
                            haffine_matrix[row_index, column_index - k.step] + extension_score * k.step,
                            dp_matrix[row_index, column_index - k.step] + k.score
                        )
                    else
                        dp_matrix[row_index,column_index] = min(
                            dp_matrix[row_index,column_index],
                            dp_matrix[row_index, column_index - k.step] + k.score
                        )
                    end
                end
            end

            # find overall best move
            dp_matrix[row_index, column_index] = min(
                dp_matrix[row_index, column_index], 
                haffine_matrix[row_index, column_index],
                vaffine_matrix[row_index, column_index]
            )
            # allows us to end in extension state (add constant to expression to adjust sensitivity of extension)
            # TODO make sure this is done in Backtracking
            if column_index == column_boundary && edge_extension_end
                vaffine_matrix[row_index,column_index] = dp_matrix[row_index,column_index] + end_extension_biase
            elseif row_index == row_boundary && edge_extension_end
                haffine_matrix[row_index,column_index] = dp_matrix[row_index,column_index] + end_extension_biase
            end
        end
    end

    #display(dp_matrix[1:end,1:end])
    # Backtracking
    res_A = LongDNA{4}("")
    res_B = LongDNA{4}("")
    # Start at the final cell
    x = column_boundary
    y = row_boundary
    # Flags for affine moves
    must_move_ver = false
    must_move_hor = false
    while x > column_offset || y > row_offset
        top_sequence_pos = x-column_offset
        left_sequence_pos = y-row_offset
        #println("col ", x)
        #println("row ", y)
        if x == column_offset # first row
            push!(res_A, DNA_Gap)
            push!(res_B, B2[y - 1])
            y -= 1
        elseif y == row_offset # first col
            push!(res_A, A2[x - 1])
            push!(res_B, DNA_Gap)
            x -= 1
        else
            # record previous position
            px = x
            py = y

            if !must_move_hor
                
                for k ∈ vgap_moves
                    # TODO check where infinity and such
                    !( ((top_sequence_pos) % k.vertical_stride == k.vertical_phase) &&
                       ((left_sequence_pos-k.step) % k.horizontal_stride == k.horizontal_phase) ) ? continue : 
                    # check if the move leads to the current cell
                    if k.extensionAble
                        current_score = must_move_ver ? vaffine_matrix[y, x] : dp_matrix[y, x]
                        can_move_affine = (current_score == vaffine_matrix[y-k.step, x] + extension_score * k.step)
                        can_move_regular = (current_score == dp_matrix[y-k.step, x] + k.score)
                        #gap extension
                    else
                        current_score = dp_matrix[y,x]
                        can_move_affine = (false)
                        can_move_regular = (current_score == dp_matrix[y-k.step,x] + k.score)
                    end
                    
                    if can_move_affine || can_move_regular
                        for i ∈ 1 : k.step
                            push!(res_A, DNA_Gap)
                            push!(res_B, B2[y - i])
                        end
                        y -= k.step

                        # constrain next move
                        if !(y == row_boundary && edge_extension_end)
                            must_move_ver = !can_move_regular
                        end
                        break
                    end
                end
            end

            if !must_move_ver

                # iterate through horizontal Move moves
                for k ∈ hgap_moves
                    !((top_sequence_pos-k.step) % k.vertical_stride == k.vertical_phase && 
                      (left_sequence_pos) % k.horizontal_stride == k.horizontal_phase) ? continue :
                    # check if the move leads to the current cell
                    if k.extensionAble
                        current_score = must_move_hor ? haffine_matrix[y, x] : dp_matrix[y, x]
                        can_move_affine = (current_score == haffine_matrix[y, x-k.step] + extension_score * k.step)
                        can_move_regular = (current_score == dp_matrix[y,x-k.step] + k.score)
                    else
                        current_score = dp_matrix[y, x]
                        can_move_affine = (false)
                        can_move_regular = (current_score == dp_matrix[y,x-k.step] + k.score)
                    end

                    if can_move_affine || can_move_regular
                        
                        for i ∈ 1:k.step # TODO why is it matching B2 if move is horizontal
                            push!(res_A, A2[x - i])
                            push!(res_B, DNA_Gap)
                        end
                        x -= k.step

                        # constrain next move
                        if !(x == column_boundary && edge_extension_end)
                            must_move_hor = !can_move_regular
                        end
                        break
                    end
                end
            end

            if !must_move_hor && !must_move_ver

                # iterate through digonal match moves
                for k ∈ match_moves
                    #calculate total (mis-)match score 
                    s = sum(t -> match_score_matrix[toInt(A2[x-t]), toInt(B2[y-t])], 1 : k.step)
                    # check if the move leads to the current cell
                    if dp_matrix[y, x] == dp_matrix[y - k.step,x - k.step] + k.score + s
                        # record the path
                        #println("match length ", k.step)
                        for i ∈ 1:k.step
                            push!(res_A, A2[x - i])
                            push!(res_B, B2[y - i])
                        end
                        x -= k.step
                        y -= k.step
                        break
                    end
                end
            end

            # if no move was found
            if px == x && py == y
                error("Backtracking failed")
            end
        end
    end
    return reverse(res_A), reverse(res_B)
end

#A = LongDNA{2}("CCGACCCCGATTCCCGTTA")
#B = LongDNA{2}("GACCTTCACCGTTA")
#                 TCACCG---GTCTA---GCCCCCTACAAAAGGCGACATCTGCCCTGGCCG---ATTGGCTACCAGACGA
#A = LongDNA{4}("ACTGCT")
#                 TCACCGCTGATCGACTGGCACCCTACAAAA---GACATCTGCCCTGGCCGCTGGTTGGCTACCAGACGA
#B = LongDNA{4}("GCCCTAC")
#match_moves = [Move(1, 0.0,1,0), Move(3, 0.0,1,0)]
#gap_moves = [Move(3, 1.0,1,0), Move(1, 2.0,1,0)]
#alignment = nw_align(A, B, .0, 0.5, match_moves, gap_moves, gap_moves, 0.5)
#println(alignment[1])
#println(alignment[2])