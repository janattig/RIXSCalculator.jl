#########################
# PERMUTATION FUNCTIONS #
#########################

# permutations of a list and sign associated to sorting this list
function permutation_sign(numbers :: Vector{Int64})
    return permutation_sign!(deepcopy(numbers))
end

function permutation_sign!(numbers :: Vector{Int64})
    # define a sign
    s = 1
    # get N = length of list
    N = length(numbers)
    # permutate until everything is in ascending order
    p = true
    @inbounds while p==true
        p = false
        for e in 1:N-1
            if numbers[e] > numbers[e+1]
                numbers[e], numbers[e+1] = numbers[e+1], numbers[e]
                p = true
                s = -s
            elseif numbers[e] == numbers[e+1]
                s = 0
                return s
            end
        end
    end
    # return the sign
    return s
end
function permutated(numbers :: Vector{Int64})
    # get the permutation sign and new permutation
    numbers_new = deepcopy(numbers)
    s = permutation_sign!(numbers_new)
    # return sign and new list
    return (s, numbers_new)
end

# n element permutations of 1:x
function getPositivePermutations(n :: Integer, x :: Integer) :: Vector{Vector{Int64}}
    # return combinations
    return collect(combinations(1:x, n))
end
