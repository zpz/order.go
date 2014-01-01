package stats




func rank(data OrderInterface, r []int, stable bool) []int {
    n := data.Len()

    idx := make([]int, n)
    if stable {
        StableOrder_(data, idx)
    } else {
        Order_(data, idx)
    }

    if r == nil {
        r = make([]int, n)
    }
    order_to_rank(idx, r)

    return r[0:n]
}




func order_to_rank(in, out []int) {
    for i := range in {
        out[in[i]] = i
    }
}


// Rank_ returns the sample ranks of the elements in the collection data.
// For example, conceptually, Rank_([]int{3, 5, 2, 6}) returns []int{1, 2, 0, 3}.
// Note that the elements of data are re-ordered in this function.
func Rank_(data OrderInterface, r []int) []int {
    return rank(data, r, false)
}




func StableRank_(data OrderInterface, r []int) []int {
    return rank(data, r, true)
}





func Rank(x []float64, r []int) []int {
    return Rank_(make_float64_index_slice(x), r)
}




func StableRank(x []float64, r []int) []int {
    return StableRank_(make_float64_index_slice(x), r)
}




func IntRank(x []int, r []int) []int {
    return Rank_(make_int_index_slice(x), r)
}




func IntStableRank(x []int, r []int) []int {
    return StableRank_(make_int_index_slice(x), r)
}




func StringRank(x []string, r []int) []int {
    return Rank_(make_string_index_slice(x), r)
}




func StringStableRank(x []string, r []int) []int {
    return StableRank_(make_string_index_slice(x), r)
}


