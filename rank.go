package stats

func rank(data OrderInterface, r []int, stable bool) []int {
	n := data.Len()

	idx := make([]int, n)
	if stable {
		StableOrder_(data, idx)
	} else {
		Order_(data, idx)
	}

	r = use_int_slice(r, n)

	order_to_rank(idx, r)

	return r
}

func order_to_rank(in, out []int) {
	for i := range in {
		out[in[i]] = i
	}
}

// Rank_ returns the sample ranks of the elements in data.
// For example, the ranks of {3, 5, 2, 6} would be {1, 2, 0, 3}.
// If there are identical elements in data, their ranks may not follow
// the order of their appearances in data.
// Note that the elements of data are re-ordered in this function.
// If r is nil, a new slice will be allocated and used;
// otherwise r must be of the correct length.
func Rank_(data OrderInterface, r []int) []int {
	return rank(data, r, false)
}

// StableRank_ is similar to Rank_, except that when there are identical
// elements in data, their orders are preserved.
// For example, the ranks of {3, 5, 2, 2, 6} would be {2, 3, 0, 1, 4}.
func StableRank_(data OrderInterface, r []int) []int {
	return rank(data, r, true)
}

// Rank is similar to Rank_, except that it is specialized for float64
// input data.
func Rank(x []float64, r []int) []int {
	return Rank_(make_float64_index_slice(x), r)
}

// StableRank is similar to StableRank_, except that it is specialized
// for float64 input data.
func StableRank(x []float64, r []int) []int {
	return StableRank_(make_float64_index_slice(x), r)
}

// IntRank is similar to Rank, except that it is for an int slice.
func IntRank(x []int, r []int) []int {
	return Rank_(make_int_index_slice(x), r)
}

// IntStableRank is similar to StableRank, except that it is for an int slice.
func IntStableRank(x []int, r []int) []int {
	return StableRank_(make_int_index_slice(x), r)
}

// StringRank is similar to Rank, except that it is for a string slice.
func StringRank(x []string, r []int) []int {
	return Rank_(make_string_index_slice(x), r)
}

// StringStableRank is similar to StableRank, except that it is for a string slice.
func StringStableRank(x []string, r []int) []int {
	return StableRank_(make_string_index_slice(x), r)
}
