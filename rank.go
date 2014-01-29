package stats

func rank(data OrderInterface, r []int, stable bool) []int {
	n := data.Len()

	idx := make([]int, n)
	if stable {
		StableOrder(data, idx)
	} else {
		Order(data, idx)
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

// Rank returns the sample ranks of the elements in data.
// For example, the ranks of {3, 5, 2, 6} would be {1, 2, 0, 3}.
// If there are identical elements in data, their ranks may not follow
// the order of their appearances in data.
// Note that the elements of data are re-ordered in this function.
// If r is nil, a new slice will be allocated and used;
// otherwise r must be of the correct length.
func Rank(data OrderInterface, r []int) []int {
	return rank(data, r, false)
}

// StableRank is similar to Rank_, except that when there are identical
// elements in data, their orders are preserved.
// For example, the ranks of {3, 5, 2, 2, 6} would be {2, 3, 0, 1, 4}.
func StableRank(data OrderInterface, r []int) []int {
	return rank(data, r, true)
}

// FloatRank is similar to Rank, except that it is specialized for float64
// input data.
func FloatRank(x []float64, r []int) []int {
	return Rank(make_float64_index_slice(x), r)
}

// FloatStableRank is similar to StableRank, except that it is specialized
// for float64 input data.
func FloatStableRank(x []float64, r []int) []int {
	return StableRank(make_float64_index_slice(x), r)
}

// IntRank is similar to Rank, except that it is for an int slice.
func IntRank(x []int, r []int) []int {
	return Rank(make_int_index_slice(x), r)
}

// IntStableRank is similar to StableRank, except that it is for an int slice.
func IntStableRank(x []int, r []int) []int {
	return StableRank(make_int_index_slice(x), r)
}

// StringRank is similar to Rank, except that it is for a string slice.
func StringRank(x []string, r []int) []int {
	return Rank(make_string_index_slice(x), r)
}

// StringStableRank is similar to StableRank, except that it is for a string slice.
func StringStableRank(x []string, r []int) []int {
	return StableRank(make_string_index_slice(x), r)
}
