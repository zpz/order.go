package stats

func IntSeq(from, to int, out []int) []int {
	n := to - from
	if out == nil {
		out = make([]int, n)
	} else {
		if len(out) != n {
			panic("wrong length for output")
		}
	}

	for i := 0; i < n; i++ {
		out[i] = from
		from++
	}
	return out
}
