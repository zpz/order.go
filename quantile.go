package stats

import (
	"math"
	"sort"
)

// Median returns the median of data.
// Note that elements in data will be re-ordered in this operation.
func Median(data []float64) float64 {
	n := len(data)
	if n < 1 {
		return math.NaN()
	}
	if n == 1 {
		return data[0]
	}

	mid := n / 2

	if mid*2 == n {
		PartialSort(sort.Float64Slice(data), mid)
		PartialSort(sort.Float64Slice(data[:mid]), mid-1)
		return (data[mid] + data[mid-1]) * 0.5
	}

	PartialSort(sort.Float64Slice(data), mid)
	return data[mid]
}
