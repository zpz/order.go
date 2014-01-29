package stats

import (
	"math"
	"sort"
)

// FloatMedian returns the median of data.
func FloatMedian(data []float64) float64 {
	n := len(data)
	if n < 1 {
		return math.NaN()
	}
	if n == 1 {
		return data[0]
	}
	if n == 2 {
		return (data[0] + data[1]) * 0.5
	}

	return float_median(float_clone(data))
}

// float_median returns the median of data.
// Note that elements in data will be re-ordered in this operation.
func float_median(data []float64) float64 {
	n := len(data)
	mid := n / 2
	if mid*2 == n {
		PartialSort(sort.Float64Slice(data), mid)
		PartialSort(sort.Float64Slice(data[:mid]), mid-1)
		return (data[mid] + data[mid-1]) * 0.5
	}
	PartialSort(sort.Float64Slice(data), mid)
	return data[mid]
}

// IntMedian returns the median of data.
func IntMedian(data []int) float64 {
	n := len(data)
	if n < 1 {
		return math.NaN()
	}
	if n == 1 {
		return float64(data[0])
	}
	if n == 2 {
		return float64(data[0]+data[1]) * 0.5
	}

	return int_median(int_clone(data))
}

// int_median returns the median of data.
// Note that elements in data will be re-ordered in this operation.
func int_median(data []int) float64 {
	n := len(data)
	mid := n / 2
	if mid*2 == n {
		PartialSort(sort.IntSlice(data), mid)
		PartialSort(sort.IntSlice(data[:mid]), mid-1)
		return float64(data[mid]+data[mid-1]) * 0.5
	}
	PartialSort(sort.IntSlice(data), mid)
	return float64(data[mid])
}
