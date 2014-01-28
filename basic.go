package stats

import (
	"math"
)

func FloatSum(x []float64) float64 {
	sum := 0.0
	for _, v := range x {
		sum += v
	}
	return sum
}

func FloatProd(x []float64) float64 {
	prod := 1.0
	for _, v := range x {
		prod *= v
	}
	return prod
}

func FloatMin(x []float64) (float64, int) {
	if len(x) < 1 {
		return math.NaN(), -1
	}
	min := x[0]
	min_idx := 0
	for i, v := range x {
		if v < min {
			min = v
			min_idx = i
		}
	}
	return min, min_idx
}

func FloatMax(x []float64) (float64, int) {
	if len(x) < 1 {
		return math.NaN(), -1
	}
	max := x[0]
	max_idx := 0
	for i, v := range x {
		if v > max {
			max = v
			max_idx = i
		}
	}
	return max, max_idx
}

func FloatMinMax(x []float64) (float64, float64, int, int) {
	if len(x) < 1 {
		return math.NaN(), math.NaN(), -1, -1
	}
	min := x[0]
	max := x[0]
	min_idx := 0
	max_idx := 0
	for i, v := range x {
		if v < min {
			min, min_idx = v, i
		}
		if v > max {
			max, max_idx = v, i
		}
	}
	return min, max, min_idx, max_idx
}

func FloatRange(x []float64) (float64, float64) {
	min, max, _, _ := FloatMinMax(x)
	return min, max
}

func FloatMean(x []float64) float64 {
	m := 0.0
	for i, v := range x {
		m += (v - m) / float64(i+1.0)
	}
	return m
}

func FloatAverage(x []float64, f func(float64) float64) float64 {
	res := 0.0
	for i, v := range x {
		res += (f(v) - res) / float64(i+1.0)
	}
	return res
}

// FloatWeightedMean computes the weighted mean of x with relative
// weights wt. The weights are not assumed to sum to 1;
// their relative values are used.
// x and wt must have the same length;
// elements of wt must all be non-negative;
// these two requirements are not checked.
func FloatWeightedMean(x, wt []float64) float64 {
	return FloatDot(x, wt) / FloatSum(wt)
}

// FloatMSD computes mean squared deviation about a specified value.
func FloatMSD(x []float64, center float64) float64 {
	return FloatAverage(x, func(x float64) float64 {
		return (x - center) * (x - center)
	})
}

func FloatVar(x []float64) (float64, float64) {
	mean, mean_xx := 0.0, 0.0
	for i, v := range x {
		mean += (v - mean) / float64(i+1.0)
		mean_xx += (v*v - mean_xx) / float64(i+1.0)
	}
	n := len(x)
	variance := (mean_xx - mean*mean) * float64(n) / float64(n-1)
	return variance, mean
}

func FloatSd(x []float64) (float64, float64) {
	v, m := FloatVar(x)
	return math.Sqrt(v), m
}

// LogSumExp computes the log of the sum of the exponentials of the
// values in x, in a numerically stable way.
func LogSumExp(x []float64) float64 {
	maxval, _ := FloatMax(x)
	if math.IsInf(maxval, 0) {
		return maxval
	}
	var lse float64
	for _, val := range x {
		lse += math.Exp(val - maxval)
	}
	return math.Log(lse) + maxval
}
