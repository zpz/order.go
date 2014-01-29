package stats

import (
	"github.com/zpz/matrix.go/dense"
	"math"
)

// FloatCov computes the covariance matrix between member slices of data.
// All member slices of data must have the same length;
// this is not checked.
func FloatCov(data [][]float64, out *dense.Dense) *dense.Dense {
	p := len(data)
	if p < 1 {
		return nil
	}

	n := len(data[0])

	out = use_dense(out, p, p)

	means := make([]float64, p)
	for i := 0; i < p; i++ {
		v, m := FloatVar(data[i])
		means[i] = m
		out.Set(i, i, v)
	}

	factor := float64(n) / float64(n-1.0)

	for i := 0; i < p; i++ {
		for j := i + 1; j < p; j++ {
			v := 0.0
			for k := 0; k < n; k++ {
				vv := (data[i][k] - means[i]) * (data[j][k] - means[j])
				v += (vv - v) / float64(k+1.0)
			}
			v *= factor
			out.Set(i, j, v)
			out.Set(j, i, v)
		}
	}

	return out
}

// x is the output of FloatCov.
// out can be the same as x, amounting to in-place conversion.
func Cov2Cor(x, out *dense.Dense) *dense.Dense {
	p := x.Rows()
	out = use_dense(out, p, p)

	sd := x.GetDiag(nil)
	FloatTransform(sd, math.Sqrt, sd)

	out.FillDiag(1.0)
	for i := 0; i < p; i++ {
		for j := i + 1; j < p; j++ {
			v := x.Get(i, j) / sd[i] / sd[j]
			out.Set(i, j, v)
			out.Set(j, i, v)
		}
	}
	return out
}

func FloatCor(data [][]float64, out *dense.Dense) *dense.Dense {
	out = FloatCov(data, out)
	Cov2Cor(out, out)
	return out
}

// FloatWeightedCov computes covariance matrix for weighted data.
func FloatWeightedCov(
	data [][]float64,
	wt []float64,
	out *dense.Dense) *dense.Dense {

	p := len(data)
	if p < 1 {
		return nil
	}
	n := len(data[0])

	// Make a copy of the weights and normalize.
	w := FloatScale(wt, 1.0/FloatSum(wt), nil)

	out = use_dense(out, p, p)

	means := make([]float64, p)
	for i := 0; i < p; i++ {
		means[i] = FloatDot(data[i], w)
	}

	factor := 1.0 / (1.0 - FloatDot(w, w))

	for i := 0; i < p; i++ {
		for j := i; j < p; j++ {
			// sum(w_i * (x_i - m_x) * (y_i - m_y))
			// = sum(w_i * x_i * y_i) - sum(w_i * x_i * m_y)
			//   - sum(w_i * m_x * y_i) + sum(w_i * m_x * m_y)
			// = sum(w_i * x_i * y_i) - m_x * m_y
			v := 0.0
			for k := 0; k < n; k++ {
				v += w[k] * data[i][k] * data[j][k]
			}
			v *= factor
			out.Set(i, j, v)
			if i != j {
				out.Set(j, i, v)
			}
		}
	}

	return out
}
