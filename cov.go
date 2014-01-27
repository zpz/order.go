package stats

import (
	"github.com/zpz/matrix.go/dense"
	"math"
)

// FloatCov computes the covariance between x and y.
func FloatCov(x, y []float64) float64 {
	n := len(x)
	assert(len(y) == n, "lengths of x and y differ")

	m_x, m_y, m_xy := 0.0, 0.0, 0.0
	for i := 0; i < n; i++ {
		r := 1.0 / float64(i+1.0)
		vx := x[i]
		vy := y[i]
		m_x += (vx - m_x) * r
		m_y += (vy - m_y) * r
		m_xy += (vx*vy - m_xy) * r
	}

	cov := m_xy - m_x*m_y
	cov *= float64(n) / float64(n-1.0)
	return cov
}

// FloatCor computes the correlation coefficient between
// x and y, returning the correlation coef and covariance between x and
// y, as well as the standard deviations of x and y.
func FloatCor(x, y []float64) (cor, cov, sd_x, sd_y float64) {
	n := len(x)
	assert(len(y) == n, "lengths of x and y differ")

	m_x, m_y, m_xx, m_yy, m_xy := 0.0, 0.0, 0.0, 0.0, 0.0
	for i := 0; i < n; i++ {
		r := 1.0 / float64(i+1.0)
		vx := x[i]
		vy := y[i]
		m_x += (vx - m_x) * r
		m_y += (vy - m_y) * r
		m_xx += (vx*vx - m_xx) * r
		m_yy += (vy*vy - m_yy) * r
		m_xy += (vx*vy - m_xy) * r
	}

	r := float64(n) / float64(n-1)
	cov = (m_xy - m_x*m_y) * r
	sd_x = math.Sqrt((m_xx - m_x*m_x) * r)
	sd_y = math.Sqrt((m_yy - m_y*m_y) * r)
	cor = cov / (sd_x * sd_y)

	return cor, cov, sd_x, sd_y
}

func FloatCovMatrix(
	// All member slices of data must have the same length;
	// this is not checked.
	data [][]float64,
	out *dense.Dense) *dense.Dense {

	p := len(data)
	if p < 1 {
		return nil
	}

	n := len(data[0])

	out = use_matrix(out, p, p)

	means := make([]float64, p)
	for i := 0; i < p; i++ {
		means[i] = FloatMean(data[i])
	}

	factor := float64(n) / float64(n-1.0)

	for i := 0; i < p; i++ {
		for j := i; j < p; j++ {
			v := FloatDot(data[i], data[j])/float64(n) -
				means[i]*means[j]
			v *= factor
			out.Set(i, j, v)
			if j > i {
				out.Set(j, i, v)
			}
		}
	}

	return out
}

/*
func CovMatrixWt(
    x *dense.Dense,
    wt []float64
) *dense.Dense {
    return &dense.Dense{}
    // TODO
    xx <- .colSums(wt * x, nrow(x), ncol(x), FALSE)
    xx <- x - rep(xx, each = nrow(x))
        # Each dimension (col) is shifted to being relative
        # to the corresponding weighted mean.
    crossprod(sqrt(wt) * xx) / (1 - sum(wt * wt))
        # FIXME: this line is a time consumer,
        #        taking 2/3 of the time in this function.
        # or
        # (1 - exp(logsumlog(2 * logwt.cov)))
        # if 'logwt.cov' is available.
        # But is it beneficial?
}
*/
