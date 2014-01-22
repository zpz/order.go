package stats

import (
	"github.com/zpz/matrix.go/dense"
	"math"
)

// Cov_ computes the covariance between Numeric slices x and y.
func Cov_(x, y Numeric) float64 {
	n := x.Len()
	assert(y.Len() == n, "lengths of x and y differ")

	m_x, m_y, m_xy := 0.0, 0.0, 0.0
	for i := 0; i < n; i++ {
		r := 1.0 / float64(i+1.0)
		vx := x.Get(i)
		vy := y.Get(i)
		m_x += (vx - m_x) * r
		m_y += (vy - m_y) * r
		m_xy += (vx*vy - m_xy) * r
	}

	cov := m_xy - m_x*m_y
	cov *= float64(n) / float64(n-1.0)
	return cov
}

// Cov_ computes the covariance between float64 slices x and y.
func Cov(x, y []float64) float64 {
	return Cov_(Float64Slice(x), Float64Slice(y))
}

// Cor_ computes the correlation coefficient between the Numeric slices
// x and y, returning the correlation coef and covariance between x and
// y, as well as the standard deviations of x and y.
func Cor_(x, y Numeric) (cor, cov, sd_x, sd_y float64) {
	n := x.Len()
	assert(y.Len() == n, "lengths of x and y differ")

	m_x, m_y, m_xx, m_yy, m_xy := 0.0, 0.0, 0.0, 0.0, 0.0
	for i := 0; i < n; i++ {
		r := 1.0 / float64(i+1.0)
		vx := x.Get(i)
		vy := y.Get(i)
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

// Cor_ computes the correlation coefficient between the float64 slices
// x and y, returning the correlation coef and covariance between x and
// y, as well as the standard deviations of x and y.
func Cor(x, y []float64) (cor, cov, sd_x, sd_y float64) {
	return Cor_(Float64Slice(x), Float64Slice(y))
}

// Extract columns of a matrix, and center each column.
func col_center(x *dense.Dense, sd, out []float64) {
	nrow, ncol := x.Dims()
	for icol, k := 0, 0; icol < ncol; icol++ {
		for irow := 0; irow < nrow; irow++ {
			out[k] = x.Get(irow, icol)
			k++
		}
		if sd == nil {
			Center(out[k-nrow:k], out[k-nrow:k])
		} else {
			s, m := Sd(out[k-nrow : k])
			sd[icol] = s
			Shift(out[k-nrow:k], -m, out[k-nrow:k])
		}
	}
}

// CovSlice returns the lower-triangular elements of the cov matrix
// between the columns of x.
func CovSlice(
	x *dense.Dense,
	out *LowerTri,
	workspace []float64) *LowerTri {

	nrow, ncol := x.Dims()

	if out == nil {
		out = NewLowerTri(ncol, nil)
	} else {
		assert(out.Dim() == ncol,
			"Provided 'out' is of wrong dimensionality")
	}

	workspace = use_slice(workspace, nrow*ncol)

	col_center(x, nil, workspace)

	for icol := 0; icol < ncol; icol++ {
		for irow := icol; irow < nrow; irow++ {
			out.Set(irow, icol,
				Dot(workspace[irow*nrow:(irow+1)*nrow],
					workspace[icol*nrow:(icol+1)*nrow])/
					float64(nrow-1))
		}
	}

	return out
}

// CorSlice returns the lower-triangular elements of the cor matrix
// between the columns of x.
func CorSlice(
	x *dense.Dense,
	out *LowerTri,
	workspace []float64) *LowerTri {

	nrow, ncol := x.Dims()

	if out == nil {
		out = NewLowerTri(ncol, nil)
	} else {
		assert(out.Dim() == ncol,
			"Provided 'out' is of wrong dimensionality")
	}

	workspace = use_slice(workspace, nrow*ncol)

	sd := make([]float64, ncol)

	col_center(x, sd, workspace)

	for icol := 0; icol < ncol; icol++ {
		for irow := icol; irow < nrow; irow++ {
			out.Set(irow, icol,
				Dot(workspace[irow*nrow:(irow+1)*nrow],
					workspace[icol*nrow:(icol+1)*nrow])/
					(float64(nrow-1)*sd[irow]*sd[icol]))
		}
	}

	return out
}

// CovMatrix returns the cov matrix between columns of x and columns of
// y. Element (i, j) of the cov matrix is the covariance between
// column i of x and column j of y.
func CovMatrix(
	x, y *dense.Dense,
	out *dense.Dense,
	workspace []float64) *dense.Dense {

	nrow_x, ncol_x := x.Dims()
	nrow_y, ncol_y := y.Dims()

	assert(nrow_x == nrow_y, "Number of rows of x and y mismatch")

	out = use_matrix(out, ncol_x, ncol_y)

	workspace = use_slice(workspace, nrow_x*ncol_x+nrow_y*ncol_y)
	work_x := workspace[0 : nrow_x*ncol_x]
	work_y := workspace[nrow_x*ncol_x : (nrow_x*ncol_x + nrow_y*ncol_y)]

	col_center(x, nil, work_x)
	col_center(y, nil, work_y)

	for ix := 0; ix < ncol_x; ix++ {
		for iy := 0; iy < ncol_y; iy++ {
			out.Set(ix, iy,
				Dot(work_x[ix*nrow_x:(ix+1)*nrow_x],
					work_y[iy*nrow_y:(iy+1)*nrow_y])/
					float64(nrow_x-1))
		}
	}

	return out
}

// CorMatrix returns the cor matrix between columns of x and columns of
// y. Element (i, j) of the cor matrix is the correlation between
// column i of x and column j of y.
func CorMatrix(
	x, y *dense.Dense,
	out *dense.Dense,
	workspace []float64) *dense.Dense {

	nrow_x, ncol_x := x.Dims()
	nrow_y, ncol_y := y.Dims()

	assert(nrow_x == nrow_y, "Number of rows of x and y mismatch")

	out = use_matrix(out, ncol_x, ncol_y)

	workspace = use_slice(workspace, nrow_x*ncol_x+nrow_y*ncol_y)
	work_x := workspace[0 : nrow_x*ncol_x]
	work_y := workspace[nrow_x*ncol_x : (nrow_x*ncol_x + nrow_y*ncol_y)]

	sd := make([]float64, ncol_x+ncol_y)
	sd_x := sd[:ncol_x]
	sd_y := sd[ncol_x:]

	col_center(x, sd_x, work_x)
	col_center(y, sd_y, work_y)

	for ix := 0; ix < ncol_x; ix++ {
		for iy := 0; iy < ncol_y; iy++ {
			out.Set(ix, iy,
				Dot(work_x[ix*nrow_x:(ix+1)*nrow_x],
					work_y[iy*nrow_y:(iy+1)*nrow_y])/
					(float64(nrow_x-1)*sd_x[ix]*sd_y[iy]))
		}
	}

	return out
}
