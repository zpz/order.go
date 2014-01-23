package stats

import (
	"github.com/pmylund/sortutil"
	"github.com/zpz/matrix.go/dense"
	"log"
	"math"
	"sort"
)

type NormixKind int

const (
	UniCovMix NormixKind = iota
	// All mixture components have the same cov matrix.
	ScaledCovMix
	// The cov matrices of all mixture components
	// are scaled versions of the same cov matrix.
	// Each mixture component has its own scaling factor.
	FreeCovMix
	// Each mixture component has a distinct cov matrix.
)

// Normix defines a normal mixture distribution.
type Normix struct {
	kind NormixKind

	// Log weight of each mixture component.
	// Length is the number of mixture components.
	logweight []float64

	// Each row is the mean of one mixture component.
	// Each col is one dimension.
	mean *dense.Dense

	// Each row is the lower-triangular (including diagonal) elements
	// of the cov matrix of one mixture component,
	// ordered in column major.
	// One row per mixture component.
	// The number of columns is equal to
	// (ndim * ndim + ndim)/2
	cov *dense.Dense

	// If kind is UniCovMix, cov has 1 row, cov_scale is nil.
	// If kind is ScaledCovMix, cov has 1 row, cov_scale contains
	// scaling factor of the cov matrix for each mixture component.
	// If kind is FreeCovMix, cov has as many rows as there are mixture
	// components, and cov_scale is nil.
	cov_scale []float64
}

// NewNormix creates a Normix object and returns a pointer to it.
// It allocates memory for all fields,
// which will be populated later.
// The allocated memory ensures that the Normix object
// holds its own data---its slice members do not point
// to data outside of the object.
func NewNormix(n_dim, n_mix int, kind NormixKind) *Normix {
	if n_dim < 1 || n_mix < 1 {
		panic("NewNormix: positive arguments expected for n_dim and n_mix")
	}

	var mix Normix

	mix.kind = kind
	mix.logweight = make([]float64, n_mix)
	mix.mean = dense.NewDense(n_mix, n_dim)
	mix.cov_scale = nil

	switch kind {
	case UniCovMix:
		mix.cov = dense.NewDense(1, n_dim*(n_dim+1)/2)
	case ScaledCovMix:
		// Usually this is not used when n_dim is 1,
		// but it is allowed.
		mix.cov = dense.NewDense(1, n_dim*(n_dim+1)/2)
		mix.cov_scale = make([]float64, n_mix)
	case FreeCovMix:
		mix.cov = dense.NewDense(n_mix, (n_dim*(n_dim+1))/2)
	default:
		panic("NewNormix: unrecognized value for kind")
	}

	return &mix
}

// Kind returns the kind of the normal mixture.
func (mix *Normix) Kind() NormixKind {
	return mix.kind
}

// Dim returns the number of dimensions of the normal mixture distribution.
func (mix *Normix) Dim() int {
	return mix.mean.Cols()
}

// Size returns the number of mixture components in the normal mixture
// distribution.
func (mix *Normix) Size() int {
	return len(mix.logweight)
}

/*
// LogweightSlice returns the slice of log-weights
// of mixture components in the distribution.
// One use case is to set the weights via the returned slice.
func (mix *Normix) LogweightSlice() []float64 {
    return mix.logweight
}



// MeanSlice returns the slice of the mean vector
// for the mixture component at the specified index.
// One use case is to set the mean vector via
// the returned slice.
func (mix *Normix) MeanSlice(i_mix int) []float64 {
    if i_mix < 0 || i_mix >= mix.NMix() {
        panic(fmt.Errorf("Normix.MeanSlice: index out of range"))
    }
    return mix.mean.RowView(i_mix)
}



// CovSlice returns the slice of the covariances.
// Note that the meaning of this slice differs according
// to the kind of the norm mixture.
func (mix *Normix) CovSlice(i_mix int) []float64 {
    if i_mix < 0 || i_mix >= mix.NMix() {
        panic(fmt.Errorf("Normix.CovSlice: i_mix out of range"))
    }
    if mix.kind == UniCovMix || mix.kind == ScaledCovMix {
        return mix.cov.RowView(0)      // ignoring 'i_mix'
    }
    return mix.cov.RowView(i_mix)
}



// CovScaleSlice returns the cov_scale member of
// a Normix. This is allowed only when the kind
// is ScaledCovMix.
func (mix *Normix) CovScaleSlice() []float64 {
    if mix.kind != ScaledCovMix {
        panic(fmt.Errorf("the kind of the normal mixture is not ScaledCovMIx"))
    }
    return mix.cov_scale
}


// Not complete. Just output some info for now.
func (mix *Normix) String() string {
	ndim := mix.Dim()
	nmix := mix.Size()
	s := fmt.Sprintf(
		"Normal mixture density\n"+
			"  dimensionality: %d\n"+
			"  number of mixture components: %d\n",
		ndim, nmix)
	return s
}
*/

// Pooled mean and variance of a random vector represented by a
// mixture---
//   Let there be 'n' mixtures for a k-vector
//   with weights a_i and component means m_i
//   and component cov matrices v_i,
//   then the overall mean and cov are
//     m = a_1 m_1 + ... + a_n m_n
//     v = a_1(v_i + m_i m_i^T) + ... + a_n(v_n + m_n m_n^T) - m m^T
//
// This function does not assume normality of the mixtures.

/*
summary.normmix <- function(object, ...)
{
    n.mix <- length(object$logweight)
    p <- ncol(object$mean)

    if (p > 1)
    {
        if (nrow(object$mean) > 1)
            mm <- .colSums(
                        object$mean * exp(object$logweight),
                        nrow(object$mean),
                        ncol(object$mean),
                        FALSE
                        )
        else
            mm <- object$mean[1, ]

        if (nrow(object$cov) > 1)
        {
            lower.tri.idx <- lower.tri(matrix(0, p, p), diag = TRUE)
            my.mean <- object$mean[1, ]

            cc <- rep(0, ncol(object$cov))
            for (i in 1 : n.mix)
            {
                if (nrow(object$mean) > 1)
                    my.mean <- object$mean[i, ]
                cc <- cc + exp(object$logweight[i]) * (object$cov[i, ] +
                    tcrossprod(my.mean, my.mean)[lower.tri.idx])
            }

            cc <- matrix(cc[inv.lower.tri(p)], p)
            cc <- cc - mm %*% t(mm)
        } else
        {
            cc <- matrix(object$cov[1, inv.lower.tri(p)], p)
            if (nrow(object$mean) > 1)
            {
                for (i in 1 : nrow(object$mean))
                    cc <- cc + exp(object$logweight[i]) *
                        tcrossprod(object$mean[i, ], object$mean[i, ])
                cc <- cc - mm %*% t(mm)
            }
        }
    } else
    {
        if (length(object$mean) > 1)
            mm <- sum(c(object$mean) * exp(object$logweight))
        else
            mm <- c(object$mean)

        cc <- sum(exp(object$logweight) *
                (c(object$cov) + c(object$mean)^2))
        cc <- cc - mm^2
    }

    list(mean = mm, cov = cc)
}
*/

// Density computes the pdf of each row of x in the Normix mix.
func (mix *Normix) Density(x *dense.Dense, out []float64) []float64 {
	ndim := mix.Dim()
	assert(x.Cols() == ndim, "Wrong shape for input x")

	nx := x.Rows()

	zz := mix.density_stats(x, nil)

	nmix := mix.Size()

	if nmix == 1 {
		out = zz.GetData(out)
	} else {
		for imix := 0; imix < nmix; imix++ {
			Shift(zz.RowView(imix), mix.logweight[imix], zz.RowView(imix))
		}
		lw := make([]float64, nmix)
		for ix := 0; ix < nx; ix++ {
			out[ix] = LogSumExp(zz.GetCol(ix, lw))
		}
	}

	return Transform(out, math.Exp, out)
}

// Random generates n random samples from the normal mixture
// distribution and returns the sample in a slice,
// one case after another.
func (mix *Normix) Random(n int, out *dense.Dense) *dense.Dense {
	ndim := mix.Dim()
	out = use_matrix(out, n, ndim)

	mixidx := LogweightedSample(mix.logweight, n, nil)
	sort.Ints(mixidx)

	cov_mat := dense.NewDense(ndim, ndim)

	switch mix.kind {
	case UniCovMix:
		mvn := NewMvnormal(make([]float64, ndim),
			symm_fill(mix.cov.RowView(0), cov_mat))

		mvn.Random(n, &RNG{}, out)

		for i, idx := range mixidx {
			z := out.RowView(i)
			Add(z, mix.mean.RowView(idx), z)
		}

	case ScaledCovMix:
		mvn := NewMvnormal(make([]float64, ndim),
			symm_fill(mix.cov.RowView(0), cov_mat))

		mvn.Random(n, &RNG{}, out)

		for i, idx := range mixidx {
			z := out.RowView(i)
			Scale(z, math.Sqrt(mix.cov_scale[idx]), z)
			Add(z, mix.mean.RowView(idx), z)
		}

	case FreeCovMix:
		for iz := 0; iz < n; {
			imix := mixidx[iz]
			nmix := 1
			iz++
			for ; iz < n && mixidx[iz] == imix; iz++ {
				nmix++
			}
			// When the weights of the mixture components are
			// highly uneven, it's common that a small number
			// of the high-weight components are repeatedly
			// used in generating samples.
			// We find the number of samples each mixture
			// component generates, so that these samples
			// are generated at once.

			mvn := NewMvnormal(mix.mean.RowView(imix),
				symm_fill(mix.cov.RowView(imix), cov_mat))
			z := out.SubmatrixView(iz-nmix, 0, nmix, ndim)
			mvn.Random(nmix, &RNG{}, z)
		}
	}

	return out
}

// Marginal returns the marginal distribution, as a Normix,
// of the dimensions specified by dims.
func (mix *Normix) Marginal(dims []int) *Normix {
	ndim := len(dims)
	nmix := mix.Size()
	kind := mix.kind

	assert(ndim > mix.Dim(), "too many elements in dims")

	// TODO: check that dims does not contain duplicate elements.

	z := NewNormix(ndim, nmix, kind)

	// Get z.logweight.
	copy(z.logweight, mix.logweight)

	// Get z.mean.
	for imix := 0; imix < nmix; imix++ {
		pick_float64s(mix.mean.RowView(imix), dims, z.mean.RowView(imix))
	}

	// Get z.cov.
	covidx := lowertri_subsetter(mix.Dim(), dims)
	if kind == FreeCovMix {
		for imix := 0; imix < nmix; imix++ {
			pick_float64s(mix.cov.RowView(imix), covidx, z.cov.RowView(imix))
		}
	} else {
		pick_float64s(mix.cov.RowView(0), covidx, z.cov.RowView(0))
	}

	// Get z.cov_scale.
	if kind == ScaledCovMix {
		copy(z.cov_scale, mix.cov_scale)
	}

	return z
}

// Conditional derives conditional density,
// given observed values for some dimensions in a normal mixture,
// of the other dimensions as another normal mixture.
func (mix *Normix) Conditional(
	data []float64,
	// Data vector.
	dims []int,
	// Values in data correspond to these dimensions in mix.
	wt_tol float64,
	// In the resultant mixture density,
	// components with highest weights that collectively
	// account for (1 - wt_tol) of the total weight are kept.
	// Use 0 if you don't know better.
) *Normix {

	assert(len(data) == len(dims), "dimensionality mismatch")

	//==================================================
	// Calculate likelihoods of the data in its marginal
	// distribution.

	marginal_pdf := mix.Marginal(dims)
	loglikely := marginal_pdf.density_stats(
		dense.DenseView(data, 1, len(dims)), nil)

	//========================================
	// Update weight of each mixture component
	// to take into account the likelihoods.

	logwt := Add(loglikely.DataView(), mix.logweight, nil)
	logintlikely := LogSumExp(logwt)
	// Log integrated likelihood.
	Shift(logwt, -logintlikely, logwt)
	// Normalized; now sum(exp(logwt)) = 1.

	// Screen the mixture components and discard those
	// with negligible weights.

	var idx_keep []int

	if len(logwt) > 1 && wt_tol > 0 {
		idx_keep := lose_weight(logwt, 1-wt_tol)

		if len(idx_keep) < len(logwt) {
			logwt = pick_float64s(logwt, idx_keep, nil)
			total_wt := math.Exp(LogSumExp(logwt))

			log.Println("keeping",
				len(idx_keep), "of", mix.Size(),
				"components for a total weight of",
				total_wt)

			Shift(logwt, -LogSumExp(logwt), logwt)
			// Normalize so that weights sum to 1.
		}
	} else {
		idx_keep = make([]int, len(logwt))
		for i, _ := range idx_keep {
			idx_keep[i] = i
		}
	}

	//====================================
	// Compute conditional mean and cov.

	n_y := len(dims)
	dims_y := dims
	// 'y' indicates the conditioning dimensions and data.

	n_x := mix.Dim() - n_y
	dims_x := exclude(mix.Dim(), dims_y)
	// 'x' indicates the conditioned, i.e. target,
	// dimensions.

	mix_x := NewNormix(n_x, len(logwt), mix.kind)
	// Conditional distribution.
	copy(mix_x.logweight, logwt)

	sigma_x_covslice_idx := lowertri_subsetter(n_x+n_y, dims_x)

	sigma_y := dense.NewDense(n_y, n_y)
	// Cov matrix between dimensions dims_y.
	sigma_y_covslice_idx := lowertri_subsetter(n_x+n_y, dims_y)
	sigma_y_covslice := make([]float64, len(sigma_y_covslice_idx))
	// This holds the lower_tri elements of sigma_y.

	sigma_yx := dense.NewDense(n_y, n_x)
	// Cov matrix between dims_y and dims_x.
	sigma_yx_slice_idx := make([]int, n_x*n_y)
	// row by row.
	for k, row := 0, 0; row < n_y; row++ {
		for col := 0; col < n_x; col++ {
			sigma_yx_slice_idx[k] = lowertri_ij2idx(
				n_x+n_y, dims_y[row], dims_x[col])
			k++
		}
	}

	// workspaces
	A := dense.NewDense(n_y, n_x)

	if mix.kind == FreeCovMix {
		mu_y := make([]float64, n_y)

		for idx_new, idx_old := range idx_keep {
			mu := mix.mean.RowView(idx_old)
			sigma := mix.cov.RowView(idx_old)
			mu_x := mix_x.mean.RowView(idx_new)
			sigma_x := mix_x.cov.RowView(idx_new)

			pick_float64s(mu, dims_y, mu_y)

			pick_float64s(sigma,
				sigma_y_covslice_idx, sigma_y_covslice)
			symm_fill(sigma_y_covslice, sigma_y)

			pick_float64s(sigma,
				sigma_yx_slice_idx, sigma_yx.DataView())

			conditional_normal(
				data, mu_y,
				sigma_y, sigma_yx, A,
				mu_x, sigma_x)

			for i, idx := range dims_x {
				mu_x[i] += mu[idx]
			}
			for i, idx := range sigma_x_covslice_idx {
				sigma_x[i] += sigma[idx]
			}
		}

	} else {
		if mix_x.kind == ScaledCovMix {
			pick_float64s(mix.cov_scale, idx_keep, mix_x.cov_scale)
		}

		sigma := mix.cov.RowView(0)
		sigma_x := mix_x.cov.RowView(0)

		pick_float64s(sigma,
			sigma_y_covslice_idx, sigma_y_covslice)
		symm_fill(sigma_y_covslice, sigma_y)

		pick_float64s(sigma,
			sigma_yx_slice_idx, sigma_yx.DataView())

		mu_x_delta := make([]float64, n_x)

		conditional_normal(
			data, make([]float64, n_y),
			sigma_y, sigma_yx, A,
			mu_x_delta, sigma_x)

		for i, idx := range sigma_x_covslice_idx {
			sigma_x[i] += sigma[idx]
		}

		mu_y := make([]float64, n_y)

		for idx_new, idx_old := range idx_keep {
			mu := mix.mean.RowView(idx_old)
			mu_x := mix_x.mean.RowView(idx_new)

			pick_float64s(mu, dims_y, mu_y)
			dense.Mult(dense.DenseView(mu_y, 1, n_y), A,
				dense.DenseView(mu_x, 1, n_x))

			Subtract(mu_x_delta, mu_x, mu_x)

			for i, idx := range dims_x {
				mu_x[i] += mu[idx]
			}
		}
	}

	return mix_x
}

// density_stats computes the log-density of the data in each mixture component.
// Input x contains one 'observation' or 'case' per row.
// Upon return, row i of out contains log-densities of every row of x
// in the i-th mixture component of mix.
func (mix *Normix) density_stats(x *dense.Dense, out *dense.Dense) *dense.Dense {
	ndim, nmix := mix.Dim(), mix.Size()
	nx := x.Rows()

	assert(x.Cols() == ndim, "input x has wrong shape")
	out = use_matrix(out, nmix, nx)

	if ndim == 1 {
		var xx []float64
		if x.Contiguous() {
			xx = x.DataView()
		} else {
			xx = x.GetData(nil)
		}
		switch mix.kind {
		case UniCovMix:
			sd := math.Sqrt(mix.cov.Get(0, 0))
			for imix := 0; imix < nmix; imix++ {
				NewNormal(mix.mean.Get(imix, 0), sd).
					Density(xx, out.RowView(imix))
			}
		case ScaledCovMix:
			v := mix.cov.Get(0, 0)
			for imix := 0; imix < nmix; imix++ {
				NewNormal(mix.mean.Get(imix, 0),
					math.Sqrt(v*mix.cov_scale[imix])).
					Density(xx, out.RowView(imix))
			}
		case FreeCovMix:
			for imix := 0; imix < nmix; imix++ {
				NewNormal(mix.mean.Get(imix, 0),
					math.Sqrt(mix.cov.Get(imix, 0))).
					Density(xx, out.RowView(imix))
			}
		}
		return out
	}

	cov_mat := dense.NewDense(ndim, ndim)

	if mix.kind == FreeCovMix {
		for imix := 0; imix < nmix; imix++ {
			mvn := NewMvnormal(mix.mean.RowView(imix),
				symm_fill(mix.cov.RowView(imix), cov_mat))
			mvn.Density(x, out.RowView(imix))
		}
	} else {
		// Create a mvn with zero mean.
		mvn := NewMvnormal(make([]float64, ndim),
			symm_fill(mix.cov.RowView(0), cov_mat))

		// Subtract mean from all data so that their densities
		// are computed using the zero-mean distribution above.
		xx := dense.NewDense(nx*nmix, ndim)
		for imix := 0; imix < nmix; imix++ {
			xxview := xx.SubmatrixView(imix*nx, 0, nx, ndim)
			dense.Copy(xxview, x)
			for row := 0; row < nx; row++ {
				Subtract(xxview.RowView(row), mix.mean.RowView(imix),
					xxview.RowView(row))
			}
		}

		if mix.kind == UniCovMix {
			if out.Contiguous() {
				mvn.Density(xx, out.DataView())
			} else {
				out.SetData(mvn.Density(xx, nil))
			}
		} else { // ScaledCovMix
			symm_fill(mix.cov.RowView(0), cov_mat)
			dist := Mahalanobis(xx, make([]float64, ndim), cov_mat, nil)
			cov_det := cov_mat.Det()
			for imix := 0; imix < nmix; imix++ {
				res := out.RowView(imix)
				cov_scale := mix.cov_scale[imix]
				coef := 1.0 / math.Sqrt(
					math.Pow(2*math.Pi*cov_scale, float64(ndim))*
						cov_det)
				Scale(dist[imix*nx:(imix+1)*nx], -0.5/cov_scale, res)
				Transform(res, math.Exp, res)
				Scale(res, coef, res)
			}
		}
	}

	return out
}

// lose_weight takes a slice of log weights
// and returns the indices of the high weights such that
// their total reaches a specified fraction of the total of all weights.
// The intended use case is to discard a large number of low-weight
// elements in a highly skewed weighted sample.
// The returned indices are in order of decreasing weights except
// when keep_total == 1, in which case the return is simply
// the full list of indices.
func lose_weight(
	logweight []float64,
	// Log weights; does not need to be normalized.
	args ...float64,
	// Up to two optional arguments;
	// the first would be keep_total and the second would be keep_tol.
	// keep_total: high weights that account for at least this fraction
	// of the total weights will be kept.
	// keep_tol: keep all weights that are larger than this fraction of
	// the max weight in the input.
	// As a result, if all weights are equal, then the entire sample is
	// retained.
	// Example values for keep_total and keep_tol may be 0.99 and 0.1,
	// respectively.
) []int {

	keep_total := 0.99
	keep_tol := 0.1
	if len(args) > 0 {
		keep_total = args[0]
		if len(args) > 1 {
			keep_tol = args[1]
		}
	}

	if keep_total >= 0.999999 {
		// Simply return all.
		z := make([]int, len(logweight))
		for i := range z {
			z[i] = i
		}
		return z
	}

	// Make a copy; don't modify the input.
	logwt := append([]float64{}, logweight...)

	// Normalize, such that sum(exp(logwt)) == 1.
	Shift(logwt, -LogSumExp(logwt), logwt)

	// Get indices listing the weights in descending order.
	idx_ordered := Order(logwt, nil)
	sortutil.Reverse(idx_ordered)

	n := len(logwt)

	// Get number of entries to keep according to keep_total.
	k := 1
	sum := 0.0
	for _, idx := range idx_ordered {
		sum += math.Exp(logwt[idx])
		if sum < keep_total {
			k++
		} else {
			break
		}
	}
	if k > n {
		k = n
	}

	// Get additional number of entries to keep according to keep_tol.
	if k < n {
		bar := math.Log(keep_tol) + logwt[idx_ordered[0]]
		for k < n {
			if logwt[idx_ordered[k]] > bar {
				k++
			} else {
				break
			}
		}
	}

	if k > n/2 {
		return idx_ordered[:k]
	}

	idx := make([]int, k)
	copy(idx, idx_ordered[:k])
	return idx
}

// exclude returns integers between 0 (inclusive) and
// n (exclusive), excluding those in index.
// Does not assume elements in include are sorted
// (otherwise the code can be more efficient).
func exclude(n int, index []int) []int {
	idx := make([]int, n)
	for _, i := range index {
		idx[i] = 1
	}
	k := 0
	for j := 0; j < n; j++ {
		if idx[j] == 0 {
			idx[k] = j
			k++
		}
	}
	return idx[:k]
}

// pick_float64s returns a subslice with the elements
// at the specified indices.
func pick_float64s(x []float64, index []int, y []float64) []float64 {
	y = use_slice(y, len(index))
	for i, j := range index {
		y[i] = x[j]
	}
	return y
}

// lowertri_subsetter takes the dimensionality of a square matrix,
// and outputs the indices into its col-major on- and below-diagonal
// element slice that would produce a similar lower-triangular
// slice for the specified dimensions.
func lowertri_subsetter(ndim_in int, dims []int) []int {
	ndim_out := len(dims)
	n_out := (ndim_out * (ndim_out + 1)) / 2
	out := make([]int, n_out)

	k := 0 // index in out to write next
	for col := 0; col < ndim_out; col++ {
		for row := col; row < ndim_out; row++ {
			out[k] = lowertri_ij2idx(ndim_in, dims[row], dims[col])
			k++
		}
	}

	return out
}

// lowertri_ij2idx returns the index of the element
// in a lower-triangular slice that represents the element (row, col)
// in the full matrix.
func lowertri_ij2idx(ndim, row, col int) int {
	// Total number of elements up to, but not including,
	// column j in the slice:
	//        n + (n-1) +...+ (n-j+1)
	//      = (n + n-j+1)/2 * (n - (n-j+1) + 1)
	//      = j * (n + n - j + 1) / 2
	//      = j * n - j * (j-1) / 2
	//
	// Check:
	//   j = 0: --> 0
	//   j = 1: --> n
	//   j = 2: --> n + n - 1
	//   j = 3: --> 3 * (n - 1) = n + (n-1) + (n-2)

	if row < col {
		row, col = col, row
	}

	return col*ndim - (col*(col-1))/2 + row - col
}

// symm_fill fills the symmetric matrix mat with data, which contains
// mat's on- and below-diagonal elements in col major.
func symm_fill(data []float64, mat *dense.Dense) *dense.Dense {
	// len(data) = (n*n + n)/2 = n * (n+1) / 2
	n := int(math.Sqrt(2.0 * float64(len(data))))
	mat = use_matrix(mat, n, n)
	for k, i := n, 0; k > 0; k-- {
		copy(mat.RowView(i)[i:], data[:k])
		for j := 1; j < k; j++ {
			mat.Set(n-k+j, i, data[j])
		}
		i++
		data = data[k:]
	}
	return mat
}

func conditional_normal(
	y []float64, // data
	mu_y []float64, // mean
	sigma_yy, sigma_yx *dense.Dense, // cov_yy, cov_yx
	A *dense.Dense, // workspace
	mu_x_delta, sigma_xx_delta []float64,
	// Conditional mean of x is original mean of x
	// plus mu_x_delta;
	// conditional cov of x (the lower-tri elements)
	// is original cov of x plus sigma_xx_delta.
) {
	// mu_{x|y} = mu_x + S_{xy} Inv(S_{yy}) (y - mu_y)
	// S_{x|y} = S_{xx} - S_{xy} Inv(S_{yy}) S_{yx}

	dimx := sigma_yx.Cols()
	dimy := len(mu_y)
	assert(
		sigma_yy.Rows() == dimy &&
			sigma_yy.Cols() == dimy &&
			sigma_yx.Rows() == dimy &&
			len(y) == dimy &&
			len(mu_x_delta) == dimx &&
			len(sigma_xx_delta) == (dimx*dimx+dimx)/2,
		"dimensions mismatch")

	A = use_matrix(A, dimy, dimx)
	dy := make([]float64, dimy)
	ddy := make([]float64, dimy)
	// workspaces

	dense.Copy(A, sigma_yx)
	// TODO: if sigma_yx is not re-used outside of this function,
	// we may avoid this copy.

	A = dense.Solve(sigma_yy, A)
	// TODO: make use of Cholesky? b/c sigma_yy is pd.

	Subtract(y, mu_y, dy)
	dense.Mult(dense.DenseView(dy, 1, dimy), A,
		dense.DenseView(mu_x_delta, 1, dimx))

	k := 0
	for col := 0; col < dimx; col++ {
		A.GetCol(col, dy) // re-use the space 'dy'
		for row := col; row < dimx; row++ {
			// Need the element (row, col) of sigma_xy %*% A,
			// which is the dot product of the col-th column of A
			// and the row-th column of sigma_yx.
			sigma_yx.GetCol(row, ddy)
			sigma_xx_delta[k] = -Dot(dy, ddy)
		}
	}
}
