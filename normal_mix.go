package stats

import (
    "sort"
    "fmt"
    "log"
    "math"
    "github.com/gonum/floats"
    "github.com/gonum/matrix/mat64"
    "github.com/pmylund/sortutil"
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
// len(logweight) reflects the number of mixture components.
// len(mean) / len(logweight) reflects the dimensionality.
type Normix struct {
    kind        NormixKind
    logweight   []float64
        // nmix = len(logweight)
    mean        []float64
        // ndim = len(mean) / nmix
    cov         []float64
        // Lower-triangular (including diagonal) elements,
        // ordered in column major, of the covariance matrix
        // of each mixture component, one component after another.
        // For each mixture component, there are
        // n = (ndim * ndim + ndim)/2 elements.
    cov_scale   []float64
        // If kind is UniCovMix,
        //      len(cov) is n, len(cov_scale) is 0;
        // if kind is ScaledCovMix,
        //      len(cov) is n, len(cov_scale) is nmix;
        // if kind is FreeCovMix,
        //      len(cov) is n * ndim, len(cov_scale) is 0.
}



// NewNormix creates a Normix.
// It allocates memory for all fields,
// which will be populated later.
// The allocated memory ensures that the Normix object
// holds its own data---its slice members do not point
// to data outside of the object.
func NewNormix(n_dim, n_mix int, kind NormixKind) *Normix {
    if n_dim < 1 || n_mix < 1 {
        panic(fmt.Errorf("NewNormix: positive arguments expected for n_dim and n_mix"))
    }

    mix := Normix{
        kind:       kind,
        logweight:  make([]float64, n_mix),
        mean:       make([]float64, n_dim * n_mix)}

    switch kind {
    case UniCovMix:
        mix.cov = make([]float64, n_dim * (n_dim + 1) / 2)
    case ScaledCovMix:
            // Usually this is not used when n_dim is 1,
            // but it is allowed.
        mix.cov = make([]float64, n_dim * (n_dim + 1) / 2)
        mix.cov_scale = make([]float64, n_mix)
    case FreeCovMix:
        mix.cov = make([]float64, n_mix * (n_dim * (n_dim + 1)) / 2)
    default:
        panic(fmt.Errorf("Unrecognized value for kind"))
    }

    return &mix
}




// Kind returns the kind of the normal mixture.
func (mix *Normix) Kind() NormixKind {
    return mix.kind
}



// NDim returns the number of dimensions of the distribution.
func (mix *Normix) NDim() int {
    return len(mix.mean) / mix.NMix()
}



// NMix returns the number of mixture components.
func (mix *Normix) NMix() int {
    return len(mix.logweight)
}




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
    ndim := mix.NDim()
    return mix.mean[ndim * i_mix : ndim * (i_mix + 1)]
}



// CovSlice returns the slice of the covariances.
// Note that the meaning of this slice differs according
// to the kind of the norm mixture.
func (mix *Normix) CovSlice(i_mix int) []float64 {
    if i_mix < 0 || i_mix >= mix.NMix() {
        panic(fmt.Errorf("Normix.CovSlice: i_mix out of range"))
    }

    if mix.kind == UniCovMix || mix.kind == ScaledCovMix {
        return mix.cov      // ignoring 'i_mix'
    }

    ndim := mix.NDim()
    n := ndim * (ndim + 1) / 2
    return mix.cov[n * i_mix : n * (i_mix + 1)]
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
    ndim := mix.NDim()
    nmix := mix.NMix()
    s := fmt.Sprintf(
        "Normal mixture density\n" +
        "  dimensionality: %d\n" +
        "  number of mixture components: %d\n",
        ndim, nmix)
    if nmix > 1 {
        if nmix <= 10 {
            w := append([]float64{}, mix.logweight...)
            floats.Apply(math.Exp, w)
            s += fmt.Sprintf(
                "  component weights: %v\n",
                w)
        } else {
            //s += fmt.Sprintf(
            //    "  entropy of component weights: %v\n",
            //    weight_entropy(mix.logweight, true, "entropy"))
                // FIXME: implement weight_entropy
        }
    }
    return s
}


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




// Density computes the pdf of x in the Normix mix.
// Each row of x is an 'observation'.
// Return density for each observation.
func (mix *Normix) Density(x *mat64.Dense, out []float64) []float64 {
    ndim := mix.NDim()
    assert(x.Ncols() == ndim, "Wrong shape for input x")

    nx := x.Nrows()

    zz := mix.density_stats(x, nil)

    nmix := mix.NMix()

    if nmix == 1 {
        copy(out, zz)
    } else {
        lw := make([]float64, nmix)
        for ix := 0; ix < nx; ix++ {
            for imix := 0; imix < nmix; imix++ {
                lw[imix] = zz[imix * nx + ix] + mix.logweight[imix]
            }
            out[ix] = floats.LogSumExp(lw)
        }
    }

    return Transform(out, math.Exp, out)
}




// Random generates n random samples from the normal mixture
// distribution and returns the sample in a slice,
// one case after another.
func (mix *Normix) Random(n int, out *mat64.Dense) *mat64.Dense {
    ndim := mix.NDim()
    out = use_matrix(out, n, ndim)

    mixidx := LogweightedSample(mix.logweight, n, nil)
    sort.Ints(mixidx)

    switch mix.kind {
    case UniCovMix:
        mvn := NewMvnormal(make([]float64, ndim),
            NewLowerTri(ndim, mix.cov).Expand(nil))

        mvn.Random(n, &RNG{}, out)

        for i, idx := range mixidx {
            z := out.RowView(i)
            Add(z, mix.MeanSlice(idx), z)
        }

    case ScaledCovMix:
        mvn := NewMvnormal(make([]float64, ndim),
            NewLowerTri(ndim, mix.cov).Expand(nil))

        mvn.Random(n, &RNG{}, out)

        for i, idx := range mixidx {
            z := out.RowView(i)
            Scale(z, math.Sqrt(mix.cov_scale[idx]), z)
            Add(z, mix.MeanSlice(idx), z)
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

            mvn := NewMvnormal(mix.MeanSlice(imix),
                NewLowerTri(ndim, mix.CovSlice(imix)).Expand(nil))
            z := out
            z.View(iz - nmix, 0, nmix, ndim)
            mvn.Random(nmix, &RNG{}, z)
        }
    default:
        panic(fmt.Errorf("unrecognized kind for Normix"))
    }

    return out
}




// Marginal returns the marginal distribution, as a Normix,
// of the dimensions specified by dims.
func (mix *Normix) Marginal(dims []int) *Normix {
    ndim := len(dims)
    nmix := mix.NMix()
    kind := mix.kind

    z := NewNormix(ndim, nmix, kind)

    copy(z.logweight, mix.logweight)

    for imix := 0; imix < nmix; imix++ {
        pick_float64s_fill(mix.MeanSlice(imix), dims, z.MeanSlice(imix))
    }


    covidx := covslice_subset_xx_idx(mix.NDim(), dims)

    if kind == FreeCovMix {
        for imix := 0; imix < nmix; imix++ {
            pick_float64s_fill(
                mix.CovSlice(imix), covidx, z.CovSlice(imix))
        }
    } else {
        pick_float64s_fill(mix.cov, covidx, z.cov)
        if kind == ScaledCovMix {
            copy(z.cov_scale, mix.cov_scale)
        }
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
        // Indices of the elements of data in the
        // multivariate variable whose density is described by mix.
    wt_tol float64,
        // In the resultant mixture density,
        // components with highest weights that collectively
        // account for (1 - wt_tol) of the total weight are kept.
        // Use 0 if you don't know better.
    ) *Normix {


    //==================================================
    // Calculate likelihoods of the data in its marginal
    // distribution.

    marginal_pdf := mix.Marginal(dims)
    loglikely := marginal_pdf.density_stats(
        mat64.NewDense(1, len(dims), data), nil)


    //========================================
    // Update weight of each mixture component
    // to take into account the likelihoods.

    logwt := append([]float64{}, loglikely...)
    floats.Add(logwt, mix.logweight)
    logintlikely := floats.LogSumExp(logwt)
        // Log integrated likelihood.
    floats.AddConst(-logintlikely, logwt)
        // Normalized; now sum(exp(logwt)) = 1.

    // Screen the mixture components and discard those
    // with negligible weights.

    var idx_keep []int;

    if len(logwt) > 1 && wt_tol > 0 {
        idx_keep := lose_weight(logwt, 1 - wt_tol)

        if len(idx_keep) < len(logwt) {
            logwt = pick_float64s(logwt, idx_keep)
            total_wt := math.Exp(floats.LogSumExp(logwt))

            log.Println("keeping",
                len(idx_keep), "of", mix.NMix(),
                "components for a total weight of",
                total_wt)

            floats.AddConst(-floats.LogSumExp(logwt), logwt)
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

    n_x := mix.NDim() - n_y
    dims_x := exclude(mix.NDim(), dims_y)
        // 'x' indicates the conditioned, i.e. target,
        // dimensions.

    mix_x := NewNormix(n_x, len(logwt), mix.kind)
        // Conditional distribution.
    copy(mix_x.logweight, logwt)


    y_delta := make([]float64, n_y)
        // Used to hold y_data - y_mean
        // as well as some intermediate results.

    sigma_x_covslice_idx := covslice_subset_xx_idx(n_x + n_y, dims_x)
        // This holds the indices out of the xy cov slice
        // the produce the cov slice for dimensions dims_x.
    sigma_x := mat64.NewDense(n_x, n_x, nil)
        // All zeros.
        // Cov matrix between dimensions dims_x.

    sigma_y_slice := make([]float64, n_y * n_y)
        // This holds the elements for the cov matrix between
        // dimensions dims_y.
    sigma_y_covslice_idx := covslice_subset_xx_idx(n_x + n_y, dims_y)
        // This holds the indices out of the xy cov slice
        // the produce the cov slice for dimensions dims_y.
    sigma_y_covslice := make([]float64, (n_y * n_y + n_y)/2)
        // This holds the covslice that would be
        // created by sigma_y_covslice_idx and expanded into sigma_y_slice.

    sigma_yx_slice := make([]float64, n_y * n_x)
        // This holds the elements of the cov matrix
        // between dimensions dims_y and dims_x.
    sigma_yx_slice_idx := make([]int, n_x * n_y)
        // This holds the indices out of the xy cov slice
        // that would produce sigma_yx_slice.
    for n_xy, k, row := n_x + n_y, 0, 0; row < n_y; row++ {
        for col := 0; col < n_x; col++ {
            sigma_yx_slice_idx[k] = covslice_ij_idx(
                n_xy, dims_y[row], dims_x[col])
        }
    }


    if mix.kind == FreeCovMix {
        for idx_new, idx_old := range idx_keep {
            pick_float64s_fill(mix.CovSlice(idx_old),
                sigma_y_covslice_idx, sigma_y_covslice)
            covslice_expand_fill(sigma_y_covslice, sigma_y_slice)
            sigma_y := mat64.NewDense(n_y, n_y, sigma_y_slice)

            pick_float64s_fill(mix.CovSlice(idx_old),
                sigma_yx_slice_idx, sigma_yx_slice)
            sigma_yx := mat64.NewDense(n_y, n_x, sigma_yx_slice)

            A := mat64.Solve(sigma_y, sigma_yx)
                // FIXME: make use of Cholesky.

            // Get conditional cov slice.
            //
            cov_x := mix_x.CovSlice(idx_new)
            pick_float64s_fill(mix.CovSlice(idx_old), sigma_x_covslice_idx, cov_x)
            var syxt *mat64.Dense
            syxt.TCopy(sigma_yx)
            sigma_x.Mul(syxt, A)
            for k, col := 0, 0; col < n_x; col++ {
                for row := col; row < n_x; row++ {
                    cov_x[k] -= sigma_x.At(row, col)
                    k++
                }
            }

            // Get conditional mean.
            mu_old := mix.MeanSlice(idx_old)
            for i, idx := range dims_y {
                y_delta[i] = data[i] - mu_old[idx]
            }
            mu_x := mix_x.MeanSlice(idx_new)
            for i, idx := range dims_x {
                mu_x[i] = mu_old[idx] + Dot(y_delta, A.RowView(i))
            }
        }

    } else {
        if mix_x.kind == ScaledCovMix {
            copy(mix_x.cov_scale, mix.cov_scale)
        }

        pick_float64s_fill(mix.cov,
            sigma_y_covslice_idx, sigma_y_covslice)
        covslice_expand_fill(sigma_y_covslice, sigma_y_slice)
        sigma_y := mat64.NewDense(n_y, n_y, sigma_y_slice)

        pick_float64s_fill(mix.cov,
            sigma_yx_slice_idx, sigma_yx_slice)
        sigma_yx := mat64.NewDense(n_y, n_x, sigma_yx_slice)

        A := mat64.Solve(sigma_y, sigma_yx)
            // FIXME: make use of Cholesky.

        // Get conditional cov slice.
        //
        pick_float64s_fill(mix.cov, sigma_x_covslice_idx, mix_x.cov)
        var syxt *mat64.Dense
        syxt.TCopy(sigma_yx)
        sigma_x.Mul(syxt, A)
        for k, col := 0, 0; col < n_x; col++ {
            for row := col; row < n_x; row++ {
                mix_x.cov[k] -= sigma_x.At(row, col)
                k++
            }
        }

        // Get conditional mean
        for idx_new, idx_old := range idx_keep {
            mu_old := mix.MeanSlice(idx_old)
            for i, idx := range dims_y {
                y_delta[i] = data[i] - mu_old[idx]
            }
            mu_new := mix_x.MeanSlice(idx_new)
            for i, idx := range dims_x {
                mu_new[i] = mu_old[idx] + Dot(y_delta, A.RowView(i))
            }
        }
    }

    return mix_x
}





// Compute the log-density of the data in each mixture component.
// Input x contains one 'observation' or 'case' per row.
// Return a slice containing the log-densities of all the cases
// in the first mixture component, then the log-densities of all the
// cases in the second mixture component, and so on.
func (mix *Normix) density_stats(x *mat64.Dense, out []float64) []float64 {
    nx, ndim := x.Dims()
    nmix := mix.NMix()

    assert(mix.NDim() == ndim, "input x has wrong shape")
    out = use_slice(out, nx * nmix)

    if ndim == 1 {
        xx := make([]float64, x.Size())
        for i := 0; i < nx; i++ {
            xx[i] = x.At(i, 0)
        }
        if mix.kind == UniCovMix {
            sd := math.Sqrt(mix.cov[0])
            for imix := 0; imix < nmix; imix++ {
                NewNormal(mix.mean[imix], sd).Density(xx, out[nx * imix :])
            }
        } else {
            for imix := 0; imix < nmix; imix++ {
                NewNormal(mix.mean[imix], math.Sqrt(mix.cov[imix])).Density(
                    xx, out[nx * imix :])
            }
        }
        return out
    }


    switch mix.kind {
    case FreeCovMix:
        for imix := 0; imix < nmix; imix++ {
            mvn := NewMvnormal(mix.MeanSlice(imix),
                NewLowerTri(ndim, mix.CovSlice(imix)).Expand(nil))
            mvn.Density(x, out[nx * imix :])
        }
    case UniCovMix:
        xx := mat64.NewDense(nx, ndim, nil)

        mvn := NewMvnormal(make([]float64, ndim),
            NewLowerTri(ndim, mix.cov).Expand(nil))

        for imix := 0; imix < nmix; imix++ {
            xx.Copy(x)
            for row := 0; row < nx; row++ {
                z := xx.RowView(row)
                Subtract(z, mix.MeanSlice(imix), z)
            }

            mvn.Density(xx, out[nx*imix :])
        }

    case ScaledCovMix:
        panic(fmt.Errorf("not implemented"))
        // FIXME: to be done
    default:
        panic(fmt.Errorf("unrecognized kind of Normix"))
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
    floats.AddConst(-floats.LogSumExp(logwt), logwt)

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
func pick_float64s(x []float64, index []int) []float64 {
    y := make([]float64, len(index))
    pick_float64s_fill(x, index, y)
    return y
}




func pick_float64s_fill(x []float64, index []int, y []float64) {
    for i, j := range index {
        y[i] = x[j]
    }
}

