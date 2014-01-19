package stats

// TODO: follow the example of the 'random' standard library
// of C++11 to achieve control and reproducibility of RNG.

import (
	"math/rand"
)

type RandomNumberGenerator interface {
	Seed(int64)
	Next() float64
}

// TODO: temporary solution
type RNG struct{}

func (r *RNG) Seed(seed int64) {
	rand.Seed(seed)
}

func (r *RNG) Next() float64 {
	return rand.Float64()
}
