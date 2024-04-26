package mathematical_func

import "math"

// sigmoid(-x)函数
func SigmoidReversed(x float64) float64 {
	x = -x
	if x > 0 {
		return 1.0 / (1.0 + math.Exp(-x))
	} else {
		return math.Exp(x) / (1.0 + math.Exp(x))
	}
}

// sigmoid(x)函数
func Sigmoid(x float64) float64 {
	if x > 0 {
		return 1.0 / (1.0 + math.Exp(-x))
	} else {
		return math.Exp(x) / (1.0 + math.Exp(x))
	}
}
