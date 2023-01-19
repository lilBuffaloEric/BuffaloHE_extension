package piecewisefunc

import (
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

/*
Approximate polynomial of the sign function is constructed with the composition of several component polynomials.
Denote the Approximation of sign function as P_a(x), where a represents the precision parameter (the larger it is, the more precise
result we obtain) and is usually between (7~14). Denote the component polynomials as p_a0(x), p_a1(x) ,..., p_ak(x).
Below stores the component polynomials computed by Remez algorithm in different setting of a. Their even terms coeffecients computed
this way are ignored because they are so close to zero.
Reference:
Precise approximation of convolutional neural networks for homomorphically encrypted data. arXiv preprint, abs/2105.10879, 2021. http://arxiv.org/abs/2105.10879.
*/

type SignComposition struct {
	Precision    int
	ComponentNum int
	Components   [3][29]complex128
	TruncateAt   [3]int
}

var (
	SignComponentsA7 = [3][29]complex128{
		{complex(0, 0), complex(7.30445164958251, 0), complex(0, 0), complex(-3.46825871108659e1, 0),
			complex(0, 0), complex(5.98596518298826e1, 0), complex(0, 0), complex(-3.18755225906466e1, 0)},

		{complex(0, 0), complex(2.40085652217597, 0), complex(0, 0), complex(-2.63125454261783, 0),
			complex(0, 0), complex(1.54912674773593, 0), complex(0, 0), complex(-3.31172956504304e-1, 0)},
	}
	TruncateAtA7 = [3]int{8, 8, 0}

	SignComponentsA10 = [3][29]complex128{
		{complex(0, 0), complex(1.08541842577442e1, 0), complex(0, 0), complex(-6.22833925211098e1, 0),
			complex(0, 0), complex(1.14369227820443e2, 0), complex(0, 0), complex(-6.28023496973074e1, 0)},

		{complex(0, 0), complex(4.13976170985111, 0), complex(0, 0), complex(-5.84997640211679, 0),
			complex(0, 0), complex(2.94376255659280, 0), complex(0, 0), complex(-4.54530437460152e-1, 0)},

		{complex(0, 0), complex(3.29956739043733, 0), complex(0, 0), complex(-7.84227260291355, 0),
			complex(0, 0), complex(1.28907764115564e1, 0), complex(0, 0), complex(-1.24917112584486e1, 0),
			complex(0, 0), complex(6.94167991428074, 0), complex(0, 0), complex(-2.04298067399942, 0),
			complex(0, 0), complex(2.46407138926031e-1, 0)},
	}
	TruncateAtA10 = [3]int{8, 8, 14}
)

var (
	Sign7  = SignComposition{Precision: 7, ComponentNum: 2, Components: SignComponentsA7, TruncateAt: TruncateAtA7}
	Sign10 = SignComposition{Precision: 10, ComponentNum: 3, Components: SignComponentsA10, TruncateAt: TruncateAtA10}
)

func Sign_evaluate(params ckks.Parameters, rlk *rlwe.RelinearizationKey, ctIn *rlwe.Ciphertext, signfunc SignComposition) (ctOut *rlwe.Ciphertext, err error) {
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	ctOut = ctIn
	for i := 0; i < signfunc.ComponentNum; i++ {
		Components := signfunc.Components[i][0:signfunc.TruncateAt[i]]
		poly := ckks.NewPoly(Components)
		ctOut, err = evaluator.EvaluatePoly(ctOut, poly, ctOut.Scale)
		if err != nil {
			return nil, err
		}
	}
	return
}

func ReLU_evaluate(params ckks.Parameters, rlk *rlwe.RelinearizationKey, ctIn *rlwe.Ciphertext, signfunc SignComposition) (ctOut *rlwe.Ciphertext, err error) {
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	var ctSign *rlwe.Ciphertext
	ctSign, err = Sign_evaluate(params, rlk, ctIn, signfunc)
	if err != nil {
		return nil, err
	}
	ctInNeg := evaluator.MulRelinNew(ctSign, ctIn)
	evaluator.Rescale(ctInNeg, ctIn.Scale, ctInNeg)
	ctOut = evaluator.AddNew(ctIn, ctInNeg)
	ctOut = evaluator.MultByConstNew(ctOut, 0.5)
	return
}

func Compare_evaluate(params ckks.Parameters, rlk *rlwe.RelinearizationKey, ctA *rlwe.Ciphertext, ctB *rlwe.Ciphertext, signfunc SignComposition) (ctOut *rlwe.Ciphertext, err error) {
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	var ctSign *rlwe.Ciphertext
	ctSign = evaluator.SubNew(ctA, ctB)
	ctSign, err = Sign_evaluate(params, rlk, ctSign, signfunc)
	if err != nil {
		return nil, err
	}
	ctOut = evaluator.AddConstNew(ctSign, 1)
	ctOut = evaluator.MultByConstNew(ctOut, 0.5)
	return
}
