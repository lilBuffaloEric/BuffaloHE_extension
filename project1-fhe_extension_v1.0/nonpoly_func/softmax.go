package nonpolyfunc

import (
	"errors"
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

var (
	Exp_LeastSquare12 = [13]complex128{
		complex(0.082783231973869, 0), complex(-1.075804435198659e-04, 0), complex(6.291247875189161e-08, 0),
	}
)

// this function is only for experiment purpose, since it will internally decrypt a ciphertext
func InvByGoldSchmidt(params ckks.Parameters, sk *rlwe.SecretKey, pk *rlwe.PublicKey, rlk *rlwe.RelinearizationKey, ctIn *rlwe.Ciphertext, interval float64, iterNum int) (ctOut *rlwe.Ciphertext, err error) {
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	decryptor := ckks.NewDecryptor(params, sk)
	encryptor := ckks.NewEncryptor(params, pk)
	encoder := ckks.NewEncoder(params)
	// first shrink the input
	if ctIn.Level() <= 2 { // prevent dropping to one (some bugs may be triggered)
		ctOut = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctIn), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
	} else {
		ctOut = ctIn.CopyNew()
	}
	ctOut = evaluator.MultByConstNew(ctOut, -2/interval)
	err = evaluator.Rescale(ctOut, ctIn.Scale, ctOut)
	if err != nil {
		return nil, err
	}
	// initialization
	cta0 := evaluator.AddConstNew(ctOut, 2)
	ctb0 := evaluator.AddConstNew(ctOut, 1)
	// Into the iteration.
	// Two levels consumed for each iteration.
	for i := 0; i < iterNum; i++ {
		if ctb0.Level() < 1 { // Do one Recryption.
			ctb0 = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctb0), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
		}
		if cta0.Level() < 2 {
			cta0 = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(cta0), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
		}
		ctb0 = evaluator.MulRelinNew(ctb0, ctb0)
		err = evaluator.Rescale(ctb0, ctIn.Scale, ctb0)
		if err != nil {
			return nil, err
		}
		cta0 = evaluator.MulRelinNew(cta0, evaluator.AddConstNew(ctb0, 1))
		err = evaluator.Rescale(cta0, ctIn.Scale, cta0)
		if err != nil {
			return nil, err
		}
	}
	ctOut = cta0
	if ctOut.Level() < 1 {
		ctOut = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctOut), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
	}
	ctOut = evaluator.MultByConstNew(ctOut, 2/interval)
	err = evaluator.Rescale(ctOut, ctIn.Scale, ctOut)
	fmt.Printf("Levels consumed by InvByGoldSchmidt: %d\n", ctIn.Level()-ctOut.Level())
	return
}

func ExpByLeastSquare(params ckks.Parameters, rlk *rlwe.RelinearizationKey, ctIn *rlwe.Ciphertext, interval int, apprxpoly []complex128) (ctOut *rlwe.Ciphertext, err error) {
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	// first shrink the input by 1/interval.
	if interval < 0 {
		return nil, errors.New("invalid input : interval <0")
	}
	if interval > 1 {
		ctOut = evaluator.MultByConstNew(ctIn, 1/interval)
		err = evaluator.Rescale(ctOut, ctIn.Scale, ctOut)
		if err != nil {
			return nil, err
		}
	}
	poly := ckks.NewPoly(apprxpoly)
	ctOut, err = evaluator.EvaluatePoly(ctOut, poly, ctIn.Scale)
	// Then exponentially amplify by interval
	if interval > 1 {
		Ilen := int(math.Floor(math.Log2(float64(interval)) + 1))
		ctPow := ctOut
		for j := Ilen - 2; j >= 0; j-- {
			ctPow = evaluator.MulRelinNew(ctPow, ctPow)
			err = evaluator.Rescale(ctPow, ctIn.Scale, ctPow)
			if err != nil {
				return nil, err
			}
			if (1<<j)&interval == 1 {
				ctOut = evaluator.MulRelinNew(ctOut, ctPow)
				err = evaluator.Rescale(ctOut, ctIn.Scale, ctOut)
				if err != nil {
					return nil, err
				}
			} else {
				ctOut = ctPow
			}
		}
	}
	return
}
