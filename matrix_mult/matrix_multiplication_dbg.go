package matrix_mult

import (
	"errors"
	"fmt"
	"math"
	auxio "project1-fhe_extension/auxiliary_io"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

func SquareMatrix_Mult_dbg(params ckks.Parameters, sk *rlwe.SecretKey, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctMatrixA *rlwe.Ciphertext, ctMatrixB *rlwe.Ciphertext, d int) (ctOut *rlwe.Ciphertext, err error) {
	if float64(d*d) != math.Pow(2, float64(params.LogSlots())) {
		return nil, errors.New("d^2 != 2^logSlots")
	}
	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	LTA := make([]ckks.LinearTransform, d)
	LTB := make([]ckks.LinearTransform, d)
	var LTtest ckks.LinearTransform // testing

	var U_sigma map[int][]float64
	var U_tao map[int][]float64
	var U_colShift_k map[int][]float64
	var U_rowShift_k map[int][]float64
	U_sigma, err = Gen_sigma_diagonalVecotrs(d)
	if err != nil {
		return nil, err
	}
	U_tao, err = Gen_tao_diagonalVectors(d)
	if err != nil {
		return nil, err
	}

	sigmaALT := ckks.GenLinearTransform(encoder, U_sigma, ctMatrixA.Level(), ctMatrixA.Scale, params.LogSlots())
	taoBLT := ckks.GenLinearTransform(encoder, U_tao, ctMatrixB.Level(), ctMatrixB.Scale, params.LogSlots())
	ctSigmaA_list := evaluator.LinearTransformNew(ctMatrixA, sigmaALT)
	ctTaoB_list := evaluator.LinearTransformNew(ctMatrixB, taoBLT)
	var ctSigmaA *rlwe.Ciphertext
	var ctTaoB *rlwe.Ciphertext
	err = evaluator.Rescale(ctSigmaA_list[0], ctMatrixA.Scale, ctSigmaA_list[0])
	if err != nil {
		return nil, err
	}
	ctSigmaA = ctSigmaA_list[0]
	err = evaluator.Rescale(ctTaoB_list[0], ctMatrixB.Scale, ctTaoB_list[0])
	if err != nil {
		return nil, err
	}
	ctTaoB = ctTaoB_list[0]
	auxio.Quick_check_matrix_full(params, sk, ctSigmaA, d, d)
	auxio.Quick_check_matrix_full(params, sk, ctTaoB, d, d)

	for k := 0; k < d; k++ {
		U_colShift_k, err = Gen_colShift_diagonalVectors(d, k)
		if err != nil {
			return nil, err
		}
		U_rowShift_k, err = Gen_rowShift_diagonalVectors(d, k)
		if err != nil {
			return nil, err
		}
		LTA[k] = ckks.GenLinearTransform(encoder, U_colShift_k, ctSigmaA.Level(), ctSigmaA.Scale, params.LogSlots())
		LTB[k] = ckks.GenLinearTransform(encoder, U_rowShift_k, ctTaoB.Level(), ctTaoB.Scale, params.LogSlots())
		if k == 0 { // testing
			LTtest = ckks.GenLinearTransform(encoder, U_rowShift_k, ctTaoB.Level(), ctTaoB.Scale, params.LogSlots())
		}
	}
	ctMatrixALT := evaluator.LinearTransformNew(ctSigmaA, LTA)
	ctMatrixBLT := evaluator.LinearTransformNew(ctTaoB, LTB)
	ctTestB0 := evaluator.LinearTransformNew(ctTaoB, LTB[0])         // testing
	auxio.Quick_check_matrix_full(params, sk, ctTestB0[0], d, d)     // testing
	ctTestSingle := evaluator.LinearTransformNew(ctTaoB, LTtest)     // testing
	auxio.Quick_check_matrix_full(params, sk, ctTestSingle[0], d, d) // testing
	for k := 0; k < d; k++ {
		err = evaluator.Rescale(ctMatrixALT[k], ctSigmaA.Scale, ctMatrixALT[k])
		if err != nil {
			return nil, err
		}
		err = evaluator.Rescale(ctMatrixBLT[k], ctTaoB.Scale, ctMatrixBLT[k])
		if err != nil {
			return nil, err
		}
		auxio.Quick_check_matrix_full(params, sk, ctMatrixALT[k], d, d)
		auxio.Quick_check_matrix_full(params, sk, ctMatrixBLT[k], d, d) // FIXME: Unknow error occurs in the first result of ctMatrixBLT

		if k == 0 {
			ctOut = evaluator.MulRelinNew(ctSigmaA, ctTaoB)
			auxio.Quick_check_matrix_full(params, sk, ctOut, d, d)

		} else {
			ctOut = evaluator.AddNew(ctOut, evaluator.MulRelinNew(ctMatrixALT[k], ctMatrixBLT[k]))
			auxio.Quick_check_matrix_full(params, sk, ctOut, d, d)
		}
	}
	// err = evaluator.Rescale(ctOut, Scale, ctOut)
	if err != nil {
		return nil, err
	}
	fmt.Printf("%d levels consumed for Matrix Multiplication\n", utils.MinInt(ctMatrixA.Level(), ctMatrixB.Level())-ctOut.Level())
	return
}

func Matrix_multirowReplicate_withColOrder_dbg(params ckks.Parameters, sk *rlwe.SecretKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, m0 int, m1 int, f int) (ctOut *rlwe.Ciphertext, err error) {
	if m1 <= 0 || m0 <= 0 {
		return nil, errors.New("invalid input: m1 <= 0 or m0 <= 0 ")
	}
	if m1*m0 != (1 << params.LogSlots()) {
		return nil, errors.New("plaintext matrix's size dose not match the total number of slots")
	}
	T := m0 / f
	if T*f != m0 {
		return nil, errors.New("plaintext matrix's column size cannot be divided by f multirows")
	}
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: galk})
	var ctTemp *rlwe.Ciphertext
	ctOut = ctIn
	Tlen := int(math.Floor(math.Log2(float64(T)) + 1))
	e := 1
	for j := Tlen - 2; j >= 0; j-- {
		ctTemp = evaluator.RotateNew(ctOut, -e*f)
		auxio.Quick_check_matrix_full(params, sk, ctTemp, m0, m1) // testing
		ctOut = evaluator.AddNew(ctTemp, ctOut)
		auxio.Quick_check_matrix_full(params, sk, ctOut, m0, m1) // testing
		e = 2 * e
		if (1<<j)&T == 1 {
			ctTemp = evaluator.RotateNew(ctOut, -f)
			evaluator.Add(ctIn, ctTemp, ctOut)
			e = e + 1
		}
	}
	return
}

func Matrix_rowRotation_withColOrder_dbg(params ckks.Parameters, sk *rlwe.SecretKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, m0 int, m1 int, i int) (ctOut *rlwe.Ciphertext, err error) {
	if m1 <= 0 || m0 <= 0 {
		return nil, errors.New("invalid input: m1 <= 0 or m0 <= 0 ")
	}
	if m1*m0 != (1 << params.LogSlots()) {
		return nil, errors.New("plaintext matrix's size dose not match the total number of slots")
	}
	if i == 0 {
		ctOut = ctIn
		return
	}
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: galk})
	encoder := ckks.NewEncoder(params)
	i %= m0
	if i < 0 {
		i += m0
	}
	RMask := make([][]float64, m0)
	LMask := make([][]float64, m0)
	var LMaskvec []float64
	var RMaskvec []float64
	Scale := ctIn.Scale
	for j := 0; j < m0; j++ {
		RMask[j] = make([]float64, m1)
		LMask[j] = make([]float64, m1)
	}
	for j := m0 - i; j < m0; j++ {
		for k := 0; k < m1; k++ {
			RMask[j][k] = 1
		}
	}
	for j := 0; j < m0-i; j++ {
		for k := 0; k < m1; k++ {
			LMask[j][k] = 1
		}
	}
	RMaskvec, err = Col_orderingInv(RMask)
	if err != nil {
		return nil, err
	}
	LMaskvec, err = Col_orderingInv(LMask)
	if err != nil {
		return nil, err
	}
	ptLMask := encoder.EncodeNew(LMaskvec, params.MaxLevel(), Scale, params.LogSlots())
	ptRMask := encoder.EncodeNew(RMaskvec, params.MaxLevel(), Scale, params.LogSlots())
	ctOut = ctIn
	ctRRotated := evaluator.RotateNew(ctOut, i-m0)
	ctLRotated := evaluator.RotateNew(ctOut, i)
	auxio.Quick_check_matrix_full(params, sk, evaluator.MulNew(ctRRotated, ptRMask), m0, m1)
	auxio.Quick_check_matrix_full(params, sk, evaluator.MulNew(ctLRotated, ptLMask), m0, m1)
	ctOut = evaluator.AddNew(evaluator.MulNew(ctRRotated, ptRMask), evaluator.MulNew(ctLRotated, ptLMask))
	err = evaluator.Rescale(ctOut, Scale, ctOut)
	return
}

func RectangleMatrix_Mult_withColOrder_dbg(params ckks.Parameters, sk *rlwe.SecretKey, galk *rlwe.RotationKeySet, rlk *rlwe.RelinearizationKey, ctMatrixA *rlwe.Ciphertext, ctMatrixB *rlwe.Ciphertext, m0 int, m1 int, m int, l int, n int) (ctOut *rlwe.Ciphertext, err error) {
	var ctMaskedA *rlwe.Ciphertext
	var ctReplicatedA *rlwe.Ciphertext
	var ctRotatedB *rlwe.Ciphertext
	var ctReplicatedB *rlwe.Ciphertext
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: galk, Rlk: rlk})
	if m <= l {
		for i := 0; i < l; i++ {
			ctMaskedA, err = Matrix_masking_withColOrder(params, rlk, ctMatrixA, m0, m1, l, i)
			if err != nil {
				return nil, err
			}
			ctReplicatedA, err = Matrix_rowTotalSum(params, galk, ctMaskedA, m0, m1)
			if err != nil {
				return nil, err
			}
			ctRotatedB, err = Matrix_rowRotation_withColOrder_dbg(params, sk, galk, ctMatrixB, m0, m1, i)
			if err != nil {
				return nil, err
			}
			if i == 0 {
				ctOut = evaluator.MulRelinNew(ctReplicatedA, ctRotatedB)

			} else {
				ctOut = evaluator.AddNew(ctOut, evaluator.MulRelinNew(ctReplicatedA, ctRotatedB))
			}
		}
	} else {
		ctReplicatedB, err = Matrix_multirowReplicate_withColOrder_dbg(params, sk, galk, ctMatrixB, m0, m1, l)
		for i := 0; i < l; i++ {
			ctMaskedA, err = Matrix_masking_withColOrder(params, rlk, ctMatrixA, m0, m1, l, i)
			if err != nil {
				return nil, err
			}
			ctReplicatedA, err = Matrix_rowTotalSum(params, galk, ctMaskedA, m0, m1)
			if err != nil {
				return nil, err
			}
			ctRotatedB, err = Matrix_rowRotation_withColOrder_dbg(params, sk, galk, ctReplicatedB, m0, m1, i)
			if i == 0 {
				ctOut = evaluator.MulRelinNew(ctReplicatedA, ctRotatedB)

			} else {
				ctOut = evaluator.AddNew(ctOut, evaluator.MulRelinNew(ctReplicatedA, ctRotatedB))
			}
		}
	}
	// evaluator.Rescale(ctOut, params.DefaultScale(), ctOut)
	fmt.Printf("%d levels consumed for Matrix Multiplication\n", utils.MinInt(ctMatrixA.Level(), ctMatrixB.Level())-ctOut.Level())
	return
}
