package matrix_mult

import (
	"errors"
	"fmt"
	"math"
	"time"

	auxio "project1-fhe_extension_v1.0/auxiliary_io"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

func Sigma_linearTransform(params ckks.Parameters, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, d int) (ctOut *rlwe.Ciphertext, err error) {
	if d <= 0 || float64(d*d) != math.Pow(2, float64(params.LogSlots())) {
		return nil, errors.New("d<0 or d^2 != 2^logSlots")
	}
	var U_sigma map[int][]float64
	U_sigma, err = Gen_sigma_diagonalVecotrs(d)
	if err != nil {
		return nil, err
	}
	Scale := ctIn.Scale
	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	sigmaLT := ckks.GenLinearTransform(encoder, U_sigma, ctIn.Level(), Scale, params.LogSlots())
	ctOut_list := evaluator.LinearTransformNew(ctIn, sigmaLT)
	ctOut = ctOut_list[0]
	err = evaluator.Rescale(ctOut, ctIn.Scale, ctOut)
	if err != nil {
		return nil, err
	}
	fmt.Printf("%d levels consumed for Sigma LinearTransform\n", ctIn.Level()-ctOut.Level())
	return
}

func Sigma_linearTransformBSGS(params ckks.Parameters, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, d int, BSGSRatio float64) (ctOut *rlwe.Ciphertext, err error) {
	if d <= 0 || float64(d*d) != math.Pow(2, float64(params.LogSlots())) {
		return nil, errors.New("d<0 or d^2 != 2^logSlots")
	}
	var U_sigma map[int][]float64
	U_sigma, err = Gen_sigma_diagonalVecotrs(d)
	if err != nil {
		return nil, err
	}
	Scale := ctIn.Scale
	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	sigmaLT := ckks.GenLinearTransformBSGS(encoder, U_sigma, ctIn.Level(), Scale, BSGSRatio, params.LogSlots())
	ctOut_list := evaluator.LinearTransformNew(ctIn, sigmaLT)
	ctOut = ctOut_list[0]
	err = evaluator.Rescale(ctOut, ctIn.Scale, ctOut)
	if err != nil {
		return nil, err
	}
	fmt.Printf("%d levels consumed for Sigma LinearTransformBSGS\n", ctIn.Level()-ctOut.Level())
	return
}

func Tao_linearTransform(params ckks.Parameters, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, d int) (ctOut *rlwe.Ciphertext, err error) {
	if d <= 0 || float64(d*d) != math.Pow(2, float64(params.LogSlots())) {
		return nil, errors.New("d<0 or d^2 != 2^logSlots")
	}
	var U_tao map[int][]float64
	U_tao, err = Gen_tao_diagonalVectors(d)
	if err != nil {
		return nil, err
	}
	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	taoLT := ckks.GenLinearTransform(encoder, U_tao, ctIn.Level(), ctIn.Scale, params.LogSlots())
	ctOut_list := evaluator.LinearTransformNew(ctIn, taoLT)
	ctOut = ctOut_list[0]
	err = evaluator.Rescale(ctOut, ctIn.Scale, ctOut)
	if err != nil {
		return nil, err
	}
	fmt.Printf("%d levels consumed for Tao LinearTransform\n", ctIn.Level()-ctOut.Level())
	return
}

func Tao_linearTransformBSGS(params ckks.Parameters, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, d int, BSGSRatio float64) (ctOut *rlwe.Ciphertext, err error) {
	if d <= 0 || float64(d*d) != math.Pow(2, float64(params.LogSlots())) {
		return nil, errors.New("d<0 or d^2 != 2^logSlots")
	}
	var U_tao map[int][]float64
	U_tao, err = Gen_tao_diagonalVectors(d)
	if err != nil {
		return nil, err
	}
	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	taoLT := ckks.GenLinearTransformBSGS(encoder, U_tao, ctIn.Level(), ctIn.Scale, BSGSRatio, params.LogSlots())
	ctOut_list := evaluator.LinearTransformNew(ctIn, taoLT)
	ctOut = ctOut_list[0]
	err = evaluator.Rescale(ctOut, ctIn.Scale, ctOut)
	if err != nil {
		return nil, err
	}
	fmt.Printf("%d levels consumed for Tao LinearTransformBSGS\n", ctIn.Level()-ctOut.Level())
	return
}

func ColShift_linearTransform(params ckks.Parameters, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, d int, k int) (ctOut *rlwe.Ciphertext, err error) {
	if d <= 0 || float64(d*d) != math.Pow(2, float64(params.LogSlots())) {
		return nil, errors.New("d<0 or d^2 != 2^logSlots")
	}
	if k <= -d || k >= d {
		return nil, errors.New("k <= -d, k >= d")
	}
	if k < 0 {
		k += d
	}
	var U_colShift map[int][]float64
	U_colShift, err = Gen_colShift_diagonalVectors(d, k)
	if err != nil {
		return nil, err
	}

	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	colShiftLT := ckks.GenLinearTransform(encoder, U_colShift, ctIn.Level(), ctIn.Scale, params.LogSlots())
	ctOut_list := evaluator.LinearTransformNew(ctIn, colShiftLT)
	ctOut = ctOut_list[0]
	err = evaluator.Rescale(ctOut, ctIn.Scale, ctOut)
	if err != nil {
		return nil, err
	}
	fmt.Printf("%d levels consumed for ColShift LinearTransform\n", ctIn.Level()-ctOut.Level())
	return
}

func RowShift_linearTransform(params ckks.Parameters, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, d int, k int) (ctOut *rlwe.Ciphertext, err error) {
	if d <= 0 || float64(d*d) != math.Pow(2, float64(params.LogSlots())) {
		return nil, errors.New("d<0 or d^2 != 2^logSlots")
	}
	if k <= -d || k >= d {
		return nil, errors.New("k <= -d, k >= d")
	}
	if k < 0 {
		k += d
	}
	var U_rowShift map[int][]float64
	U_rowShift, err = Gen_rowShift_diagonalVectors(d, k)
	if err != nil {
		return nil, err
	}

	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	rowShiftLT := ckks.GenLinearTransform(encoder, U_rowShift, ctIn.Level(), ctIn.Scale, params.LogSlots())
	ctOut_list := evaluator.LinearTransformNew(ctIn, rowShiftLT)
	ctOut = ctOut_list[0]
	err = evaluator.Rescale(ctOut, ctIn.Scale, ctOut)
	if err != nil {
		return nil, err
	}
	fmt.Printf("%d levels consumed for RowShift LinearTransform\n", ctIn.Level()-ctOut.Level())
	return
}

func Transpose_linearTransform(params ckks.Parameters, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, d int) (ctOut *rlwe.Ciphertext, err error) {
	if d <= 0 || float64(d*d) != math.Pow(2, float64(params.LogSlots())) {
		return nil, errors.New("d<0 or d^2 != 2^logSlots")
	}

	var U_transpose map[int][]float64
	U_transpose, err = Gen_transpose_diagonalVectors(d)
	if err != nil {
		return nil, err
	}

	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	transposeLT := ckks.GenLinearTransform(encoder, U_transpose, ctIn.Level(), ctIn.Scale, params.LogSlots())
	ctOut_list := evaluator.LinearTransformNew(ctIn, transposeLT)
	ctOut = ctOut_list[0]
	err = evaluator.Rescale(ctOut, ctIn.Scale, ctOut)
	if err != nil {
		return nil, err
	}
	fmt.Printf("%d levels consumed for Transpose LinearTransform\n", ctIn.Level()-ctOut.Level())
	return
}

func TransposeCtao_linearTransform(params ckks.Parameters, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, d int) (ctOut *rlwe.Ciphertext, err error) {
	if d <= 0 || float64(d*d) != math.Pow(2, float64(params.LogSlots())) {
		return nil, errors.New("d<0 or d^2 != 2^logSlots")
	}
	var U_transCtao map[int][]float64
	U_transCtao, err = Gen_trans_C_tao_diagonalVectors(d)
	if err != nil {
		return nil, err
	}
	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	transCtaoLT := ckks.GenLinearTransform(encoder, U_transCtao, ctIn.Level(), ctIn.Scale, params.LogSlots())
	ctOut_list := evaluator.LinearTransformNew(ctIn, transCtaoLT)
	ctOut = ctOut_list[0]
	err = evaluator.Rescale(ctOut, ctIn.Scale, ctOut)
	if err != nil {
		return nil, err
	}
	fmt.Printf("%d levels consumed for LinearTransform\n", ctIn.Level()-ctOut.Level())
	return
}

func SquareMatrix_Mult(params ckks.Parameters, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctMatrixA *rlwe.Ciphertext, ctMatrixB *rlwe.Ciphertext, d int) (ctOut *rlwe.Ciphertext, err error) {
	if float64(d*d) != math.Pow(2, float64(params.LogSlots())) {
		return nil, errors.New("d^2 != 2^logSlots")
	}
	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	LTA := make([]ckks.LinearTransform, d)
	LTB := make([]ckks.LinearTransform, d)
	// var LTtest ckks.LinearTransform // dbg testing
	now := time.Now() // dbg testing
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
	// now := time.Now() // dbg testing
	ctSigmaA_list := evaluator.LinearTransformNew(ctMatrixA, sigmaALT)
	ctTaoB_list := evaluator.LinearTransformNew(ctMatrixB, taoBLT)
	// fmt.Printf("normal LT consume %s\n", time.Since(now)) // dbg testing
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
	// auxio.Quick_check_matrix_full(params, sk, ctSigmaA, d, d)
	// auxio.Quick_check_matrix_full(params, sk, ctTaoB, d, d)

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
		/*
			if k == 0 { // testing
				LTtest = ckks.GenLinearTransform(encoder, U_rowShift_k, ctTaoB.Level(), ctTaoB.Scale, params.LogSlots())
			}
		*/
	}
	ctMatrixALT := evaluator.LinearTransformNew(ctSigmaA, LTA)
	ctMatrixBLT := evaluator.LinearTransformNew(ctTaoB, LTB)
	// ctTestB0 := evaluator.LinearTransformNew(ctTaoB, LTB[0])         // dbg testing
	// auxio.Quick_check_matrix_full(params, sk, ctTestB0[0], d, d)     // dbg testing
	// ctTestSingle := evaluator.LinearTransformNew(ctTaoB, LTtest)     // dbg testing
	// auxio.Quick_check_matrix_full(params, sk, ctTestSingle[0], d, d) // dbg testing
	for k := 0; k < d; k++ {
		err = evaluator.Rescale(ctMatrixALT[k], ctSigmaA.Scale, ctMatrixALT[k])
		if err != nil {
			return nil, err
		}
		err = evaluator.Rescale(ctMatrixBLT[k], ctTaoB.Scale, ctMatrixBLT[k])
		if err != nil {
			return nil, err
		}
		// auxio.Quick_check_matrix_full(params, sk, ctMatrixALT[k], d, d) // dbg testing
		// auxio.Quick_check_matrix_full(params, sk, ctMatrixBLT[k], d, d) // FIXME: Unknow error occurs in the first result of ctMatrixBLT
		if k == 0 {
			var ctTemp *rlwe.Ciphertext                                             // dbg testing
			ctTemp, err = RowShift_linearTransform(params, rlk, galk, ctTaoB, d, 0) // dbg testing
			ctOut = evaluator.MulNew(ctMatrixALT[k], ctTemp)                        // dbg testing
			// auxio.Quick_check_matrix_full(params, sk, ctOut, d, d) // dbg testing

		} else {
			evaluator.MulAndAdd(ctMatrixALT[k], ctMatrixBLT[k], ctOut) // dbg testing
			// auxio.Quick_check_matrix_full(params, sk, ctOut, d, d) // dbg testing
		}
	}
	// err = evaluator.Rescale(ctOut, Scale, ctOut)
	evaluator.Relinearize(ctOut, ctOut)
	if err != nil {
		return nil, err
	}

	fmt.Printf("%d levels, %s consumed for Matrix Multiplication\n", utils.MinInt(ctMatrixA.Level(), ctMatrixB.Level())-ctOut.Level(), time.Since(now))
	return
}

// diagonal masking only
// this function has used params.MaxLevel() setting the plaintext's level for internal use, this should be inspect carefully while using bootstrap scheme.
func Matrix_masking_withRowOrder(params ckks.Parameters, rlk *rlwe.RelinearizationKey, ctIn *rlwe.Ciphertext, m0 int, m1 int, l int, i int) (ctOut *rlwe.Ciphertext, err error) {
	if l <= 0 || m0 <= 0 || m1 <= 0 {
		return nil, errors.New("l<=0 || m0 <= 0 || m1 <= 0")
	}
	if m0*m1 > (1 << params.LogSlots()) {
		return nil, errors.New("hypercube volumn larger than number slots available")
	}
	if m0 < l {
		return nil, errors.New("m0 < l")
	}
	if i < 0 || i >= l {
		return nil, errors.New("i<0 || i>= l")
	}

	Ui := make([][]float64, m0)
	for I := 0; I < m0; I++ {
		Ui[I] = make([]float64, m1)
		for J := 0; J < m1; J++ {
			if (J-(i+I))%l == 0 {
				Ui[I][J] = 1
			} else {
				Ui[I][J] = 0
			}
		}
	}
	var ui []float64
	ui, err = Row_orderingInv(Ui)
	if err != nil {
		return nil, err
	}
	Scale := ctIn.Scale
	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	ptUi := encoder.EncodeNew(ui, params.MaxLevel(), Scale, params.LogSlots())
	// auxio.Quick_decode_vector(params, ptUi)
	ctOut = evaluator.MulNew(ctIn, ptUi)
	err = evaluator.Rescale(ctOut, Scale, ctOut)
	fmt.Printf("%d levels consumed for Matrix masking withRowOrder\n", ctIn.Level()-ctOut.Level())
	return
}

// diagonal masking only
// this function has used params.MaxLevel() setting the plaintext's level for internal use, this should be inspect carefully while using bootstrap scheme.
func Matrix_masking_withColOrder(params ckks.Parameters, rlk *rlwe.RelinearizationKey, ctIn *rlwe.Ciphertext, m0 int, m1 int, l int, i int) (ctOut *rlwe.Ciphertext, err error) {
	if l <= 0 || m0 <= 0 || m1 <= 0 {
		return nil, errors.New("l<=0 || m0 <= 0 || m1 <= 0")
	}
	if m0*m1 > (1 << params.LogSlots()) {
		return nil, errors.New("matrix's volumn larger than number slots available")
	}
	if m0 < l {
		return nil, errors.New("m0 < l")
	}
	if i < 0 || i >= l {
		return nil, errors.New("i<0 || i>= l")
	}

	Ui := make([][]float64, m0)
	for I := 0; I < m0; I++ {
		Ui[I] = make([]float64, m1)
		for J := 0; J < m1; J++ {
			if (J-(i+I))%l == 0 {
				Ui[I][J] = 1
			} else {
				Ui[I][J] = 0
			}
		}
	}
	var ui []float64
	ui, err = Col_orderingInv(Ui)
	if err != nil {
		return nil, err
	}
	Scale := ctIn.Scale
	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	ptUi := encoder.EncodeNew(ui, params.MaxLevel(), Scale, params.LogSlots())
	auxio.Quick_decode_vector(params, ptUi)
	ctOut = evaluator.MulNew(ctIn, ptUi)
	err = evaluator.Rescale(ctOut, Scale, ctOut)
	fmt.Printf("%d levels consumed for Matrix masking withRowOrder\n", ctIn.Level()-ctOut.Level())
	return
}

// input Ciphertext must be some encryption of a plaintext space that can be represented by a Matrix with size m0 x m1,
// which means m0 x m1 = 2^{logslots}.
// Moreover, If the plaintext matrix is column-ordering encoded, then for each row of the plain matrix,
// all the elements (in the row) will become the sum of themselves.
// If the plaintext matrix is row-ordering encoded, then for each column of the plain matrix,
// all the elements (in the column) will become the sum of themselves.
func Matrix_rowTotalSum_withColOrder(params ckks.Parameters, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, m0 int, m1 int) (ctOut *rlwe.Ciphertext, err error) {
	if m1 <= 0 || m0 <= 0 {
		return nil, errors.New("invalid input: m1 <= 0 or m0 <= 0 ")
	}
	if m1*m0 != (1 << params.LogSlots()) {
		return nil, errors.New("plaintext matrix's size dose not match the total number of slots")
	}
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: galk})
	ctOut = ctIn
	for i := m0; i < m1*m0; i = (i << 1) {
		ctTemp := evaluator.RotateNew(ctOut, i)
		evaluator.Add(ctOut, ctTemp, ctOut)
	}
	return
}

// input Ciphertext must be some encryption of a plaintext space that can be represented by a Matrix with size m0 x m1,
// which means m0 x m1 = 2^{logslots}. Moreover, the plaintext matrix should be column ordering encoded.
func Matrix_multirowReplicate_withColOrder(params ckks.Parameters, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, m0 int, m1 int, f int) (ctOut *rlwe.Ciphertext, err error) {
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
		// auxio.Quick_check_matrix_full(params, sk, ctTemp, m0, m1) // dbg testing
		ctOut = evaluator.AddNew(ctTemp, ctOut)
		// auxio.Quick_check_matrix_full(params, sk, ctOut, m0, m1) // dbg testing
		e = 2 * e
		if (1<<j)&T == 1 {
			ctTemp = evaluator.RotateNew(ctOut, -f)
			evaluator.Add(ctIn, ctTemp, ctOut)
			e = e + 1
		}
	}
	return
}

// this function has used params.MaxLevel() setting the plaintext's level for internal use, this should be inspect carefully while using bootstrap scheme.
func Matrix_rowRotation_withColOrder(params ckks.Parameters, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, m0 int, m1 int, i int) (ctOut *rlwe.Ciphertext, err error) {
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
	// auxio.Quick_check_matrix_full(params, sk, evaluator.MulNew(ctRRotated, ptRMask), m0, m1) // dbg testing
	// auxio.Quick_check_matrix_full(params, sk, evaluator.MulNew(ctLRotated, ptLMask), m0, m1) // dbg testing
	ctOut = evaluator.AddNew(evaluator.MulNew(ctRRotated, ptRMask), evaluator.MulNew(ctLRotated, ptLMask))
	err = evaluator.Rescale(ctOut, Scale, ctOut)
	return
}

func RectangleMatrix_Mult_withColOrder(params ckks.Parameters, galk *rlwe.RotationKeySet, rlk *rlwe.RelinearizationKey, ctMatrixA *rlwe.Ciphertext, ctMatrixB *rlwe.Ciphertext, m0 int, m1 int, m int, l int, n int) (ctOut *rlwe.Ciphertext, err error) {
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
			ctReplicatedA, err = Matrix_rowTotalSum_withColOrder(params, galk, ctMaskedA, m0, m1)
			if err != nil {
				return nil, err
			}
			ctRotatedB, err = Matrix_rowRotation_withColOrder(params, galk, ctMatrixB, m0, m1, i)
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
		ctReplicatedB, err = Matrix_multirowReplicate_withColOrder(params, galk, ctMatrixB, m0, m1, l)
		for i := 0; i < l; i++ {
			ctMaskedA, err = Matrix_masking_withColOrder(params, rlk, ctMatrixA, m0, m1, l, i)
			if err != nil {
				return nil, err
			}
			ctReplicatedA, err = Matrix_rowTotalSum_withColOrder(params, galk, ctMaskedA, m0, m1)
			if err != nil {
				return nil, err
			}
			ctRotatedB, err = Matrix_rowRotation_withColOrder(params, galk, ctReplicatedB, m0, m1, i)
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

// this function has used params.MaxLevel() setting the plaintext's level for internal use, this should be inspect carefully while using bootstrap scheme.
func Matrix_Aggregate_withRowOrder(params ckks.Parameters, sk *rlwe.SecretKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, row int, col int, axis int) (ctOut *rlwe.Ciphertext, err error) {
	if row <= 0 || col <= 0 {
		return nil, errors.New("invalid input: m1 <= 0 or m0 <= 0 ")
	}
	if row*col != (1 << params.LogSlots()) {
		return nil, errors.New("matrix's size dose not match the total number of slots")
	}
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: galk})
	ctOut = ctIn.CopyNew()
	if axis == 0 {
		// Aggregate in Row order. The ctOut will have replicated rows repeated col times.
		for i := col; i < row*col; i = (i << 1) {
			// fmt.Printf("before Rotation:\n")
			// auxio.Quick_check_matrix_full(params, sk, ctOut, row, col)
			ctTemp := evaluator.RotateNew(ctOut, i)
			evaluator.Add(ctOut, ctTemp, ctOut)
			// fmt.Printf("After Rotation and add:\n")
			// auxio.Quick_check_matrix_full(params, sk, ctOut, row, col)
		}
	}
	if axis == 1 {
		encoder := ckks.NewEncoder(params)
		MaskVec := make([]float64, 1<<params.LogSlots())
		for i := 0; i < len(MaskVec); i += col {
			MaskVec[i] = 1
		}
		ptMask := encoder.EncodeNew(MaskVec, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		// Aggregate in Col order. The ctOut will have replicated cols repeated row times.
		for i := 1; i < col; i = (i << 1) {
			ctTemp := evaluator.RotateNew(ctOut, i)
			evaluator.Add(ctOut, ctTemp, ctOut)
		}
		evaluator.Mul(ctOut, ptMask, ctOut)
		err = evaluator.Rescale(ctOut, params.DefaultScale(), ctOut)
		if err != nil {
			panic(err)
		}
		for i := -1; i > -col; i = (i << 1) {
			ctTemp := evaluator.RotateNew(ctOut, i)
			evaluator.Add(ctOut, ctTemp, ctOut)
		}
	}
	return
}

// axis is the current axis of the Vector.
func ReplicateVec_Switch_Axis(params ckks.Parameters, sk *rlwe.SecretKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, d int, axis int) (ctOut *rlwe.Ciphertext, err error) {
	if d <= 0 {
		return nil, errors.New("invalid input:d <= 0 ")
	}
	if d*d != (1 << params.LogSlots()) {
		return nil, errors.New("the 2-power of vector dimension is not equal to the total number of slots")
	}
	targetAxis := (axis + 1) & (2 - 1)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: galk})
	encoder := ckks.NewEncoder(params)
	Mask := make([]float64, d*d)
	for i := 0; i < d; i++ {
		Mask[i*d+i] = 1
	}
	ptMask := encoder.EncodeNew(Mask, ctIn.Level(), params.DefaultScale(), params.LogSlots())
	ctTemp := evaluator.MulNew(ctIn, ptMask)
	err = evaluator.Rescale(ctTemp, params.DefaultScale(), ctTemp)
	if err != nil {
		panic(err)
	}
	ctOut, err = Matrix_Aggregate_withRowOrder(params, sk, galk, ctTemp, d, d, targetAxis)
	return

}

func ReplicatedVec_Combine_dbg(params ckks.Parameters, sk *rlwe.SecretKey, ctIn []*rlwe.Ciphertext, d int, axis int) (ctOut *rlwe.Ciphertext, err error) {
	if d <= 0 {
		return nil, errors.New("invalid input:d <= 0 ")
	}
	if d*d != (1 << params.LogSlots()) {
		return nil, errors.New("the 2-power of vector dimension is not equal to the total number of slots")
	}
	if len(ctIn) > d {
		return nil, errors.New("cannot combine more than d vectors")
	}
	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{})
	Mask := make([]float64, d*d)
	Scale := rlwe.NewScale(float64(params.RingQ().Modulus[ctIn[0].Level()])).Mul(params.DefaultScale()).Div(ctIn[0].Scale)
	if axis == 1 {
		for i := 0; i < len(ctIn); i++ {
			for j := i; j < i+(d*d); j += d {
				Mask[j] = 1
			}
			Maskpt := encoder.EncodeNew(Mask, ctIn[i].Level(), Scale, params.LogSlots())
			if i == 0 {
				ctOut = evaluator.MulNew(ctIn[i], Maskpt)
			} else {
				evaluator.MulAndAdd(ctIn[i], Maskpt, ctOut)
			}
			for j := i; j < i+(d*d); j += d {
				Mask[j] = 0
			}
		}
	} else if axis == 0 {
		for i := 0; i < len(ctIn); i++ {
			for j := i * d; j < (i+1)*d; j++ {
				Mask[j] = 1
			}
			Maskpt := encoder.EncodeNew(Mask, ctIn[i].Level(), Scale, params.LogSlots())
			if i == 0 {
				ctOut = evaluator.MulNew(ctIn[i], Maskpt)
			} else {
				evaluator.MulAndAdd(ctIn[i], Maskpt, ctOut)
			}
			for j := i * d; j < (i+1)*d; j++ {
				Mask[j] = 0
			}
		}
	}
	return
}

func ReplicateVec_Decompose_dbg(params ckks.Parameters, sk *rlwe.SecretKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, d int, axis int, vecNum int) (ctOut []*rlwe.Ciphertext, err error) {
	ctOut = make([]*rlwe.Ciphertext, vecNum)
	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: galk})
	Mask := make([]float64, d*d)
	if axis == 1 {
		for i := 0; i < len(ctOut); i++ {
			for j := i; j < i+(d*d); j += d {
				Mask[j] = 1
			}
			Maskpt := encoder.EncodeNew(Mask, ctIn.Level(), params.DefaultScale(), params.LogSlots())
			ctOut[i] = evaluator.MulNew(ctIn, Maskpt)
			err = evaluator.Rescale(ctOut[i], params.DefaultScale(), ctOut[i])
			if err != nil {
				return nil, err
			}
			ctOut[i], err = Matrix_Aggregate_withRowOrder(params, sk, galk, ctOut[i], d, d, axis)
			if err != nil {
				return nil, err
			}
			for j := i; j < i+(d*d); j += d {
				Mask[j] = 0
			}
		}
	} else if axis == 0 {
		for i := 0; i < len(ctOut); i++ {
			for j := i * d; j < (i+1)*d; j++ {
				Mask[j] = 1
			}
			Maskpt := encoder.EncodeNew(Mask, ctIn.Level(), params.DefaultScale(), params.LogSlots())
			ctOut[i] = evaluator.MulNew(ctIn, Maskpt)
			err = evaluator.Rescale(ctOut[i], params.DefaultScale(), ctOut[i])
			if err != nil {
				return nil, err
			}
			ctOut[i], err = Matrix_Aggregate_withRowOrder(params, sk, galk, ctOut[i], d, d, axis)
			if err != nil {
				return nil, err
			}
			for j := i * d; j < (i+1)*d; j++ {
				Mask[j] = 0
			}
		}
	}
	return
}
