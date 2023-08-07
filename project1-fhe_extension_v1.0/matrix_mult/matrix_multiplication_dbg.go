package matrix_mult

import (
	"errors"
	"fmt"
	"math"
	"sort"
	"sync"
	"time"

	auxio "project1-fhe_extension_v1.0/auxiliary_io"

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

func SquareMatrix_MultBSGS_dbg(params ckks.Parameters, sk *rlwe.SecretKey, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctMatrixA *rlwe.Ciphertext, ctMatrixB *rlwe.Ciphertext, d int, sigmaBSGSRatio float64, taoBSGSRatio float64) (ctOut *rlwe.Ciphertext, err error) {
	if float64(d*d) != math.Pow(2, float64(params.LogSlots())) {
		return nil, errors.New("d^2 != 2^logSlots")
	}
	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	LTA := make([]ckks.LinearTransform, d)
	LTB := make([]ckks.LinearTransform, d) // option 1
	//StepsB := make([]int, d) // option 2
	var U_sigma map[int][]float64
	var U_tao map[int][]float64
	var U_colShift_k map[int][]float64
	var U_rowShift_k map[int][]float64 // option 1
	U_sigma, err = Gen_sigma_diagonalVecotrs(d)
	if err != nil {
		return nil, err
	}
	U_tao, err = Gen_tao_diagonalVectors(d)
	if err != nil {
		return nil, err
	}

	sigmaALT := ckks.GenLinearTransformBSGS(encoder, U_sigma, ctMatrixA.Level(), ctMatrixA.Scale, sigmaBSGSRatio, params.LogSlots())
	taoBLT := ckks.GenLinearTransformBSGS(encoder, U_tao, ctMatrixB.Level(), ctMatrixB.Scale, taoBSGSRatio, params.LogSlots())
	now := time.Now() // dbg testing
	ctSigmaA_list := evaluator.LinearTransformNew(ctMatrixA, sigmaALT)
	ctTaoB_list := evaluator.LinearTransformNew(ctMatrixB, taoBLT)
	fmt.Printf("BSGS LT consume %s\n", time.Since(now)) // dbg testing
	var ctSigmaA *rlwe.Ciphertext
	var ctTaoB *rlwe.Ciphertext
	var ctMatrixALT []*rlwe.Ciphertext
	//var ctMatrixBLT map[int]*rlwe.Ciphertext // option 2
	var ctMatrixBLT []*rlwe.Ciphertext // option 1
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
	for k := 0; k < d; k++ {
		U_colShift_k, err = Gen_colShift_diagonalVectors(d, k)
		if err != nil {
			return nil, err
		}
		LTA[k] = ckks.GenLinearTransform(encoder, U_colShift_k, ctSigmaA.Level(), ctSigmaA.Scale, params.LogSlots())
		//StepsB[k] = k * d //option 2
		U_rowShift_k, err = Gen_rowShift_diagonalVectors(d, k)                                                   // option 1
		LTB[k] = ckks.GenLinearTransform(encoder, U_rowShift_k, ctTaoB.Level(), ctTaoB.Scale, params.LogSlots()) // option 1
	}

	// ColShiftRotationsNeeded :=
	// fmt.Printf("Rotation needed for ColumnShift linearTransformation:") // dbg testing.
	// for i:=0;

	ctMatrixALT = evaluator.LinearTransformNew(ctSigmaA, LTA)
	//ctMatrixBLT = evaluator.RotateHoistedNew(ctTaoB, StepsB) // option 2
	ctMatrixBLT = evaluator.LinearTransformNew(ctTaoB, LTB) // option 1

	for k := 0; k < d; k++ {
		err = evaluator.Rescale(ctMatrixALT[k], ctSigmaA.Scale, ctMatrixALT[k])
		if err != nil {
			return nil, err
		}
		err = evaluator.Rescale(ctMatrixBLT[k], ctTaoB.Scale, ctMatrixBLT[k]) // option 1
		if k == 0 {
			//ctOut = evaluator.MulRelinNew(ctMatrixALT[k], ctMatrixBLT[k*d]) // option 2
			var ctTemp *rlwe.Ciphertext                                             // dbg testing
			ctTemp, err = RowShift_linearTransform(params, rlk, galk, ctTaoB, d, 0) // dbg testing
			// auxio.Quick_check_matrix(params, sk, ctTemp, d, d)                      // dbg testing
			ctOut = evaluator.MulNew(ctMatrixALT[k], ctTemp) // option 1 dbg testing

		} else {
			//ctOut = evaluator.AddNew(ctOut, evaluator.MulRelinNew(ctMatrixALT[k], ctMatrixBLT[k*d])) // option 2
			//ctTemp := evaluator.MulRelinNew(ctMatrixALT[k], ctMatrixBLT[k]) // dbg testing
			//fmt.Printf("middle mult of iteration no.%d checking:\n", k)     // dbg testing
			//auxio.Quick_check_matrix(params, sk, ctTemp, d, d)              // dbg testing
			evaluator.MulAndAdd(ctMatrixALT[k], ctMatrixBLT[k], ctOut) // option 1
			// ctOut = evaluator.AddNew(ctOut, evaluator.MulNew(ctMatrixALT[k], ctMatrixBLT[k])) // option 1 (redundant)
		}
		//fmt.Printf("middle result of iteration no.%d checking:\n", k) // dbg testing
		//auxio.Quick_check_matrix(params, sk, ctOut, d, d)             // dbg testing
	}
	evaluator.Relinearize(ctOut, ctOut)
	if err != nil {
		return nil, err
	}
	fmt.Printf("%d levels consumed for Matrix Multiplication\n", utils.MinInt(ctMatrixA.Level(), ctMatrixB.Level())-ctOut.Level())
	return
}

// this function has used params.MaxLevel() setting the plaintext's level for internal use, this should be inspect carefully while using bootstrap scheme.
func SquareMatrix_MultBSGS_DecomposeVer_dbg(params ckks.Parameters, sk *rlwe.SecretKey, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctMatrixA *rlwe.Ciphertext, ctMatrixB *rlwe.Ciphertext, SigmaMaps []map[int][]float64, TauMaps []map[int][]float64, sigmaBSGSRatio float64, taoBSGSRatio float64) (ctOut *rlwe.Ciphertext, err error) {
	d := int(math.Sqrt(float64(len(SigmaMaps[0][0]))))
	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	SigmaLTs := make([]ckks.LinearTransform, 2)
	TauLTs := make([]ckks.LinearTransform, 2)
	var ColShiftMap map[int][]float64
	// var RowShiftMap map[int][]float64
	var ctSigmaA_duplicate *rlwe.Ciphertext
	var ctSigmaA_colShift *rlwe.Ciphertext
	// var ctTauB_rowShift *rlwe.Ciphertext

	SigmaLTs[0] = ckks.GenLinearTransformBSGS(encoder, SigmaMaps[0], ctMatrixA.Level(), ctMatrixA.Scale, sigmaBSGSRatio, params.LogSlots())
	SigmaLTs[1] = ckks.GenLinearTransformBSGS(encoder, SigmaMaps[1], ctMatrixA.Level(), ctMatrixA.Scale, sigmaBSGSRatio, params.LogSlots())
	TauLTs[0] = ckks.GenLinearTransformBSGS(encoder, TauMaps[0], ctMatrixB.Level(), ctMatrixB.Scale, taoBSGSRatio, params.LogSlots())
	TauLTs[1] = ckks.GenLinearTransformBSGS(encoder, TauMaps[1], ctMatrixB.Level(), ctMatrixB.Scale, taoBSGSRatio, params.LogSlots())

	now := time.Now() // dbg testing
	ctSigmaA := evaluator.LinearTransformNew(ctMatrixA, SigmaLTs[0])[0]
	ctSigmaA = evaluator.LinearTransformNew(ctSigmaA, SigmaLTs[1])[0]
	evaluator.Rescale(ctSigmaA, ctMatrixA.Scale, ctSigmaA)
	// auxio.Quick_check_matrix(params, sk, ctSigmaA, d, d) // dbg testing
	// ctSigmaA_duplicate = ctSigmaA.CopyNew()
	// auxio.Quick_check_matrix(params, sk, ctSigmaA_duplicate, d, d) // dbg testing
	ctTauB := evaluator.LinearTransformNew(ctMatrixB, TauLTs[0])[0]
	ctTauB = evaluator.LinearTransformNew(ctTauB, TauLTs[1])[0]
	evaluator.Rescale(ctTauB, ctMatrixB.Scale, ctTauB)
	fmt.Printf("BSGS Decomposed LT consume %s \n", time.Since(now))
	// auxio.Quick_check_matrix(params, sk, ctTauB, d, d) // dbg testing

	ctSum := evaluator.MulNew(ctSigmaA, ctTauB)
	// auxio.Quick_check_matrix(params, sk, ctSum, d, d) // dbg testing

	for k := 1; k < d; k++ { // Method 1
		ColShiftMap, _ = Gen_colShift_diagonalVectors(d, k)
		// RowShiftMap, _ = Gen_rowShift_diagonalVectors(d, k)
		ColShiftDiagVec_k := encoder.EncodeNew(ColShiftMap[k], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		ColShiftDiagVec_kmd := encoder.EncodeNew(ColShiftMap[k-d], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		// RowShiftDiagVec := encoder.EncodeNew(RowShiftMap[k], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		// --------------
		ctSigmaA = evaluator.RotateNew(ctSigmaA, 1)                    // option 1
		ctSigmaA_duplicate = evaluator.RotateNew(ctSigmaA, (-d + d*d)) // option 1 // Should have better solution for rotation -d
		ctSigmaA_colShift = evaluator.MulNew(ctSigmaA, ColShiftDiagVec_k)
		// auxio.Quick_check_matrix(params, sk, ctSigmaA_colShift, d, d) // dbg testing
		evaluator.MulAndAdd(ctSigmaA_duplicate, ColShiftDiagVec_kmd, ctSigmaA_colShift)
		// auxio.Quick_check_matrix(params, sk, ctSigmaA_colShift, d, d) // dbg testing
		// --------------
		ctTauB = evaluator.RotateNew(ctTauB, d)
		// ctTauB_rowShift = evaluator.MulNew(ctTauB, RowShiftDiagVec)
		// auxio.Quick_check_matrix(params, sk, ctTauB, d, d) // dbg testing
		// --------------
		evaluator.MulAndAdd(ctSigmaA_colShift, ctTauB, ctSum)
		// auxio.Quick_check_matrix(params, sk, ctSum, d, d) // dbg testing
	}

	evaluator.Relinearize(ctSum, ctSum)
	err = evaluator.Rescale(ctSum, params.DefaultScale(), ctSum)
	if err != nil {
		return nil, err
	}
	ctOut = ctSum
	fmt.Printf("%d levels consumed for Matrix Multiplication\n", utils.MinInt(ctMatrixA.Level(), ctMatrixB.Level())-ctOut.Level())
	return
}

func SquareMatrix_MultBSGS_DecomposeVer_NeedLTs_dbg(params ckks.Parameters, sk *rlwe.SecretKey, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctMatrixA *rlwe.Ciphertext, ctMatrixB *rlwe.Ciphertext, SigmaLTs []ckks.LinearTransform, TauLTs []ckks.LinearTransform, ColShiftPTs map[int][]*rlwe.Plaintext, d int) (ctOut *rlwe.Ciphertext, err error) {

	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})

	ctSigmaA := evaluator.LinearTransform4ArithmeticSeqNew(ctMatrixA, SigmaLTs[0])[0]
	ctSigmaA = evaluator.LinearTransform4ArithmeticSeqNew(ctSigmaA, SigmaLTs[1])[0]
	evaluator.Rescale(ctSigmaA, ctMatrixA.Scale, ctSigmaA)
	// fmt.Printf("ctSigmaA has scale %f, level %d\n", math.Log2(ctSigmaA.Scale.Float64()), ctSigmaA.Level())
	// auxio.Quick_check_matrix(params, sk, ctSigmaA, d, d) // dbg testing
	// ctSigmaA_duplicate = ctSigmaA.CopyNew()
	// auxio.Quick_check_matrix(params, sk, ctSigmaA_duplicate, d, d) // dbg testing
	ctTauB := evaluator.LinearTransform4ArithmeticSeqNew(ctMatrixB, TauLTs[0])[0]
	ctTauB = evaluator.LinearTransform4ArithmeticSeqNew(ctTauB, TauLTs[1])[0]
	evaluator.Rescale(ctTauB, ctMatrixB.Scale, ctTauB)
	// fmt.Printf("ctTauB has scale %f, level %d\n", math.Log2(ctTauB.Scale.Float64()), ctTauB.Level())
	// auxio.Quick_check_matrix(params, sk, ctTauB, d, d)

	ctSum := evaluator.MulNew(ctSigmaA, ctTauB)
	// fmt.Printf("ctSum has scale %f, level %d\n", math.Log2(ctSum.Scale.Float64()), ctSum.Level())
	// auxio.Quick_check_matrix(params, sk, ctSum, d, d) // dbg testing

	/*
		// Method 1: require only one ~ two keys.
		encoder := ckks.NewEncoder(params)
		var ColShiftMap map[int][]float64
		// var RowShiftMap map[int][]float64
		var ctSigmaA_duplicate *rlwe.Ciphertext
		var ctSigmaA_colShift *rlwe.Ciphertext
		for k := 1; k < d; k++ {
			ColShiftMap, _ = Gen_colShift_diagonalVectors(d, k)
			// RowShiftMap, _ = Gen_rowShift_diagonalVectors(d, k)
			ColShiftDiagVec_k := encoder.EncodeNew(ColShiftMap[k], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
			ColShiftDiagVec_kmd := encoder.EncodeNew(ColShiftMap[k-d], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
			// RowShiftDiagVec := encoder.EncodeNew(RowShiftMap[k], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
			// --------------
			ctSigmaA = evaluator.RotateNew(ctSigmaA, 1)                    // option 1
			ctSigmaA_duplicate = evaluator.RotateNew(ctSigmaA, (-d + d*d)) // option 1 // Should have better solution for rotation -d
			ctSigmaA_colShift = evaluator.MulNew(ctSigmaA, ColShiftDiagVec_k)
			// auxio.Quick_check_matrix(params, sk, ctSigmaA_colShift, d, d) // dbg testing
			evaluator.MulAndAdd(ctSigmaA_duplicate, ColShiftDiagVec_kmd, ctSigmaA_colShift)
			// auxio.Quick_check_matrix(params, sk, ctSigmaA_colShift, d, d) // dbg testing
			// --------------
			ctTauB = evaluator.RotateNew(ctTauB, d)
			// ctTauB_rowShift = evaluator.MulNew(ctTauB, RowShiftDiagVec)
			// auxio.Quick_check_matrix(params, sk, ctTauB, d, d) // dbg testing
			// --------------
			evaluator.MulAndAdd(ctSigmaA_colShift, ctTauB, ctSum)

			auxio.Quick_check_matrix(params, sk, ctSum, d, d) // dbg testing
		}
	*/

	// Method2 : using Hoisting
	// Create Rotations for ColShifts
	Rotations4Colshift := make([]int, SigmaLTs[0].N1) // we will do N1 rotations in one hoistedRotation.
	for i := range Rotations4Colshift {
		Rotations4Colshift[i] = i + 1 // generate Rotations 1~N1
	}
	Rotations4Colshift = append(Rotations4Colshift, -d+d*d) // generate an extra Rotation -d+d*d for rotation 1-d ~ N1-d mod d*d
	// Create Rotations for RowShifts
	Rotations4Rowshift := make([]int, TauLTs[0].N1) // we will do N1 rotations in one hoistedRotation
	for i := range Rotations4Rowshift {
		Rotations4Rowshift[i] = (i + 1) * d // generate Rotations 1*d ~ N1*d
	}
	if SigmaLTs[0].N1 >= TauLTs[0].N1 {
		var Rotations4Rowshift_temp []int
		for k := 0; k < d; k += SigmaLTs[0].N1 {
			if k+SigmaLTs[0].N1 >= d {
				Rotations4Colshift = append(Rotations4Colshift[:d-k-1], Rotations4Colshift[len(Rotations4Colshift)-1]) // Rotations4Colshift[d-k-2] contains rotation step d-k-1
				// this will constraint Rotations4ColShift into 1~d-k-1, where k+d-k-1 = d-1 is the last rotation.
			}
			/*
				fmt.Printf("\nCurrent Col Rots when k = %d: ", k)
				for _, rot := range Rotations4Colshift[:(len(Rotations4Colshift) - 1)] {
					fmt.Printf("%d ", rot)
				} // dbg testing
			*/
			ctSigmaA_colShift_list1 := evaluator.RotateHoistedNew(ctSigmaA, Rotations4Colshift)                                                                                           // Generate ctSigmaA_colShift_list1 = Rotate(ctSigmaA;1,2,3,...,N1-1,-d+d*d)
			ctSigmaA_colShift_list2 := evaluator.RotateHoistedNew(ctSigmaA_colShift_list1[Rotations4Colshift[len(Rotations4Colshift)-1]], Rotations4Colshift[:len(Rotations4Colshift)-1]) // Generate ctSigma_colShift_list2 = Rotate(Rotate(ctSigma;-d+d*d);1,2,3,...,N-1)
			// first store the new ctSigmaA = Rotate(ctSigmaA;N1)
			ctSigmaA = ctSigmaA_colShift_list1[Rotations4Colshift[(len(Rotations4Colshift)-1)-1]]

			// fmt.Printf("\nStart computing row shift, using ColShiftPTs: ")
			for _, rot := range Rotations4Colshift[:len(Rotations4Colshift)-1] {
				ctSigmaA_colShift_list1[rot] = evaluator.AddNew(evaluator.MulNew(ctSigmaA_colShift_list1[rot], ColShiftPTs[k+rot][0]), evaluator.MulNew(ctSigmaA_colShift_list2[rot], ColShiftPTs[k+rot][1]))
				err = evaluator.Rescale(ctSigmaA_colShift_list1[rot], params.DefaultScale(), ctSigmaA_colShift_list1[rot])
				if err != nil {
					panic(err)
				}
				// fmt.Printf("%d ", k+rot)
			}
			for j := 0; j < len(Rotations4Colshift)-1; j += TauLTs[0].N1 {
				if j+TauLTs[0].N1 >= len(Rotations4Colshift)-1 { // Since SigmaLTs[0].N1 >= TauLTs[0].N1 then j > 0
					Rotations4Rowshift_temp = Rotations4Rowshift[:(len(Rotations4Colshift)-1)-j] // Rotations4Rowshift[len(Rotations4Colshift)-j-2] contains rotation step len(Rotations4Colshift)-j-1
				} else {
					Rotations4Rowshift_temp = Rotations4Rowshift
				}
				/*
					fmt.Printf("\nCurrent Row Rots /d when k = %d, j= %d: ", k, j)
					for _, rot := range Rotations4Rowshift_temp {
						fmt.Printf("%d ", rot/d)
					} // dbg testing
				*/

				ctTauB_rowShift_list := evaluator.RotateHoistedNew(ctTauB, Rotations4Rowshift_temp) // Generate ctTauB_rowShift_list = Rotate(ctTauB;1*d,2*d,...,(N1-1)*d)

				/*
					fmt.Printf("\nDo MulandAdd: ")
					ctTauB_rowShift_list_idx := make([]int, 0)
					for i, _ := range ctTauB_rowShift_list {
						ctTauB_rowShift_list_idx = append(ctTauB_rowShift_list_idx, i)
					}
					sort.Ints(ctTauB_rowShift_list_idx)
					for _, i := range ctTauB_rowShift_list_idx {
						// fmt.Printf("This round's ctSigmaAColShift:\n")
						// auxio.Quick_check_matrix(params, sk, ctSigmaA_colShift_list1[i/d], d, d)
						// fmt.Printf("This round's ctTauBRowShift:\n")
						// auxio.Quick_check_matrix(params, sk, ctTauB_rowShift_list[i], d, d)
						evaluator.MulAndAdd(ctTauB_rowShift_list[i], ctSigmaA_colShift_list1[i/d], ctSum)
						fmt.Printf("Get ctSum: \n")
						auxio.Quick_check_matrix_full(params, sk, ctSum, d, d)
					}
				*/

				for i, ctTauB_rowShift := range ctTauB_rowShift_list {
					evaluator.SetScale(ctTauB_rowShift, ctSigmaA_colShift_list1[i/d].Scale)
					ctTemp := evaluator.MulNew(ctTauB_rowShift, ctSigmaA_colShift_list1[i/d])
					// fmt.Printf("ctTemp has scale %f, level %d\n", math.Log2(ctTemp.Scale.Float64()), ctTemp.Level())
					if k == 0 {

						evaluator.SetScale(ctSum, ctTemp.Scale)
						// fmt.Printf("For the first time ctSum will be reset to scale %f, level %d\n", math.Log2(ctSum.Scale.Float64()), ctSum.Level())
					}
					evaluator.Add(ctSum, ctTemp, ctSum)
					// evaluator.MulAndAdd(ctTauB_rowShift, ctSigmaA_colShift_list1[i/d], ctSum)
					// fmt.Printf("%d: [%d,%d] ", count, i, i/d)

				}

				// stotre the new ctTauB = Rotate(ctTauB;N1*d)
				ctTauB = ctTauB_rowShift_list[Rotations4Rowshift_temp[len(Rotations4Rowshift_temp)-1]]
			}
		}
	} else { // Unfinished...
		for j := 0; j < d; j += TauLTs[0].N1 {
			if j+TauLTs[0].N1 >= d {
				Rotations4Rowshift = Rotations4Rowshift[:d-j]
			}
			ctTauB_rowShift_list := evaluator.RotateHoistedNew(ctTauB, Rotations4Rowshift)
			for k := 0; k < len(Rotations4Rowshift); k += SigmaLTs[0].N1 {
				if k+SigmaLTs[0].N1 >= len(Rotations4Rowshift) {
					Rotations4Colshift = append(Rotations4Colshift[:len(Rotations4Rowshift)-k], Rotations4Colshift[len(Rotations4Colshift)-1])
				}
				ctSigmaA_colShift_list1 := evaluator.RotateHoistedNew(ctSigmaA, Rotations4Colshift)                                                                       // Generate ctSigmaA_colShift_list1 = Rotate(ctSigmaA;1,2,3,...,N1-1,-d+d*d)
				ctSigmaA_colShift_list2 := evaluator.RotateHoistedNew(ctSigmaA_colShift_list1[len(Rotations4Colshift)-1], Rotations4Colshift[:len(Rotations4Colshift)-1]) // Generate ctSigma_colShift_list2 = Rotate(Rotate(ctSigma;-d+d*d);1,2,3,...,N-1)
				for _, rot := range Rotations4Colshift[:len(Rotations4Colshift)-1] {
					ctSigmaA_colShift_list1[rot] = evaluator.AddNew(evaluator.MulNew(ctSigmaA_colShift_list1[rot], ColShiftPTs[(j+k)+rot][0]), evaluator.MulNew(ctSigmaA_colShift_list2[rot], ColShiftPTs[(j+k)+rot][1]))
					evaluator.MulAndAdd(ctTauB_rowShift_list[(k+rot)*d], ctSigmaA_colShift_list1[rot], ctSum)
				}
			}
		}
	}

	evaluator.Relinearize(ctSum, ctSum)
	err = evaluator.Rescale(ctSum, params.DefaultScale(), ctSum)
	if err != nil {
		return nil, err
	}
	ctOut = ctSum
	fmt.Printf("%d levels consumed for Matrix Multiplication\n", utils.MinInt(ctMatrixA.Level(), ctMatrixB.Level())-ctOut.Level())
	return
}

// incorrect output...
func SquareMatrix_MultBSGS_DecomposeVer2_NeedLTs_dbg(params ckks.Parameters, sk *rlwe.SecretKey, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctMatrixA *rlwe.Ciphertext, ctMatrixB []*rlwe.Ciphertext, DSigmaLTs1 []ckks.LinearTransform, DSigmaLTs2 []ckks.LinearTransform, DTaulTs []ckks.LinearTransform, ColShiftPTs interface{}, d int, goroutineNum ...int) (ctOut []*rlwe.Ciphertext, err error) {

	var ctDSigmaA []*rlwe.Ciphertext
	ctOut = make([]*rlwe.Ciphertext, len(ctMatrixB))

	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	// Computing SigmaA:
	LTstemp := make([]ckks.LinearTransform, 2)
	LTstemp[0] = DSigmaLTs1[len(DSigmaLTs1)-1]
	LTstemp[1] = DSigmaLTs2[len(DSigmaLTs2)-1]

	//auxio.Quick_check_matrix(params, sk, ctMatrixA, d, d)
	ctDSigmaA = evaluator.LinearTransform4ArithmeticSeqNew(ctMatrixA, LTstemp, true)
	err = evaluator.Rescale(ctDSigmaA[0], params.DefaultScale(), ctDSigmaA[0])
	if err != nil {
		return nil, err
	}
	err = evaluator.Rescale(ctDSigmaA[1], params.DefaultScale(), ctDSigmaA[1])
	if err != nil {
		return nil, err
	}
	for i := len(DSigmaLTs1) - 2; i >= 0; i-- {
		ctDSigmaA[0] = evaluator.LinearTransform4ArithmeticSeqNew(ctDSigmaA[0], DSigmaLTs1[i])[0]
		err = evaluator.Rescale(ctDSigmaA[0], params.DefaultScale(), ctDSigmaA[0])
	}
	for i := len(DSigmaLTs2) - 2; i >= 0; i-- {
		ctDSigmaA[1] = evaluator.LinearTransform4ArithmeticSeqNew(ctDSigmaA[1], DSigmaLTs2[i])[0]
		err = evaluator.Rescale(ctDSigmaA[1], params.DefaultScale(), ctDSigmaA[1])
	}
	evaluator.SetScale(ctDSigmaA[1], ctDSigmaA[0].Scale)
	evaluator.Add(ctDSigmaA[0], ctDSigmaA[1], ctDSigmaA[0])
	// auxio.Quick_check_matrix(params, sk, ctDSigmaA[0], d, d)

	var ctDSigmaA_ColShifted []*rlwe.Ciphertext
	switch ColShiftMask := ColShiftPTs.(type) {
	case []ckks.LinearTransform:
		ctDSigmaA_ColShifted = evaluator.ColShiftRestricted(ctDSigmaA[0], ColShiftMask, LTstemp[0].N1, d)

		for i := range ctDSigmaA_ColShifted {
			err = evaluator.Rescale(ctDSigmaA_ColShifted[i], params.DefaultScale(), ctDSigmaA_ColShifted[i])
			if err != nil {
				return nil, err
			}
			//auxio.Quick_check_matrix(params, sk, ctDSigmaA_ColShifted[i], d, d)
		}

	}

	// Computing TauBs:
	DTauN1 := DTaulTs[len(DTaulTs)-1].N1
	rots := make([]int, DTauN1)
	for i := range rots {
		rots[i] = (i + 1) * d
	}
	//var ctDTauB_RowShifted map[int]*rlwe.Ciphertext
	//var ctSum *rlwe.Ciphertext
	var wg sync.WaitGroup
	var maxSubRoutes int
	if len(goroutineNum) > 0 {
		maxSubRoutes = goroutineNum[0]
	} else {
		maxSubRoutes = 1
	}
	var subrouteNum = maxSubRoutes
	for i := 0; i < len(ctMatrixB); i++ {
		wg.Add(1)
		subrouteNum--
		go func(i int) (err error) {
			defer wg.Done()
			var ctDTauB_RowShifted map[int]*rlwe.Ciphertext
			evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: galk, Rlk: rlk})
			ctDTauB := evaluator.LinearTransform4ArithmeticSeqNew(ctMatrixB[i], DTaulTs[len(DTaulTs)-1])
			err = evaluator.Rescale(ctDTauB[0], params.DefaultScale(), ctDTauB[0])
			if err != nil {
				return err
			}
			for j := len(DTaulTs) - 2; j >= 0; j-- {
				ctDTauB[0] = evaluator.LinearTransform4ArithmeticSeqNew(ctDTauB[0], DTaulTs[j])[0]
				err = evaluator.Rescale(ctDTauB[0], params.DefaultScale(), ctDTauB[0])
				if err != nil {
					return err
				}
			}
			// auxio.Quick_check_matrix(params, sk, ctDTauB[0], d, d)
			ctSum := evaluator.MulNew(ctDSigmaA[0], ctDTauB[0])
			auxio.Quick_check_infos(ctSum, "ctSum")
			auxio.Quick_check_matrix(params, sk, ctSum, d, d)
			for k := 0; k < d-1; k += DTauN1 {
				if k+DTauN1 > d-1 {
					ctDTauB_RowShifted = evaluator.RotateHoistedNew(ctDTauB[0], rots[0:d-1-k])
				} else {
					ctDTauB_RowShifted = evaluator.RotateHoistedNew(ctDTauB[0], rots)
				}
				// Sorted way

				ctDTauB_RowShifted_list := make([]int, 0)
				for i, _ := range ctDTauB_RowShifted {
					ctDTauB_RowShifted_list = append(ctDTauB_RowShifted_list, i)
				}
				sort.Ints(ctDTauB_RowShifted_list)
				for _, i := range ctDTauB_RowShifted_list {
					key := i / d
					auxio.Quick_check_infos(ctDSigmaA_ColShifted[k+key-1], "ctA")
					auxio.Quick_check_matrix(params, sk, ctDSigmaA_ColShifted[k+key-1], d, d)
					auxio.Quick_check_infos(ctDTauB_RowShifted[i], "ctB")
					auxio.Quick_check_matrix(params, sk, ctDTauB_RowShifted[i], d, d)
					ctTemp := evaluator.MulNew(ctDTauB_RowShifted[i], ctDSigmaA_ColShifted[k+key-1])
					auxio.Quick_check_infos(ctTemp, "ctTemp")
					auxio.Quick_check_matrix(params, sk, ctTemp, d, d)
					if k == 0 {
						ctSum.SetScale(ctTemp.Scale)
					}
					evaluator.Add(ctSum, ctTemp, ctSum)
					auxio.Quick_check_infos(ctSum, "ctSum")
					auxio.Quick_check_matrix(params, sk, ctSum, d, d)
				}

				// Map way
				/*
					for key, ct := range ctDTauB_RowShifted {
						key = key / d
						auxio.Quick_check_infos(ctDSigmaA_ColShifted[k+key-1], "ctA")
						auxio.Quick_check_matrix(params, sk, ctDSigmaA_ColShifted[k+key-1], d, d)
						auxio.Quick_check_infos(ct, "ctB")
						auxio.Quick_check_matrix(params, sk, ct, d, d)
						ctTemp := evaluator.MulNew(ctDSigmaA_ColShifted[k+key-1], ct)
						auxio.Quick_check_infos(ctTemp, "ctTemp")
						auxio.Quick_check_matrix(params, sk, ctTemp, d, d)
						if k == 0 {
							ctSum.SetScale(ctTemp.Scale)
						}
						evaluator.Add(ctSum, ctTemp, ctSum)
						auxio.Quick_check_infos(ctSum, "ctSum")
						auxio.Quick_check_matrix(params, sk, ctSum, d, d)

					}
				*/
				ctDTauB[0] = ctDTauB_RowShifted[d*len(ctDTauB_RowShifted)].CopyNew()
			}
			evaluator.Relinearize(ctSum, ctSum)
			err = evaluator.Rescale(ctSum, params.DefaultScale(), ctSum)
			if err != nil {
				return err
			}
			ctOut[i] = ctSum.CopyNew()
			fmt.Printf("subroutine %d-th for matrix multiplication Done.\n", i)
			return nil
		}(i)
		if subrouteNum <= 0 {
			wg.Wait()
			subrouteNum = maxSubRoutes
		}
	}
	if subrouteNum != maxSubRoutes {
		wg.Wait()
	}

	return
}

func SquareMatrix_MultBSGS_DecomposeVer3_NeedLTs_dbg(params ckks.Parameters, sk *rlwe.SecretKey, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctMatrixAdecomp []*rlwe.Ciphertext, ctMatrixB []*rlwe.Ciphertext, DTaulTs []ckks.LinearTransform, d int, goroutineNum ...int) (ctOut []*rlwe.Ciphertext, err error) {
	ctOut = make([]*rlwe.Ciphertext, len(ctMatrixB))
	// Computing TauBs:
	DTauN1 := DTaulTs[len(DTaulTs)-1].N1
	rots := make([]int, DTauN1)
	for i := range rots {
		rots[i] = (i + 1) * d
	}
	//var ctDTauB_RowShifted map[int]*rlwe.Ciphertext
	//var ctSum *rlwe.Ciphertext
	var wg sync.WaitGroup
	var maxSubRoutes int
	if len(goroutineNum) > 0 {
		maxSubRoutes = goroutineNum[0]
	} else {
		maxSubRoutes = 1
	}
	var subrouteNum = maxSubRoutes
	for i := 0; i < len(ctMatrixB); i++ {
		wg.Add(1)
		subrouteNum--
		go func(i int) (err error) {
			defer wg.Done()
			var ctDTauB_RowShifted map[int]*rlwe.Ciphertext
			evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: galk, Rlk: rlk})
			ctDTauB := evaluator.LinearTransform4ArithmeticSeqNew(ctMatrixB[i], DTaulTs[len(DTaulTs)-1])
			err = evaluator.Rescale(ctDTauB[0], params.DefaultScale(), ctDTauB[0])
			if err != nil {
				return err
			}
			for j := len(DTaulTs) - 2; j >= 0; j-- {
				ctDTauB[0] = evaluator.LinearTransform4ArithmeticSeqNew(ctDTauB[0], DTaulTs[j])[0]
				err = evaluator.Rescale(ctDTauB[0], params.DefaultScale(), ctDTauB[0])
				if err != nil {
					return err
				}
			}
			// auxio.Quick_check_matrix(params, sk, ctDTauB[0], d, d)
			ctSum := evaluator.MulNew(ctMatrixAdecomp[0], ctDTauB[0])
			// auxio.Quick_check_matrix(params, sk, ctSum, d, d)
			for k := 0; k < d-1; k += DTauN1 {
				if k+DTauN1 > d-1 {
					ctDTauB_RowShifted = evaluator.RotateHoistedNew(ctDTauB[0], rots[0:d-1-k])
				} else {
					ctDTauB_RowShifted = evaluator.RotateHoistedNew(ctDTauB[0], rots)
				}
				// Sorted way
				/*
					ctDTauB_RowShifted_list := make([]int, 0)
					for i, _ := range ctDTauB_RowShifted {
						ctDTauB_RowShifted_list = append(ctDTauB_RowShifted_list, i)
					}
					sort.Ints(ctDTauB_RowShifted_list)
					for _, i := range ctDTauB_RowShifted_list {
						key := i / d
						ctTemp := evaluator.MulNew(ctDTauB_RowShifted[i], ctDSigmaA_ColShifted[k+key-1])
						if k == 0 {
							ctSum.SetScale(ctTemp.Scale)
						}
						evaluator.Add(ctSum, ctTemp, ctSum)
						auxio.Quick_check_matrix(params, sk, ctSum, d, d)
					}
				*/
				// Map way

				for key, ct := range ctDTauB_RowShifted {
					key = key / d
					ctTemp := evaluator.MulNew(ctMatrixAdecomp[k+key], ct) // ctTemp := evaluator.MulNew(ctMatrixAdecomp[k+key-1], ct)
					/*
						if k == 0 {
							ctSum.SetScale(ctTemp.Scale)
						}
					*/
					evaluator.Add(ctSum, ctTemp, ctSum)
				}
				ctDTauB[0] = ctDTauB_RowShifted[d*len(ctDTauB_RowShifted)].CopyNew()
			}
			evaluator.Relinearize(ctSum, ctSum)
			err = evaluator.Rescale(ctSum, params.DefaultScale(), ctSum)
			if err != nil {
				return err
			}
			ctOut[i] = ctSum.CopyNew()
			fmt.Printf("subroutine %d-th for matrix multiplication Done.\n", i)
			return nil
		}(i)
		if subrouteNum <= 0 {
			wg.Wait()
			subrouteNum = maxSubRoutes
		}
	}
	if subrouteNum != maxSubRoutes {
		wg.Wait()
	}
	return
}

func SquareMatrix_MultBSGS_DecomposeVer4_NeedLTs_dbg(params ckks.Parameters, sk *rlwe.SecretKey, rlk *rlwe.RelinearizationKey, ctMatrixAdecomp []*rlwe.Ciphertext, ctMatrixBdecomp [][]*rlwe.Ciphertext, d int, goroutineNum ...int) (ctOut []*rlwe.Ciphertext, err error) {
	ctOut = make([]*rlwe.Ciphertext, len(ctMatrixBdecomp))
	var wg sync.WaitGroup
	var maxSubRoutes int
	if len(goroutineNum) > 0 {
		maxSubRoutes = goroutineNum[0]
	} else {
		maxSubRoutes = 1
	}
	var subrouteNum = maxSubRoutes

	for i := 0; i < len(ctMatrixBdecomp); i++ {
		wg.Add(1)
		subrouteNum--
		go func(i int) {
			defer wg.Done()
			evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
			var ctSum *rlwe.Ciphertext
			for key, ct := range ctMatrixBdecomp[i] {
				//auxio.Quick_check_infos(ctMatrixAdecomp[key], "ctA")
				//auxio.Quick_check_matrix(params, sk, ctMatrixAdecomp[key], d, d)
				//auxio.Quick_check_infos(ct, "ctB")
				//auxio.Quick_check_matrix(params, sk, ct, d, d)
				ctTemp := evaluator.MulNew(ctMatrixAdecomp[key], ct)
				//auxio.Quick_check_infos(ctTemp, "ctTemp")
				//auxio.Quick_check_matrix(params, sk, ctTemp, d, d)
				if key == 0 {
					ctSum = ctTemp.CopyNew()
				} else {
					evaluator.Add(ctSum, ctTemp, ctSum)
				}
				//auxio.Quick_check_infos(ctSum, "ctSum")
				//auxio.Quick_check_matrix(params, sk, ctSum, d, d)
			}
			evaluator.Relinearize(ctSum, ctSum)
			evaluator.Rescale(ctSum, params.DefaultScale(), ctSum)
			ctOut[i] = ctSum.CopyNew()
			fmt.Printf("subroutine %d-th for matrix multiplication Done.\n", i)
		}(i)
		if subrouteNum <= 0 {
			wg.Wait()
			subrouteNum = maxSubRoutes
		}
	}
	if subrouteNum != maxSubRoutes {
		wg.Wait()
	}
	return
}

func Sigma_linearTransform_MultiThread_dbg(params ckks.Parameters, sk *rlwe.SecretKey, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctMatrixA []*rlwe.Ciphertext, DSigmaLTs1 []ckks.LinearTransform, DSigmaLTs2 []ckks.LinearTransform, ColShiftPTs interface{}, d int, goroutineNum ...int) (ctOut [][]*rlwe.Ciphertext, err error) {

	var wg sync.WaitGroup
	var maxSubRoutes int
	if len(goroutineNum) > 0 {
		maxSubRoutes = goroutineNum[0]
	} else {
		maxSubRoutes = 1
	}
	ctOut = make([][]*rlwe.Ciphertext, len(ctMatrixA))
	LTstemp := make([]ckks.LinearTransform, 2)
	LTstemp[0] = DSigmaLTs1[len(DSigmaLTs1)-1]
	LTstemp[1] = DSigmaLTs2[len(DSigmaLTs2)-1]
	var subrouteNum = maxSubRoutes
	for i := 0; i < len(ctMatrixA); i++ {
		wg.Add(1)
		subrouteNum--
		go func(i int) {
			defer wg.Done()
			evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
			ctDSigmaA := evaluator.LinearTransform4ArithmeticSeqNew(ctMatrixA[i], LTstemp, true)
			err = evaluator.Rescale(ctDSigmaA[0], params.DefaultScale(), ctDSigmaA[0])
			err = evaluator.Rescale(ctDSigmaA[1], params.DefaultScale(), ctDSigmaA[1])
			for i := len(DSigmaLTs1) - 2; i >= 0; i-- {
				ctDSigmaA[0] = evaluator.LinearTransform4ArithmeticSeqNew(ctDSigmaA[0], DSigmaLTs1[i])[0]
				err = evaluator.Rescale(ctDSigmaA[0], params.DefaultScale(), ctDSigmaA[0])
			}
			for i := len(DSigmaLTs2) - 2; i >= 0; i-- {
				ctDSigmaA[1] = evaluator.LinearTransform4ArithmeticSeqNew(ctDSigmaA[1], DSigmaLTs2[i])[0]
				err = evaluator.Rescale(ctDSigmaA[1], params.DefaultScale(), ctDSigmaA[1])
			}
			evaluator.SetScale(ctDSigmaA[1], ctDSigmaA[0].Scale)
			evaluator.Add(ctDSigmaA[0], ctDSigmaA[1], ctDSigmaA[0])
			// auxio.Quick_check_matrix(params, sk, ctDSigmaA[0], d, d)
			switch ColShiftMask := ColShiftPTs.(type) {
			case []ckks.LinearTransform:
				ctOut[i] = evaluator.ColShiftRestricted(ctDSigmaA[0], ColShiftMask, LTstemp[0].N1, d)
				for j := range ctOut[i] {
					err = evaluator.Rescale(ctOut[i][j], params.DefaultScale(), ctOut[i][j])
					//auxio.Quick_check_matrix(params, sk, ctDSigmaA_ColShifted[i], d, d)
				}
				evaluator.SetScale(ctDSigmaA[0], ctOut[i][0].Scale)
				ctOut[i] = append([]*rlwe.Ciphertext{ctDSigmaA[0]}, ctOut[i]...) // insert the non-shifted SigmaA as the first element of ctOut
			}
		}(i)
		if subrouteNum <= 0 {
			wg.Wait()
			subrouteNum = maxSubRoutes
		}
	}
	wg.Wait()

	return
}

func Tau_linearTransform_MultiThread_dbg(params ckks.Parameters, sk *rlwe.SecretKey, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctMatrixB []*rlwe.Ciphertext, DTauLTs []ckks.LinearTransform, d int, goroutineNum ...int) (ctOut [][]*rlwe.Ciphertext, err error) {
	DTauN1 := DTauLTs[len(DTauLTs)-1].N1
	rots := make([]int, DTauN1)
	for i := range rots {
		rots[i] = (i + 1) * d
	}
	var wg sync.WaitGroup
	var maxSubRoutes int
	if len(goroutineNum) > 0 {
		maxSubRoutes = goroutineNum[0]
	} else {
		maxSubRoutes = 1
	}
	ctOut = make([][]*rlwe.Ciphertext, len(ctMatrixB))
	for i := 0; i < len(ctMatrixB); i++ {
		ctOut[i] = make([]*rlwe.Ciphertext, d)
	}
	var subrouteNum = maxSubRoutes
	for i := 0; i < len(ctMatrixB); i++ {
		wg.Add(1)
		subrouteNum--
		go func(i int) {
			defer wg.Done()
			var ctDTauB_RowShifted map[int]*rlwe.Ciphertext
			evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
			ctDTauB := evaluator.LinearTransform4ArithmeticSeqNew(ctMatrixB[i], DTauLTs[len(DTauLTs)-1])
			err = evaluator.Rescale(ctDTauB[0], params.DefaultScale(), ctDTauB[0])
			for j := len(DTauLTs) - 2; j >= 0; j-- {
				ctDTauB[0] = evaluator.LinearTransform4ArithmeticSeqNew(ctDTauB[0], DTauLTs[j])[0]
				err = evaluator.Rescale(ctDTauB[0], params.DefaultScale(), ctDTauB[0])
			}
			ctOut[i][0] = ctDTauB[0].CopyNew()

			for k := 0; k < d-1; k += DTauN1 {
				if k+DTauN1 > d-1 {
					ctDTauB_RowShifted = evaluator.RotateHoistedNew(ctDTauB[0], rots[0:d-1-k])
				} else {
					ctDTauB_RowShifted = evaluator.RotateHoistedNew(ctDTauB[0], rots)
				}
				for key := range ctDTauB_RowShifted {
					ctOut[i][k+key/d] = ctDTauB_RowShifted[key].CopyNew()
				}
				ctDTauB[0] = ctDTauB_RowShifted[d*len(ctDTauB_RowShifted)].CopyNew()
			}
		}(i)
		if subrouteNum <= 0 {
			wg.Wait()
			subrouteNum = maxSubRoutes
		}
	}
	wg.Wait()

	/*
		for i := range ctOut {
			fmt.Printf("ctOut[%d]:\n", i)
			for j := range ctOut[i] {
				auxio.Quick_check_matrix(params, sk, ctOut[i][j], d, d)
			}
		}
	*/
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

// this function has used params.MaxLevel() setting the plaintext's level for internal use, this should be inspect carefully while using bootstrap scheme.
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
			ctReplicatedA, err = Matrix_rowTotalSum_withColOrder(params, galk, ctMaskedA, m0, m1)
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
			ctReplicatedA, err = Matrix_rowTotalSum_withColOrder(params, galk, ctMaskedA, m0, m1)
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

func Sigma_linearTransformBSGS_dbg(params ckks.Parameters, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, d int, BSGSRatio float64) (ctOut *rlwe.Ciphertext, err error) {
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
	RotationsNeeded := sigmaLT.Rotations()
	fmt.Printf("Rotation Steps Needed for SigmaLT with BSGSRation %f: ", BSGSRatio)
	for _, i := range RotationsNeeded {
		fmt.Printf("%d ", i)
	}
	fmt.Printf("\n")
	ctOut_list := evaluator.LinearTransformNew(ctIn, sigmaLT)
	ctOut = ctOut_list[0]
	err = evaluator.Rescale(ctOut, ctIn.Scale, ctOut)
	if err != nil {
		return nil, err
	}
	fmt.Printf("%d levels consumed for Sigma LinearTransformBSGS\n", ctIn.Level()-ctOut.Level())
	return
}

func Sigma_linearTransform_dbg(params ckks.Parameters, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, d int) (ctOut *rlwe.Ciphertext, err error) {
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

	RotationsNeeded := sigmaLT.Rotations()
	fmt.Printf("Rotation Steps Needed for SigmaLT: ")
	for _, i := range RotationsNeeded {
		fmt.Printf("%d ", i)
	}
	fmt.Printf("\n")

	ctOut_list := evaluator.LinearTransformNew(ctIn, sigmaLT)
	ctOut = ctOut_list[0]
	err = evaluator.Rescale(ctOut, ctIn.Scale, ctOut)
	if err != nil {
		return nil, err
	}
	fmt.Printf("%d levels consumed for Sigma LinearTransform\n", ctIn.Level()-ctOut.Level())
	return
}

func General_linearTransform_dbg(params ckks.Parameters, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, U map[int][]float64, BSGSRatio ...float64) (ctOut *rlwe.Ciphertext, err error) {
	var LT ckks.LinearTransform
	Scale := ctIn.Scale
	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	now := time.Now()
	if len(BSGSRatio) == 0 {
		LT = ckks.GenLinearTransform(encoder, U, ctIn.Level(), Scale, params.LogSlots())
	} else {
		LT = ckks.GenLinearTransformBSGS(encoder, U, ctIn.Level(), Scale, BSGSRatio[0], params.LogSlots())
	}
	fmt.Printf("LinearTransformGeneration takes %s\n", time.Since(now))

	RotationsNeeded := LT.Rotations()
	fmt.Printf("Rotation Steps Needed for LinearTransform: ")
	for _, i := range RotationsNeeded {
		fmt.Printf("%d ", i)
	}
	if len(BSGSRatio) != 0 {
		fmt.Printf(" (BSGSRatio: %f)", BSGSRatio[0])
	}
	fmt.Printf("\n")

	now = time.Now()
	ctOut_list := evaluator.LinearTransformNew(ctIn, LT)
	ctOut = ctOut_list[0]
	err = evaluator.Rescale(ctOut, ctIn.Scale, ctOut)
	if err != nil {
		return nil, err
	}
	fmt.Printf("%d levels consumed for General LinearTransform in %s \n", ctIn.Level()-ctOut.Level(), time.Since(now))
	return
}
