package check

import (
	"fmt"
	"math"
	"sync"
	"time"

	auxio "project1-fhe_extension_v1.0/auxiliary_io"
	mtrxmult "project1-fhe_extension_v1.0/matrix_mult"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func LinTransX_check() {
	var err error
	//	var now time.Time

	d := 128
	Sigma_BSGSN1 := 16
	Tau_BSGSN1 := 16
	DSigma_BSGSN1 := 16
	DTau_BSGSN1 := 16
	Trans_BSGSN1 := 8

	var sk *rlwe.SecretKey
	var pk *rlwe.PublicKey
	var rlk *rlwe.RelinearizationKey
	var galk *rlwe.RotationKeySet

	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: 15,
		Q: []uint64{0x4000000120001, 0x10000140001, 0xffffe80001, // 50 + 17 x 40 =
			0x10000290001, 0xffffc40001, 0x100003e0001,
			0x10000470001, 0x100004b0001, 0xffffb20001,
			//0x10000500001, 0x10000650001, 0xffff940001,
			//0xffff8a0001, 0xffff820001, 0xffff780001,
			//0x10000890001, 0xffff750001, 0x10000960001
		},
		P:            []uint64{0x40000001b0001, 0x3ffffffdf0001, 0x4000000270001}, // 50, 50, 50
		LogSlots:     int(math.Log2(float64(d * d))),
		DefaultScale: 1 << 40,
	})

	/*
		params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
			LogN: 16,
			Q: []uint64{0x80000000080001, 0x2000000a0001, 0x2000000e0001, 0x1fffffc20001, // 55 + 5 x 45
				0x200000440001, 0x200000500001},
			P:            []uint64{0x80000000440001, 0x7fffffffba0001, 0x80000000500001, 0x7fffffffaa0001}, // 4 x 55
			LogSlots:     int(math.Log2(float64(d * d))),
			DefaultScale: 1 << 45,
		})
	*/

	if err != nil {
		panic(err)
	}
	kgen := ckks.NewKeyGenerator(params)
	sk, pk = kgen.GenKeyPair()
	rlk = kgen.GenRelinearizationKey(sk, 1)

	// create a Matrix A
	MatrixA := make([][]float64, d)
	for i := 0; i < d; i++ {
		MatrixA[i] = make([]float64, d)
		for j := 0; j < d; j++ {
			MatrixA[i][j] = float64(i*d + j)
		}
	}
	// auxio.Print_matrix_f64_2d(MatrixA, d, d)
	// encode matrixA into a vector using row ordering
	var a []float64
	a, err = mtrxmult.Row_orderingInv(MatrixA)
	if err != nil {
		panic(err)
	}

	// encode into plaintext and encrypt :
	encoder := ckks.NewEncoder(params)
	encryptor := ckks.NewEncryptor(params, pk)
	ptA := encoder.EncodeNew(a, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
	ctA := encryptor.EncryptNew(ptA)

	/*
		ctAs_shifted := make([]*rlwe.Ciphertext, 0)
		ctAs_shifted = append(ctAs_shifted, ctA)
		for i := 0; i < DSigma_BSGSN1-1; i++ {
			a = append(a[1:], a[0])
			ctAs_shifted = append(ctAs_shifted, encryptor.EncryptNew(encoder.EncodeNew(a, params.MaxLevel(), params.DefaultScale(), params.LogSlots())))
		}
	*/

	// create the Map of Sigma LinearTransformation
	var SigmaDiagonalMap map[int][]float64
	SigmaDiagonalMap, err = mtrxmult.Gen_sigma_diagonalVecotrs(d)
	if err != nil {
		panic(err)
	}
	// create the Map of Tau LinearTransformation
	var TauDiagonalMap map[int][]float64
	TauDiagonalMap, err = mtrxmult.Gen_tao_diagonalVectors(d)
	if err != nil {
		panic(err)
	}

	// create the Decomposed Maps of Sigma LinearTransformation
	SigmaMatrix := mtrxmult.DiagonalVectors2Matrix(SigmaDiagonalMap, d*d)
	DSigmaDiagonalMaps, err := mtrxmult.Converge2DiagonalDecompose_Sigma(SigmaMatrix)
	if err != nil {
		panic(err)
	}
	// create the Decomposed Maps of Tau LinearTransformation

	TauMatrix := mtrxmult.DiagonalVectors2Matrix(TauDiagonalMap, d*d)
	DTauDiagonalMaps, err := mtrxmult.Converge2DiagonalDecompose_Tao(TauMatrix)
	if err != nil {
		panic(err)
	}

	idx, _, _ := ckks.BsgsIndex4ArithmeticSeq(SigmaDiagonalMap, params.Slots(), Sigma_BSGSN1, 1)
	var rots = make(map[int]bool)
	for outidx, inidces := range idx {
		rots[outidx] = true
		for _, i := range inidces {
			rots[i] = true
		}
	}
	fmt.Printf("SigmaLT %d Rotations: \n", len(rots))

	idx, _, _ = ckks.BsgsIndex4ArithmeticSeq(DSigmaDiagonalMaps[0], params.Slots(), Sigma_BSGSN1, 1)
	rots = make(map[int]bool)
	for outidx, inidces := range idx {
		rots[outidx] = true
		for _, i := range inidces {
			rots[i] = true
		}
	}
	fmt.Printf("DSigmaLT %d Rotations: \n", len(rots))

	idx, _, _ = ckks.BsgsIndex4ArithmeticSeq(TauDiagonalMap, params.Slots(), Tau_BSGSN1, d)
	rots = make(map[int]bool)
	for outidx, inidces := range idx {
		rots[outidx] = true
		for _, i := range inidces {
			rots[i] = true
		}
	}
	fmt.Printf("TauLT %d Rotations : \n", len(rots))

	idx, _, _ = ckks.BsgsIndex4ArithmeticSeq(DTauDiagonalMaps[0], params.Slots(), Tau_BSGSN1, d)
	rots = make(map[int]bool)
	for outidx, inidces := range idx {
		rots[outidx] = true
		for _, i := range inidces {
			rots[i] = true
		}
	}
	fmt.Printf("DTauLT %d Rotations : \n", len(rots))

	// create the Map of Transpose LinearTransformation
	/*
		var TransDiagonalMap map[int][]float64
		TransDiagonalMap, err = mtrxmult.Gen_transpose_diagonalVectors(d)
		if err != nil {
			panic(err)
		}
	*/

	// generate LinearTransform
	SigmaLT := ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, SigmaDiagonalMap, params.MaxLevel(), params.DefaultScale(), Sigma_BSGSN1, 1, params.LogSlots())
	TauLT := ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, TauDiagonalMap, params.MaxLevel(), params.DefaultScale(), Tau_BSGSN1, d, params.LogSlots())
	// TransLT := ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, TransDiagonalMap, params.MaxLevel(), params.DefaultScale(), Trans_BSGSN1, (d - 1), params.LogSlots())
	DSigmaLTs := make([]ckks.LinearTransform, 2)
	DTauLTs := make([]ckks.LinearTransform, 2)
	DSigmaLTs[0] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DSigmaDiagonalMaps[0], params.MaxLevel(), params.DefaultScale(), DSigma_BSGSN1, 1, params.LogSlots())
	DSigmaLTs[1] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DSigmaDiagonalMaps[1], params.MaxLevel(), params.DefaultScale(), DSigma_BSGSN1, 1, params.LogSlots())
	DTauLTs[0] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DTauDiagonalMaps[0], params.MaxLevel(), params.DefaultScale(), DTau_BSGSN1, d, params.LogSlots())
	DTauLTs[1] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DTauDiagonalMaps[1], params.MaxLevel(), params.DefaultScale(), DTau_BSGSN1, d, params.LogSlots())

	fmt.Printf("Test configuration: d: %d, Sigma'sBSGSN1: %d, Tau'sBSGSN1: %d, DSigma'sBSGSN1:%d, DTau'sBSGSN1: %d, Tran'sBSGSN1: %d\n", d, Sigma_BSGSN1, Tau_BSGSN1, DSigma_BSGSN1, DTau_BSGSN1, Trans_BSGSN1)

	// generate Rotations

	fmt.Printf("\nSigmaLT %d Rotations: ", len(SigmaLT.Rotations4ArithmeticSeq()))
	for _, i := range SigmaLT.Rotations4ArithmeticSeq() {
		fmt.Printf("%d ", i)
	}
	fmt.Printf("\nDSigmaLT %d Rotations: ", len(DSigmaLTs[0].Rotations4ArithmeticSeq()))
	for _, i := range DSigmaLTs[0].Rotations4ArithmeticSeq() {
		fmt.Printf("%d ", i)
	}

	fmt.Printf("\nTauLT %d Rotations : ", len(TauLT.Rotations4ArithmeticSeq()))
	for _, i := range TauLT.Rotations4ArithmeticSeq() {
		fmt.Printf("%d ", i)
	}

	fmt.Printf("\nDTauLT %d Rotations : ", len(DTauLTs[0].Rotations4ArithmeticSeq()))
	for _, i := range DTauLTs[0].Rotations4ArithmeticSeq() {
		fmt.Printf("%d ", i)
	}

	/*
		fmt.Printf("\nTransLT %d Rotations: ", len(TransLT.Rotations4ArithmeticSeq()))
		for _, i := range TransLT.Rotations4ArithmeticSeq() {
			fmt.Printf("%d ", i)
		}
	*/

	// generate Keys.
	now := time.Now()
	SigmaRotKeys := kgen.GenRotationKeysForRotations(SigmaLT.Rotations4ArithmeticSeq(), false, sk)
	fmt.Printf("\nSigma Rot keys generation consumes %s\n", time.Since(now))
	now = time.Now()
	TauRotKeys := kgen.GenRotationKeysForRotations(TauLT.Rotations4ArithmeticSeq(), false, sk)
	fmt.Printf("\nTau Rot keys generation consumes %s\n", time.Since(now))

	now = time.Now()
	DSigmaRotKeys := kgen.GenRotationKeysForRotations(DSigmaLTs[0].Rotations4ArithmeticSeq(), false, sk)
	fmt.Printf("\nDSigma Rot keys generation consumes %s\n", time.Since(now))

	now = time.Now()
	DTauRotKeys := kgen.GenRotationKeysForRotations(DTauLTs[0].Rotations4ArithmeticSeq(), false, sk)
	fmt.Printf("\nDTau Rot keys generation consumes %s\n", time.Since(now))

	//TransRotKeys := kgen.GenRotationKeysForRotations(TransLT.Rotations4ArithmeticSeq(), false, sk)

	// create evaluatorX
	evaluatorX4Sigma := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: SigmaRotKeys, Rlk: rlk})
	evaluatorX4Tau := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: TauRotKeys, Rlk: rlk})
	// evaluatorX4Trans := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: TransRotKeys, Rlk: rlk})
	// Evaluate LinearTransformation:

	times := 5
	var sumtime time.Duration
	var ctDSigmaA, ctDTauA []*rlwe.Ciphertext
	var ctSigmaA, ctTauA []*rlwe.Ciphertext

	/*
		fmt.Printf("----------------------------- Experiment Check -----------------------------\n")
		for i:=0; i<DSigma_BSGSN1;i++{
			ctDSigmaA = evaluatorX4Sigma.LinearTransform4ArithmeticSeqNew(ctAs_shifted[i], DSigmaLTs[0])
			auxio.Quick_check_matrix(params,sk,ctDSigmaA[0],d,d)
		}
		for i:=0;
	*/

	fmt.Printf("----------------------------- Single Thread Mod -----------------------------\n")

	sumtime = 0
	for i := 0; i < times; i++ {
		now = time.Now()
		ctDSigmaA = evaluatorX4Sigma.LinearTransform4ArithmeticSeqNew(ctA, DSigmaLTs[0])
		elapsed1 := time.Since(now)
		fmt.Printf("iter %dth: DSigmaLT1 consumes %s \n", i, elapsed1)
		now = time.Now()
		ctDSigmaA = evaluatorX4Sigma.LinearTransform4ArithmeticSeqNew(ctDSigmaA[0], DSigmaLTs[1])
		elpased2 := time.Since(now)
		fmt.Printf("iter %dth: DSigmaLT2 consumes %s \n", i, elpased2)
		if i > 0 {
			sumtime += (elapsed1 + elpased2)
		}
	}
	fmt.Printf("DSigmaLT consumes %s in average\n", sumtime/time.Duration(times-1))

	sumtime = 0
	for i := 0; i < times; i++ {
		now = time.Now()
		ctSigmaA = evaluatorX4Sigma.LinearTransform4ArithmeticSeqNew(ctA, SigmaLT)
		elapsed1 := time.Since(now)
		fmt.Printf("iter %dth: SigmaLT consumes %s \n", i, elapsed1)
		if i > 0 {
			sumtime += elapsed1
		}
	}
	fmt.Printf("SigmaLT consumes %s in average\n", sumtime/time.Duration(times-1))

	sumtime = 0
	for i := 0; i < times; i++ {
		now = time.Now()
		ctDTauA = evaluatorX4Tau.LinearTransform4ArithmeticSeqNew(ctA, DTauLTs[0])
		elapsed1 := time.Since(now)
		fmt.Printf("iter %dth: DTauLT1 consumes %s \n", i, elapsed1)
		now = time.Now()
		ctDTauA = evaluatorX4Tau.LinearTransform4ArithmeticSeqNew(ctDTauA[0], DTauLTs[1])
		elpased2 := time.Since(now)
		fmt.Printf("iter %dth: DTauLT2 consumes %s \n", i, elpased2)
		if i > 0 {
			sumtime += (elapsed1 + elpased2)
		}
	}
	fmt.Printf("DTauLT consumes %s in average\n", sumtime/time.Duration(times-1))

	sumtime = 0
	for i := 0; i < times; i++ {
		now = time.Now()
		ctTauA = evaluatorX4Tau.LinearTransform4ArithmeticSeqNew(ctA, TauLT)
		elapsed1 := time.Since(now)
		fmt.Printf("iter %dth: TauLT consumes %s \n", i, elapsed1)
		if i > 0 {
			sumtime += elapsed1
		}
	}
	fmt.Printf("TauLT consumes %s in average\n", sumtime/time.Duration(times-1))

	fmt.Printf("----------------------------- Multi-Thread Md -----------------------------\n")

	var wg sync.WaitGroup
	var threadNum = 5

	sumtime = 0
	for i := 0; i < times; i++ {
		now = time.Now()
		for j := 0; j < threadNum; j++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
				evaluator_subroute := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: SigmaRotKeys, Rlk: rlk})
				evaluator_subroute.LinearTransform4ArithmeticSeqNew(ctA, DSigmaLTs[0])
				evaluator_subroute.LinearTransform4ArithmeticSeqNew(ctA, DSigmaLTs[1])
			}()
		}
		wg.Wait()
		if i > 0 {
			sumtime += (time.Since(now) / time.Duration(threadNum))
		}
	}
	fmt.Printf("DSigmaLT consumes %s in average in %d-threads mod\n", sumtime/time.Duration(times-1), threadNum)

	sumtime = 0
	for i := 0; i < times; i++ {
		now = time.Now()
		for j := 0; j < threadNum; j++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
				evaluator_subroute := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: SigmaRotKeys, Rlk: rlk})
				evaluator_subroute.LinearTransform4ArithmeticSeqNew(ctA, SigmaLT)
			}()
		}
		wg.Wait()
		if i > 0 {
			sumtime += (time.Since(now) / time.Duration(threadNum))
		}
	}
	fmt.Printf("SigmaLT consumes %s in average in %d-threads mod\n", sumtime/time.Duration(times-1), threadNum)

	sumtime = 0
	for i := 0; i < times; i++ {
		now = time.Now()
		for j := 0; j < threadNum; j++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
				evaluator_subroute := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: TauRotKeys, Rlk: rlk})
				evaluator_subroute.LinearTransform4ArithmeticSeqNew(ctA, DTauLTs[0])
				evaluator_subroute.LinearTransform4ArithmeticSeqNew(ctA, DTauLTs[1])
			}()
		}
		wg.Wait()
		if i > 0 {
			sumtime += (time.Since(now) / time.Duration(threadNum))
		}
	}
	fmt.Printf("DTauLT consumes %s in average in %d-threads mod\n", sumtime/time.Duration(times-1), threadNum)

	sumtime = 0
	for i := 0; i < times; i++ {
		now = time.Now()
		for j := 0; j < threadNum; j++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
				evaluator_subroute := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: TauRotKeys, Rlk: rlk})
				evaluator_subroute.LinearTransform4ArithmeticSeqNew(ctA, TauLT)
			}()
		}
		wg.Wait()
		if i > 0 {
			sumtime += (time.Since(now) / time.Duration(threadNum))
		}
	}
	fmt.Printf("TauLT consumes %s in average in %d-threads mod\n", sumtime/time.Duration(times-1), threadNum)

	// ctTransA := evaluatorX4Trans.LinearTransform4ArithmeticSeqNew(ctA, TransLT)

	// Check the result:
	// auxio.Quick_check_matrix(params, sk, ctA, d, d)
	auxio.Quick_check_matrix(params, sk, ctSigmaA[0], d, d)
	auxio.Quick_check_matrix(params, sk, ctDSigmaA[0], d, d)
	auxio.Quick_check_matrix(params, sk, ctTauA[0], d, d)
	auxio.Quick_check_matrix(params, sk, ctDTauA[0], d, d)
	// auxio.Quick_check_matrix(params, sk, ctTransA[0], d, d)

	print(pk)
	print(rlk)
	print(galk)
	print(a)
	print(SigmaDiagonalMap)
	print(Sigma_BSGSN1)
	print(ctDSigmaA)
	print(ctDTauA)
	print(DSigmaRotKeys)
	print(DTauRotKeys)
	// print(ctTransA)
	print(ctSigmaA)
	print(ctTauA)

}

func LinTransX2_check() {
	var err error
	//	var now time.Time

	d := 128
	Sigma_BSGSN1 := 16
	Tau_BSGSN1 := 16
	DSigma_BSGSN1 := 16
	DTau_BSGSN1 := 16
	Trans_BSGSN1 := 16
	SigmaMaxDiagNo := d / 4
	TauMaxDiagNO := d / 4 * d

	var sk *rlwe.SecretKey
	var pk *rlwe.PublicKey
	var rlk *rlwe.RelinearizationKey
	var galk *rlwe.RotationKeySet

	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: 15,
		Q: []uint64{0x4000000120001, 0x10000140001, 0xffffe80001, // 50 + 17 x 40 =
			0x10000290001, 0xffffc40001, 0x100003e0001,
			0x10000470001, 0x100004b0001, 0xffffb20001,
			0x10000500001, 0x10000650001, // 0xffff940001,
			// 0xffff8a0001, 0xffff820001, 0xffff780001,
			// 0x10000890001, 0xffff750001, 0x10000960001,
		},
		P:            []uint64{0x40000001b0001, 0x3ffffffdf0001, 0x4000000270001}, // 50, 50, 50
		LogSlots:     int(math.Log2(float64(d * d))),
		DefaultScale: 1 << 40,
	})
	if err != nil {
		panic(err)
	}
	Maxlevel := params.MaxLevel()

	/*
		params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
			LogN:  16,
			Sigma: rlwe.DefaultSigma,
			H:     32768,
			Q: []uint64{
				0x10000000006e0001, // 60 Q0
				0x2000000a0001,     // 45
				0x2000000e0001,     // 45
				0x1fffffc20001,     // 45
				0x200000440001,     // 45
				0x200000500001,     // 45
				0x200000620001,     // 45
				0x1fffff980001,     // 45
				0x2000006a0001,     // 45
				0x1fffff7e0001,     // 45
				0x3ffffe80001,      // 42 StC
				0x3ffffd20001,      // 42 StC
				0x3ffffca0001,      // 42 StC
				0xffffffffffc0001,  // 60 ArcSine
				0xfffffffff240001,  // 60 ArcSine
				0x1000000000f00001, // 60 ArcSine
				0xfffffffff840001,  // 60 Double angle
				0x1000000000860001, // 60 Double angle
				0xfffffffff6a0001,  // 60 Sine
				0x1000000000980001, // 60 Sine
				0xfffffffff5a0001,  // 60 Sine
				0x1000000000b00001, // 60 Sine
				0x1000000000ce0001, // 60 Sine
				0xfffffffff2a0001,  // 60 Sine
				0x400000000360001,  // 58 CtS
				0x3ffffffffbe0001,  // 58 CtS
				0x400000000660001,  // 58 CtS
				0x4000000008a0001,  // 58 CtS
			},
			P: []uint64{
				0x1fffffffffe00001, // Pi 61
				0x1fffffffffc80001, // Pi 61
				0x1fffffffffb40001, // Pi 61
				0x1fffffffff500001, // Pi 61
				0x1fffffffff420001, // Pi 61
			},
			LogSlots:     14,
			DefaultScale: 1 << 45,
		})
		if err != nil {
			panic(err)
		}
		Maxlevel := 9
	*/

	kgen := ckks.NewKeyGenerator(params)
	sk, pk = kgen.GenKeyPair()
	rlk = kgen.GenRelinearizationKey(sk, 1)
	// create a Matrix A
	MatrixA := make([][]float64, d)
	for i := 0; i < d; i++ {
		MatrixA[i] = make([]float64, d)
		for j := 0; j < d; j++ {
			MatrixA[i][j] = float64(i*d + j)
		}
	}
	// auxio.Print_matrix_f64_2d(MatrixA, d, d)
	// encode matrixA into a vector using row ordering
	var a []float64
	a, err = mtrxmult.Row_orderingInv(MatrixA)
	if err != nil {
		panic(err)
	}

	// encode into plaintext and encrypt :
	encoder := ckks.NewEncoder(params)
	encryptor := ckks.NewEncryptor(params, pk)
	ptA := encoder.EncodeNew(a, Maxlevel, params.DefaultScale(), params.LogSlots())
	ctA := encryptor.EncryptNew(ptA)
	fmt.Printf("Single ciphertext has: %d bytes", ctA.MarshalBinarySize())

	// create the Map of Sigma LinearTransformation
	var SigmaDiagonalMap map[int][]float64
	SigmaDiagonalMap, err = mtrxmult.Gen_sigma_diagonalVecotrs(d)
	if err != nil {
		panic(err)
	}
	// create the Map of Tau LinearTransformation
	var TauDiagonalMap map[int][]float64
	TauDiagonalMap, err = mtrxmult.Gen_tao_diagonalVectors(d)
	if err != nil {
		panic(err)
	}

	// create the Maps of Decomposed Sigma LinearTransformation
	DSigmaDiagonalMaps1, DSigmaDiagonalMaps2, err := mtrxmult.GenSigmaDiagnalDecomposeMatrices(d, SigmaMaxDiagNo)
	if err != nil {
		panic(err)
	}

	// create the Maps of Decomposed Tau LinearTransformation
	DTauDiagonalMaps, err := mtrxmult.GenTauDiagonalDecomposeMatrices(d, TauMaxDiagNO)
	if err != nil {
		panic(err)
	}

	// generate LinearTransform
	SigmaLT := ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, SigmaDiagonalMap, Maxlevel, params.DefaultScale(), Sigma_BSGSN1, 1, params.LogSlots())
	TauLT := ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, TauDiagonalMap, Maxlevel, params.DefaultScale(), Tau_BSGSN1, d, params.LogSlots())
	DSigmaLTs1 := make([]ckks.LinearTransform, len(DSigmaDiagonalMaps1))
	DSigmaLTs2 := make([]ckks.LinearTransform, len(DSigmaDiagonalMaps2))
	DTauLTs := make([]ckks.LinearTransform, len(DTauDiagonalMaps))
	for i := range DSigmaLTs1 {
		if i == len(DSigmaLTs1)-1 {
			DSigmaLTs1[i] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DSigmaDiagonalMaps1[i], Maxlevel, params.DefaultScale(), DSigma_BSGSN1, 1, params.LogSlots())
		} else {
			DSigmaLTs1[i] = ckks.GenLinearTransform(encoder, DSigmaDiagonalMaps1[i], Maxlevel, params.DefaultScale(), params.LogSlots())
		}
	}
	for i := range DSigmaLTs2 {
		if i == len(DSigmaLTs2)-1 {
			DSigmaLTs2[i] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DSigmaDiagonalMaps2[i], Maxlevel, params.DefaultScale(), DSigma_BSGSN1, 1, params.LogSlots())
		} else {
			DSigmaLTs2[i] = ckks.GenLinearTransform(encoder, DSigmaDiagonalMaps2[i], Maxlevel, params.DefaultScale(), params.LogSlots())
		}
	}
	for i := range DTauLTs {
		if i == len(DTauLTs)-1 {
			DTauLTs[i] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DTauDiagonalMaps[i], Maxlevel, params.DefaultScale(), DTau_BSGSN1, d, params.LogSlots())
		} else {
			DTauLTs[i] = ckks.GenLinearTransform(encoder, DTauDiagonalMaps[i], Maxlevel, params.DefaultScale(), params.LogSlots())
		}
	}

	fmt.Printf("Test configuration: level:%d, d: %d, Sigma'sBSGSN1: %d, Tau'sBSGSN1: %d, DSigma'sBSGSN1:%d, DTau'sBSGSN1: %d, Tran'sBSGSN1: %d, Target SigmaMaxDiagNo: %d, Target TauMaxDiagNo: %d \n", params.MaxLevel(), d, Sigma_BSGSN1, Tau_BSGSN1, DSigma_BSGSN1, DTau_BSGSN1, Trans_BSGSN1, SigmaMaxDiagNo, TauMaxDiagNO)

	// generate Rotations
	fmt.Printf("\nSigmaLT %d Rotations: ", len(SigmaLT.Rotations4ArithmeticSeq()))
	for _, i := range SigmaLT.Rotations4ArithmeticSeq() {
		fmt.Printf("%d ", i)
	}

	var DSigmaRotset = make(map[int]bool)
	for k := range DSigmaLTs1 {
		for _, i := range DSigmaLTs1[k].Rotations4ArithmeticSeq() {
			if !DSigmaRotset[i] {
				DSigmaRotset[i] = true
			}
		}
	}
	for k := range DSigmaLTs2 {
		for _, i := range DSigmaLTs2[k].Rotations4ArithmeticSeq() {
			if !DSigmaRotset[i] {
				DSigmaRotset[i] = true
			}
		}
	}
	fmt.Printf("\nDSigmaLT %d Rotations: ", len(DSigmaRotset))
	for rot := range DSigmaRotset {
		fmt.Printf("%d ", rot)
	}

	var DTauRotset = make(map[int]bool)
	fmt.Printf("\nTauLT %d Rotations : ", len(TauLT.Rotations4ArithmeticSeq()))
	for _, i := range TauLT.Rotations4ArithmeticSeq() {
		fmt.Printf("%d ", i)
	}

	for k := range DTauLTs {
		for _, i := range DTauLTs[k].Rotations4ArithmeticSeq() {
			if !DTauRotset[i] {
				DTauRotset[i] = true
			}
		}
	}
	fmt.Printf("\nDTauLT %d Rotations: ", len(DTauRotset))
	for rot := range DTauRotset {
		fmt.Printf("%d ", rot)
	}
	//-------------------------------------------------- MtrxMult Test 1------------------------------------------------------
	fmt.Printf("MtrxMult Ver2 test\n")
	MtrxMultRotset := make([]int, 0)
	for rot := range DSigmaRotset {
		MtrxMultRotset = append(MtrxMultRotset, rot)
	}
	for rot := range DTauRotset {
		MtrxMultRotset = append(MtrxMultRotset, rot)
	}
	MtrxMultRotset = append(MtrxMultRotset, -d+d*d)
	MtrxMultKeys := kgen.GenRotationKeysForRotations(MtrxMultRotset, false, sk)

	// Gen ColShift LTs:
	ColShiftLTs := make([]ckks.LinearTransform, 0)
	for k := 1; k < d; k++ {
		U, err := mtrxmult.Gen_colShift_diagonalVectors(d, k)
		if err != nil {
			panic(err)
		}
		LT := ckks.GenLinearTransform(encoder, U, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		ColShiftLTs = append(ColShiftLTs, LT)
	}

	ctBs := make([]*rlwe.Ciphertext, 1)
	for i := range ctBs {
		ctBs[i] = ctA.CopyNew()
	}
	/*
		ctRslt, err := mtrxmult.SquareMatrix_MultBSGS_DecomposeVer2_NeedLTs_dbg(params, sk, rlk, MtrxMultKeys, ctA, ctBs, DSigmaLTs1, DSigmaLTs2, DTauLTs, ColShiftLTs, d, 1)
		if err != nil {
			panic(err)
		}
		for i := range ctRslt {
			auxio.Quick_check_matrix(params, sk, ctRslt[i], d, d)
		}
	*/

	//-------------------------------------------------- MtrxMult Test 2------------------------------------------------------
	fmt.Printf("MtrxMtul Ver4 test\n")
	/*
		MtrxMultRotset := make([]int, 0)
		for rot := range DSigmaRotset {
			MtrxMultRotset = append(MtrxMultRotset, rot)
		}
		for rot := range DTauRotset {
			MtrxMultRotset = append(MtrxMultRotset, rot)
		}
		MtrxMultRotset = append(MtrxMultRotset, -d+d*d)
		MtrxMultKeys := kgen.GenRotationKeysForRotations(MtrxMultRotset, false, sk)

		// Gen ColShift LTs:
		ColShiftLTs := make([]ckks.LinearTransform, 0)
		for k := 1; k < d; k++ {
			U, err := mtrxmult.Gen_colShift_diagonalVectors(d, k)
			if err != nil {
				panic(err)
			}
			LT := ckks.GenLinearTransform(encoder, U, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
			ColShiftLTs = append(ColShiftLTs, LT)
		}

		ctBs := make([]*rlwe.Ciphertext, 3)
		for i := range ctBs {
			ctBs[i] = ctA.CopyNew()
		}
	*/

	ctAs_dmp, _ := mtrxmult.Sigma_linearTransform_MultiThread_dbg(params, sk, rlk, MtrxMultKeys, ctBs, DSigmaLTs1, DSigmaLTs2, ColShiftLTs, d, 1)
	ctBs_dmp, _ := mtrxmult.Tau_linearTransform_MultiThread_dbg(params, sk, rlk, MtrxMultKeys, ctBs, DTauLTs, d, 1)
	ctCs, _ := mtrxmult.SquareMatrix_MultBSGS_DecomposeVer4_NeedLTs_dbg(params, sk, rlk, ctAs_dmp[0], ctBs_dmp, d, 1)
	for _, ct := range ctCs {
		auxio.Quick_check_matrix(params, sk, ct, d, d)
	}
	/*
		if 1 == 1 {
			return
		}
	*/

	//-------------------------------------------------- LinTrans Test ------------------------------------------------------

	// generate Keys.
	now := time.Now()
	SigmaRotKeys := kgen.GenRotationKeysForRotations(SigmaLT.Rotations4ArithmeticSeq(), false, sk)
	fmt.Printf("\nSigma Rot keys generation consumes %s and %d MB\n", time.Since(now), SigmaRotKeys.MarshalBinarySize()/(1024*1024))
	now = time.Now()
	TauRotKeys := kgen.GenRotationKeysForRotations(TauLT.Rotations4ArithmeticSeq(), false, sk)
	fmt.Printf("\nTau Rot keys generation consumes %s and %d MB\n", time.Since(now), TauRotKeys.MarshalBinarySize()/(1024*1024))

	DSigmaRots := make([]int, 0)
	for rot := range DSigmaRotset {
		DSigmaRots = append(DSigmaRots, rot)
	}
	now = time.Now()
	DSigmaRotKeys := kgen.GenRotationKeysForRotations(DSigmaRots, false, sk)
	fmt.Printf("\nDSigma Rot keys generation consumes %s and %d MB\n", time.Since(now), DSigmaRotKeys.MarshalBinarySize()/(1024*1024))

	DTauRots := make([]int, 0)
	for rot := range DTauRotset {
		DTauRots = append(DTauRots, rot)
	}
	now = time.Now()
	DTauRotKeys := kgen.GenRotationKeysForRotations(DTauRots, false, sk)
	fmt.Printf("\nDTau Rot keys generation consumes %s and %d MB\n", time.Since(now), DTauRotKeys.MarshalBinarySize()/(1024*1024))

	// create evaluatorX
	evaluatorX4Sigma := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: SigmaRotKeys, Rlk: rlk})
	evaluatorX4Tau := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: TauRotKeys, Rlk: rlk})

	times := 5
	var sumtime time.Duration
	var ctDSigmaA, ctDTauA []*rlwe.Ciphertext
	var ctSigmaA, ctTauA []*rlwe.Ciphertext

	fmt.Printf("----------------------------- Single Thread Mod -----------------------------\n")
	sumtime = 0
	LTstemp := make([]ckks.LinearTransform, 2)
	LTstemp[0] = DSigmaLTs1[len(DSigmaLTs1)-1]
	LTstemp[1] = DSigmaLTs2[len(DSigmaLTs2)-1]
	for i := 0; i < times; i++ {
		now = time.Now()
		ctDSigmaA = evaluatorX4Sigma.LinearTransform4ArithmeticSeqNew(ctA, LTstemp, true)
		//auxio.Quick_check_matrix(params, sk, ctDSigmaA[0], d, d)
		// auxio.Quick_check_matrix(params, sk, ctDSigmaA[1], d, d)
		elapsed := time.Since(now)
		for i := len(DSigmaLTs1) - 2; i >= 0; i-- {
			now = time.Now()
			ctDSigmaA[0] = evaluatorX4Sigma.LinearTransform4ArithmeticSeqNew(ctDSigmaA[0], DSigmaLTs1[i])[0]
			evaluatorX4Sigma.Rescale(ctDSigmaA[0], params.DefaultScale(), ctDSigmaA[0])
			//auxio.Quick_check_matrix(params, sk, ctDSigmaA[0], d, d)
			elapsed += time.Since(now)
		}
		for i := len(DSigmaLTs2) - 2; i >= 0; i-- {
			now = time.Now()
			ctDSigmaA[1] = evaluatorX4Sigma.LinearTransform4ArithmeticSeqNew(ctDSigmaA[1], DSigmaLTs2[i])[0]
			evaluatorX4Sigma.Rescale(ctDSigmaA[1], params.DefaultScale(), ctDSigmaA[1])
			//auxio.Quick_check_matrix(params, sk, ctDSigmaA[1], d, d)
			elapsed += time.Since(now)
		}
		now = time.Now()
		evaluatorX4Sigma.SetScale(ctDSigmaA[1], ctDSigmaA[0].Scale)
		evaluatorX4Sigma.Add(ctDSigmaA[0], ctDSigmaA[1], ctDSigmaA[0])
		// auxio.Quick_check_matrix(params, sk, ctDSigmaA[0], d, d)
		elapsed += time.Since(now)
		if i > 0 {
			sumtime += elapsed
		}
	}
	fmt.Printf("DSigmaLT consumes %s in average\n", sumtime/time.Duration(times-1))

	sumtime = 0
	for i := 0; i < times; i++ {
		now = time.Now()
		ctSigmaA = evaluatorX4Sigma.LinearTransform4ArithmeticSeqNew(ctA, SigmaLT)
		elapsed1 := time.Since(now)
		fmt.Printf("iter %dth: SigmaLT consumes %s \n", i, elapsed1)
		if i > 0 {
			sumtime += elapsed1
		}
	}
	fmt.Printf("SigmaLT consumes %s in average\n", sumtime/time.Duration(times-1))

	sumtime = 0
	for i := 0; i < times; i++ {
		now = time.Now()
		ctDTauA = evaluatorX4Tau.LinearTransform4ArithmeticSeqNew(ctA, DTauLTs[len(DTauLTs)-1])
		elapsed := time.Since(now)
		for i := len(DTauLTs) - 2; i >= 0; i-- {
			now = time.Now()
			ctDTauA[0] = evaluatorX4Tau.LinearTransform4ArithmeticSeqNew(ctDTauA[0], DTauLTs[i])[0]
			elapsed += time.Since(now)
		}
		if i > 0 {
			sumtime += elapsed
		}
	}
	fmt.Printf("DTauLT consumes %s in average\n", sumtime/time.Duration(times-1))

	sumtime = 0
	for i := 0; i < times; i++ {
		now = time.Now()
		ctTauA = evaluatorX4Tau.LinearTransform4ArithmeticSeqNew(ctA, TauLT)
		elapsed1 := time.Since(now)
		fmt.Printf("iter %dth: TauLT consumes %s \n", i, elapsed1)
		if i > 0 {
			sumtime += elapsed1
		}
	}
	fmt.Printf("TauLT consumes %s in average\n", sumtime/time.Duration(times-1))

	// Check the result:
	// auxio.Quick_check_matrix(params, sk, ctA, d, d)
	auxio.Quick_check_infos(ctSigmaA[0], "ctSigmaA")
	auxio.Quick_check_matrix(params, sk, ctSigmaA[0], d, d)
	auxio.Quick_check_infos(ctDSigmaA[0], "ctDSigmaA")
	auxio.Quick_check_matrix(params, sk, ctDSigmaA[0], d, d)
	auxio.Quick_check_infos(ctTauA[0], "ctTauA")
	auxio.Quick_check_matrix(params, sk, ctTauA[0], d, d)
	auxio.Quick_check_infos(ctDTauA[0], "ctDTauA")
	auxio.Quick_check_matrix(params, sk, ctDTauA[0], d, d)

	print(pk)
	print(rlk)
	print(galk)
	print(a)
	print(SigmaDiagonalMap)
	print(Sigma_BSGSN1)
	print(ctDSigmaA)
	print(ctDTauA)
	print(DSigmaRotKeys)
	print(DTauRotKeys)
	// print(ctTransA)
	print(ctSigmaA)
	print(ctTauA)

}

func LinTransX3_check() {
	var err error
	//	var now time.Time

	d := 128
	dmpFactor := 4
	Sigma_BSGSN1 := 8
	Tau_BSGSN1 := 8
	DSigma_BSGSN1 := 8
	DTau_BSGSN1 := 8
	Trans_BSGSN1 := 8
	SigmaMaxDiagNo := d / dmpFactor
	TauMaxDiagNO := d / dmpFactor * d

	var sk *rlwe.SecretKey
	var pk *rlwe.PublicKey
	var rlk *rlwe.RelinearizationKey
	var galk *rlwe.RotationKeySet

	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: 15,
		Q: []uint64{0x4000000120001, 0x10000140001, 0xffffe80001, // 50 + 17 x 40 =
			0x10000290001, 0xffffc40001, 0x100003e0001,
			0x10000470001, 0x100004b0001, 0xffffb20001,
			0x10000500001, 0x10000650001, 0xffff940001,
			0xffff8a0001, 0xffff820001, 0xffff780001,
			// 0x10000890001, 0xffff750001, 0x10000960001,
		},
		P:            []uint64{0x40000001b0001, 0x3ffffffdf0001, 0x4000000270001}, // 50, 50, 50
		LogSlots:     int(math.Log2(float64(d * d))),
		DefaultScale: 1 << 40,
	})
	if err != nil {
		panic(err)
	}
	Minlevel := params.MaxLevel() - 4 // dmpFacotr should not greater than 4.
	Maxlevel := Minlevel + dmpFactor

	/*
		params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
			LogN:  16,
			Sigma: rlwe.DefaultSigma,
			H:     32768,
			Q: []uint64{
				0x10000000006e0001, // 60 Q0
				0x2000000a0001,     // 45
				0x2000000e0001,     // 45
				0x1fffffc20001,     // 45
				0x200000440001,     // 45
				0x200000500001,     // 45
				0x200000620001,     // 45
				0x1fffff980001,     // 45
				0x2000006a0001,     // 45
				0x1fffff7e0001,     // 45
				0x3ffffe80001,      // 42 StC
				0x3ffffd20001,      // 42 StC
				0x3ffffca0001,      // 42 StC
				0xffffffffffc0001,  // 60 ArcSine
				0xfffffffff240001,  // 60 ArcSine
				0x1000000000f00001, // 60 ArcSine
				0xfffffffff840001,  // 60 Double angle
				0x1000000000860001, // 60 Double angle
				0xfffffffff6a0001,  // 60 Sine
				0x1000000000980001, // 60 Sine
				0xfffffffff5a0001,  // 60 Sine
				0x1000000000b00001, // 60 Sine
				0x1000000000ce0001, // 60 Sine
				0xfffffffff2a0001,  // 60 Sine
				0x400000000360001,  // 58 CtS
				0x3ffffffffbe0001,  // 58 CtS
				0x400000000660001,  // 58 CtS
				0x4000000008a0001,  // 58 CtS
			},
			P: []uint64{
				0x1fffffffffe00001, // Pi 61
				0x1fffffffffc80001, // Pi 61
				0x1fffffffffb40001, // Pi 61
				0x1fffffffff500001, // Pi 61
				0x1fffffffff420001, // Pi 61
			},
			LogSlots:     14,
			DefaultScale: 1 << 45,
		})
		if err != nil {
			panic(err)
		}
		Maxlevel := 9
	*/

	kgen := ckks.NewKeyGenerator(params)
	sk, pk = kgen.GenKeyPair()
	rlk = kgen.GenRelinearizationKey(sk, 1)
	// create a Matrix A
	MatrixA := make([][]float64, d)
	for i := 0; i < d; i++ {
		MatrixA[i] = make([]float64, d)
		for j := 0; j < d; j++ {
			MatrixA[i][j] = float64(i*d + j)
		}
	}
	// auxio.Print_matrix_f64_2d(MatrixA, d, d)
	// encode matrixA into a vector using row ordering
	var a []float64
	a, err = mtrxmult.Row_orderingInv(MatrixA)
	if err != nil {
		panic(err)
	}

	// encode into plaintext and encrypt :
	encoder := ckks.NewEncoder(params)
	encryptor := ckks.NewEncryptor(params, pk)
	ptA := encoder.EncodeNew(a, Maxlevel, params.DefaultScale(), params.LogSlots())
	ctA := encryptor.EncryptNew(ptA)
	fmt.Printf("Single ciphertext has: %d bytes", ctA.MarshalBinarySize())

	// create the Map of Sigma LinearTransformation
	var SigmaDiagonalMap map[int][]float64
	SigmaDiagonalMap, err = mtrxmult.Gen_sigma_diagonalVecotrs(d)
	if err != nil {
		panic(err)
	}
	// create the Map of Tau LinearTransformation
	var TauDiagonalMap map[int][]float64
	TauDiagonalMap, err = mtrxmult.Gen_tao_diagonalVectors(d)
	if err != nil {
		panic(err)
	}

	// create the Maps of Decomposed Sigma LinearTransformation
	DSigmaDiagonalMaps1, DSigmaDiagonalMaps2, err := mtrxmult.GenSigmaDiagnalDecomposeMatrices_Ver2(d, SigmaMaxDiagNo)
	if err != nil {
		panic(err)
	}

	// create the Maps of Decomposed Tau LinearTransformation
	DTauDiagonalMaps, err := mtrxmult.GenTauDiagonalDecomposeMatrices_Ver2(d, TauMaxDiagNO)
	if err != nil {
		panic(err)
	}

	// generate LinearTransform
	SigmaLT := ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, SigmaDiagonalMap, Maxlevel, params.DefaultScale(), Sigma_BSGSN1, 1, params.LogSlots())
	TauLT := ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, TauDiagonalMap, Maxlevel, params.DefaultScale(), Tau_BSGSN1, d, params.LogSlots())
	DSigmaLTs1 := make([]ckks.LinearTransform, len(DSigmaDiagonalMaps1))
	DSigmaLTs2 := make([]ckks.LinearTransform, len(DSigmaDiagonalMaps2))
	DTauLTs := make([]ckks.LinearTransform, len(DTauDiagonalMaps))
	for i := range DSigmaLTs1 {
		if i == 0 {
			DSigmaLTs1[i] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DSigmaDiagonalMaps1[i], Maxlevel, params.DefaultScale(), DSigma_BSGSN1, 1, params.LogSlots())
		} else {
			DSigmaLTs1[i] = ckks.GenLinearTransform(encoder, DSigmaDiagonalMaps1[i], Maxlevel, params.DefaultScale(), params.LogSlots())
		}
	}
	for i := range DSigmaLTs2 {
		if i == 0 {
			DSigmaLTs2[i] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DSigmaDiagonalMaps2[i], Maxlevel, params.DefaultScale(), DSigma_BSGSN1, 1, params.LogSlots())
		} else {
			DSigmaLTs2[i] = ckks.GenLinearTransform(encoder, DSigmaDiagonalMaps2[i], Maxlevel, params.DefaultScale(), params.LogSlots())
		}
	}
	for i := range DTauLTs {
		if i == 0 {
			DTauLTs[i] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DTauDiagonalMaps[i], Maxlevel, params.DefaultScale(), DTau_BSGSN1, d, params.LogSlots())
		} else {
			DTauLTs[i] = ckks.GenLinearTransform(encoder, DTauDiagonalMaps[i], Maxlevel, params.DefaultScale(), params.LogSlots())
		}
	}

	fmt.Printf("Test configuration: level:%d, d: %d, Sigma'sBSGSN1: %d, Tau'sBSGSN1: %d, DSigma'sBSGSN1:%d, DTau'sBSGSN1: %d, Tran'sBSGSN1: %d, Target SigmaMaxDiagNo: %d, Target TauMaxDiagNo: %d \n", params.MaxLevel(), d, Sigma_BSGSN1, Tau_BSGSN1, DSigma_BSGSN1, DTau_BSGSN1, Trans_BSGSN1, SigmaMaxDiagNo, TauMaxDiagNO)

	// generate Rotations
	fmt.Printf("\nSigmaLT %d Rotations: ", len(SigmaLT.Rotations4ArithmeticSeq()))
	for _, i := range SigmaLT.Rotations4ArithmeticSeq() {
		fmt.Printf("%d ", i)
	}

	var DSigmaRotset = make(map[int]bool)
	for k := range DSigmaLTs1 {
		for _, i := range DSigmaLTs1[k].Rotations4ArithmeticSeq() {
			if !DSigmaRotset[i] {
				DSigmaRotset[i] = true
			}
		}
	}
	for k := range DSigmaLTs2 {
		for _, i := range DSigmaLTs2[k].Rotations4ArithmeticSeq() {
			if !DSigmaRotset[i] {
				DSigmaRotset[i] = true
			}
		}
	}
	fmt.Printf("\nDSigmaLT %d Rotations: ", len(DSigmaRotset))
	for rot := range DSigmaRotset {
		fmt.Printf("%d ", rot)
	}

	var DTauRotset = make(map[int]bool)
	fmt.Printf("\nTauLT %d Rotations : ", len(TauLT.Rotations4ArithmeticSeq()))
	for _, i := range TauLT.Rotations4ArithmeticSeq() {
		fmt.Printf("%d ", i)
	}

	for k := range DTauLTs {
		for _, i := range DTauLTs[k].Rotations4ArithmeticSeq() {
			if !DTauRotset[i] {
				DTauRotset[i] = true
			}
		}
	}
	fmt.Printf("\nDTauLT %d Rotations: ", len(DTauRotset))
	for rot := range DTauRotset {
		fmt.Printf("%d ", rot)
	}

	//-------------------------------------------------- LinTrans Test ------------------------------------------------------

	// generate Keys.
	now := time.Now()
	SigmaRotKeys := kgen.GenRotationKeysForRotations(SigmaLT.Rotations4ArithmeticSeq(), false, sk)
	fmt.Printf("\nSigma Rot keys generation consumes %s and %d MB\n", time.Since(now), SigmaRotKeys.MarshalBinarySize()/(1024*1024))
	now = time.Now()
	TauRotKeys := kgen.GenRotationKeysForRotations(TauLT.Rotations4ArithmeticSeq(), false, sk)
	fmt.Printf("\nTau Rot keys generation consumes %s and %d MB\n", time.Since(now), TauRotKeys.MarshalBinarySize()/(1024*1024))

	DSigmaRots := make([]int, 0)
	for rot := range DSigmaRotset {
		DSigmaRots = append(DSigmaRots, rot)
	}
	now = time.Now()
	DSigmaRotKeys := kgen.GenRotationKeysForRotations(DSigmaRots, false, sk)
	fmt.Printf("\nDSigma Rot keys generation consumes %s and %d MB\n", time.Since(now), DSigmaRotKeys.MarshalBinarySize()/(1024*1024))

	DTauRots := make([]int, 0)
	for rot := range DTauRotset {
		DTauRots = append(DTauRots, rot)
	}
	now = time.Now()
	DTauRotKeys := kgen.GenRotationKeysForRotations(DTauRots, false, sk)
	fmt.Printf("\nDTau Rot keys generation consumes %s and %d MB\n", time.Since(now), DTauRotKeys.MarshalBinarySize()/(1024*1024))

	// create evaluatorX
	evaluatorX4Sigma := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: SigmaRotKeys, Rlk: rlk})
	evaluatorX4Tau := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rtks: TauRotKeys, Rlk: rlk})

	times := 5
	var sumtime time.Duration
	var ctDSigmaA, ctDTauA []*rlwe.Ciphertext
	var ctSigmaA, ctTauA []*rlwe.Ciphertext

	fmt.Printf("ctA in level %d has %d bytes\n", ctA.Level(), ctA.MarshalBinarySize())
	fmt.Printf("ctA in level %d has %d bytes\n", ctA.Level()-(dmpFactor-1), evaluatorX4Sigma.DropLevelNew(ctA, (dmpFactor-1)).MarshalBinarySize())

	fmt.Printf("----------------------------- Single Thread Mod -----------------------------\n")
	sumtime = 0
	LTstemp := make([]ckks.LinearTransform, 2)
	LTstemp[0] = DSigmaLTs1[len(DSigmaLTs1)-1]
	LTstemp[1] = DSigmaLTs2[len(DSigmaLTs2)-1]
	for i := 0; i < times; i++ {
		now = time.Now()
		//auxio.Quick_check_matrix(params, sk, ctDSigmaA[0], d, d)
		// auxio.Quick_check_matrix(params, sk, ctDSigmaA[1], d, d)
		ctDSigmaA = evaluatorX4Sigma.LinearTransform4ArithmeticSeqNew(ctA, DSigmaLTs1[len(DSigmaLTs1)-1])
		ctDSigmaA = append(ctDSigmaA, evaluatorX4Sigma.LinearTransform4ArithmeticSeqNew(ctA, DSigmaLTs2[len(DSigmaLTs2)-1])[0])
		elapsed := time.Since(now)
		for i := len(DSigmaLTs1) - 2; i >= 0; i-- {
			now = time.Now()
			evaluatorX4Sigma.Rescale(ctDSigmaA[0], params.DefaultScale(), ctDSigmaA[0])
			ctDSigmaA[0] = evaluatorX4Sigma.LinearTransform4ArithmeticSeqNew(ctDSigmaA[0], DSigmaLTs1[i])[0]
			//auxio.Quick_check_matrix(params, sk, ctDSigmaA[0], d, d)
			elapsed += time.Since(now)
		}
		for i := len(DSigmaLTs2) - 2; i >= 0; i-- {
			now = time.Now()
			evaluatorX4Sigma.Rescale(ctDSigmaA[1], params.DefaultScale(), ctDSigmaA[1])
			ctDSigmaA[1] = evaluatorX4Sigma.LinearTransform4ArithmeticSeqNew(ctDSigmaA[1], DSigmaLTs2[i])[0]
			//auxio.Quick_check_matrix(params, sk, ctDSigmaA[1], d, d)
			elapsed += time.Since(now)
		}
		now = time.Now()
		//evaluatorX4Sigma.SetScale(ctDSigmaA[1], ctDSigmaA[0].Scale)
		evaluatorX4Sigma.Add(ctDSigmaA[0], ctDSigmaA[1], ctDSigmaA[0])
		// auxio.Quick_check_matrix(params, sk, ctDSigmaA[0], d, d)
		elapsed += time.Since(now)
		if i > 0 {
			sumtime += elapsed
		}
	}
	fmt.Printf("DSigmaLT consumes %s in average\n", sumtime/time.Duration(times-1))

	sumtime = 0
	for i := 0; i < times; i++ {
		now = time.Now()
		ctSigmaA = evaluatorX4Sigma.LinearTransform4ArithmeticSeqNew(evaluatorX4Sigma.DropLevelNew(ctA, dmpFactor-1), SigmaLT)
		elapsed1 := time.Since(now)
		fmt.Printf("iter %dth: SigmaLT consumes %s \n", i, elapsed1)
		if i > 0 {
			sumtime += elapsed1
		}
	}
	fmt.Printf("SigmaLT consumes %s in average\n", sumtime/time.Duration(times-1))

	sumtime = 0
	for i := 0; i < times; i++ {
		now = time.Now()
		ctDTauA = evaluatorX4Tau.LinearTransform4ArithmeticSeqNew(ctA, DTauLTs[len(DTauLTs)-1])
		elapsed := time.Since(now)
		for i := len(DTauLTs) - 2; i >= 0; i-- {
			now = time.Now()
			evaluatorX4Tau.Rescale(ctDTauA[0], params.DefaultScale(), ctDTauA[0])
			ctDTauA[0] = evaluatorX4Tau.LinearTransform4ArithmeticSeqNew(ctDTauA[0], DTauLTs[i])[0]
			elapsed += time.Since(now)
		}
		if i > 0 {
			sumtime += elapsed
		}
	}
	fmt.Printf("DTauLT consumes %s in average\n", sumtime/time.Duration(times-1))

	sumtime = 0
	for i := 0; i < times; i++ {
		now = time.Now()
		ctTauA = evaluatorX4Tau.LinearTransform4ArithmeticSeqNew(evaluatorX4Tau.DropLevelNew(ctA, dmpFactor-1), TauLT)
		elapsed1 := time.Since(now)
		fmt.Printf("iter %dth: TauLT consumes %s \n", i, elapsed1)
		if i > 0 {
			sumtime += elapsed1
		}
	}
	fmt.Printf("TauLT consumes %s in average\n", sumtime/time.Duration(times-1))

	// Check the result:
	// auxio.Quick_check_matrix(params, sk, ctA, d, d)
	auxio.Quick_check_infos(ctSigmaA[0], "ctSigmaA")
	auxio.Quick_check_matrix(params, sk, ctSigmaA[0], d, d)
	auxio.Quick_check_infos(ctDSigmaA[0], "ctDSigmaA")
	auxio.Quick_check_matrix(params, sk, ctDSigmaA[0], d, d)
	auxio.Quick_check_infos(ctTauA[0], "ctTauA")
	auxio.Quick_check_matrix(params, sk, ctTauA[0], d, d)
	auxio.Quick_check_infos(ctDTauA[0], "ctDTauA")
	auxio.Quick_check_matrix(params, sk, ctDTauA[0], d, d)

	print(pk)
	print(rlk)
	print(galk)
	print(a)
	print(SigmaDiagonalMap)
	print(Sigma_BSGSN1)
	print(ctDSigmaA)
	print(ctDTauA)
	print(DSigmaRotKeys)
	print(DTauRotKeys)
	// print(ctTransA)
	print(ctSigmaA)
	print(ctTauA)

}

func MatrixMultX_check(loadkeysFromDisk bool) {
	var err error
	//	var now time.Time

	// N1
	N1Sigma := 16
	N1Tau := 16
	N1Trans := 16

	// initialize encryption scheme params.
	var params ckks.Parameters
	var sk *rlwe.SecretKey
	var pk *rlwe.PublicKey
	var rlk *rlwe.RelinearizationKey
	var kgen rlwe.KeyGenerator

	// generate Keys, we can load keys from disk or generate it.
	if loadkeysFromDisk {
		params, sk, pk, rlk, err = LoadKeys_check()
		if err != nil {
			panic(err)
		}
	} else {
		params, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
			LogN: 15,
			Q: []uint64{1152921504598720513,
				1099490000897, 1099498258433,
				1099499175937, 1099499569153,
				1099500617729}, // 60 + 40*5
			P:            []uint64{1152921504606584833}, // 60
			LogSlots:     14,
			DefaultScale: 1 << 40,
		})
		if err != nil {
			panic(err)
		}
		kgen = ckks.NewKeyGenerator(params)
		sk, pk = kgen.GenKeyPair()
		// Relinearization key
		rlk = kgen.GenRelinearizationKey(sk, 1)

		// Print Keys' size:
		fmt.Printf("Encryption Scheme Params occupy %d bytes\n", params.MarshalBinarySize())
		fmt.Printf("Secret Key:%d bytes, Public Key:%d bytes\n", sk.MarshalBinarySize(), pk.MarshalBinarySize())
		fmt.Printf("Relinearization Key: %d bytes\n", rlk.MarshalBinarySize())

		// Store this set of keys, for convenience:
		keysTracker := auxio.NewTracker("Keys", -1)
		n := 0
		n, err = keysTracker.StoreUpdateOne(&params)
		if err != nil {
			panic(err)
		} else {
			fmt.Printf("Store params %d bytes\n", n)
		}
		n, err = keysTracker.StoreUpdateOne(sk)
		if err != nil {
			panic(err)
		} else {
			fmt.Printf("Store sk %d bytes\n", n)
		}
		n, err = keysTracker.StoreUpdateOne(pk)
		if err != nil {
			panic(err)
		} else {
			fmt.Printf("Store pk %d bytes\n", n)
		}
		n, err = keysTracker.StoreUpdateOne(rlk)
		if err != nil {
			panic(err)
		} else {
			fmt.Printf("Store rlk %d bytes\n", n)
		}
		n, err = keysTracker.StoreFinish()
		if err != nil {
			panic(err)
		} else {
			fmt.Printf("Store filetracker %d bytes\n", n)
		}
	}

	d := int(math.Sqrt(float64(int(1) << params.LogSlots())))

	// Encryptor
	encryptor := ckks.NewEncryptor(params, pk)

	// encoder
	encoder := ckks.NewEncoder(params)

	// evaluator
	// evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})

	// create the Map and LinearTransform object of Transpose LinearTransformation
	fmt.Printf("create the Map and LinearTransform object of Transpose LinearTransformation\n")
	var TransposeDiagonalMap map[int][]float64
	TransposeDiagonalMap, err = mtrxmult.Gen_transpose_diagonalVectors(d)
	if err != nil {
		panic(err)
	}
	TransposeLT := ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, TransposeDiagonalMap, params.MaxLevel(), params.DefaultScale(), N1Trans, (d - 1), params.LogSlots())

	// create the Map of Sigma LinearTransformation
	fmt.Printf("create the Map of Sigma LinearTransformation\n")
	var SigmaDiagonalMap map[int][]float64
	SigmaDiagonalMap, err = mtrxmult.Gen_sigma_diagonalVecotrs(d)
	if err != nil {
		panic(err)
	}

	// create the Map of Tau LinearTransformation
	fmt.Printf("create the Map of Tau LinearTransformation\n")
	var TauDiagonalMap map[int][]float64
	TauDiagonalMap, err = mtrxmult.Gen_tao_diagonalVectors(d)
	if err != nil {
		panic(err)
	}

	// create the Decomposed Maps and LinearTransform object of Sigma LinearTransformation
	fmt.Printf("create the Decomposed Maps and LinearTransform object of Sigma LinearTransformation\n")
	SigmaMatrix := mtrxmult.DiagonalVectors2Matrix(SigmaDiagonalMap, d*d)
	DSigmaLTs := make([]ckks.LinearTransform, 2)
	DSigmaDiagonalMaps, err := mtrxmult.Converge2DiagonalDecompose_Sigma(SigmaMatrix)
	if err != nil {
		panic(err)
	}
	DSigmaLTs[0] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DSigmaDiagonalMaps[0], params.MaxLevel(), params.DefaultScale(), N1Sigma, 1, params.LogSlots())
	DSigmaLTs[1] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DSigmaDiagonalMaps[1], params.MaxLevel(), params.DefaultScale(), N1Sigma, 1, params.LogSlots())

	// create the Decomposed Maps of Tau LinearTransformation
	fmt.Printf("create the Decomposed Maps of Tau LinearTransformation\n")
	TauMatrix := mtrxmult.DiagonalVectors2Matrix(TauDiagonalMap, d*d)
	DTauLTs := make([]ckks.LinearTransform, 2)
	DTauDiagonalMaps, err := mtrxmult.Converge2DiagonalDecompose_Tao(TauMatrix)
	if err != nil {
		panic(err)
	}
	DTauLTs[0] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DTauDiagonalMaps[0], params.MaxLevel(), params.DefaultScale(), N1Tau, d, params.LogSlots())
	DTauLTs[1] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DTauDiagonalMaps[1], params.MaxLevel(), params.DefaultScale(), N1Tau, d, params.LogSlots())

	// create the Map of Plaintext for ColShift Diagonals
	ColShiftPTs := make(map[int][]*rlwe.Plaintext)
	var ColShiftMap map[int][]float64
	for k := 1; k < d; k++ {
		ColShiftMap, err = mtrxmult.Gen_colShift_diagonalVectors(d, k)
		if err != nil {
			panic(err)
		}
		ColShiftPTs[k] = make([]*rlwe.Plaintext, 2)
		ColShiftPTs[k][0] = encoder.EncodeNew(ColShiftMap[k], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		ColShiftPTs[k][1] = encoder.EncodeNew(ColShiftMap[k-d], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
	}

	// Create Rotation keys correspond to Transposition, MatrixMult and RowTotalSum
	var galk4Transpose *rlwe.RotationKeySet
	var galk4MtrxMult *rlwe.RotationKeySet
	var galk4RowtotalSum *rlwe.RotationKeySet

	// Create File Tracker to Store or load these keys
	fmt.Printf("create File Tracker to Store or Load the Rotation keys\n")
	var RtkTracker4Transpose *auxio.Filetracker
	var RtkTracker4MtrxMult *auxio.Filetracker
	var RtkTracker4RowtotalSum *auxio.Filetracker

	if loadkeysFromDisk {
		// we are not going to imediately get these keys from disk when they already exists in disk.
		fmt.Printf("we only Init the Tracker for rotationKeySets, but are not going to imediately get these keys from disk.\n")
		RtkTracker4Transpose = auxio.NewTracker4File("Rtk4Transpose")
		RtkTracker4MtrxMult = auxio.NewTracker4File("Rtk4MtrxMult")
		RtkTracker4RowtotalSum = auxio.NewTracker4File("Rtk4RowtotalSum")
	} else {
		// Create File Tracker to Store or load these keys
		fmt.Printf("Init File Tracker to Store or Load the Rotation keys\n")
		RtkTracker4Transpose = auxio.NewTracker("Rtk4Transpose", -1)
		RtkTracker4MtrxMult = auxio.NewTracker("Rtk4MtrxMult", -1)
		RtkTracker4RowtotalSum = auxio.NewTracker("Rtk4RowtotalSum", -1)

		// Create Rotation keys correspond to Transposition, MatrixMult and RowTotalSum
		fmt.Printf("create Rotation keys correspond to Transposition, MatrixMult and RowTotalSum\n")
		Rotation4Transpose := TransposeLT.Rotations4ArithmeticSeq()
		Rotation4MtrxMult := DSigmaLTs[0].Rotations4ArithmeticSeq()
		Rotation4MtrxMult = append(Rotation4MtrxMult, DTauLTs[0].Rotations()...)
		Rotation4MtrxMult = append(Rotation4MtrxMult, -d+d*d)
		galk4Transpose = kgen.GenRotationKeysForRotations(Rotation4Transpose, false, sk)
		galk4MtrxMult = kgen.GenRotationKeysForRotations(Rotation4MtrxMult, false, sk)
		Rotation4RowTotalSum := make([]int, d)
		for i := d; i < d*d; i = (i << 1) {
			Rotation4RowTotalSum = append(Rotation4RowTotalSum, i)
		}
		// Rotation4RowTotalSum = append(Rotation4RowTotalSum, -d+d*d)
		galk4RowtotalSum = kgen.GenRotationKeysForRotations(Rotation4RowTotalSum, false, sk)

		// Store the rotation keys if they are freshly generated.
		_, err = RtkTracker4Transpose.StoreUpdateOne(galk4Transpose)
		if err != nil {
			panic(err)
		} else {
			fmt.Printf("Store galk4Transpose %d bytes\n", galk4Transpose.MarshalBinarySize())
		}
		_, err = RtkTracker4Transpose.StoreFinish()
		if err != nil {
			panic(err)
		}

		_, err = RtkTracker4MtrxMult.StoreUpdateOne(galk4MtrxMult)
		if err != nil {
			panic(err)
		} else {
			fmt.Printf("Store galk4MtrxMult %d bytes\n", galk4MtrxMult.MarshalBinarySize())
		}
		_, err = RtkTracker4MtrxMult.StoreFinish()
		if err != nil {
			panic(err)
		}

		_, err = RtkTracker4RowtotalSum.StoreUpdateOne(galk4RowtotalSum)
		if err != nil {
			panic(err)
		} else {
			fmt.Printf("Store galk4RowtotalSum %d bytes\n", galk4RowtotalSum.MarshalBinarySize())
		}
		_, err = RtkTracker4RowtotalSum.StoreFinish()
		if err != nil {
			panic(err)
		}
	}

	// Free the memory:
	galk4Transpose = nil
	galk4MtrxMult = nil
	galk4RowtotalSum = nil

	// create a Matrix A
	MatrixA := make([][]float64, d)
	for i := 0; i < d; i++ {
		MatrixA[i] = make([]float64, d)
		for j := 0; j < d; j++ {
			MatrixA[i][j] = float64(i*d+j) * 0.0001
		}
	}
	auxio.Print_matrix_f64_2d(MatrixA, d, d)
	// encode matrixA into a vector using row ordering
	var a []float64
	a, err = mtrxmult.Row_orderingInv(MatrixA)
	if err != nil {
		panic(err)
	}

	// encode into plaintext and encrypt :
	ptA := encoder.EncodeNew(a, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
	ctA := encryptor.EncryptNew(ptA)

	// Retrieve the Rotation keys for MatrixMult
	RtkTracker4MtrxMult = auxio.NewTracker4File("Rtk4MtrxMult")
	galk4MtrxMult = new(rlwe.RotationKeySet)
	_, err = RtkTracker4MtrxMult.ReadUpdateOne(galk4MtrxMult)
	if err != nil {
		panic(err)
	}

	// Compute MatrixMultX:
	fmt.Printf("Start Computing SquareMatrixMult\n")
	now := time.Now()
	var ctAA *rlwe.Ciphertext
	ctAA, err = mtrxmult.SquareMatrix_MultBSGS_DecomposeVer_NeedLTs_dbg(params, sk, rlk, galk4MtrxMult, ctA, ctA, DSigmaLTs, DTauLTs, ColShiftPTs, d)
	if err != nil {
		panic(err)
	}
	fmt.Printf("SqaureMatrix Mult completed in %s", time.Since(now))
	auxio.Quick_check_matrix(params, sk, ctAA, d, d)

	fmt.Printf("Compare with exact result of SquareMatrix Mult:\n")
	checkA, _ := mtrxmult.SquareMatrix_product_permute_version(MatrixA, MatrixA)
	if err != nil {
		panic(err)
	}
	auxio.Print_matrix_f64_2d(checkA, d, d)
}

func ErrorCtrlX_check(loadkeysFromDisk bool) {
	var err error
	// initialize encryption scheme params.
	var params ckks.Parameters
	var sk *rlwe.SecretKey
	var pk *rlwe.PublicKey
	var rlk *rlwe.RelinearizationKey
	var kgen rlwe.KeyGenerator

	// generate Keys, we can load keys from disk or generate it.
	if loadkeysFromDisk {
		params, sk, pk, rlk, err = LoadKeys_check()
		if err != nil {
			panic(err)
		}
	} else {
		params, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
			LogN: 15,
			Q: []uint64{1152921504598720513,
				1099490000897, 1099498258433,
				1099499175937, 1099499569153,
				1099500617729}, // 60 + 40*5
			P:            []uint64{1152921504606584833}, // 60
			LogSlots:     14,
			DefaultScale: 1 << 40,
		})
		if err != nil {
			panic(err)
		}
		kgen = ckks.NewKeyGenerator(params)
		sk, pk = kgen.GenKeyPair()
		// Relinearization key
		rlk = kgen.GenRelinearizationKey(sk, 1)

		// Print Keys' size:
		fmt.Printf("Encryption Scheme Params occupy %d bytes\n", params.MarshalBinarySize())
		fmt.Printf("Secret Key:%d bytes, Public Key:%d bytes\n", sk.MarshalBinarySize(), pk.MarshalBinarySize())
		fmt.Printf("Relinearization Key: %d bytes\n", rlk.MarshalBinarySize())
	}

	// Encryptor
	encryptor := ckks.NewEncryptor(params, pk)

	// encoder
	encoder := ckks.NewEncoder(params)

	// evaluator
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})

	// Scales
	scales := make([]rlwe.Scale, params.MaxLevel()+1)
	for i := params.MaxLevel(); i >= 0; i-- {
		if i == params.MaxLevel() {
			scales[i] = rlwe.NewScale(float64(params.RingQ().Modulus[params.MaxLevel()]))
		} else {
			scale_qip1 := rlwe.NewScale(float64(params.RingQ().Modulus[i+1]))
			scales[i] = (scales[i+1].Mul(scales[i+1])).Div(scale_qip1)
		}
		scales[i] = rlwe.NewScale(1 << 40)
	}

	// create a Vec a
	// encode  a vector using row ordering
	var a = make([]float64, 1<<params.LogSlots())
	for i := range a {
		a[i] = float64(i % 128)
	}
	// encode into plaintext and encrypt :
	ptA := encoder.EncodeNew(a, params.MaxLevel(), scales[params.MaxLevel()], params.LogSlots())
	ctA := encryptor.EncryptNew(ptA)
	auxio.Quick_check_vector(params, sk, ctA) // Check the original Precision.
	ctApow2 := evaluator.MulRelinNew(ctA, ctA)
	err = evaluator.Rescale(ctApow2, scales[ctA.Level()-1], ctApow2)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Check the Precision after one Multiplication at Maximum Level.\n")
	auxio.Quick_check_vector(params, sk, ctApow2) // Check the Precision after one Multiplication at Maximum Level.

	ctApow4 := evaluator.MulRelinNew(ctApow2, ctApow2)
	err = evaluator.Rescale(ctApow4, scales[ctApow2.Level()-1], ctApow4)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Check the Precision after one Multiplication at Maximum Level -1.\n")
	auxio.Quick_check_vector(params, sk, ctApow4) // Check the Precision after one Multiplication at Maximum Level -1 .

	ctA_duplicate := ctA.CopyNew()
	ctA_duplicate.Resize(ctA_duplicate.Degree(), ctApow4.Level()+1)
	evaluator.SetScale(ctA_duplicate, scales[ctApow4.Level()])
	ctApow5 := evaluator.MulRelinNew(ctApow4, ctA_duplicate)
	err = evaluator.Rescale(ctApow5, scales[ctApow4.Level()-1], ctApow5)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Check the Precision after a Multiplication with ciphertexts at different levels.\n")
	auxio.Quick_check_vector(params, sk, ctApow5) // Check the Precision after a Multiplication with ciphertexts at different levels.

}
