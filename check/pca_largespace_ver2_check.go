package check

import (
	"fmt"
	"io/ioutil"
	"math"
	"path"
	"runtime/debug"
	"sync"
	"time"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	auxio "project1-fhe_extension_v1.0/auxiliary_io"
	mtrxmult "project1-fhe_extension_v1.0/matrix_mult"
)

func PCA_largeSpace_Ver2_check(loadkeysFromDisk bool) {
	var err error
	matrixCols := 128
	matrixRows := 128
	rows := 200
	cols := 257

	// BSGSRatio:
	/*
		// BGSGRatio4Transpose := 2
		BSGSRatio4Sigma := 2
		BSGSRatio4Tau := 2
	*/
	// Maximum Routines:
	maxSubRoutes := 1

	// N1
	DSigma_BSGSN1 := 16
	DTau_BSGSN1 := 16
	N1Trans := 16

	// runtime.GOMAXPROCS(6)
	d := 1 << 7

	// MaxDiagVec
	SigmaMaxDiagVec := d / 4
	TauMaxDiagVec := (d / 4) * d

	// necessary file paths:
	var KeysPath = "Keys"
	var Rtk4Transpose_path = "Rtk4Transpose"
	var Rtk4MtrxMult_path = "Rtk4MtrxMult"
	var Rtk4RowtotalSum_path = "Rtk4RowtotalSum"
	var ctCov_path = "ctCov"
	var ctMean_path = "ctMean"
	var csvfilepath = "mnist.csv"
	var ctDtau_path = "ctDtau" // will become a very large file...

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
		params, err = ckks.NewParametersFromLiteral(PN15QP720S14)
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
		keysTracker := auxio.NewTracker(KeysPath, -1)
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

	inputLevel := 9

	// Encryptor
	encryptor := ckks.NewEncryptor(params, pk)

	// Decryptor
	// decryptor := ckks.NewDecryptor(params, sk)

	// encoder
	encoder := ckks.NewEncoder(params)

	// evaluator
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})

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
	DSigmaDiagonalMaps1, DSigmaDiagonalMaps2, err := mtrxmult.GenSigmaDiagnalDecomposeMatrices(d, SigmaMaxDiagVec)
	if err != nil {
		panic(err)
	}

	// create the Map of Tau LinearTransformation
	fmt.Printf("create the Map of Tau LinearTransformation\n")
	DTauDiagonalMaps, err := mtrxmult.GenTauDiagonalDecomposeMatrices(d, TauMaxDiagVec)
	if err != nil {
		panic(err)
	}

	// create the Decomposed Maps and LinearTransform object of Sigma LinearTransformation
	fmt.Printf("create the Decomposed Maps and LinearTransform object of Sigma LinearTransformation\n")
	DSigmaLTs1 := make([]ckks.LinearTransform, len(DSigmaDiagonalMaps1))
	DSigmaLTs2 := make([]ckks.LinearTransform, len(DSigmaDiagonalMaps2))
	for i := range DSigmaLTs1 {
		if i == len(DSigmaLTs1)-1 {
			DSigmaLTs1[i] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DSigmaDiagonalMaps1[i], params.MaxLevel(), params.DefaultScale(), DSigma_BSGSN1, 1, params.LogSlots())
		} else {
			DSigmaLTs1[i] = ckks.GenLinearTransform(encoder, DSigmaDiagonalMaps1[i], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		}
	}
	for i := range DSigmaLTs2 {
		if i == len(DSigmaLTs2)-1 {
			DSigmaLTs2[i] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DSigmaDiagonalMaps2[i], params.MaxLevel(), params.DefaultScale(), DSigma_BSGSN1, 1, params.LogSlots())
		} else {
			DSigmaLTs2[i] = ckks.GenLinearTransform(encoder, DSigmaDiagonalMaps2[i], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		}
	}

	// create the Decomposed Maps of Tau LinearTransformation
	fmt.Printf("create the Decomposed Maps of Tau LinearTransformation\n")
	DTauLTs := make([]ckks.LinearTransform, len(DTauDiagonalMaps))
	for i := range DTauLTs {
		if i == len(DTauLTs)-1 {
			DTauLTs[i] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DTauDiagonalMaps[i], params.MaxLevel(), params.DefaultScale(), DTau_BSGSN1, d, params.LogSlots())
		} else {
			DTauLTs[i] = ckks.GenLinearTransform(encoder, DTauDiagonalMaps[i], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		}
	}

	// create the Map of Plaintext for ColShift Diagonals
	fmt.Printf("create the Map of Plaintext for ColShift Diagonals\n")
	ColShiftLTs := make([]ckks.LinearTransform, 0)
	for k := 1; k < d; k++ {
		U, err := mtrxmult.Gen_colShift_diagonalVectors(d, k)
		if err != nil {
			panic(err)
		}
		LT := ckks.GenLinearTransform(encoder, U, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		ColShiftLTs = append(ColShiftLTs, LT)
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
		RtkTracker4Transpose = auxio.NewTracker4File(Rtk4Transpose_path)
		RtkTracker4MtrxMult = auxio.NewTracker4File(Rtk4MtrxMult_path)
		RtkTracker4RowtotalSum = auxio.NewTracker4File(Rtk4RowtotalSum_path)
	} else {
		// Create File Tracker to Store or load these keys
		fmt.Printf("Init File Tracker to Store or Load the Rotation keys\n")
		RtkTracker4Transpose = auxio.NewTracker(Rtk4Transpose_path, -1)
		RtkTracker4MtrxMult = auxio.NewTracker(Rtk4MtrxMult_path, -1)
		RtkTracker4RowtotalSum = auxio.NewTracker(Rtk4RowtotalSum_path, -1)

		// kgen = ckks.NewKeyGenerator(params)

		// Create Rotation keys correspond to Transposition, MatrixMult and RowTotalSum
		fmt.Printf("create Rotation keys correspond to Transposition, MatrixMult and RowTotalSum\n")
		fmt.Printf("creating Rotation keys for Transposition\n")
		Rotation4Transpose := TransposeLT.Rotations4ArithmeticSeq()
		galk4Transpose = kgen.GenRotationKeysForRotations(Rotation4Transpose, false, sk)

		fmt.Printf("creating Rotation keys for MatrixMult\n")
		var MtrxMultRotmap = make(map[int]bool)
		for k := range DSigmaLTs1 {
			for _, i := range DSigmaLTs1[k].Rotations4ArithmeticSeq() {
				if !MtrxMultRotmap[i] {
					MtrxMultRotmap[i] = true
				}
			}
		}
		for k := range DSigmaLTs2 {
			for _, i := range DSigmaLTs2[k].Rotations4ArithmeticSeq() {
				if !MtrxMultRotmap[i] {
					MtrxMultRotmap[i] = true
				}
			}
		}
		for k := range DTauLTs {
			for _, i := range DTauLTs[k].Rotations4ArithmeticSeq() {
				if !MtrxMultRotmap[i] {
					MtrxMultRotmap[i] = true
				}
			}
		}
		var Rotation4MtrxMult = make([]int, 0)
		for rot := range MtrxMultRotmap {
			Rotation4MtrxMult = append(Rotation4MtrxMult, rot)
		}
		Rotation4MtrxMult = append(Rotation4MtrxMult, -d+d*d)
		galk4MtrxMult = kgen.GenRotationKeysForRotations(Rotation4MtrxMult, false, sk)

		fmt.Printf("creating Rotation keys for RowtotalSum\n")
		Rotation4RowTotalSum := make([]int, d)
		for i := d; i < d*d; i = (i << 1) {
			Rotation4RowTotalSum = append(Rotation4RowTotalSum, i)
		}
		// Rotation4RowTotalSum = append(Rotation4RowTotalSum, -d+d*d)
		galk4RowtotalSum = kgen.GenRotationKeysForRotations(Rotation4RowTotalSum, false, sk)

		// Store the rotation keys if they are freshly generated.
		fmt.Printf("Store the rotation keys since they are freshly generated.\n")
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
	debug.FreeOSMemory()

	// Construct cols-1/matrixCols x cols-1/matrixCols ciphertexts:
	Cov := make([][]*rlwe.Ciphertext, (cols-1)/matrixCols)
	for i := 0; i < (cols-1)/matrixCols; i++ {
		Cov[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	}
	// Construct [rows/d] matrices with size d x d
	ORSB1 := make([][][]float64, int(math.Ceil(float64(rows)/float64(d))))
	for i := 0; i < int(math.Ceil(float64(rows)/float64(d))); i++ {
		ORSB1[i] = make([][]float64, d)
		for j := 0; j < d; j++ {
			ORSB1[i][j] = make([]float64, d)
		}
	}

	// Construct ((cols-1)/matrixCols)^2 Ciphertexts representing a cols x cols Covariance Matrix
	ctCov := make([][]*rlwe.Ciphertext, (cols-1)/matrixCols)
	ctMean := make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	ctMMT := make([][]*rlwe.Ciphertext, (cols-1)/matrixCols) // ciphertext for MeanVec · MeanVec^T
	for i := 0; i < (cols-1)/matrixCols; i++ {
		ctCov[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
		ctMMT[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	}
	// We need some more temp Ciphertext to help...
	// var ctTemp *rlwe.Ciphertext
	// var ctSum *rlwe.Ciphertext
	// ctSum := encryptor.EncryptNew(auxio.Encode_single_float64(params, float64(0), params.MaxLevel(), rlwe.NewScale(1)))
	// ctSum := encryptor.EncryptZeroNew(params.MaxLevel())                                   // ciphertext for inner product of two Rows of Subblock
	ctORSBT := make([]*rlwe.Ciphertext, int(math.Ceil(float64(rows)/float64(matrixRows))))               // ciphetexts for One Row of Subblock Transposed
	ctORSB := make([]*rlwe.Ciphertext, int(math.Ceil(float64(rows)/float64(matrixRows))))                // ciphetexts for One Row of Subblock.
	ctORSBT_SigmaDecomp := make([][]*rlwe.Ciphertext, int(math.Ceil(float64(rows)/float64(matrixRows)))) // Sigma Decomposed ciphetexts for One Row of Subblock Transposed.
	ctARSB_TauDecomp := make([][][]*rlwe.Ciphertext, int(math.Ceil(float64(rows)/float64(matrixRows))))
	// ctSubRbuff := make([]*rlwe.Ciphertext, maxSubRoutes)                                   // subroutine buffer.
	// ctOCSB := make([]*rlwe.Ciphertext, int(math.Ceil(float64(cols)/float64(matrixCols))))
	var ctOCSBSum []*rlwe.Ciphertext
	// ctOCSBSum := make([]*rlwe.Ciphertext, int(math.Ceil(float64(cols)/float64(matrixCols))))
	IdentityVec := make([]float64, 1<<params.LogSlots())
	for i := 0; i < len(IdentityVec); i++ {
		IdentityVec[i] = 0.5
	}
	// ptIdentity := encoder.EncodeNew(IdentityVec, params.MaxLevel(), params.DefaultScale(), params.LogSlots())

	// Need to Store some of the ciphertexts into disk.
	var Covtracker = auxio.NewTracker(ctCov_path, -1)
	var Meantracker = auxio.NewTracker(ctMean_path, -1)
	var Dtautracker = auxio.NewTracker(ctDtau_path, -1)
	// var Datatracker = auxio.NewTracker("ctData", -1)

	// Start Computing X*X^T and Mean
	fmt.Printf("Start Computing X*X^T and Mean\n")
	var elapsed time.Duration
	var now time.Time
	var wg sync.WaitGroup
	for k := 1; k < cols; k += matrixCols {
		now = time.Now()
		fmt.Printf("Extract the %d th OneRowOfSubblock, and by the way compute their Tau transformation result in all row shift\n", (k-1)/matrixCols)
		subrouteNum := maxSubRoutes
		for j := 1; j < rows+1; j += matrixRows {
			wg.Add(1)
			subrouteNum--
			go func(d int, j int) {
				now := time.Now()
				defer wg.Done()
				var sb [][]string            // One subblock in string type
				var sbf64 [][]float64        // One subblick in float64 type
				var MatrixSB []float64       // Row odering Vector for one subblock
				var ptMtrxSB *rlwe.Plaintext // Plaintext for one subblock
				var ctMtrxSB *rlwe.Ciphertext
				encoder_subroute := ckks.NewEncoder(params)
				encryptor_subroute := ckks.NewEncryptor(params, pk)
				if j+matrixRows > rows+1 {
					sb, err = auxio.GetOneCSVSubblock(csvfilepath, j, k, rows-j+1, matrixCols)
				} else {
					sb, err = auxio.GetOneCSVSubblock(csvfilepath, j, k, matrixRows, matrixCols)
				}
				sbf64, err = auxio.Switch2d_str2f64(sb)
				if err != nil {
					panic(err)
				}
				MatrixSB, err = mtrxmult.Row_orderingInvZeroPad(sbf64, d)
				ptMtrxSB = encoder_subroute.EncodeNew(MatrixSB, inputLevel, params.DefaultScale(), params.LogSlots())
				ctMtrxSB = encryptor_subroute.EncryptNew(ptMtrxSB)
				ctORSB[j/matrixRows] = ctMtrxSB.CopyNew()
				fmt.Printf("the %d th routine done in %s \n", j, time.Since(now))
			}(d, j)
			if subrouteNum <= 0 {
				wg.Wait()
				subrouteNum = maxSubRoutes
			}
		}
		wg.Wait()
		fmt.Printf("Extraction complete in %s with one ciphertext: %d bytes, take up space %d MB in total\n", time.Since(now), ctORSB[0].MarshalBinarySize(), ctORSB[0].MarshalBinarySize()*len(ctORSBT)/(1024*1024)) // dbg testing
		elapsed += time.Since(now)

		// Retrieve the Rotation keys for MatrixMult
		RtkTracker4MtrxMult = auxio.NewTracker4File(Rtk4MtrxMult_path)
		galk4MtrxMult = new(rlwe.RotationKeySet)
		_, err = RtkTracker4MtrxMult.ReadUpdateOne(galk4MtrxMult)
		if err != nil {
			panic(err)
		}
		fmt.Printf("Retrieve the Rotation keys with %d MB for MatrixMult\n", galk4MtrxMult.MarshalBinarySize()/(1024*1024))

		fmt.Printf("Start to compute and Store one ORSB's TauDecomposition into file\n")
		now = time.Now()
		ctARSB_TauDecomp[(k-1)/matrixCols], err = mtrxmult.Tau_linearTransform_MultiThread_dbg(params, sk, rlk, galk4MtrxMult, ctORSB, DTauLTs, d, maxSubRoutes)
		/*
			for j := 0; j < len(ctORSB); j += maxSubRoutes {
				var ctDTaus [][]*rlwe.Ciphertext
				if j+maxSubRoutes > len(ctORSB) {
					ctDTaus, err = mtrxmult.Tau_linearTransform_MultiThread_dbg(params, sk, rlk, galk4MtrxMult, ctORSB[j:len(ctORSB)], DTauLTs, d, maxSubRoutes)
				} else {
					ctDTaus, err = mtrxmult.Tau_linearTransform_MultiThread_dbg(params, sk, rlk, galk4MtrxMult, ctORSB[j:j+maxSubRoutes], DTauLTs, d, maxSubRoutes)
				}
				if err != nil {
					panic(err)
				}
				for _, cts := range ctDTaus {
					for _, ct := range cts {
						Dtautracker.StoreUpdateOne(ct)
					}
				}
			}
		*/
		fmt.Printf("Compute and Store one ORSB's TauDecomposition into disk in %s\n", time.Since(now))
		elapsed += time.Since(now)
	}
	Dtautracker.StoreFinish()
	galk4MtrxMult = nil
	debug.FreeOSMemory()

	for k := 1; k < cols; k += matrixCols {
		// Load the Rotation keys for Transposition
		RtkTracker4Transpose = auxio.NewTracker4File(Rtk4Transpose_path)
		galk4Transpose = new(rlwe.RotationKeySet)
		_, err = RtkTracker4Transpose.ReadUpdateOne(galk4Transpose)
		if err != nil {
			panic(err)
		}
		fmt.Printf("Retrieve the Rotation keys with %d MB for Transposition\n", galk4Transpose.MarshalBinarySize()/(1024*1024))

		// Generate evaluation key for transposition and other operation:
		elk := rlwe.EvaluationKey{Rlk: rlk, Rtks: galk4Transpose}
		// Extract the Transposition of OneRowOfSubblock.
		now = time.Now()
		fmt.Printf("Extract the Transposition of %d th OneRowOfSubblock, and by the way compute their Sigma transformation result in all column shift\n", (k-1)/matrixCols)
		subrouteNum := maxSubRoutes
		for j := 1; j < rows+1; j += matrixRows {
			// extract a matrixRows x matrixCols matrix
			// MultiThread Mode
			wg.Add(1)
			subrouteNum--
			go func(d int, ORSB1 [][][]float64, j int) {
				now := time.Now()
				defer wg.Done()
				var sb [][]string            // One subblock in string type
				var sbf64 [][]float64        // One subblick in float64 type
				var MatrixSB []float64       // Row odering Vector for one subblock
				var ptMtrxSB *rlwe.Plaintext // Plaintext for one subblock
				var ctMtrxSB *rlwe.Ciphertext
				encoder_subroute := ckks.NewEncoder(params)
				encryptor_subroute := ckks.NewEncryptor(params, pk)
				evaluator_subroute := ckks.NewEvaluator(params, elk)
				if j+matrixRows > rows+1 {
					sb, err = auxio.GetOneCSVSubblock(csvfilepath, j, k, rows-j+1, matrixCols)
				} else {
					sb, err = auxio.GetOneCSVSubblock(csvfilepath, j, k, matrixRows, matrixCols)
				}
				sbf64, err = auxio.Switch2d_str2f64(sb)
				if err != nil {
					panic(err)
				}
				MatrixSB, err = mtrxmult.Row_orderingInvZeroPad(sbf64, d)
				ptMtrxSB = encoder_subroute.EncodeNew(MatrixSB, inputLevel, params.DefaultScale(), params.LogSlots())
				ctMtrxSB = encryptor_subroute.EncryptNew(ptMtrxSB)
				// auxio.Quick_decode_matrix_full(params, ptMtrxSB, d, d)
				ctORSB[j/matrixRows] = ctMtrxSB                                                                       // dbg testing
				ctORSBT[j/matrixRows] = evaluator_subroute.LinearTransform4ArithmeticSeqNew(ctMtrxSB, TransposeLT)[0] // encrypt the Transposed Plaintext for one subblock into the OneRowOfSubblock ciphertext list.
				if err != nil {
					panic(err)
				}
				err = evaluator_subroute.Rescale(ctORSBT[j/matrixRows], params.DefaultScale(), ctORSBT[j/matrixRows])
				if err != nil {
					panic(err)
				}
				// copy(ORSB1[j/matrixRows], sbf64)
				fmt.Printf("the %d th routine done in %s \n", j, time.Since(now))
				fmt.Printf("ctORSBT at level %d, Scale %f\n", ctORSBT[j/matrixRows].Level(), math.Log2(ctORSBT[j/matrixRows].Scale.Float64()))
				// auxio.Quick_check_matrix_full(params, sk, ctORSBT[j/matrixRows], d, d)
				// return
			}(d, ORSB1, j)
			if subrouteNum <= 0 {
				wg.Wait()
				subrouteNum = maxSubRoutes
			}
		}
		wg.Wait()
		fmt.Printf("Extractiong & Transposition complete in %s with one ciphertext: %d bytes, take up space %d MB in total\n", time.Since(now), ctORSBT[0].MarshalBinarySize(), ctORSBT[0].MarshalBinarySize()*len(ctORSBT)*2/(1024*1024)) // dbg testing
		elapsed += time.Since(now)
		// Free the Rotation keys for Transposition
		galk4Transpose = nil

		// We then compute the Sum of One Row of Subblocks to prepare for the Aggregation of Mean vector computation
		now = time.Now()
		fmt.Printf("We then compute the Sum of One Row of Subblocks to prepare for the Aggregation of Mean vector computation\n")
		ctMean[(k-1)/matrixCols] = encryptor.EncryptZeroNew(params.MaxLevel()) // dbg testing
		for j := 1; j < rows+1; j += matrixRows {
			evaluator.Add(ctMean[(k-1)/matrixCols], ctORSB[(j-1)/matrixRows], ctMean[(k-1)/matrixCols])
			ctORSB[(j-1)/matrixRows] = nil // free the space.
		}
		fmt.Printf("Sum of Subblocks complete in %s\n", time.Since(now))
		elapsed += time.Since(now)
		debug.FreeOSMemory()

		// Retrieve the Rotation keys for MatrixMult
		RtkTracker4MtrxMult = auxio.NewTracker4File(Rtk4MtrxMult_path)
		galk4MtrxMult = new(rlwe.RotationKeySet)
		_, err = RtkTracker4MtrxMult.ReadUpdateOne(galk4MtrxMult)
		if err != nil {
			panic(err)
		}
		fmt.Printf("Retrieve the Rotation keys with %d MB for MatrixMult", galk4MtrxMult.MarshalBinarySize()/(1024*1024))

		// Compute the Sigma Decomposition of the ciphertexts in ctORSBT:
		//now = time.Now()
		//fmt.Printf("Compute the Sigma Decomposition of the ciphertexts in ctORSBT\n")
		// ctORSBT_SigmaDecomp, err = mtrxmult.Sigma_linearTransform_MultiThread_dbg(params, sk, rlk, galk4MtrxMult, ctORSBT, DSigmaLTs1, DSigmaLTs2, ColShiftLTs, d, len(ctORSBT))
		//fmt.Printf("Sigma Decomposition of the ciphertexts in ctORSBT Done in %s\n", time.Since(now))
		//elapsed += time.Since(now)

		/* dbg testing
		for i := 0; i < len(ctORSBT); i++ {
			auxio.Quick_check_matrix(params, sk, ctORSBT[i], d, d)
		}
		*/

		// Start Computing Matrix Mult for each element in ctORSBT, each of them will be multiplied by elements in ctOCSB in paralle.
		now = time.Now()
		fmt.Printf("Start Computing Matrix Mult for each element in ctORSBT, each of them will be multiplied by elements in ctOCSB in paralle.\n")
		ctOCSBSum = make([]*rlwe.Ciphertext, (k-1)/matrixCols+1)

		subrouteNum_outer := maxSubRoutes / ((k-1)/matrixCols + 1)
		cond := sync.NewCond(&sync.Mutex{})
		j := 1
		for idx := range ctORSBT_SigmaDecomp {
			if ctORSBT_SigmaDecomp[idx] == nil {
				wg.Wait()
				nowtmp := time.Now()
				var ctORSBT_SigmaDecompTmp [][]*rlwe.Ciphertext
				if idx+maxSubRoutes > len(ctORSBT) {
					fmt.Printf("Computing SigmaLinTrans for ctORSBT[%d:%d]\n", idx, len(ctORSBT))
					ctORSBT_SigmaDecompTmp, err = mtrxmult.Sigma_linearTransform_MultiThread_dbg(params, sk, rlk, galk4MtrxMult, ctORSBT[idx:], DSigmaLTs1, DSigmaLTs2, ColShiftLTs, d, maxSubRoutes)
				} else {
					fmt.Printf("Computing SigmaLinTrans for ctORSBT[%d:%d]\n", idx, idx+maxSubRoutes)
					ctORSBT_SigmaDecompTmp, err = mtrxmult.Sigma_linearTransform_MultiThread_dbg(params, sk, rlk, galk4MtrxMult, ctORSBT[idx:idx+maxSubRoutes], DSigmaLTs1, DSigmaLTs2, ColShiftLTs, d, maxSubRoutes)
				}
				if err != nil {
					panic(err)
				}
				fmt.Printf("ctORSBT_SigmaDecompTmp occupying %d MB\n", ctORSBT_SigmaDecompTmp[0][0].MarshalBinarySize()*len(ctORSBT_SigmaDecompTmp)*len(ctORSBT_SigmaDecompTmp[0])/(1024*1024))
				copy(ctORSBT_SigmaDecomp[idx:idx+len(ctORSBT_SigmaDecompTmp)], ctORSBT_SigmaDecompTmp)
				fmt.Printf("SigmaLinTrans complete in %s\n", time.Since(nowtmp))
				for i := 0; i < idx; i++ {
					ctORSBT_SigmaDecomp[i] = nil
				}
				nowtmp = time.Now()
				debug.FreeOSMemory()
				elapsed -= time.Since(nowtmp)
			}
			wg.Add(1)
			subrouteNum_outer--
			go func(idx int, j int) {
				defer wg.Done()
				// ctOCSB := make([]*rlwe.Ciphertext, (k-1)/matrixCols+1)
				ctOCSB_TauDmp := make([][]*rlwe.Ciphertext, (k-1)/matrixCols+1)
				var wg_subroutine sync.WaitGroup
				// Retrieve OCSB of one row
				subrouteNum_inner := ((k-1)/matrixCols + 1)
				for i := 1; i <= k; i += matrixCols {
					fmt.Printf("Retrieve dataset[%d][%d], and load it into ctOCSB\n", idx, (i-1)/matrixCols)
					ctOCSB_TauDmp[(i-1)/matrixCols] = make([]*rlwe.Ciphertext, d)
					wg_subroutine.Add(1)
					subrouteNum_inner--
					go func(d int, i int, j int) {
						defer wg_subroutine.Done()
						/*
							Dtautracker_subroute := auxio.NewTracker4File(ctDtau_path)
							Dtautracker_subroute.ReadInit()
							Dtautracker_subroute.Locate(((i-1)/matrixCols)*len(ctORSBT_SigmaDecomp)*d + ((j-1)/matrixRows)*d)
							for idx := range ctOCSB_TauDmp[(i-1)/matrixCols] {
								ctOCSB_TauDmp[(i-1)/matrixCols][idx] = ckks.NewCiphertext(params, 1, params.MaxLevel())
								Dtautracker_subroute.ReadUpdateOne(ctOCSB_TauDmp[(i-1)/matrixCols][idx])
							}
						*/
						ctOCSB_TauDmp[(i-1)/matrixCols] = ctARSB_TauDecomp[(i-1)/matrixCols][(j-1)/matrixRows]
					}(d, i, j)
					if subrouteNum_inner <= 0 {
						wg_subroutine.Wait()
						subrouteNum_inner = ((k-1)/matrixCols + 1)
					}
				}
				wg_subroutine.Wait()
				/*
					for i := range ctOCSB_TauDmp {
						fmt.Printf("ctOCSB_TauDmp[%d]:\n", i)
						for j := range ctOCSB_TauDmp[i] {
							auxio.Quick_check_matrix(params, sk, ctOCSB_TauDmp[i][j], d, d)
						}
					}
				*/
				// Compute MtrxMult between element in ctORSBT and elements in ctOCSB:
				fmt.Printf("Compute MtrxMult between ctORSBT[%d] and elements in %d-th ctOCSB, and add them to the previous result in ctOCSBSum:\n", idx, idx)
				ctOCSBTemp, err := mtrxmult.SquareMatrix_MultBSGS_DecomposeVer4_NeedLTs_dbg(params, sk, rlk, ctORSBT_SigmaDecomp[idx], ctOCSB_TauDmp, d, len(ctOCSB_TauDmp))
				// ctOCSBTemp, err := mtrxmult.SquareMatrix_MultBSGS_DecomposeVer3_NeedLTs_dbg(params, sk, rlk, galk4MtrxMult, ctORSBT_SigmaDecomp[idx], ctOCSB, DTauLTs, d, len(ctOCSB))
				// ctOCSBTemp, err := mtrxmult.SquareMatrix_MultBSGS_DecomposeVer2_NeedLTs_dbg(params, sk, rlk, galk4MtrxMult, ct, ctOCSB, DSigmaLTs1, DSigmaLTs2, DTauLTs, ColShiftLTs, d, len(ctOCSB))
				if err != nil {
					panic(err)
				}
				cond.L.Lock()
				if ctOCSBSum[0] == nil {
					fmt.Printf("Let ctOCSBTemp[%d] as the first value of ctOCBSum\n", idx)
					for i := range ctOCSBTemp {
						ctOCSBSum[i] = ctOCSBTemp[i].CopyNew()
						//fmt.Printf("First ctOCSBSum[%d] value:\n", i)
						//auxio.Quick_check_matrix(params, sk, ctOCSBSum[i], d, d)
					}
				} else {
					fmt.Printf("Adding ctOCSBTemp[%d] into ctOCBSum\n", idx)
					for i := range ctOCSBTemp {
						//fmt.Printf("New value that will be added to ctOCSBSum[%d]:\n", i)
						//auxio.Quick_check_matrix(params, sk, ctOCSBTemp[i], d, d)
						evaluator.Add(ctOCSBTemp[i], ctOCSBSum[i], ctOCSBSum[i])
						//fmt.Printf("After adding a new one, ctOCSBSum[%d] value:\n", i)
						//auxio.Quick_check_matrix(params, sk, ctOCSBSum[i], d, d)
					}
				}
				cond.L.Unlock()
			}(idx, j)
			if subrouteNum_outer <= 0 {
				wg.Wait()
				subrouteNum_outer = maxSubRoutes / ((k-1)/matrixCols + 1)
			}
			j += matrixRows
		}
		wg.Wait()
		fmt.Printf("Scale the result in ctOCSBSum, and load them into ctCov:\n")
		for i := range ctOCSBSum {
			err = evaluator.Rescale(ctOCSBSum[i], params.DefaultScale(), ctOCSBSum[i])
			if err != nil {
				panic(err)
			}
			//fmt.Printf("ctCov[%d][%d] before scaling has scale %f, level %d\n", (k-1)/matrixCols, i, math.Log2(ctOCSBSum[i].Scale.Float64()), ctOCSBSum[i].Level())
			// auxio.Quick_check_matrix(params, sk, ctOCSBSum[i], d, d)
			// ctOCSBSum[i] = encryptor.EncryptNew(encoder.EncodeNew(auxio.DecryptDecode(params, sk, ctOCSBSum[i]), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
			evaluator.MultByConst(ctOCSBSum[i], 1.0/float64(rows), ctOCSBSum[i])
			err = evaluator.Rescale(ctOCSBSum[i], params.DefaultScale(), ctOCSBSum[i])
			if err != nil {
				panic(err)
			}
			fmt.Printf("ctCov[%d][%d] now has scale %f, level %d\n", (k-1)/matrixCols, i, math.Log2(ctOCSBSum[i].Scale.Float64()), ctOCSBSum[i].Level())
			ctCov[(k-1)/matrixCols][i] = ctOCSBSum[i].CopyNew()
			auxio.Quick_check_matrix(params, sk, ctCov[(k-1)/matrixCols][i], d, d)
		}

		galk4MtrxMult = nil
		fmt.Printf("X*X^T %d/7 process Done in %s\n", (k-1)/matrixCols+1, time.Since(now))
		elapsed += time.Since(now)
		debug.FreeOSMemory()

		// At the same time we computie 1/7 of the Mean vector by Aggregating the OneRowOfSubblock ciphertexts and RowTotalSum.
		now = time.Now()
		fmt.Printf("At the same time we computie 1/7 of the Mean vector by Aggregating the OneRowOfSubblock ciphertexts and RowTotalSum.\n")
		// First Retrieve the RotationKeys for RowtotalSum :
		RtkTracker4RowtotalSum = auxio.NewTracker4File(Rtk4RowtotalSum_path)
		galk4RowtotalSum = new(rlwe.RotationKeySet)
		_, err = RtkTracker4RowtotalSum.ReadUpdateOne(galk4RowtotalSum)
		if err != nil {
			panic(err)
		}
		ctMean[(k-1)/matrixCols], err = mtrxmult.Matrix_Aggregate_withRowOrder(params, sk, galk4RowtotalSum, ctMean[(k-1)/matrixCols], d, d, 0)
		if err != nil {
			panic(err)
		}
		// auxio.Quick_check_matrix_full(params, sk, ctMean[(k-1)/matrixCols], d, d) // dbg testing
		evaluator.MultByConst(ctMean[(k-1)/matrixCols], 1.0/float64(rows), ctMean[(k-1)/matrixCols])
		fmt.Printf("Check the MeanVec\n")
		auxio.Quick_check_matrix(params, sk, ctMean[(k-1)/matrixCols], d, d) // dbg testing
		fmt.Printf("Mean %d/7 process Done in %s\n", (k-1)/matrixCols+1, time.Since(now))
		elapsed += time.Since(now)
		// We want to store the ctMean into disk, to reduce the memory occupation...
		// Store the mean:
		Meantracker.StoreUpdateOne(ctMean[(k-1)/matrixCols])
		ctMean[(k-1)/matrixCols] = nil                  // release the memory.
		galk4RowtotalSum = nil                          // free the Rotation key for RowtotalSum
		for i := 0; i < len(ctORSBT_SigmaDecomp); i++ { // free the decomposition space
			ctORSBT_SigmaDecomp[i] = nil
		}
		debug.FreeOSMemory()

	}

	// Load the Rotation keys for Transposition
	fmt.Printf("Retrieve the Rotation keys for Transposition\n")
	RtkTracker4Transpose = auxio.NewTracker4File(Rtk4Transpose_path)
	galk4Transpose = new(rlwe.RotationKeySet)
	_, err = RtkTracker4Transpose.ReadUpdateOne(galk4Transpose)
	if err != nil {
		panic(err)
	}
	// Generate evaluation key for transposition and other operation:
	elk := rlwe.EvaluationKey{Rlk: rlk, Rtks: galk4Transpose}
	evaluator = ckks.NewEvaluator(params, elk)

	// Using the symmetric property of Covariance to compute the rest of its part.
	now = time.Now()
	fmt.Printf("Using the symmetric property of Covariance to compute the rest of its part.\n")
	for i := 0; i < len(ctCov); i++ {
		for j := 0; j < len(ctCov); j++ {
			if i < j { // in the previous loop we've got ctCov[k][i] where i<k
				err = evaluator.Rescale(ctCov[j][i], params.DefaultScale(), ctCov[j][i])
				fmt.Printf("The ctCov[%d][%d], it will be used to compute ctCov[%d][%d]\n", j, i, i, j)
				auxio.Quick_check_matrix(params, sk, ctCov[j][i], d, d)
				if err != nil {
					panic(err)
				}
				ctCov[i][j] = evaluator.LinearTransform4ArithmeticSeqNew(ctCov[j][i], TransposeLT)[0]
				err = evaluator.Rescale(ctCov[i][j], params.DefaultScale(), ctCov[i][j])
				if err != nil {
					panic(err)
				}
				fmt.Printf("Syncronize the level and scale of ctCov[%d][%d] and ctCov[%d][%d]\n", i, j, j, i)
				evaluator.SetScale(ctCov[j][i], ctCov[i][j].Scale) // Syncronize the level and scale of ctCov[i][j] and ctCov[j][i]
				fmt.Printf("ctCov[%d][%d] now has scale %f, level %d, same as ctCov[%d][%d] scale %f, level %d\n", j, i, math.Log2(ctCov[j][i].Scale.Float64()), ctCov[j][i].Level(), i, j, math.Log2(ctCov[i][j].Scale.Float64()), ctCov[i][j].Level())
			}
		}
	}
	for i := 0; i < len(ctCov); i++ {
		err = evaluator.Rescale(ctCov[i][i], params.DefaultScale(), ctCov[i][i])
		if err != nil {
			panic(err)
		}
		if len(ctCov) > 1 {
			evaluator.SetScale(ctCov[i][i], ctCov[0][1].Scale) // Syncronize the level and scale of ctCov[i][i] and arbitary ctCov[i][j] where i < j.
			fmt.Printf("ctCov[%d][%d] is synchronized to scale %f, level %d\n", i, i, math.Log2(ctCov[i][i].Scale.Float64()), ctCov[i][i].Level())
		}
		for j := 0; j < len(ctCov); j++ {
			fmt.Printf("ctCov[%d][%d] stored into disk:\n", i, j)
			auxio.Quick_check_matrix(params, sk, ctCov[i][j], d, d)
			_, err = Covtracker.StoreUpdateOne(ctCov[i][j])
			if err != nil {
				panic(err)
			}
		}
	}
	fmt.Printf("Rest subblocks of Covariance computed and synchronized in %s\n", time.Since(now))
	elapsed += time.Since(now)

	fmt.Printf("Covariance Process Done in %s\n", elapsed)

	// Finish all the storing process:
	_, err = Covtracker.StoreFinish()
	if err != nil {
		panic(err)
	}
	_, err = Meantracker.StoreFinish()
	if err != nil {
		panic(err)
	}

	// Release the ctORSBT's memory, since all of the subblock ciphertexts are stored in disk.
	for i := 0; i < len(ctORSBT); i++ {
		ctORSBT[i] = nil
	}

}

func PCA_Aggregate() {

	targetPath := "CovSum"
	// initialize encryption scheme params.
	var err error
	matrixCols := 112
	// matrixRows := 128
	// rows := 256
	cols := 785

	// Size of the Square Matrix.
	// d := 1 << 7

	var params ckks.Parameters
	// var sk *rlwe.SecretKey
	// var pk *rlwe.PublicKey
	var rlk *rlwe.RelinearizationKey
	// var kgen rlwe.KeyGenerator

	params, _, _, rlk, err = LoadKeys_check()
	if err != nil {
		panic(err)
	}

	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})

	filePath := "E:\\user-eric\\Go\\project1-fhe_extension_v1.0"
	CovPaths := make([]string, 0)
	result := GetAllFile(filePath)
	for i := 0; i < len(result); i++ {
		suffixName := GetSuffixName(result[i])
		if suffixName == "cov" {
			CovPaths = append(CovPaths, result[i])
		}
	}
	ctCov := make([][]*rlwe.Ciphertext, (cols-1)/cols)
	ctCov_tmp := make([][]*rlwe.Ciphertext, (cols-1)/cols)
	for i := range ctCov {
		ctCov[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
		ctCov_tmp[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
		for j := range ctCov[i] {
			ctCov[i][j] = new(rlwe.Ciphertext)
			ctCov_tmp[i][j] = new(rlwe.Ciphertext)
		}
	}
	for idx, p := range CovPaths {
		Covtracker := auxio.NewTracker4File(p)
		if idx == 0 {
			for i := 0; i < len(ctCov); i++ {
				//_, err = auxio.ReadCtVec4file(Covtracker.Fp, Covtracker.Hierarchy[0], ctCov[i], i*Covtracker.Hierarchy[0]*Covtracker.Hierarchy[1])
				for j := 0; j < len(ctCov); j++ {
					Covtracker.ReadUpdateOne(ctCov[i][j])
				}
			}
		} else {
			for i := 0; i < len(ctCov); i++ {
				//_, err = auxio.ReadCtVec4file(Covtracker.Fp, Covtracker.Hierarchy[0], ctCov[i], i*Covtracker.Hierarchy[0]*Covtracker.Hierarchy[1])
				for j := 0; j < len(ctCov); j++ {
					Covtracker.ReadUpdateOne(ctCov_tmp[i][j])
					evaluator.Add(ctCov[i][j], ctCov_tmp[i][j], ctCov[i][j])
				}
			}
		}
	}

	Covtracker := auxio.NewTracker(targetPath, -1)
	for i := 0; i < len(ctCov); i++ {
		//_, err = auxio.ReadCtVec4file(Covtracker.Fp, Covtracker.Hierarchy[0], ctCov[i], i*Covtracker.Hierarchy[0]*Covtracker.Hierarchy[1])
		for j := 0; j < len(ctCov); j++ {
			_, err = Covtracker.StoreUpdateOne(ctCov[i][j])
			if err != nil {
				panic(err)
			}
		}
	}
	Covtracker.StoreFinish()

}

// Get all files in the directory
func GetAllFile(pathname string) []string {
	rd, _ := ioutil.ReadDir(pathname)
	var filePath []string
	for _, fi := range rd {
		if fi.IsDir() {
			filePath = append(filePath, pathname+"/"+fi.Name())
			GetAllFile(pathname + fi.Name() + "\\")
		} else {
			filePath = append(filePath, pathname+"/"+fi.Name())
		}
	}
	return filePath
}

// Get suffix.
func GetSuffixName(filePath string) string {
	var fileSuffix string
	fileSuffix = path.Ext(filePath)
	return fileSuffix
}
