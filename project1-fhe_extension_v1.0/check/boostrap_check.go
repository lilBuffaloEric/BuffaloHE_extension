package check

import (
	"fmt"
	"math"
	"math/rand"
	"runtime/debug"
	"sync"
	"time"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	auxio "project1-fhe_extension_v1.0/auxiliary_io"
	mtrxmult "project1-fhe_extension_v1.0/matrix_mult"
	nonpolyfunc "project1-fhe_extension_v1.0/nonpoly_func"
)

func PCA_btpVersion_check(loadkeysFromDisk bool, loadPCAkeysFromDisk bool) {
	var err error
	matrixCols := 128
	matrixRows := 128
	rows := 60000
	cols := 257

	// BSGSRatio:
	/*
		// BGSGRatio4Transpose := 2
		BSGSRatio4Sigma := 2
		BSGSRatio4Tau := 2
	*/
	// Maximum Routines:
	maxSubRoutes := 5

	// N1
	DSigma_BSGSN1 := 16
	DTau_BSGSN1 := 16
	N1Trans := 16

	// runtime.GOMAXPROCS(6)
	d := 1 << 7

	// MaxDiagVec
	SigmaMaxDiagVec := d / 2
	TauMaxDiagVec := (d / 2) * d

	// necessary file paths:
	var KeysPath_btpVersion = "Keys_btpVer"
	var Rtk4Transpose_path = "Rtk4Transpose_btpVer"
	var Rtk4MtrxMult_path = "Rtk4MtrxMult_btpVer"
	var Rtk4RowtotalSum_path = "Rtk4RowtotalSum_btpVer"
	var ctCov_path = "ctCov_btpVer"
	var ctMean_path = "ctMean_btpVer"
	var csvfilepath = "mnist.csv"

	// initialize encryption scheme params.

	var ckksParams ckks.ParametersLiteral
	var btpParams bootstrapping.Parameters
	var params ckks.Parameters

	var sk *rlwe.SecretKey
	var pk *rlwe.PublicKey
	var rlk *rlwe.RelinearizationKey
	var btpevk bootstrapping.EvaluationKeys
	var kgen rlwe.KeyGenerator

	// generate Keys, we can load keys from disk or generate them.
	if loadkeysFromDisk {
		// In this function, we do not need to perform boostrapping, so we only pick up the keys we need in this section.
		params, btpParams, sk, pk, btpevk.Rlk, _, _, _, err = LoadKeys_btpVersion(KeysPath_btpVersion)
		if err != nil {
			panic(err)
		}
	} else {
		paramSet := bootstrapping.N16QP1788H32768H32
		ckksParams = paramSet.SchemeParams
		btpParams = paramSet.BootstrappingParams
		ckksParams.LogSlots = 14
		params, err = ckks.NewParametersFromLiteral(ckksParams)
		if err != nil {
			panic(err)
		}
		fmt.Println()
		fmt.Printf("CKKS parameters: logN = %d, logSlots = %d, H(%d; %d), logQP = %d, levels = %d, scale= 2^%f, sigma = %f \n", params.LogN(), params.LogSlots(), params.HammingWeight(), btpParams.EphemeralSecretWeight, params.LogQP(), params.QCount(), math.Log2(params.DefaultScale().Float64()), params.Sigma())

		kgen = ckks.NewKeyGenerator(params)
		sk, pk = kgen.GenKeyPair()
		fmt.Println()

		fmt.Println("Generating bootstrapping keys...")
		btpevk = bootstrapping.GenEvaluationKeys(btpParams, params, sk)
		fmt.Println("Done")

		// Print Keys' size:
		fmt.Printf("Encryption Scheme btpParams occupy %d bytes\n", btpParams.MarshalBinarySize())
		fmt.Printf("Encryption Scheme params    occupy %d bytes\n", params.MarshalBinarySize())
		fmt.Printf("Encryption Scheme secretkey occupy %d bytes\n", sk.MarshalBinarySize())
		fmt.Printf("Encryption Scheme publickey occupy %d bytes\n", pk.MarshalBinarySize())
		fmt.Printf("(btp) Relinearization Key   occupy %d bytes\n", btpevk.Rlk.MarshalBinarySize())
		fmt.Printf("(btp) Rotation Keys         occupy %d bytes\n", btpevk.Rtks.MarshalBinarySize())
		fmt.Printf("(btp) DtS Switching Key     occupy %d bytes\n", btpevk.SwkDtS.MarshalBinarySize())
		fmt.Printf("(btp) StD Switching Key     occupy %d bytes\n", btpevk.SwkStD.MarshalBinarySize())

		_, err = StoreKeys_btpVersion(KeysPath_btpVersion, params, btpParams, sk, pk, btpevk.Rlk, btpevk.Rtks, btpevk.SwkDtS, btpevk.SwkStD)
		if err != nil {
			panic(err)
		}
		btpevk.Rtks = nil
		btpevk.SwkDtS = nil
		btpevk.SwkStD = nil
		debug.FreeOSMemory()
	}

	// let rlk be a copy of btpevk.Rlk, they share the same memory address.
	rlk = btpevk.Rlk // History lefting problems... rlk should be deleted, but...

	// Encryptor
	encryptor := ckks.NewEncryptor(params, pk)

	// Decryptor
	// decryptor := ckks.NewDecryptor(params, sk)

	// encoder
	encoder := ckks.NewEncoder(params)

	// evaluator
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})

	// In the Boostrapping Version, we do not begin the performance of operations from the top level,
	// instead we will begin at the first level which is not included in bootstrap procedure, the number of this level can be
	// computed by:
	MaxLevel := btpParams.SlotsToCoeffsParameters.LevelStart - len(btpParams.SlotsToCoeffsParameters.ScalingFactor)
	inputLevel := 7

	// create the Map and LinearTransform object of Transpose LinearTransformation
	fmt.Printf("create the Map and LinearTransform object of Transpose LinearTransformation\n")
	var TransposeDiagonalMap map[int][]float64
	TransposeDiagonalMap, err = mtrxmult.Gen_transpose_diagonalVectors(d)
	if err != nil {
		panic(err)
	}
	TransposeLT := ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, TransposeDiagonalMap, MaxLevel, params.DefaultScale(), N1Trans, (d - 1), params.LogSlots())

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
			DSigmaLTs1[i] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DSigmaDiagonalMaps1[i], MaxLevel, params.DefaultScale(), DSigma_BSGSN1, 1, params.LogSlots())
		} else {
			DSigmaLTs1[i] = ckks.GenLinearTransform(encoder, DSigmaDiagonalMaps1[i], MaxLevel, params.DefaultScale(), params.LogSlots())
		}
	}
	for i := range DSigmaLTs2 {
		if i == len(DSigmaLTs2)-1 {
			DSigmaLTs2[i] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DSigmaDiagonalMaps2[i], MaxLevel, params.DefaultScale(), DSigma_BSGSN1, 1, params.LogSlots())
		} else {
			DSigmaLTs2[i] = ckks.GenLinearTransform(encoder, DSigmaDiagonalMaps2[i], MaxLevel, params.DefaultScale(), params.LogSlots())
		}
	}

	// create the Decomposed Maps of Tau LinearTransformation
	fmt.Printf("create the Decomposed Maps of Tau LinearTransformation\n")
	DTauLTs := make([]ckks.LinearTransform, len(DTauDiagonalMaps))
	for i := range DTauLTs {
		if i == len(DTauLTs)-1 {
			DTauLTs[i] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DTauDiagonalMaps[i], MaxLevel, params.DefaultScale(), DTau_BSGSN1, d, params.LogSlots())
		} else {
			DTauLTs[i] = ckks.GenLinearTransform(encoder, DTauDiagonalMaps[i], MaxLevel, params.DefaultScale(), params.LogSlots())
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
		LT := ckks.GenLinearTransform(encoder, U, MaxLevel, params.DefaultScale(), params.LogSlots())
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

	if loadPCAkeysFromDisk {
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

		kgen = ckks.NewKeyGenerator(params)

		// Create Rotation keys correspond to Transposition, MatrixMult and RowTotalSum
		fmt.Printf("create Rotation keys correspond to Transposition, MatrixMult and RowTotalSum\n")
		fmt.Printf("creating Rotation keys for Transposition\n")
		Rotation4Transpose := TransposeLT.Rotations4ArithmeticSeq()
		galk4Transpose = kgen.GenRotationKeysForRotations(Rotation4Transpose, false, sk)
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
		galk4Transpose = nil
		debug.FreeOSMemory()

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
		galk4MtrxMult = nil
		debug.FreeOSMemory()

		fmt.Printf("creating Rotation keys for RowtotalSum\n")
		Rotation4RowTotalSum := make([]int, d)
		for i := d; i < d*d; i = (i << 1) {
			Rotation4RowTotalSum = append(Rotation4RowTotalSum, i)
		}
		// Rotation4RowTotalSum = append(Rotation4RowTotalSum, -d+d*d)
		galk4RowtotalSum = kgen.GenRotationKeysForRotations(Rotation4RowTotalSum, false, sk)
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
	ctORSBT := make([]*rlwe.Ciphertext, int(math.Ceil(float64(rows)/float64(matrixRows))))                  // ciphetexts for One Row of Subblock Transposed
	ctORSB := make([]*rlwe.Ciphertext, int(math.Ceil(float64(rows)/float64(matrixRows))))                   // ciphetexts for One Row of Subblock.
	var ctORSBT_SigmaDecomp = make([][]*rlwe.Ciphertext, int(math.Ceil(float64(rows)/float64(matrixRows)))) // Sigma Decomposed ciphetexts for One Row of Subblock Transposed.
	// ctSubRbuff := make([]*rlwe.Ciphertext, maxSubRoutes)                                   // subroutine buffer.
	// ctOCSB := make([]*rlwe.Ciphertext, int(math.Ceil(float64(cols)/float64(matrixCols))))
	var ctOCSBSum []*rlwe.Ciphertext
	// ctOCSBSum := make([]*rlwe.Ciphertext, int(math.Ceil(float64(cols)/float64(matrixCols))))
	IdentityVec := make([]float64, 1<<params.LogSlots())
	for i := 0; i < len(IdentityVec); i++ {
		IdentityVec[i] = 0.5
	}
	// ptIdentity := encoder.EncodeNew(IdentityVec, MaxLevel, params.DefaultScale(), params.LogSlots())

	// Need to Store some of the ciphertexts into disk.
	var Covtracker = auxio.NewTracker(ctCov_path, -1)
	var Meantracker = auxio.NewTracker(ctMean_path, -1)
	// var Datatracker = auxio.NewTracker("ctData", -1)

	// Start Computing X*X^T and Mean
	fmt.Printf("Start Computing X*X^T and Mean\n")
	var elapsed time.Duration
	var now time.Time
	var wg sync.WaitGroup
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
				// //auxio.Quick_check_matrix_full(params, sk, ctORSBT[j/matrixRows], d, d)
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
		ctMean[(k-1)/matrixCols] = encryptor.EncryptZeroNew(MaxLevel) // dbg testing
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
			//auxio.Quick_check_matrix(params, sk, ctORSBT[i], d, d)
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
				ctOCSB := make([]*rlwe.Ciphertext, (k-1)/matrixCols+1)
				var wg_subroutine sync.WaitGroup
				// Retrieve OCSB of one row
				subrouteNum_inner := ((k-1)/matrixCols + 1)
				for i := 1; i <= k; i += matrixCols {
					fmt.Printf("Retrieve dataset[%d][%d], and load it into ctOCSB\n", idx, (i-1)/matrixCols)
					wg_subroutine.Add(1)
					subrouteNum_inner--
					go func(d int, i int, j int) {
						defer wg_subroutine.Done()
						var sb [][]string
						var sbf64 [][]float64
						var MatrixSB []float64
						var ptMtrxSB *rlwe.Plaintext
						encryptor_subroute := ckks.NewEncryptor(params, pk)
						encoder_subroute := ckks.NewEncoder(params)
						// evaluator_subroute := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
						if j+matrixRows > rows+1 {
							sb, err = auxio.GetOneCSVSubblock(csvfilepath, j, i, rows-j+1, matrixCols)
						} else {
							sb, err = auxio.GetOneCSVSubblock(csvfilepath, j, i, matrixRows, matrixCols)
						}
						sbf64, err = auxio.Switch2d_str2f64(sb)
						if err != nil {
							panic(err)
						}
						MatrixSB, err = mtrxmult.Row_orderingInvZeroPad(sbf64, d)
						ptMtrxSB = encoder_subroute.EncodeNew(MatrixSB, inputLevel, params.DefaultScale(), params.LogSlots())
						ctOCSB[(i-1)/matrixCols] = encryptor_subroute.EncryptNew(ptMtrxSB)
					}(d, i, j)
					if subrouteNum_inner <= 0 {
						wg_subroutine.Wait()
						subrouteNum_inner = ((k-1)/matrixCols + 1)
					}
				}
				wg_subroutine.Wait()

				// Compute MtrxMult between element in ctORSBT and elements in ctOCSB:
				fmt.Printf("Compute MtrxMult between ctORSBT[%d] and elements in %d-th ctOCSB, and add them to the previous result in ctOCSBSum:\n", idx, idx)
				ctOCSBTemp, err := mtrxmult.SquareMatrix_MultBSGS_DecomposeVer3_NeedLTs_dbg(params, sk, rlk, galk4MtrxMult, ctORSBT_SigmaDecomp[idx], ctOCSB, DTauLTs, d, len(ctOCSB))
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
						////auxio.Quick_check_matrix(params, sk, ctOCSBSum[i], d, d)
					}
				} else {
					fmt.Printf("Adding ctOCSBTemp[%d] into ctOCBSum\n", idx)
					for i := range ctOCSBTemp {
						//fmt.Printf("New value that will be added to ctOCSBSum[%d]:\n", i)
						////auxio.Quick_check_matrix(params, sk, ctOCSBTemp[i], d, d)
						evaluator.Add(ctOCSBTemp[i], ctOCSBSum[i], ctOCSBSum[i])
						//fmt.Printf("After adding a new one, ctOCSBSum[%d] value:\n", i)
						////auxio.Quick_check_matrix(params, sk, ctOCSBSum[i], d, d)
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
			// //auxio.Quick_check_matrix(params, sk, ctOCSBSum[i], d, d)
			// ctOCSBSum[i] = encryptor.EncryptNew(encoder.EncodeNew(auxio.DecryptDecode(params, sk, ctOCSBSum[i]), MaxLevel, params.DefaultScale(), params.LogSlots()))
			evaluator.MultByConst(ctOCSBSum[i], 1.0/float64(rows), ctOCSBSum[i])
			err = evaluator.Rescale(ctOCSBSum[i], params.DefaultScale(), ctOCSBSum[i])
			if err != nil {
				panic(err)
			}
			fmt.Printf("ctCov[%d][%d] now has scale %f, level %d\n", (k-1)/matrixCols, i, math.Log2(ctOCSBSum[i].Scale.Float64()), ctOCSBSum[i].Level())
			ctCov[(k-1)/matrixCols][i] = ctOCSBSum[i].CopyNew()
			//auxio.Quick_check_matrix(params, sk, ctCov[(k-1)/matrixCols][i], d, d)
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
		// //auxio.Quick_check_matrix_full(params, sk, ctMean[(k-1)/matrixCols], d, d) // dbg testing
		evaluator.MultByConst(ctMean[(k-1)/matrixCols], 1.0/float64(rows), ctMean[(k-1)/matrixCols])
		fmt.Printf("Check the MeanVec\n")
		//auxio.Quick_check_matrix(params, sk, ctMean[(k-1)/matrixCols], d, d) // dbg testing
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
				//auxio.Quick_check_matrix(params, sk, ctCov[j][i], d, d)
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
			//auxio.Quick_check_matrix(params, sk, ctCov[i][j], d, d)
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

func CovMatrix_btpVersion_check(btpparamsOn bool, reencryptionMode bool) {
	var err error
	matrixCols := 128
	// matrixRows := 128
	// rows := 256
	cols := 257

	// Size of the Square Matrix.
	d := 1 << 7

	// N1
	N1Trans := 16

	// Target Number of EigenVectors.
	TargetEigVecNum := 4

	// necessary file paths:
	var KeysPath_btpVersion string
	var Rtk4Transpose_path string
	var Rtk4RowtotalSum_path string
	var Rtk4ColtotalSum_path string
	var ctCov_path string
	var ctMean_path string
	if btpparamsOn {

		KeysPath_btpVersion = "Keys_btpVer"
		Rtk4Transpose_path = "Rtk4Transpose_btpVer"
		Rtk4RowtotalSum_path = "Rtk4RowtotalSum_btpVer"
		Rtk4ColtotalSum_path = "Rtk4ColtotalSum_btpVer"
		ctCov_path = "ctCov_btpVer"
		ctMean_path = "ctMean_btpVer"
	} else {
		KeysPath_btpVersion = "Keys"
		Rtk4Transpose_path = "Rtk4Transpose"
		Rtk4RowtotalSum_path = "Rtk4RowtotalSum"
		Rtk4ColtotalSum_path = "Rtk4ColtotalSum"
		ctCov_path = "ctCov"
		ctMean_path = "ctMean"
	}

	// initialize encryption scheme params.

	//var ckksParams ckks.ParametersLiteral
	var btpParams bootstrapping.Parameters
	var params ckks.Parameters
	var sk *rlwe.SecretKey
	var pk *rlwe.PublicKey
	var rlk *rlwe.RelinearizationKey
	var btpevk bootstrapping.EvaluationKeys
	var bootstrapper *bootstrapping.Bootstrapper
	// var kgen rlwe.KeyGenerator

	if btpparamsOn {
		if reencryptionMode {
			params, btpParams, sk, pk, btpevk.Rlk, btpevk.Rtks, btpevk.SwkDtS, btpevk.SwkStD, err = LoadKeys_btpVersion(KeysPath_btpVersion)
			if err != nil {
				panic(err)
			}
			// let rlk be a copy of btpevk.Rlk, they share the same memory address.
			rlk = btpevk.Rlk // History lefting problems... rlk should be deleted, but...

			// bootstrapper
			if bootstrapper, err = bootstrapping.NewBootstrapper(params, btpParams, btpevk); err != nil {
				panic(err)
			}
		} else {
			params, btpParams, sk, pk, rlk, _, _, _, err = LoadKeys_btpVersion(KeysPath_btpVersion)
		}
	} else {
		params, sk, pk, rlk, err = LoadKeys_check()
		if err != nil {
			panic(err)
		}
	}

	// Encryptor
	encryptor := ckks.NewEncryptor(params, pk)

	// Decryptor
	decryptor := ckks.NewDecryptor(params, sk)

	// encoder
	encoder := ckks.NewEncoder(params)

	// In the Boostrapping Version, we do not begin the performance of operations from the top level,
	// instead we will begin at the first level which is not included in bootstrap procedure, the number of this level can be
	// computed by:
	var MaxLevel int
	if btpparamsOn {
		MaxLevel = btpParams.SlotsToCoeffsParameters.LevelStart - len(btpParams.SlotsToCoeffsParameters.ScalingFactor)
	} else {
		MaxLevel = params.MaxLevel()
	}

	// create the Map and LinearTransform object of Transpose LinearTransformation
	fmt.Printf("create the Map and LinearTransform object of Transpose LinearTransformation\n")
	var TransposeDiagonalMap map[int][]float64
	TransposeDiagonalMap, err = mtrxmult.Gen_transpose_diagonalVectors(d)
	if err != nil {
		panic(err)
	}
	TransposeLT := ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, TransposeDiagonalMap, MaxLevel, params.DefaultScale(), N1Trans, (d - 1), params.LogSlots())

	// create FileTrackers
	var Meantracker *auxio.Filetracker
	var Covtracker *auxio.Filetracker
	var RtkTracker4Transpose *auxio.Filetracker
	var RtkTracker4RowtotalSum *auxio.Filetracker
	var galk4Transpose *rlwe.RotationKeySet
	var galk4RowtotalSum *rlwe.RotationKeySet
	var galk4ColtotalSum *rlwe.RotationKeySet
	var RtkTracker4ColtotalSum *auxio.Filetracker

	// Construct ((cols-1)/matrixCols)^2 Ciphertexts representing a cols x cols Covariance Matrix
	ctCov := make([][]*rlwe.Ciphertext, (cols-1)/matrixCols)
	ctMean := make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	ctMMT := make([][]*rlwe.Ciphertext, (cols-1)/matrixCols) // ciphertext for MeanVec · MeanVec^T
	ctEigVec := make([][]*rlwe.Ciphertext, TargetEigVecNum)
	ctEigVal := make([]*rlwe.Ciphertext, TargetEigVecNum)
	for i := 0; i < (cols-1)/matrixCols; i++ {
		ctCov[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
		ctMMT[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	}
	for i := 0; i < len(ctEigVec); i++ {
		ctEigVec[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	}
	var ctTemp *rlwe.Ciphertext
	IdentityVec := make([]float64, 1<<params.LogSlots())
	for i := 0; i < len(IdentityVec); i++ {
		IdentityVec[i] = 0.5
	}
	// ptIdentity := encoder.EncodeNew(IdentityVec, params.MaxLevel(), params.DefaultScale(), params.LogSlots())

	// Compute Mean*Mean^T, first retrieving the Mean from disk.
	//_, err = auxio.ReadCtVec4file(Meantracker.Fp, Meantracker.Hierarchy[0], ctMean, 0)
	Meantracker = auxio.NewTracker4File(ctMean_path)
	RtkTracker4Transpose = auxio.NewTracker4File(Rtk4Transpose_path)
	RtkTracker4RowtotalSum = auxio.NewTracker4File(Rtk4RowtotalSum_path)
	// RtkTracker4ColtotalSum = auxio.NewTracker4File(Rtk4ColtotalSum_path)

	galk4RowtotalSum = new(rlwe.RotationKeySet)
	galk4ColtotalSum = new(rlwe.RotationKeySet)
	RtkTracker4RowtotalSum.ReadUpdateOne(galk4RowtotalSum)
	//RtkTracker4ColtotalSum.ReadUpdateOne(galk4ColtotalSum)

	// Create RotationKeys for ColumnTotalSum
	kgen := ckks.NewKeyGenerator(params)
	RtkTracker4ColtotalSum = auxio.NewTracker(Rtk4ColtotalSum_path, -1)
	Rots4ColtotalSum := make([]int, 0)
	for i := 1; i < d; i = (i << 1) {
		Rots4ColtotalSum = append(Rots4ColtotalSum, i) // wairing to check
		Rots4ColtotalSum = append(Rots4ColtotalSum, -i)
	}
	galk4ColtotalSum = kgen.GenRotationKeysForRotations(Rots4ColtotalSum, false, sk)
	_, err = RtkTracker4ColtotalSum.StoreUpdateOne(galk4ColtotalSum)
	if err != nil {
		panic(err)
	}
	_, err = RtkTracker4ColtotalSum.StoreFinish()
	if err != nil {
		panic(err)
	}

	// timer
	var now time.Time
	var elapsed time.Duration
	cnt := 0
	var btpnow time.Time
	var btpelapsed time.Duration

	// Retrieve the galk4Transpose
	fmt.Printf("Retrieve the galk4Transpose\n")
	galk4Transpose = new(rlwe.RotationKeySet)
	RtkTracker4Transpose.ReadUpdateOne(galk4Transpose)

	// evaluator
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk4Transpose})

	for i := 0; i < len(ctMean); i++ {
		ctMean[i] = new(rlwe.Ciphertext)
		Meantracker.ReadUpdateOne(ctMean[i])
		evaluator.Rescale(ctMean[i], params.DefaultScale(), ctMean[i])
		fmt.Printf("The %d-th ctMean with scale %f, level %d: \n", i, math.Log2(ctMean[i].Scale.Float64()), ctMean[i].Level())
	}

	fmt.Printf("Compute the Transpose of Mean vector\n")
	now = time.Now()
	for i := 0; i < len(ctMean); i++ {

		ctTemp = evaluator.LinearTransform4ArithmeticSeqNew(ctMean[i], TransposeLT)[0]
		fmt.Printf("The %d-th ctMean^T with scale %f, level %d: \n", i, math.Log2(ctTemp.Scale.Float64()), ctTemp.Level())
		//auxio.Quick_check_matrix(params, sk, ctTemp, d, d)
		if err != nil {
			panic(err)
		}
		for j := 0; j < len(ctMean); j++ {
			// fmt.Printf("The i:%d, j:%d ctMMT: \n", i, j)
			ctMMT[i][j] = evaluator.MulRelinNew(ctTemp, ctMean[j])
			evaluator.Rescale(ctMMT[i][j], params.DefaultScale(), ctMMT[i][j])
			fmt.Printf("The ctMMT[%d][%d] with scale %f, level %d: \n", i, j, math.Log2(ctMMT[i][j].Scale.Float64()), ctMMT[i][j].Level())
		}
	}
	fmt.Printf("Transposition of Mean vector complete in %s\n", time.Since(now))
	elapsed += time.Since(now)
	// free the Transpose Rotation keys:
	galk4Transpose = nil
	debug.FreeOSMemory()

	// Retrieve the X*X^T from disk, and do Cov = X*X^T-Mean*Mean^T
	fmt.Printf("Retrieve the X*X^T from disk, and do Cov = X*X^T-Mean*Mean^T\n")
	now = time.Now()
	Covtracker = auxio.NewTracker4File(ctCov_path)
	for i := 0; i < len(ctCov); i++ {
		//_, err = auxio.ReadCtVec4file(Covtracker.Fp, Covtracker.Hierarchy[0], ctCov[i], i*Covtracker.Hierarchy[0]*Covtracker.Hierarchy[1])
		for j := 0; j < len(ctCov); j++ {
			fmt.Printf("The i:%d, j:%d ctMM^T with scale %f, level %d: \n", i, j, math.Log2(ctMMT[i][j].Scale.Float64()), ctMMT[i][j].Level())
			//auxio.Quick_check_matrix(params, sk, ctMMT[i][j], d, d)
			ctCov[i][j] = new(rlwe.Ciphertext)
			Covtracker.ReadUpdateOne(ctCov[i][j])
			fmt.Printf("The Original i:%d, j:%d ctX^TX has scale %f, level %d: \n", i, j, math.Log2(ctCov[i][j].Scale.Float64()), ctCov[i][j].Level())

			// Synchronizing scale Method 1:
			/*
				for k := 0; k <= ctCov[i][j].Level()-ctMMT[i][j].Level(); k++ {
					evaluator.Mul(ctCov[i][j], ptIdentity, ctCov[i][j])
					evaluator.MultByConst(ctCov[i][j], 2, ctCov[i][j])
					err = evaluator.Rescale(ctCov[i][j], params.DefaultScale(), ctCov[i][j])
					if err != nil {
						panic(err)
					}
					fmt.Printf("Drop ctX^TX to Scale %f, level %d\n", math.Log2(ctCov[i][j].Scale.Float64()), ctCov[i][j].Level())
				}
			*/
			// Synchronizing scale Method 2:
			/*
				err = evaluator.Rescale(ctCov[i][j], params.DefaultScale(), ctCov[i][j])
				if err != nil {
					panic(err)
				}
			*/
			if ctMMT[i][j].Level() >= ctCov[i][j].Level() {
				evaluator.SetScale(ctMMT[i][j], ctCov[i][j].Scale)
			} else {
				evaluator.SetScale(ctCov[i][j], ctMMT[i][j].Scale)
			}

			// evaluator.Rescale(ctCov[i][j], params.DefaultScale(), ctCov[i][j])
			fmt.Printf("The i:%d, j:%d ctX^TX with scale %f, level %d: \n", i, j, math.Log2(ctCov[i][j].Scale.Float64()), ctCov[i][j].Level())
			//auxio.Quick_check_matrix(params, sk, ctCov[i][j], d, d)
			evaluator.Sub(ctCov[i][j], ctMMT[i][j], ctCov[i][j])
			fmt.Printf("The i:%d, j:%d ctCov: with scale %f, level %d: \n", i, j, math.Log2(ctCov[i][j].Scale.Float64()), ctCov[i][j].Level())
			// //auxio.Quick_check_matrix_full(params, sk, ctCov[i][j], d, d)
			// do one recryption,this is for the following PowerMethod.
			// ctCov[i][j] = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctCov[i][j]), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
			//auxio.Quick_check_matrix(params, sk, ctCov[i][j], d, d)

			if ctCov[i][j].Level() < 6 {
				if reencryptionMode {
					// ctCov[i][j] = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctCov[i][j]), params.LogSlots()), MaxLevel, params.DefaultScale(), params.LogSlots()))
					btpnow = time.Now()
					fmt.Printf("Before Bootstrap, ctCov[%d][%d] :\n", i, j)
					//auxio.Quick_check_matrix(params, sk, ctCov[i][j], d, d)

					/*
						fmt.Printf("Debug scenario: \n")
						ctTemp := bootstrapper.Bootstrap_dbg(ctCov[i][j])
						//auxio.Quick_check_infos(ctTemp, "ctTemp original")
						//auxio.Quick_check_matrix(params, sk, ctTemp, d, d)
						ctCov[i][j] = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctCov[i][j]), params.LogSlots()), MaxLevel, params.DefaultScale(), params.LogSlots()))
						ctTemp = bootstrapper.Bootstrap_dbg(ctCov[i][j])
						//auxio.Quick_check_infos(ctTemp, "ctTemp after manual re-encrypt")
						//auxio.Quick_check_matrix(params, sk, ctTemp, d, d)
					*/

					ctCov[i][j] = bootstrapper.Bootstrap(ctCov[i][j])
					fmt.Printf("After Bootstrap in %s, ctCov[%d][%d]:\n", time.Since(btpnow), i, j)
					//auxio.Quick_check_matrix(params, sk, ctCov[i][j], d, d)
					elapsed += time.Since(btpnow)
				} else {
					ctCov[i][j] = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctCov[i][j]), params.LogSlots()), MaxLevel, params.DefaultScale(), params.LogSlots()))
				}
				fmt.Printf("ctCov is brought back to level %d\n", ctCov[i][j].Level())
				cnt++
			}
		}
	}
	elapsed += time.Since(now)
	fmt.Printf("Cov complete in %s\n", time.Since(now))

	// Compute the k most dominant eigenVectors and the corresponding eigenValues, We will skip the Boostrapping procedure, and use recryption instead.
	// create a Vi
	ptVi := encoder.EncodeNew(generateRandomArray(params.Slots()), MaxLevel, params.DefaultScale(), params.LogSlots())
	// ptVi := auxio.Encode_single_float64(params, 1.0, MaxLevel, params.DefaultScale())
	ctVi := make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	ctVip1 := make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	ctViT := make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	/*
		for i := 0; i < len(ctVi); i++ {
			ctVi[i] = encryptor.EncryptNew(ptVi)
		}
	*/

	// PowerMethod Iteration, iteratively compute: Vi <- Vi^T * Cov , notice that Cov is a semmetric matrix.
	var ctSum *rlwe.Ciphertext
	// var ctSign *rlwe.Ciphertext

	IterNum := 15
	// NewtonIterNum := 12
	VecMod := 1 // We will have two Vector Mod: ColVector Mod (1) & RowVector Mod (0)
	var InvSRTiters = [14]int{7, 14, 13, 14, 13, 14, 13, 14, 14, 13, 14, 13, 14, 23}
	if IterNum-1 == 4 {
		InvSRTiters[3] = 24
	}

	fmt.Printf("Entering PowerMethod\n")
	now = time.Now()
	for k := 0; k < TargetEigVecNum; k++ {
		VecMod = 1
		for i := 0; i < len(ctVi); i++ {
			vec := make([]float64, params.Slots())
			random_arr := generateRandomArray(matrixCols)
			if VecMod == 1 {
				for j := range random_arr {
					for m := 0; m < matrixCols; m++ {
						vec[j*matrixCols+m] = random_arr[j]
					}
				}
			} else {
				for j := 0; i < params.Slots()/matrixCols; i++ {
					copy(vec[j*matrixCols:(j+1)*matrixCols], random_arr)
				}
			}
			ptVi = encoder.EncodeNew(vec, MaxLevel, params.DefaultScale(), params.LogSlots())
			ctVi[i] = encryptor.EncryptNew(ptVi)
		}
		for t := 0; t < IterNum; t++ {
			fmt.Printf("------------------ this is the %dth iteration ---------------\n", t)
			fmt.Printf("Begin Cov LT\n")
			// Check if a re-encryption is needed.
			if ctVi[0].Level() < 4 || t == IterNum-1 {
				var ctTemp *rlwe.Ciphertext
				if VecMod == 0 {
					/*
						fmt.Printf("Before Combination:\n")
						for _, ct := range ctVi {
							//auxio.Quick_check_matrix(params, sk, ct, d, d)
						}
					*/
					ctTemp, err = mtrxmult.ReplicatedVec_Combine_dbg(params, sk, ctVi, d, VecMod)
					if err != nil {
						panic(err)
					}
					// Decide bootstrap or re-encrypt
					if reencryptionMode {
						btpnow = time.Now()
						fmt.Printf("Before Bootstrap, ctTemp:\n")
						//auxio.Quick_check_matrix(params, sk, ctTemp, d, d)
						ctTemp = bootstrapper.Bootstrap(ctTemp)
						fmt.Printf("After Bootstrap in %s, ctTemp:\n", time.Since(btpnow))
						//auxio.Quick_check_matrix(params, sk, ctTemp, d, d)
						elapsed += time.Since(btpnow)
					} else {
						ctTemp = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctTemp), params.LogSlots()), MaxLevel, params.DefaultScale(), params.LogSlots()))
					}
					ctVi, err = mtrxmult.ReplicateVec_Decompose_dbg(params, sk, galk4RowtotalSum, ctTemp, d, VecMod, len(ctVi))
					/*
						fmt.Printf("After Combination:\n")
						for _, ct := range ctVi {
							//auxio.Quick_check_matrix(params, sk, ct, d, d)
						}
					*/
					if err != nil {
						panic(err)
					}
					fmt.Printf("ctVip1 is previously at level %d, we will bring it back to level %d\n", ctTemp.Level(), ctVip1[0].Level())
					cnt++
				} else if VecMod == 1 && (ctVi[0].Level() < 3 || t == IterNum-1) {
					ctTemp, err = mtrxmult.ReplicatedVec_Combine_dbg(params, sk, ctVi, d, VecMod)
					if err != nil {
						panic(err)
					}
					// Decide bootstrap or re-encrypt
					if reencryptionMode {
						btpnow = time.Now()
						fmt.Printf("Before Bootstrap, ctTemp:\n")
						//auxio.Quick_check_matrix(params, sk, ctTemp, d, d)
						ctTemp = bootstrapper.Bootstrap(ctTemp)
						fmt.Printf("After Bootstrap in %s, ctTemp:\n", time.Since(btpnow))
						//auxio.Quick_check_matrix(params, sk, ctTemp, d, d)
						elapsed += time.Since(btpnow)
					} else {
						ctTemp = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctTemp), params.LogSlots()), MaxLevel, params.DefaultScale(), params.LogSlots()))
					}
					ctVi, err = mtrxmult.ReplicateVec_Decompose_dbg(params, sk, galk4ColtotalSum, ctTemp, d, VecMod, len(ctVi))
					if err != nil {
						panic(err)
					}
					fmt.Printf("ctVip1 is previously at level %d, we will bring it back to level %d\n", ctTemp.Level(), ctVip1[0].Level())
					cnt++
				}
			}
			if VecMod == 0 {
				for i := 0; i < len(ctCov); i++ {
					// We will do ctCov[j][i] * ctVi[j] , but will actually do ColAggregate(\sum(ctCov[i][j] * ctVi[j])) to complete this task, where ctCov[i][j]^T = ctCov[j][i]
					fmt.Printf("ctVi[%d]:\n", 0)
					//auxio.Quick_check_infos(ctVi[0], "ctVi")
					ctSum = evaluator.MulNew(ctVi[0], ctCov[i][0]) // Consume one level
					fmt.Printf("ctCov[%d][%d]:\n", i, 0)
					//auxio.Quick_check_infos(ctCov[i][0], "ctCov")
					// //auxio.Quick_check_matrix(params, sk, ctCov[i][0], d, d)
					// //auxio.Quick_check_matrix(params, sk, ctVi[0], d, d)
					fmt.Printf("MulAndAddResult:\n")
					//auxio.Quick_check_infos(ctSum, "ctSum")
					// //auxio.Quick_check_matrix(params, sk, ctSum, d, d)
					for j := 1; j < len(ctCov); j++ {
						fmt.Printf("ctCov[%d][%d]:\n", i, j)
						////auxio.Quick_check_matrix(params, sk, ctCov[i][j], d, d)
						fmt.Printf("ctVi[%d]:\n", j)
						// //auxio.Quick_check_matrix(params, sk, ctVi[j], d, d)
						evaluator.MulAndAdd(ctVi[j], ctCov[i][j], ctSum)
						fmt.Printf("MulAndAddResult:\n")
						////auxio.Quick_check_matrix(params, sk, ctSum, d, d)
					}
					evaluator.Relinearize(ctSum, ctSum)
					//fmt.Printf("Before Aggregation:\n")
					////auxio.Quick_check_matrix(params, sk, ctSum, d, d)
					err = evaluator.Rescale(ctSum, params.DefaultScale(), ctSum)
					if err != nil {
						panic(err)
					}
					ctSum, err = mtrxmult.Matrix_Aggregate_withRowOrder(params, sk, galk4ColtotalSum, ctSum, d, d, 1) // Consume one level
					if err != nil {
						panic(err)
					}
					fmt.Printf("The %dth ctSum Aggregated with scale %f level %d: \n", i, math.Log2(ctSum.Scale.Float64()), ctSum.Level())
					////auxio.Quick_check_matrix_full(params, sk, ctSum, d, d)

					ctVip1[i] = ctSum.CopyNew()
					fmt.Printf("The %dth Vip1 with scale %f level %d: \n", i, math.Log2(ctVip1[i].Scale.Float64()), ctVip1[i].Level())
					////auxio.Quick_check_matrix_full(params, sk, ctVip1[i], d, d)
					ctSum = nil

				}
				// now the ctVi has become the ColVector mod, we switch the sign:
				VecMod = 1

			} else if VecMod == 1 {
				for i := 0; i < len(ctCov); i++ {
					// we will do ctCov[i][j] * ctVi[j], but will actually do RowAggregate(\sum(ctCov[j][i] * ctVi[j])) to complete this task, where ctCov[i][j]^T = ctCov[j][i]
					ctSum = evaluator.MulNew(ctVi[0], ctCov[0][i]) // Consume one level
					fmt.Printf("ctCov[%d][%d]:\n", 0, i)
					////auxio.Quick_check_matrix(params, sk, ctCov[0][i], d, d)
					fmt.Printf("ctVi[%d]:\n", 0)
					////auxio.Quick_check_matrix(params, sk, ctVi[0], d, d)
					fmt.Printf("MulAndAddResult:\n")
					////auxio.Quick_check_matrix(params, sk, ctSum, d, d)
					for j := 1; j < len(ctCov); j++ {
						fmt.Printf("ctCov[%d][%d]:\n", j, i)
						////auxio.Quick_check_matrix(params, sk, ctCov[j][i], d, d)
						fmt.Printf("ctVi[%d]:\n", j)
						////auxio.Quick_check_matrix(params, sk, ctVi[j], d, d)
						evaluator.MulAndAdd(ctVi[j], ctCov[j][i], ctSum)
						fmt.Printf("MulAndAddResult:\n")
						////auxio.Quick_check_matrix(params, sk, ctSum, d, d)
					}
					evaluator.Relinearize(ctSum, ctSum)
					//fmt.Printf("Before Aggregation:\n")
					////auxio.Quick_check_matrix(params, sk, ctSum, d, d)
					err = evaluator.Rescale(ctSum, params.DefaultScale(), ctSum)
					if err != nil {
						panic(err)
					}
					ctSum, err = mtrxmult.Matrix_Aggregate_withRowOrder(params, sk, galk4RowtotalSum, ctSum, d, d, 0)
					if err != nil {
						panic(err)
					}
					fmt.Printf("The %dth ctSum Aggregated with scale %f level %d: \n", i, math.Log2(ctSum.Scale.Float64()), ctSum.Level())
					//auxio.Quick_check_matrix(params, sk, ctSum, d, d)

					ctVip1[i] = ctSum.CopyNew()
					fmt.Printf("The %dth Vip1 with scale %f level %d: \n", i, math.Log2(ctVip1[i].Scale.Float64()), ctVip1[i].Level())
					//auxio.Quick_check_matrix(params, sk, ctVip1[i], d, d)
					ctSum = nil
					// ctSum = encryptor.EncryptZeroNew(params.MaxLevel())
				}
				// now the ctVi has become the RowVecotr mod, we switch the sign:
				VecMod = 0
			}
			fmt.Printf("Cov LT complete in %s\n", time.Since(now))

			// if this is not the last turn, we will do the normalisation. The last turn is only for computing eigen value and do not update eigen vector.
			if t < IterNum-1 {
				// Computing Inner product.
				ctSum = evaluator.MulNew(ctVip1[0], ctVip1[0])
				for i := 1; i < len(ctVip1); i++ { // Sum ctVip1
					evaluator.MulAndAdd(ctVip1[i], ctVip1[i], ctSum)
				}
				ctSum = evaluator.RelinearizeNew(ctSum)
				err = evaluator.Rescale(ctSum, params.DefaultScale(), ctSum) // Consume one level
				if err != nil {
					panic(err)
				}

				// Aggregate ctSum
				fmt.Printf("Before Inner Product:\n")
				//auxio.Quick_check_matrix(params, sk, ctSum, d, d)
				// VecMod Now represents the ctVip1's form.
				if VecMod == 1 {
					ctSum, err = mtrxmult.Matrix_Aggregate_withRowOrder(params, sk, galk4RowtotalSum, ctSum, d, d, 0)
				} else if VecMod == 0 {
					ctSum, err = mtrxmult.Matrix_Aggregate_withRowOrder(params, sk, galk4ColtotalSum, ctSum, d, d, 1) // Consume one level
				}
				if err != nil {
					panic(err)
				}
				fmt.Printf("Inner Product with Scale %f, Level %d:\n", math.Log2(ctSum.Scale.Float64()), ctSum.Level())
				//auxio.Quick_check_matrix(params, sk, ctSum, d, d)

				// Taylor Init Guess:
				if ctSum.Level() <= 2 {
					fmt.Printf("ctSum is now at level %d, we will bring it back to level %d\n", ctSum.Level(), MaxLevel)
					if reencryptionMode {
						btpnow = time.Now()
						fmt.Printf("Before Bootstrap, ctSum:\n")
						//auxio.Quick_check_matrix(params, sk, ctSum, d, d)
						ctSum = bootstrapper.Bootstrap(ctSum)
						fmt.Printf("After Bootstrap in %s, ctSum:\n", time.Since(btpnow))
						//auxio.Quick_check_matrix(params, sk, ctSum, d, d)
						elapsed += time.Since(btpnow)
					} else {
						ctSum = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctSum), params.LogSlots()), MaxLevel, params.DefaultScale(), params.LogSlots()))
					}
					cnt++
				}
				ctGuess, err := nonpolyfunc.TaylorInitNew(params, rlk, ctSum, nonpolyfunc.Inv_sqrt_taypor1_0to2pow16[:]) // Consume one level *2
				if err != nil {
					panic(err)
				}
				fmt.Printf("Taylor Guess with Scale %f, Level %d:\n", math.Log2(ctGuess.Scale.Float64()), ctGuess.Level())
				//auxio.Quick_check_matrix(params, sk, ctGuess, d, d)

				// Do one Recryption:
				// ctGuess = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctGuess), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
				// ctSum = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctSum), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))

				// InvSqrt by Newton:
				var ctInvsqrt *rlwe.Ciphertext
				var btptimes int
				//var NewtonIterNumtmp = NewtonIterNum
				/*
					if t == IterNum-2 {
						NewtonIterNumtmp = NewtonIterNum + 5
					}
				*/
				var NewtonIterNumtmp = InvSRTiters[t]
				fmt.Printf("NewtonIter times: %d\n", NewtonIterNumtmp)
				if reencryptionMode {
					ctInvsqrt, btptimes, err = nonpolyfunc.InvSqrtByNewton_btpVer3_dbg(params, rlk, pk, sk, ctSum, ctGuess, NewtonIterNumtmp, MaxLevel, reencryptionMode, bootstrapper) // Consume one level * 3 * NewtonIterNum
				} else {
					ctInvsqrt, btptimes, err = nonpolyfunc.InvSqrtByNewton_btpVer3_dbg(params, rlk, pk, sk, ctSum, ctGuess, NewtonIterNumtmp, MaxLevel, reencryptionMode) // Consume one level * 3 * NewtonIterNum
				}
				if err != nil {
					panic(err)
				}
				cnt += btptimes
				fmt.Printf("InvSqrt with Scale %f, Level %d:\n", math.Log2(ctInvsqrt.Scale.Float64()), ctInvsqrt.Level())
				//auxio.Quick_check_matrix(params, sk, ctInvsqrt, d, d)

				// Do one Recryption:
				// ctInvsqrt = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctInvsqrt), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
				// fmt.Printf("InvSqrt after recryption with Scale %f, Level %d:\n", math.Log2(ctInvsqrt.Scale.Float64()), ctInvsqrt.Level())

				// Normalisation
				if ctInvsqrt.Level() <= 1 {
					fmt.Printf("ctInvsqrt is now at level %d, we will bring it back to level %d\n", ctInvsqrt.Level(), MaxLevel)
					if reencryptionMode {
						btpnow = time.Now()
						fmt.Printf("Before Bootstrap, ctInvsqrt:\n")
						//auxio.Quick_check_matrix(params, sk, ctInvsqrt, d, d)
						ctInvsqrt = bootstrapper.Bootstrap(ctInvsqrt)
						fmt.Printf("After Bootstrap in %s, ctInvsqrt:\n", time.Since(btpnow))
						//auxio.Quick_check_matrix(params, sk, ctInvsqrt, d, d)
						elapsed += time.Since(btpnow)
					} else {
						ctInvsqrt = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctInvsqrt), params.LogSlots()), MaxLevel, params.DefaultScale(), params.LogSlots()))
					}
					cnt++
				}
				if ctVip1[0].Level() <= 1 {
					var ctTemp *rlwe.Ciphertext
					ctTemp, err = mtrxmult.ReplicatedVec_Combine_dbg(params, sk, ctVip1, d, VecMod)
					if err != nil {
						panic(err)
					}
					if reencryptionMode {
						btpnow = time.Now()
						fmt.Printf("Before Bootstrap, ctTemp:\n")
						//auxio.Quick_check_matrix(params, sk, ctTemp, d, d)
						ctTemp = bootstrapper.Bootstrap(ctTemp)
						fmt.Printf("After Bootstrap in %s, ctTemp:\n", time.Since(btpnow))
						//auxio.Quick_check_matrix(params, sk, ctTemp, d, d)
						elapsed += time.Since(btpnow)
					} else {
						ctTemp = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctTemp), params.LogSlots()), MaxLevel, params.DefaultScale(), params.LogSlots()))
					}
					ctVip1, err = mtrxmult.ReplicateVec_Decompose_dbg(params, sk, galk4ColtotalSum, ctTemp, d, VecMod, len(ctVip1))
					if err != nil {
						panic(err)
					}
					cnt++
					fmt.Printf("ctVip1 is previously at level %d, we will bring it back to level %d\n", ctTemp.Level(), ctVip1[0].Level())
				}
				for i := 0; i < len(ctVip1); i++ {
					fmt.Printf("ctVip1[%d] before normalisation:\n", i)
					//auxio.Quick_check_matrix(params, sk, ctVip1[i], d, d)
					evaluator.MulRelin(ctVip1[i], ctInvsqrt, ctVip1[i])
					err = evaluator.Rescale(ctVip1[i], params.DefaultScale(), ctVip1[i]) // Consume one level
					if err != nil {
						panic(err)
					}
					fmt.Printf("ctVip1[%d] after normalisation:\n", i)
					//auxio.Quick_check_matrix(params, sk, ctVip1[i], d, d)

				}
				fmt.Printf("ctVip1 after normalisation has Scale %f, Level %d, this will be used to update ctVi\n", math.Log2(ctVip1[0].Scale.Float64()), ctVip1[0].Level())

				// Update ctVi with ctVip1
				for i := 0; i < len(ctVip1); i++ {
					ctVi[i] = ctVip1[i].CopyNew()
					fmt.Printf("The %dth iteration's Vip1 %d result:\n", t, i)
					// //auxio.Quick_check_matrix_full(params, sk, ctVi[i], d, d)
				}

			} else { // if this is the last turn, compute the eigen value corresponding to the eigen vector.
				// Recall that at this point, ctVip1 should contain Cov @ ctVi, the eigen value we want can be computed by :
				// eig_val = <(Cov @ ctVi), ctVi> / <ctVi, ctVi> = <ctVip1, ctVi> / <ctVi,ctVi>, if we have confidence in the
				// normalisation procedure in previous turn, than this equation can be reduced to eig_val = <ctVip1,ctVi> / 1
				// After this, we will have to compute the Shifted Covariance Matrix: ShiftedCov = Cov - eigval * eigvec^T @ eigvec

				fmt.Printf("This is the last iteration of %dth EigenVector computation\n", k)
				// So we will Compute the Inner product <ctVip1,ctVi>
				// VecMod Now represents the ctVip1's form.
				if ((VecMod + 1) & (2 - 1)) == 1 { // if ctVi is now in column mod, then we create a replication of its row mod version
					// this will cost one level.
					for i := 0; i < len(ctVip1); i++ {
						ctViT[i], err = mtrxmult.ReplicateVec_Switch_Axis(params, sk, galk4RowtotalSum, ctVi[i], d, (VecMod+1)&(2-1)) // Consume one level
						fmt.Printf("The ctViT with Scale %f, level %d:\n", math.Log2(ctViT[i].Scale.Float64()), ctViT[i].Level())
						//auxio.Quick_check_matrix(params, sk, ctViT[i], d, d)
						fmt.Printf("The ctVi  with Scale %f, level %d:\n", math.Log2(ctVi[i].Scale.Float64()), ctVi[i].Level())
						//auxio.Quick_check_matrix(params, sk, ctVi[i], d, d)
						if err != nil {
							panic(err)
						}
					}
				} else if ((VecMod + 1) & (2 - 1)) == 0 { // if ctVi is now in row mod, then we create a replication of its column mod version
					// this will cost one level.
					for i := 0; i < len(ctVip1); i++ {
						ctViT[i], err = mtrxmult.ReplicateVec_Switch_Axis(params, sk, galk4ColtotalSum, ctVi[i], d, (VecMod+1)&(2-1)) // Consume one level *2
						if err != nil {
							panic(err)
						}
					}
				}

				ctSum = evaluator.MulNew(ctVip1[0], ctViT[0])
				for i := 1; i < len(ctVip1); i++ { // Sum ctVip1
					evaluator.MulAndAdd(ctVip1[i], ctViT[i], ctSum)
				}
				ctSum = evaluator.RelinearizeNew(ctSum)
				err = evaluator.Rescale(ctSum, params.DefaultScale(), ctSum) // Consume one level
				if err != nil {
					panic(err)
				}

				// Aggregate ctSum then we get the eig_val.
				// fmt.Printf("Before Inner Product:\n")
				// //auxio.Quick_check_matrix(params, sk, ctSum, d, d)
				// VecMod Now represents the ctVip1's form.
				if VecMod == 1 { // since ctVip1 is in column mod, then ctSum is the multiple of ctVip1 and ctViT that is in Column mod
					ctSum, err = mtrxmult.Matrix_Aggregate_withRowOrder(params, sk, galk4RowtotalSum, ctSum, d, d, 0)
				} else if VecMod == 0 { // since ctVip1 is in Row mod, then ctSum is the multiple of ctVip1 and ctViT that is in row mod
					ctSum, err = mtrxmult.Matrix_Aggregate_withRowOrder(params, sk, galk4ColtotalSum, ctSum, d, d, 1) // Consume one level
				}
				if err != nil {
					panic(err)
				}
				fmt.Printf("Aggregate ctSum then we get the eig_val with Scale %f, level %d:\n", math.Log2(ctSum.Scale.Float64()), ctSum.Level())
				// DO one Recryption.
				//auxio.Quick_check_matrix(params, sk, ctSum, d, d)

				// Store the eig_val:
				ctEigVal[k] = ctSum.CopyNew()
				// Store the eig_Vec:
				for i := 0; i < len(ctVi); i++ {
					ctEigVec[k][i] = ctVi[i].CopyNew()
				}
				if ctEigVal[k].Level() <= ctCov[0][0].Level() {
					fmt.Printf("ctEigVal is now at level %d, we will bring it back to level %d\n", ctEigVal[k].Level(), MaxLevel)
					if reencryptionMode {
						btpnow = time.Now()
						fmt.Printf("Before Bootstrap, ctEigVal[k]:\n")
						//auxio.Quick_check_matrix(params, sk, ctEigVal[k], d, d)
						ctEigVal[k] = bootstrapper.Bootstrap(ctEigVal[k])
						fmt.Printf("After Bootstrap in %s, ctEigVal[k]:\n", time.Since(btpnow))
						//auxio.Quick_check_matrix(params, sk, ctEigVal[k], d, d)
						elapsed += time.Since(btpnow)
					} else {
						ctEigVal[k] = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctEigVal[k]), params.LogSlots()), MaxLevel, params.DefaultScale(), params.LogSlots()))
					}
					cnt++
				}

				// Update Cov as the Shifted Cov, We will use container ctMMT to temporarily store the eigvec^T @ eigvec
				for i := 0; i < len(ctEigVec[k]); i++ {
					fmt.Printf("Check the ctEigVec[%d][%d] to see wether it is really of form %d:\n", k, i, (VecMod+1)&(2-1))
					// //auxio.Quick_check_matrix_full(params, sk, ctEigVec[k][i], d, d)
					ctTemp = ctViT[i]
					fmt.Printf("The ctEigVec^T[%d][%d] is with scale %f, level %d: \n", k, i, math.Log2(ctTemp.Scale.Float64()), ctTemp.Level())
					//auxio.Quick_check_matrix(params, sk, ctTemp, d, d)
					if err != nil {
						panic(err)
					}
					for j := 0; j < len(ctEigVec[k]); j++ {
						// fmt.Printf("The i:%d, j:%d ctMMT: \n", i, j)
						if ((VecMod + 1) & (2 - 1)) == 1 {
							fmt.Printf("Internal turn: The ctEigVec^T[%d][%d], it will be multiplied with ctEigVec[%d][%d] to become ctMMT[%d][%d]\n:", k, i, k, j, j, i)
							ctMMT[j][i] = evaluator.MulRelinNew(ctTemp, ctEigVec[k][j])
							evaluator.Rescale(ctMMT[j][i], params.DefaultScale(), ctMMT[j][i]) // Consume one level
							ctMMT[j][i] = evaluator.MulRelinNew(ctMMT[j][i], ctEigVal[k])
							evaluator.Rescale(ctMMT[j][i], params.DefaultScale(), ctMMT[j][i]) // Consume one level
							fmt.Printf("The eigvec^T @ eigvec * eigval[%d][%d] with scale %f, level %d: \n", i, j, math.Log2(ctMMT[j][i].Scale.Float64()), ctMMT[j][i].Level())
							if ctMMT[j][i].Level() >= ctCov[j][i].Level() {
								ctMMT[j][i].SetScale(ctCov[j][i].Scale)
							} else {
								ctCov[j][i].SetScale(ctMMT[j][i].Scale)
							}
							//auxio.Quick_check_matrix(params, sk, ctMMT[j][i], d, d)

							fmt.Printf("The Original ctCov[%d][%d] with scale %f, level %d: \n", i, j, math.Log2(ctCov[j][i].Scale.Float64()), ctCov[j][i].Level())
							//auxio.Quick_check_matrix(params, sk, ctCov[j][i], d, d)
							ctCov[j][i] = evaluator.SubNew(ctCov[j][i], ctMMT[j][i])
							fmt.Printf("The Shifted ctCov[%d][%d] with scale %f, level %d: \n", i, j, math.Log2(ctCov[j][i].Scale.Float64()), ctCov[j][i].Level())
							//auxio.Quick_check_matrix(params, sk, ctCov[j][i], d, d)
						} else if ((VecMod + 1) & (2 - 1)) == 0 { // idx should be flipped here .
							ctMMT[i][j] = evaluator.MulRelinNew(ctTemp, ctEigVec[k][j])
							evaluator.Rescale(ctMMT[i][j], params.DefaultScale(), ctMMT[i][j]) // Consume one level
							ctMMT[i][j] = evaluator.MulRelinNew(ctMMT[i][j], ctEigVal[k])
							evaluator.Rescale(ctMMT[i][j], params.DefaultScale(), ctMMT[i][j]) // Consume one level
							fmt.Printf("The eigvec^T @ eigvec * eigval[%d][%d] with scale %f, level %d: \n", i, j, math.Log2(ctMMT[i][j].Scale.Float64()), ctMMT[i][j].Level())
							if ctMMT[i][j].Level() >= ctCov[i][j].Level() {
								ctMMT[i][j].SetScale(ctCov[i][j].Scale)
							} else {
								ctCov[i][j].SetScale(ctMMT[i][j].Scale)
							}
							//auxio.Quick_check_matrix(params, sk, ctMMT[i][j], d, d)

							fmt.Printf("The Original ctCov[%d][%d] with scale %f, level %d: \n", i, j, math.Log2(ctCov[j][i].Scale.Float64()), ctCov[i][j].Level())
							//auxio.Quick_check_matrix(params, sk, ctCov[i][j], d, d)
							ctCov[i][j] = evaluator.SubNew(ctCov[i][j], ctMMT[i][j])
							fmt.Printf("The Shifted ctCov[%d][%d] with scale %f, level %d: \n", i, j, math.Log2(ctCov[i][j].Scale.Float64()), ctCov[i][j].Level())
							//auxio.Quick_check_matrix(params, sk, ctCov[i][j], d, d)
						}

					}
				}
			}

		}
	}
	elapsed += time.Since(now)
	fmt.Printf("PowerMethod Done with btp times %d with total btptime %s, and total time %s\n", cnt, btpelapsed, time.Since(now))
	ctEigenVectors := make([][]float64, TargetEigVecNum)
	EigVecs := make([][]float64, TargetEigVecNum)
	for i := range ctEigenVectors {
		EigVec := make([]complex128, 0)
		for j := range ctEigVec[i] {
			fmt.Printf("Now extract partial EigenVector from ctEigVec[%d][%d]\n", i, j)
			//auxio.Quick_check_matrix_full(params, sk, ctEigVec[i][j], d, d)
			partialEigVec_dup := auxio.DecryptDecode(params, sk, ctEigVec[i][j])
			partialEigVec := make([]complex128, matrixCols)
			for k := 0; k < matrixCols; k++ {
				partialEigVec[k] = partialEigVec_dup[k*d+1] // Do not take the first column since it's precision is worse.
			}
			EigVec = append(EigVec, partialEigVec...)
		}
		EigVecs[i] = auxio.Complex2Float(EigVec)
	}
	auxio.ExportToCSV(EigVecs, "EigenVectors")

}

func LoadKeys_btpVersion(path string) (params ckks.Parameters, btpParams bootstrapping.Parameters, sk *rlwe.SecretKey, pk *rlwe.PublicKey, btprlk *rlwe.RelinearizationKey, btpgalks *rlwe.RotationKeySet, swkDtS *rlwe.SwitchingKey, swkStD *rlwe.SwitchingKey, err error) {
	// Allocate Memory for keys
	sk = new(rlwe.SecretKey)
	pk = new(rlwe.PublicKey)
	btprlk = new(rlwe.RelinearizationKey)
	btpgalks = new(rlwe.RotationKeySet)
	swkDtS = new(rlwe.SwitchingKey)
	swkStD = new(rlwe.SwitchingKey)

	// create Keys' File Tracker.
	var keysTracker = auxio.NewTracker4File(path)
	var n int

	// Load sk,pk, etc...
	n, err = keysTracker.ReadUpdateOne(&params)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Read params %d bytes\n", n)
	}

	n, err = keysTracker.ReadUpdateOne(&btpParams)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Read btpParams %d bytes\n", n)
	}

	n, err = keysTracker.ReadUpdateOne(sk)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Read sk %d bytes\n", n)
	}

	n, err = keysTracker.ReadUpdateOne(pk)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Read pk %d bytes\n", n)
	}

	n, err = keysTracker.ReadUpdateOne(btprlk)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Read btprlk %d bytes\n", n)
	}

	n, err = keysTracker.ReadUpdateOne(btpgalks)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Read btpgalks %d bytes\n", n)
	}

	n, err = keysTracker.ReadUpdateOne(swkDtS)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Read swkDtS %d bytes\n", n)
	}

	n, err = keysTracker.ReadUpdateOne(swkStD)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Read swkStD %d bytes\n", n)
	}

	return

}

func StoreKeys_btpVersion(path string, params ckks.Parameters, btpParams bootstrapping.Parameters, sk *rlwe.SecretKey, pk *rlwe.PublicKey, btprlk *rlwe.RelinearizationKey, btpgalks *rlwe.RotationKeySet, swkDtS *rlwe.SwitchingKey, swkStD *rlwe.SwitchingKey) (n int, err error) {
	keysTracker := auxio.NewTracker(path, -1)
	n = 0
	n, err = keysTracker.StoreUpdateOne(&params)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Store params %d bytes\n", n)
	}

	n, err = keysTracker.StoreUpdateOne(&btpParams)
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

	n, err = keysTracker.StoreUpdateOne(btprlk)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Store btprlk %d bytes\n", n)
	}

	n, err = keysTracker.StoreUpdateOne(btpgalks)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Store btpgalks %d bytes\n", n)
	}

	n, err = keysTracker.StoreUpdateOne(swkDtS)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Store swkDtS %d bytes\n", n)
	}

	n, err = keysTracker.StoreUpdateOne(swkStD)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Store swkStD %d bytes\n", n)
	}

	n, err = keysTracker.StoreFinish()
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Store filetracker %d bytes\n", n)
	}
	return
}

func StoreAndLoadKeys_btpVersion_check(loadkeysFromDisk bool) {
	var err error

	var KeysPath_btpVersion = "Keys_btpVer"

	// initialize encryption scheme params.

	var ckksParams ckks.ParametersLiteral
	var btpParams bootstrapping.Parameters
	var params ckks.Parameters

	var sk *rlwe.SecretKey
	var pk *rlwe.PublicKey
	// var rlk *rlwe.RelinearizationKey
	var btpevk bootstrapping.EvaluationKeys
	var kgen rlwe.KeyGenerator

	// generate Keys, we can load keys from disk or generate them.
	if loadkeysFromDisk {
		// In this function, we do not need to perform boostrapping, so we only pick up the keys we need in this section.
		params, btpParams, sk, pk, btpevk.Rlk, btpevk.Rtks, btpevk.SwkDtS, btpevk.SwkStD, err = LoadKeys_btpVersion(KeysPath_btpVersion)
		if err != nil {
			panic(err)
		}
	} else {
		paramSet := bootstrapping.N16QP1788H32768H32
		ckksParams = paramSet.SchemeParams
		btpParams = paramSet.BootstrappingParams
		ckksParams.LogSlots = 14
		params, err = ckks.NewParametersFromLiteral(ckksParams)
		if err != nil {
			panic(err)
		}
		fmt.Println()
		fmt.Printf("CKKS parameters: logN = %d, logSlots = %d, H(%d; %d), logQP = %d, levels = %d, scale= 2^%f, sigma = %f \n", params.LogN(), params.LogSlots(), params.HammingWeight(), btpParams.EphemeralSecretWeight, params.LogQP(), params.QCount(), math.Log2(params.DefaultScale().Float64()), params.Sigma())

		kgen = ckks.NewKeyGenerator(params)
		sk, pk = kgen.GenKeyPair()
		fmt.Println()
		fmt.Println("Generating bootstrapping keys...")
		btpevk = bootstrapping.GenEvaluationKeys(btpParams, params, sk)
		fmt.Println("Done")

		// Print Keys' size:
		fmt.Printf("Encryption Scheme btpParams occupy %d bytes\n", btpParams.MarshalBinarySize())
		fmt.Printf("Encryption Scheme params    occupy %d bytes\n", params.MarshalBinarySize())
		fmt.Printf("Encryption Scheme secretkey occupy %d bytes\n", sk.MarshalBinarySize())
		fmt.Printf("Encryption Scheme publickey occupy %d bytes\n", pk.MarshalBinarySize())
		fmt.Printf("(btp) Relinearization Key   occupy %d bytes\n", btpevk.Rlk.MarshalBinarySize())
		fmt.Printf("(btp) Rotation Keys         occupy %d bytes\n", btpevk.Rtks.MarshalBinarySize())
		fmt.Printf("(btp) DtS Switching Key     occupy %d bytes\n", btpevk.SwkDtS.MarshalBinarySize())
		fmt.Printf("(btp) StD Switching Key     occupy %d bytes\n", btpevk.SwkStD.MarshalBinarySize())

		_, err = StoreKeys_btpVersion(KeysPath_btpVersion, params, btpParams, sk, pk, btpevk.Rlk, btpevk.Rtks, btpevk.SwkDtS, btpevk.SwkStD)
		if err != nil {
			panic(err)
		}
	}

	// let rlk be a copy of btpevk.Rlk, they share the same memory address.
	//rlk = btpevk.Rlk // History lefting problems... rlk should be deleted, but...

	// bootstrapper
	var btp *bootstrapping.Bootstrapper
	if btp, err = bootstrapping.NewBootstrapper(params, btpParams, btpevk); err != nil {
		panic(err)
	}

	// Encryptor
	encryptor := ckks.NewEncryptor(params, pk)

	// Decryptor
	decryptor := ckks.NewDecryptor(params, sk)

	// encoder
	encoder := ckks.NewEncoder(params)

	// evaluator
	// evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})

	// In the Boostrapping Version, we do not begin the performance of operations from the top level,
	// instead we will begin at the first level which is not included in bootstrap procedure, the number of this level can be
	// computed by:
	MaxLevel := btpParams.SlotsToCoeffsParameters.LevelStart - len(btpParams.SlotsToCoeffsParameters.ScalingFactor)

	d := 128

	MatrixA := make([][]float64, d)
	for i := 0; i < d; i++ {
		MatrixA[i] = make([]float64, d)
		for j := 0; j < d; j++ {
			MatrixA[i][j] = float64(i*d + j)
		}
	}

	// Generate a random plaintext
	valuesWant := make([]complex128, params.Slots())
	valuesWant2 := make([]complex128, params.Slots())
	for i := range valuesWant {
		valuesWant[i] = complex(utils.RandFloat64(-10, 10), 0)
		valuesWant2[i] = valuesWant[i] * valuesWant[i]
	}

	plaintext := encoder.EncodeNew(valuesWant, MaxLevel, rlwe.NewScale(float64(params.RingQ().Modulus[MaxLevel])), params.LogSlots())

	// Encrypt
	ciphertext1 := encryptor.EncryptNew(plaintext)
	btp.MulRelin(ciphertext1, ciphertext1, ciphertext1)
	err = btp.Rescale(ciphertext1, params.DefaultScale(), ciphertext1)
	btp.DropLevel(ciphertext1, MaxLevel-2)
	btp.SetScale(ciphertext1, params.DefaultScale())
	if err != nil {
		panic(err)
	}

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of values vs. ciphertext")
	valuesTest1 := printDebug(params, ciphertext1, valuesWant2, decryptor, encoder)

	// Bootstrap the ciphertext (homomorphic re-encryption)
	// It takes a ciphertext at level 0 (if not at level 0, then it will reduce it to level 0)
	// and returns a ciphertext at level MaxLevel - k, where k is the depth of the bootstrapping circuit.
	// CAUTION: the scale of the ciphertext MUST be equal (or very close) to params.DefaultScale()
	// To equalize the scale, the function evaluator.SetScale(ciphertext, parameters.DefaultScale()) can be used at the expense of one level.
	fmt.Println()

	fmt.Printf("cihphtext now at level %d, with scale %f\n", ciphertext1.Level(), math.Log2(ciphertext1.Scale.Float64()))
	fmt.Println("Bootstrapping...")
	now := time.Now()
	ciphertext2 := btp.Bootstrap(ciphertext1)
	fmt.Printf("Done in %s\n", time.Since(now))

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of ciphertext vs. Bootstrapp(ciphertext)")
	printDebug(params, ciphertext2, valuesTest1, decryptor, encoder)

}

func printDebug(params ckks.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []complex128, decryptor rlwe.Decryptor, encoder ckks.Encoder) (valuesTest []complex128) {

	valuesTest = encoder.Decode(decryptor.DecryptNew(ciphertext), params.LogSlots())

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))

	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale.Float64()))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])

	precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, params.LogSlots(), 0)

	fmt.Println(precStats.String())
	fmt.Println()

	return
}

func generateRandomArray(t int) []float64 {
	rand.Seed(time.Now().UnixNano()) // 使用当前时间作为随机种子

	arr := make([]float64, t) // 创建长度为t的切片

	for i := 0; i < t; i++ {
		r := rand.Float64()
		arr[i] = r // 生成0~1之间的随机浮点数
	}
	fmt.Printf("Random: %f,%f,%f,%f", arr[0], arr[1], arr[2], arr[3])

	return arr
}