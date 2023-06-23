package ckks

import (
	"fmt"
	"runtime"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// TraceNew maps X -> sum((-1)^i * X^{i*n+1}) for 0 <= i < N and returns the result on a new ciphertext.
// For log(n) = logSlots.
func (eval *evaluator) TraceNew(ctIn *rlwe.Ciphertext, logSlots int) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ctIn.Level())
	eval.Trace(ctIn, logSlots, ctOut)
	return
}

// Average returns the average of vectors of batchSize elements.
// The operation assumes that ctIn encrypts SlotCount/'batchSize' sub-vectors of size 'batchSize'.
// It then replaces all values of those sub-vectors by the component-wise average between all the sub-vectors.
// Example for batchSize=4 and slots=8: [{a, b, c, d}, {e, f, g, h}] -> [0.5*{a+e, b+f, c+g, d+h}, 0.5*{a+e, b+f, c+g, d+h}]
// Operation requires log2(SlotCout/'batchSize') rotations.
// Required rotation keys can be generated with 'RotationsForInnerSumLog(batchSize, SlotCount/batchSize)”
func (eval *evaluator) Average(ctIn *rlwe.Ciphertext, logBatchSize int, ctOut *rlwe.Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("ctIn.Degree() != 1 or ctOut.Degree() != 1")
	}

	if logBatchSize > eval.params.LogSlots() {
		panic("cannot Average: batchSize must be smaller or equal to the number of slots")
	}

	ringQ := eval.params.RingQ()

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	n := eval.params.Slots() / (1 << logBatchSize)

	// pre-multiplication by n^-1
	for i := 0; i < level+1; i++ {
		Q := ringQ.Modulus[i]
		bredParams := ringQ.BredParams[i]
		mredparams := ringQ.MredParams[i]
		invN := ring.ModExp(uint64(n), Q-2, Q)
		invN = ring.MForm(invN, Q, bredParams)

		ring.MulScalarMontgomeryVec(ctIn.Value[0].Coeffs[i], ctOut.Value[0].Coeffs[i], invN, Q, mredparams)
		ring.MulScalarMontgomeryVec(ctIn.Value[1].Coeffs[i], ctOut.Value[1].Coeffs[i], invN, Q, mredparams)
	}

	eval.InnerSum(ctOut, 1<<logBatchSize, n, ctOut)
}

// RotateHoistedNew takes an input Ciphertext and a list of rotations and returns a map of Ciphertext, where each element of the map is the input Ciphertext
// rotation by one element of the list. It is much faster than sequential calls to Rotate.
func (eval *evaluator) RotateHoistedNew(ctIn *rlwe.Ciphertext, rotations []int) (ctOut map[int]*rlwe.Ciphertext) {
	ctOut = make(map[int]*rlwe.Ciphertext)
	for _, i := range rotations {
		ctOut[i] = NewCiphertext(eval.params, 1, ctIn.Level())
	}
	eval.RotateHoisted(ctIn, rotations, ctOut)
	return
}

// RotateHoisted takes an input Ciphertext and a list of rotations and populates a map of pre-allocated Ciphertexts,
// where each element of the map is the input Ciphertext rotation by one element of the list.
// It is much faster than sequential calls to Rotate.
func (eval *evaluator) RotateHoisted(ctIn *rlwe.Ciphertext, rotations []int, ctOut map[int]*rlwe.Ciphertext) {
	levelQ := ctIn.Level()
	eval.DecomposeNTT(levelQ, eval.params.PCount()-1, eval.params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)
	for _, i := range rotations {
		eval.AutomorphismHoisted(levelQ, ctIn, eval.BuffDecompQP, eval.params.GaloisElementForColumnRotationBy(i), ctOut[i])
	}
}

// LinearTransform is a type for linear transformations on ciphertexts.
// It stores a plaintext matrix diagonalized in diagonal form and
// can be evaluated on a ciphertext by using the evaluator.LinearTransform method.
type LinearTransform struct {
	LogSlots   int                 // Log of the number of slots of the plaintext (needed to compute the appropriate rotation keys)
	N1         int                 // N1 is the number of inner loops of the baby-step giant-step algorithm used in the evaluation (if N1 == 0, BSGS is not used).
	Level      int                 // Level is the level at which the matrix is encoded (can be circuit dependent)
	Scale      rlwe.Scale          // Scale is the scale at which the matrix is encoded (can be circuit dependent)
	Vec        map[int]ringqp.Poly // Vec is the matrix, in diagonal form, where each entry of vec is an indexed non-zero diagonal.
	CommonDiff int                 // the Common Difference of the Diagonalvectors' Arithmetic Sequence.
}

// NewLinearTransform allocates a new LinearTransform with zero plaintexts at the specified level.
// If BSGSRatio == 0, the LinearTransform is set to not use the BSGS approach.
// Method will panic if BSGSRatio < 0.
func NewLinearTransform(params Parameters, nonZeroDiags []int, level, logSlots int, BSGSRatio float64) LinearTransform {
	vec := make(map[int]ringqp.Poly)
	slots := 1 << logSlots
	levelQ := level
	levelP := params.PCount() - 1
	var N1 int
	if BSGSRatio == 0 {
		N1 = 0
		for _, i := range nonZeroDiags {
			idx := i
			if idx < 0 {
				idx += slots
			}
			vec[idx] = params.RingQP().NewPolyLvl(levelQ, levelP)
		}
	} else if BSGSRatio > 0 {
		N1 = FindBestBSGSSplit(nonZeroDiags, slots, BSGSRatio)
		index, _, _ := BsgsIndex(nonZeroDiags, slots, N1)
		for j := range index {
			for _, i := range index[j] {
				vec[j+i] = params.RingQP().NewPolyLvl(levelQ, levelP)
			}
		}
	} else {
		panic("cannot NewLinearTransform: BSGS ratio cannot be negative")
	}

	return LinearTransform{LogSlots: logSlots, N1: N1, Level: level, Vec: vec}
}

// Rotations returns the list of rotations needed for the evaluation
// of the linear transform.
func (LT *LinearTransform) Rotations() (rotations []int) {
	slots := 1 << LT.LogSlots

	rotIndex := make(map[int]bool)

	var index int

	N1 := LT.N1

	if LT.N1 == 0 {

		for j := range LT.Vec {
			rotIndex[j] = true
		}

	} else {

		for j := range LT.Vec {

			index = ((j / N1) * N1) & (slots - 1)
			rotIndex[index] = true

			index = j & (N1 - 1)
			rotIndex[index] = true
		}
	}

	rotations = make([]int, len(rotIndex))
	var i int
	for j := range rotIndex {
		rotations[i] = j
		i++
	}

	return rotations
}

// Encode encodes on a pre-allocated LinearTransform the linear transforms' matrix in diagonal form `value`.
// values.(type) can be either map[int][]complex128 or map[int][]float64.
// User must ensure that 1 <= len([]complex128\[]float64) <= 2^logSlots < 2^logN.
// It can then be evaluated on a ciphertext using evaluator.LinearTransform.
// Evaluation will use the naive approach (single hoisting and no baby-step giant-step).
// Faster if there is only a few non-zero diagonals but uses more keys.
func (LT *LinearTransform) Encode(encoder Encoder, value interface{}, scale rlwe.Scale) {

	enc, ok := encoder.(*encoderComplex128)
	if !ok {
		panic("cannot Encode: encoder should be an encoderComplex128")
	}

	dMat := interfaceMapToMapOfInterface(value)
	slots := 1 << LT.LogSlots
	N1 := LT.N1

	if N1 == 0 {
		for i := range dMat {
			idx := i
			if idx < 0 {
				idx += slots
			}

			if _, ok := LT.Vec[idx]; !ok {
				panic("cannot Encode: error encoding on LinearTransform: input does not match the same non-zero diagonals")
			}

			enc.Embed(dMat[i], LT.LogSlots, scale, true, LT.Vec[idx])
		}
	} else {

		index, _, _ := BsgsIndex(value, slots, N1)

		var values interface{}
		switch value.(type) {
		case map[int][]complex128:
			values = make([]complex128, slots)
		case map[int][]float64:
			values = make([]float64, slots)
		}

		for j := range index {

			rot := -j & (slots - 1)

			for _, i := range index[j] {
				// manages inputs that have rotation between 0 and slots-1 or between -slots/2 and slots/2-1
				v, ok := dMat[j+i]
				if !ok {
					v = dMat[j+i-slots]
				}

				if _, ok := LT.Vec[j+i]; !ok {
					panic("cannot Encode: error encoding on LinearTransform BSGS: input does not match the same non-zero diagonals")
				}

				copyRotInterface(values, v, rot)

				enc.Embed(values, LT.LogSlots, scale, true, LT.Vec[j+i])
			}
		}
	}

	LT.Scale = scale
}

// GenLinearTransform allocates and encode a new LinearTransform struct from the linear transforms' matrix in diagonal form `value`.
// values.(type) can be either map[int][]complex128 or map[int][]float64.
// User must ensure that 1 <= len([]complex128\[]float64) <= 2^logSlots < 2^logN.
// It can then be evaluated on a ciphertext using evaluator.LinearTransform.
// Evaluation will use the naive approach (single hoisting and no baby-step giant-step).
// Faster if there is only a few non-zero diagonals but uses more keys.
func GenLinearTransform(encoder Encoder, value interface{}, level int, scale rlwe.Scale, logslots int) LinearTransform {

	enc, ok := encoder.(*encoderComplex128)
	if !ok {
		panic("cannot GenLinearTransform: encoder should be an encoderComplex128")
	}

	params := enc.params
	dMat := interfaceMapToMapOfInterface(value)
	vec := make(map[int]ringqp.Poly)
	slots := 1 << logslots
	levelQ := level
	levelP := params.PCount() - 1
	for i := range dMat {

		idx := i
		if idx < 0 {
			idx += slots
		}
		vec[idx] = params.RingQP().NewPolyLvl(levelQ, levelP)
		enc.Embed(dMat[i], logslots, scale, true, vec[idx])
	}

	return LinearTransform{LogSlots: logslots, N1: 0, Vec: vec, Level: level, Scale: scale}
}

// GenLinearTransformBSGS allocates and encodes a new LinearTransform struct from the linear transforms' matrix in diagonal form `value` for evaluation with a baby-step giant-step approach.
// values.(type) can be either map[int][]complex128 or map[int][]float64.
// User must ensure that 1 <= len([]complex128\[]float64) <= 2^logSlots < 2^logN.
// LinearTransform types can be be evaluated on a ciphertext using evaluator.LinearTransform.
// Evaluation will use the optimized approach (double hoisting and baby-step giant-step).
// Faster if there is more than a few non-zero diagonals.
// BSGSRatio is the maximum ratio between the inner and outer loop of the baby-step giant-step algorithm used in evaluator.LinearTransform.
// Optimal BSGSRatio value is between 4 and 16 depending on the sparsity of the matrix.
func GenLinearTransformBSGS(encoder Encoder, value interface{}, level int, scale rlwe.Scale, BSGSRatio float64, logSlots int) (LT LinearTransform) {

	enc, ok := encoder.(*encoderComplex128)
	if !ok {
		panic("cannot GenLinearTransformBSGS: encoder should be an encoderComplex128")
	}

	params := enc.params

	slots := 1 << logSlots

	// N1*N2 = N
	N1 := FindBestBSGSSplit(value, slots, BSGSRatio)

	index, _, _ := BsgsIndex(value, slots, N1)

	vec := make(map[int]ringqp.Poly)

	dMat := interfaceMapToMapOfInterface(value)
	levelQ := level
	levelP := params.PCount() - 1

	var values interface{}
	switch value.(type) {
	case map[int][]complex128:
		values = make([]complex128, slots)
	case map[int][]float64:
		values = make([]float64, slots)
	}

	for j := range index {

		rot := -j & (slots - 1)

		for _, i := range index[j] {

			// manages inputs that have rotation between 0 and slots-1 or between -slots/2 and slots/2-1
			v, ok := dMat[j+i]
			if !ok {
				v = dMat[j+i-slots]
			}
			vec[j+i] = params.RingQP().NewPolyLvl(levelQ, levelP)

			copyRotInterface(values, v, rot)

			enc.Embed(values, logSlots, scale, true, vec[j+i])
		}
	}

	return LinearTransform{LogSlots: logSlots, N1: N1, Vec: vec, Level: level, Scale: scale}
}

func copyRotInterface(a, b interface{}, rot int) {
	switch a.(type) {
	case []complex128:

		ac128 := a.([]complex128)
		bc128 := b.([]complex128)

		n := len(ac128)

		if len(bc128) >= rot {
			copy(ac128[:n-rot], bc128[rot:])
			copy(ac128[n-rot:], bc128[:rot])
		} else {
			copy(ac128[n-rot:], bc128)
		}
	case []float64:

		af64 := a.([]float64)
		bf64 := b.([]float64)

		n := len(af64)

		if len(bf64) >= rot {
			copy(af64[:n-rot], bf64[rot:])
			copy(af64[n-rot:], bf64[:rot])
		} else {
			copy(af64[n-rot:], bf64)
		}
	}
}

// BsgsIndex returns the index map and needed rotation for the BSGS matrix-vector multiplication algorithm.
func BsgsIndex(el interface{}, slots, N1 int) (index map[int][]int, rotN1, rotN2 []int) {
	index = make(map[int][]int)
	rotN1Map := make(map[int]bool)
	rotN2Map := make(map[int]bool)
	var nonZeroDiags []int
	switch element := el.(type) {
	case map[int][]complex128:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case map[int][]float64:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case map[int]bool:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case map[int]ringqp.Poly:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case []int:
		nonZeroDiags = element
	}

	for _, rot := range nonZeroDiags {
		rot &= (slots - 1)
		idxN1 := ((rot / N1) * N1) & (slots - 1)
		idxN2 := rot & (N1 - 1)
		if index[idxN1] == nil {
			index[idxN1] = []int{idxN2}
		} else {
			index[idxN1] = append(index[idxN1], idxN2)
		}
		rotN1Map[idxN1] = true
		rotN2Map[idxN2] = true
	}

	rotN1 = []int{}
	for i := range rotN1Map {
		rotN1 = append(rotN1, i)
	}

	rotN2 = []int{}
	for i := range rotN2Map {
		rotN2 = append(rotN2, i)
	}

	return
}

func interfaceMapToMapOfInterface(m interface{}) map[int]interface{} {
	d := make(map[int]interface{})
	switch el := m.(type) {
	case map[int][]complex128:
		for i := range el {
			d[i] = el[i]
		}
	case map[int][]float64:
		for i := range el {
			d[i] = el[i]
		}
	default:
		panic("cannot interfaceMapToMapOfInterface: invalid input, must be map[int][]complex128 or map[int][]float64")
	}
	return d
}

// FindBestBSGSSplit finds the best N1*N2 = N for the baby-step giant-step algorithm for matrix multiplication.
func FindBestBSGSSplit(diagMatrix interface{}, maxN int, maxRatio float64) (minN int) {

	for N1 := 1; N1 < maxN; N1 <<= 1 {

		_, rotN1, rotN2 := BsgsIndex(diagMatrix, maxN, N1)

		nbN1, nbN2 := len(rotN1)-1, len(rotN2)-1

		if float64(nbN2)/float64(nbN1) == maxRatio {
			return N1
		}

		if float64(nbN2)/float64(nbN1) > maxRatio {
			return N1 / 2
		}
	}

	return 1
}

// LinearTransformNew evaluates a linear transform on the Ciphertext "ctIn" and returns the result on a new Ciphertext.
// The linearTransform can either be an (ordered) list of PtDiagMatrix or a single PtDiagMatrix.
// In either case, a list of Ciphertext is returned (the second case returning a list
// containing a single Ciphertext). A PtDiagMatrix is a diagonalized plaintext matrix constructed with an Encoder using
// the method encoder.EncodeDiagMatrixAtLvl(*).
func (eval *evaluator) LinearTransformNew(ctIn *rlwe.Ciphertext, linearTransform interface{}) (ctOut []*rlwe.Ciphertext) {

	switch LTs := linearTransform.(type) {
	case []LinearTransform:
		ctOut = make([]*rlwe.Ciphertext, len(LTs))

		var maxLevel int
		for _, LT := range LTs {
			maxLevel = utils.MaxInt(maxLevel, LT.Level)
		}

		minLevel := utils.MinInt(maxLevel, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)

		for i, LT := range LTs {
			ctOut[i] = NewCiphertext(eval.params, 1, minLevel)

			if LT.N1 == 0 {
				eval.MultiplyByDiagMatrix(ctIn, LT, eval.BuffDecompQP, ctOut[i])
			} else {
				eval.MultiplyByDiagMatrixBSGS(ctIn, LT, eval.BuffDecompQP, ctOut[i])
			}

			ctOut[i].MetaData = ctIn.MetaData
			ctOut[i].Scale = ctIn.Scale.Mul(LT.Scale)
		}

	case LinearTransform:

		minLevel := utils.MinInt(LTs.Level, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)

		ctOut = []*rlwe.Ciphertext{NewCiphertext(eval.params, 1, minLevel)}

		if LTs.N1 == 0 {
			eval.MultiplyByDiagMatrix(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		} else {
			eval.MultiplyByDiagMatrixBSGS(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		}

		ctOut[0].MetaData = ctIn.MetaData
		ctOut[0].Scale = ctIn.Scale.Mul(LTs.Scale)
	}
	return
}

// LinearTransform evaluates a linear transform on the pre-allocated Ciphertexts.
// The linearTransform can either be an (ordered) list of PtDiagMatrix or a single PtDiagMatrix.
// In either case a list of Ciphertext is returned (the second case returning a list
// containing a single Ciphertext). A PtDiagMatrix is a diagonalized plaintext matrix constructed with an Encoder using
// the method encoder.EncodeDiagMatrixAtLvl(*).
func (eval *evaluator) LinearTransform(ctIn *rlwe.Ciphertext, linearTransform interface{}, ctOut []*rlwe.Ciphertext) {

	switch LTs := linearTransform.(type) {
	case []LinearTransform:
		var maxLevel int
		for _, LT := range LTs {
			maxLevel = utils.MaxInt(maxLevel, LT.Level)
		}

		minLevel := utils.MinInt(maxLevel, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)

		for i, LT := range LTs {
			if LT.N1 == 0 {
				eval.MultiplyByDiagMatrix(ctIn, LT, eval.BuffDecompQP, ctOut[i])
			} else {
				eval.MultiplyByDiagMatrixBSGS(ctIn, LT, eval.BuffDecompQP, ctOut[i])
			}

			ctOut[i].MetaData = ctIn.MetaData
			ctOut[i].Scale = ctIn.Scale.Mul(LT.Scale)
		}

	case LinearTransform:
		minLevel := utils.MinInt(LTs.Level, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)
		if LTs.N1 == 0 {
			eval.MultiplyByDiagMatrix(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		} else {
			eval.MultiplyByDiagMatrixBSGS(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		}

		ctOut[0].MetaData = ctIn.MetaData
		ctOut[0].Scale = ctIn.Scale.Mul(LTs.Scale)
	}
}

// MultiplyByDiagMatrix multiplies the Ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the Ciphertext
// "ctOut". Memory buffers for the decomposed ciphertext BuffDecompQP, BuffDecompQP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The naive approach is used (single hoisting and no baby-step giant-step), which is faster than MultiplyByDiagMatrixBSGS
// for matrix of only a few non-zero diagonals but uses more keys.
func (eval *evaluator) MultiplyByDiagMatrix(ctIn *rlwe.Ciphertext, matrix LinearTransform, BuffDecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	levelQ := utils.MinInt(ctOut.Level(), utils.MinInt(ctIn.Level(), matrix.Level))
	levelP := len(ringP.Modulus) - 1

	ctOut.Resize(ctOut.Degree(), levelQ)

	QiOverF := eval.params.QiOverflowMargin(levelQ)
	PiOverF := eval.params.PiOverflowMargin(levelP)

	c0OutQP := ringqp.Poly{Q: ctOut.Value[0], P: eval.BuffQP[5].Q}
	c1OutQP := ringqp.Poly{Q: ctOut.Value[1], P: eval.BuffQP[5].P}

	ct0TimesP := eval.BuffQP[0].Q // ct0 * P mod Q
	tmp0QP := eval.BuffQP[1]
	tmp1QP := eval.BuffQP[2]
	ksRes0QP := eval.BuffQP[3]
	ksRes1QP := eval.BuffQP[4]

	ring.CopyLvl(levelQ, ctIn.Value[0], eval.buffCt.Value[0])
	ring.CopyLvl(levelQ, ctIn.Value[1], eval.buffCt.Value[1])
	ctInTmp0, ctInTmp1 := eval.buffCt.Value[0], eval.buffCt.Value[1]

	ringQ.MulScalarBigintLvl(levelQ, ctInTmp0, ringP.ModulusAtLevel[levelP], ct0TimesP) // P*c0

	var state bool
	var cnt int
	for k := range matrix.Vec {

		k &= int((ringQ.NthRoot >> 2) - 1)

		if k == 0 {
			state = true
		} else {

			galEl := eval.params.GaloisElementForColumnRotationBy(k)

			rtk, generated := eval.Rtks.Keys[galEl]
			if !generated {
				panic("cannot MultiplyByDiagMatrix: switching key not available")
			}

			index := eval.PermuteNTTIndex[galEl]

			eval.KeyswitchHoistedNoModDown(levelQ, BuffDecompQP, rtk, ksRes0QP.Q, ksRes1QP.Q, ksRes0QP.P, ksRes1QP.P)
			ringQ.AddLvl(levelQ, ksRes0QP.Q, ct0TimesP, ksRes0QP.Q)
			ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, ksRes0QP, index, tmp0QP)
			ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, ksRes1QP, index, tmp1QP)

			if cnt == 0 {
				// keyswitch(c1_Q) = (d0_QP, d1_QP)
				ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, matrix.Vec[k], tmp0QP, c0OutQP)
				ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, matrix.Vec[k], tmp1QP, c1OutQP)
			} else {
				// keyswitch(c1_Q) = (d0_QP, d1_QP)
				ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, matrix.Vec[k], tmp0QP, c0OutQP)
				ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, matrix.Vec[k], tmp1QP, c1OutQP)
			}

			if cnt%QiOverF == QiOverF-1 {
				ringQ.ReduceLvl(levelQ, c0OutQP.Q, c0OutQP.Q)
				ringQ.ReduceLvl(levelQ, c1OutQP.Q, c1OutQP.Q)
			}

			if cnt%PiOverF == PiOverF-1 {
				ringP.ReduceLvl(levelP, c0OutQP.P, c0OutQP.P)
				ringP.ReduceLvl(levelP, c1OutQP.P, c1OutQP.P)
			}

			cnt++
		}
	}

	if cnt%QiOverF == 0 {
		ringQ.ReduceLvl(levelQ, c0OutQP.Q, c0OutQP.Q)
		ringQ.ReduceLvl(levelQ, c1OutQP.Q, c1OutQP.Q)
	}

	if cnt%PiOverF == 0 {
		ringP.ReduceLvl(levelP, c0OutQP.P, c0OutQP.P)
		ringP.ReduceLvl(levelP, c1OutQP.P, c1OutQP.P)
	}

	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c0OutQP.Q, c0OutQP.P, c0OutQP.Q) // sum(phi(c0 * P + d0_QP))/P
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c1OutQP.Q, c1OutQP.P, c1OutQP.Q) // sum(phi(d1_QP))/P

	if state { // Rotation by zero
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0].Q, ctInTmp0, c0OutQP.Q) // ctOut += c0_Q * plaintext
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0].Q, ctInTmp1, c1OutQP.Q) // ctOut += c1_Q * plaintext
	}
}

// MultiplyByDiagMatrixBSGS multiplies the Ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the Ciphertext
// "ctOut". Memory buffers for the decomposed Ciphertext BuffDecompQP, BuffDecompQP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The BSGS approach is used (double hoisting with baby-step giant-step), which is faster than MultiplyByDiagMatrix
// for matrix with more than a few non-zero diagonals and uses significantly less keys.
func (eval *evaluator) MultiplyByDiagMatrixBSGS(ctIn *rlwe.Ciphertext, matrix LinearTransform, PoolDecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	levelQ := utils.MinInt(ctOut.Level(), utils.MinInt(ctIn.Level(), matrix.Level))
	levelP := len(ringP.Modulus) - 1

	ctOut.Resize(ctOut.Degree(), levelQ)

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	// Computes the N2 rotations indexes of the non-zero rows of the diagonalized DFT matrix for the baby-step giant-step algorithm
	index, _, rotN2 := BsgsIndex(matrix.Vec, 1<<matrix.LogSlots, matrix.N1)

	ring.CopyLvl(levelQ, ctIn.Value[0], eval.buffCt.Value[0])
	ring.CopyLvl(levelQ, ctIn.Value[1], eval.buffCt.Value[1])

	ctInTmp0, ctInTmp1 := eval.buffCt.Value[0], eval.buffCt.Value[1]

	// Pre-rotates ciphertext for the baby-step giant-step algorithm, does not divide by P yet
	ctInRotQP := eval.RotateHoistedNoModDownNew(levelQ, rotN2, ctInTmp0, eval.BuffDecompQP)

	// Accumulator inner loop
	tmp0QP := eval.BuffQP[1]
	tmp1QP := eval.BuffQP[2]

	// Accumulator outer loop
	cQP := rlwe.CiphertextQP{Value: [2]ringqp.Poly{eval.BuffQP[3], eval.BuffQP[4]}}
	cQP.IsNTT = true

	// Result in QP
	c0OutQP := ringqp.Poly{Q: ctOut.Value[0], P: eval.BuffQP[5].Q}
	c1OutQP := ringqp.Poly{Q: ctOut.Value[1], P: eval.BuffQP[5].P}

	ringQ.MulScalarBigintLvl(levelQ, ctInTmp0, ringP.ModulusAtLevel[levelP], ctInTmp0) // P*c0
	ringQ.MulScalarBigintLvl(levelQ, ctInTmp1, ringP.ModulusAtLevel[levelP], ctInTmp1) // P*c1

	// OUTER LOOP
	var cnt0 int
	for j := range index {

		// INNER LOOP
		var cnt1 int
		for _, i := range index[j] {
			if i == 0 {
				if cnt1 == 0 {
					ringQ.MulCoeffsMontgomeryConstantLvl(levelQ, matrix.Vec[j].Q, ctInTmp0, tmp0QP.Q)
					ringQ.MulCoeffsMontgomeryConstantLvl(levelQ, matrix.Vec[j].Q, ctInTmp1, tmp1QP.Q)
					tmp0QP.P.Zero()
					tmp1QP.P.Zero()
				} else {
					ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, matrix.Vec[j].Q, ctInTmp0, tmp0QP.Q)
					ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, matrix.Vec[j].Q, ctInTmp1, tmp1QP.Q)
				}
			} else {
				if cnt1 == 0 {
					ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, matrix.Vec[j+i], ctInRotQP[i].Value[0], tmp0QP)
					ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, matrix.Vec[j+i], ctInRotQP[i].Value[1], tmp1QP)
				} else {
					ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, matrix.Vec[j+i], ctInRotQP[i].Value[0], tmp0QP)
					ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, matrix.Vec[j+i], ctInRotQP[i].Value[1], tmp1QP)
				}
			}

			if cnt1%QiOverF == QiOverF-1 {
				ringQ.ReduceLvl(levelQ, tmp0QP.Q, tmp0QP.Q)
				ringQ.ReduceLvl(levelQ, tmp1QP.Q, tmp1QP.Q)
			}

			if cnt1%PiOverF == PiOverF-1 {
				ringP.ReduceLvl(levelP, tmp0QP.P, tmp0QP.P)
				ringP.ReduceLvl(levelP, tmp1QP.P, tmp1QP.P)
			}

			cnt1++
		}

		if cnt1%QiOverF != 0 {
			ringQ.ReduceLvl(levelQ, tmp0QP.Q, tmp0QP.Q)
			ringQ.ReduceLvl(levelQ, tmp1QP.Q, tmp1QP.Q)
		}

		if cnt1%PiOverF != 0 {
			ringP.ReduceLvl(levelP, tmp0QP.P, tmp0QP.P)
			ringP.ReduceLvl(levelP, tmp1QP.P, tmp1QP.P)
		}

		// If j != 0, then rotates ((tmp0QP.Q, tmp0QP.P), (tmp1QP.Q, tmp1QP.P)) by N1*j and adds the result on ((cQP.Value[0].Q, cQP.Value[0].P), (cQP.Value[1].Q, cQP.Value[1].P))
		if j != 0 {

			// Hoisting of the ModDown of sum(sum(phi(d1) * plaintext))
			eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, tmp1QP.Q, tmp1QP.P, tmp1QP.Q) // c1 * plaintext + sum(phi(d1) * plaintext) + phi(c1) * plaintext mod Q

			galEl := eval.params.GaloisElementForColumnRotationBy(j)

			rtk, generated := eval.Rtks.Keys[galEl]
			if !generated {
				panic("cannot MultiplyByDiagMatrixBSGS: switching key not available")
			}

			rotIndex := eval.PermuteNTTIndex[galEl]

			eval.GadgetProductNoModDown(levelQ, tmp1QP.Q, rtk.GadgetCiphertext, cQP) // Switchkey(P*phi(tmpRes_1)) = (d0, d1) in base QP
			ringQP.AddLvl(levelQ, levelP, cQP.Value[0], tmp0QP, cQP.Value[0])

			// Outer loop rotations
			if cnt0 == 0 {

				ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, cQP.Value[0], rotIndex, c0OutQP)
				ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, cQP.Value[1], rotIndex, c1OutQP)
			} else {
				ringQP.PermuteNTTWithIndexAndAddNoModLvl(levelQ, levelP, cQP.Value[0], rotIndex, c0OutQP)
				ringQP.PermuteNTTWithIndexAndAddNoModLvl(levelQ, levelP, cQP.Value[1], rotIndex, c1OutQP)
			}

			// Else directly adds on ((cQP.Value[0].Q, cQP.Value[0].P), (cQP.Value[1].Q, cQP.Value[1].P))
		} else {
			if cnt0 == 0 {
				ringQP.CopyLvl(levelQ, levelP, tmp0QP, c0OutQP)
				ringQP.CopyLvl(levelQ, levelP, tmp1QP, c1OutQP)
			} else {
				ringQP.AddNoModLvl(levelQ, levelP, c0OutQP, tmp0QP, c0OutQP)
				ringQP.AddNoModLvl(levelQ, levelP, c1OutQP, tmp1QP, c1OutQP)
			}
		}

		if cnt0%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(levelQ, ctOut.Value[0], ctOut.Value[0])
			ringQ.ReduceLvl(levelQ, ctOut.Value[1], ctOut.Value[1])
		}

		if cnt0%PiOverF == PiOverF-1 {
			ringP.ReduceLvl(levelP, c0OutQP.P, c0OutQP.P)
			ringP.ReduceLvl(levelP, c1OutQP.P, c1OutQP.P)
		}

		cnt0++
	}

	if cnt0%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, ctOut.Value[0], ctOut.Value[0])
		ringQ.ReduceLvl(levelQ, ctOut.Value[1], ctOut.Value[1])
	}

	if cnt0%PiOverF != 0 {
		ringP.ReduceLvl(levelP, c0OutQP.P, c0OutQP.P)
		ringP.ReduceLvl(levelP, c1OutQP.P, c1OutQP.P)
	}

	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ctOut.Value[0], c0OutQP.P, ctOut.Value[0]) // sum(phi(c0 * P + d0_QP))/P
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ctOut.Value[1], c1OutQP.P, ctOut.Value[1]) // sum(phi(d1_QP))/P

	ctInRotQP = nil
	runtime.GC()
}

func BsgsIndex4ArithmeticSeq(el interface{}, slots, N1 int, commonDiff int) (index map[int][]int, rotN1, rotN2 []int) {
	index = make(map[int][]int)
	rotN1Map := make(map[int]bool)
	rotN2Map := make(map[int]bool)
	var nonZeroDiags []int
	switch element := el.(type) {
	case map[int][]complex128:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case map[int][]float64:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case map[int]bool:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case map[int]ringqp.Poly:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case []int:
		nonZeroDiags = element
	}
	for _, rot := range nonZeroDiags {
		// The nonZeroDiags may contains rotations that has form -k*(commonDiff) mod slots for some positive integer k, scale 1/commonDiff cant directly apply on them.
		// So we should check these rotations first.
		var rotDivByCD int
		if rot%commonDiff != 0 { // check if rotation is divisible by commonDiff, panic when commonDiff is zero
			rotDivByCD = (rot - slots)           // turn the rotations with form -k*(commonDiff) mod slots back to -k*(commonDiff) < -slots.
			rotDivByCD = rotDivByCD / commonDiff // then we can get rotDivByCD = -k
		} else {
			rotDivByCD = rot / commonDiff // rotations with form k*(commonDiff) mod slots for some integer k such that -slots < k*(commonDiff) < slots can be scaled by 1/commonDiff directly.
		}
		rotDivByCD &= (slots - 1)
		idxN1 := (((rotDivByCD / N1) * N1) * commonDiff) & (slots - 1) // outer loop index
		idxN2 := ((rotDivByCD & (N1 - 1)) * commonDiff) & (slots - 1)  // inner loop index
		if index[idxN1] == nil {
			index[idxN1] = []int{idxN2}
		} else {
			index[idxN1] = append(index[idxN1], idxN2)
		}
		rotN1Map[idxN1] = true
		rotN2Map[idxN2] = true
	}

	rotN1 = []int{}
	for i := range rotN1Map {
		rotN1 = append(rotN1, i)
	}
	rotN2 = []int{}
	for i := range rotN2Map {
		rotN2 = append(rotN2, i)
	}

	return

}

func (eval *evaluator) MultiplyByDiagMatrixBSGS4ArithmeticSeq(ctIn *rlwe.Ciphertext, matrix LinearTransform, PoolDecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext) {
	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	levelQ := utils.MinInt(ctOut.Level(), utils.MinInt(ctIn.Level(), matrix.Level))
	levelP := len(ringP.Modulus) - 1

	ctOut.Resize(ctOut.Degree(), levelQ)

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	// Computes the N2 rotations indexes of the non-zero rows of the diagonalized DFT matrix for the baby-step giant-step algorithm
	index, _, rotN2 := BsgsIndex4ArithmeticSeq(matrix.Vec, 1<<matrix.LogSlots, matrix.N1, matrix.CommonDiff)

	ring.CopyLvl(levelQ, ctIn.Value[0], eval.buffCt.Value[0])
	ring.CopyLvl(levelQ, ctIn.Value[1], eval.buffCt.Value[1])

	ctInTmp0, ctInTmp1 := eval.buffCt.Value[0], eval.buffCt.Value[1]

	// Pre-rotates ciphertext for the baby-step giant-step algorithm, does not divide by P yet
	//now := time.Now() // dbg testing
	ctInRotQP := eval.RotateHoistedNoModDownNew(levelQ, rotN2, ctInTmp0, eval.BuffDecompQP)
	//fmt.Printf("%s consumed for Pre-rotating %d ciphertexts for the baby-step giant-step algorithm, does not divide by P yet\n", time.Since(now), len(ctInRotQP)) // dbg testing

	// dbg testing:

	var bytes int
	for _, ct := range ctInRotQP {
		bytes += ct.MarshalBinarySize()
	}
	fmt.Printf("\nmiddle ciphertexts take up %d bytes\n", bytes)

	// Accumulator inner loop
	tmp0QP := eval.BuffQP[1]
	tmp1QP := eval.BuffQP[2]

	// Accumulator outer loop
	cQP := rlwe.CiphertextQP{Value: [2]ringqp.Poly{eval.BuffQP[3], eval.BuffQP[4]}}
	cQP.IsNTT = true

	// Result in QP
	c0OutQP := ringqp.Poly{Q: ctOut.Value[0], P: eval.BuffQP[5].Q}
	c1OutQP := ringqp.Poly{Q: ctOut.Value[1], P: eval.BuffQP[5].P}

	ringQ.MulScalarBigintLvl(levelQ, ctInTmp0, ringP.ModulusAtLevel[levelP], ctInTmp0) // P*c0
	ringQ.MulScalarBigintLvl(levelQ, ctInTmp1, ringP.ModulusAtLevel[levelP], ctInTmp1) // P*c1

	// OUTER LOOP
	var cnt0 int
	for j := range index {
		//now = time.Now()
		// INNER LOOP
		var cnt1 int
		for _, i := range index[j] {
			if i == 0 {
				if cnt1 == 0 {
					ringQ.MulCoeffsMontgomeryConstantLvl(levelQ, matrix.Vec[j].Q, ctInTmp0, tmp0QP.Q)
					ringQ.MulCoeffsMontgomeryConstantLvl(levelQ, matrix.Vec[j].Q, ctInTmp1, tmp1QP.Q)
					tmp0QP.P.Zero()
					tmp1QP.P.Zero()
				} else {
					ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, matrix.Vec[j].Q, ctInTmp0, tmp0QP.Q)
					ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, matrix.Vec[j].Q, ctInTmp1, tmp1QP.Q)
				}
			} else {
				if cnt1 == 0 {
					ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, matrix.Vec[j+i], ctInRotQP[i].Value[0], tmp0QP)
					ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, matrix.Vec[j+i], ctInRotQP[i].Value[1], tmp1QP)
				} else {
					ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, matrix.Vec[j+i], ctInRotQP[i].Value[0], tmp0QP)
					ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, matrix.Vec[j+i], ctInRotQP[i].Value[1], tmp1QP)
				}
			}

			if cnt1%QiOverF == QiOverF-1 {
				ringQ.ReduceLvl(levelQ, tmp0QP.Q, tmp0QP.Q)
				ringQ.ReduceLvl(levelQ, tmp1QP.Q, tmp1QP.Q)
			}

			if cnt1%PiOverF == PiOverF-1 {
				ringP.ReduceLvl(levelP, tmp0QP.P, tmp0QP.P)
				ringP.ReduceLvl(levelP, tmp1QP.P, tmp1QP.P)
			}

			cnt1++
			//fmt.Printf("	%s consumed for %dth outer loop's %dth inner loop\n", time.Since(now), cnt0, cnt1-1)
		}

		if cnt1%QiOverF != 0 {
			ringQ.ReduceLvl(levelQ, tmp0QP.Q, tmp0QP.Q)
			ringQ.ReduceLvl(levelQ, tmp1QP.Q, tmp1QP.Q)
		}

		if cnt1%PiOverF != 0 {
			ringP.ReduceLvl(levelP, tmp0QP.P, tmp0QP.P)
			ringP.ReduceLvl(levelP, tmp1QP.P, tmp1QP.P)
		}

		// If j != 0, then rotates ((tmp0QP.Q, tmp0QP.P), (tmp1QP.Q, tmp1QP.P)) by N1*j and adds the result on ((cQP.Value[0].Q, cQP.Value[0].P), (cQP.Value[1].Q, cQP.Value[1].P))
		if j != 0 {

			// Hoisting of the ModDown of sum(sum(phi(d1) * plaintext))
			eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, tmp1QP.Q, tmp1QP.P, tmp1QP.Q) // c1 * plaintext + sum(phi(d1) * plaintext) + phi(c1) * plaintext mod Q

			galEl := eval.params.GaloisElementForColumnRotationBy(j)

			rtk, generated := eval.Rtks.Keys[galEl]
			if !generated {
				panic("cannot MultiplyByDiagMatrixBSGS: switching key not available")
			}

			rotIndex := eval.PermuteNTTIndex[galEl]

			eval.GadgetProductNoModDown(levelQ, tmp1QP.Q, rtk.GadgetCiphertext, cQP) // Switchkey(P*phi(tmpRes_1)) = (d0, d1) in base QP
			ringQP.AddLvl(levelQ, levelP, cQP.Value[0], tmp0QP, cQP.Value[0])

			// Outer loop rotations
			if cnt0 == 0 {

				ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, cQP.Value[0], rotIndex, c0OutQP)
				ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, cQP.Value[1], rotIndex, c1OutQP)
			} else {
				ringQP.PermuteNTTWithIndexAndAddNoModLvl(levelQ, levelP, cQP.Value[0], rotIndex, c0OutQP)
				ringQP.PermuteNTTWithIndexAndAddNoModLvl(levelQ, levelP, cQP.Value[1], rotIndex, c1OutQP)
			}

			// Else directly adds on ((cQP.Value[0].Q, cQP.Value[0].P), (cQP.Value[1].Q, cQP.Value[1].P))
		} else {
			if cnt0 == 0 {
				ringQP.CopyLvl(levelQ, levelP, tmp0QP, c0OutQP)
				ringQP.CopyLvl(levelQ, levelP, tmp1QP, c1OutQP)
			} else {
				ringQP.AddNoModLvl(levelQ, levelP, c0OutQP, tmp0QP, c0OutQP)
				ringQP.AddNoModLvl(levelQ, levelP, c1OutQP, tmp1QP, c1OutQP)
			}
		}

		if cnt0%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(levelQ, ctOut.Value[0], ctOut.Value[0])
			ringQ.ReduceLvl(levelQ, ctOut.Value[1], ctOut.Value[1])
		}

		if cnt0%PiOverF == PiOverF-1 {
			ringP.ReduceLvl(levelP, c0OutQP.P, c0OutQP.P)
			ringP.ReduceLvl(levelP, c1OutQP.P, c1OutQP.P)
		}

		cnt0++
		//fmt.Printf("%s consumed for %dth outer loop \n", time.Since(now), cnt0-1) // dbg testing.

	}

	if cnt0%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, ctOut.Value[0], ctOut.Value[0])
		ringQ.ReduceLvl(levelQ, ctOut.Value[1], ctOut.Value[1])
	}

	if cnt0%PiOverF != 0 {
		ringP.ReduceLvl(levelP, c0OutQP.P, c0OutQP.P)
		ringP.ReduceLvl(levelP, c1OutQP.P, c1OutQP.P)
	}

	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ctOut.Value[0], c0OutQP.P, ctOut.Value[0]) // sum(phi(c0 * P + d0_QP))/P
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ctOut.Value[1], c1OutQP.P, ctOut.Value[1]) // sum(phi(d1_QP))/P

	ctInRotQP = nil
	runtime.GC()

}

// Same functionality as MultiplyByDiagMatrixBSGS4ArithmeticSeq, but this routine can accepts multiple LinearTrnsformation. User must ensure that
// All the input LinearTransform instances can share the group of inner loop pre-rotated ciphertexts. Also, user must ensure that all the LinearTransform instances has the same level and scale.
func (eval *evaluator) MultiplyByDiagMatricesBSGS4ArithmeticSeq(ctIn *rlwe.Ciphertext, matrices []LinearTransform, PoolDecompQP []ringqp.Poly, ctOut []*rlwe.Ciphertext) {
	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	levelQ := utils.MinInt(ctOut[0].Level(), utils.MinInt(ctIn.Level(), matrices[0].Level))
	levelP := len(ringP.Modulus) - 1

	for i := range ctOut {
		ctOut[i].Resize(ctOut[i].Degree(), levelQ)
	}

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	// Computes the N2 rotations indexes of the non-zero rows of the diagonalized DFT matrix for the baby-step giant-step algorithm
	indexes := make([]map[int][]int, len(matrices))
	rotN2set := make(map[int]bool)
	var rotN2 []int
	for i, matrix := range matrices {
		indexes[i], _, rotN2 = BsgsIndex4ArithmeticSeq(matrix.Vec, 1<<matrix.LogSlots, matrix.N1, matrix.CommonDiff)
		for rot := range rotN2 {
			rotN2set[rot] = true
		}
	}
	rotN2 = make([]int, 0, len(rotN2set))
	for key := range rotN2set {
		rotN2 = append(rotN2, key)
	}

	ring.CopyLvl(levelQ, ctIn.Value[0], eval.buffCt.Value[0])
	ring.CopyLvl(levelQ, ctIn.Value[1], eval.buffCt.Value[1])

	ctInTmp0, ctInTmp1 := eval.buffCt.Value[0], eval.buffCt.Value[1]

	// Pre-rotates ciphertext for the baby-step giant-step algorithm, does not divide by P yet
	ctInRotQP := eval.RotateHoistedNoModDownNew(levelQ, rotN2, ctInTmp0, eval.BuffDecompQP)

	// dbg testing:

	var bytes int
	for _, ct := range ctInRotQP {
		bytes += ct.MarshalBinarySize()
	}
	fmt.Printf("\nmiddle ciphertexts take up %d bytes\n", bytes)

	ringQ.MulScalarBigintLvl(levelQ, ctInTmp0, ringP.ModulusAtLevel[levelP], ctInTmp0) // P*c0
	ringQ.MulScalarBigintLvl(levelQ, ctInTmp1, ringP.ModulusAtLevel[levelP], ctInTmp1) // P*c1

	for k, matrix := range matrices {

		// Accumulator inner loop
		tmp0QP := eval.BuffQP[1]
		tmp1QP := eval.BuffQP[2]

		// Accumulator outer loop
		cQP := rlwe.CiphertextQP{Value: [2]ringqp.Poly{eval.BuffQP[3], eval.BuffQP[4]}}
		cQP.IsNTT = true

		// Result in QP
		c0OutQP := ringqp.Poly{Q: ctOut[k].Value[0], P: eval.BuffQP[5].Q}
		c1OutQP := ringqp.Poly{Q: ctOut[k].Value[1], P: eval.BuffQP[5].P}

		// OUTER LOOP
		var cnt0 int
		for j := range indexes[k] {
			//now = time.Now()
			// INNER LOOP
			var cnt1 int
			for _, i := range indexes[k][j] {
				if i == 0 {
					if cnt1 == 0 {
						ringQ.MulCoeffsMontgomeryConstantLvl(levelQ, matrix.Vec[j].Q, ctInTmp0, tmp0QP.Q)
						ringQ.MulCoeffsMontgomeryConstantLvl(levelQ, matrix.Vec[j].Q, ctInTmp1, tmp1QP.Q)
						tmp0QP.P.Zero()
						tmp1QP.P.Zero()
					} else {
						ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, matrix.Vec[j].Q, ctInTmp0, tmp0QP.Q)
						ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, matrix.Vec[j].Q, ctInTmp1, tmp1QP.Q)
					}
				} else {
					if cnt1 == 0 {
						ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, matrix.Vec[j+i], ctInRotQP[i].Value[0], tmp0QP)
						ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, matrix.Vec[j+i], ctInRotQP[i].Value[1], tmp1QP)
					} else {
						ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, matrix.Vec[j+i], ctInRotQP[i].Value[0], tmp0QP)
						ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, matrix.Vec[j+i], ctInRotQP[i].Value[1], tmp1QP)
					}
				}

				if cnt1%QiOverF == QiOverF-1 {
					ringQ.ReduceLvl(levelQ, tmp0QP.Q, tmp0QP.Q)
					ringQ.ReduceLvl(levelQ, tmp1QP.Q, tmp1QP.Q)
				}

				if cnt1%PiOverF == PiOverF-1 {
					ringP.ReduceLvl(levelP, tmp0QP.P, tmp0QP.P)
					ringP.ReduceLvl(levelP, tmp1QP.P, tmp1QP.P)
				}

				cnt1++
				//fmt.Printf("	%s consumed for %dth outer loop's %dth inner loop\n", time.Since(now), cnt0, cnt1-1)
			}

			if cnt1%QiOverF != 0 {
				ringQ.ReduceLvl(levelQ, tmp0QP.Q, tmp0QP.Q)
				ringQ.ReduceLvl(levelQ, tmp1QP.Q, tmp1QP.Q)
			}

			if cnt1%PiOverF != 0 {
				ringP.ReduceLvl(levelP, tmp0QP.P, tmp0QP.P)
				ringP.ReduceLvl(levelP, tmp1QP.P, tmp1QP.P)
			}

			// If j != 0, then rotates ((tmp0QP.Q, tmp0QP.P), (tmp1QP.Q, tmp1QP.P)) by N1*j and adds the result on ((cQP.Value[0].Q, cQP.Value[0].P), (cQP.Value[1].Q, cQP.Value[1].P))
			if j != 0 {

				// Hoisting of the ModDown of sum(sum(phi(d1) * plaintext))
				eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, tmp1QP.Q, tmp1QP.P, tmp1QP.Q) // c1 * plaintext + sum(phi(d1) * plaintext) + phi(c1) * plaintext mod Q

				galEl := eval.params.GaloisElementForColumnRotationBy(j)

				rtk, generated := eval.Rtks.Keys[galEl]
				if !generated {
					panic("cannot MultiplyByDiagMatrixBSGS: switching key not available")
				}

				rotIndex := eval.PermuteNTTIndex[galEl]

				eval.GadgetProductNoModDown(levelQ, tmp1QP.Q, rtk.GadgetCiphertext, cQP) // Switchkey(P*phi(tmpRes_1)) = (d0, d1) in base QP
				ringQP.AddLvl(levelQ, levelP, cQP.Value[0], tmp0QP, cQP.Value[0])

				// Outer loop rotations
				if cnt0 == 0 {

					ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, cQP.Value[0], rotIndex, c0OutQP)
					ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, cQP.Value[1], rotIndex, c1OutQP)
				} else {
					ringQP.PermuteNTTWithIndexAndAddNoModLvl(levelQ, levelP, cQP.Value[0], rotIndex, c0OutQP)
					ringQP.PermuteNTTWithIndexAndAddNoModLvl(levelQ, levelP, cQP.Value[1], rotIndex, c1OutQP)
				}

				// Else directly adds on ((cQP.Value[0].Q, cQP.Value[0].P), (cQP.Value[1].Q, cQP.Value[1].P))
			} else {
				if cnt0 == 0 {
					ringQP.CopyLvl(levelQ, levelP, tmp0QP, c0OutQP)
					ringQP.CopyLvl(levelQ, levelP, tmp1QP, c1OutQP)
				} else {
					ringQP.AddNoModLvl(levelQ, levelP, c0OutQP, tmp0QP, c0OutQP)
					ringQP.AddNoModLvl(levelQ, levelP, c1OutQP, tmp1QP, c1OutQP)
				}
			}

			if cnt0%QiOverF == QiOverF-1 {
				ringQ.ReduceLvl(levelQ, ctOut[k].Value[0], ctOut[k].Value[0])
				ringQ.ReduceLvl(levelQ, ctOut[k].Value[1], ctOut[k].Value[1])
			}

			if cnt0%PiOverF == PiOverF-1 {
				ringP.ReduceLvl(levelP, c0OutQP.P, c0OutQP.P)
				ringP.ReduceLvl(levelP, c1OutQP.P, c1OutQP.P)
			}

			cnt0++
			//fmt.Printf("%s consumed for %dth outer loop \n", time.Since(now), cnt0-1) // dbg testing.

		}

		if cnt0%QiOverF != 0 {
			ringQ.ReduceLvl(levelQ, ctOut[k].Value[0], ctOut[k].Value[0])
			ringQ.ReduceLvl(levelQ, ctOut[k].Value[1], ctOut[k].Value[1])
		}

		if cnt0%PiOverF != 0 {
			ringP.ReduceLvl(levelP, c0OutQP.P, c0OutQP.P)
			ringP.ReduceLvl(levelP, c1OutQP.P, c1OutQP.P)
		}

		eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ctOut[k].Value[0], c0OutQP.P, ctOut[k].Value[0]) // sum(phi(c0 * P + d0_QP))/P
		eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ctOut[k].Value[1], c1OutQP.P, ctOut[k].Value[1]) // sum(phi(d1_QP))/P
	}
	ctInRotQP = nil
	runtime.GC()

}

// if shareInnerLoop is set to true, then the LinearTransform instances will share the pre-rotated ciphertexts for their innerLoop, this require the LinearTransform has same level and scale for correctness concern
// the shareInnerLoop mode is now only available for LinearTransform instances in ArithmeticSeq mode.
func (eval *evaluator) LinearTransform4ArithmeticSeqNew(ctIn *rlwe.Ciphertext, linearTransform interface{}, shareInnerLoop ...bool) (ctOut []*rlwe.Ciphertext) {
	switch LTs := linearTransform.(type) {
	case []LinearTransform:
		ctOut = make([]*rlwe.Ciphertext, len(LTs))

		var maxLevel int
		for _, LT := range LTs {
			maxLevel = utils.MaxInt(maxLevel, LT.Level)
		}

		minLevel := utils.MinInt(maxLevel, ctIn.Level())

		// now := time.Now() // dbg testing
		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)
		// fmt.Printf("Ciphertext Decomposition consumes %s\n", time.Since(now)) // dbg testing

		if len(shareInnerLoop) > 0 && shareInnerLoop[0] {
			for i := range LTs {
				ctOut[i] = NewCiphertext(eval.params, 1, minLevel)
			}
			eval.MultiplyByDiagMatricesBSGS4ArithmeticSeq(ctIn, LTs, eval.BuffDecompQP, ctOut)
			for i, LT := range LTs {
				ctOut[i].Scale = ctIn.Scale.Mul(LT.Scale)
			}
		} else {

			for i, LT := range LTs {
				ctOut[i] = NewCiphertext(eval.params, 1, minLevel)

				if LT.N1 == 0 {
					eval.MultiplyByDiagMatrix(ctIn, LT, eval.BuffDecompQP, ctOut[i])
				} else if LT.N1 != 0 && LT.CommonDiff == 0 {
					eval.MultiplyByDiagMatrixBSGS(ctIn, LT, eval.BuffDecompQP, ctOut[i])
				} else if LT.N1 != 0 && LT.CommonDiff != 0 {
					eval.MultiplyByDiagMatrixBSGS4ArithmeticSeq(ctIn, LT, eval.BuffDecompQP, ctOut[i])
				}

				ctOut[i].MetaData = ctIn.MetaData
				ctOut[i].Scale = ctIn.Scale.Mul(LT.Scale)
			}

		}

	case LinearTransform:

		minLevel := utils.MinInt(LTs.Level, ctIn.Level())

		// now := time.Now() // dbg testing
		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)
		// fmt.Printf("Ciphertext Decomposition consumes %s\n", time.Since(now)) // dbg testing

		ctOut = []*rlwe.Ciphertext{NewCiphertext(eval.params, 1, minLevel)}

		if LTs.N1 == 0 {
			eval.MultiplyByDiagMatrix(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		} else if LTs.N1 != 0 && LTs.CommonDiff == 0 {
			eval.MultiplyByDiagMatrixBSGS(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		} else if LTs.N1 != 0 && LTs.CommonDiff != 0 {
			eval.MultiplyByDiagMatrixBSGS4ArithmeticSeq(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		}

		ctOut[0].MetaData = ctIn.MetaData
		ctOut[0].Scale = ctIn.Scale.Mul(LTs.Scale)
	}
	return

}

func GenLinearTransformBSGS4ArithmeticSeq(encoder Encoder, value interface{}, level int, scale rlwe.Scale, N1 int, commonDiff int, logSlots int) (LT LinearTransform) {

	enc, ok := encoder.(*encoderComplex128)
	if !ok {
		panic("cannot GenLinearTransformBSGS: encoder should be an encoderComplex128")
	}

	params := enc.params

	slots := 1 << logSlots

	index, _, _ := BsgsIndex4ArithmeticSeq(value, slots, N1, commonDiff)

	vec := make(map[int]ringqp.Poly)

	dMat := interfaceMapToMapOfInterface(value)
	levelQ := level
	levelP := params.PCount() - 1

	var values interface{}
	switch value.(type) {
	case map[int][]complex128:
		values = make([]complex128, slots)
	case map[int][]float64:
		values = make([]float64, slots)
	}

	for j := range index {

		rot := -j & (slots - 1)

		for _, i := range index[j] {

			// manages inputs that have rotation between 0 and slots-1 or between -slots/2 and slots/2-1
			v, ok := dMat[j+i]
			if !ok {
				v = dMat[j+i-slots]
			}
			vec[j+i] = params.RingQP().NewPolyLvl(levelQ, levelP)

			copyRotInterface(values, v, rot)

			enc.Embed(values, logSlots, scale, true, vec[j+i])
		}
	}

	return LinearTransform{LogSlots: logSlots, N1: N1, Vec: vec, Level: level, Scale: scale, CommonDiff: commonDiff}

}

func (LT *LinearTransform) Rotations4ArithmeticSeq() (rotations []int) {
	slots := 1 << LT.LogSlots

	rotIndex := make(map[int]bool)

	N1 := LT.N1

	if LT.N1 == 0 {

		for j := range LT.Vec {
			rotIndex[j] = true
		}

	} else {
		for j := range LT.Vec {
			// see BsgsIndexArithmeticSeq for specific details.
			// The nonZeroDiags may contains rotations that has form -k*(commonDiff) mod slots for some positive integer k, scale 1/commonDiff cant directly apply on them.
			// So we should check these rotations first.
			var jDivByCD int
			if j%LT.CommonDiff != 0 { // check if rotation is divisible by commonDiff, panic when commonDiff is zero
				jDivByCD = (j - slots)              // turn the rotations with form -k*(commonDiff) mod slots back to -k*(commonDiff) < -slots.
				jDivByCD = jDivByCD / LT.CommonDiff // then we can get rotDivByCD = -k
			} else {
				jDivByCD = j / LT.CommonDiff // rotations with form k*(commonDiff) mod slots for some integer k such that -slots < k*(commonDiff) < slots can be scaled by 1/commonDiff directly.
			}
			jDivByCD &= (slots - 1)
			idxN1 := (((jDivByCD / N1) * N1) * LT.CommonDiff) & (slots - 1) // outer loop index
			idxN2 := ((jDivByCD & (N1 - 1)) * LT.CommonDiff) & (slots - 1)  // inner loop index
			rotIndex[idxN1] = true
			rotIndex[idxN2] = true

			/* original method, have some problems in some specific circumstances.
			index = ((((j / LTX.CommonDiff) / N1) * N1) * LTX.CommonDiff) & (slots - 1)
			rotIndex[index] = true

			index = (((j / LTX.CommonDiff) & (N1 - 1)) * LTX.CommonDiff) & (slots - 1)
			rotIndex[index] = true
			*/
		}
	}

	rotations = make([]int, len(rotIndex))
	var i int
	for j := range rotIndex {
		rotations[i] = j
		i++
	}

	return rotations

}

func MergeAndDistinct(a []int, b []int) (c []int) {
	set := make(map[int]bool)
	for _, value := range a {
		set[value] = true
	}
	for _, value := range b {
		set[value] = true
	}
	c = make([]int, 0, len(set))
	for key := range set {
		c = append(c, key)
	}
	return
}

// Do LTs and hoisted rotations on the same cihphertext, this might save a lot of time doing NTTDecomposition.
// Notice that this function only support naive type of LinearTransformation implementation and do not support BSGS optimization.
// Also, the minLevel among the ctIn and LTs should be the level of ctIn.
func (eval *evaluator) ColShiftRestricted(ctIn *rlwe.Ciphertext, LTs []LinearTransform, AvailableRots int, dimension int) (ctOut []*rlwe.Ciphertext) {
	AvailableKeys := AvailableRots
	TargetRots := len(LTs)
	d := dimension

	ctOut = make([]*rlwe.Ciphertext, len(LTs))

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	levelQ := ctIn.Level() // this should be the minlevel, otherwise this might go wrong.
	levelP := len(ringP.Modulus) - 1

	for i := range ctOut {
		ctOut[i] = NewCiphertext(eval.params, 1, levelQ)
		ctOut[i].Resize(ctOut[i].Degree(), levelQ)
	}

	QiOverF := eval.params.QiOverflowMargin(levelQ)
	PiOverF := eval.params.PiOverflowMargin(levelP)

	ct0TimesP := eval.BuffQP[0].Q // ct0 * P mod Q
	tmp0QP := eval.BuffQP[1]
	tmp1QP := eval.BuffQP[2]
	ksRes0QP := eval.BuffQP[3]
	ksRes1QP := eval.BuffQP[4]

	tmp1BuffDecompQP := make([]ringqp.Poly, len(eval.BuffDecompQP))
	tmp2BuffDecompQP := make([]ringqp.Poly, len(eval.BuffDecompQP))
	for i, poly := range eval.BuffDecompQP {
		tmp1BuffDecompQP[i] = poly.CopyNew()
		tmp2BuffDecompQP[i] = poly.CopyNew()
	}
	var tmpBuffDecompQP []ringqp.Poly

	var ctInBase *rlwe.Ciphertext
	var ctInNegBase *rlwe.Ciphertext
	ctInBase = ctIn.CopyNew()
	ctInNegBase = ctIn.CopyNew()

	eval.DecomposeNTT(levelQ, eval.params.PCount()-1, eval.params.PCount(), ctInBase.Value[1], ctInBase.IsNTT, tmp1BuffDecompQP)
	eval.AutomorphismHoisted(levelQ, ctInBase, tmp1BuffDecompQP, eval.params.GaloisElementForColumnRotationBy(-d+d*d), ctInNegBase)
	eval.DecomposeNTT(levelQ, eval.params.PCount()-1, eval.params.PCount(), ctInNegBase.Value[1], ctInNegBase.IsNTT, tmp2BuffDecompQP)

	var state bool
	var cnt int

	// AvailableRots should be divisible by TargetRots. (since this is a incomplete version)
	for t := 0; t < TargetRots; t += AvailableKeys {
		var end int
		if t+AvailableKeys > TargetRots {
			end = TargetRots - t
		} else {
			end = AvailableKeys
		}

		for i, matrix := range LTs[t : t+end] {
			c0OutQP := ringqp.Poly{Q: ctOut[i+t].Value[0], P: eval.BuffQP[5].Q}
			c1OutQP := ringqp.Poly{Q: ctOut[i+t].Value[1], P: eval.BuffQP[5].P}

			var ctInTmp0, ctInTmp1 *ring.Poly
			state = false
			cnt = 0

			for k := range matrix.Vec {
				if 0 <= k && k < d {
					ring.CopyLvl(levelQ, ctInBase.Value[0], eval.buffCt.Value[0])
					ring.CopyLvl(levelQ, ctInBase.Value[1], eval.buffCt.Value[1])
					ctInTmp0, ctInTmp1 = eval.buffCt.Value[0], eval.buffCt.Value[1]
					tmpBuffDecompQP = tmp1BuffDecompQP
					//eval.BuffDecompQP = tmp1BuffDecompQP
				} else {
					ring.CopyLvl(levelQ, ctInNegBase.Value[0], eval.buffCt.Value[0])
					ring.CopyLvl(levelQ, ctInNegBase.Value[1], eval.buffCt.Value[1])
					ctInTmp0, ctInTmp1 = eval.buffCt.Value[0], eval.buffCt.Value[1]
					tmpBuffDecompQP = tmp2BuffDecompQP // testing
					//eval.BuffDecompQP = tmp2BuffDecompQP // testing
				}

				ringQ.MulScalarBigintLvl(levelQ, ctInTmp0, ringP.ModulusAtLevel[levelP], ct0TimesP) // P*c0

				if k == 0 {
					state = true
				} else {
					var tmpk int
					if 0 <= k && k < d {
						tmpk = k % AvailableKeys
						if tmpk == 0 && k != 0 {
							tmpk = AvailableKeys
						}
					} else {
						tmpk = ((k + d + d*d) % (d * d))
						if tmpk%AvailableKeys == 0 && tmpk != 0 {
							tmpk = AvailableKeys
						} else {
							tmpk = tmpk % AvailableKeys
						}
					}
					galEl := eval.params.GaloisElementForColumnRotationBy(tmpk)
					rtk, generated := eval.Rtks.Keys[galEl]
					if !generated {
						panic("cannot MultiplyByDiagMatrix: switching key not available")
					}
					index := eval.PermuteNTTIndex[galEl]

					eval.KeyswitchHoistedNoModDown(levelQ, tmpBuffDecompQP, rtk, ksRes0QP.Q, ksRes1QP.Q, ksRes0QP.P, ksRes1QP.P)
					ringQ.AddLvl(levelQ, ksRes0QP.Q, ct0TimesP, ksRes0QP.Q)
					ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, ksRes0QP, index, tmp0QP)
					ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, ksRes1QP, index, tmp1QP)

					if cnt == 0 {
						// keyswitch(c1_Q) = (d0_QP, d1_QP)
						ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, matrix.Vec[k], tmp0QP, c0OutQP)
						ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, matrix.Vec[k], tmp1QP, c1OutQP)
					} else {
						// keyswitch(c1_Q) = (d0_QP, d1_QP)
						ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, matrix.Vec[k], tmp0QP, c0OutQP)
						ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, matrix.Vec[k], tmp1QP, c1OutQP)
					}

					if cnt%QiOverF == QiOverF-1 {
						ringQ.ReduceLvl(levelQ, c0OutQP.Q, c0OutQP.Q)
						ringQ.ReduceLvl(levelQ, c1OutQP.Q, c1OutQP.Q)
					}

					if cnt%PiOverF == PiOverF-1 {
						ringP.ReduceLvl(levelP, c0OutQP.P, c0OutQP.P)
						ringP.ReduceLvl(levelP, c1OutQP.P, c1OutQP.P)
					}

					cnt++

				}

			}
			if cnt%QiOverF == 0 {
				ringQ.ReduceLvl(levelQ, c0OutQP.Q, c0OutQP.Q)
				ringQ.ReduceLvl(levelQ, c1OutQP.Q, c1OutQP.Q)
			}

			if cnt%PiOverF == 0 {
				ringP.ReduceLvl(levelP, c0OutQP.P, c0OutQP.P)
				ringP.ReduceLvl(levelP, c1OutQP.P, c1OutQP.P)
			}

			eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c0OutQP.Q, c0OutQP.P, c0OutQP.Q) // sum(phi(c0 * P + d0_QP))/P
			eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c1OutQP.Q, c1OutQP.P, c1OutQP.Q) // sum(phi(d1_QP))/P

			if state { // Rotation by zero
				ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0].Q, ctInTmp0, c0OutQP.Q) // ctOut += c0_Q * plaintext
				ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0].Q, ctInTmp1, c1OutQP.Q) // ctOut += c1_Q * plaintext
			}
			ctOut[i+t].MetaData = ctIn.MetaData
			ctOut[i+t].Scale = ctIn.Scale.Mul(matrix.Scale)

		}

		eval.AutomorphismHoisted(levelQ, ctInBase, tmp1BuffDecompQP, eval.params.GaloisElementForColumnRotationBy(AvailableKeys), ctInBase)
		eval.AutomorphismHoisted(levelQ, ctInNegBase, tmp2BuffDecompQP, eval.params.GaloisElementForColumnRotationBy(AvailableKeys), ctInNegBase)
		eval.DecomposeNTT(levelQ, eval.params.PCount()-1, eval.params.PCount(), ctInBase.Value[1], ctInBase.IsNTT, tmp1BuffDecompQP)
		eval.DecomposeNTT(levelQ, eval.params.PCount()-1, eval.params.PCount(), ctInNegBase.Value[1], ctInNegBase.IsNTT, tmp2BuffDecompQP)
	}

	return
}
