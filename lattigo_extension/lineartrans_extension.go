package lattigoextension

import (
	"runtime"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

type CKKSevaluatorX struct {
	Eval   ckks.Evaluator
	Params ckks.Parameters
}

type CKKSLinearTransformX struct {
	LT         ckks.LinearTransform
	CommonDiff int
}

func (evaluatorX *CKKSevaluatorX) MultiplyByDiagMatrixBSGS4ArithmeticSeq(ctIn *rlwe.Ciphertext, matrix CKKSLinearTransformX, PoolDecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext) {
	params := evaluatorX.Eval.GetRLWEEvaluator().Parameters()
	eval := evaluatorX.Eval.GetRLWEEvaluator()

	ringQ := params.RingQ()
	ringP := params.RingP()
	ringQP := params.RingQP()

	levelQ := utils.MinInt(ctOut.Level(), utils.MinInt(ctIn.Level(), matrix.LT.Level))
	levelP := len(ringP.Modulus) - 1

	ctOut.Resize(ctOut.Degree(), levelQ)

	QiOverF := params.QiOverflowMargin(levelQ) >> 1
	PiOverF := params.PiOverflowMargin(levelP) >> 1

	// Computes the N2 rotations indexes of the non-zero rows of the diagonalized DFT matrix for the baby-step giant-step algorithm
	index, _, rotN2 := BsgsIndex4ArithmeticSeq(matrix.LT.Vec, 1<<matrix.LT.LogSlots, matrix.LT.N1, matrix.CommonDiff)

	ring.CopyLvl(levelQ, ctIn.Value[0], evaluatorX.Eval.BuffCt().Value[0])
	ring.CopyLvl(levelQ, ctIn.Value[1], evaluatorX.Eval.BuffCt().Value[1])

	ctInTmp0, ctInTmp1 := eval.BuffCt.Value[0], eval.BuffCt.Value[1]

	// Pre-rotates ciphertext for the baby-step giant-step algorithm, does not divide by P yet
	ctInRotQP := evaluatorX.Eval.RotateHoistedNoModDownNew(levelQ, rotN2, ctInTmp0, eval.BuffDecompQP)

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
					ringQ.MulCoeffsMontgomeryConstantLvl(levelQ, matrix.LT.Vec[j].Q, ctInTmp0, tmp0QP.Q)
					ringQ.MulCoeffsMontgomeryConstantLvl(levelQ, matrix.LT.Vec[j].Q, ctInTmp1, tmp1QP.Q)
					tmp0QP.P.Zero()
					tmp1QP.P.Zero()
				} else {
					ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, matrix.LT.Vec[j].Q, ctInTmp0, tmp0QP.Q)
					ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, matrix.LT.Vec[j].Q, ctInTmp1, tmp1QP.Q)
				}
			} else {
				if cnt1 == 0 {
					ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, matrix.LT.Vec[j+i], ctInRotQP[i].Value[0], tmp0QP)
					ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, matrix.LT.Vec[j+i], ctInRotQP[i].Value[1], tmp1QP)
				} else {
					ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, matrix.LT.Vec[j+i], ctInRotQP[i].Value[0], tmp0QP)
					ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, matrix.LT.Vec[j+i], ctInRotQP[i].Value[1], tmp1QP)
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

			galEl := params.GaloisElementForColumnRotationBy(j)

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

func (evaluatorX *CKKSevaluatorX) LinearTransformNew(ctIn *rlwe.Ciphertext, linearTransform interface{}) (ctOut []*rlwe.Ciphertext) {
	eval := evaluatorX.Eval.GetRLWEEvaluator()
	params := evaluatorX.Eval.GetRLWEEvaluator().Parameters()
	params_ckks := evaluatorX.Params
	switch LTs := linearTransform.(type) {
	/*
		case []ckks.LinearTransform:

			ctOut = make([]*rlwe.Ciphertext, len(LTs))

			var maxLevel int
			for _, LT := range LTs {
				maxLevel = utils.MaxInt(maxLevel, LT.Level)
			}

			minLevel := utils.MinInt(maxLevel, ctIn.Level())
			eval.DecomposeNTT(minLevel, params.PCount()-1, params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)

			for i, LT := range LTs {
				ctOut[i] = ckks.NewCiphertext(params_ckks, 1, minLevel)

				if LT.N1 == 0 {
					evaluatorX.Eval.MultiplyByDiagMatrix(ctIn, LT, eval.BuffDecompQP, ctOut[i])
				} else {
					evaluatorX.Eval.MultiplyByDiagMatrixBSGS(ctIn, LT, eval.BuffDecompQP, ctOut[i])
				}

				ctOut[i].MetaData = ctIn.MetaData
				ctOut[i].Scale = ctIn.Scale.Mul(LT.Scale)
			}

		case ckks.LinearTransform:

			minLevel := utils.MinInt(LTs.Level, ctIn.Level())
			eval.DecomposeNTT(minLevel, params.PCount()-1, params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)

			ctOut = []*rlwe.Ciphertext{ckks.NewCiphertext(params_ckks, 1, minLevel)}

			if LTs.N1 == 0 {
				evaluatorX.Eval.MultiplyByDiagMatrix(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
			} else {
				evaluatorX.Eval.MultiplyByDiagMatrixBSGS(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
			}

			ctOut[0].MetaData = ctIn.MetaData
			ctOut[0].Scale = ctIn.Scale.Mul(LTs.Scale)
	*/
	case []ckks.LinearTransform:
		ctOut = evaluatorX.Eval.LinearTransformNew(ctIn, LTs)
	case ckks.LinearTransform:
		ctOut = evaluatorX.Eval.LinearTransformNew(ctIn, LTs)
	case []CKKSLinearTransformX:
		ctOut = make([]*rlwe.Ciphertext, len(LTs))

		var maxLevel int
		for _, LTX := range LTs {
			maxLevel = utils.MaxInt(maxLevel, LTX.LT.Level)
		}

		minLevel := utils.MinInt(maxLevel, ctIn.Level())
		eval.DecomposeNTT(minLevel, params.PCount()-1, params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)

		for i, LTX := range LTs {
			ctOut[i] = ckks.NewCiphertext(params_ckks, 1, minLevel)

			if LTX.LT.N1 == 0 && LTX.CommonDiff == 0 {
				evaluatorX.Eval.MultiplyByDiagMatrix(ctIn, LTX.LT, eval.BuffDecompQP, ctOut[i])
			} else if LTX.LT.N1 != 0 && LTX.CommonDiff == 0 {
				evaluatorX.Eval.MultiplyByDiagMatrixBSGS(ctIn, LTX.LT, eval.BuffDecompQP, ctOut[i])
			} else if LTX.LT.N1 != 0 && LTX.CommonDiff != 0 {
				evaluatorX.MultiplyByDiagMatrixBSGS4ArithmeticSeq(ctIn, LTX, eval.BuffDecompQP, ctOut[i])
			} else {
				// constructing....
			}

			ctOut[i].MetaData = ctIn.MetaData
			ctOut[i].Scale = ctIn.Scale.Mul(LTX.LT.Scale)
		}

	case CKKSLinearTransformX:
		minLevel := utils.MinInt(LTs.LT.Level, ctIn.Level())
		eval.DecomposeNTT(minLevel, params.PCount()-1, params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)

		ctOut = []*rlwe.Ciphertext{ckks.NewCiphertext(params_ckks, 1, minLevel)}

		if LTs.LT.N1 == 0 && LTs.CommonDiff == 0 {
			evaluatorX.Eval.MultiplyByDiagMatrix(ctIn, LTs.LT, eval.BuffDecompQP, ctOut[0])
		} else if LTs.LT.N1 != 0 && LTs.CommonDiff == 0 {
			evaluatorX.Eval.MultiplyByDiagMatrixBSGS(ctIn, LTs.LT, eval.BuffDecompQP, ctOut[0])
		} else if LTs.LT.N1 != 0 && LTs.CommonDiff != 0 {
			evaluatorX.MultiplyByDiagMatrixBSGS4ArithmeticSeq(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		} else {
			// constructing...
		}

		ctOut[0].MetaData = ctIn.MetaData
		ctOut[0].Scale = ctIn.Scale.Mul(LTs.LT.Scale)

	}

	return
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

func GenLinearTransformBSGS4ArithmeticSeq(params ckks.Parameters, value interface{}, level int, scale rlwe.Scale, N1 int, commonDiff int, logSlots int) (LT CKKSLinearTransformX) {

	encoder := ckks.NewEncoder(params)

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

			encoder.Embed(values, logSlots, scale, true, vec[j+i])
		}
	}

	return CKKSLinearTransformX{LT: ckks.LinearTransform{LogSlots: logSlots, N1: N1, Vec: vec, Level: level, Scale: scale}, CommonDiff: commonDiff}

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

func (LTX *CKKSLinearTransformX) Rotations() (rotations []int) {
	slots := 1 << LTX.LT.LogSlots

	rotIndex := make(map[int]bool)

	// var index int

	N1 := LTX.LT.N1

	if LTX.LT.N1 == 0 {

		for j := range LTX.LT.Vec {
			rotIndex[j] = true
		}

	} else {

		for j := range LTX.LT.Vec {
			// see BsgsIndexArithmeticSeq for specific details.
			// The nonZeroDiags may contains rotations that has form -k*(commonDiff) mod slots for some positive integer k, scale 1/commonDiff cant directly apply on them.
			// So we should check these rotations first.
			var jDivByCD int
			if j%LTX.CommonDiff != 0 { // check if rotation is divisible by commonDiff, panic when commonDiff is zero
				jDivByCD = (j - slots)               // turn the rotations with form -k*(commonDiff) mod slots back to -k*(commonDiff) < -slots.
				jDivByCD = jDivByCD / LTX.CommonDiff // then we can get rotDivByCD = -k
			} else {
				jDivByCD = j / LTX.CommonDiff // rotations with form k*(commonDiff) mod slots for some integer k such that -slots < k*(commonDiff) < slots can be scaled by 1/commonDiff directly.
			}
			jDivByCD &= (slots - 1)
			idxN1 := (((jDivByCD / N1) * N1) * LTX.CommonDiff) & (slots - 1) // outer loop index
			idxN2 := ((jDivByCD & (N1 - 1)) * LTX.CommonDiff) & (slots - 1)  // inner loop index
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
