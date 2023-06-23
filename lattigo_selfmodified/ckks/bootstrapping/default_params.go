package bootstrapping

import (
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/advanced"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

type defaultParametersLiteral struct {
	SchemeParams        ckks.ParametersLiteral
	BootstrappingParams Parameters
}

// The parameters provided hereunder are the parameters used in the paper
// Bootstrapping for Approximate Homomorphic Encryption with Negligible
// Failure-Probability by Using Sparse-Secret Encapsulation,
// https://eprint.iacr.org/2022/024

// DefaultParametersSparse is a set of default bootstrapping parameters with H=192 as main secret and H=32 as ephemeral secret.
var DefaultParametersSparse = []defaultParametersLiteral{N16QP1546H192H32, N16QP1547H192H32, N16QP1553H192H32, N15QP768H192H32}

// DefaultParametersDense is a set of default bootstrapping parameters with H=N/2 as main secret and H=32 as ephemeral secret.
var DefaultParametersDense = []defaultParametersLiteral{N16QP1767H32768H32, N16QP1788H32768H32, N16QP1793H32768H32, N15QP880H16384H32}

var (
	// N16QP1546H192H32 is a default bootstrapping parameters for a main secret with H=192 and an ephemeral secret with H=32.
	// Residual Q : 420 bits.
	// Precision : 26.6 bits for 2^{15} slots.
	// Failure : 2^{-138.7} for 2^{15} slots.
	N16QP1546H192H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:  16,
			Sigma: rlwe.DefaultSigma,
			H:     192,
			Q: []uint64{
				0x10000000006e0001, // 60 Q0
				0x10000140001,      // 40
				0xffffe80001,       // 40
				0xffffc40001,       // 40
				0x100003e0001,      // 40
				0xffffb20001,       // 40
				0x10000500001,      // 40
				0xffff940001,       // 40
				0xffff8a0001,       // 40
				0xffff820001,       // 40
				0x7fffe60001,       // 39 StC
				0x7fffe40001,       // 39 StC
				0x7fffe00001,       // 39 StC
				0xfffffffff840001,  // 60 Sine (double angle)
				0x1000000000860001, // 60 Sine (double angle)
				0xfffffffff6a0001,  // 60 Sine
				0x1000000000980001, // 60 Sine
				0xfffffffff5a0001,  // 60 Sine
				0x1000000000b00001, // 60 Sine
				0x1000000000ce0001, // 60 Sine
				0xfffffffff2a0001,  // 60 Sine
				0x100000000060001,  // 56 CtS
				0xfffffffff00001,   // 56 CtS
				0xffffffffd80001,   // 56 CtS
				0x1000000002a0001,  // 56 CtS
			},
			P: []uint64{
				0x1fffffffffe00001, // Pi 61
				0x1fffffffffc80001, // Pi 61
				0x1fffffffffb40001, // Pi 61
				0x1fffffffff500001, // Pi 61
				0x1fffffffff420001, // Pi 61
			},
			LogSlots:     15,
			DefaultScale: 1 << 40,
		},

		Parameters{
			EphemeralSecretWeight: 32,
			SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.SlotsToCoeffs,
				RepackImag2Real:     true,
				LevelStart:          12,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{0x7fffe60001},
					{0x7fffe40001},
					{0x7fffe00001},
				},
			},
			EvalModParameters: advanced.EvalModLiteral{
				Q:             0x10000000006e0001,
				LevelStart:    20,
				SineType:      advanced.Cos1,
				MessageRatio:  256.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ArcSineDeg:    0,
				ScalingFactor: 1 << 60,
			},
			CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.CoeffsToSlots,
				RepackImag2Real:     true,
				LevelStart:          24,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{0x100000000060001},
					{0xfffffffff00001},
					{0xffffffffd80001},
					{0x1000000002a0001},
				},
			},
		},
	}

	// N16QP1547H192H32 is a default bootstrapping parameters for a main secret with H=192 and an ephemeral secret with H=32.
	// Residual Q : 285 bits.
	// Precision : 32.1 bits for 2^{15} slots.
	// Failure : 2^{-138.7} for 2^{15} slots.
	N16QP1547H192H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:  16,
			Sigma: rlwe.DefaultSigma,
			H:     192,
			Q: []uint64{
				0x10000000006e0001, // 60 Q0
				0x2000000a0001,     // 45
				0x2000000e0001,     // 45
				0x1fffffc20001,     // 45
				0x200000440001,     // 45
				0x200000500001,     // 45
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
			},
			LogSlots:     15,
			DefaultScale: 1 << 45,
		},

		Parameters{
			EphemeralSecretWeight: 32,
			SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.SlotsToCoeffs,
				RepackImag2Real:     true,
				LevelStart:          8,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{0x3ffffe80001},
					{0x3ffffd20001},
					{0x3ffffca0001},
				},
			},
			EvalModParameters: advanced.EvalModLiteral{
				Q:             0x10000000006e0001,
				LevelStart:    19,
				SineType:      advanced.Cos1,
				MessageRatio:  4.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ArcSineDeg:    7,
				ScalingFactor: 1 << 60,
			},
			CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.CoeffsToSlots,
				RepackImag2Real:     true,
				LevelStart:          23,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{0x400000000360001},
					{0x3ffffffffbe0001},
					{0x400000000660001},
					{0x4000000008a0001},
				},
			},
		},
	}

	// N16QP1553H192H32 is a default bootstrapping parameters for a main secret with H=192 and an ephemeral secret with H=32.
	// Residual Q : 505 bits.
	// Precision : 19.1 bits for 2^{15} slots.
	// Failure : 2^{-138.7} for 2^{15} slots.
	N16QP1553H192H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:  16,
			Sigma: rlwe.DefaultSigma,
			H:     192,
			Q: []uint64{
				0x80000000080001,   // 55 Q0
				0xffffffffffc0001,  // 60
				0x10000000006e0001, // 60
				0xfffffffff840001,  // 60
				0x1000000000860001, // 60
				0xfffffffff6a0001,  // 60
				0x1000000000980001, // 60
				0xfffffffff5a0001,  // 60
				0x1000000000b00001, // 60 StC  (30)
				0x1000000000ce0001, // 60 StC  (30+30)
				0x80000000440001,   // 55 Sine (double angle)
				0x7fffffffba0001,   // 55 Sine (double angle)
				0x80000000500001,   // 55 Sine
				0x7fffffffaa0001,   // 55 Sine
				0x800000005e0001,   // 55 Sine
				0x7fffffff7e0001,   // 55 Sine
				0x7fffffff380001,   // 55 Sine
				0x80000000ca0001,   // 55 Sine
				0x200000000e0001,   // 53 CtS
				0x20000000140001,   // 53 CtS
				0x20000000280001,   // 53 CtS
				0x1fffffffd80001,   // 53 CtS
			},
			P: []uint64{
				0x1fffffffffe00001, // Pi 61
				0x1fffffffffc80001, // Pi 61
				0x1fffffffffb40001, // Pi 61
				0x1fffffffff500001, // Pi 61
				0x1fffffffff420001, // Pi 61
			},
			LogSlots:     15,
			DefaultScale: 1 << 30,
		},

		Parameters{
			EphemeralSecretWeight: 32,
			SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.SlotsToCoeffs,
				RepackImag2Real:     true,
				LevelStart:          9,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{1073741824.0},
					{1073741824.0062866, 1073741824.0062866},
				},
			},
			EvalModParameters: advanced.EvalModLiteral{
				Q:             0x80000000080001,
				LevelStart:    17,
				SineType:      advanced.Cos1,
				MessageRatio:  256.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ArcSineDeg:    0,
				ScalingFactor: 1 << 55,
			},
			CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.CoeffsToSlots,
				RepackImag2Real:     true,
				LevelStart:          21,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{0x200000000e0001},
					{0x20000000140001},
					{0x20000000280001},
					{0x1fffffffd80001},
				},
			},
		},
	}

	// N15QP768H192H32 is a default bootstrapping parameters for a main secret with H=192 and an ephemeral secret with H=32.
	// Residual Q : 110 bits.
	// Precision : 15.4 bits for 2^{14} slots.
	// Failure : 2^{-139.7} for 2^{14} slots.
	N15QP768H192H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:  15,
			Sigma: rlwe.DefaultSigma,
			H:     192,
			Q: []uint64{
				0x1fff90001,       // 32 Q0
				0x4000000420001,   // 50
				0x1fc0001,         // 25
				0xffffffffffc0001, // 60 StC (30+30)
				0x4000000120001,   // 50 Sine
				0x40000001b0001,   // 50 Sine
				0x3ffffffdf0001,   // 50 Sine
				0x4000000270001,   // 50 Sine
				0x3ffffffd20001,   // 50 Sine
				0x3ffffffcd0001,   // 50 Sine
				0x4000000350001,   // 50 Sine
				0x3ffffffc70001,   // 50 Sine
				0x1fffffff50001,   // 49 CtS
				0x1ffffffea0001,   // 49 CtS
			},
			P: []uint64{
				0x7fffffffe0001, // 51
				0x8000000110001, // 51
			},
			LogSlots:     14,
			DefaultScale: 1 << 25,
		},

		Parameters{
			EphemeralSecretWeight: 32,
			SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.SlotsToCoeffs,
				RepackImag2Real:     true,
				LevelStart:          3,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{1073741823.9998779, 1073741823.9998779},
				},
			},
			EvalModParameters: advanced.EvalModLiteral{
				Q:             0x1fff90001,
				LevelStart:    11,
				SineType:      advanced.Cos1,
				MessageRatio:  256.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ArcSineDeg:    0,
				ScalingFactor: 1 << 50,
			},
			CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.CoeffsToSlots,
				RepackImag2Real:     true,
				LevelStart:          13,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{0x1fffffff50001},
					{0x1ffffffea0001},
				},
			},
		},
	}

	// N16QP1767H32768H32 is a default bootstrapping parameters for a main secret with H=32768 and an ephemeral secret with H=32.
	// Residual Q : 580 bits.
	// Precision : 23.0 bits for 2^{15} slots.
	// Failure : 2^{-138.7} for 2^{15} slots.
	N16QP1767H32768H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:  16,
			Sigma: rlwe.DefaultSigma,
			H:     32768,
			Q: []uint64{
				0x10000000006e0001, // 60 Q0
				0x10000140001,      // 40
				0xffffe80001,       // 40
				0xffffc40001,       // 40
				0x100003e0001,      // 40
				0xffffb20001,       // 40
				0x10000500001,      // 40
				0xffff940001,       // 40
				0xffff8a0001,       // 40
				0xffff820001,       // 40
				0xffff780001,       // 40
				0x10000960001,      // 40
				0x10000a40001,      // 40
				0xffff580001,       // 40
				0x7fffe60001,       // 39 StC
				0x7fffe40001,       // 39 StC
				0x7fffe00001,       // 39 StC
				0xfffffffff840001,  // 60 Sine (double angle)
				0x1000000000860001, // 60 Sine (double angle)
				0xfffffffff6a0001,  // 60 Sine
				0x1000000000980001, // 60 Sine
				0xfffffffff5a0001,  // 60 Sine
				0x1000000000b00001, // 60 Sine
				0x1000000000ce0001, // 60 Sine
				0xfffffffff2a0001,  // 60 Sine
				0x100000000060001,  // 56 CtS
				0xfffffffff00001,   // 56 CtS
				0xffffffffd80001,   // 56 CtS
				0x1000000002a0001,  // 56 CtS
			},
			P: []uint64{
				0x1fffffffffe00001, // Pi 61
				0x1fffffffffc80001, // Pi 61
				0x1fffffffffb40001, // Pi 61
				0x1fffffffff500001, // Pi 61
				0x1fffffffff420001, // Pi 61
				0x1fffffffff380001, // Pi 61
			},
			LogSlots:     15,
			DefaultScale: 1 << 40,
		},

		Parameters{
			EphemeralSecretWeight: 32,
			SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.SlotsToCoeffs,
				RepackImag2Real:     true,
				LevelStart:          16,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{0x7fffe60001},
					{0x7fffe40001},
					{0x7fffe00001},
				},
			},
			EvalModParameters: advanced.EvalModLiteral{
				Q:             0x10000000006e0001,
				LevelStart:    24,
				SineType:      advanced.Cos1,
				MessageRatio:  256.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ArcSineDeg:    0,
				ScalingFactor: 1 << 60,
			},
			CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.CoeffsToSlots,
				RepackImag2Real:     true,
				LevelStart:          28,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{0x100000000060001},
					{0xfffffffff00001},
					{0xffffffffd80001},
					{0x1000000002a0001},
				},
			},
		},
	}

	// N16QP1788H32768H32 is a default bootstrapping parameters for a main secret with H=32768 and an ephemeral secret with H=32.
	// Residual Q : 465 bits.
	// Precision : 29.0 bits for 2^{15} slots.
	// Failure : 2^{-138.7} for 2^{15} slots.
	N16QP1788H32768H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
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
			LogSlots:     15,
			DefaultScale: 1 << 45,
		},

		Parameters{
			EphemeralSecretWeight: 32,
			SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.SlotsToCoeffs,
				RepackImag2Real:     true,
				LevelStart:          12,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{0x3ffffe80001},
					{0x3ffffd20001},
					{0x3ffffca0001},
				},
			},
			EvalModParameters: advanced.EvalModLiteral{
				Q:             0x10000000006e0001,
				LevelStart:    23,
				SineType:      advanced.Cos1,
				MessageRatio:  4.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ArcSineDeg:    7,
				ScalingFactor: 1 << 60,
			},
			CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.CoeffsToSlots,
				RepackImag2Real:     true,
				LevelStart:          27,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{0x400000000360001},
					{0x3ffffffffbe0001},
					{0x400000000660001},
					{0x4000000008a0001},
				},
			},
		},
	}

	// N16QP1793H32768H32 is a default bootstrapping parameters for a main secret with H=32768 and an ephemeral secret with H=32.
	// Residual Q : 745 bits.
	// Precision : 17.8 bits for 2^{15} slots.
	// Failure : 2^{-138.7} for 2^{15} slots.
	N16QP1793H32768H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:  16,
			Sigma: rlwe.DefaultSigma,
			H:     32768,
			Q: []uint64{
				0x80000000080001,   // 55 Q0
				0xffffffffffc0001,  // 60
				0x10000000006e0001, // 60
				0xfffffffff840001,  // 60
				0x1000000000860001, // 60
				0xfffffffff6a0001,  // 60
				0x1000000000980001, // 60
				0xfffffffff5a0001,  // 60
				0xfffffffff2a0001,  // 60
				0xfffffffff240001,  // 60
				0x1000000000f00001, // 60
				0xffffffffefe0001,  // 60
				0x1000000000b00001, // 60 StC  (30)
				0x1000000000ce0001, // 60 StC  (30+30)
				0x80000000440001,   // 55 Sine (double angle)
				0x7fffffffba0001,   // 55 Sine (double angle)
				0x80000000500001,   // 55 Sine
				0x7fffffffaa0001,   // 55 Sine
				0x800000005e0001,   // 55 Sine
				0x7fffffff7e0001,   // 55 Sine
				0x7fffffff380001,   // 55 Sine
				0x80000000ca0001,   // 55 Sine
				0x200000000e0001,   // 53 CtS
				0x20000000140001,   // 53 CtS
				0x20000000280001,   // 53 CtS
				0x1fffffffd80001,   // 53 CtS
			},
			P: []uint64{
				0x1fffffffffe00001, // Pi 61
				0x1fffffffffc80001, // Pi 61
				0x1fffffffffb40001, // Pi 61
				0x1fffffffff500001, // Pi 61
				0x1fffffffff420001, // Pi 61
			},
			LogSlots:     15,
			DefaultScale: 1 << 30,
		},

		Parameters{
			EphemeralSecretWeight: 32,
			SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.SlotsToCoeffs,
				RepackImag2Real:     true,
				LevelStart:          13,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{1073741824.0},
					{1073741824.0062866, 1073741824.0062866},
				},
			},
			EvalModParameters: advanced.EvalModLiteral{
				Q:             0x80000000080001,
				LevelStart:    21,
				SineType:      advanced.Cos1,
				MessageRatio:  256.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ArcSineDeg:    0,
				ScalingFactor: 1 << 55,
			},
			CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.CoeffsToSlots,
				RepackImag2Real:     true,
				LevelStart:          25,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{0x200000000e0001},
					{0x20000000140001},
					{0x20000000280001},
					{0x1fffffffd80001},
				},
			},
		},
	}

	// N15QP880H16384H32 is a default bootstrapping parameters for a main secret with H=16384 and an ephemeral secret with H=32.
	// Residual Q : 166 bits.
	// Precision : 17.3 bits for 2^{14} slots.
	// Failure : 2^{-139.7} for 2^{14} slots.
	N15QP880H16384H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:  15,
			Sigma: rlwe.DefaultSigma,
			H:     16384,
			Q: []uint64{
				0x10000140001,      // 40 Q0
				0x7ffe0001,         // 31
				0x7ff80001,         // 31
				0x80140001,         // 31
				0x7fea0001,         // 31
				0x1000000000ce0001, // 60 StC  (30+30)
				0x80000000080001,   // 55 Sine (double angle)
				0x80000000440001,   // 55 Sine (double angle)
				0x7fffffffba0001,   // 55 Sine
				0x80000000500001,   // 55 Sine
				0x7fffffffaa0001,   // 55 Sine
				0x800000005e0001,   // 55 Sine
				0x7fffffff7e0001,   // 55 Sine
				0x7fffffff380001,   // 55 Sine
				0x10000000060001,   // 52 CtS
				0xffffffff00001,    // 52 CtS
			},
			P: []uint64{
				0x100000000060001, // 56
				0x1000000002a0001, // 56
			},
			LogSlots:     14,
			DefaultScale: 1 << 31,
		},

		Parameters{
			EphemeralSecretWeight: 32,
			SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.SlotsToCoeffs,
				RepackImag2Real:     true,
				LevelStart:          5,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{1073741824.0062866, 1073741824.0062866},
				},
			},
			EvalModParameters: advanced.EvalModLiteral{
				Q:             0x10000140001,
				LevelStart:    13,
				SineType:      advanced.Cos1,
				MessageRatio:  256.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ArcSineDeg:    0,
				ScalingFactor: 1 << 55,
			},
			CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.CoeffsToSlots,
				RepackImag2Real:     true,
				LevelStart:          15,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{0x10000000060001},
					{0xffffffff00001},
				},
			},
		},
	}
)
