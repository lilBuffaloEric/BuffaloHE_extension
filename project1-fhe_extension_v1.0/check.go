package main

import "project1-fhe_extension_v1.0/check"

func main() {
	// ------------------------------ Internal Components Check (Not all available) ------------------------------
	// check.CovMatrix_check()
	// check.MatrixMultX_check(true)
	// check.LinTransX_check()
	// check.LinTransX2_check()
	// check.LinTransX3_check()
	// check.Transepose_check()
	// check.OnlyMatrixMult_check()
	// check.Invsqrt_check() // checked
	// check.SquareMatrix_product_check() // checked
	// check.SquareMatrix_encode_check() // checked
	// check.SquareMatrix_permute_check() // checked
	// check.Piecewise_check() // checked
	// check.Softmax_check() // checked
	// check.Csvio_check()
	// check.PCA_check(false) // 202356: update N1 as 32 for all Lintrans.
	// check.PCA_largeSpace_check(true)
	// check.CiphertextIO_check()
	// check.Matrixvisualize_check()
	// check.DiskIO_check()
	// check.Correctness_check()
	// check.ErrorCtrlX_check(true)
	// check.StoreAndLoadKeys_btpVersion_check(true)

	// ------------------------------ Available for testing ------------------------------
	// check.LinTransX3_check() // Routine Testing the Diagonal Convergence Decomposition of the Lintrans in homomorphic matrix multiplication.
	// check.PCA_btpVersion_check(false, false) // Routine Computing Cov' = 1/N*X^T*X and U = u*u^T
	check.PCA_largeSpace_Ver2_check(false) // Routine Computing Cov' = 1/N*X^T*X and U = u*u^T (theoratically faster, but need more space to cache intermediate homomorphic matrix-mult results in memory)
	// check.CovMatrix_btpVersion_check(false, false) // Routine Computing Cov = Cov' - U and the approximate eigenvectors of Cov. Use bootstrpping if reencryptionMode is set to "Ture", but correctness cannot be guaranteed.
	// check.PandaPCA_check(true) // Routine simulating the Panda's scheme (2021). This can only be reguarded as a speed testing, and no correctness can be guaranteed.
}


