# BuffaloHE_extension
A personal extension of Lattigo Library. ¯\\_(ツ)_/¯


Built based on the library Lattigo, this extension aims to construct more complex and interesting homomorphic functions using api
provided by Lattigo. The extension is now still being developed and is expected to contain functionality such
as (but not limited to): Evaluation of non-polynomial function inverse square root and softmax, remote sampling on a cloud server, Matrix 
multiplication for different shapes, and shallow neural network training routine.

We now mainly focus on a privacy-preserving PCA protocol for cloud computing scenarios. To run the codes, one can do the following:
- `cd project1-fhe_extension_v1.0`
- `mk vendor`
- `go mod tidy` (draw all dependencies, they will be automatically loaded into `vendor`)
- Replace the specific files in `vendor/github.com/tuneinsight/lattigo/v4` with the files in `~/BuffaloHE_extension/lattigo_selfmodified`
- Unzip the `tailored_datasets.rar` and put the `.csv` files under the `project1-fhe_extension_v1.0` directory.
- Uncomment routines in `check.go`, do `go build check.go` and `go run check.go` to execute expected tests.

More detailed comments are coming soon...
