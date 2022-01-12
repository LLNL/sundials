
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_memory.h>

#ifdef SUNDIALS_SUNMATRIX_MAGMADENSE
#include <sunmatrix/sunmatrix_magmadense.h>
#endif
#ifdef BUILD_SUNMATRIX_ONEMKLDENSE
#include <sunmatrix/sunmatrix_onemkldense.h>
#endif

namespace sundials
{

struct _SUNExecQueue
{ };

typedef struct _SUNExecQueue *SUNExecQueue;

template<SUNMatrix_ID ID>
class BlockMatrix
{
public:
  BlockMatrix(SUNContext sunctx = nullptr)
    :  A_(SUNMatNewEmpty(sunctx)),
       isOwner_(true)
  { }

  BlockMatrix(SUNMatrix A) : A_(A), isOwner_(false) { }

  virtual ~BlockMatrix() { if (isOwner_) SUNMatDestroy(A_); }

  operator SUNMatrix()
  {
    return A_;
  }

  virtual sunindextype BlockRows() = 0;
  virtual sunindextype BlockCols() = 0;
  virtual sunindextype Rows() = 0;
  virtual sunindextype Cols() = 0;
  virtual sunindextype NumBlocks() = 0;
  virtual realtype* Data() = 0;
  virtual realtype** BlockData() = 0;


  // SUNMatrix unary operations
  int Zero() { return SUNMatZero(A_); }

  // SUNMatrix binary operations
  static int Copy(BlockMatrix& A, BlockMatrix& B) { return SUNMatCopy(A, B); }
  static int ScaleAdd(realtype c, BlockMatrix& A, BlockMatrix& B) { return SUNMatScaleAdd(c, A, B); }

protected:
  bool isOwner_;
  SUNMatrix A_;

};

template<SUNMatrix_ID ID>
class BlockDenseMatrix : public BlockMatrix<ID>
{
public:
  BlockDenseMatrix(SUNContext sunctx = nullptr) : BlockMatrix<ID>(sunctx) { }

  BlockDenseMatrix(SUNMatrix A) : BlockMatrix<ID>(A) { }

  BlockDenseMatrix(sunindextype M, sunindextype N, SUNMemoryType memType,
                   SUNMemoryHelper mhelp, SUNExecQueue Q, SUNContext sunctx);

  BlockDenseMatrix(sunindextype nBlocks, sunindextype M, sunindextype N,
                   SUNMemoryType memType, SUNMemoryHelper mhelp, SUNExecQueue Q,
                   SUNContext sunctx);

  sunindextype BlockRows();
  sunindextype BlockCols();
  sunindextype Rows();
  sunindextype Cols();
  sunindextype NumBlocks();
  realtype* Data();
  realtype** BlockData();

  sunindextype LDim();
  sunindextype BlockLDim();

  realtype* Block(sunindextype k);
  realtype* Column(sunindextype j);
  realtype* BlockColumn(sunindextype k, sunindextype j);
};

#ifdef SUNDIALS_SUNMATRIX_MAGMADENSE
template<>
BlockDenseMatrix<SUNMATRIX_MAGMADENSE>::BlockDenseMatrix(sunindextype M, sunindextype N,
  SUNMemoryType memType, SUNMemoryHelper mhelp, SUNExecQueue Q, SUNContext sunctx)
{
  A_ = SUNMatrix_MagmaDense(M, N, memType, mhelp, Q, sunctx);
}

template<>
BlockDenseMatrix<SUNMATRIX_MAGMADENSE>::BlockDenseMatrix(sunindextype nBlocks, sunindextype M,
  sunindextype N, SUNMemoryType memType, SUNMemoryHelper mhelp, SUNExecQueue Q, SUNContext sunctx)
{
  A_ = SUNMatrix_MagmaDenseBlock(nBlocks, M, N, memType, mhelp, Q, sunctx);
}

template<>
sunindextype BlockDenseMatrix<SUNMATRIX_MAGMADENSE>::BlockRows()
{ return SUNMatrix_MagmaDense_BlockRows(A_); }

template<>
sunindextype BlockDenseMatrix<SUNMATRIX_MAGMADENSE>::BlockCols()
{ return SUNMatrix_MagmaDense_BlockColumns(A_); }

template<>
sunindextype BlockDenseMatrix<SUNMATRIX_MAGMADENSE>::Rows()
{ return SUNMatrix_MagmaDense_Rows(A_); }

template<>
sunindextype BlockDenseMatrix<SUNMATRIX_MAGMADENSE>::Cols()
{ return SUNMatrix_MagmaDense_Columns(A_); }

template<>
sunindextype BlockDenseMatrix<SUNMATRIX_MAGMADENSE>::NumBlocks()
{ return SUNMatrix_MagmaDense_NumBlocks(A_); }

template<>
sunindextype BlockDenseMatrix<SUNMATRIX_MAGMADENSE>::LDim()
{ return SUNMatrix_MagmaDense_LData(A_); }

template<>
sunindextype BlockDenseMatrix<SUNMATRIX_MAGMADENSE>::BlockLDim()
{ return SUNMatrix_MagmaDense_BlockLData(A_); }

template<>
realtype* BlockDenseMatrix<SUNMATRIX_MAGMADENSE>::Data()
{ return SUNMatrix_MagmaDense_Data(A_); }

template<>
realtype** BlockDenseMatrix<SUNMATRIX_MAGMADENSE>::BlockData()
{ return SUNMatrix_MagmaDense_BlockData(A_); }

template<>
realtype* BlockDenseMatrix<SUNMATRIX_MAGMADENSE>::Block(sunindextype k)
{ return SUNMatrix_MagmaDense_Block(A_, k); }

template<>
realtype* BlockDenseMatrix<SUNMATRIX_MAGMADENSE>::Column(sunindextype j)
{ return SUNMatrix_MagmaDense_Column(A_, j); }

template<>
realtype* BlockDenseMatrix<SUNMATRIX_MAGMADENSE>::BlockColumn(sunindextype k, sunindextype j)
{ return SUNMatrix_MagmaDense_BlockColumn(A_, k, j); }
#endif

#ifdef SUNDIALS_SUNMATRIX_ONEMKLDENSE
template<>
BlockDenseMatrix<SUNMATRIX_ONEMKLDENSE>::BlockDenseMatrix(sunindextype M, sunindextype N,
  SUNMemoryType memType, SUNMemoryHelper mhelp, SUNExecQueue Q, SUNContext sunctx)
{
  A_ = SUNMatrix_MagmaDense(M, N, memType, mhelp, Q, sunctx);
}

template<>
BlockDenseMatrix<SUNMATRIX_ONEMKLDENSE>::BlockDenseMatrix(sunindextype nBlocks, sunindextype M,
  sunindextype N, SUNMemoryType memType, SUNMemoryHelper mhelp, SUNExecQueue Q, SUNContext sunctx)
{
  A_ = SUNMatrix_OneMklDenseBlock(nBlocks, M, N, memType, mhelp, Q, sunctx);
}

template<>
sunindextype BlockDenseMatrix<SUNMATRIX_ONEMKLDENSE>::BlockRows()
{ return SUNMatrix_OneMklDense_BlockRows(A_); }

template<>
sunindextype BlockDenseMatrix<SUNMATRIX_ONEMKLDENSE>::BlockCols()
{ return SUNMatrix_OneMklDense_BlockColumns(A_); }

template<>
sunindextype BlockDenseMatrix<SUNMATRIX_ONEMKLDENSE>::Rows()
{ return SUNMatrix_OneMklDense_Rows(A_); }

template<>
sunindextype BlockDenseMatrix<SUNMATRIX_ONEMKLDENSE>::Cols()
{ return SUNMatrix_OneMklDense_Columns(A_); }

template<>
sunindextype BlockDenseMatrix<SUNMATRIX_ONEMKLDENSE>::NumBlocks()
{ return SUNMatrix_OneMklDense_NumBlocks(A_); }

template<>
sunindextype BlockDenseMatrix<SUNMATRIX_ONEMKLDENSE>::LDim()
{ return SUNMatrix_OneMklDense_LData(A_); }

template<>
sunindextype BlockDenseMatrix<SUNMATRIX_ONEMKLDENSE>::BlockLDim()
{ return SUNMatrix_OneMklDense_BlockLData(A_); }

template<>
realtype* BlockDenseMatrix<SUNMATRIX_ONEMKLDENSE>::Data()
{ return SUNMatrix_OneMklDense_Data(A_); }

template<>
realtype** BlockDenseMatrix<SUNMATRIX_ONEMKLDENSE>::BlockData()
{ return SUNMatrix_OneMklDense_BlockData(A_); }

template<>
realtype* BlockDenseMatrix<SUNMATRIX_ONEMKLDENSE>::Block(sunindextype k)
{ return SUNMatrix_OneMklDense_Block(A_, k); }

template<>
realtype* BlockDenseMatrix<SUNMATRIX_ONEMKLDENSE>::Column(sunindextype j)
{ return SUNMatrix_OneMklDense_Column(A_, j); }

template<>
realtype* BlockDenseMatrix<SUNMATRIX_ONEMKLDENSE>::BlockColumn(sunindextype k, sunindextype j)
{ return SUNMatrix_OneMklDense_BlockColumn(A_, k, j); }
#endif

template<SUNMatrix_ID ID>
class BlockSparseMatrix : public BlockMatrix<ID>
{
public:
  BlockSparseMatrix(SUNContext sunctx = nullptr) : BlockMatrix<ID>(sunctx) { }

  BlockSparseMatrix(SUNMatrix A) : BlockMatrix<ID>(A) { }

  BlockSparseMatrix(sunindextype M, sunindextype N, sunindextype NNZ, SUNMemoryType memType,
                    SUNMemoryHelper mhelp, SUNExecQueue Q, SUNContext sunctx);

  BlockSparseMatrix(sunindextype nBlocks, sunindextype M, sunindextype N, sunindextype NNZ,
                    SUNMemoryType memType, SUNMemoryHelper mhelp, SUNExecQueue Q,
                    SUNContext sunctx);

  sunindextype BlockRows();
  sunindextype BlockCols();
  sunindextype Rows();
  sunindextype Cols();
  sunindextype NumBlocks();
  realtype* Data();
  realtype** BlockData();

  realtype* Block(sunindextype k);
  realtype* Column(sunindextype j);
  realtype* BlockColumn(sunindextype k, sunindextype j);
};

#ifdef SUNDIALS_SUNMATRIX_ONEMKLDENSE
template<>
BlockSparseMatrix<SUNMATRIX_CUSPARSE>::BlockSparseMatrix(sunindextyCUSPARSEdextype N,
  SUNMemoryType memType, SUNMemoryHelper mhelp, SUNExecQueue Q, SUNContext sunctx)
{
  A_ = SUNMatrix_MagmaDense(M, N, memType, mhelp, Q, sunctx);
}

template<>
BlockSparseMatrix<SUNMATRIX_CUSPARSE>::BlockSparseMatrix(sunindextyCUSPARSE sunindextype M,
  sunindextype N, SUNMemoryType memType, SUNMemoryHelper mhelp, SUNExecQueue Q, SUNContext sunctx)
{
  A_ = SUNMatrix_OneMklDenseBlock(nBlocks, M, N, memType, mhelp, Q, sunctx);
}

template<>
sunindextype BlockSparseMatrix<SUNMATRIX_CUSPARSE>::BlockRows()
{ return SUNMatrix_OneMklDense_BlockRows(A_); }

template<>
sunindextype BlockSparseMatrix<SUNMATRIX_CUSPARSE>::BlockCols()
{ return SUNMatrix_OneMklDense_BlockColumns(A_); }

template<>
sunindextype BlockSparseMatrix<SUNMATRIX_CUSPARSE>::Rows()
{ return SUNMatrix_OneMklDense_Rows(A_); }

template<>
sunindextype BlockSparseMatrix<SUNMATRIX_CUSPARSE>::Cols()
{ return SUNMatrix_OneMklDense_Columns(A_); }

template<>
sunindextype BlockSparseMatrix<SUNMATRIX_CUSPARSE>::NumBlocks()
{ return SUNMatrix_OneMklDense_NumBlocks(A_); }

template<>
sunindextype BlockSparseMatrix<SUNMATRIX_CUSPARSE>::LDim()
{ return SUNMatrix_OneMklDense_LData(A_); }

template<>
sunindextype BlockSparseMatrix<SUNMATRIX_CUSPARSE>::BlockLDim()
{ return SUNMatrix_OneMklDense_BlockLData(A_); }

template<>
realtype* BlockSparseMatrix<SUNMATRIX_CUSPARSE>::Data()
{ return SUNMatrix_OneMklDense_Data(A_); }

template<>
realtype** BlockSparseMatrix<SUNMATRIX_CUSPARSE>::BlockData()
{ return SUNMatrix_OneMklDense_BlockData(A_); }

template<>
realtype* BlockSparseMatrix<SUNMATRIX_CUSPARSE>::Block(sunindextype k)
{ return SUNMatrix_OneMklDense_Block(A_, k); }

template<>
realtype* BlockSparseMatrix<SUNMATRIX_CUSPARSE>::Column(sunindextype j)
{ return SUNMatrix_OneMklDense_Column(A_, j); }

template<>
realtype* BlockSparseMatrix<SUNMATRIX_CUSPARSE>::BlockColumn(sunindextype k, sunindextype j)
{ return SUNMatrix_OneMklDense_BlockColumn(A_, k, j); }
#endif


}
