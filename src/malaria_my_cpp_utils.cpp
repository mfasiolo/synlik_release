#include "synlik.h"

using namespace Rcpp;

/*
 * Prints the entries of a NumericVector
 */ 
inline void print(const NumericMatrix & Mat)
{
  for(int iRow = 0; iRow < Mat.nrow(); iRow++)
  {
    for(int iCol = 0; iCol < Mat.ncol(); iCol++)
     std::cout << Mat(iRow, iCol) << " ";
     
     std::cout << std::endl;
  }
}



/* 
 * Naive factorial calculation
 */
inline unsigned int myFactorial(int n)
{
    unsigned int ret = 1;
    for(unsigned int i = 1; i <= n; ++i)
      ret *= i;
    return ret;
}



/*
 * Takes the rowIndex-th row of Mat and stores it in vector Vett
 */ 
inline void RowToVector(const NumericMatrix & Mat, 
                              NumericVector & Vett, 
                        const int rowIndex, 
                              int colStart = -1,
                              int colEnd   = -1)
{
  
if(colStart == -1 && colEnd == -1)
{
  colStart = 0;
  colEnd = Mat.ncol();
} else {
        if( (colEnd - colStart) != Vett.size() ) stop("RowToVector: (colEnd - colStart) != Vett.size()");
        }

int iVett = 0;
for(int iCol = colStart; iCol != colEnd; iCol++, iVett++)
{
Vett[iVett] = Mat(rowIndex, iCol);
}
}

  
/*
 *  Equivalent to R seq(a, b, by)
 */ 
inline NumericVector seqCpp(const double & from, 
                            const double & to, 
                            const double & by )
{
  if(from > to && by > 0.0) stop("seqCpp: if from > to then you need by < 0");
  
  NumericVector outSeq( (to + 1e-10*(to - from) - from) / by  + 1 );
  NumericVector::iterator theEnd = outSeq.end();
  double tmp = from;
  
  for(NumericVector::iterator iter = outSeq.begin(); iter != theEnd; iter++, tmp += by)
  {
  *iter = tmp;
  }
  
  return outSeq;
}
  
  
