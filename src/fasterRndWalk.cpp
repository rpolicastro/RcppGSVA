#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double fasterRndWalk(IntegerVector gSetIdx, IntegerVector geneRanking, int j, NumericMatrix Ra) {
  int n = geneRanking.size();
  int k = gSetIdx.size();
  IntegerVector idxs = match(gSetIdx, geneRanking) - 1;
  
  double sum1 = 0;
  double sum2 = 0;
  for (int i = 0; i < k; ++i) {
    int idx = idxs[i];
    double value = Ra(geneRanking[idx] - 1, j - 1);
    sum1 += value * (n - idx);
    sum2 += value;
  }
  
  double stepCDFinGeneSet2 = sum1 / sum2;
  
  double stepCDFoutGeneSet2 = (n * (n + 1) / 2 - sum(n - idxs + 1)) / (n - k);
  
  double walkStat = stepCDFinGeneSet2 - stepCDFoutGeneSet2;
  
  return walkStat;
}
