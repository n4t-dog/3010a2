import java.util.Arrays;

public class Matrix{
  public static void main(String[] args){
    double[][] coeff1 = {{1.0,2.0,3.0},{2.0,3.0,1.0},{3.0,1.0,2.0}};
    double[][] coeff2 = {{1.0,2.0,3.0},{2.0,3.0,1.0},{3.0,1.0,2.0}};

    double[] cons1 = {6.0,6.0,6.0};
    double[] cons2 = {6.0,6.0,6.0};

    double[] naiveSol = naiveGaussian(coeff1,cons1);
    double[] SPPSol = SPPGaussian(coeff2,cons2);

    printArray("naive solution", naiveSol);
    printArray("spp solution", SPPSol);
  }

  public static void fwdElimination(double[][] coeff, double[] cons){
    int n = coeff.length;
    for(int k=0;k<n-1;k++){
      for(int i=k+1;i<n;i++){
        double mult = coeff[i][k] / coeff[k][k];
        for(int j=k;j<n;j++)
          coeff[i][j] -= mult * coeff[k][j];
        cons[i] -= mult * cons[k];
      }
    }
  }

  public static void backSubst(double[][] coeff, double[] cons, double[] sol){
    int n = coeff.length;
    sol[n-1] = cons[n-1] / coeff[n-1][n-1];
    for(int i=n-2;i>=0;i--){
      double sum = cons[i];
      for(int j=i+1;j<n;j++){
        sum -= coeff[i][j] * sol[j];
      }
      sol[i] = sum / coeff[i][i];
    }
  }

  public static double[] naiveGaussian(double[][] coeff, double[] cons){
    double[] sol = new double[coeff.length];
    fwdElimination(coeff,cons);
    backSubst(coeff,cons,sol);
    return sol;
  }

  public static void SPPFwdElimination(double[][] coeff, double[] cons, int[] ind){
    int n = coeff.length;
    double[] scaling = new double[n];
    for(int i=0;i<n;i++){
      double smax = 0;
      for(int j=1;j<n;j++){
        smax = Math.max(smax,Math.abs(coeff[i][j]));
      }
      scaling[i] = smax;
    }
    for(int k=0;k<n-1;k++){
      double rmax = 0;
      int maxInd = k;
      for(int i=k;i<n;i++){
        double r = Math.abs(coeff[ind[i]][k] / scaling[ind[i]]);
        if(r>rmax){
          rmax = r;
          maxInd = i;
        }
      }
      int temp = ind[maxInd];
      ind[maxInd] = ind[k];
      ind[k] = temp;
      for(int i=k+1;i<n;i++){
        double mult = coeff[ind[i]][k] / coeff[ind[k]][k];
        for(int j=k+1;j<n;j++){
          coeff[ind[i]][j] -= mult * coeff[ind[k]][j];
        }
        cons[ind[i]] -= mult * cons[ind[k]];
      }
    }
  }

  public static void SPPBackSubst(double[][] coeff, double[] cons, double[] sol, int[] ind){
    int n = coeff.length;
    sol[n-1] = cons[ind[n-1]] / coeff[ind[n-1]][n-1];
    for(int i=n-2;i>=0;i--){
      double sum = cons[ind[i]];
      for(int j=i+1;j<n;j++){
        sum -= coeff[ind[i]][j] * sol[j];
      }
      sol[i] = sum / coeff[ind[i]][i];
    }
  }

  public static double[] SPPGaussian(double[][] coeff, double[] cons){
    int n = coeff.length;
    double[] sol = new double[n];
    int[] ind = new int[n];
    for(int i=0;i<n;i++){
      ind[i] = i;
    }
    SPPFwdElimination(coeff,cons,ind);
    SPPBackSubst(coeff,cons,sol,ind);
    return sol;
  }


  public static void printArray(String tag, double[] input){
    String output = tag + ": {";
    for(int i=0;i<input.length;i++)
      output += input[i] + ((i<input.length-1) ? "," : "");
    output += "}";
    System.out.println(output);
  }

  public static void printMatrix(String tag, double[][] input){
    for(int i=0;i<input.length;i++){
      printArray(tag+i,input[i]);
    }
  }
}
