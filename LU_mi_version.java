import java.io.*;

public class LU_mi_version{
	static int n = 5;
	static double [][] a;
	static double [][] l;
	static double [][] u;
	static double [] b;
	static double [] x;
	static double [] y;
	public static void main(String[] args){
		Generar();
		Doolitle();
		Imprimir(a);
		Imprimir(l);
		Imprimir(u);
		Imprimir(AxB(l,u));
		y=SustitucionProgresiva();
		x=SustitucionRegresiva(y);
		Imprimir(x);
	}

	public static void Generar(){
		Random rnd =new Random();
		a= new double[n][n];
		l=new double[n][n];
		u=new double[n][n];
		b=new double[n];
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				a[i][j] = rnd.nextDouble() * 100;
				if(i==j){
					l[i][j]=1;
				}else{
					l[i][j]=0;
				}
				u[i][j]=0;
			}
			b[i] = rnd.nextDouble() * 100;
		}
	}
	public static void Doolitle(){
		for(int i=0;i<n;i++){
			for(int j=0;j<i;j++){
				double suma = 0;
				for(int k=0;k<j;k++){
					suma+=(l[i][k]*u[k][j]);
				}
				l[i][j]=(a[i][j]-suma)/u[j][j];
			}
			for(int j=i;j<n;j++){
				double suma = 0;
				for(int k=0;k<i;k++){
					suma+=(l[i][k]*u[k][j]);
				}
				u[i][j]=(a[i][j]-suma)/l[i][i];
			}
		}
	}
	public static double[] SustitucionProgresiva(){
		y = new double[n];
		for(int i=0;i<n;i++){
			double suma = 0;
			for(int j=0;j<i;j++){
				suma+=(l[i][j]*y[j]);
			}
			y[i] = (b[i]-suma)/l[i][i];
		}
		return y;
	}
	public static double[] SustitucionRegresiva(double[] y){
		x= new double[n];
		for(int i=n-1;i>=0;i--){
			double suma = 0;
			for(int j=i+1;j<n;j++){
				suma+=(u[i][j]*x[j]);
			}
			x[i] =(y[i]-suma)/u[i][i];
		}
		return x;
	} 
	public static void Imprimir(double[][] m){
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				System.out.printf("%12.2f",m[i][j]);
			}
			System.out.println();
		}
		System.out.println();
	}
	public static void Imprimir(double[] x){
		for(int i=0;i<n;i++){
			System.out.printf("%12.2f",x[i]);
		}
		System.out.println();		
	}
	public static double[][] AxB(double[][]A,double[][]B) {
        double[][] C=new double[n][n];
        for(int i=0;i<n;i++) {
            for(int j=0;j<n;j++) {
                C[i][j] = 0;
                for(int k=0;k<n;k++) {
                    C[i][j] = C[i][j] + A[i][k]*B[k][j];
                }
            }
        }
        return C;
    }
}