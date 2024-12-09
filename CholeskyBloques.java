import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.RandomAccessFile;
import java.io.IOException;
import java.io.InterruptedIOException;
import java.io.PrintWriter;
import java.util.LinkedList;
import java.util.Random;
import java.util.Scanner;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.io.File;
import java.io.BufferedReader;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;
import java.util.concurrent.RecursiveTask;
import java.util.concurrent.ForkJoinPool;
//====================================================================================================================================
//-factorizacion cholesky por bloques, tal como lo indica el paper el tamaño del bloque NB
//-influye en la eficiencia de la paralelizacion ,se prefieren NB grandes 64,128,256
//--con el uso de pool de hilos se consigue poco costo de sincronizacion y buen balance de carga====================================================================================================================================
//-tal como menciona el paper, la arquitectura del ordenador tambien influye
public class CholeskyBloques {
    private static int BLOCK = 3; //tamaño de datos BLOCK -1=3, el ultimo 'digito' es " ",al generar la simetrica el tamaño se duplica
    private static String FILENAMEMATRIZ = "DATACholeskyBloques.TXT"; 
    private static int N_global = 640; // no necesariamente multiplo de NB , PAR o IMPAR , se crean k bloques si es impar con un bloque mayor 
    private static double [][] A_Global =new double[N_global][N_global];
    private static double [][]LGlobal = new double[N_global][N_global]; 
    private static int NB = 64;//necesariamente de gran tamaño para evitar granularidad
    private static int k_global=N_global/NB;
    private static double[][][][] ABloques = new double[k_global][k_global][NB][NB];
    private static double[][][][] LBloques = new double[k_global][k_global][NB][NB];
    //------------------------------------------------------------------------------------------
//-inicio main----------------------------------------------------------------------------------------- 
    public static void main(String[]args){  
        DataSetCholesky dataGlobal = new DataSetCholesky(FILENAMEMATRIZ, N_global, N_global, BLOCK-1, "FileWriter");
        A_Global = dataGlobal.ReadDataRAF();
        //Imprimir("AGlobal", A_Global);
        long inicio = System.nanoTime();
        double[][] LSerial =CholeskySerial(Copiar(A_Global));   
        long finSerial = System.nanoTime()-inicio; 
        inicio = System.nanoTime();
        ObtenerABloquesHilosPool();
        CholeskybloquesHilosPool();
        ObtenerLGlobalHilosPool();
        long finParalelo = System.nanoTime()-inicio;
        //Imprimir("SERIAL", LSerial);
        //Imprimir("PARALELO", LGlobal);
        System.out.println("tiempo serial: "+ finSerial/10000000+"\ttiempo paralelo: "+finParalelo/10000000+"\n");
}//--fin main
//---------------------------------CHOLESKY POR BLOQUES-------------------------------------------------------------------------------------------------------      
//---------------Cholesky serial, es la base para posibles paralelizaciones-------------------------------------------------------------
    public static void Choleskybloques(){
        for(int j=0;j<k_global;j++){
            LBloques[j][j] =CholeskySerial(ABloques[j][j]);
            //Imprimir("LBLOQUE[j][j]",LBloques[j][j]);
            //trsm L21=A21 L11^-t          
            //trsm L31 =A31 L11^-t       
            for(int i=j+1;i<k_global;i++){
                //System.out.printf("hacia abajo ↓ %d\t",i); 
                LBloques[i][j] = TRSM(ABloques[i][j],LBloques[j][j]); //
            }
            //syrk A22=A22 -  L21L21^t,    paralelizar
            //syrk A32=A32 - L31 L21^t
            //syrk A33=A33 -  L31L31^t
            for(int i=j+1;i<k_global;i++){
                //System.out.printf("en las diagonales y ↓ %d\t",i);
                for(int ii=i;ii<k_global;ii++){
                    ABloques[ii][i] = SYRK(ABloques[ii][i],LBloques[ii][j],LBloques[i][j]);
                }
            }
        }    
    }
//-------------Paralelizando los grupos Trsm y Syrk con un array de Thread para cada grupo de tareas, las tareas de actualizaciones usan un hilo por tarea------------------------------------------------------------------------------------------------------------------------------
    public static void CholeskybloquesHilos(){   //
        for(int j=0;j<k_global;j++){ 
            LBloques[j][j] =CholeskySerial(ABloques[j][j]);
            //Imprimir("LBLOQUE[j][j]",LBloques[j][j]);
            //trsm L21=A21 L11^-t          
            //trsm L31 =A31 L11^-t      
            final int jj= j;
            Thread[]Trsm = new Thread[k_global-(jj+1)];
            for(int hil =0;hil<Trsm.length;hil++){
                final int hill = hil;
                double [][]LL=Copiar(LBloques[jj][jj]); //los hilos acceden a la misma matriz entonces se crea una copia
                Trsm[hill] = new Thread(new Runnable() {             //                     for(int i=j+1;i<k_global;i++){
                    public void run(){                              //                          LBloques[i][j] = TRSM(ABloques[i][j],LBloques[j][j]);                     
                            LBloques[(jj+1)+hill][jj] = TRSM(ABloques[(jj+1)+hill][jj],LL);               //                 i←jj+1+hill  
                        }
                    });
            }
            for (Thread thread : Trsm) {
                thread.start();
            }
            for (Thread thread : Trsm) {
                try{
                    thread.join();
                }catch(InterruptedException e){
                    Thread.currentThread().interrupt();
                }
            }  
            //syrk A22=A22 -  L21L21^t,    
            //syrk A32=A32 - L31 L21^t
            //syrk A33=A33 -  L31L31^t
            Thread[]Syrk = new Thread[k_global-(jj+1)];
            for(int hil =0 ;hil<Syrk.length;hil++){
                final int hill=hil;
                Syrk[hill]=new Thread(new Runnable() {
                    public void run(){ 
                            for(int ii=jj+1+hill;ii<k_global;ii++){
                                ABloques[ii][jj+1+hill] = SYRK(ABloques[ii][jj+1+hill],LBloques[ii][jj],LBloques[jj+1+hill][jj]);
                            }
                        }
                    });
                }
                for (Thread thread : Syrk) {
                    thread.start();
                }
                for (Thread thread : Syrk) {
                    try{
                        thread.join();
                    }catch(InterruptedException e){
                        Thread.currentThread().interrupt();
                    }
                } 
        }    
    }
//------------------------primera tentativa de optimizacion de los bloques trsm y syrk , definiendo explicitamente las dependencias entre ellas----------------------------------------------------
//-------------------------sin los arrays de hilos para los grupos trsm y syrk --evitando la sobrecarga de gestion 
    public static void CholeskybloquesCompletableFuture(){
        for(int j=0;j<k_global;j++){ 
            LBloques[j][j] =CholeskySerial(ABloques[j][j]);
            //Imprimir("LBLOQUE[j][j]",LBloques[j][j]);
            //trsm L21=A21 L11^-t          
            //trsm L31 =A31 L11^-t      
            final int jj= j;
            CompletableFuture<Void> trsm =CompletableFuture.runAsync(()->{
                for(int i=jj+1;i<k_global;i++){ 
                    LBloques[i][jj] = TRSM(ABloques[i][jj],LBloques[jj][jj]); //
                } 
            });  
            //syrk A22=A22 -  L21L21^t,    
            //syrk A32=A32 - L31 L21^t
            //syrk A33=A33 -  L31L31^t
            CompletableFuture<Void> syrk = trsm.thenRunAsync(()->{
                for(int i=jj+1;i<k_global;i++){ 
                    for(int ii=i;ii<k_global;ii++){
                        ABloques[ii][i] = SYRK(ABloques[ii][i],LBloques[ii][jj],LBloques[i][jj]);
                    }
                }
            });
            try{
                CompletableFuture.allOf(trsm,syrk).get(); 
            }catch(InterruptedException | ExecutionException e){
                e.printStackTrace();
            }
        }    
    }
//--------------------------------dentro de los ambitos de las tareas de añaden los arrays de thread-------------------------------------------------------------------
    public static void CholeskybloquesCompletableFutureHilos(){   
        for(int j=0;j<k_global;j++){ 
            LBloques[j][j] =CholeskySerial(ABloques[j][j]);
            //Imprimir("LBLOQUE[j][j]",LBloques[j][j]);
            //trsm L21=A21 L11^-t         
            //trsm L31 =A31 L11^-t       
            final int jj= j;
            CompletableFuture<Void> trsm =CompletableFuture.runAsync(()->{
                Thread[]Trsm = new Thread[k_global-(jj+1)];
                for(int hil =0;hil<Trsm.length;hil++){
                    final int hill = hil;
                    double [][]LL=Copiar(LBloques[jj][jj]); //los hilos acceden a la misma matriz entonces se crea una copia
                    Trsm[hill] = new Thread(new Runnable() {             //                     for(int i=j+1;i<k_global;i++){
                        public void run(){                              //                          LBloques[i][j] = TRSM(ABloques[i][j],LBloques[j][j]); 
                            LBloques[(jj+1)+hill][jj] = TRSM(ABloques[(jj+1)+hill][jj],LL);         //}          i←jj+1+hill                      //paralelizar               
                        }
                    });
                }
                for (Thread thread : Trsm) {
                    thread.start();
                }
                for (Thread thread : Trsm) {
                    try{
                        thread.join();
                    }catch(InterruptedException e){
                        Thread.currentThread().interrupt();
                    }
                } 
            });  
            //syrk A22=A22 -  L21L21^t,    
            //syrk A32=A32 - L31 L21^t
            //syrk A33=A33 -  L31L31^t
            CompletableFuture<Void> syrk = trsm.thenRunAsync(()->{
                Thread[]Syrk = new Thread[k_global-(jj+1)];
                for(int hil =0 ;hil<Syrk.length;hil++){
                    final int hill=hil;
                    Syrk[hill]=new Thread(new Runnable() {
                        public void run(){ 
                            for(int ii=jj+1+hill;ii<k_global;ii++){
                                ABloques[ii][jj+1+hill] = SYRK(ABloques[ii][jj+1+hill],LBloques[ii][jj],LBloques[jj+1+hill][jj]);
                            }
                        }
                    });
                }
                for (Thread thread : Syrk) {
                    thread.start();
                }
                for (Thread thread : Syrk) {
                    try{
                        thread.join();
                    }catch(InterruptedException e){
                        Thread.currentThread().interrupt();
                    }
                }
            });
            try{
                CompletableFuture.allOf(trsm,syrk).get(); 
            }catch(InterruptedException | ExecutionException e){
                e.printStackTrace();
            }
        }    
    }
//-------------Reemplazando los array de Threads por Instream, que gestionaran la creacion y destruccion de la misma cantidad de  k hilos----------------------------------------------------------------------------------------
    public static void CholeskybloquesIntStream() {
        for (int j = 0; j < k_global; j++) {
            LBloques[j][j] = CholeskySerial(ABloques[j][j]);
            final int jj=j;
            // TRSM 
            IntStream.range(jj + 1, k_global).parallel().forEach(i -> {
                LBloques[i][jj] = TRSM(ABloques[i][jj], LBloques[jj][jj]);
            });
    
            // SYRK 
            IntStream.range(jj + 1, k_global).parallel().forEach(i -> {
                for (int ii = i; ii < k_global; ii++) {
                    ABloques[ii][i] = SYRK(ABloques[ii][i], LBloques[ii][jj], LBloques[i][jj]);
                }
            });
        }
    }
//---------------Ubicando ambos grupos de tareas(hilos) dentro de tareas sincronizadas------------------------------------------------------------------------------
    public static void CholeskybloquesConcurrentIntStream(){
        for(int j=0;j<k_global;j++){ 
            LBloques[j][j] =CholeskySerial(ABloques[j][j]);
            //Imprimir("LBLOQUE[j][j]",LBloques[j][j]);
            //trsm L21=A21 L11^-t          
            //trsm L31 =A31 L11^-t       
            final int jj= j;
            CompletableFuture<Void> trsm = CompletableFuture.runAsync(() -> {
                IntStream.range(jj + 1, k_global).parallel().forEach(i -> {
                    LBloques[i][jj] = TRSM(ABloques[i][jj], LBloques[jj][jj]);
                });
            });  
            //syrk A22=A22 -  L21L21^t,    paralelizar
            //syrk A32=A32 - L31 L21^t
            //syrk A33=A33 -  L31L31^t
            CompletableFuture<Void> syrk = trsm.thenRunAsync(() -> {
                IntStream.range(jj + 1, k_global).forEach(i -> {
                    IntStream.range(i, k_global).parallel().forEach(ii -> {
                        ABloques[ii][i] = SYRK(ABloques[ii][i], LBloques[ii][jj], LBloques[i][jj]);
                    });
                });
            });
            
            try{
                CompletableFuture.allOf(trsm,syrk).get(); 
            }catch(InterruptedException | ExecutionException e){
                e.printStackTrace();
            }
        }    
    }
//------------------Esta la opcion optima (EL METODO QUE USARA MAIN), uso de herramientas optimizadas y poca sobrecarga de operaciones------------------------------------------------------------------------------
    public static void CholeskybloquesHilosPool(){   //
        for(int j=0;j<k_global;j++){ 
            LBloques[j][j] =CholeskySerial(ABloques[j][j]);
            //Imprimir("LBLOQUE[j][j]",LBloques[j][j]);
            //trsm L21=A21 L11^-t          
            //trsm L31 =A31 L11^-t      
            final int jj= j;
            ExecutorService poolTrsm =Executors.newFixedThreadPool(4); //reemplaza al array de Thread para el grupo trsm
            for(int i=jj+1;i<k_global;i++){
                    double [][]LL=Copiar(LBloques[jj][jj]); //los hilos acceden a la misma matriz entonces se crea una copia
                    final int ii=i;
                    poolTrsm.submit(()->{             //                        
                        LBloques[ii][jj] = TRSM(ABloques[ii][jj],LL);           
                    });
            }
            poolTrsm.shutdown();
            try{
                poolTrsm.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }      
            catch(InterruptedException e){
                Thread.currentThread().interrupt();
            }
            ExecutorService poolSyrk=Executors.newFixedThreadPool(4); //pool de 4 hilos para las tareas syrk
            //syrk A22=A22 -  L21L21^t,    
            //syrk A32=A32 - L31 L21^t
            //syrk A33=A33 -  L31L31^t
            for(int col =j+1 ;col<k_global;col++){
                final int coll=col;
                poolSyrk.submit(() ->{
                    for(int ii=coll;ii<k_global;ii++){
                        ABloques[ii][coll] = SYRK(ABloques[ii][coll],LBloques[ii][jj],LBloques[coll][jj]);
                    }
                });
            }
            poolSyrk.shutdown();
            try{
                poolSyrk.awaitTermination(Long.MAX_VALUE,TimeUnit.NANOSECONDS);
            }catch(InterruptedException e){
                Thread.currentThread().interrupt();
            } 
        }    
    }
//-------------------------------Lj+1j = Aj+1j * Ljj^-t------------------------------------------------------------------
    public static double[][] TRSM(double[][]Ab,double[][]L){
        double[][]U=Transponer(L);
         
        double[][]INV =InversaTriangularSuperior(U);
    
        return ProdParalelo(Ab, INV); 
    }
//--------------------------------A j+1j+1 = Aj+1j+1 - Lj+1jLj+1j^-t------------------------------------------------
    public static double[][]SYRK(double[][]M,double[][]L,double[][]Lt){
        double[][]LT = Transponer(Lt);

        double[][] LLT = ProdParalelo(L, LT); 
        
        return Resta(M,LLT);
    }
//---------------------OBTENCION DE LOS SUBBLOQUES Aij (ABLOQUES)----------------------------------------------
//---------------------se obtienen los subbloques a partir del A_Global, se itera en las columnas----------------------------------------------
public static void ObtenerABloquesSerial(){
    for (int i = 0; i < k_global; i++) {
        for (int j = 0; j < k_global; j++) {
            ABloques[i][j] = new double[NB][NB];  // Creando las  submatrices de tamaño NB x NB
        }
    } 
    //dividir en k bloques de tamaño NB(modo SERIAL)
    int k = k_global;
    for(int columna=0;columna<k;columna++){
        System.out.println("columna: "+columna);
        for(int abajo=columna;abajo<k;abajo++){ //maneja los bloques(numero del bloque) hacia abajo
            System.out.println("abajo: "+abajo);
            for(int i=abajo*NB;i<(abajo+1)*NB;i++){//las filas en cada bloque-abajo
                System.out.println("i: "+i);                                         
                System.arraycopy(A_Global[i], columna*NB, ABloques[abajo][columna][i-abajo*NB],0, NB);
            }                    //la fila i hasta NB filas                    //  ↑avanza en las filas hasta NB filas para evitar usar otro                                                                                          
        }
    }    
    /*  ejemplo de como funciona el llenado de una columna de bloques
        double [][] A11=new double[NB][NB];   
        for(int i=(1-1)*NB;i<1*NB;i++){
            System.arraycopy(A[i], 0, A11[i], 0, NB);
        }
        double [][] A21=new double[NB][NB];
        for(int i=(2-1)*NB;i<2*NB;i++){
            System.arraycopy(A[i], 0, A21[i], 0, NB);
        }
        double [][] A31=new double[NB][NB];
        for(int i=(3-1)*NB;i<3*NB;i++){
            System.arraycopy(A[i], 0, A31[i], 0, NB);
    }*/  
}
//--------------------------paralelizacion de las k columnas , 1 hilo por columna---------------------------------------------------
public static void ObtenerABloquesParalela(){
    for (int i = 0; i < k_global; i++) {
        for (int j = 0; j < k_global; j++) {
            ABloques[i][j] = new double[NB][NB];  
        }
    }     
    Thread [] columnas =new Thread[k_global];
    for (int hil=0;hil<k_global;hil++) {       //un hilo por columna ,se obtienen k(k+1)/2 bloques ,solo se procesa esa cantidad pues A es simetrica
        final int  hilo = hil;
        columnas[hilo]=new Thread(new Runnable(){
            public void run(){
                for(int abajo=hilo;abajo<k_global;abajo++){
                    for(int i=abajo*NB;i<(abajo+1)*NB;i++){
                        System.arraycopy(A_Global[i], hilo*NB, ABloques[abajo][hilo][i-abajo*NB], 0, NB);
                    }
                }
            }
        });
    }
    for (Thread thread : columnas) {
        thread.start();
    }
    for (Thread thread : columnas) {
        try{
            thread.join();
        }catch(InterruptedException e){
            Thread.currentThread().interrupt();
        }
    } 
}
//----------esta es la opcion Optima,se usa en MAIN, reemplazando el array de Thread por un pool de hilos de 4 hilos exactamente-------------------------------------------------------
public static void ObtenerABloquesHilosPool(){
    for (int i = 0; i < k_global; i++) {
        for (int j = 0; j < k_global; j++) {
            ABloques[i][j] = new double[NB][NB];  // Asignando submatrices de tamaño NB x NB
        }
    } 
    ExecutorService poolBloques = Executors.newFixedThreadPool(4);
    for (int hil=0;hil<k_global;hil++) {       //un hilo por columna se obtienen k(k+1)/2 bloques ,pues A es simetrica
        final int  hilo = hil;
        poolBloques.submit(()->{
                for(int abajo=hilo;abajo<k_global;abajo++){
                    for(int i=abajo*NB;i<(abajo+1)*NB;i++){
                        System.arraycopy(A_Global[i], hilo*NB, ABloques[abajo][hilo][i-abajo*NB], 0, NB);
                    } 
                }
        });
    }
    poolBloques.shutdown();
    try{
            poolBloques.awaitTermination(Long.MAX_VALUE,TimeUnit.NANOSECONDS);
    }catch(InterruptedException e){
            Thread.currentThread().interrupt();
    }    
}

//-----------------OBTENER L (ensamblar los Lij) L final--------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------
public static void ObtenerLGlobalSerial(){
    for(int k=0;k<k_global;k++){        //modo serial para la obtencion de los bloques
        for(int filas=k;filas<k_global;filas++){
            for(int fil=0;fil<NB;fil++){                    //        ↓determina la fila  ↓determina la columna
                   System.arraycopy(LBloques[filas][k][fil], 0, LGlobal[filas*NB+fil],k*NB, NB); 
            }
        } 
    }
}
//--------------------Se ensambla los bloques Lij obtenidos con k(1 por columna) hilos creados --------------------------------------------------------
public static void ObtenerLGlobalParalela(){
    /*     for(int k=0;k<k_global;k++){        //modo serial para la obtencion de los bloques
           for(int filas=k;filas<k_global;filas++){
                for(int fil=0;fil<NB;fil++){                      ↓determina la fila  ↓determina la columna
                   System.arraycopy(LBloques[filas][k][fil], 0, LGlobal[filas*NB+fil],k*NB, NB); 
                }
           } 
        }
    */  //EL METODO ObtenerLGlobal sigue el mismo planteamiento , se trabaja por columnas k(no hay dependencias entre ellas), hacia abajo  y las filas de los bloques
        Thread [] columnas = new Thread[k_global];
        for(int hil=0;hil<columnas.length;hil++){
            final int  hilo =hil;
            columnas[hilo] = new Thread(new Runnable(){
                public void run(){
                    for(int filas=hilo;filas<k_global;filas++){
                        for(int fil=0;fil<NB;fil++){
                            System.arraycopy(LBloques[filas][hilo][fil], 0, LGlobal[NB*filas+fil], hilo*NB, NB);
                        }
                    }
                }
            });
        }
        for (Thread thread : columnas) {
            thread.start();
        }
        for (Thread thread : columnas) {
            try{
                thread.join();
            }catch(InterruptedException e){
                Thread.currentThread().interrupt();
            }
        }   
    }
//--------opcion Optima, este metodo se usa en MAIN,Se ensambla los bloques obtenidos pero resulta una mejora usando un pool de 4 hilos-------------------------------------------------------------------------------------------------------------------
public static void ObtenerLGlobalHilosPool(){
   //aqui sencillamente se usa herramientas ya optimizadas de java.util.concurrent, reemplazando el array
    //de columnas por un pool de 4 hilos , hace exactamente lo mismo , pero (se entiende) esta optmizado
    ExecutorService poolLGlobal = Executors.newFixedThreadPool(4);
    for(int hil=0;hil<k_global;hil++){
        final int  hilo =hil;
        poolLGlobal.submit(()->{
                for(int filas=hilo;filas<k_global;filas++){
                    for(int fil=0;fil<NB;fil++){
                        System.arraycopy(LBloques[filas][hilo][fil], 0, LGlobal[NB*filas+fil], hilo*NB, NB);
                    }
                }
        });
    }
    poolLGlobal.shutdown();
    try{
        poolLGlobal.awaitTermination(Long.MAX_VALUE,TimeUnit.NANOSECONDS);
    }catch(InterruptedException e){
        Thread.currentThread().interrupt();
    }  
} 
//----------------------------------------------------------------------------- ---------------------------------------------------------------
//-------------------este metodo solo se hizo para el testeo manual , no es usado por los metodos de cholesky por bloques-----------------------------------------------------------
    public static void CopiarBloque(double[][] destino, double[][] bloque, int filaInicio, int columnaInicio) {
        int rows = bloque.length;
        int cols = bloque[0].length;
        for (int i = 0; i < rows; i++) {
            System.arraycopy(bloque[i], 0, destino[filaInicio + i], columnaInicio, cols);
        }
    }
//------------------------------------------------------------------------------------------------------------------------------------------
//----------------CHOLESKY serial , usando para la otencion del LL^t de los bloques diagonales Ajj antes de las tareas trsm y syrk---------------------------------------------------------------------------------------------------------------------------
    public static double[][] CholeskySerial(double[][] M){
        int n = M.length;
        double [][]G =new double[M.length][M.length];
        double[][]MM=Copiar(M);
        //Imprimir(G);
        for(int j =0;j<n;j++){
            double suma = 0;
            for(int k =0;k<j;k++){
                suma+=(G[j][k]*G[j][k]);
            }
            //System.out.printf("\nM[j][j]-suma:%1.2f\n",M[j][j]-suma);
            G[j][j] = Math.pow(MM[j][j]-suma,0.5);
            for(int i =j+1;i<n;i++){
                double sumai = 0;
                for(int m =0;m<j;m++){
                    sumai+=(G[i][m]*G[j][m]);       
                }
                G[i][j] = (MM[i][j]-sumai)/G[j][j];
            }
        }
        return G;
    }   
//-----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------
public static double[][] Resta(double [][]A,double[][]B){
    int n = A.length;
    double [][] resta = new double[n][n]; 
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            resta[i][j] = A[i][j]-B[i][j];
        }
    }
    return resta;
}
//----------------------------------------------------------------------------------------------------------------------------------- 
//--------------------metodo producto paralelo usado por trsm y syrk --------------------------------------------------------------------------------------------------------------
 public static  double[][] ProdParalelo(double[][]M1,double[][]M2){
    double [][]PROD=new double[M1.length][M2[0].length];
    int N=M1.length;
    double[][] temp1;
    double[][] temp2;
    double[][] temp3;			  // este analisis maneja si N es par o impar
    double[][] temp4;			  // mitad = (int)n/2 + r={0,1}   y esto en los indices pero era menos legible
    int mitad = (int)N/2 + (N%2); //para n par n/2    para  n impar n/2 + 1
    if(N%2==0){					  // mitad =2  si n=4    mitad=3 si n=5
        temp1 = new double[mitad][mitad];
        temp2 = new double[mitad][mitad];
        temp3 = new double[mitad][mitad];
        temp4 = new double[mitad][mitad];
    }
    else{									//ejemplo n=5
        temp1 = new double[mitad][mitad];     // [3][3] superior izq
        temp2 = new double[mitad][mitad-1];	  // [3][2] superior der	
        temp3 = new double[mitad-1][mitad];	  // [2][3] inferior izq
        temp4 = new double[mitad-1][mitad-1]; // [2][2] inferior der
    }
    ExecutorService poolProducto = Executors.newFixedThreadPool(4);
    poolProducto.submit(()->{   
        for(int i=0;i<mitad;i++){
            for(int j=0;j<mitad;j++){
                double suma=0;
                for(int k=0;k<N;k++){
                    suma+=M1[i][k]*M2[k][j];
                }
                temp1[i][j]=suma;
                }
            }
        });
    poolProducto.submit(()->{
        for(int i=0;i<mitad;i++){
            for(int j=mitad;j<N;j++){
                double suma=0;
                for(int k=0;k<N;k++){
                    suma+=M1[i][k]*M2[k][j];
                }
                temp2[i][j-mitad]=suma;
            }
        }
    });
    poolProducto.submit(()->{
        for(int i=mitad;i<N;i++){
            for(int j=0;j<mitad;j++){
                double suma=0;
                for(int k=0;k<N;k++){
                    suma+=M1[i][k]*M2[k][j];
                }
                temp3[i-mitad][j]=suma;
            }
        }
    });
    poolProducto.submit(()->{
        for(int i=mitad;i<N;i++){
            for(int j=mitad;j<N;j++){
                double suma=0;
                for(int k=0;k<N;k++){
                    suma+=M1[i][k]*M2[k][j];
                }
                temp4[i-mitad][j-mitad]=suma;
            }
        }
    });
    poolProducto.shutdown();
    try{
        poolProducto.awaitTermination(Long.MAX_VALUE,TimeUnit.NANOSECONDS);
    }catch(InterruptedException e){
        Thread.currentThread().interrupt();
    }						
    //ahora reemzamblando a PROD, java proporciona este metodo corroborado
    for (int i = 0; i < mitad; i++) {										//ejemplo 	
        System.arraycopy(temp1[i], 0, PROD[i], 0, mitad); //N par n=4 mitad=2
        System.arraycopy(temp2[i], 0, PROD[i], mitad, N/2);		//N=5 impar mitad=3
    }																	//[3][3]    [3][2]
    for (int i = 0; i < N / 2; i++) {									//[2][3]	[2][2]	
        System.arraycopy(temp3[i], 0, PROD[i + mitad], 0, mitad);
        System.arraycopy(temp4[i], 0, PROD[i + mitad], mitad, N / 2);
    }
    return PROD;
}
//---------------------------------------------------------------------------------------------------------------------------------
//-------------usado para evitar modificar en las referencias en  trsm-------------------------------------------------------------
    public static double[][] Copiar(double[][]A){
		int n=A.length;
		double [][]COPIA=new double[n][n];
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				COPIA[i][j]=A[i][j];
			}
		}
		return COPIA;
	}
//--------------------------------------------------------------------------------------------------------------------------------------------
//------------------------usado por TRSM luego de transponer Ljj se obtiene una triangular superior, del cual se requiere su inversa------------------------------------------------------
    public static double[][] InversaTriangularSuperior(double[][]U){
        double [][] INV = new double[U.length][U.length];
        double factor = Math.pow(10,15);
        for(int j=U.length-1;j>=0;j--){
            if(U[j][j]==0){
                throw new ArithmeticException("no existe inversa");
            }
            INV[j][j] =1/U[j][j];
            for(int i=j-1;i>=0;i--){
                double suma = 0;
                for(int k=i+1;k<=j;k++){
                    suma +=(U[i][k]*INV[k][j]);    
                }
                INV[i][j] = -Math.round(suma*factor/U[i][i])/factor;
            }
        }
        return INV;
    }
//------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------
    public static double[][] InversaTriangularInferior(double[][]L){
        double factor = Math.pow(10,15);
        double[][]INV = new double[L.length][L.length];
        for(int j=0;j<L.length;j++){
            if(L[j][j]==0){
                throw new ArithmeticException("La matriz no es invertible");
            }
            INV[j][j]=Math.round(factor/L[j][j])/factor;
            //System.out.println(" "+INV[j][j]);
            for(int i=j+1;i<L.length;i++){
                double suma = 0;
                for(int k=j;k<=i-1;k++){
                    suma += (L[i][k] * INV[k][j]);  
                }
                INV[i][j] = -Math.round(suma*factor/L[i][i])/factor; 
            }
        }
        return INV;
    }
//--------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------     
    public static synchronized void Imprimir(String cad,double[][] M){
        System.out.println(cad);
        for(int i=0;i<M.length;i++){
            for(int j=0;j<M[0].length;j++){
                System.out.printf("%1.5f ",M[i][j]);
            }
            System.out.println();
        }
        System.out.println();
    } 
//----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------
    public static double [][] Transponer(double[][]M){
        double[][]T=new double  [M[0].length][M.length];
        for(int i=0;i<M.length;i++){
            for(int j=0;j<M[0].length;j++){
                T[j][i]=M[i][j];
            }
        }
        return T;
    }   
//----------------------------------------------------------------------------------------------------------------------------------------
}
//============================================================================================================================================
//============================================================================================================================================ 
//---Esta es la clase para crear los datos. Se definen dos formas de escribir los datos ------------------------------------------------- 
//---mediante PrintWriter que a su vez usa Scanner para la lectura
//---mediante FileWriter que usa RandomAccessFile para la lectura
//---estas dos maneras se reciben como argumento del constructor mediante el atributo "mode"
class DataSetCholesky{
    private static int filas;
    private static int columnas;
    private static String name;
    private static int BLOCK;
    private static String mode;
    DataSetCholesky(String nam,int fil,int col,int block,String mod){
        this.filas = fil;
        this.columnas = col;
        this.name =nam;
        this.BLOCK=block;
        this.mode =mod;
        if(mod.equals("PrintWriter")){
            CreateDataPrintWriter();
        }else{
            CreateDataFileWriter();
        }
        
    }
    //--------------------------------------------------------------------------------------------------------------------------------
    public static void CreateDataPrintWriter(){
        double num;
        Random ran = new Random();
        PrintWriter pw = null;
        try{
            pw = new PrintWriter(name);
            for(int i = 0;i<filas;i++){
                for(int j=0;j<columnas;j++){
                    double min = Math.pow(10,BLOCK-2);
                    double max= Math.pow(10,BLOCK-1);
                    System.out.printf("max= %f , min =%f\n",max,min);
                    num = min + ran.nextDouble()*(max - min);
                    pw.printf("%1.1f",num);
                }
                pw.println("");
            }
            pw.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
        }
    }
    //----------------------------------------------------------------------------------------------------------------------------------
    public static void CreateDataFileWriter(){
        double num;
        Random ran = new Random();
        FileWriter pw = null;
        try{
            pw = new FileWriter(name);
            for(int i = 0;i<filas;i++){
                for(int j=0;j<columnas;j++){
                    double min = Math.pow(10,BLOCK-2);
                    double max= Math.pow(10,BLOCK-1);
                    //System.out.printf("max= %f , min =%f\n",max,min);
                    num = min + ran.nextDouble()*(max - min);
                    pw.write((long)num+" ");
                }
            }
            pw.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
        }
    }
    //-----------------------------------------------------------------------------------------------------------------------------------
    public static double[][] ReadDataScanner(){
        if(mode.equals("PrintWriter")){
            double [][]M=new double[filas][columnas];
        Scanner scanner = null;
        double num;
        Path ruta = Paths.get(name);
        try{
            scanner = new Scanner(Files.newInputStream(ruta));
            scanner.useDelimiter("[;\\s]+");
            int i =0;
            int j =0;
            while(scanner.hasNext()){
                if(scanner.hasNextDouble()){
                    M[i][j] = scanner.nextDouble();
                    System.out.printf("%1.1f\t",M[i][j]);
                    j++;
                    if(j == columnas){
                        j = 0;
                        i++;
                        System.out.println();
                    }
                }else{
                    scanner.next();
                }
            }
            scanner.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
        }
        return M;
        }else return null;
    }
    //--------------------------------------------------------------------------------------------------------------------------------------
    public static void WriteData(String FILE){
        if(mode.equals("PrintWriter")){
            PrintWriter pw = null;
        double[][] A= ReadDataScanner();
        try{    
            pw = new PrintWriter(FILE);
            for(int i =0;i<A.length;i++){
                for(int j=0;j<A[0].length;j++){
                    //System.out.printf("write%1.1f\t",A[i][j]);
                    pw.printf("%1.1f;",A[i][j]);
                }
                System.out.println();
                pw.println("");
            }
            pw.close();
        }catch(FileNotFoundException e){
            System.out.println(e.getMessage());
        }
        }
    }
    //-----------------------------------------------------------------------------------------------------------------------------------
    public static void WriteData(String FILE,double[][]M){
        if(mode.equals("PrintWriter")){
            PrintWriter pw = null;
        try{
            pw = new PrintWriter(FILE);
            for(int i=0;i<M.length;i++){
                for(int j=0;j<M[0].length;j++){
                    //System.out.printf("write%1.1f\t",M[i][j]);
                    pw.printf("%1.1f;",M[i][j]);
                }
                //System.out.println();
                pw.println("");
            }
            pw.close();
        }catch(FileNotFoundException e){
            System.out.println(e.getMessage());
        }
        finally{
            pw.close();
        }
        }
    }
    //----------------------------------------------------------------------------------------------------------------------------------------
    public static String BufferToString(byte[]BUFFER){
        String cad = "";
        for(int i=0;i<BUFFER.length;i++){
            cad = cad + (char)BUFFER[i];
        }
        return cad;
    }
    //--------------------------------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------------------------------
    public static double[][] ReadDataRAF(String FILE){
        if(mode.equals("FileWriter")){
            double [][]M=new double[filas][columnas];
        RandomAccessFile RAF =null;
        try{
            RAF = new RandomAccessFile(FILE,"r");
            int W = (int)RAF.length()/(filas*columnas);
           // System.out.println("dentro de readRAF"+W);
            byte []BUFFER = new byte[W];
            for(int i =0;i<filas;i++){
                for(int j=0;j<columnas;j++){
                    RAF.seek(i*columnas*W + W*j);       //buena forma de direccionar
                    RAF.read(BUFFER);
                    String cad = BufferToString(BUFFER).trim();
                    M[i][j] = Double.parseDouble(cad);
                    System.out.printf("%1.1f\t", M[i][j]);
                    if(j==columnas - 1){
                        System.out.println();
                    }
                }
            }
            RAF.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
        }
        return M;
        }
        else{return null;}
    }
    //------------------------------------------------------------------------------------------------------------------------------------------
    public static double[][] ReadDataRAF(){
        if(mode.equals("FileWriter")){
            double [][]M=new double[filas][columnas];
        RandomAccessFile RAF =null;
        try{
            RAF = new RandomAccessFile(name,"r");
            int W = (int)RAF.length()/(filas*columnas);
           // System.out.println("dentro de readRAF"+W);
            byte []BUFFER = new byte[W];
            for(int i =0;i<filas;i++){
                for(int j=0;j<columnas;j++){
                    RAF.seek(i*columnas*W + W*j);       //buena forma de direccionar
                    RAF.read(BUFFER);
                    String cad = BufferToString(BUFFER).trim();
                    M[i][j] = Double.parseDouble(cad);
                    //System.out.printf("%1.1f\t", M[i][j]);
                    if(j==columnas - 1){
                        //System.out.println();
                    }
                }
            }
            RAF.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
        }
        return M = ProdParalelo(M, Transponer(M));
        }
        else{return null;}
    }
    //--------------------------------------------------------------------------------------------------------------------------------------------
    public static void WriteDataRAF(String FILE){
        if(mode.equals("FileWriter")){
            double [][]M=new double[filas][columnas];
        RandomAccessFile RAF =null;
        FileWriter fw = null;
        try{
            RAF = new RandomAccessFile(name,"r");
            fw = new FileWriter(FILE);
            int W = (int)RAF.length()/(filas*columnas);
           // System.out.println("dentro de readRAF"+W);
            byte []BUFFER = new byte[W];
            for(int i =0;i<filas;i++){
                for(int j=0;j<columnas;j++){
                    RAF.seek(i*columnas*W + W*j);
                    RAF.read(BUFFER);
                    String cad = BufferToString(BUFFER).trim();
                    M[i][j] = Double.parseDouble(cad);
                    fw.write((long)M[i][j]+" ");
                }
            }
            RAF.close();
            fw.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
            }
        }
    }
    //---------------------------------------------------------------------------------------------------------------------------------------
    public static double[][] ProdParalelo(double[][]M1,double[][]M2){
        double [][]PROD=new double[M1.length][M2[0].length];
        int N=M1.length;
        double[][] temp1;
        double[][] temp2;
        double[][] temp3;			  // este analisis maneja si N es par o impar
        double[][] temp4;			  // mitad = (int)n/2 + r={0,1}   y esto en los indices pero era menos legible
        int mitad = (int)N/2 + (N%2); //para n par n/2    para  n impar n/2 + 1
        if(N%2==0){					  // mitad =2  si n=4    mitad=3 si n=5
            temp1 = new double[mitad][mitad];
            temp2 = new double[mitad][mitad];
            temp3 = new double[mitad][mitad];
            temp4 = new double[mitad][mitad];
        }
        else{									//ejemplo n=5
            temp1 = new double[mitad][mitad];     // [3][3] superior izq
            temp2 = new double[mitad][mitad-1];	  // [3][2] superior der	
            temp3 = new double[mitad-1][mitad];	  // [2][3] inferior izq
            temp4 = new double[mitad-1][mitad-1]; // [2][2] inferior der
        }
        ExecutorService poolProducto = Executors.newFixedThreadPool(4);
        poolProducto.submit(()->{
            for(int i=0;i<mitad;i++){
                for(int j=0;j<mitad;j++){
                    double suma=0;
                    for(int k=0;k<N;k++){
                        suma+=M1[i][k]*M2[k][j];
                    }
                    temp1[i][j]=suma;
                    }
                }
        });
        poolProducto.submit(()->{
            for(int i=0;i<mitad;i++){
                for(int j=mitad;j<N;j++){
                    double suma=0;
                    for(int k=0;k<N;k++){
                        suma+=M1[i][k]*M2[k][j];
                    }
                    temp2[i][j-mitad]=suma;
                }
            }
        });
        poolProducto.submit(()->{
            for(int i=mitad;i<N;i++){
                for(int j=0;j<mitad;j++){
                    double suma=0;
                    for(int k=0;k<N;k++){
                        suma+=M1[i][k]*M2[k][j];
                    }
                    temp3[i-mitad][j]=suma;
                }
            }
        });
        poolProducto.submit(()->{
            for(int i=mitad;i<N;i++){
                for(int j=mitad;j<N;j++){
                    double suma=0;
                    for(int k=0;k<N;k++){
                        suma+=M1[i][k]*M2[k][j];
                    }
                    temp4[i-mitad][j-mitad]=suma;
                }
            }
        });
        poolProducto.shutdown();
        try{
                poolProducto.awaitTermination(Long.MAX_VALUE,TimeUnit.NANOSECONDS);
        }catch(InterruptedException e){
                e.printStackTrace();
        }						
        		//ahora reemzamblando a PROD, java proporciona este metodo corroborado
        for (int i = 0; i < mitad; i++) {										//ejemplo 	
            System.arraycopy(temp1[i], 0, PROD[i], 0, mitad); //N par n=4 mitad=2
            System.arraycopy(temp2[i], 0, PROD[i], mitad, N/2);		//N=5 impar mitad=3
        }																	//[3][3]    [3][2]
        for (int i = 0; i < N / 2; i++) {									//[2][3]	[2][2]	
            System.arraycopy(temp3[i], 0, PROD[i + mitad], 0, mitad);
            System.arraycopy(temp4[i], 0, PROD[i + mitad], mitad, N / 2);
        }
        return PROD;
    }
    //----------------------------------------------------------------------------------------------------------------------------
    public static double [][] Transponer(double[][]M){
        double[][]T=new double  [M[0].length][M.length];
        for(int i=0;i<M.length;i++){
            for(int j=0;j<M[0].length;j++){
                T[j][i]=M[i][j];
            }
        }
        return T;
    }
    //------------------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------------
}
