import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.Random;
//=============================================================ssssss====
class AccesoDisco2 {
    private static String FILENAME="DATADISCO.TXT";
    private static int filas= 8;
    private static int columnas = 8;
    private static int  []anchos; 
    AccesoDisco2(){
    }
    //-------------------------------------------------------
    //--ejemplo de uso--------
    public static void main(String[]args){
        CrearData();
        double[][]A=ReadData();
        WriteInFile("DATADISCOINFILE.TXT", A);
        long inicio = System.nanoTime();
        double productoSerial = ProductoInterno(1,4,FILENAME,FILENAME);
        long tiempoSerial = System.nanoTime()-inicio;
        inicio = System.nanoTime();
        double productoParalelo = ProductoInternoPool(1, 4, FILENAME, FILENAME);
        long tiempoParalela = System.nanoTime() - inicio;
        System.out.println("serial: "+productoSerial+" tiempo="+tiempoSerial/1000000+" prodParalela: "+productoParalelo+"tiempo="+tiempoParalela/1000000);
        double[][]prod= ProductoSerial(FILENAME, FILENAME);
    }
    //---------------------------------------------------------------
    public static void CrearData(){
        anchos =new int[columnas];
        Random ran  =new Random();
        for(int i=0;i<anchos.length;i++){
            anchos[i] = ran.nextInt(4)+1;
            //System.out.println(anchos[i]+"\t");
        }
        System.out.println();
        FileWriter fw;
        try{
            fw = new FileWriter(FILENAME);
            double max;
            double min;
            for(int i= 0;i<filas;i++){
                for(int j=0;j<columnas;j++){
                    max = Math.pow(10,anchos[j]);
                    min = Math.pow(10,anchos[j]-1);
                    double num = ran.nextDouble()*(max-min)+min;
                    fw.write((long)num +"");
                }
            }
            fw.close();
        }catch(IOException e){
            e.printStackTrace();
        }
    }
    //---------------------------------------------------------------    
    public static double[][] ReadData(){
        RandomAccessFile RAF ;
        double[][] M = new double[filas][columnas];
        try{ 
            RAF = new RandomAccessFile(FILENAME,"rw");
            int anchoColumnas = PosicionarSeekColumna(columnas);
            for(int i= 0;i<filas;i++){
                for(int j=0;j<columnas;j++){
                    byte[]record = new byte[anchos[j]];
                    RAF.seek(i*anchoColumnas+PosicionarSeekColumna(j));
                    RAF.read(record);
                    String cad = ConvertirRecord(record);
                    M[i][j] =Double.parseDouble(cad); 
                } 
            }
        }catch(IOException e){
            e.printStackTrace();
        }
        return M;
    }
    //---------------------------------------------
    //---------------------------------------------
    public static void WriteInFile(String FILE,double[][]M){
        FileWriter fw;
        try{    
            fw = new FileWriter(FILE);
            for(int i=0;i<M.length;i++){
                for(int j=0;j<M[0].length;j++){
                    fw.write(M[i][j]+" ");
                }
                fw.write("\n");
            }
            fw.close();
        }catch(IOException e){
            e.printStackTrace();
        }
    }
    //-------------------------------------------------
    public static String ConvertirRecord(byte[]record){
        String cadena = "";
        for(int i=0;i<record.length;i++){
            cadena =cadena + (char)record[i];
        }
        return cadena.trim();
    }
    //--------------------------------------------------
    private static int PosicionarSeekColumna(int j){
        int p=0;
        for(int i=0;i<j;i++){
            p+=anchos[i];
        }
        return  p;
    }
    //-------------------------------------------------------------------
    public static double ProductoInterno(int fil, int colum,String FILE1,String FILE2){
        int fila = fil-1;
        int columna = colum -1;
        RandomAccessFile Rfila ;
        RandomAccessFile Rcolumna;
        double suma=0.0;
        int anchoColumnas =PosicionarSeekColumna(columnas);
        try{
            Rfila = new RandomAccessFile(FILE1, "r");
            Rcolumna = new RandomAccessFile(FILE2, "r");
            //celdaFila * celdaColumna
            for(int k=0;k<columnas;k++){
                double eFila =ObtenerCelda(Rfila,fila*anchoColumnas+PosicionarSeekColumna(k),anchos[k]);
                double eColumna=ObtenerCelda(Rcolumna,k*anchoColumnas+PosicionarSeekColumna(columna),anchos[columna]);
                suma+=(eFila+eColumna);
                /* esto lo hace obtener celda 
                Rfila.seek(fila*anchoColumnas+PosicionarSeekColumna(k));
                byte [] recordF  = new byte[anchos[k]];
                Rfila.read(recordF);
                Double eFila = Double.parseDouble(ConvertirRecord(recordF));
                Rcolumna.seek(k*anchoColumnas+PosicionarSeekColumna(columna));
                byte [] recordC = new byte[anchos[columna]];
                Rcolumna.read(recordC);
                Double eColumna = Double.parseDouble(ConvertirRecord(recordC));
                suma+=(eFila*eColumna);
                */
            }
            Rfila.close();
            Rcolumna.close();
        }catch(IOException e){
            e.printStackTrace();
        }
        return suma;
    } 
    //-----------------------------------------------------------------------------
    public static double ProductoInternoPool(int fil, int colum,String FILE1,String FILE2){
        int fila = fil;
        int columna = colum;
        double []suma=new double[]{0,0,0,0};
        int anchoColumnas =PosicionarSeekColumna(columnas);
        ExecutorService poolProductoInterno = Executors.newFixedThreadPool(4);
        poolProductoInterno.submit(()->{
            try{
                RandomAccessFile Rfila = new RandomAccessFile(FILE1, "r");
                RandomAccessFile Rcolumna = new RandomAccessFile(FILE2, "r");
                //celdaFila * celdaColumna
                for(int k=0;k<(int)(1*columnas/4);k++){
                    double eFila =ObtenerCelda(Rfila,fila*anchoColumnas+PosicionarSeekColumna(k),anchos[k]);
                    double eColumna=ObtenerCelda(Rcolumna,k*anchoColumnas+PosicionarSeekColumna(columna),anchos[columna]);
                    suma[0]+=(eFila+eColumna);
                }
                Rfila.close();
                Rcolumna.close();
            }catch(IOException e){
                e.printStackTrace();
            }
        });
        poolProductoInterno.submit(()->{
            try{
            RandomAccessFile Rfila = new RandomAccessFile(FILE1, "r");
            RandomAccessFile Rcolumna = new RandomAccessFile(FILE2, "r");
            //celdaFila * celdaColumna
            for(int k=(int)(columnas/4);k<(int)(2*columnas/4);k++){
                double eFila =ObtenerCelda(Rfila,fila*anchoColumnas+PosicionarSeekColumna(k),anchos[k]);
                double eColumna=ObtenerCelda(Rcolumna,k*anchoColumnas+PosicionarSeekColumna(columna),anchos[columna]);
                suma[1]+=(eFila+eColumna);
            }
            Rfila.close();
            Rcolumna.close();
        }catch(IOException e){
            e.printStackTrace();
        }});
        poolProductoInterno.submit(()->{
            try{
                RandomAccessFile Rfila = new RandomAccessFile(FILE1, "r");
                RandomAccessFile Rcolumna = new RandomAccessFile(FILE2, "r");
                //celdaFila * celdaColumna
                for(int k=(int)(2*columnas/4);k<(int)(3*columnas/4);k++){
                    double eFila =ObtenerCelda(Rfila,fila*anchoColumnas+PosicionarSeekColumna(k),anchos[k]);
                    double eColumna=ObtenerCelda(Rcolumna,k*anchoColumnas+PosicionarSeekColumna(columna),anchos[columna]);
                    suma[2]+=(eFila+eColumna);
                }
                Rfila.close();
                Rcolumna.close();
            }catch(IOException e){
                e.printStackTrace();
            }
        });
        poolProductoInterno.submit( ()->{
            try{
                RandomAccessFile Rfila = new RandomAccessFile(FILE1, "r");
                RandomAccessFile Rcolumna = new RandomAccessFile(FILE2, "r");
                //celdaFila * celdaColumna
                for(int k=(int)(3*columnas/4);k<columnas;k++){
                    double eFila =ObtenerCelda(Rfila,fila*anchoColumnas+PosicionarSeekColumna(k),anchos[k]);
                    double eColumna=ObtenerCelda(Rcolumna,k*anchoColumnas+PosicionarSeekColumna(columna),anchos[columna]);
                    suma[3]+=(eFila+eColumna);
                }
                Rfila.close();
                Rcolumna.close();
            }catch(IOException e){
                e.printStackTrace();
            }
        });
        poolProductoInterno.shutdown();
        try{
            poolProductoInterno.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        }catch(InterruptedException e){
            Thread.currentThread().interrupt();
        }
        return (suma[0]+suma[1]+suma[2]+suma[3]);
    } 
    //------------------------------------------------------------------------------------------------
    public static double ObtenerCelda(RandomAccessFile RAF,int P,int W) throws IOException {
        RAF.seek(P);
        byte [] recordF  = new byte[W];
        RAF.read(recordF);
        Double e = Double.parseDouble(ConvertirRecord(recordF));
        return e;
    }
    //-----------------------------------------------------------------------------------------------------
    public static double[][] ProductoSerial(String FILE1,String FILE2){
        double[][]PROD = new double[filas][columnas];
        FileWriter fw ;
        try{
            fw = new FileWriter("PRODUCTO.TXT");
            for(int f=0;f<filas;f++){
                for(int c=0;c<columnas;c++){
                    PROD[f][c] =  ProductoInternoPool(f, c, FILE1, FILE2);
                    fw.write((double)PROD[f][c]+" ");
                }
                fw.write("\n");
        }
        fw.close();
    }catch(IOException e){
        e.printStackTrace();
    }
        return PROD;
    } 
    //------------------------------------------------------------------------------------------------------
}
