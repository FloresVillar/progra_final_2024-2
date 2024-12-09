import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Scanner;
import java.util.Random;
//=================================================================
public class AccesoDisco {
    private static int BLOCK =3;
    private static String FILENAME="DATADISCO.TXT";
    private static int filas= 6;
    private static int columnas = 6;
    private static int  []anchos; 
    AccesoDisco(){
        anchos =new int[columnas];
        Random ran  =new Random();
        for(int i=0;i<anchos.length;i++){
            anchos[i] = ran.nextInt(5)+1;
        }
    }
    //-------------------------------------------------------
    public static void main(String[]args){
        CrearData();
        double[][]A=ReadData();
        WriteInFile("DATADISCOINFILE.TXT", A);
    }
    //---------------------------------------------------------------
    public static void CrearData(){
        FileWriter fw;
        Random ran = new Random();
        try{
            fw = new FileWriter(FILENAME);
            double max = Math.pow(10,BLOCK-1);
            double min = Math.pow(10,BLOCK-2);
            for(int i= 0;i<filas;i++){
                for(int j=0;j<columnas;j++){
                    double num = ran.nextDouble()*(max-min)+min;
                    fw.write((long)num +" ");
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
            int anchoDato = (int)(RAF.length()/(filas*columnas));
            System.out.println(anchoDato);
            byte[]record = new byte[anchoDato];
            for(int i= 0;i<filas;i++){
                for(int j=0;j<columnas;j++){
                    RAF.seek(i*anchoDato*columnas+anchoDato*j);
                    RAF.read(record);
                    String cad = leerRecord(record);
                    M[i][j] =Double.parseDouble(cad);
                    System.out.print(M[i][j]+"\t");
                }
                System.out.println();
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
    public static String leerRecord(byte[]record){
        String cadena = " ";
        for(int i=0;i<record.length;i++){
            cadena =cadena + (char)record[i];
        }
        return cadena.trim();
    }
    //--------------------------------------------------
    private static int PosicionarSeekColumna(int j){
        int p=0;
        for(int i=0;i<=j;i++){
            p+=anchos[i];
        }
        return  p;
    }
    //-------------------------------------------------------------------
    public static double ProductoInterno(int fila, int columna,double[][]M){
        RandomAccessFile Rfila ;
        RandomAccessFile Rcolumna;
        int anchoColumnas =PosicionarSeekColumna(columnas);
        try{
            Rfila = new RandomAccessFile(FILENAME, "r");
            Rcolumna = new RandomAccessFile(FILENAME, "r");
            //celdaFila * celdaColumna
            double suma =0.0;
            for(int k=0;k<M[0].length;k++){
                Rfila.seek(fila*anchoColumnas+PosicionarSeekColumna(k));
                byte [] record  = new byte[anchos[k]];
                Rcolumna.seek(k*anchoColumnas+PosicionarSeekColumna(columna));
                //suma+=(celdaFila*celdaColumna);
            }
        }catch(IOException e){
            e.printStackTrace();
        }
        return 1.0;
    }    
}
