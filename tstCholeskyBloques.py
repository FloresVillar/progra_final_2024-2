import numpy as np

if __name__ == "__main__":
    L11 =np.array([[22802.67337, 0.00000, 0.00000 ,0.00000],
[17871.10999, 10086.02968, 0.00000 ,0.00000],
[13851.43965 ,8399.94470, 6795.69980 ,0.00000],
[16205.15437, 2135.01879, 6051.63983, 13238.08015]])
    A10=np.array([[368509378.00000 ,334329452.00000 ,262652770.00000 ,380720332.00000],
[397114485.00000, 350432462.00000 ,316364739.00000 ,375965118.00000], 
[466365863.00000, 426471740.00000 ,345866718.00000 ,432874355.00000],
[371443738.00000 ,371697971.00000 ,337512163.00000 ,386924395.00000]])
    L11t=np.transpose(L11)
    L11ti=np.linalg.inv(L11t)
    print(L11)
    print(L11t)
    print(L11ti)
    print(A10@L11ti)
    A = np.array([[1,2,4],[3,2,4],[6,7,8]])
    U = np.array([[3,5,4],[0,4,5],[0,0,9]])
    print(A)
    print(U)
    print(A@U)
    M =np.array([[-6395205761.41770800 ,-8066020962.38804100, -5984246602.65509000, -5401668591.06543600],
[-8066020962.38804100, -10055691338.84700800 ,-7453983359.02241200, -6761615613.37243500]
,[-5984246602.65509000 ,-7453983359.02241200 ,-5568834084.26806200, -5018381694.00170500]
,[-5401668591.06543600, -6761615613.37243500 ,-5018381694.00170500 ,-4544204069.29161600]
])
    print(M)
    l,v =np.linalg.eig(M)
    print(l)
    print(v)
"""//Imprimir("ABloques[1][0]",ABloques[1][0]);
        //Imprimir("ABloques[2][0]",ABloques[2][0]);
        //Imprimir("Transpuesta de L11",Transponer(L11));
       // Imprimir("L11 L11^t= ABloques[0][0] prodParalelo", ProdParalelo(L11, Transponer(L11)));
       // Imprimir("L11 L11^t= ABloques[0][0] prodTriangular", ProductoTriangular(L11, Transponer(L11)));
        //Imprimir("ABloques[0][0]", ABloques[0][0]);
        
        //Imprimir("L11^t L11^-t",ProdParalelo(Transponer(L11),InversaTriangularSuperior(Transponer(L11))));
        //Imprimir("L11^t L11^-t",ProductoTriangular(Transponer(L11),InversaTriangularSuperior(Transponer(L11))));
       // Imprimir("L11^t",Transponer(L11));
        //Imprimir("L11^-t",InversaTriangularSuperior(Transponer(L11)));
        Imprimir("ABloques[2][1]", ABloques[1][0]);
        //Imprimir("ABloques[1][0] L11^-t productoTriangular",ProductoTriangular(ABloques[1][0], InversaTriangularSuperior(Transponer(L11))));
        //Imprimir("ABloques[1][0] L11^-t productoParalelo",ProdParalelo(ABloques[1][0], InversaTriangularSuperior(Transponer(L11))));
        //Imprimir("ABloques[1][0] L11^-t productoSerial",ProdSerial(ABloques[1][0], InversaTriangularSuperior(Transponer(L11))));
        
        LGlobal manual
224.20303 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
173.84689 110.51361 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
165.86306 59.26169 88.26947 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
188.49433 96.63830 -41.37551 42.33197 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
152.40650 58.83079 9.13090 -17.84168 79.82166 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 
146.66171 60.70128 3.60238 -20.59716 -28.96060 67.13986 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
166.39828 70.87069 -31.08940 1.05606 -19.89239 -2.03376 77.41090 0.00000 0.00000 0.00000 0.00000 0.00000
156.75078 80.30110 -34.92648 36.21625 11.28833 30.89030 24.30188 76.57175 0.00000 0.00000 0.00000 0.00000
167.09854 60.98288 33.07573 -37.29767 29.56542 4.90258 -27.27030 22.57641 52.15873 0.00000 0.00000 0.00000
154.24858 85.16022 12.96168 -27.69598 26.22064 42.72453 -59.50533 31.91825 25.19706 40.60272 0.00000 0.00000
198.38269 87.58909 39.43565 11.22552 -29.40991 1.64595 31.38385 12.76084 -15.36405 -27.98947 43.15526 0.00000
193.79756 50.02911 1.09138 -38.51254 -5.70463 58.38705 4.28480 -6.08239 13.51346 13.11122 -28.29168 15.07672

factorizacion cholesky Serial
224.20303 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
173.84689 110.51361 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 
165.86306 59.26169 88.26947 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
188.49433 96.63830 -41.37551 42.33197 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
152.40650 58.83079 9.13090 -17.84168 79.82166 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
146.66171 60.70128 3.60238 -20.59716 -28.96060 67.13986 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
166.39828 70.87069 -31.08940 1.05606 -19.89239 -2.03376 77.41090 0.00000 0.00000 0.00000 0.00000 0.00000
156.75078 80.30110 -34.92648 36.21625 11.28833 30.89030 24.30188 76.57175 0.00000 0.00000 0.00000 0.00000
167.09854 60.98288 33.07573 -37.29767 29.56542 4.90258 -27.27030 22.57641 52.15873 0.00000 0.00000 0.00000
154.24858 85.16022 12.96168 -27.69598 26.22064 42.72453 -59.50533 31.91825 25.19706 40.60272 0.00000 0.00000
198.38269 87.58909 39.43565 11.22552 -29.40991 1.64595 31.38385 12.76084 -15.36405 -27.98947 43.15526 0.00000
193.79756 50.02911 1.09138 -38.51254 -5.70463 58.38705 4.28480 -6.08239 13.51346 13.11122 -28.29168 15.07672


















LGlobal manual
236.79316 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
164.65847 123.54994 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 
205.04393 48.98246 87.82770 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
228.40187 31.06191 -30.26272 97.01501 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
177.33620 27.79517 44.80149 40.00952 74.49405 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
155.29165 20.84917 -42.39954 -3.28448 75.06149 60.20864 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
194.49042 121.26679 -2.37994 57.74595 40.86788 13.28884 34.98295 0.00000 0.00000 0.00000 0.00000 0.00000
152.43684 113.74335 -19.46626 23.10936 2.14027 18.17459 -22.22747 45.12791 0.00000 0.00000 0.00000 0.00000
239.18765 47.81650 29.61313 21.92439 -2.67935 23.05199 -35.15035 15.92032 48.35964 0.00000 0.00000 0.00000
151.50776 4.74273 20.52819 -4.93464 -27.77264 0.84809 -10.38770 15.42197 0.18824 60.65751 0.00000 0.00000
162.20486 46.69040 -35.14986 -34.36495 -19.70452 32.24332 -15.97056 13.08230 7.82976 35.84936 27.93818 0.00000 
180.23747 63.67768 45.35370 15.30692 -2.60285 -13.84147 -0.34623 7.59684 2.29330 -50.49023 -12.67626 10.53301

factorizacion cholesky Serial
236.79316 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
164.65847 123.54994 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
205.04393 48.98246 87.82770 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
228.40187 31.06191 -30.26272 97.01501 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 
177.33620 27.79517 44.80149 40.00952 74.49405 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
155.29165 20.84917 -42.39954 -3.28448 75.06149 60.20864 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
194.49042 121.26679 -2.37994 57.74595 40.86788 13.28884 34.98295 0.00000 0.00000 0.00000 0.00000 0.00000
152.43684 113.74335 -19.46626 23.10936 2.14027 18.17459 -22.22747 45.12791 0.00000 0.00000 0.00000 0.00000
239.18765 47.81650 29.61313 21.92439 -2.67935 23.05199 -35.15035 15.92032 48.35964 0.00000 0.00000 0.00000
151.50776 4.74273 20.52819 -4.93464 -27.77264 0.84809 -10.38770 15.42197 0.18824 60.65751 0.00000 0.00000
162.20486 46.69040 -35.14986 -34.36495 -19.70452 32.24332 -15.97056 13.08230 7.82976 35.84936 27.93818 0.00000
180.23747 63.67768 45.35370 15.30692 -2.60285 -13.84147 -0.34623 7.59684 2.29330 -50.49023 -12.67626 10.53301 


LGlobal iteracion
236.79316 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 
164.65847 123.54994 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
205.04393 48.98246 87.82770 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
228.40187 31.06191 -30.26272 97.01501 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
177.33620 27.79517 44.80149 40.00952 74.49405 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
155.29165 20.84917 -42.39954 -3.28448 75.06149 60.20864 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 
194.49042 121.26679 -2.37994 57.74595 40.86788 13.28884 34.98295 0.00000 0.00000 0.00000 0.00000 0.00000
152.43684 113.74335 -19.46626 23.10936 2.14027 18.17459 -22.22747 45.12791 0.00000 0.00000 0.00000 0.00000
239.18765 47.81650 29.61313 21.92439 -2.67935 23.05199 -35.15035 15.92032 48.35964 0.00000 0.00000 0.00000
151.50776 4.74273 20.52819 -4.93464 -27.77264 0.84809 -10.38770 15.42197 0.18824 60.65751 0.00000 0.00000
162.20486 46.69040 -35.14986 -34.36495 -19.70452 32.24332 -15.97056 13.08230 7.82976 35.84936 27.93818 0.00000
180.23747 63.67768 45.35370 15.30692 -2.60285 -13.84147 -0.34623 7.59684 2.29330 -50.49023 -12.67626 10.53301

        
        """
        """//-------METODOS PARA DETERMINAR SI UNA MATRIZ SIMETRICA ES DEFINIDA POSITIVA; no se llega a un buen resultado----------------------------------------------------------------
    public static ClassPositiva BuscarPositiva(double[][]M){
        ClassQR qr;
        double[][] M1 = new double[M.length][M[0].length];
        double epsilon=10000;
        int iter=0;
        while(iter<100&&epsilon>0.000001){
            qr = QR(M);
            M1 = ProdParalelo(qr.R, qr.Q);
            //norma de probenius o infinita M1 - M
            epsilon = NormaFrobenius(M,M1);
            //Imprimir(M1);
            System.out.println(epsilon);
            M = M1;
            iter++;
        }
        System.out.println("dentro de BuscarPositiva");
        Imprimir("M",M);
        //los elemetos de la diagonal son positivos?
        boolean flag=true;
        for(int i=0;i<M.length;i++){
            if(M[i][i]<0){
                flag =false;
            }
        }  
        return new ClassPositiva(M, flag);
    }
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    public static double NormaFrobenius(double[][]M,double[][]P){
        double norma = 0;
        for(int i=0;i<M.length;i++){
            double suma =0;
            for(int j=0;j<M[0].length;j++){
                if(i!=j){
                    suma = Math.pow(P[i][j] -M[i][j],2);
                }
            }
            norma +=suma;
        }
        return Math.sqrt(norma);
    }
    //---------------------------------------------------------------------------
    
    //---------------------------------------------------------------------------------
    public static ClassQR QR(double[][]M){
        int filas=M.length;
        int columnas = M[0].length;
        double x;
        double [][]Q =new double[filas][columnas];
        double[][]R=new double[filas][columnas];
        for (int i = 0; i < columnas; i++) {
            double res=0;
            for(int k=0;k<filas;k++){    //x =Math.sqrt(A.prodEsc(i, i));
                res+=(M[k][i]* M[k][i]);
            }
            x = Math.sqrt(res);
            R[i][i] =x;
            for (int k = 0; k < filas; k++) {
                    x = M[k][i] / R[i][i];
                    Q[k][i] = x;
            }
            for (int j = i + 1; j < columnas; j++) {//filas?gpt
                    //primero completar los R[i][j]
                    double aux=0;
                    for(int m=0;m<columnas;m++){
                            aux+=(Q[m][i]*M[m][j]);
                    }
                    R[i][j]=aux;
                    for (int k = 0; k < filas; k++) {
                            x = M[k][j] - R[i][j] * Q[k][i];
                            M[k] [j]= x;
                    }
            }
        }
        return new ClassQR(Q, R);
    }

     //-----------------------------------------------------------------------
    public static double ProductoTriple(double[][]X,double[][]M){
        double [][]Xt = Transponer(X);
        double[][]PROD=new double[X.length][M[0].length];
        for(int i=0;i<X.length;i++){
            for(int j=0;j<M[0].length;j++){
                double sum=0;
                for(int k=0;k<M.length;k++){
                    sum += X[i][k]*M[k][j];
                }
                PROD[i][j] = sum;
            }
        }
        double XAXt=0;
        for(int i=0;i<PROD.length;i++){
            for(int j=0;j<Xt[0].length;j++){
                double sum=0;
                for(int k=0;k<Xt.length;k++){
                    sum += X[i][k]*M[k][j];
                }
                XAXt = sum;
            }
        }
        return XAXt;
    }
    //----------------------------------------------------------------------
    //--------------------------------------------------------------------------
    
    //-------------------------------------------------------------------------
    //----------------------------------------------------------------------
    public static void esPositiva(){
        boolean flag=false;
        while(!flag){
            int n=50;
            GenerarMatriz(FILENAMEMATRIZ);
            for(int  i=0;i<n;i++){
                double [][]x = GenerarX(FILENAMEX);
                double p = ProductoTriple(x, A_Global);
                if(p<0){
                    flag = false;
                }else{
                    flag = true;
                }
            }
        }
    }
    public static double[][] GenerarX(String FILENAME){
        WriteDataX(FILENAME);
        double num;
        int k,n;
        double [][]X =new double[1][N_global];
        try{
            RandomAccessFile RAF = new RandomAccessFile(FILENAME, "r");
            n=(int)RAF.length()/BLOCK;
            k = 0;
            while(k<n){
                RAF.seek(k*BLOCK);
                RAF.read(RECORD);
                num = convertir();
                X[0][k] =num;
                k++;                
            }
            RAF.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
        }
        //Imprimir(X);
        return X;
    }
    MEETODOS PARA LA ESCRITURA Y LECTURA DENTRO DE LA CLASE PRINCIPAL REEMPLAZADOS POR DATASET
       public static void WriteDataMatriz(String FILENAME){
        double num;
        try{
            FileWriter fw = new FileWriter(FILENAME);
            //generando el aletario
            double max= Math.pow(10,BLOCK-1);
            double min = Math.pow(10,BLOCK-2);
            //int l = (int)(N*(N+1)/2);   //para una matriz simetrica
            for (int i=0;i<N_global*N_global;i++){
                num = Math.random()*(max - min ) + min;
                fw.write((long)num + " ");
            }
            fw.close();
        }catch(IOException e){
            e.printStackTrace();
        }
    }
    //----------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------
    public static void GenerarMatriz(String FILENAME){
        WriteDataMatriz(FILENAME);
        int n,k,i=0,j=0;
        double num;
        try{
            RandomAccessFile RAF = new RandomAccessFile(FILENAME,"r");
            n = (int)RAF.length()/BLOCK;
            k=0;
            //System.out.println("n: "+ n +" k"+ k);
            while(k<n){
                RAF.seek(k*BLOCK);      //posicionando el puntero
                RAF.read(RECORD);       //almacenando den RECORD "123 ", "432 "....
                num = convertir();      //llamando a convertir
                //almacenar en matriz
                A_Global[i][j]=A_Global[j][i]= num;      //pues sera simetrica
                if(j == i){           // j hasta j==i (la diagonal)
                    j = 0;
                    i++;
                    if(i ==N_global){
                        break;  //fila N=4 , no posible
                    }
                }else{
                    j++;
                } 
                k++;
            }
            RAF.close();
        }catch(IOException e){
            e.printStackTrace();
        }
       //Imprimir(A);
    }
    //-----------------------------------------------------------------------
    //----------------------------------------------------------------------
    public static void GenerarMatrizSimetrica(String FILENAME){
        WriteDataMatriz(FILENAME);
        int n,k,i=0,j=0;
        double num;
        try{
            RandomAccessFile RAF = new RandomAccessFile(FILENAME,"r");
            n = (int)RAF.length()/BLOCK;
            k=0;
            //System.out.println("n: "+ n +" k"+ k);
            while(k<n){
                RAF.seek(k*BLOCK);      //posicionando el puntero
                RAF.read(RECORD);       //almacenando den RECORD "123 ", "432 "....
                num = convertir();      //llamando a convertir
                //almacenar en matriz
                A_Global[i][j]= num;      //pues sera simetrica
                if(j == N_global-1){           // j hasta j==i (la diagonal)
                    j = 0;
                    i++;
                }else{
                    j++;
                } 
                k++;
            }
            RAF.close();
        }catch(IOException e){
            e.printStackTrace();
        }
       A_Global=ProdParalelo(A_Global, Transponer(A_Global));
    }
    //-----------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------
    public static double convertir(){
        String CAD = " ";
        for(int i=0;i<BLOCK-1;i++){
            CAD = CAD + (char)RECORD[i];
        }
        return Double.parseDouble(CAD);
    }
    //-----------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------
   
    //------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
    public static void WriteDataX(String FILENAME){
        double num;
        try{
            FileWriter fw = new FileWriter(FILENAME);
            for(int i=0;i<N_global;i++){
                double max = Math.pow(10, BLOCK-1);
                double min = Math.pow(10,BLOCK-2);
                num = Math.random()*(max - min)+min;
                fw.write((long)num + " ");
            }
            fw.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
        }
    }
    TESTEOS
    //GenerarMatriz(FILENAMEMATRIZ);
        //esPositiva();
        //Imprimir(A);   
         //TESTEANDO Clase QR para determinar si una matriz es definida positiva via A'=RQ
        /*
        GenerarMatriz(FILENAMEMATRIZ);
        ClassQR qr = QR(A);
        Imprimir(qr.Q);
        Imprimir(qr.R);
        Imprimir(ProdParalelo(qr.Q, qr.R));
        double t = NormaFrobenius(A, qr.Q);
        System.out.println(t);
        System.out.println(BuscarPositiva(A));
        //para buscar positiva dado que ya es simetrica:boolean flag =false;
        ClassPositiva positiva=new ClassPositiva(A,true);
        while(!flag){
            GenerarMatriz(FILENAMEMATRIZ);
            positiva = BuscarPositiva(A);
            flag = pos.flag;
        }
        A =positiva.M;
        */
        /*  //En lugar de usar QR se hizo A= AA^t y se tiene una simetrica
        GenerarMatrizSimetrica(FILENAMEMATRIZ);
        Imprimir(A);
        double [][] G = CholeskySerial(A);
        Imprimir(G);
            //TESTEANDO si cholesky serial realemente descompuso A en G (triangular inferior)
        /*
        Imprimir(ProdParalelo(G, Transponer(G)));
        */
        // TESTEANDO la obtencion de los bloques para el tratamiento paralelo, se otendran solo desde la diagonal hacia↓
        // es decir (k)(k+1)/2  BLOQUES
        /*
        ObtenerABloques();
        for(int i=0;i<ABloques.length;i++){
            for(int j=0;j<ABloques[0].length;j++){
                
                Imprimir(ABloques[i][j]); //
            }
        } 
        */
        //TESTEANDO DataSetCholesky
        /* 
        DataSetCholesky data = new DataSetCholesky(FILENAMEMATRIZ, N,N,BLOCK,"FileWriter");//2modos de CrearData PrintWriter(Scanner) o FileWriter(RandomAccessFile)
        A = DataSetCholesky.ReadDataRAF(FILENAMEMATRIZ);
        DataSetCholesky.WriteDataRAF("DATACholeskyBloquesCOPIA.TXT");
        Imprimir(DataSetCholesky.ReadDataRAF("DATACholeskyBloquesCOPIA.TXT"));
        Imprimir(A);
        */
        //Imprimir(A);
        //double[][] X=GenerarX(FILENAMEX);
        //double nTest=ProductoTriple(X, A);
        //System.out.println(nTest); se estaba testeando estos metodos
       // TESTEANDO CholeskyBloques        falta revisar
       
       /*
       for(int i=0;i<LBloques.length;i++){
            for(int j=0;j<LBloques[0].length;j++){
                Imprimir(LBloques[i][j]); //
            }
        }  
        */
            // TESTEANDO las triangulares 
        /* 
        ObtenerABloques();
        double [][] L11= CholeskySerial(ABloques[0][0]);
        Imprimir(L11);
        ImprimirInversa(InversaTriangularInferior(L11));
        Imprimir(ProdParalelo(L11, InversaTriangularInferior(L11)));
        Imprimir(Transponer(L11));
        ImprimirInversa(InversaTriangularSuperior(Transponer(L11)));
        Imprimir(ProdParalelo(Transponer(L11), InversaTriangularSuperior(Transponer(L11))));
        */
        /*     FUNCIONA OBTENIENDO BLOQUES AFUERA
        GenerarMatrizSimetrica(FILENAMEMATRIZ);
        Imprimir("A",A_Global);
        ObtenerABloques();
       for(int i=0;i<ABloques.length;i++){
         for(int j=0;j<ABloques[0].length;j++){
                
            Imprimir("ABloques",ABloques[i][j]); 
        }
        }
    */
    /* 
        DataSetCholesky dataGlobal = new DataSetCholesky("DATASETCHOLESKY.TXT", N_global, N_global, BLOCK-1, "FileWriter");
        A_Global = dataGlobal.ReadDataRAF();
        Imprimir("A con DataSetCholseky", A_Global);
        ObtenerABloques();
        for(int i=0;i<ABloques.length;i++){
            for(int j=0;j<ABloques[0].length;j++){
                Imprimir("ABloques",ABloques[i][j]); 
            }
        }
        Imprimir("ABloques[0][0]", ABloques[0][0]);
        double [][] L11= CholeskySerial(ABloques[0][0]);
        Imprimir("L11",L11);
        double[][]L21=TRSM(ABloques[1][0], L11);
        Imprimir("L21 = TRSM(ABloques[2][1] L11)",L21);
        Imprimir("ABloques[3][1]", ABloques[2][0]);
        double[][]L31=TRSM(ABloques[2][0], L11);
        Imprimir("L31 =TRSM ABloques[3][1]", L31);
        Imprimir("ABloques[2][2]",ABloques[1][1]);
        double[][]A22p=SYRK(ABloques[1][1], L21,L21);
        Imprimir("SYRK A22=A22 -L21L21^t",A22p);
        double[][]A32p=SYRK(ABloques[2][1], L31, L21);
        Imprimir("SYRK A32=A32 -L31 L21^t", A32p);
        double[][]A33p =SYRK(ABloques[2][2],L31,L31);
        Imprimir("SYRK A33 =A33 -L31 L31^t", A33p);
        double[][]L22=CholeskySerial(A22p);
        Imprimir("L22 ", L22);
        double[][]L32=TRSM(A32p, L22);
        Imprimir("L32", L32);
        double[][]A33pp =SYRK(A33p, L32,L32);
        Imprimir("A33pp", A33pp);
        double[][]L33 = CholeskySerial(A33pp);
        Imprimir("L33", L33);
        double [][]LBloquesF =new double[N_global][N_global];
        CopiarBloque(LBloquesF, L11, 0, 0);  // L11 -> posición (0, 0)
        CopiarBloque(LBloquesF, L21, 4, 0);  // L21 -> posición (4, 0)
        CopiarBloque(LBloquesF, L31, 8, 0);  // L31 -> posición (8, 0)
        CopiarBloque(LBloquesF, L22, 4, 4);  // L22 -> posición (4, 4)
        CopiarBloque(LBloquesF, L32, 8, 4);  // L32 -> posición (8, 4)
        CopiarBloque(LBloquesF, L33, 8, 8);
        Imprimir("LGlobal manual ", LBloquesF);
        //corroborando con choleskySerial
        double[][]ASerial = Copiar(A_Global);
        Imprimir("factorizacion cholesky Serial",CholeskySerial(ASerial));
    */        
        //LOS RESULTADOS INDICAN QUE TODO ES CORRECTO LO CUAL ES GENIAL GENIAL GENIAL !!!
            //EL TESTEO ANTERIOR SE HIZO DE FORMA NO ITERADA PASO A PASO y RESULTA EN LA OBTENCION DE LOS LBLOQUES
            //AHORA SE INTENTA LA FORMA ITERADA y RESULTA UN ERROR en el manejo de los indices seguramente FALTA CORREGIR ELLO
            //UNICAMENTE SE REQUIERE CORREGIR EL MANEJO DE LOS INDICES DENTRO DE LOS FOR DE Choleskybloques()
    
        //DataSetCholesky dataGlobal = new DataSetCholesky("DATASETCHOLESKY.TXT", 12, 12, 3, "FileWriter");
    /*  A_Global = dataGlobal.ReadDataRAF();
        Imprimir("A con DataSetCholseky", A_Global);
        ObtenerABloques();
       for(int i=0;i<ABloques.length;i++){
            for(int j=0;j<ABloques[0].length;j++){
                Imprimir("ABloques",ABloques[i][j]); 
            }
        } 
    */
          //CORREGIDO , funciona perfectamente, solo se trataba de un tema de indices, se corroboran que los 3 resultados
        // Cholesky manual*   cholesky serial   choleskyBloques  DAN LA MISMA TRIAANGULAR INFERIOR
    /*  for(int i=0;i<ABloques.length;i++){
            for(int j=0;j<ABloques[0].length;j++){
                Imprimir("LBloques",LBloques[i][j]); 
            }
        }
    */
    """