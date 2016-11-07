public class HMM {
    public  double[][] A;
    public  double[][] B;
    private  double[] Pi;
    public static int N = 5;
    public int M = 9;
    private static double[] scalee;
    private static double[] scale;
    HMM() {
        A = new double[N][N];
        B = new double[N][M];
        Pi = new double[N];

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if(i==j)
                    A[i][j] = (double)1.0;
                else
                    A[i][j]=(double) 0.001;
            }
            Pi[i]= (double) 1.0/N;
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                B[i][j] = 0.001;
            }
        }
        B[1][0]=0.5;
        B[1][2]=0.5;
        B[2][7]=1.0/3;
        B[2][3]=1.0/3;
        B[2][5]=1.0/3;


        for (int i = 0; i < M; i++) {
            B[0][i]=1.0/M;
            B[4][i]=(double)1.0/M;
            B[3][i]=(double)1.0/M;
        }
    }

    public void train(int[] O) {
        if (O.length == 1) {
            return;
        }
        double lamda = 100;
       
        double delta = 100;
        int seqcount = O.length;
        scalee = new double[seqcount];
        do {

            //init and calc alpha
            scalee =new double[seqcount];
            double[][] alpha = new double[seqcount][N];
            for (int i = 0; i < N; i++) {
                alpha[0][i] = Pi[i] * B[i][O[0]];
                scalee[0] += alpha[0][i];
            }
            for (int i = 0; i < N; i++) {
                alpha[0][i] = alpha[0][i] / scalee[0];
            }
            for (int j = 1; j < seqcount; j++) {
                for (int i = 0; i < N; i++) {
                    alpha[j][i] = calcalpha(i, j, A, B, O, alpha, N);
                    scalee[j] += alpha[j][i];
                }
                for (int i = 0; i < N; i++) {
                    alpha[j][i] = alpha[j][i] / scalee[j];
                }

            }

            //init and calc beta
            double[][] beta = new double[seqcount][N];

            for (int i = 0; i < N; i++) {
                beta[seqcount - 1][i] = 1 / scalee[seqcount - 1];
            }
            for (int j = seqcount - 2; j > -1; j--) {
                for (int i = 0; i < N; i++) {
                    beta[j][i] = calcbeta(i, j, A, B, O, beta, N) / scalee[j];
                }
            }


            //init and calc gamma
            double[][][] digamma = new double[seqcount][N][N];
            for (int k = 0; k < seqcount - 2; k++) {
                double below = 0;
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        digamma[k][i][j] = alpha[k][i] * beta[k + 1][j] * A[i][j] * B[j][O[k + 1]];
                        below += digamma[k][i][j];
                    }
                }
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        digamma[k][i][j] = digamma[k][i][j] / below;
                    }
                }
            }

            //init and calc digamma
            double[][] gamma = new double[seqcount][N];
            for (int i = 0; i < seqcount - 1; i++) {
                double below = 0;
                for (int j = 0; j < N; j++) {
                    for (int k = 0; k < N; k++) {
                        below += alpha[i][j] * beta[i + 1][k] * A[j][k] * B[k][O[i + 1]];
                    }
                }
                for (int j = 0; j < N; j++) {
                    gamma[i][j] = 0;
                    for (int k = 0; k < N; k++) {
                        digamma[i][j][k] = alpha[i][j] * A[j][k] * B[k][O[i + 1]] * beta[i + 1][k] / below;
                        gamma[i][j] += digamma[i][j][k];
                    }
                }
            }
            double lamda2 = 0;

            A = updatea(digamma, gamma, seqcount, N, N);
            B = updateb(digamma, gamma, seqcount, O, M, N);
            Pi = updatepi(gamma, N);

            lamda2 = calclamda(alpha, seqcount, N);
            delta = Math.abs(lamda2 - lamda);
            lamda = lamda2;
            
            

        }
        while (delta > 0.01);
    }

    public double[] predict(int[] O) {                     

        int seqcount= O.length;
        scale = new double[seqcount];
        double [] result=new double[M];
        double[][] alpha = new double[seqcount][N];
        for (int i = 0; i < N; i++) {
            alpha[0][i] = Pi[i] * B[i][O[0]];
            scale[0] += alpha[0][i];
        }
        for (int i = 0; i < N; i++) {
            alpha[0][i] = alpha[0][i] / scale[0];
        }
        for (int j = 1; j < seqcount; j++) {
            for (int i = 0; i < N; i++) {
                alpha[j][i] = calcalpha(i, j, A, B, O, alpha, N);
                scale[j] += alpha[j][i];
            }
            for (int i = 0; i < N; i++) {
                alpha[j][i] = alpha[j][i] / scale[j];
            }

        }
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                double temp=0;
                for (int k = 0; k < N; k++) {
                    temp+=A[k][j]*alpha[seqcount-1][k];
                }
                temp*=B[j][i];
                result[i]+=temp;
            }
        }
        double sumtemp=0;
        for (int i = 0; i < M; i++) {
            sumtemp+=result[i];
        }
        for (int i = 0; i < M; i++) {
            result[i]/=sumtemp;
        }
        return result;
    }

    public  double Probability_HMM_Obsequence(int[] O){   
    	int seqcount=O.length;
        double result=(float) 0.0;

        double []preprob=new double[N];
        for (int i = 0; i < N; i++) {
            preprob[i]=Pi[i]*B[i][O[0]];
        }
        for (int i = 1; i < seqcount; i++) {
            double []probtmp=new double[N];
            for (int ii = 0; ii < N; ii++) {
                probtmp[ii]= (double) 0.0;
            }
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < N; k++) {
                    probtmp[j]+=preprob[k]*A[k][j];
                }
                probtmp[j]=probtmp[j]*B[j][O[i]];
            }
            preprob=probtmp;

        }

        for (int i = 0; i < N; i++) {
            result+=preprob[i];
        }
        return result;

    }

    private static double calclamda(double[][] alpha, int seqcount , int  N) {
        double result=0;
        for (int i = 0; i < seqcount; i++) {
            result+=Math.log(scalee[i]);
        }
        return result;
    }

    private static double[] updatepi(double[][] gamma,int N) {
        double [] result=new double[N];
        for (int i = 0; i < N; i++) {
            result[i] = (gamma[0][i]*.99999+.00001);
        }
        return result;
    }

    private static double[][] updateb(double[][][] digamma, double[][] gamma, int seqcount, int[] O, int Bcol, int Brow) {
        double [][]result=new double[Brow][Bcol];
        for (int i = 0; i < Brow; i++) {
            for (int j = 0; j < Bcol; j++) {
                double above=0;
                double below=0;
                for (int k = 0; k < seqcount-1; k++) {
                    above+=gamma[k][i]*(O[k]==j?1:0);
                    below+=gamma[k][i];
                }
                result[i][j]= (double) ((above/below)*.99999+.00001);
            }
        }
        return result;
    }

    private static double[][] updatea(double[][][] digamma, double[][] gamma, int seqcount, int Acol, int Arow) {
        double [][]result=new double[Arow][Acol];
        for (int i = 0; i < Arow; i++) {
            for (int j = 0; j < Acol; j++) {
                double above=0;
                double below=0;
                for (int k = 0; k < seqcount-1; k++) {
                    above+=digamma[k][i][j];
                    below+=gamma[k][i];
                }
                result[i][j]= ((above/below)*.99999+.00001);
            }
        }
        return result;
    }

    private static double calcalpha(int i, int t, double[][] A, double[][] B, int[] O, double[][] alpha, int N) {
        double result=0;

        for (int j = 0; j < N; j++) {
            result+=A[j][i]*alpha[t-1][j];
        }
        result*=B[i][O[t]];
        return result;
    }

    private static double calcbeta(int i, int t, double[][] A, double[][] B, int[] O,double[][] beta,int N) {
        double result=0;
        for (int j = 0; j < N; j++) {
            result+=beta[t+1][j]*B[j][O[t+1]]*A[i][j];
        }
        return result;

    }


}

