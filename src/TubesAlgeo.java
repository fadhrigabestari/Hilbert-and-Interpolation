import java.util.*;
import java.text.*;
import java.lang.*;
import java.io.*;

public class TubesAlgeo{
	//Inisialisasi cara input ke program
	private static Scanner in = new Scanner(System.in);
	//Format penulisan nilai 3 di belakang koma
	private static DecimalFormat df = new DecimalFormat("#.###");
	//Inisialisasi nilai epsilon dalam program
	private static double EPSILON = 1e-24;

	//Mengecek input N
  public static boolean IsNValid(int n) {
    return(n >= 2); //True saat n untuk interpolasi valid, atau lebih dari 1
  }

	//Mengecek input antara Gap
  public static boolean IsGapValid(float start, float end) {
    return(end > start); //True saat nilai end > start untuk interpolasi
  }

	//Mengecek input antara min max
  public static boolean IsInValid(int i, int min, int max) {
    return(i >= min && i <= max); //True saat nilai input diantara min dan max
  }

	//Mengecek index
  public static boolean IsIdxValid(int m, int n) {
    return(m > 0 && n > 0); //True saat nilai index tidak negatif dan tidak 0
  }

	//Fungsi pangkat
  public static float sqr(float x, int n) {
		int i;
		float c = 1;
		for(i = 0; i < n; i++) {
			c = x * c; //Mengalikan nilai x sebanyak n kali yang diletakkan di c
		}
		return c; //Mengembalikan nilai c
	}

	//Membaca matrix dengan baris n dan kolom m
  public static float[][] BacaMatriks(int n, int m) {
		int i, j;
		float[][] M;
		M = new float[n][m];
		for (i=0; i<n; i++) {
			for (j=0; j<m; j++) {
				M[i][j] = in.nextFloat();
			}
		}
		return M;
	}

  //Membaca array dengan index n
	public static float[] BacaArray(int n) {
		int i;
		float[] T = new float[n];
		for (i=0; i<n; i++) {
				T[i] = in.nextFloat();
		}
		return T;
	}

	//Menulis matrix biasa
	public static void TulisMatriks(float[][] M) {
		int i, j;

		for (i=0; i<M.length; i++) {
			for (j=0 ; j<M[i].length; j++) {
				System.out.print(df.format(M[i][j]) + "\t");
			}
			System.out.println();
		}
	}

	//Menulis augmented matrix
	public static void TulisMatriksDiperbesar (float[][] M, float[] T) {
		int i,j;
		for (i=0; i<M.length; i++) {
			for (j=0 ; j<M[i].length; j++) {
				System.out.print(df.format(M[i][j]) + "\t");
			}
			System.out.print("|\t" + df.format(T[i]));
			System.out.println();
		}
	}

 	//Mengcopy matrix dan array ke matrix dan array baru
  public static void CopyMatrix(float[][] Msrc, float[][] Mdest, float[] Tsrc, float[] Tdest) {
		for (int i=0; i<Msrc.length;i++) {
			for (int j=0; j<Msrc[0].length; j++) {
				Mdest[i][j] = Msrc[i][j];
			}
			Tdest[i] = Tsrc[i];
		}
	}

	//Melaksanakan fungsi gauss
	public static void Gauss(float[][] M, float[] T) {
		int i, j, k, l, m, n, o;
		int max;
		float[] temparr;
		float tempfloat;

		m = M.length; n = M[0].length;

		if (m>n) {o=n;} else {o=m;};
		l=0;
		for (k=0; k<o; k++) {
			//MENCARI PIVOT
			max = k;
			for (i=k; i<m; i++) {
				if (Math.abs(M[i][k])>Math.abs(M[max][k])) {
					max = i;
				}
			}

			//MENUKAR BARIS DENGAN PIVOT
			temparr = M[k];
			M[k] = M[max];
			M[max] = temparr;

			tempfloat = T[k];
			T[k] = T[max];
			T[max] = tempfloat;

			l = k;
			while (l<n && M[k][l]==0) {
				l++;
			}
			if (l<n) {
				//PENYEDERHANAAN
				tempfloat = M[k][l];
				if (Math.abs(tempfloat) > EPSILON) {
					T[k] = T[k] /  tempfloat;
					for (j=k; j<n; j++) {
						M[k][j] = M[k][j] / tempfloat;
					}
				}
				//ELIMINASI
				for (i = k+1; i < m; i++) {
					tempfloat = M[i][l] / M[k][l];
					T[i] = T[i] - (tempfloat * T[k]);
					for (j=0; j<n; j++) {
						M[i][j] = M[i][j] - (tempfloat * M[k][j]);
					}
				}
			}
		}
	}

	//Melaksanakan fungsi gauss jordan
	public static void GaussJordan(float[][] M, float[] T) {
		int i, j, k, l, m, n, o;
		int max;
		float[] temparr;
		float tempfloat;

		m = M.length; n = M[0].length;

		if (m>n) {o=n;} else {o=m;};
		l=0;
		for (k=0; k<o; k++) {
			//MENCARI PIVOT
			max = k;
			for (i=k; i<m; i++) {
				if (Math.abs(M[i][k])>Math.abs(M[max][k])) {
					max = i;
				}
			}

			//MENUKAR BARIS DENGAN PIVOT
			temparr = M[k];
			M[k] = M[max];
			M[max] = temparr;

			tempfloat = T[k];
			T[k] = T[max];
			T[max] = tempfloat;

			l = k;
			while (l<n && M[k][l]==0) {
				l++;
			}
			if (l<n) {
				//PENYEDERHANAAN
				tempfloat = M[k][l];
				if (Math.abs(tempfloat) > EPSILON) {
					T[k] = T[k] /  tempfloat;
					for (j=k; j<n; j++) {
						M[k][j] = M[k][j] / tempfloat;
					}
				}
				//ELIMINASI
				for (i = k+1; i < m; i++) {
					tempfloat = M[i][l] / M[k][l];
					T[i] = T[i] - (tempfloat * T[k]);
					for (j=0; j<n; j++) {
						M[i][j] = M[i][j] - (tempfloat * M[k][j]);
					}
				}
			}
		}

		//ELIMINASI GAUSS JORDAN
		boolean isfound;
		for (k=0; k<n; k++) {
			i=m-1;
			while (i>=0 && M[i][k]!=1) {
				i--;
			}

			if (i!=0 && i>=0) {
				for (l=0; l<i; l++) {
					tempfloat = M[l][k] / M[i][k];
					T[l] = T[l] - (tempfloat * T[i]);
					for (j=0; j<n; j++) {
						M[l][j] = M[l][j] - (tempfloat * M[i][j]);
					}
				}
			}
		}
	}

	//Menyimpan jawaban hasil gauss di file external sebagai pilihan
  public static void SimpanJawabanG(float[][] Mgauss, float[] Tgauss, float[][] Mori, float[] Tori, String[] jawaban) {
		int i, j;
		char N;

		System.out.println("Apakah anda ingin menyimpan kalkulasi ke file? (Y/N)");
		N = in.next().charAt(0);

		if (N=='Y' || N=='y') {
			System.out.println("Masukan nama file dimana kalkulasi ingin disimpan : ");
			in.nextLine();
			String filename = in.nextLine() ;
			try {
				PrintWriter writer = new PrintWriter(filename, "UTF-8");
				//MENULIS AUGMENTED MATRIX
				writer.println("Augmented Matrix");
				for (i=0; i<Mori.length; i++) {
					for ( j=0; j<Mori[0].length; j++) {
						writer.print(df.format(Mori[i][j]) + "\t");
					}
					writer.print("|\t" + df.format(Tori[i]));
					writer.println("");
				}
				//MENULIS HASIL GAUSS
				writer.println("Echelon Form Matrix");
				for (i=0; i<Mgauss.length; i++) {
					for (j=0; j<Mgauss[0].length; j++) {
						writer.print(df.format(Mgauss[i][j]) + "\t");
					}
					writer.print("|\t" + df.format(Tgauss[i]));
					writer.println("");
				}
				//MENULIS JAWABAN
				if (jawaban.length==1) {
					writer.println(jawaban[0]);
				} else {
					writer.println("Solusi : ");
					for (i=0; i<jawaban.length; i++) {
						writer.println("X" + i + " = " + jawaban[i]);
					}
				}
				writer.close();
			} catch (IOException e) {
			}
      System.out.println("File berhasil disimpan");
		}
	}

	//Melaksanakan fungsi gauss jordan dan menuliskannya sebagai output
  public static void JordanG(float[][] M, float[] T, float[][] Mori, float[] Tori) {
		GaussJordan(M,T); //Melaksanakan procedure gauss jordan
		System.out.println("\nReduced Echelon Form Matrix\n");
		TulisMatriksDiperbesar(M,T);
		System.out.println();
	}

	//Melaksanakan fungsi gauss dan menuliskannya sebagai output
	public static void OrdG(float[][] M, float[] T, float[][] Mori, float[] Tori) {
		Gauss(M,T);
		System.out.println("Echelon Form Matrix\n");
		TulisMatriksDiperbesar(M,T);
		System.out.println();
	}

	//Membaca file external dengan format matrix untuk gauss dan gauss jordan
  public static float[][] BacaFileMatriks() {
		System.out.println("Masukan Nama File yang Berisi Augmented Matrix : ");
    in.nextLine();
    String filename = in.nextLine();
		String str;
		int row, col, i, j;
		float[][] failed = new float[0][0];
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			try {
				List<String> list = new ArrayList<String>();
				while((str = br.readLine()) != null) {
					if (!str.equals("")) {
						list.add(str);
					}
				}

				String[] temparr = list.toArray(new String[0]);
				String[] parts = temparr[0].split("\t");
				row = temparr.length;
				col = parts.length-1;

				float[][] M = new float[row][col];

				for (i=0; i<row; i++){
		      parts = temparr[i].split("\t");
		      for (j=0; j<parts.length-1; j++) {
		        if (parts[j].equals("|")) {
		          M[i][j] = Float.parseFloat(parts[j+1]);
		        } else {
		          M[i][j] = Float.parseFloat(parts[j]);
		        }
		      }
		    }
				return M;
			} catch (IOException e) { return failed;}
		} catch (FileNotFoundException ex) {
      System.out.println("File tidak ditemukan!");
      return BacaFileMatriks();
    }
	}

	//Proses penyelesaian soal gauss, disertakan outputnya
  public static void SolusiG(float[][] M, float[] T, float[][] Mori, float[] Tori) {
		boolean iszero = true;
		int i, j, k, l;
		float[] X = new float[M[0].length];
		String[] Y = new String[M[0].length];
		float sum, tempfloat;
		int countnonzero;
		String Z = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
		int N = Z.length();

		Random r = new Random();
    GaussJordan(M,T);
		for (j=0; j<M[0].length; j++) {
			if (!(Math.abs(M[M.length-1][j]) < EPSILON)) {
				iszero = false;
			}
		}

		if (iszero) {
			if (!((Math.abs(T[T.length-1])) < EPSILON)) {
				String[] gagal = new String[1];
				gagal[0] = "Tidak ada Solusi";
				System.out.println("Tidak ada Solusi");
				SimpanJawabanG(M, T, Mori, Tori, gagal);
			} else {
				//Mengisi semua tabel solusi parameter dengen empty string
				for (i=0; i<Y.length;i++) {
					Y[i] = "";
				}
				//Mengisi semua tabel solusi dengan angka NaN
				for (i=0; i<X.length; i++) {
					X[i] = Float.NaN;
				}
				//Mendapatkan solusi nonparameter (jika ada)
				for (i = M.length-1; i>=0; i--) {
					sum = 0;
					l = 0;
					while (l<M[0].length && M[i][l]==0) {
						l++;
					}
					if (l<M[0].length) {
						for (j = l+1; j < M[0].length; j++) {
							if (M[i][j] != 0) {
								sum = sum + (M[i][j] * X[l]) ;
							}
						}
						X[l] = (T[i] - sum) / M[i][l];
					}
				}
				//Mencari solusi parameter;
				for (i=M.length-1; i>=0; i--){
					sum = 0;
					l = 0;
					while (l<M[0].length && M[i][l]==0) {
						l++;
					}
					if (l<M[0].length) {
						if (!Float.isNaN(X[l])) {
							Y[l] = String.format("%.3f" ,(X[l]));
						} else {
							for (j=M[0].length-1; j>l;j--) {
								if (M[i][j] != 0) {
									if (Float.isNaN(X[j])) {
										if (Y[j]=="") {
											Y[j] = Character.toString(Z.charAt(r.nextInt(N)));
										}
										if (M[i][j]>0) {
											if (Y[l]=="") {
												Y[l] = ("-"+ M[i][j]+ "" +Y[j]);
											} else {
												Y[l] += " - "+ M[i][j]+ "" + Y[j];
											}
										} else {
											if (Y[l]=="") {
												Y[l] = "" + (-1*M[i][j]) + "" +Y[j];
											} else {
												Y[l] += " + "+ (-1*M[i][j])+ ""  + Y[j];
											}
										}
									} else {
										sum = sum + X[j];
									}
								}
							}
							X[l] = (T[i] - sum) / M[i][l];
							if (X[l] != 0 && !Float.isNaN(X[l])) {
								if (X[l] > 0)
									Y[l] += " + " + X[l];
								else
									Y[l] += " - " + (-1*X[l]);
							}
						}
					}
				}
				for (i=0; i<Y.length; i++) {
					System.out.println("X" + i + " = " + Y[i]);
				}
				SimpanJawabanG(M, T, Mori, Tori, Y);
			}
		} else {
			for (i = M.length-1; i>=0; i--) {
				sum = 0;
				l = 0;
				while (l<M[0].length && M[i][l]==0) {
					l++;
				}
				if (l<M[0].length) {
					for (j = l+1; j < M[0].length; j++) {
						if (M[i][j] != 0) {
							sum = sum + (M[i][j] * X[j]) ;
						}
					}
					X[l] = (T[i] - sum) / M[i][l];
				}
			}

			for (i=0; i<T.length; i++) {
				Y[i] = String.format("%.3f" ,(X[i]));
				System.out.println("X" + i + " = " + Y[i]);
			}
			SimpanJawabanG(M, T, Mori, Tori, Y);
		}
	}

	//Memisahkan augmented matrix saat baca file external
  public static void SplitAugmented(float[][] Maug, float[][] M, float[] T) {
    int i,j;
    for (i=0; i<Maug.length; i++) {
      for (j=0; j<Maug[0].length; j++) {
        if (j==Maug[0].length-1) {
          T[i] = Maug[i][j];
        } else {
          M[i][j] = Maug[i][j];
        }
      }
    }
  }

	//Menu untuk gauss, berisikan semua fungsi dan prosedur gauss
  public static void GaussMenu() {
		int n,m,i;
    float[][] M;
    float[] T;
		System.out.println("Tipe pembacaan data : ");
		System.out.println("1. Pembacaan data dari pengguna");
		System.out.println("2. Pembacaan data dari file external");
		System.out.print("Tipe pembacaan data yang ada inginkan : ");
		do {
			i = in.nextInt();
		}while(!IsInValid(i,1,2));

		if (i == 1) {
			do {
				System.out.print("Masukkan jumlah variabel : ");
				n = in.nextInt();
				System.out.print("Masukkan jumlah persamaan : ");
				m = in.nextInt();
			}while(!IsIdxValid(m,n)); //Agar n dan m diatas 0

			System.out.println("Masukkan koefisien persamaan : ");
			M = BacaMatriks(m, n); //Membaca Matriks dengan m baris dan n kolom
			System.out.println("Masukkan solusi persamaan : ");
			T = BacaArray(m); //Membaca solusi persamaan
		} else {
      float[][] Maug = BacaFileMatriks();
      M = new float[Maug.length][Maug[0].length-1];
  		T = new float[Maug.length];
      SplitAugmented(Maug, M, T);
		}

    float[][] Mori = new float[M.length][M[0].length];
    float[] Tori = new float[M.length];
		CopyMatrix(M, Mori, T, Tori);

		System.out.println();
		System.out.println("Augmented Matrix\n");
		TulisMatriksDiperbesar(M,T); //Menulis matriks
		System.out.println();

		System.out.println("Tipe Gauss : ");
		System.out.println("1. Penyederhanaan SPL dengan metode Gauss");
		System.out.println("2. Penyederhanaan SPL dengan metode Gauss Jordan");
		System.out.print("Tipe Gauss yang ingin digunakan : ");
		do {
			i = in.nextInt();
		}while(!IsInValid(i,1,2)); //Agar input diantara 1 dan 2

		if(i == 1) {
			OrdG(M,T,Mori,Tori); //Gauss biasa untuk input 1
		}
		else {
			JordanG(M,T,Mori,Tori); //Gauss Jordan untuk input 2
		}
		SolusiG(M,T,Mori,Tori);
		System.out.println();
	}

	//Mengisi matrix dengan input interval [a..b] dengan n
  public static float[][] IsiMIntvl(float start, float end, int n) {
		int i,j;
		float h = (end - start) / n; //h sebagai nilai jarak antar x
		float[][] M = new float[n+1][n+1]; //Jumlah kolom dan baris adalah + 1 dari n (mulai dari 0)
		for(i = 0; i <= n; i++) {
			for(j = 0; j <= n; j++) {
				M[i][j] = sqr(start,j); //Memangkatkan angka start sebanyak j kali (sesuai kolom)
			}
			start = start + h; //Menambahkan nilai start dengan selisih antar x
		}
		return M; //Mengembalikan nilai M
	}

	//Mengisi matrix dengan x inputan manual dari user
  public static float[][] IsiMInterpolasi(float[] T, int n) {
		float[][] M = new float[n][n]; //Jumlah kolom dan baris adalah n + 1
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < n; j++) {
				M[i][j] = sqr(T[i],j); //Isi dengan pangkat j dari Array T (nilai x)
			}
		}
		return M; //Mengembalikan M
	}

	//Menulis bentuk polinomial interpolasi setelah pendekatan f(x)
  public static void TulisPolinom(float[] T) {
		System.out.print("p(x) = " + df.format(T[0]) + " + ");
		for(int i = 1; i < T.length; i++) {
			if(T[i] != Float.NaN) {
				System.out.print(df.format(T[i]) + "x^" + i);
				if(i < T.length - 1) {
					System.out.print(" + ");
				}
			}
		}
		System.out.println("\n");
	}

	public static float xfungsi (float x) {
    float result;
    result = (1 + (float)Math.sqrt(x) + sqr(x,2));
    if(x >= 0) {
      for(int i = 0; i < Math.abs(x); i++) {
        result = result / 2.71828f;
      }
    } else {
      for(int i = 0; i < Math.abs(x); i++) {
        result = result * 2.71828f;
      }
    }
    return result;
  }

	//Memasukkan nilai x ke dalam f(x)
  public static double IsiX(float[] T, float X) {
	 	double result = 0;
		for(int i = 0; i < T.length; i++) {
			result = result + (double)(T[i] * sqr(X,i));
		}
		return result; //Mengembalikan result
	}

	//Menulis nilai IsiX ke layar sebagai output
  public static void TulisXFx(float[] T, float X) {
		System.out.print("Nilai f(" + X + ") : ");
		double result = IsiX(T,X);
		System.out.println(df.format(result)); //Ditulis 3 angka dibelakang koma
		System.out.println();
	}

	//Memberi opsi kepada user untuk menyimpan data ke dalam file external
  public static void SimpanJawabanI(float[][] M, float[] T, float[][] Mori, float[] Tori, float X, int type) {
		int i, j;
		char N;

		System.out.println("Tipe penyimpanan data : ");
		System.out.println("1. Menyimpan jawaban ke file external");
		System.out.println("2. Tidak menyimpan jawaban ke file external");
		System.out.print("Tipe penyimpanan data yang diinginkan : ");
		do {
			i = in.nextInt();
		}while(!IsInValid(i,1,2));

		if(i == 1) {
			System.out.print("Masukkan nama file dimana kalkulasi ingin disimpan : ");
			in.nextLine();
			String filename = in.nextLine();
			try {
				PrintWriter writer = new PrintWriter(filename, "UTF-8");
				//MENULIS AUGMENTED MATRIX
				writer.println("Augmented Matrix");
				writer.println();
				for(i = 0; i < Mori.length; i++) {
					for(j = 0; j < Mori[0].length; j++) {
						writer.print(df.format(Mori[i][j]) + "\t");
					}
					writer.print("|\t" + df.format(Tori[i]));
					writer.println("\n");
				}
				//MENULIS HASIL GAUSS
				writer.println("Echelon Form Matrix");
				writer.println();
				for (i=0; i<M.length; i++) {
					for (j=0; j<M[0].length; j++) {
						writer.print(df.format(M[i][j]) + "\t");
					}
					writer.print("|\t" + df.format(T[i]));
					writer.println("\n");
				}
				//MENULIS POLINOM INTERPOLASI
				writer.print("p(x) = " + df.format(T[0]) + " + ");
				for(i = 1; i < T.length; i++) {
					writer.print(df.format(T[i]) + "x^" + i);
					if(i < T.length - 1) {
						writer.print(" + ");
					}
				}
				writer.println();
				if(type == 1) {
					writer.print("Nilai f(" + X + ") : ");
					double result = IsiX(T,X);
					writer.println(df.format(result));
				}
				writer.close();
			} catch (IOException e) {}
		} else {
			System.out.println();
		}
	}

	//Menyelesaikan kalkulasi estimasi x dalam f(x)
  public static int ProcessXFx(float [][] M, float[] T, float[][] Mori, float[] Tori, float X) {
		int c;
		System.out.println("Estimasi nilai x dalam f(x) :");
		System.out.println("1. Ya, estimasi nilai x dalam f(x)");
		System.out.println("2. Tidak, sudahi estimasi nilai x");
		System.out.print("Pilihan yang Anda inginkan : ");
		do {
			c = in.nextInt();
		}while(!IsInValid(c,1,2));
		System.out.println();
		if(c == 1) {
			System.out.print("Masukkan nilai X yang ingin di-estimasi : ");
			X = in.nextFloat();
			TulisXFx(T,X);
			SimpanJawabanI(M,T,Mori,Tori,X,1);
		}
		return c;
	}

	//Menyelesaikan kalkulasi interpolasi secara menyeluruh
  public static void SolusiI(float[][] M, float[]T, float[][] Mori, float[] Tori, float X) {
		int i;
		System.out.println("\nAugmented Matrix\n");
    TulisMatriksDiperbesar(M,T); //Menulis bentuk augmented matriks
    System.out.println("\n");
    GaussJordan(M,T); //Menggunakan fungsi gauss jordan
    System.out.println("Echelon Form Matrix\n");
    TulisMatriksDiperbesar(M,T);
    System.out.println();
    TulisPolinom(T); //Menulis polinom jawaban
		if(ProcessXFx(M,T,Mori,Tori,X) == 1) {
			do {
				i = ProcessXFx(M,T,Mori,Tori,X);
			} while(i == 1);
		}
		else {
			SimpanJawabanI(M,T,Mori,Tori,X,2);
		}
  }

	//Menu dengan isi interpolasi untuk input [a..b] dengan n
	public static void GapInterpolate() {
		float start, end; int n; float X = 0; int i;
		System.out.print("Masukan selang antar titik : ");
		do {
			start = in.nextFloat();
			end = in.nextFloat();
		}while(!IsGapValid(start,end)); //Agar start < end

		System.out.print("Masukkan jumlah n : ");
		do {
			n = in.nextInt();
		}while(!IsNValid(n)); //Agar n >= 2
		float h = (end-start) / n;
		float[][] M = IsiMIntvl(start,end,n); //Matriks diisi oleh nilai
		System.out.println("Masukkan solusi persamaan : ");
		System.out.println("1. Masukkan solusi fx secara manual");
		System.out.println("2. Masukkan fx ke dalam fx yang ada");
		System.out.print("Tipe pemasukkan yang diinginkan : ");
		do {
		i = in.nextInt();
	} while (!IsInValid(i,1,2));
		if(i == 1) {
			float[] T = BacaArray(n + 1); //Mengisi array solusi
			float[][] Mori = new float[n + 1][n + 1];
			float[] Tori = new float[n + 1];
			CopyMatrix(M, Mori, T, Tori);
			SolusiI(M,T,Mori,Tori,X);
		} else {
			float[] T = new float[M.length];
			float[][] Mori = new float[n + 1][n + 1];
			float[] Tori = new float[n + 1];
			for(i = 0; i < M.length; i++) {
				T[i] = xfungsi(start);
				start = start + h;
			}
			CopyMatrix(M, Mori, T, Tori);
			SolusiI(M,T,Mori,Tori,X);
		}
	}

	//Menu dengan isi interpolasi untuk input x secara manual
  public static void ManualInterpolate() {
		float X = 0; int n;
		System.out.print("Masukkan jumlah nilai x yang ada : ");
		do {
			n = in.nextInt();
		}while(!IsNValid(n)); //Agar n >= 2

		float[] TI = new float[n]; //Array nilai x yang ada
		System.out.print("Masukkan nilai x : ");
		for(int i = 0; i < TI.length; i++) {
			TI[i] = in.nextFloat();
		}
		float[][] M = IsiMInterpolasi(TI,n); //Matriks diisi oleh nilai
		float[][] Mori = new float[n][n];
		System.out.print("Masukkan solusi persamaan f(x) : ");
		float[] T = BacaArray(n); //Membaca solusi ke array
		float[] Tori = new float[n];
		CopyMatrix(M,Mori,T,Tori);

		SolusiI(M,T,Mori,Tori,X);
	}

	//Membaca file external dan memasukkannya ke matrix untuk interpolasi
  public static float[][] BacaFileI() {
		System.out.print("Masukkan nama file untuk kalkulasi interpolasi : ");
		in.nextLine();
		String filename = in.nextLine();
		String str;
		int row, col, i, j;
		float[][] failed = new float[0][0];
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			try {
				List<String> list = new ArrayList<String>();
				while((str = br.readLine()) != null) {
					if (!str.equals("")) {
						list.add(str);
					}
				}

				String[] temparr = list.toArray(new String[0]);
				String[] parts = temparr[0].split(" ");
				row = temparr.length;
				col = parts.length;

				float[][] M = new float[row][col];

				for (i = 0; i < row; i++){
		      parts = temparr[i].split(" ");
		      for (j = 0; j < col; j++) {
		          M[i][j] = Float.parseFloat(parts[j]);
		      }
		    }
				return M;
		  } catch (IOException e) { return failed;}
		} catch (FileNotFoundException ex) {
        System.out.println("File tidak ditemukan!");
        return BacaFileI() ;
    }
	}

	//Memasukkan baris pertama matrix ke dalam array nilai x
	public static float[] XFromFile(float[][] M) {
		float[] T = new float[M[0].length];
		for(int i = 0; i < M[0].length; i++) {
			T[i] = M[0][i];
		}
		return T;
	}

	//Memasukkan baris kedua matrix ke dalam array nilai f(x)
	public static float[] FxFromFile(float[][] M) {
		float X = 0;
		float[] T = new float[M[0].length];
		for(int i = 0; i < M[0].length; i++) {
			T[i] = M[1][i];
		}
		return T;
	}

	//Menu untuk interpolasi dari file external
  public static void ExtInterpolate() {
		float X = 0;
		float[][] Mex = BacaFileI();
		float[] Tx = XFromFile(Mex);
		float[] TFx = FxFromFile(Mex);

		float[][] M = IsiMInterpolasi(Tx,Tx.length); //Matriks diisi oleh nilai
		float[][] Mori = new float[M.length][M[0].length];
		float[] Tori = new float[TFx.length];
		CopyMatrix(M,Mori,TFx,Tori);

		SolusiI(M,TFx,Mori,Tori,X);
	}

	//Menu yang berisi semua jenis kalkulasi interpolasi
  public static void InterpolateMenu() {
		int i;
		System.out.println("Tipe pembacaan data : ");
		System.out.println("1. Pembacaan data dari pengguna");
		System.out.println("2. Pembacaan data dari file external");
		System.out.print("Tipe pembacaan data yang ada inginkan : ");
		do {
			i = in.nextInt();
		}while(!IsInValid(i,1,2));
		System.out.println();
		if (i == 1) {
			System.out.println("Tipe Interpolasi : "); //Membuat menu untuk pilihan
			System.out.println("1. Memasukkan data selang [a..b] dengan n");
			System.out.println("2. Memasukkan nilai x secara manual");
			System.out.print("Tipe interpolasi yang ingin Anda lakukan : ");
			do {
				i = in.nextInt();
			}while(!IsInValid(i,1,2)); //Agar input diantara 1 dan 2
			System.out.println();
			if(i == 1) {
				GapInterpolate(); //Menjalankan penghampiran fungsi
			}
			else {
				ManualInterpolate(); //Menjalankan pencarian x dari fungsi
			}
		} else {
			ExtInterpolate();
		}
	}

	//Menu yang berisi interpolate menu dan gauss menu
  public static void MainMenu() {
    boolean repeat = true;
    while(repeat) {
      int i;
			System.out.println("  _____       _                 ____    ____	");
			System.out.println(" |           / \\     |      |  |    |  |    | ");
			System.out.println(" |   __     /   \\    |      |  |____   |____	");
			System.out.println(" |     }   /-----\\   |      |       |       | ");
			System.out.println(" |_____|  /       \\  |______|  |____|  |____| ");
			System.out.println("==============================================");
      System.out.println("      _         _____    _____   ____   ____    _____               _       ____    ");
      System.out.println("  |  |  \\   |     |     |       |    | |    |  |     |  |          / \\     |    |  | ");
      System.out.println("  |  |   \\  |     |     |_____  |____| |____|  |     |  |         /   \\    |____   |");
      System.out.println("  |  |    \\ |     |     |       | \\    |       |     |  |        /-----\\        |  |");
      System.out.println("  |  |     \\|     |     |_____  |  \\   |       |_____|  |_____  /       \\  |____|  | ");
      System.out.println(" =====================================================================================");
      System.out.println("Selamat datang!");
      System.out.println("Materi pada program ini :");
      System.out.println("1. Penyederhanaan SPL dengan metode Gauss");
      System.out.println("2. Interpolasi");
      System.out.print("Materi yang ingin di akses : ");
      do {
        i = in.nextInt();
      }while(!IsInValid(i,1,2));
      System.out.println();
      if(i == 1) {
        GaussMenu();
      }
      else {
        InterpolateMenu();
      }

      System.out.println("Pilihan untuk pengulangan program :");
      System.out.println("1. Ya, ulangi program");
      System.out.println("2. Tidak, sudahi program");
      System.out.print("Pilihan yang Anda inginkan : ");
      do {
        i = in.nextInt();
      }while(!IsInValid(i,1,2));

      if(i == 2) {
        repeat = false;
      }
    }
  }

	//Main program yang berisikan main menu
  public static void main(String[] args) {
    MainMenu();
  }

}
