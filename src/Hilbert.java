import java.util.*;
import java.text.*;
import java.io.*;

public class Hilbert {
	private static DecimalFormat df = new DecimalFormat("#.###");
	private static Scanner in = new Scanner (System.in);
	
	public static void main(String[] args) {
		System.out.println("Masukan n Matriks Hilbert : ");
		int n = in.nextInt();
		float[][] M =  new float[n][n];
		float[] T = new float[n];

		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				M[i][j] = 1.0f / ((i+1.0f) + (j+1.0f) - 1.0f);
			}
			T[i] = 1.0f;
		}
		
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				System.out.print(df.format(M[i][j]) + "\t");
			}
			System.out.println("\n");
		}
		
		System.out.println("Masukan Nama File untuk disimpan : ");
		in.nextLine();
		String filename = in.nextLine() ;
		try{
			PrintWriter writer = new PrintWriter(filename, "UTF-8");
			for (int i=0; i<n; i++) {
				for (int j=0; j<n; j++) {
					writer.print((M[i][j]) + "\t");
				}
			writer.print("|\t" + T[i]);
			writer.println("\n");
			}
			writer.close();
			System.out.println("Pembuatan Matrix Hilbert Berhasil");
		} catch (IOException e) {
			System.out.println("Pembuatan Matrix Hilbert Gagal");
		}
	}
}
