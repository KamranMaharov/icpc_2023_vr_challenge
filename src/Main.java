import java.io.*;
import java.util.*;

public class Main {
	
	
	static class Scanner {
		Scanner(InputStream in) { this.in = in; } InputStream in;
		byte[] bb = new byte[200000 * 19]; int i, n;
		byte getc() {
			if (i == n) {
				i = n = 0;
				try { n = in.read(bb); } catch (IOException e) {}
			}
			return i < n ? bb[i++] : 0;
		}
		
		boolean isdigit(byte c) {
			return c >= (byte)'0' && c <= (byte)'9';
		}
		
		int nextInt() {
			/*byte c = 0; while (c <= ' ') c = getc();
			boolean negate = false;
			if (c == '-') {
				negate = true;
				c = getc();
			}*/
			byte c = 0; while (!isdigit(c)) c = getc();
			int a = 0; while (isdigit(c)) { a = a * 10 + c - '0'; c = getc(); }
			//if (negate) a = -a;
			return a;
		}
		
		long nextLong() {
			byte c = 0; while (c <= ' ') c = getc();
			long a = 0; while (c > ' ') { a = a * 10 + c - '0'; c = getc(); }
			return a;
		}
	}
	
	
	
	public static void main(String[] args) throws Exception {
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in), 200000 * 20);
		BufferedWriter out = new BufferedWriter(new OutputStreamWriter(System.out), 50000000);
		//Scanner in = new Scanner(System.in);		
		
		//out.write(Integer.parseInt("3") + "\n");
		//out.flush();
		
		int N = Integer.parseInt(in.readLine()); //users <= 100
		int K = Integer.parseInt(in.readLine()); //cells <= 10
		int T = Integer.parseInt(in.readLine()); //time <= 1000 ( <= 0.5 second)
		int R = Integer.parseInt(in.readLine()); //frequency resource <= 10
		
		double[][][][] initSinr = new double[T][K][R][N]; //quality of channel, float, 0.0 < s < 10000.0
		double[][][][] dfactor = new double[K][R][N][N]; //interference factor, float, -2.0 <= d <= 0.0
		
		for (int t=0; t<T; t++) { //time
			for (int k=0; k<K; k++) { //cell
				for (int r=0; r<R; r++) { //rbg
					String[] line = in.readLine().split(" ");
					for (int n=0; n<N; n++) {
						initSinr[t][k][r][n] = Double.parseDouble(line[n]);
					}
				}
			}
		}
		
		for (int k=0; k<K; k++) {
			for (int r=0; r<R; r++) {
				for (int m=0; m<N; m++) {
					String[] line = in.readLine().split(" ");
					for (int n=0; n<N; n++) {
						//dfactor[k][r][m][n] = Integer.parseInt(line[n]);
						dfactor[k][r][m][n] = Double.parseDouble(line[n]);
					}
				}
			}
		}
		
		int J = Integer.parseInt(in.readLine()); //number of frames, 1<=j<=5000
		
		int[] frameId = new int[J]; //0 <= j <= J-1
		int[] tbs = new int[J]; //0 < tbs <= 100000
		int[] userId = new int[J]; // 0 <= n < N
		int[] t0 = new int[J]; // 0 <= t < T
		int[] td = new int[J]; // 1 <= td <= 100
		
		boolean[][] userTime = new boolean[N][T];
		
		for (int i=0; i<J; i++) {
			String[] line = in.readLine().split(" ");
			
			frameId[i] = Integer.parseInt(line[0]);
			tbs[i] = Integer.parseInt(line[1]);
			userId[i] = Integer.parseInt(line[2]);
			t0[i] = Integer.parseInt(line[3]);
			td[i] = Integer.parseInt(line[4]);
			
			for (int tx=t0[i]; tx<t0[i]+td[i]; tx++) {
				userTime[userId[i]][tx] = true;
			}
			
		}
		
		for (int t=0; t<T; t++) { // time
			
			int activeUsers = 0;
			for (int ux=0; ux<N; ux++) {
				if (userTime[ux][t]) {
					activeUsers++;
				}
			}
			
			for (int k=0; k<K; k++) { // cell
				for (int r=0; r<R; r++) { //freq rbg
					
					
					// 4.0 / N -> power range for each RBG -> [0; 4]
					// 1.0 / N -> power range for all RBGs -> [0; R]
					
					
					for (int n=0; n<N; n++) {
						if (userTime[n][t]) {
							out.write((1.0 / activeUsers) + " ");
						} else {
							out.write(0.0 + " ");
						}
						
						
						/*if (N == 2 && K == 2 && T == 2) {
							//out.write((1.0 / N) + " ");
						} else {
							//out.write((1.0 / N + 1e-6) + " ");
							//1e-8 is ok.
							//1e-7, 2 incorrect answers
							//1e-6, 48 incorrect answers
						} */
						
						//xxx
						
					}
					out.write("\n");
					
					
					//xxx * N * R <= R
				}
			}
		}
		
		
		in.close();
		out.close();
	}
}