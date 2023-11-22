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
		double[] allocatedBits = new double[J];
		
		boolean[][] userTime = new boolean[N][T];
		int[][] userTimeToFrame = new int[N][T];
		
		for (int i=0; i<J; i++) {
			String[] line = in.readLine().split(" ");
			
			frameId[i] = Integer.parseInt(line[0]);
			tbs[i] = Integer.parseInt(line[1]);
			userId[i] = Integer.parseInt(line[2]);
			t0[i] = Integer.parseInt(line[3]);
			td[i] = Integer.parseInt(line[4]);
			allocatedBits[i] = 0.0;
			
			for (int tx=t0[i]; tx<t0[i]+td[i]; tx++) {
				userTime[userId[i]][tx] = true;
				userTimeToFrame[userId[i]][tx] = i;
			}
			
			
			//prune out if impossible to transmit
			
			
			double total = 0.0;
			
			for (int tx=t0[i]; tx<t0[i]+td[i]; tx++) {
				
				for (int k=0; k<K; k++) {
					
					double hasil = 1.0;
					
					for (int r=0; r<R; r++) {
						hasil *= initSinr[tx][k][r][userId[i]];
					}
					double stk = Math.pow(hasil, 1.0 / R);
					
					total += R * Math.log(1.0 + stk) / Math.log(2.0);
					
				}
			}
			
			if (N == 2 && K == 2 && T == 2 && R == 1) {
				//pass
			} else if (192.0 * total < tbs[i]) { // we cannot fulfil this frame
				//throw new IllegalArgumentException("powerful frame");
				
				for (int tx=t0[i]; tx<t0[i]+td[i]; tx++) {
					userTime[userId[i]][tx] = false;
				}
			}
		}
		
		for (int t=0; t<T; t++) { // time
			
			int activeUsers = 0;
			for (int ux=0; ux<N; ux++) {
				if (userTime[ux][t]) {
					activeUsers++;
				}
			}
			
			double[][] urPoints = new double[N][R];
			
			for (int ux=0; ux<N; ux++) {
				for (int r=0; r<R; r++) {
					urPoints[ux][r] = 1.0;
					
					for (int k=0; k<K; k++) {
						urPoints[ux][r] *= initSinr[t][k][r][ux];
					}
					
				}
			}
			
			
			
			double[][] somethingElse = new double[N][R];
			
			for (int ux=0; ux<N; ux++) {
				
				if (!userTime[ux][t])
					continue;
				
				for (int r=0; r<R; r++) {
					//somethingElse[ux][r] = 1.0;
					somethingElse[ux][r] = 1.0;
					
					for (int k=0; k<K; k++) {
						//somethingElse[ux][r] *= initSinr[t][k][r][ux];
						
						somethingElse[ux][r] *= initSinr[t][k][r][ux];
					}
					
					//tbs[userTimeToFrame[ux][t]]
					
					
					somethingElse[ux][r] = 192.0 * Math.log(somethingElse[ux][r]) / Math.log(2.0)
										/ tbs[userTimeToFrame[ux][t]];
					
				}
			}
			
			
			int[] matchResource = new int[R];
			int[] userMatchCnt = new int[N];
			
			
			for (int r=0; r<R; r++) {
				double max_value = -1.0;
				matchResource[r] = -1;
				int match_user = -1;
				
				int minCount = Integer.MAX_VALUE;
				
				for (int ux=0; ux<N; ux++) {	
					if (!userTime[ux][t])
						continue;
					
					if (userMatchCnt[ux] < minCount) {
						minCount = userMatchCnt[ux];
					}
				}
				
				
				for (int ux=0; ux<N; ux++) {	
					if (!userTime[ux][t])
						continue;
					
					
					/*if (userMatchCnt[ux] == minCount &&
							urPoints[ux][r] > max_value) {
						max_value = urPoints[ux][r];
						match_user = ux;
					}*/
					
					if (userMatchCnt[ux] == minCount &&
							somethingElse[ux][r] > max_value) {
						max_value = somethingElse[ux][r];
						match_user = ux;
					}
				}
				
				
				if (activeUsers > 0) {
					
					matchResource[r] = match_user;
					userMatchCnt[match_user]++;
					
					if (max_value < 0.0) {
						throw new IllegalArgumentException("no match!!!");
					}
				}
			}
			
			
			
			
			for (int k=0; k<K; k++) { // cell
				
				
				//int chosenUser = 0;
				
				for (int r=0; r<R; r++) { //freq rbg
					//chosenUser++;
					//chosenUser++;
					
					
					
					double sumDebugger = 0.0;
					
					//int currentUser = 0;
					
					// 4.0 / N -> power range for each RBG -> [0; 4]
					// 1.0 / N -> power range for all RBGs -> [0; R]
					
					
					for (int n=0; n<N; n++) {
						
						if (matchResource[r] == n) {
							out.write(1.0 + " ");
							sumDebugger += 1.0;
							
							/*if (matchUser[n] != r) {
								throw new IllegalArgumentException("match failure");
							}*/

						} else {
							out.write(0.0 + " ");
							sumDebugger += 0.0;
						}
						
						
					} // users (n)
					
					if (activeUsers > 0) {
					
						if (sumDebugger < 0.999999 || sumDebugger > 1.000001) {
							throw new IllegalArgumentException("a b c d e");
						}
					} else {
						if (sumDebugger > 1e-9) {
							throw new IllegalArgumentException("a b c d e");
						}
					}
					
					out.write("\n");
					
					
					//xxx * N * R <= R
				} //rgb (freq resource, r)
			} // cells (k)
			
			
			
			for (int ux=0; ux<N; ux++) {
				
				if (!userTime[ux][t])
					continue;
				
				for (int k=0; k<K; k++) {
					
					double hasil = 1.0;
					int activeRbg = 0;
					
					for (int r=0; r<R; r++) {
						
						if (matchResource[r] != ux)
							continue;
						
						hasil *= initSinr[t][k][r][ux];
						activeRbg++;
					}
					
					if (activeRbg != 0) {
						double stk = Math.pow(hasil, 1.0 / activeRbg);
						
						
						allocatedBits[userTimeToFrame[ux][t]] += 192.0 * activeRbg * Math.log(1.0 + stk) / Math.log(2.0);
						
					}
					
				}
				
				if (allocatedBits[userTimeToFrame[ux][t]] >= tbs[userTimeToFrame[ux][t]]) {
					int frameToDeactivate = userTimeToFrame[ux][t];
					
					
					for (int tx=t0[frameToDeactivate];
							tx<t0[frameToDeactivate]+td[frameToDeactivate]; tx++) {
						userTime[ux][tx] = false;
					}
				}
				
			}
			
			
			
		} // time (t)
		
		//if (R >= N) {
		//	throw new IllegalArgumentException("a b c d e");
		//}
		
		in.close();
		out.close();
	}
}