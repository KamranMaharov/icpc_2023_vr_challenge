import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;


import java.util.Comparator;
import java.util.HashMap;

public class Main {
	
	private static class MaskGetter {
		int count;
		List<Integer> masks;
		int maskIdx;
		
		void recurse(int object, int number, int mask) {
			if (number == 7)
				return;
			
			mask ^= (1 << object);
			masks.add(mask);
			
			for (int nextObject = object + 1; nextObject < count; nextObject++) {
				recurse(nextObject, number + 1, mask);
			}
		}
		
		MaskGetter(int count) {
			this.count = count;
			this.masks = new ArrayList<>();
			
			for (int start = 0; start < this.count; start++) {
				recurse(start, 1, 0);
			}
			
			this.maskIdx = 0;
		}
		
		int next() {
			if (maskIdx < masks.size()) {
				return masks.get(maskIdx++);
			}
			return -1;
		}
	}
	
	
	private static class ResourceCell implements Comparable<ResourceCell> {
		int r;
		int k;
		double signal;
		
		ResourceCell(
				int r,
				int k,
				double signal) {
			this.r = r;
			this.k = k;
			this.signal = signal;
		}

		@Override
		public int compareTo(Main.ResourceCell other) {
			if (this.signal == other.signal)
				return 0;
			else if (this.signal > other.signal)
				return -1;
			return 1;
		}
		
		
	}
	
	
	public static void main(String[] args) throws Exception {
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in), 200000 * 20);
		BufferedWriter out = new BufferedWriter(new OutputStreamWriter(System.out), 50000000);
		
		
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
		
		
		int initDiscarded = 0;
		
		boolean oneTimed = false;
		
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
			
			if (td[i] == 1)
				oneTimed = true;
			
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
				initDiscarded++;
				
				for (int tx=t0[i]; tx<t0[i]+td[i]; tx++) {
					//userTime[userId[i]][tx] = false;
				}
			}
		}
		
		//following 2 var's are in sync
		int fulfilledV1 = 0;
		int fulfilledV2 = 0;
		
		int totalShareables = 0;
		int totalBestCount = 0;
		int totalExtraAdded = 0;
		
		int totalAllocs = 0;
		
		for (int t=0; t<T; t++) {
			//out.write("TIME " + t + "\n");
			
			

			double[] cellPower = new double[K];
			double[][] cellResourcePower = new double[K][R];
			double[][][] userCellResourcePower = new double[N][K][R];
			
			Set<Integer> bestResources = new HashSet<>();
			
			
			
			
			for (int iter = 1; iter <= 3; iter++) {
			
			
			List<Integer>[][] shared = new List[R][K];
			
			for (int r=0; r<R; r++) {
				for (int k=0; k<K; k++) {
					shared[r][k] = new ArrayList<>();
				}
			}
			
			for (int ux=0; ux<N; ux++) {
				if (!userTime[ux][t]) {
					continue;
				}
				
				for (int r=0; r<R; r++) {
					
					
					if (bestResources.contains(r)) {
						continue;
					}
					
					
					for (int k=0; k<K; k++) {
						
						double canBeGiven = Math.min(4.0 - cellResourcePower[k][r], R - cellPower[k]);
						
						if (tbs[userTimeToFrame[ux][t]] < allocatedBits[userTimeToFrame[ux][t]])
							throw new IllegalArgumentException("get out!!!");
						if (t0[userTimeToFrame[ux][t]] + td[userTimeToFrame[ux][t]] <= t)
							throw new IllegalArgumentException("get out!!!");
						
						if ( 1.0 * (tbs[userTimeToFrame[ux][t]] - allocatedBits[userTimeToFrame[ux][t]]) 
								<= 192 * Math.log(1.0 + canBeGiven * initSinr[t][k][r][ux]) / Math.log(2.0) ) {
							//out.write("time " + t + " STRONG " + ux + " " + r + " " + k + " " +
							//		tbs[userTimeToFrame[ux][t]] + " " + 192 * Math.log(1.0 + initSinr[t][k][r][ux]) / Math.log(2.0) +	"\n");
							shared[r][k].add(ux);
						}
					}
				}
			}
			
			
			int bestUserCount = 1;
			Map<Integer, Double> bestPowerDist = new HashMap<>();
			int bestResource = -1;
			int bestCell = -1;
			
			for (int r=0; r<R; r++) {
				
				for (int k=0; k<K; k++) {
					if (shared[r][k].size() == 0)
						continue;
					
					//out.write("SH " + r + ", " + k + " ::::::  ");
					//for (Integer ux : shared[r][k]) {
					//	out.write(ux + " ");
					//}
					//out.write("\n");
					
					
					
					

					int num = shared[r][k].size();
					
					MaskGetter mg = new MaskGetter(shared[r][k].size());
					

					for (int mask=mg.next(); mask != -1; mask = mg.next()) {
						
						if (Integer.bitCount(mask) <= bestUserCount ||
								Integer.bitCount(mask) > 6) { // best user found has more users
							continue;
						}
						
						List<Integer> subUsers = new ArrayList<>();
						for (int idx=0; idx<num; idx++) {
							if ((mask & (1<<idx)) == 0)
								continue;
							
							subUsers.add(shared[r][k].get(idx));
						}
						
						if (Integer.bitCount(mask) != subUsers.size())
							throw new IllegalArgumentException("no way!!!");
						
						double totalPower = 0.0;
						Map<Integer, Double> powerDist = new HashMap<>();
						
						for (int idx1=0; idx1<subUsers.size(); idx1++) {
							
							
							double interference = 1.0;
							
							
							for (int idx2=0; idx2<subUsers.size(); idx2++) {
								if (idx1 == idx2)
									continue;
								
								interference *= Math.pow(Math.E,
										dfactor[k][r][subUsers.get(idx1)][subUsers.get(idx2)]);
									
							}
							
							int user1 = subUsers.get(idx1);
							
							if (tbs[userTimeToFrame[user1][t]] < allocatedBits[userTimeToFrame[user1][t]])
								throw new IllegalArgumentException("get out!!!");
							if (t0[userTimeToFrame[user1][t]] + td[userTimeToFrame[user1][t]] <= t)
								throw new IllegalArgumentException("get out!!!");
							
							double power = ( Math.pow(2.0, 1.0 * (tbs[userTimeToFrame[user1][t]] - allocatedBits[userTimeToFrame[user1][t]])
									 / 192.0)  - 1.0 ) /
												(initSinr[t][k][r][user1] * interference);
							if (power < 0.0) {
								throw new IllegalArgumentException("fdsaff");
								//out.write("POWER " + " " + power + "\n");
							}
							totalPower += power;
							
							if (totalPower > Math.min(4.0 - cellResourcePower[k][r], R - cellPower[k])) {
								break;
							}
							
							
							powerDist.put(user1, power);
							//out.write("POWER " + power + "\n");
						}
						
						if (totalPower <= Math.min(4.0 - cellResourcePower[k][r], R - cellPower[k])) {
							bestUserCount = subUsers.size();
							bestPowerDist = powerDist;
							bestResource = r;
							bestCell = k;
						}
						
						//out.write(null);
						if (totalPower <= 4.0) {
							//out.write("    TOTAL " + totalPower + " " + mask + " :::: " + "  , USERS: ");
							//for (final Integer user : subUsers) {
							//	out.write(user + " ");
							//}
							//out.write("\n");
						}
					}

				}
			}
			
			
			
			
			
			
			
			if (bestCell >= 0) {
				double alreadyPowered = 0.0;
				for (final Double pw : bestPowerDist.values()) {
					alreadyPowered += pw;
				}
				
				cellPower[bestCell] += alreadyPowered;
				cellResourcePower[bestCell][bestResource] += alreadyPowered;
				bestResources.add(bestResource);
				
				//out.write("TIME " + t + " ::: USERS :: ");
				
				for (final Integer ux : bestPowerDist.keySet()) {
					//out.write(ux + " ");
					userCellResourcePower[ux][bestCell][bestResource] = bestPowerDist.get(ux);
					
					int frame = userTimeToFrame[ux][t];
					
					if (1.0 * (tbs[frame] - allocatedBits[frame])  <= 0.0) {
						throw new IllegalArgumentException("fdafafa");
					}
					
					allocatedBits[frame] += 1.0 * (tbs[frame] - allocatedBits[frame]);
					
					
					if (allocatedBits[frame] >= tbs[frame]) {
						//out.write("full_allocation " + ux + " during subset force" + "\n");
						
						
						for (int tt=t0[frame];
								tt<t0[frame] + td[frame];
								tt++) {
							userTime[ux][tt] = false;
						}
						
						
						totalAllocs++;
						
					}
					
				}
				
				//out.write("\n");
				//out.write("resource _ cell : " + bestResource + " " + bestCell + "\n");
				
			}
			
			
			
			}
			
			int reman = 0;
			for (int r=0; r<R; r++) {	
				if (bestResources.contains(r))
					continue;
				reman++;
			}
			
			boolean[] reserved = new boolean[R];
			
			for (int r=0; r<R; r++) {
				
				if (bestResources.contains(r))
					continue;
				
				boolean breakable = false;
				
				for (int ux=0; ux<N; ux++) {
					
					if (!userTime[ux][t]) {
						continue;
					}
					
					for (double maxPower = 0.05; maxPower <= 4.0; maxPower += 0.05) {
					
					double total = 0.0;
					
					for (int k=0; k<K; k++) {
						
						// give power in bestCell
						
						if (cellPower[k] > R || cellResourcePower[k][r] > 4.0) {
							throw new IllegalArgumentException("no way!!!");
						}
						
						double power = Math.min( Math.min(maxPower, 1.0 * (R - cellPower[k]) / reman),
								Math.min(R - cellPower[k], 4 - cellResourcePower[k][r]) );
						
						
						total += Math.log(1.0 + power * initSinr[t][k][r][ux]) / Math.log(2);
						
					}
					
					if (allocatedBits[userTimeToFrame[ux][t]] + total * 192.0 >=
							tbs[userTimeToFrame[ux][t]]) {
						
						allocatedBits[userTimeToFrame[ux][t]] += total * 192.0;
						
						for (int k=0; k<K; k++) {
							
							// give power in bestCell
							
							if (cellPower[k] > R || cellResourcePower[k][r] > 4.0) {
								throw new IllegalArgumentException("no way!!!");
							}
							
							double power = Math.min(Math.min(maxPower, 1.0 * (R - cellPower[k]) / reman),
									Math.min(R - cellPower[k], 4 - cellResourcePower[k][r]) );
							
							userCellResourcePower[ux][k][r] = power;
							cellPower[k] += power;
							cellResourcePower[k][r] += power;	
						}
						
						reserved[r] = true;
						
						for (int tt=t0[userTimeToFrame[ux][t]];
								tt<t0[userTimeToFrame[ux][t]] + td[userTimeToFrame[ux][t]];
								tt++) {
							userTime[ux][tt] = false;
						}
						
						
						totalAllocs++;
						breakable = true;
						break;
						
					}
					
					
					}
					
					if (breakable)
						break;
					
				}
				reman--;
			}
			
			
			/*
			for (int ux=0; ux<N; ux++) {
				if (!userTime[ux][t]) {
					continue;
				}
				
				
				for (int r=0; r<R; r++) {
					if (bestResources.contains(r) || reserved[r])
						continue;
					
					List<ResourceCell> rcs = new ArrayList<>();
					
					for (int k=0; k<K; k++) {
						rcs.add(new ResourceCell(r, k, initSinr[t][k][r][ux]));
					}
					
					Collections.sort(rcs);
					
					double total = 0.0;
					
					for (int ki=0; ki<K; ki++) {
						int k = rcs.get(ki).k;
						
						double powerLeft = Math.min(4.0,
								Math.min(R - cellPower[k], 4.0 - cellResourcePower[k][r]));
						
						total += Math.log(1.0 + initSinr[t][k][r][ux] * powerLeft) / Math.log(2.0);
					}
					
					
					
					
					if (allocatedBits[userTimeToFrame[ux][t]] + total * 192.0 >=
							tbs[userTimeToFrame[ux][t]]) {
						
						allocatedBits[userTimeToFrame[ux][t]] += total * 192.0;
						
						for (int ki=0; ki<K; ki++) {
							int k = rcs.get(ki).k;
							// give power in bestCell
							
							if (cellPower[k] > R || cellResourcePower[k][r] > 4.0) {
								throw new IllegalArgumentException("no way!!!");
							}
							
							double power = Math.min(4.0,
									Math.min(R - cellPower[k], 4 - cellResourcePower[k][r]) );
							
							userCellResourcePower[ux][k][r] = power;
							cellPower[k] += power;
							cellResourcePower[k][r] += power;	
						}
						
						reserved[r] = true;
						
						for (int tt=t0[userTimeToFrame[ux][t]];
								tt<t0[userTimeToFrame[ux][t]] + td[userTimeToFrame[ux][t]];
								tt++) {
							userTime[ux][tt] = false;
						}
						
						
						totalAllocs++;
						break;
						
					}
					
					
					
					
					
					
				}
				
			}*/
			
			
			
			
			
			
			double[] currentAlloc = new double[N];
			boolean[][] choice = new boolean[R][N];
			
			int resRemaining = 0;
			
			for (int r=0; r<R; r++) {	
				if (bestResources.contains(r) || reserved[r])
					continue;
				resRemaining++;
			}
			
			for (int r=0; r<R; r++) {
				
				if (bestResources.contains(r) || reserved[r])
					continue;
				
				double maxdiff = 0.0;
				double okpower = 0.0;
				double nextAlloc = -1.0;
				int user = -1;
				
				for (int ux=0; ux<N; ux++) {
					
					if (!userTime[ux][t]) {
						continue;
					}
					
					
					for (double maxPower = 0.1; maxPower <= 4.0; maxPower += 0.1) {
					
					choice[r][ux] = true;
					
					double total = 0.0;
					
					for (int k=0; k<K; k++) {
						
						int bCount = 0;
						
						double hasil = 1.0;
						
						for (int ri=0; ri<R; ri++) {
							if (choice[ri][ux]) {
								bCount++;
								
								double power = -1;
								if (ri == r) {
									power = Math.min(
												Math.min(1.0 * (R-cellPower[k]) / resRemaining, maxPower),
												Math.min(R - cellPower[k], 4.0 - cellResourcePower[k][ri]));
								} else {
									power = userCellResourcePower[ux][k][ri];
								}
								
								hasil *= power * initSinr[t][k][ri][ux];
							}
						}
						
						double snt = Math.pow(hasil, 1.0 / bCount);
						total += bCount * Math.log(1.0 + snt) / Math.log(2.0);
					}
					
					double diff = Math.min(total * 192.0 - currentAlloc[ux], tbs[userTimeToFrame[ux][t]] - allocatedBits[userTimeToFrame[ux][t]]) /
							(tbs[userTimeToFrame[ux][t]] - allocatedBits[userTimeToFrame[ux][t]]);
					
					if (diff > maxdiff) {
						maxdiff = diff ;
						okpower = maxPower;
						nextAlloc = total * 192.0;
						user = ux;
					}
					
					choice[r][ux] = false;
					
					
					
					}
					
					
					
				}
				
				if (user >= 0) {
					choice[r][user] = true;
					allocatedBits[userTimeToFrame[user][t]] += (nextAlloc - currentAlloc[user]);
					currentAlloc[user] = nextAlloc;
					
					for (int k=0; k<K; k++) {
						double power = Math.min(
								Math.min(1.0 * (R-cellPower[k]) / resRemaining, okpower),
								Math.min(R - cellPower[k], 4.0 - cellResourcePower[k][r]));
						cellPower[k] += power;
						cellResourcePower[k][r] += power;
						userCellResourcePower[user][k][r] = power;
					}
					
					//out.write("allocating resource " + r + " to user " + user + "\n");
					
					if (allocatedBits[userTimeToFrame[user][t]] >= tbs[userTimeToFrame[user][t]]) {
						//out.write("full_allocation " + user + "\n");
						
						
						for (int tt=t0[userTimeToFrame[user][t]];
								tt<t0[userTimeToFrame[user][t]] + td[userTimeToFrame[user][t]];
								tt++) {
							userTime[user][tt] = false;
						}
						
						
						totalAllocs++;
						
					}
					
				}
				resRemaining--;
				
			}
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			

			for (int k=0; k<K; k++) {	
				
				for (int r=0; r<R; r++) {
					
					
					
					for (int ux=0; ux<N; ux++) {
					
					
						out.write(userCellResourcePower[ux][k][r] + " ");
						
						
					}
					
					out.write("\n");
				}
				
			}
			
			
			
		}
		
		
		//out.write(totalAllocs + "\n");
		
		//out.write(totalxxBestCount + "\n");
		//out.write(totalShareables + "\n");
		//out.write(totalExtraAdded + "\n");
		//if (fulfilledV1 != fulfilledV2)
	//		throw new IllegalArgumentException("inconsistency");
		
	//	if (fulfilledV1 < J * 0.4 && false) {
	//		throw new IllegalArgumentException("too low signals");
	//	}
		
		//out.write(J + " " + fulfilledV1 + " " + fulfilledV2 + "\n");
		//out.write(J + " " + initDixscarded + " " + fulfilledV1 + " " + discarded + "\n");
		
		in.close();
		out.close();
	}
}