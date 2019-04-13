package jmetal.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import jmetal.core.SolutionSet;
import jmetal.util.Distance;

public class ArchiveTruncation {

	public ArchiveTruncation() {}
	
	/* 通过标记某个点是否被删除来降低时间复杂度	 */
	/*
	public SolutionSet ArchiveTruncation1(SolutionSet archive, int size) {
		
		int archiveSize = archive.size();
		if(archiveSize <= size) {
			return archive;
		}
		
		
		double dis[][] = new double[archiveSize][archiveSize];
		List<Edge> edgeList = new LinkedList<Edge>();
		Distance distance = new Distance();
		boolean isRemove[] = new boolean[archiveSize];
		
		for(int i=0; i<archiveSize; i++) {
			for(int j=0; j<i; j++) {
				dis[i][j] = distance.distanceBetweenObjectives(archive.get(i), archive.get(j));
				dis[j][i] = dis[i][j];
				edgeList.add(new Edge(i, j, dis[i][j]));
			}
			isRemove[i] = false;
		}
		Collections.sort(edgeList, new EdgeComparator());

		Edge edge;
		int pa, pb, toRemove = 0;
		double length, min, max;
		List<Double> a = new ArrayList<Double>();
		List<Double> b = new ArrayList<Double>();
		List<Integer> remove = new ArrayList<Integer>();
		
		for(int k=size; k<archiveSize; k++) {
			edge = edgeList.get(0);
			pa = edge.getPa(); 
			pb = edge.getPb();
			while(isRemove[pa] == true || isRemove[pb] == true) {
				edgeList.remove(0);
				edge = edgeList.get(0);
				pa = edge.getPa(); 
				pb = edge.getPb();
			}
			length = edge.getLength();
			a.clear();
			b.clear();
			for(int i=0; i<archiveSize; i++) {
				if(i != pa && i != pb) {
					a.add(dis[pa][i]);
					b.add(dis[pb][i]);
				}
			}
			Collections.sort(a);
			Collections.sort(b);
			
			for(int i=0; i<a.size(); i++) {
				if(a.get(i) != b.get(i)) {
					if(a.get(i) < b.get(i)) {
						toRemove = pa;
					} else {
						toRemove = pb;
					}
					break;
				}
			}
			isRemove[toRemove] = true;
			remove.add(toRemove);
			for(int i=0; i<archiveSize; i++) {
				dis[i][toRemove] = Double.MAX_VALUE;
			}
		}
		Collections.sort(remove);
		while(remove.size() > 0) {
			int i = remove.get(remove.size()-1);
			archive.remove(i);
			remove.remove(remove.size()-1);
		}
		
		return archive;
	}
 */
		
	/*根据所有点到该点距离平方的倒数和来判断拥挤程度 */
	/*
	public SolutionSet ArchiveTruncation2(SolutionSet archive, int size) {
		int archiveSize = archive.size();
		if(archiveSize <= size) {
			return archive;
		}
		
		double max = 0, min = Double.MAX_VALUE;
		double dis[][] = new double[archive.size()][archive.size()];
		double L[][] = new double[archive.size()][archive.size()];
		Distance distance = new Distance();
		
		for(int i=0; i<archiveSize; i++) {
			for(int j=i+1; j<archiveSize; j++) {
				dis[i][j] = distance.distanceBetweenObjectives(archive.get(i), archive.get(j));
				dis[j][i] = dis[i][j];
				if(dis[i][j] > max) {
					max = dis[i][j];
				}
				if(dis[i][j] < min) {
					min = dis[i][j];
				}
			}
		}
		
		int TL = archiveSize * archiveSize;
		List<Node> density = new LinkedList<Node>();
		double d = 0;
		for(int i=0; i<archiveSize; i++) {
			d = 0;
			for(int j=0; j<archiveSize; j++) {
				if(i == j) continue;
				if(dis[i][j] == min) {
					L[i][j] = L[j][i] = 1;
				} else {
					L[i][j] = L[j][i] = (int) Math.ceil(TL * (dis[i][j]-min) / (max-min));
				}
				d += (1.0 / Math.pow(L[i][j], 2));
			}
			density.add(new Node(i, d));
		}
		
		int toRemove, l=0, lr =0;
		while(archive.size() > size) {
			archiveSize = archive.size();
			max = 0;
			toRemove = 0;
			int k = 0;
			for(Node node : density) {
				if(node.getD() > max) {
					max = node.getD();
					toRemove = k;
					lr = node.getLocation();
				}
				k++;
			}
			for(int i=0; i<archiveSize; i++) {
				if(i != toRemove) {
					l = density.get(i).location;
					d = density.get(i).getD();
					density.get(i).setD(d-(1.0 / Math.pow(L[l][lr], 2)));
				}
			}
			archive.remove(toRemove);
			density.remove(toRemove);
		}
		
		return archive;
	}
	public class Node {
		int location;
		double d;
		public Node(int location, double d) {
			this.location = location;
			this.d = d;
		}
		public int getLocation() {
			return location;
		}
		public void setLocation(int location) {
			this.location = location;
		}
		public double getD() {
			return d;
		}
		public void setD(double d) {
			this.d = d;
		}
	}
*/
	
	/* 所有边组成一个堆 */
	public SolutionSet ArchiveTruncationByHeap(SolutionSet archive, int size) {
		int archiveSize = archive.size();
		if(archiveSize <= size) {
			return archive;
		}
		
		double dis[][] = new double[archiveSize][archiveSize];
		int location[][] = new int[archiveSize][archiveSize];
		Edge edges[] = new Edge[archiveSize*archiveSize];
		boolean isRemove[] = new boolean[archiveSize];
		Distance distance = new Distance();
		
		int k = 1;
		for(int i=0; i<archiveSize; i++) {
			for(int j=i+1; j<archiveSize; j++) {
				dis[i][j] = distance.distanceBetweenObjectives(archive.get(i), archive.get(j));
				dis[j][i] = dis[i][j];
				location[i][j] = location[j][i] = k;
				edges[k++] = new Edge(i, j, dis[i][j]);
			}
			isRemove[i] = false;
		}
		
		int heapSize = k - 1;
		MinHeap edgeHeap = new MinHeap();
		edgeHeap.bulid(edges, location, heapSize);
		
		double disToA[], disToB[];
		int heapASize, heapBSize;
		MinHeapDouble heapA = null, heapB = null;
		List<Integer> needToRemove = new ArrayList<Integer>();
		int pa, pb, toRemove = 0;
		Edge shortest = null;
		for(k=size; k<archiveSize; k++) {
			shortest = edgeHeap.getAndDelete(edges, location, 1);
			pa = shortest.getPa();
			pb = shortest.getPb();
			
			while(isRemove[pa] == true || isRemove[pb] == true) {
				shortest = edgeHeap.getAndDelete(edges, location, 1);
				pa = shortest.getPa();
				pb = shortest.getPb();
			}
			
			disToA = new double[archiveSize];
			disToB = new double[archiveSize];

			int n = 1;
//			double da = 0, db = 0;
			
			for(int i=0; i<archiveSize; i++) {
				if(i != pa && i != pb) {
					disToA[n] = dis[pa][i];
					disToB[n] = dis[pb][i];
					n ++;
//					da += 1.0 / Math.pow(dis[pa][i], 2);
//					db += 1.0 / Math.pow(dis[pb][i], 2);
				}
			}
//			if(da > db) {
//				toRemove = pa;
//			} else {
//				toRemove = pb;
//			}
			heapASize = heapBSize = n-1;
			heapA = new MinHeapDouble();
			heapB = new MinHeapDouble();
			heapA.build(disToA, heapASize);
			heapB.build(disToB, heapBSize);
			
			while(heapA.getHeapSize() > 0 && heapB.getHeapSize() > 0) {
				double sa = heapA.getAndDelete(disToA, 1);
				double sb = heapB.getAndDelete(disToB, 1);
				if(sa != sb) {
					if(sa < sb) {
						toRemove = pa;
					} else {
						toRemove = pb;
					}
					break;
				}
			}
			
			isRemove[toRemove] = true;
			needToRemove.add(toRemove);
			for(int i=0; i<archiveSize; i++) {
				dis[i][toRemove] = Double.MAX_VALUE;
			}
		}
		
		Collections.sort(needToRemove);
		while(needToRemove.size() > 0) {
			int i = needToRemove.get(needToRemove.size()-1);
			archive.remove(i);
			needToRemove.remove(needToRemove.size()-1);
		}
		
		return archive;
	}
	
	/* 每个点到其他点的距离组成一个堆 */
	public SolutionSet ArchiveTruncationByMultiHeaps(SolutionSet archive, int size) {
		int archiveSize = archive.size();
		if(archiveSize <= size) {
			return archive;
		}
		double dis[][] = new double[archiveSize][archiveSize];
		int location[][] = new int[archiveSize][archiveSize];
		Edge edges[][] = new Edge[archiveSize][archiveSize];
		int heapSize[] = new int[archiveSize];
		MinHeap h[] = new MinHeap[archiveSize];
		boolean isRemove[] = new boolean[archiveSize];
		Distance distance = new Distance();
		
		for(int i=0; i<archiveSize; i++) {
			h[i] = new MinHeap();
			int k = 1;
			for(int j=0; j<archiveSize; j++) {
				if(i == j) continue;
				dis[i][j] = distance.distanceBetweenObjectives(archive.get(i), archive.get(j));
				location[i][j] = k;
				edges[i][k++] = new Edge(i, j, dis[i][j]);
			}
			isRemove[i] = false;
			heapSize[i] = k-1;
			h[i].bulid(edges[i], location, heapSize[i]);
		}
		
		double disToA[], disToB[];
		int heapASize, heapBSize;
		MinHeapDouble heapA = null, heapB = null;
		List<Integer> needToRemove = new ArrayList<Integer>();
		int pa=0, pb=0, toRemove = 0;
		double la, lb, l, shortest;
		for(int k=size; k<archiveSize; k++) {
			shortest = Double.MAX_VALUE;
			for(int i=0; i<archiveSize; i++) {
				if(isRemove[i] == false && edges[i][1].getLength() < shortest) {
					shortest = edges[i][1].getLength();
					pa = edges[i][1].getPa();
					pb = edges[i][1].getPb();
				}
			}
			heapASize = h[pa].getHeapSize();
			heapBSize = h[pb].getHeapSize();
			disToA = new double[heapASize+1];
			disToB = new double[heapBSize+1];
			heapA = new MinHeapDouble();
			heapB = new MinHeapDouble();
			heapA.setHeapSize(heapASize);
			heapB.setHeapSize(heapBSize);
			for(int i=1; i<=heapASize; i++) {
				disToA[i] = edges[pa][i].getLength();
				disToB[i] = edges[pb][i].getLength();
			}
			while(heapA.getHeapSize() > 0 && heapB.getHeapSize() > 0) {
				double sa = heapA.getAndDelete(disToA, 1);
				double sb = heapB.getAndDelete(disToB, 1);
				if(sa != sb) {
					if(sa < sb) {
						toRemove = pa;
					} else {
						toRemove = pb;
					}
					break;
				}
			}
				
			isRemove[toRemove] = true;
			needToRemove.add(toRemove);
			for(int i=0; i<archiveSize; i++) {
				if(isRemove[i] == false) {
					h[i].getAndDelete(edges[i], location, location[i][toRemove]);
				}
			}
			
		}
		
		Collections.sort(needToRemove);
		while(needToRemove.size() > 0) {
			int i = needToRemove.get(needToRemove.size()-1);
			archive.remove(i);
			needToRemove.remove(needToRemove.size()-1);
		}
		return archive;
	}
	
	
}
