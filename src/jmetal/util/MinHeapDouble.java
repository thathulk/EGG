package jmetal.util;

public class MinHeapDouble {
	
	private int heapSize;
	public MinHeapDouble() {}
	
	public int getHeapSize() {
		return heapSize;
	}
	public void setHeapSize(int heapSize) {
		this.heapSize = heapSize;
	}
	
	public void build(double heap[], int heapSize) {
		this.setHeapSize(heapSize);
		for(int i=heapSize/2; i>=1; i--) {
			heapIFY(heap, i);
		}
	}
	
	public double getAndDelete(double heap[], int k) {
		double dis = heap[k];
		heap[k] = heap[heapSize];
		this.setHeapSize(heapSize-1);
		heapIFY(heap, k);
		return dis;
	}
	
	public void heapIFY(double heap[], int k) {
		int l = k * 2;
		int r = l + 1;
		int min = 0;
		if(l <= heapSize && heap[l] < heap[k]) {
			min = l;
		} else {
			min = k;
		}
		if(r <= heapSize && heap[r] < heap[min]) {
			min = r;
		} 
		if(min != k) {
			double temp = heap[k];
			heap[k] = heap[min];
			heap[min] = temp;
			
			heapIFY(heap, min);
		}
	}
}
