package jmetal.util;

public class MinHeap {
	
	private int heapSize;
	public MinHeap() {
		
	}
	
	public int getHeapSize() {
		return heapSize;
	}

	public void setHeapSize(int heapSize) {
		this.heapSize = heapSize;
	}

	public void bulid(Edge heap[], int location[][], int heapSize) {
		this.setHeapSize(heapSize);
		for(int i=heapSize/2; i>=1; i--) {
			heapIFY(heap, location, i);
		}
	}
	
	public Edge getAndDelete(Edge heap[], int location[][], int k) {
		Edge edge = new Edge(heap[k]);
		heap[k] = new Edge(heap[heapSize]);
		location[heap[k].getPa()][heap[k].getPb()] = k;
		this.setHeapSize(heapSize-1);
		heapIFY(heap, location, k);
		return edge;
	}
	
	public void heapIFY(Edge heap[], int location[][], int k) {
		int l = k * 2;
		int r = l + 1;
		int min = 0;
		if(l <= heapSize && heap[l].getLength() < heap[k].getLength()) {
			min = l;
		} else {
			min = k;
		}
		if(r <= heapSize && heap[r].getLength() < heap[min].getLength()) {
			min = r;
		} 
		if(min != k) {
			Edge temp = new Edge(heap[k]);
			heap[k] = new Edge(heap[min]);
			heap[min] = new Edge(temp);
			
			location[heap[k].getPa()][heap[k].getPb()] = k;
			location[heap[min].getPa()][heap[min].getPb()] = min;
			
			heapIFY(heap, location, min);
		}
	}

}
