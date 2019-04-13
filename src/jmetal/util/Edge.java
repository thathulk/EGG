package jmetal.util;

public class Edge {
	
	private int pa;
	private int pb;
	private double length;
	
	public Edge() {}
	
	public Edge(Edge edge) {
		this.pa = edge.getPa();
		this.pb = edge.getPb();
		this.length = edge.getLength();
	}
	
	public Edge(int a, int b, double length) {
		this.pa = a;
		this.pb = b;
		this.length = length;
	}
	
	public int getPa() {
		return pa;
	}
	public void setPa(int pa) {
		this.pa = pa;
	}
	public int getPb() {
		return pb;
	}
	public void setPb(int pb) {
		this.pb = pb;
	}
	public double getLength() {
		return length;
	}
	public void setLength(double length) {
		this.length = length;
	}
	
	
	
}
