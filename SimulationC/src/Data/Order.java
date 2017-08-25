package Data;

public class Order {
	public long time;
	public double first, mid, volume;

	public Order(long time, double first, double mid, double volume) {
		this.time = time;
		this.first = first;
		this.mid = mid;
		this.volume = volume;
	}

}
