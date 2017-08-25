package Data;

import java.util.ArrayList;
import java.util.List;

public class LstNode {
	public double price;
	public List<Order> lo = new ArrayList<Order>();

	public LstNode(double val) {
		this.price = val;
		lo.clear();
	}

	public void add(long time, double first, double mid, double volume) {
		if (lo.size() == 0) {
			lo.add(new Order(time, first, mid, volume));
			return;
		}
		if (Math.abs(volume - lo.get(lo.size() - 1).volume) < 1e-6)
			return;
		lo.add(new Order(time, first, mid, volume));
	}

}
