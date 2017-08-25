package Data;

import java.util.ArrayList;
import java.util.List;

public class TradeSecond {
	public List<Trade> lt = new ArrayList<Trade>();
	public int max = 0, min = 1000000;

	public TradeSecond() {
		lt.clear();
	}

	public void add(Trade t) {
		lt.add(t);
		min = Math.min(t.lprice, min);
		max = Math.max(t.lprice, max);
	}
}
