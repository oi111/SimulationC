package Data;

public class TCor {
	double ax, ay, xy, ax2, ay2;
	public int tot;

	public TCor() {
		tot = 0;
		ax = ay = xy = ax2 = ay2 = 0;
	}

	public void init() {
		tot = 0;
		ax = ay = xy = ax2 = ay2 = 0;
	}

	public void add(double x, double y) {
		ax += x;
		ay += y;
		ax2 += x * x;
		ay2 += y * y;
		xy += x * y;
		tot++;
	}

	public double calCorrelation() {
		if (tot * ax2 - ax * ax == 0)
			return 0;
		return (tot * xy - ax * ay) / Math.sqrt((tot * ax2 - ax * ax) * (tot * ay2 - ay * ay));
	}
}
