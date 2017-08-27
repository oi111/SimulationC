package Data;

public class CalK {
	double x2, xy, sumx, sumy;
	public int n = 0;
	public double k = 0, b = 0;

	public CalK() {
		x2 = xy = sumx = sumy = n = 0;
	}

	public void add(double x, double y) {
		n++;
		xy += x * y;
		sumx += x;
		sumy += y;
		x2 += x * x;
	}

	public void cal() {
		k = (n * xy - sumx * sumy) / (n * x2 - sumx * sumx + 1e-8);
		b = sumy / (n + 1e-8) - k * sumx / (n + 1e-8);
	}
}
