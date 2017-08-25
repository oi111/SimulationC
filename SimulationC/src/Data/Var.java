package Data;

public class Var {
	public double n;
	double x, x2;

	public Var() {
		x = 0;
		x2 = 0;
		n = 0;
	}

	public void add(double tmp) {
		n++;
		x += tmp;
		x2 += tmp * tmp;
	}

	public double cal() {
		return (x2 / (n + 1e-12) - x * x / (n * n + 1e-12)) * n;
	}

	public double cal2() {
		return (x2 / (n + 1e-12) - x * x / (n * n + 1e-12));
	}
}
