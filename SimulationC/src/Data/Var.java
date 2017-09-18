package Data;

public class Var {
	public double n;
	double x, x2;
	public double n1, n2, n3;

	public Var() {
		x = 0;
		x2 = 0;
		n = 0;
		n1 = 0;
	}

	public void add(double tmp) {
		n++;
		x += tmp;
		x2 += tmp * tmp;
		if (tmp > 0.02)
			n1++;
		if (tmp > 0.1)
			n2++;
		if (tmp > 0.5)
			n3++;
	}

	public double getPro() {
		return n1 / (n + 1e-8);
	}

	public double cal() {
		return (x2 / (n + 1e-12) - x * x / (n * n + 1e-12)) * n;
	}

	public double cal2() {
		return (x2 / (n + 1e-12) - x * x / (n * n + 1e-12));
	}
}
