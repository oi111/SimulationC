package Data;

public class PowerLaw {
	double t1, tmin;
	int num;

	public PowerLaw() {
		t1 = 0;
		tmin = 1;
		num = 0;
	}

	public void add(double val) {
		t1 += Math.log(val);
		num++;
	}

	public double cal() {
		return 1 + num / (t1 - num * Math.log(tmin));
	}
}
