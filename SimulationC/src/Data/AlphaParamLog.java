package Data;

public class AlphaParamLog {
	double sumx, x2;
	public int num;
	public double alpha, gamma, beta;
	public double sin3, sin4, cos3, cos4, u, sigma;
	double t3 = 0.3, t4 = 0.4;

	public AlphaParamLog() {
		init();
	}

	void init() {
		sumx = 0;
		x2 = 0;
		num = 0;
		sigma = beta = alpha = gamma = u = 0;
		sin3 = sin4 = cos3 = cos4 = 0;
		alpha = beta = 0;
	}

	public void add(double val) {
		if (val == 0)
			return;
		// add(val);
		addLog(val);
		sin3 += Math.sin(t3 * val);
		sin4 += Math.sin(t4 * val);
		cos3 += Math.cos(t3 * val);
		cos4 += Math.cos(t4 * val);

	}

	public void cal() {
		calLog();
		calU();
	}

	void calU() {
		u = (Math.pow(t4, alpha) * Math.atan(sin3 / cos3) - Math.pow(t3, alpha) * Math.atan(sin4 / cos4))
				/ (t3 * Math.pow(t4, alpha) - t4 * Math.pow(t3, alpha));
	}

	public void addLog(double val) {

		val = Math.log(Math.abs(val));
		sumx += val;
		x2 += val * val;
		num++;
	}

	public void calLog() {
		if (num <= 1)
			return;
		double euler = 0.57721566;
		double average = sumx / num;
		double delt = (x2 - sumx * sumx / num) / (num - 1.0);
		alpha = 1.0 / Math.sqrt(delt * 6.0 / (Math.PI * Math.PI) - 0.5);
		gamma = Math.exp((average - euler * (1.0 / alpha - 1.0)) * alpha);
		sigma = Math.pow(gamma, 1.0 / alpha);
	}

}
