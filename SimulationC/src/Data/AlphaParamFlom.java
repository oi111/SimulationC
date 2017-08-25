package Data;

import java.util.ArrayList;
import java.util.List;

import Uti.Gamma;

public class AlphaParamFlom {
	double sumx, x2;
	public int num;
	public double alpha, gamma, beta;
	public double p;
	public double ap, a_p, sp, s_p;
	public double sin3, sin4, cos3, cos4, u, sigma;
	public double sin1, sin2, cos1, cos2;
	double t3 = 0.3, t4 = 0.4;
	Gamma ga = new Gamma();
	List<Double> ld = new ArrayList<Double>();

	public AlphaParamFlom() {

	}

	void init() {
		sumx = 0;
		x2 = 0;
		num = 0;
		sigma = beta = alpha = gamma = u = 0;
		p = 0.05;
		sin3 = sin4 = cos3 = cos4 = 0;
	}

	public void add(double val) {
		// if (val == 0)
		// return;
		if (val == 0)
			return;
		addFlom(val);
		sin3 += Math.sin(t3 * val);
		sin4 += Math.sin(t4 * val);
		cos3 += Math.cos(t3 * val);
		cos4 += Math.cos(t4 * val);
	}

	public void addFlom(double val) {
		ap += Math.pow(Math.abs(val), p);
		a_p += Math.pow(Math.abs(val), -p);
		sp += (val > 0 ? 1 : -1) * Math.pow(Math.abs(val), p);
		s_p += (val > 0 ? 1 : -1) * Math.pow(Math.abs(val), -p);
		sumx += val;
		x2 += val * val;
		num++;
	}

	public void cal() {
		calFlom();
		calU();
	}

	void calU() {
		u = (Math.pow(t4, alpha) * Math.atan(sin3 / cos3) - Math.pow(t3, alpha) * Math.atan(sin4 / cos4))
				/ (t3 * Math.pow(t4, alpha) - t4 * Math.pow(t3, alpha));
	}

	public void calFlom() {
		ap /= num;
		a_p /= num;
		sp /= num;
		s_p /= num;
		double q = p * Math.PI / 2.0;
		alpha = 1.0 / (q * (ap * a_p / Math.tan(q) + sp * s_p * Math.tan(q)));
		alpha = p * Math.PI / findValue(alpha);
		double theta = Math.asin(ap * s_p * alpha * Math.sin(p * Math.PI / alpha)) * alpha / 2.0 / p;
		beta = Math.tan(theta) / Math.tan(alpha * Math.PI / 2.0);

		double tmp = ga.gamma(1 - p) * Math.cos(p * Math.PI / 2.0) / ga.gamma(1 - p / alpha)
				* Math.cos(p * theta / alpha) * ap;
		gamma = Math.abs(Math.cos(theta)) * Math.pow(tmp, alpha / p);
		sigma = Math.pow(gamma, 1.0 / alpha);
	}

	double findValue(double val) {
		double left, right, mid = 0, tmp;
		left = 0;
		right = Math.PI;
		while (Math.abs(left - right) > 1e-8) {
			mid = (left + right) / 2.0;
			tmp = Math.sin(mid) / mid;
			if (tmp < val)
				right = mid;
			else
				left = mid;
		}
		return (left + right) / 2.0;
	}

}
