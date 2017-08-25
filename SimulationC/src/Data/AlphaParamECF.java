package Data;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class AlphaParamECF {
	double sumx, x2;
	public int num;
	public double alpha, gamma, beta;
	public double sin3, sin4, cos3, cos4, u, u0, u1, sigma0, sigma, sigma1;
	public double sin1, sin2, cos1, cos2;
	double t3 = 0.3, t4 = 0.4;

	List<Double> ld = new ArrayList<Double>();

	public AlphaParamECF() {
		init();
	}

	public void init() {
		sumx = 0;
		x2 = 0;
		num = 0;
		sigma = beta = alpha = gamma = u = 0;
		sin3 = sin4 = cos3 = cos4 = 0;
		u1 = sigma1 = 0;
	}

	public void add(double val) {
		ld.add(val);
	}

	public void cal() {
		if (ld.size() == 0)
			return;
		calInit();
		double tmp1 = sigma0, tmp2 = u0;
		u1 = u0;
		sigma1 = sigma0;
		for (int i = 1; i <= 5; i++) {
			u0 = 0;
			for (int j = 0; j < ld.size(); j++)
				ld.set(j, (ld.get(j) - u0) / sigma0);
			calK();
			calSigma();
			u1 = u1 + sigma1 * u;
			sigma1 = sigma * sigma1;
			sigma0 = sigma;
			u0 = u;
		}
		sigma = sigma1;
		u = u1;
	}

	void calInit() {
		Collections.sort(ld);
		double tot = 0, pt = 0;
		sigma0 = (ld.get((int) (ld.size() * 0.72)) - ld.get((int) (ld.size() * 0.28))) / 1.654;

		for (int i = (int) (ld.size() * 0.25); i < (int) (ld.size() * 0.75); i++, pt++) {
			tot += ld.get(i);
		}
		u0 = tot / (pt + 0.0);
	}

	void calK() {
		double t1 = 0.2;
		double t2 = 0.8;
		sin1 = sin2 = cos1 = cos2 = 0;
		for (int i = 0; i < ld.size(); i++) {
			double val = ld.get(i);
			sin1 += Math.sin(t1 * val);
			sin2 += Math.sin(t2 * val);
			cos1 += Math.cos(t1 * val);
			cos2 += Math.cos(t2 * val);
		}
		sin1 = sin1 / (ld.size() + 0.0);
		sin2 = sin2 / (ld.size() + 0.0);
		cos1 = cos1 / (ld.size() + 0.0);
		cos2 = cos2 / (ld.size() + 0.0);
		double tmp1 = Math.sqrt(sin1 * sin1 + cos1 * cos1);
		double tmp2 = Math.sqrt(sin2 * sin2 + cos2 * cos2);
		alpha = Math.log(Math.log(tmp1) / Math.log(tmp2)) / Math.log(t1 / t2);
		sigma = (Math.log(t1) * Math.log(-Math.log(tmp2)) - Math.log(t2) * Math.log(-Math.log(tmp1)))
				/ Math.log(t1 / t2);
		sigma = Math.exp(sigma);
	}

	void calSigma() {
		double t3 = 0.1;
		double t4 = 0.4;
		sin3 = sin4 = cos3 = cos4 = 0;
		for (int i = 0; i < ld.size(); i++) {
			double val = ld.get(i);
			sin3 += Math.sin(t3 * val);
			sin4 += Math.sin(t4 * val);
			cos3 += Math.cos(t3 * val);
			cos4 += Math.cos(t4 * val);
		}
		double tmp1 = Math.atan(sin3 / cos3);
		double tmp2 = Math.atan(sin4 / cos4);
		u = (Math.pow(t4, alpha) * tmp1 - Math.pow(t3, alpha) * tmp2)
				/ (t3 * Math.pow(t4, alpha) - t4 * Math.pow(t3, alpha));
		beta = (t4 * tmp1 - t3 * tmp2) / (Math.pow(sigma, alpha) * Math.tan(Math.PI / 2.0 * alpha)
				* (t4 * Math.pow(t3, alpha) - t3 * Math.pow(t4, alpha)));
	}

}
