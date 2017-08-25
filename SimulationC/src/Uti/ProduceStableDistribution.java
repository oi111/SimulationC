package Uti;

import java.util.ArrayList;
import java.util.List;

public class ProduceStableDistribution {

	int NN = 1000000;
	double alpha = 0.5, beta = 1, gamma = 1, mu = 0;
	double data[];// = new double[NN];
	double max_val = 1e4;

	public void setParam(double alpha, double beta, double gamma, double mu) {
		this.alpha = alpha;
		this.beta = beta;
		this.gamma = gamma;
	}

	void setNN(int nn) {
		this.NN = nn;
	}

	void process() {
		data = new double[NN];
		List<Double> ld = new ArrayList<Double>();
		double tot = 0;
		for (int i = 0; i < NN; i++) {
			while (true) {
				double tv = (Math.random() - 0.5) * Math.PI;
				double tw = -Math.log(Math.random());
				double tm = (1 + beta * beta * Math.tan(Math.PI / 2.0 * alpha) * Math.tan(Math.PI / 2.0 * alpha));
				tm = Math.pow(tm, 1.0 / (2 * alpha));
				double tn = -Math.atan(beta * Math.tan(Math.PI * alpha / 2.0)) / alpha;
				double t1 = Math.sin(alpha * (tv - tn)) / Math.pow(Math.cos(tv), 1.0 / alpha);
				double t2 = Math.pow(Math.cos(tv - alpha * (tv - tn)) / tw, (1 - alpha) / alpha);
				data[i] = (t1 * t2 * tm * gamma + mu);
				System.out.println(data[i]);
				// if (Math.abs(data[i]) < 10000)
				break;
			}
			ld.add(data[i]);
			tot += Math.abs(data[i]);
			// System.out.println(data[i]);
		}
		System.out.println(tot / NN);
		PdfOutput pdf = new PdfOutput();
		pdf.setParam(1000, -1000, 1000);
		pdf.openFile("xxx.txt");
		pdf.calPdf(ld);
		pdf.closeFile();

	}

	public double getNext() {
		double tv = (Math.random() - 0.5) * Math.PI;
		double tw = -Math.log(Math.random());
		double tm = (1 + beta * beta * Math.tan(Math.PI / 2.0 * alpha) * Math.tan(Math.PI / 2.0 * alpha));
		tm = Math.pow(tm, 1.0 / (2 * alpha));
		double tn = -Math.atan(beta * Math.tan(Math.PI * alpha / 2.0)) / alpha;
		double t1 = Math.sin(alpha * (tv - tn)) / Math.pow(Math.cos(tv), 1.0 / alpha);
		double t2 = Math.pow(Math.cos(tv - alpha * (tv - tn)) / tw, (1 - alpha) / alpha);
		// if (t1 * t2 * tm * gamma + mu < 0)
		// System.out.println(t1 * t2 * tm * gamma + mu);
		return (t1 * t2 * tm * gamma + mu);
	}

	public double getPosNext() {

		while (true) {
			double tmp = getNext();
			if (tmp >= 0 && tmp < max_val)
				// if (tmp >= 0)
				return tmp;
		}
	}

	public void setMax(double max) {
		this.max_val = max;
	}

}
