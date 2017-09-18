package Uti;

import java.util.List;

import Data.THR;
import Data.TK;

public class Uti {

	public static double calSigma(List<Double> ld) {
		double averager = calAverage(ld);
		double tot = 0;
		for (int i = 0; i < ld.size(); i++)
			tot += (ld.get(i) - averager) * (ld.get(i) - averager);
		tot /= (ld.size() + 1e-8);
		return Math.sqrt(tot);
	}

	public static double calAverage(List<Double> ld) {
		double tot = 0;
		for (int i = 0; i < ld.size(); i++)
			tot += ld.get(i);
		return tot / (ld.size() + 1e-8);
	}

	public static void getSoft(double b[]) {
		double ret[] = new double[b.length];
		for (int i = 0; i < b.length; i++) {
			int k = 0;
			double tot = 0;
			for (int j = -5; j < 5; j++)
				if (i + j >= 0 && i + j < b.length) {
					tot += b[i + j];
					k++;
				}
			ret[i] = tot / (k + 1e-8);
		}
		for (int i = 0; i < b.length; i++)
			b[i] = ret[i];
	}

	public static void getSoft(TK b[]) {
		double ret[] = new double[b.length];
		for (int i = 0; i < b.length; i++) {
			int k = 0;
			double tot = 0;
			for (int j = -5; j < 5; j++)
				if (i + j >= 0 && i + j < b.length) {
					tot += b[i + j].t2;
					k++;
				}
			ret[i] = tot / (k + 1e-8);
		}
		for (int i = 0; i < b.length; i++)
			b[i].t2 = ret[i];
	}

	public static void getSoft2(TK b[]) {
		double tmp[] = new double[b.length];
		for (int g = 1; g <= 10; g++) {
			for (int j = 0; j < b.length; j++) {
				int len = 0;
				if (j < 50)
					len = 1;
				else if (j < 150)
					len = 3;
				else
					len = 7;
				double tot = 0;
				int t1 = 0;
				for (int k = -len; k <= len; k++)
					if (j + k >= 0 && j + k < b.length) {
						t1++;
						tot += (b[j + k].t2);
					}
				tmp[j] = tot / t1;
			}
			for (int j = 0; j < b.length; j++)
				b[j].t2 = tmp[j];
		}
	}

	public static void getSoft3(THR b[]) {
		double tmp[] = new double[b.length];
		for (int g = 1; g <= 10; g++) {
			for (int j = 0; j < b.length; j++) {
				int len = 0;
				if (j < 50)
					len = 1;
				else if (j < 150)
					len = 3;
				else
					len = 7;
				double tot = 0;
				int t1 = 0;
				for (int k = -len; k <= len; k++)
					if (j + k >= 0 && j + k < b.length) {
						t1++;
						tot += (b[j + k].t2);
					}
				tmp[j] = tot / t1;
			}
			for (int j = 0; j < b.length; j++)
				b[j].t2 = tmp[j];
		}
	}
}
