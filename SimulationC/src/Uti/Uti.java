package Uti;

import java.util.List;

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
}
