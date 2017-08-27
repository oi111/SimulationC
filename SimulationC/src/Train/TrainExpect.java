package Train;

import Data.TK;
import OrderData.Expect;
import Uti.OutputFile;

public class TrainExpect {
	double DeltaX = 5e-5;
	OutputFile output;

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		TrainExpect t = new TrainExpect();
		// t.process(30, 400, 1000, 30, 10, 0.01, 0.00005, 0.000000001, 5000,
		// 101, 2, 0.6, 0.001, 1, 10000, 150, 50, 5000,
		// 0.0002, 2);
		t.process(Double.valueOf(args[0]), Double.valueOf(args[1]), Double.valueOf(args[2]), Double.valueOf(args[3]),
				Double.valueOf(args[4]), Double.valueOf(args[5]), Double.valueOf(args[6]), Double.valueOf(args[7]),
				Double.valueOf(args[8]), Integer.valueOf(args[9]), Integer.valueOf(args[10]), Double.valueOf(args[11]),
				Double.valueOf(args[12]), Double.valueOf(args[13]), Integer.valueOf(args[14]),
				Integer.valueOf(args[15]), Integer.valueOf(args[16]), Integer.valueOf(args[17]),
				Double.valueOf(args[18]), Double.valueOf(args[19]), Integer.valueOf(args[20]));
	}

	TK[] getFirstExpect(int nn, double sigmak) {
		double g[] = { 0, 1.394993763107846E-4, 2.7685714424764445, 0.015449367040358618, 3.0439195762885527,
				0.0017350596397439685, 21.622820373026034, 2.4235037612854265E-8 };
		TK ret[] = new TK[nn];
		for (int i = 0; i < nn; i++) {
			ret[i] = new TK(0, 0);
			ret[i].t2 = Math.exp(-g[6] * (i + 1) * DeltaX) * Math.pow((i + 1) * DeltaX + g[1], 1)
					/ Math.pow((i + 1) * DeltaX + g[5], 1.1) / 1.2 * sigmak;
			ret[i].t1 = Math.exp(-g[6] * (i + 1) * DeltaX) * Math.pow((i + 1) * DeltaX + g[1], g[2])
					* Math.pow((i + 1) * DeltaX + g[3], g[4]) / Math.pow((i + 1) * DeltaX + g[5], g[2] + g[4] + 3)
					* g[7];
		}
		return ret;
	}

	void getNextInput(TK a[], double b[]) {
		for (int i = 0; i < a.length; i++)
			a[i].t1 = (a[i].t1 * Expect.e[i] / b[i]);
	}

	void init() {
		output = new OutputFile();
		output.setFileName("Expect.txt");
		output.openFile();
	}

	void close() {
		output.closeFile();
	}

	void process(double alphad, double alphain, double alphak, double alphaout, double cc, double deltat, double deltax,
			double dk, double maxval, int nn, int num, double salpha, double sigmak, double sigmastable, int t,
			int tmax, int tn, int ts, double vk, double vsigma, int zz) {
		init();
		TK input[] = getFirstExpect(nn, sigmak);
		for (int i = 0; i < zz; i++) {
			ProDeltaNSimulationForTrainExpect p = new ProDeltaNSimulationForTrainExpect();
			p.process(alphad, alphain, alphak, alphaout, cc, deltat, deltax, dk, maxval, nn, num, salpha, sigmak,
					sigmastable, t, tmax, tn, ts, vk, vsigma, input);
			getNextInput(input, p.getOutputExpect());
			output(input);
		}
		close();
	}

	void output(TK input[]) {
		for (int i = 0; i < input.length; i++)
			output.write(input[i].t1 + " ");
		output.write("\n");
	}

}
