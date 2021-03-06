package TrainCorAndExpect;

import Data.TK;
import OrderData.CorVAndDeltaN;
import OrderData.Expect;
import Uti.InputFile;
import Uti.OutputFile;

public class TrainCorAndExpectBeginX {
	double DeltaX = 5e-5;
	double Zlimit = 1;
	OutputFile output;
	double ZRate = 0.1;

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		TrainCorAndExpectBeginX t = new TrainCorAndExpectBeginX();
		// t.process(30, 400, 1000, 30, 10, 0.01, 0.00005, 0.000000001, 5000,
		// 101, 2, 0.6, 0.001, 1, 10000, 150, 50, 5000,
		// 0.0002, 2);
		// t.process(30, 400, 1000, 30, 10, 0.01, 0.00005, 0.0000000008, 1000,
		// 100, 2, 0.6, 0.02, 1, 50000, 200, 50, 5000,
		// 0.00007, 2, 1, 1, 1.01, 1000);
		t.process(Double.valueOf(args[0]), Double.valueOf(args[1]), Double.valueOf(args[2]), Double.valueOf(args[3]),
				Double.valueOf(args[4]), Double.valueOf(args[5]), Double.valueOf(args[6]), Double.valueOf(args[7]),
				Double.valueOf(args[8]), Integer.valueOf(args[9]), Integer.valueOf(args[10]), Double.valueOf(args[11]),
				Double.valueOf(args[12]), Double.valueOf(args[13]), Integer.valueOf(args[14]),
				Integer.valueOf(args[15]), Integer.valueOf(args[16]), Integer.valueOf(args[17]),
				Double.valueOf(args[18]), Double.valueOf(args[19]), Double.valueOf(args[20]), Integer.valueOf(args[21]),
				Double.valueOf(args[22]), Integer.valueOf(args[23]));
	}

	TK[] getFirstExpect(int nn, double sigmak) {
		double g[] = { 0, 1.394993763107846E-4, 2.7685714424764445, 0.015449367040358618, 3.0439195762885527,
				0.0017350596397439685, 21.622820373026034, 2.4235037612854265E-8 };
		TK ret[] = new TK[nn];
		InputFile input = new InputFile();
		input.setFileName("Expect.txt");
		input.openFile();
		String line1 = input.read();
		String pt1[] = line1.split(" ");
		String line2 = input.read();
		String pt2[] = line2.split(" ");
		for (int i = 0; i < nn; i++) {
			ret[i] = new TK(0, 0);
			// ret[i].t2 = Double.valueOf(pt2[i]);
			ret[i].t2 = Math.exp(-g[6] * (i + 1) * DeltaX) * Math.pow((i + 1) * DeltaX + g[1], 1)
					/ Math.pow((i + 1) * DeltaX + g[5], 1.1) / 1.2 * sigmak;
			ret[i].t1 = Double.valueOf(Expect.e[i]);
			ret[i].t2 = Double.valueOf(pt2[i]);
			ret[i].t1 = Double.valueOf(pt1[i]);
		}
		input.closeFile();
		return ret;
	}

	void getNextInput(TK p[], double a[], double b[]) {
		// Uti.getSoft(b);
		for (int i = 0; i < p.length; i++) {
			if (CorVAndDeltaN.c[i] * Zlimit > b[i])
				p[i].t2 *= ZRate;
			else
				p[i].t2 /= ZRate;

			if (Expect.e[i] / a[i] >= 0.98 && Expect.e[i] / a[i] <= 1.02)
				p[i].t1 = (p[i].t1 * Expect.e[i] / a[i]);
			else {
				if (Expect.e[i] > a[i])
					p[i].t1 *= 1.02;
				else
					p[i].t1 /= 1.02;
			}
		}
	}

	void init() {
		output = new OutputFile();
		output.setFileName("CorAndExpect.txt");
		output.openFile();
	}

	void close() {
		output.closeFile();
	}

	void process(double alphad, double alphain, double alphak, double alphaout, double cc, double deltat, double deltax,
			double dk, double maxval, int nn, int num, double salpha, double sigmak, double sigmastable, int t,
			int tmax, int tn, int ts, double vk, double vsigma, double zlimit, int zn, double zrate, int zz) {
		this.ZRate = zrate;
		this.Zlimit = zlimit;
		init();
		TK input[] = getFirstExpect(nn, sigmak);
		for (int i = 0; i < zz; i++) {
			SimulationForTrainCorAndExpect p = new SimulationForTrainCorAndExpect();
			p.process(alphad, alphain, alphak, alphaout, cc, deltat, deltax, dk, maxval, nn, num, salpha, sigmak,
					sigmastable, t, tmax, tn, ts, vk, vsigma, input);
			getNextInput(input, p.getOutputExpect(), p.getOutputCor(zn));
			output(p.getOutputExpect());
			output(p.getOutputCor(zn));
			output(input);
		}
		close();
	}

	void output(TK input[]) {
		for (int i = 0; i < input.length; i++)
			output.write(input[i].t1 + " ");
		output.write("\n");
		for (int i = 0; i < input.length; i++)
			output.write(input[i].t2 + " ");
		output.write("\n");
	}

	void output(double input[]) {
		for (int i = 0; i < input.length; i++)
			output.write(input[i] + " ");
		output.write("\n");
	}

}
