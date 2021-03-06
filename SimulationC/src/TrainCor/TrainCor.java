package TrainCor;

import Data.TK;
import OrderData.CorVAndDeltaN;
import Uti.InputFile;
import Uti.OutputFile;
import Uti.Uti;

public class TrainCor {
	double DeltaX = 5e-5;
	OutputFile output;
	double ZRate = 0.1;
	double ZLimit = 0.7;

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		TrainCor t = new TrainCor();
		// t.process(30, 400, 1000, 30, 10, 0.01, 0.00005, 0.000000001, 5000,
		// 101, 2, 0.6, 0.001, 1, 10000, 150, 50, 5000,
		// 0.0002, 2);
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
		// String line3 = input.read();
		// String pt3[] = line3.split(" ");
		for (int i = 0; i < nn; i++) {
			ret[i] = new TK(0, 0);
			ret[i].t2 = Math.exp(-g[6] * (i + 1) * DeltaX) * Math.pow((i + 1) * DeltaX + g[1], 1)
					/ Math.pow((i + 1) * DeltaX + g[5], 1.1) / 1.2 * sigmak;
			ret[i].t2 = Double.valueOf(pt2[i]);
			ret[i].t1 = Double.valueOf(pt1[i]);
			// ret[i].t3 = Double.valueOf(pt3[i]);
		}
		input.closeFile();
		return ret;
	}

	void getNextInput(TK a[], double b[]) {
		Uti.getSoft(b);
		for (int i = 0; i < a.length; i++) {
			if (CorVAndDeltaN.c[i] * ZLimit > b[i])
				a[i].t2 = a[i].t2 * ZRate + 0.001;
			else
				a[i].t2 /= ZRate;
			// if (a[i].t2 < 0)
			// a[i].t2 = 0.001;
			// a[i].t2 = (a[i].t2 * Math.abs(k));
		}
		// Uti.getSoft(a);
	}

	void init() {
		output = new OutputFile();
		output.setFileName("CorVAndDeltaN.txt");
		output.openFile();
	}

	void close() {
		output.closeFile();
	}

	void process(double alphad, double alphain, double alphak, double alphaout, double cc, double deltat, double deltax,
			double dk, double maxval, int nn, int num, double salpha, double sigmak, double sigmastable, int t,
			int tmax, int tn, int ts, double vk, double vsigma, double zlimit, int zn, double zrate, int zz) {
		this.ZRate = zrate;
		this.ZLimit = zlimit;
		init();
		TK input[] = getFirstExpect(nn, sigmak);
		for (int i = 0; i < zz; i++) {
			ProDeltaNSimulationForTrainCor p = new ProDeltaNSimulationForTrainCor();
			p.process(alphad, alphain, alphak, alphaout, cc, deltat, deltax, dk, maxval, nn, num, salpha, sigmak,
					sigmastable, t, tmax, tn, ts, vk, vsigma, input);
			getNextInput(input, p.getOutputCor(zn));
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
