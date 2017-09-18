package TrainAll;

import Data.THR;
import OrderData.CorDeltaN;
import OrderData.CorVAndDeltaN;
import OrderData.Expect;
import Uti.InputFile;
import Uti.OutputFile;
import Uti.Uti;

public class TrainAllBeginX2 {
	double DeltaX = 5e-5;
	double Zlimit = 1;
	double ZlimitD = 1;
	OutputFile output;
	double ZRate = 1.05;
	double ZRateD = 1.1;

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		TrainAllBeginX2 t = new TrainAllBeginX2();
		// t.process(30, 400, 1000, 30, 10, 0.01, 0.00005, 0.000000001, 5000,
		// 101, 2, 0.6, 0.001, 1, 10000, 150, 50, 5000,
		// 0.0002, 2);
		t.process(Double.valueOf(args[0]), Double.valueOf(args[1]), Double.valueOf(args[2]), Double.valueOf(args[3]),
				Double.valueOf(args[4]), Double.valueOf(args[5]), Double.valueOf(args[6]), Double.valueOf(args[7]),
				Double.valueOf(args[8]), Integer.valueOf(args[9]), Integer.valueOf(args[10]), Double.valueOf(args[11]),
				Double.valueOf(args[12]), Double.valueOf(args[13]), Integer.valueOf(args[14]),
				Integer.valueOf(args[15]), Integer.valueOf(args[16]), Integer.valueOf(args[17]),
				Double.valueOf(args[18]), Double.valueOf(args[19]), Double.valueOf(args[20]), Double.valueOf(args[21]),
				Integer.valueOf(args[22]), Double.valueOf(args[23]), Double.valueOf(args[24]),
				Integer.valueOf(args[25]));
	}

	THR[] getFirstExpect(double dk, int nn, double sigmak) {
		double g[] = { 0, 1.394993763107846E-4, 2.7685714424764445, 0.015449367040358618, 3.0439195762885527,
				0.0017350596397439685, 21.622820373026034, 2.4235037612854265E-8 };
		THR ret[] = new THR[nn];
		InputFile input = new InputFile();
		input.setFileName("Expect.txt");
		input.openFile();
		String line1 = input.read();
		String pt1[] = line1.split(" ");
		String line2 = input.read();
		String pt2[] = line2.split(" ");
		String line3 = input.read();
		String pt3[] = line3.split(" ");
		for (int i = 0; i < nn; i++) {
			ret[i] = new THR(0, 0, 0);
			// ret[i].t2 = Double.valueOf(pt2[i]);
			ret[i].t2 = Math.exp(-g[6] * (i + 1) * DeltaX) * Math.pow((i + 1) * DeltaX + g[1], 1)
					/ Math.pow((i + 1) * DeltaX + g[5], 1.1) / 1.2 * sigmak;
			ret[i].t1 = Double.valueOf(Expect.e[i]);
			ret[i].t2 = Double.valueOf(pt2[i]);
			ret[i].t1 = Double.valueOf(pt1[i]);
			ret[i].t3 = 2.6 * Math.pow((i + 1) * DeltaX, 2) / Math.pow((i + 1) * DeltaX + 0.00001, 2)
					* Math.exp(-30 * i * DeltaX) * dk;
			ret[i].t3 = Double.valueOf(pt3[i]);

		}
		input.closeFile();
		return ret;
	}

	void getNextInputExpect(THR p[], double a[]) {
		// Uti.getSoft(b);
		for (int i = 0; i < p.length; i++) {
			p[i].t1 = (p[i].t1 * Expect.e[i] / a[i]);
		}
	}

	void getNextInputCor(THR p[], double b[]) {
		Uti.getSoft(b);
		for (int i = 0; i < p.length; i++) {
			if (CorVAndDeltaN.c[i] * Zlimit > b[i])
				p[i].t2 = p[i].t2 * ZRate + 0.001;
			else
				p[i].t2 /= ZRate;
		}
	}

	void getNextInputD(THR p[], double c[]) {
		Uti.getSoft(c);
		for (int i = 0; i < p.length; i++) {
			if (CorDeltaN.dn[i] * ZlimitD > c[i])
				p[i].t3 = p[i].t3 * ZRateD + 1e-9;
			else {
				if (CorDeltaN.dn[i] * 3 < c[i])
					p[i].t3 /= 2;
				else
					p[i].t3 /= ZRateD;
			}
		}
	}

	void init() {
		output = new OutputFile();
		output.setFileName("TrainAll.txt");
		output.openFile();
	}

	void close() {
		output.closeFile();
	}

	void process(double alphad, double alphain, double alphak, double alphaout, double cc, double deltat, double deltax,
			double dk, double maxval, int nn, int num, double salpha, double sigmak, double sigmastable, int t,
			int tmax, int tn, int ts, double vk, double vsigma, double zlimit, double zlimitd, int zn, double zrate,
			double zrated, int zz) {
		this.ZRate = zrate;
		this.Zlimit = zlimit;
		this.ZlimitD = zlimitd;
		init();
		THR input[] = getFirstExpect(dk, nn, sigmak);
		for (int i = 0; i < zz; i++) {
			SimulationForTrainCorAndExpect p = new SimulationForTrainCorAndExpect();
			// if (i % 50 < 10)
			p.process(alphad, alphain, alphak, alphaout, cc, deltat, deltax, dk, maxval, nn, num, salpha, sigmak,
					sigmastable, t, tmax, tn, ts, vk, vsigma, input);
			// else
			// p.process(alphad, alphain, alphak, alphaout, cc, deltat, deltax,
			// dk, maxval, nn, num, salpha, sigmak,
			// sigmastable, 50000, tmax, tn, ts, vk, vsigma, input);
			// if (i % 20 == 0)
			// getSoft(input);
			if (i % 20 == 10)
				Uti.getSoft3(input);
			// if (i % 50 < 10)
			// getNextInputD(input, p.getOutputCorPos());
			if (i % 20 < 10)
				getNextInputCor(input, p.getOutputCor(zn));
			if (i % 20 >= 10)
				getNextInputExpect(input, p.getOutputExpect());

			output(p.getOutputExpect());
			output(p.getOutputCor(zn));
			output(p.getOutputCorPos());
			output(input);
		}
		close();
	}

	public void getSoft(THR b[]) {
		double ret[] = new double[b.length];
		for (int i = 0; i < b.length; i++) {
			int k = 0;
			double tot = 0;
			for (int j = -5; j < 5; j++)
				if (i + j >= 0 && i + j < b.length) {
					tot += b[i + j].t3;
					k++;
				}
			ret[i] = tot / (k + 1e-8);
		}
		for (int i = 0; i < b.length; i++)
			b[i].t3 = ret[i];
	}

	void output(THR input[]) {
		for (int i = 0; i < input.length; i++)
			output.write(input[i].t1 + " ");
		output.write("\n");
		for (int i = 0; i < input.length; i++)
			output.write(input[i].t2 + " ");
		output.write("\n");
		for (int i = 0; i < input.length; i++)
			output.write(input[i].t3 + " ");
		output.write("\n");
	}

	void output(double input[]) {
		for (int i = 0; i < input.length; i++)
			output.write(input[i] + " ");
		output.write("\n");
	}

}
