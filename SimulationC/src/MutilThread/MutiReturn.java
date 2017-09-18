package MutilThread;

import java.util.ArrayList;
import java.util.List;

import Uti.OutputFile;

public class MutiReturn {
	double args[] = { 30, 400, 1000, 30, 10, 0.01, 0.00005, 0.0000000008, 5000, 11, 2, 0.6, 0.003, 9, 1000000, 200, 50,
			5000, 0.00015, 2 };
	// double args[] = { 30, 400, 1000, 30, 10, 0.1, 0.00005, 0.0000000008,
	// 1000, 11, 2, 0.6, 0.02, 1, 1000000, 200, 50,
	// 5000, 0.00007, 2 };
	int NN = 10;
	double pro[] = { 0.7, 0.8, 0.9 };
	List<List<Double>> ts = new ArrayList<List<Double>>();
	OutputFile output;
	double ans[];

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		MutiReturn m = new MutiReturn();
		// m.process(2, 0.00000005, 0.00000015, 0.00000001, 0.0002, 0.002,
		// 0.0001);
		// m.process(1, 0.000000001, 0.00000002, 0.000000001, 0.000005, 0.0001,
		// 0.000005);

		// 2 0.00000005 0.00000015 0.00000001 0.0002 0.002 0.0001
		// m.process(5, 18, 0.00001, 0.000040, 0.00001);
		m.process(Integer.valueOf(args[0]), Integer.valueOf(args[1]), Double.valueOf(args[2]), Double.valueOf(args[3]),
				Double.valueOf(args[4]));

	}

	void init() {
		output = new OutputFile();
		output.setFileName("MutiReturn.txt");
		output.openFile();
		ans = new double[4];
		ans[0] = 100000000;
	}

	void close() {
		output.closeFile();
	}

	void process(int nn, int index, double dk1, double dk2, double dk3) {

		System.out.println(nn + " " + dk1 + " " + dk2 + " " + dk3);
		this.NN = nn;
		init();

		run(index, dk1, dk2, dk3);
		// output.write(ans[0] + " " + ans[1] + " " + ans[2] + " " + ans[3] +
		// "\n");
		output.closeFile();
	}

	void run(int index, double dk1, double dk2, double dk3) {
		for (double i = dk1; i < dk2; i += dk3) {
			args[index] = i;
			runOne(i);
		}
	}

	void runOne(double val) {
		List<ProDeltaNSimulationInput> pr = new ArrayList<ProDeltaNSimulationInput>();
		pr.clear();
		for (int i = 0; i < NN; i++) {
			ProDeltaNSimulationInput p = new ProDeltaNSimulationInput();
			p.setParm(Double.valueOf(args[0]), Double.valueOf(args[1]), Double.valueOf(args[2]),
					Double.valueOf(args[3]), Double.valueOf(args[4]), Double.valueOf(args[5]), Double.valueOf(args[6]),
					Double.valueOf(args[7]), Double.valueOf(args[8]), (int) (args[9]), (int) (args[10]),
					Double.valueOf(args[11]), Double.valueOf(args[12]), Double.valueOf(args[13]), (int) (args[14]),
					(int) (args[15]), (int) (args[16]), (int) (args[17]), Double.valueOf(args[18]),
					Double.valueOf(args[19]), i);
			pr.add(p);
			Thread th = new Thread(pr.get(pr.size() - 1));
			th.start();
		}
		while (true) {
			boolean flg = true;
			for (int i = 0; i < pr.size(); i++)
				if (pr.get(i).AC == false)
					flg = false;
			if (flg)
				break;
		}
		double tot = 0;
		for (int i = 0; i < pr.size(); i++) {
			tot += pr.get(i).pdk;
			System.out.println(i + " " + val + " " + pr.get(i).pdk);
		}
		tot = tot / (pr.size() + 1e-8);
		System.out.println(pr.size() + " " + val + " " + tot);
		output.write(val + " " + tot + "\n");

		// output();
	}

}
