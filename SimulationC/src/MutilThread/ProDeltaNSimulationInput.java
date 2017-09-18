package MutilThread;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import Data.CalK;
import Data.PosInfo;
import Data.TCor;
import Data.TK;
import Data.Var;
import Uti.InputFile;
import Uti.OutputFile;
import Uti.PdfOutput;
import Uti.ProduceStableDistribution;
import cern.jet.random.StudentT;
import cern.jet.random.engine.RandomEngine;

public class ProDeltaNSimulationInput implements Runnable {
	boolean AC = false;
	double pdk, pdb;
	int INDEX;

	@Override
	public void run() {
		// TODO Auto-generated method stub
		process();
		AC = true;
	}

	int NUM_COR = 10;
	int NN = 3;
	int Num = 10;
	int T;
	int TMax = 100;
	double DeltaT;
	double DeltaX;
	double AlphaIn;
	double AlphaOut;
	double AlphaK;
	double AlphaD;
	double SAlpha;
	double VSigma;
	double MaxVal;
	double CC;
	double VK;
	double SigmaK;
	double SigmaStable;
	int BIN = 1;
	int LEN = 200;
	double DK = 1e-7;
	int BINDEX = 0;
	int TNUM = 500;
	int TS = 100;
	// int AllT;
	PosInfo ap[], ap2[];
	ProduceStableDistribution ps;
	List<List<Double>> lld = new ArrayList<List<Double>>();
	List<List<Double>> lln = new ArrayList<List<Double>>();
	List<Double> lv = new ArrayList<Double>();
	List<TK> ltk = new ArrayList<TK>();
	List<TK> p1 = new ArrayList<TK>();
	List<TK> p2 = new ArrayList<TK>();
	List<TK> p3 = new ArrayList<TK>();
	List<TK> p4 = new ArrayList<TK>();
	List<TK> r1 = new ArrayList<TK>();
	List<TK> r2 = new ArrayList<TK>();
	List<TK> r3 = new ArrayList<TK>();
	List<List<Var>> lvar = new ArrayList<List<Var>>();
	double velocity[];
	TCor mt[][];
	TCor mt2[][];
	TK llt[][];
	TK llt2[][];
	double ord[][];
	List<TCor> ltc = new ArrayList<TCor>();
	List<List<TCor>> ltc2 = new ArrayList<List<TCor>>();
	Random gr = new Random();
	StudentT st;
	OutputFile output_d, output_expect, output_p, output_corpos, output_cor, output_n, output_en, output_tn, output_v,
			output_var;
	PdfOutput pdf1, pdf2, pdf_v;
	double current_v;

	PosInfo[] init() {
		PosInfo ap[] = initOther();
		initStatic();
		initF(ap);
		initD(ap);
		initN(ap);
		return ap;
	}

	void initStatic() {
		lld.clear();
		lln.clear();
		ltk.clear();
		p1.clear();
		p2.clear();
		p3.clear();
		p4.clear();
		r1.clear();
		r2.clear();
		r3.clear();
		ltc.clear();
		ltc2.clear();
		lv.clear();
		lvar.clear();
		mt = new TCor[NN][NN];
		mt2 = new TCor[NN][NN];
		for (int i = 0; i < NN; i++) {
			lld.add(new ArrayList<Double>());
			lln.add(new ArrayList<Double>());
			ltk.add(new TK(0, 0));
			p1.add(new TK(0, 0));
			p2.add(new TK(0, 0));
			p3.add(new TK(0, 0));
			p4.add(new TK(0, 0));
			r1.add(new TK(0, 0));
			r2.add(new TK(0, 0));
			r3.add(new TK(0, 0));
			ltc.add(new TCor());

		}
		for (int i = 0; i < NUM_COR; i++) {
			List<TCor> tmp = new ArrayList<TCor>();
			for (int j = 0; j < NN; j++)
				tmp.add(new TCor());
			ltc2.add(tmp);
		}
		for (int i = 0; i < mt.length; i++)
			for (int j = 0; j < mt[i].length; j++) {
				mt[i][j] = new TCor();
				mt2[i][j] = new TCor();
			}
		llt = new TK[NN / BIN + 1][LEN];
		for (int i = 0; i < llt.length; i++)
			for (int j = 0; j < llt[i].length; j++)
				llt[i][j] = new TK(0, 0);
		llt2 = new TK[NN / BIN + 1][LEN];
		for (int i = 0; i < llt2.length; i++)
			for (int j = 0; j < llt2[i].length; j++)
				llt2[i][j] = new TK(0, 0);
		for (int i = 0; i < 200; i++) {
			List<Var> tmp = new ArrayList<Var>();
			for (int j = 0; j < NN; j++)
				tmp.add(new Var());
			lvar.add(tmp);
		}
	}

	PosInfo[] initOther() {
		ps = new ProduceStableDistribution();
		st = new StudentT(VSigma, RandomEngine.makeDefault());
		ps.setParam(SAlpha, 0, SigmaStable, 0);
		ps.setMax(MaxVal);
		PosInfo ap[];
		ap = new PosInfo[NN];
		velocity = new double[300000];
		for (int i = 0; i < NN; i++)
			ap[i] = new PosInfo();
		ord = new double[100000][NN];
		for (int i = 0; i < ord.length; i++)
			for (int j = 0; j < ord[i].length; j++)
				ord[i][j] = 0;
		return ap;
	}

	void initN(PosInfo ap[]) {
		for (int i = 0; i < ap.length; i++) {
			ap[i].n = 25;
			ap[i].h = Math.log(25);
		}
	}

	void initF(PosInfo ap[]) {
		InputFile input = new InputFile();
		input.setFileName("Expect.txt");
		input.openFile();
		double g[] = { 0, 1.394993763107846E-4, 2.7685714424764445, 0.015449367040358618, 3.0439195762885527,
				0.0017350596397439685, 21.622820373026034, 2.4235037612854265E-8 };
		String line = input.read();
		String pt[] = line.split(" ");
		for (int i = 0; i < ap.length; i++) {
			ap[i].fout = Math.exp(-g[6] * i * DeltaX) * Math.pow(i * DeltaX + g[1], 1)
					/ Math.pow(i * DeltaX + g[5], 1.1) / 1.2 * SigmaK;
			ap[i].fin = Double.valueOf(pt[i]) * ap[i].fout;

		}
		input.closeFile();
	}

	void initD(PosInfo ap[]) {
		for (int i = 0; i < ap.length; i++) {
			// ap[i].d = 1.2 * Math.pow((i + 1) * DeltaX, 4) / Math.pow((i + 1)
			// * DeltaX + 0.0006, 4)
			// * Math.exp(-AlphaD * i * DeltaX) * DK;
			ap[i].d = 2.6 * Math.pow((i + 1) * DeltaX, 2) / Math.pow((i + 1) * DeltaX + 0.00001, 2)
					* Math.exp(-AlphaD * i * DeltaX) * DK;
			// ap[i].d = 2.6 * Math.pow(i * DeltaX, 2) / Math.pow(i * DeltaX +
			// 0.00001, 2) * Math.exp(-AlphaD * i * DeltaX)
			// * DK;
			// ap[i].d = 2.6 * Math.pow((i + 1) * DeltaX, 5) / Math.pow((i + 1)
			// * DeltaX + 0.0006, 5)
			// * Math.exp(-AlphaD * i * DeltaX) * DK;
			// double x = (i + 1) * DeltaX;
			// ap[i].d = 1.0 / (1 + Math.exp(-AlphaK * (x - 5e-5 * 40))) *
			// Math.exp(-AlphaD * x) * DK;
		}
		// for (int i = 0; i < ap.length; i++)
		// System.out.print(ap[i].d + " ");
		// System.out.println();
		double a[] = new double[41];
		for (int i = 0; i <= 20; i++)
			a[i] = 0.2 / DK * 5 * 1e-8;// 0.000002;
		for (int i = 21; i < 40; i++)
			a[i] = 0.05 * (i - 20);
		// for (int i = 0; i < 40; i++)
		// ap[i].d = ap[i].d * a[i];
	}

	void openOutput(String file) {

		output_d = new OutputFile();
		output_d.setFileName(file + "Data.txt");
		output_d.openFile();

		output_expect = new OutputFile();
		output_expect.setFileName(file + "Expect.txt");
		output_expect.openFile();

		output_p = new OutputFile();
		output_p.setFileName(file + "P.txt");
		output_p.openFile();

		output_corpos = new OutputFile();
		output_corpos.setFileName(file + "CorPos.txt");
		output_corpos.openFile();

		output_cor = new OutputFile();
		output_cor.setFileName(file + "Cor.txt");
		output_cor.openFile();

		output_n = new OutputFile();
		output_n.setFileName(file + "N.txt");
		output_n.openFile();

		output_tn = new OutputFile();
		output_tn.setFileName(file + "TN.txt");
		output_tn.openFile();

		output_en = new OutputFile();
		output_en.setFileName(file + "EN.txt");
		output_en.openFile();

		output_v = new OutputFile();
		output_v.setFileName(file + "V.txt");
		output_v.openFile();

		output_var = new OutputFile();
		output_var.setFileName(file + "Var.txt");
		output_var.openFile();

		pdf1 = new PdfOutput();
		pdf1.setParam(100, -100, 100);
		pdf1.openFile(file + ".txt");
		pdf2 = new PdfOutput();
		pdf2.setParam(140, -5, 2);
		pdf2.openFile(file + "Log.txt");
		pdf_v = new PdfOutput();
		pdf_v.setParam(100, -7, -2);
		pdf_v.openFile(file + "VLog.txt");
	}

	void closeOutput() {
		output_d.closeFile();
		output_expect.closeFile();
		output_p.closeFile();
		output_corpos.closeFile();
		output_cor.closeFile();
		output_n.closeFile();
		output_tn.closeFile();
		output_en.closeFile();
		output_v.closeFile();
		output_var.closeFile();
		pdf1.closeFile();
		pdf2.closeFile();
		pdf_v.closeFile();
	}

	void setParm(double alphad, double alphain, double alphak, double alphaout, double cc, double deltat, double deltax,
			double dk, double maxval, int nn, int num, double salpha, double sigmak, double sigmastable, int t,
			int tmax, int tn, int ts, double vk, double vsigma, int index) {
		this.AlphaD = alphad;
		this.AlphaIn = alphain;
		this.AlphaOut = alphaout;
		this.AlphaK = alphak;
		this.T = t;
		this.Num = num;
		this.DeltaT = deltat;
		this.DeltaX = deltax;
		this.VSigma = vsigma;
		this.NN = nn;
		this.DK = dk;
		this.MaxVal = maxval;
		this.CC = cc;
		this.TNUM = tn;
		this.TS = ts;
		this.VK = vk;
		this.SigmaK = sigmak;
		this.TMax = tmax;
		this.SAlpha = salpha;
		this.SigmaStable = sigmastable;
		this.INDEX = index;
		ap = init();
		ap2 = init();

	}

	void process() {
		String file = "ProDeltaNSimulationVelocity_" + "AD" + (int) AlphaD + "AI" + (int) AlphaIn + "AK" + (int) AlphaK
				+ "AO" + (int) AlphaOut + "CC" + CC + "DK" + (int) (Math.log10(DK)) + "DT" + (int) (Math.log10(DeltaT))
				+ "DX" + (int) (Math.log10(DeltaX)) + "MV" + MaxVal + "NN" + NN + "NUM" + Num + "T" + T + "TM" + TMax
				+ "TN" + TNUM + "VS" + (int) Math.log10(VSigma) + "_" + INDEX;
		openOutput(file);
		outputTheoryExpect();
		for (int i = 0; i < Num; i++) {
			initStatic();
			for (int j = 0; j < T; j++) {
				processOne(j);
				// if (j % 10000 == 0)
				// System.out.println((i * T + j));
				// System.out.println(ap[0].n + "\t" + ap[1].n + "\t" + ap[2].n
				// + "\t" + ap[197].n + "\t" + ap[198].n
				// + "\t" + ap[199].n);
			}
			outputExpect();
			outputP();
			outputCorPos();
			outputCor();
			outputN();
		}
		outputCumPdf();
		outputPdfV();
		outputV();
		outputVar();
		output(file);
		outputExampleN();
		closeOutput();
	}

	void processOne(int index) {
		double tot = 0;
		tot = ap[0].n + ap2[0].n;

		// current_v = (gr.nextGaussian() * VSigma + ((ap2[1].h - ap2[0].h) -
		// (ap[1].h - ap[0].h)) / DeltaX * ap2[0].d)
		// / (tot + 1e-15);

		double t0 = 0;
		while (true) {
			t0 = st.nextDouble();
			if (Math.abs(t0) < TMax)
				break;
		}
		// t0 = 0;
		double t1 = t0 * VK;
		double t2 = ((ap2[1].n - ap2[0].n) - (ap[1].n - ap[0].n)) / DeltaX * ap2[0].d;
		current_v = (t1 + t2) / (tot + 1e-15);

		// System.out.println(tot + " " + (ap2[1].n - ap2[0].n) + " " + (ap[1].n
		// - ap[0].n) + " " + t0 + " " + t1 + " "
		// + t2 + " " + current_v);
		int pt = index % velocity.length;
		velocity[pt] = velocity[(pt + velocity.length - 1) % velocity.length] + current_v;
		lv.add((velocity[pt] - velocity[(pt + velocity.length - TNUM) % velocity.length]));
		if (current_v > 1e-3)
			output_d.write(tot + " " + (ap2[1].n - ap2[0].n) + " " + (ap[1].n - ap[0].n) + " " + t0 + " " + t1 + " "
					+ t2 + " " + current_v + "\n");
		AT();
		BT();
		Next(index);
		AP2();

	}

	void BT() {
		for (int i = 0; i < ap.length - 1; i++) {
			ap[i].bt = ps.getPosNext() * DeltaT * ap[i].fin;
		}
	}

	void AT() {

		for (int i = 0; i < ap.length - 1; i++) {

			double t1 = 0, u = 0;
			if (i != 0) {
				t1 = ((ap[i + 1].h - ap[i].h) * ap[i + 1].d + (ap[i - 1].h - ap[i].h) * ap[i].d
						+ calX2(ap[i + 1].h - ap[i].h) * ap[i].d) * DeltaT / DeltaX / DeltaX;
				u = calX2(ap[i + 1].h - ap[i].h) * ap[i].d * DeltaT / DeltaX / DeltaX;
			} else
				t1 = ((ap[i + 1].h - ap[i].h) * ap[i + 1].d + calX2(ap[i + 1].h - ap[i].h) * ap[i].d) * DeltaT / DeltaX
						/ DeltaX;
			double t2 = current_v * DeltaT / DeltaX * (ap[i + 1].h - ap[i].h);
			double t3 = -ap[i].fout * ps.getPosNext() * DeltaT;
			ap[i].p1 = t1;
			ap[i].p2 = t2;
			ap[i].p3 = t3;
			ap[i].at = t1 + t2 + t3;// t1 + t2 + t3;// t1 + t2 + t3;
			p1.get(i).t1++;
			p1.get(i).t2 += Math.abs(t1);
			p2.get(i).t1++;
			p2.get(i).t2 += Math.abs(t2);
			p3.get(i).t1++;
			p3.get(i).t2 += Math.abs(t3);
			p4.get(i).t1++;
			p4.get(i).t2 += Math.abs(u);
		}
	}

	void Next(int index) {
		double g[] = new double[ap.length];
		for (int i = 0; i < ap.length - 1; i++)
			g[i] = ap[i].n;
		for (int i = 0; i < ap.length - 1; i++) {
			double tmp = ap[i].n;
			double th = ap[i].h;
			ap[i].h = ap[i].at + Math.log(Math.exp(ap[i].h) + ap[i].bt);
			ap[i].n = Math.exp(ap[i].h);

			ltk.get(i).t1++;
			ltk.get(i).t2 += ap[i].n;
			ltc.get(i).add(current_v, ap[i].n - tmp);
			for (int j = 0; j < NUM_COR; j++) {
				int num = 1 << j;
				double delt_v = velocity[index % velocity.length]
						- velocity[(index + velocity.length - TNUM * num) % velocity.length];
				double delt_n = ap[i].n - ord[(index + ord.length - TNUM * num) % ord.length][i];
				ltc2.get(j).get(i).add(delt_v / (TNUM * num), delt_n);
			}
			// output_d.write(ap[i].n + " ");
		}
		// output_d.write("\n");
		// System.out.println(index + "\t" + ap[10].n + "\t" + ap[11].n);
		// output_d.write("\n");

		BINDEX = index % ord.length;// (BINDEX + 1) % ord.length;
		for (int i = 0; i < NN; i++)
			ord[BINDEX][i] = ap[i].n;
		int p = (BINDEX + ord.length - TNUM) % ord.length;

		ap[ap.length - 1].n = ap[ap.length - 2].n;
		ap[ap.length - 1].h = ap[ap.length - 2].h;
	}

	void AP2() {
		for (int i = 0; i < ap2.length - 1; i++) {
			ap2[i].bt = ps.getPosNext() * DeltaT * ap2[i].fin;
		}
		for (int i = 0; i < ap2.length - 1; i++) {
			double t1 = 0;
			if (i != 0)
				t1 = ((ap2[i + 1].h - ap2[i].h) * ap2[i + 1].d + (ap2[i - 1].h - ap2[i].h) * ap2[i].d
						+ calX2(ap[i + 1].h - ap[i].h) * ap2[i].d) * DeltaT / DeltaX / DeltaX;
			else
				t1 = ((ap2[i + 1].h - ap2[i].h) * ap2[i + 1].d + calX2(ap[i + 1].h - ap[i].h) * ap2[i].d) * DeltaT
						/ DeltaX / DeltaX;

			double t2 = current_v * DeltaT / DeltaX * (ap2[i + 1].h - ap2[i].h);
			double t3 = -ap2[i].fout * ps.getPosNext() * DeltaT;
			ap2[i].at = t1 - t2 + t3;// t1 - t2 + t3;
		}

		for (int i = 0; i < ap2.length - 1; i++) {
			ap2[i].h = ap2[i].at + Math.log(Math.exp(ap2[i].h) + ap2[i].bt);
			ap2[i].n = Math.exp(ap2[i].h);
		}

		ap2[ap2.length - 1].n = ap2[ap2.length - 2].n;
		ap2[ap2.length - 1].h = ap2[ap2.length - 2].h;
	}

	double calX2(double val) {

		return (1 - Math.exp(-val * val * this.CC)) / this.CC;
	}

	void output(String file) {
		for (int i = 0; i < ap.length - 1; i++) {
			pdf1.calPdf(lld.get(i));
			pdf2.calLogPdf(lld.get(i));
		}
	}

	void outputExampleN() {
		for (int i = 0; i < lln.size(); i++) {
			for (int j = 0; j < lln.get(i).size(); j++)
				output_en.write(lln.get(i).get(j) + " ");
			output_en.write("\n");
		}
	}

	void outputTheoryExpect() {
		for (int i = 0; i < ap.length - 1; i++)
			output_expect.write(ap[i].fin + " ");
		output_expect.write("\n");
		for (int i = 0; i < ap.length - 1; i++)
			output_expect.write((ap[i].fout + 1e-8) + " ");
		output_expect.write("\n");
		for (int i = 0; i < ap.length - 1; i++) {
			output_expect.write(ap[i].fin / (ap[i].fout + 1e-8) + " ");
		}
		output_expect.write("\n");
	}

	void outputExpect() {
		for (int i = 0; i < ltk.size() - 1; i++) {
			output_expect.write(ltk.get(i).t2 / (ltk.get(i).t1 + 1e-8) + " ");
		}
		output_expect.write("\n");
		for (int i = 0; i < lld.size() - 1; i++)
			output_expect.write(Uti.Uti.calAverage(lld.get(i)) + " ");
		output_expect.write("\n");
	}

	void outputP() {
		for (int i = 0; i < p1.size() - 1; i++)
			output_p.write(p1.get(i).t2 / (p1.get(i).t1 + 1e-8) + " ");
		output_p.write("\n");

		for (int i = 0; i < p2.size() - 1; i++)
			output_p.write(p2.get(i).t2 / (p2.get(i).t1 + 1e-8) + " ");
		output_p.write("\n");

		for (int i = 0; i < p3.size() - 1; i++)
			output_p.write(p3.get(i).t2 / (p3.get(i).t1 + 1e-8) + " ");
		output_p.write("\n");

		for (int i = 0; i < p4.size() - 1; i++)
			output_p.write(p4.get(i).t2 / (p4.get(i).t1 + 1e-8) + " ");
		output_p.write("\n");

		for (int i = 0; i < r1.size() - 1; i++)
			output_p.write(r1.get(i).t2 / (r1.get(i).t1 + 1e-8) + " ");
		output_p.write("\n");

		for (int i = 0; i < r2.size() - 1; i++)
			output_p.write(r2.get(i).t2 / (r2.get(i).t1 + 1e-8) + " ");
		output_p.write("\n");

	}

	void outputCorPos() {
		for (int i = 0; i < mt.length - 1; i++) {
			for (int j = 0; j < mt.length - 1; j++)
				output_corpos.write(mt[i][j].calCorrelation() + " ");
			output_corpos.write("\n");
		}
		for (int i = 0; i < mt2.length - 1; i++) {
			for (int j = 0; j < mt2.length - 1; j++)
				output_corpos.write(mt2[i][j].calCorrelation() + " ");
			output_corpos.write("\n");
		}
	}

	void outputCor() {
		for (int i = 0; i < ltc.size() - 1; i++)
			output_cor.write(ltc.get(i).calCorrelation() + " ");
		output_cor.write("\n");
		for (int j = 0; j < NUM_COR; j++) {
			for (int i = 0; i < ltc2.get(j).size() - 1; i++)
				output_cor.write(ltc2.get(j).get(i).calCorrelation() + " ");
			output_cor.write("\n");
		}
	}

	void outputN() {
		for (int i = 0; i < llt.length; i++) {
			for (int j = 0; j < llt[i].length; j++)
				output_n.write(llt[i][j].t2 / (llt[i][j].t1 + 1e-8) + " ");
			output_n.write("\n");
		}
		for (int i = 0; i < llt.length; i++) {
			for (int j = 0; j < llt[i].length; j++)
				output_n.write((llt[i][j].t1) + " ");
			output_n.write("\n");
		}
	}

	void outputV() {
		for (int i = 0; i < lv.size(); i++)
			output_v.write(lv.get(i) + " ");
		output_v.write("\n");
	}

	void outputPdfV() {
		TK tmp[] = pdf_v.calLogPdf2(lv);
		CalK ck = new CalK();
		for (int i = 0; i < tmp.length; i++)
			if (tmp[i].t2 > 1 && tmp[i].t2 < 100)
				ck.add(Math.log10(tmp[i].t1), Math.log10(tmp[i].t2));
		ck.cal();
		pdf_v.output.write(ck.k + " " + ck.b + " " + ck.n + "\n");
		pdk = ck.k;
		pdb = ck.b;
	}

	void outputCumPdf() {
		List<Double> ld = new ArrayList<Double>();
		PdfOutput pdf = new PdfOutput();
		pdf.setParam(120, -6, 0);
		pdf.openFile("out.txt");
		pdf.calCumLog(lv, ld);
		pdf.calCumLog(lv);
		pdf.closeFile();
	}

	void outputVar() {
		// for (int i = 0; i < lvar.size(); i++) {
		// for (int j = 0; j < lvar.get(i).size(); j++)
		// output_var.write(lvar.get(i).get(j).cal2() + " ");
		// output_var.write("\n");
		// }
	}
}
