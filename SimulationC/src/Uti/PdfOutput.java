package Uti;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import Data.TK;

public class PdfOutput {
	int len = 10000;
	double xmin = 0, xmax = 1000;
	double delt;
	public OutputFile output = new OutputFile();

	public PdfOutput() {
		openFile("PdfOutput.txt");
	}

	void init() {

	}

	public void setParam(int len, double xmin, double xmax) {
		this.len = len;
		this.xmax = xmax;
		this.xmin = xmin;
		this.delt = (xmax - xmin + 0.0) / len;
	}

	public void calPdf(List<Double> ld) {
		double a[] = new double[len], tot = 0;
		for (int i = 0; i < len; i++)
			a[i] = 0;
		for (int i = 0; i < ld.size(); i++) {
			tot++;
			if ((int) ((ld.get(i).doubleValue() - xmin) / delt) >= len
					|| (int) ((ld.get(i).doubleValue() - xmin) / delt) < 0)
				continue;
			a[(int) ((ld.get(i).doubleValue() - xmin) / delt)]++;
		}
		for (int i = 0; i < len; i++) {
			output.write((xmin + i * delt) + " ");
		}
		output.write("\n");
		for (int i = 0; i < len; i++) {
			output.write(a[i] / (tot + 0.0) / delt + " ");
		}
		output.write("\n");
	}

	public void calAutoPdf(List<Double> ld) {
		double a[] = new double[len], tot = 0;
		List<Double> tmp = new ArrayList<Double>();
		tmp.clear();
		if (ld.size() != 0) {
			double val = Uti.calSigma(ld);

			for (int i = 0; i < ld.size(); i++)
				tmp.add(ld.get(i) / val);
			Collections.sort(tmp);
			xmin = tmp.get((int) (tmp.size() * 0.01));
			xmax = tmp.get((int) (tmp.size() * 0.99));
			delt = (xmax - xmin + 0.0) / len;
		}
		for (int i = 0; i < len; i++)
			a[i] = 0;
		for (int i = (int) (tmp.size() * 0.01); i < (int) (tmp.size() * 0.99); i++) {
			tot++;
			if ((int) ((tmp.get(i).doubleValue() - xmin) / delt) >= len
					|| (int) ((tmp.get(i).doubleValue() - xmin) / delt) < 0)
				continue;
			a[(int) ((tmp.get(i).doubleValue() - xmin) / delt)]++;
		}
		for (int i = 0; i < len; i++) {
			output.write((xmin + i * delt) + " ");
		}
		output.write("\n");
		for (int i = 0; i < len; i++) {
			output.write(a[i] / (tot + 0.0) / delt + " ");
		}
		output.write("\n");
	}

	public void calNormalPdf(List<Double> ld) {
		double a[] = new double[len], tot = 0;
		double val = Uti.calSigma(ld);
		for (int i = 0; i < len; i++)
			a[i] = 0;

		for (int i = 0; i < ld.size(); i++) {
			tot++;
			if ((int) ((ld.get(i).doubleValue() / val - xmin) / delt) >= len
					|| (int) ((ld.get(i).doubleValue() / val - xmin) / delt) < 0)
				continue;
			a[(int) ((ld.get(i).doubleValue() / val - xmin) / delt)]++;
		}
		for (int i = 0; i < len; i++) {
			output.write((xmin + i * delt) + " ");
		}
		output.write("\n");
		for (int i = 0; i < len; i++) {
			output.write(a[i] / (tot + 0.0) / delt + " ");
		}
		output.write("\n");
	}

	public void calLogPdfToLinear(List<Double> ld) {
		double a[] = new double[len], tot = 0;
		for (int i = 0; i < len; i++)
			a[i] = 0;
		for (int i = 0; i < ld.size(); i++) {
			double tmp = Math.log(Math.abs(ld.get(i).doubleValue()));
			if ((int) ((tmp - xmin) / delt) >= len || (int) ((tmp - xmin) / delt) < 0)
				continue;
			a[(int) ((tmp - xmin) / delt)]++;
			tot++;
		}
		for (int i = 0; i < len; i++) {
			// System.out.print((xmin + i * delt) + " ");
			output.write(Math.exp(xmin + i * delt) + " ");
		}
		output.write("\n");
		// System.out.println();
		for (int i = 0; i < len; i++) {
			// System.out.print(a[i] / (tot + 0.0) + " ");
			output.write(a[i] / (tot + 0.0) / Math.exp(xmin + i * delt) + " ");
		}
		output.write("\n");
		// System.out.println();
	}

	public void calLogPdf(List<Double> ld) {
		double a[] = new double[len], b[] = new double[len], tot = 0;
		for (int i = 0; i < len; i++)
			b[i] = a[i] = 0;
		for (int i = 0; i < ld.size(); i++) {
			tot++;
			if ((int) ((Math.log10(Math.abs(ld.get(i).doubleValue())) - xmin) / delt) >= len
					|| (int) ((Math.log10(Math.abs(ld.get(i).doubleValue())) - xmin) / delt) < 0)
				continue;
			if (ld.get(i).doubleValue() > 0)
				a[(int) ((Math.log10(ld.get(i).doubleValue()) - xmin) / delt)]++;
			if (ld.get(i).doubleValue() < 0)
				b[(int) ((Math.log10(-ld.get(i).doubleValue()) - xmin) / delt)]++;
		}
		for (int i = 0; i < len; i++)
			output.write(Math.pow(10, (xmin + (i + 1) * delt)) + " ");
		output.write("\n");
		for (int i = 0; i < len; i++)
			output.write(a[i] / (tot + 0.0) / (Math.pow(10, (xmin + (i + 1) * delt)) - Math.pow(10, (xmin + i * delt)))
					+ " ");
		output.write("\n");
		for (int i = 0; i < len; i++)
			output.write(b[i] / (tot + 0.0) / (Math.pow(10, (xmin + (i + 1) * delt)) - Math.pow(10, (xmin + i * delt)))
					+ " ");
		output.write("\n");
	}

	public TK[] calLogPdf2(List<Double> ld) {
		double a[] = new double[len], b[] = new double[len], tot = 0;
		for (int i = 0; i < len; i++)
			b[i] = a[i] = 0;
		for (int i = 0; i < ld.size(); i++) {
			tot++;
			if ((int) ((Math.log10(Math.abs(ld.get(i).doubleValue())) - xmin) / delt) >= len
					|| (int) ((Math.log10(Math.abs(ld.get(i).doubleValue())) - xmin) / delt) < 0)
				continue;
			if (ld.get(i).doubleValue() > 0)
				a[(int) ((Math.log10(ld.get(i).doubleValue()) - xmin) / delt)]++;
			if (ld.get(i).doubleValue() < 0)
				b[(int) ((Math.log10(-ld.get(i).doubleValue()) - xmin) / delt)]++;
		}
		for (int i = 0; i < len; i++)
			output.write(Math.pow(10, (xmin + (i + 1) * delt)) + " ");
		output.write("\n");
		for (int i = 0; i < len; i++)
			output.write(a[i] / (tot + 0.0) / (Math.pow(10, (xmin + (i + 1) * delt)) - Math.pow(10, (xmin + i * delt)))
					+ " ");
		output.write("\n");
		TK ret[] = new TK[len];
		for (int i = 0; i < len; i++) {
			ret[i] = new TK(0, 0);
			ret[i].t1 = Math.pow(10, (xmin + (i + 1) * delt));
			ret[i].t2 = a[i] / (tot + 0.0) / (Math.pow(10, (xmin + (i + 1) * delt)) - Math.pow(10, (xmin + i * delt)));
		}
		for (int i = 0; i < len; i++)
			output.write(b[i] / (tot + 0.0) / (Math.pow(10, (xmin + (i + 1) * delt)) - Math.pow(10, (xmin + i * delt)))
					+ " ");
		output.write("\n");
		return ret;
	}

	public void calCumLog(List<Double> ld) {
		double a[] = new double[len], b[] = new double[len], tot = 0;
		for (int i = 0; i < len; i++)
			b[i] = a[i] = 0;
		for (int i = 0; i < ld.size(); i++) {
			tot++;
			if ((int) ((Math.log10(Math.abs(ld.get(i).doubleValue())) - xmin) / delt) >= len
					|| (int) ((Math.log10(Math.abs(ld.get(i).doubleValue())) - xmin) / delt) < 0)
				continue;
			if (ld.get(i).doubleValue() > 0)
				a[(int) ((Math.log10(ld.get(i).doubleValue()) - xmin) / delt)]++;
			// if (ld.get(i).doubleValue() < 0)
			// b[(int) ((Math.log10(-ld.get(i).doubleValue()) - xmin) /
			// delt)]++;
		}
		for (int i = len - 2; i >= 0; i--)
			a[i] = a[i] + a[i + 1];
		for (int i = 0; i < len; i++)
			output.write(Math.pow(10, (xmin + (i + 1) * delt)) + " ");
		output.write("\n");
		for (int i = 0; i < len; i++)
			output.write(a[i] / (tot + 0.0) + " ");
		output.write("\n");
	}

	public void calCumLog(List<Double> ld, List<Double> ret) {
		double a[] = new double[len], b[] = new double[len], tot = 0;
		for (int i = 0; i < len; i++)
			b[i] = a[i] = 0;
		for (int i = 0; i < ld.size(); i++) {
			tot++;
			if ((int) ((Math.log10(Math.abs(ld.get(i).doubleValue())) - xmin) / delt) >= len
					|| (int) ((Math.log10(Math.abs(ld.get(i).doubleValue())) - xmin) / delt) < 0)
				continue;
			if (ld.get(i).doubleValue() > 0)
				a[(int) ((Math.log10(ld.get(i).doubleValue()) - xmin) / delt)]++;
			// if (ld.get(i).doubleValue() < 0)
			// b[(int) ((Math.log10(-ld.get(i).doubleValue()) - xmin) /
			// delt)]++;
		}
		for (int i = len - 2; i >= 0; i--)
			a[i] = a[i] + a[i + 1];
		ret.clear();
		for (int i = 0; i < len; i++)
			ret.add(a[i] / (tot + 0.0));
	}

	public void openFile(String file) {
		// closeFile();
		output.setFileName(file);
		output.openFile();
	}

	public void closeFile() {
		try {
			output.closeFile();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
