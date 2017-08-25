package Uti;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;

public class OutputFile {
	public File f;
	public FileOutputStream fos;
	public OutputStreamWriter osw;
	public BufferedWriter bw;
	public String file_name;

	public void setFileName(String file) {
		file_name = file;
	}

	public void openFile() {
		try {
			f = new File(file_name);
			fos = new FileOutputStream(f);
			// osw = new OutputStreamWriter(fos, "utf-8");
			osw = new OutputStreamWriter(fos);
			bw = new BufferedWriter(osw);
		} catch (Exception e) {
			System.out.println("err in openFile in OutputFile Class");
		}
	}

	public void write(String str) {
		try {
			// System.out.println(str);
			bw.write(str);
		} catch (Exception e) {
			System.out.println("err in write in OutputFile Class");
		}
	}

	public void closeFile() {
		try {
			bw.close();
			osw.close();
			fos.close();
		} catch (Exception e) {
			System.err.println("err in closeFile in OutputFile Class");
		}
	}
}
