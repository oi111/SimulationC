package Uti;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;

public class InputFile {
	public File f;
	public FileInputStream fis;
	public InputStreamReader isw;
	public BufferedReader br;
	public String file_name;

	public void setFileName(String file) {
		file_name = file;
	}

	public void openFile() {
		try {
			f = new File(file_name);
			fis = new FileInputStream(f);
			isw = new InputStreamReader(fis, "utf-8");
			br = new BufferedReader(isw);
		} catch (Exception e) {
			System.out.println("err in openFile in InputFile Class");
		}
	}

	public String read() {
		try {
			return br.readLine();
		} catch (Exception e) {
			System.out.println("err in read in InputFile Class");
		}
		return "";
	}

	public void closeFile() {
		try {
			br.close();
			isw.close();
			fis.close();
		} catch (Exception e) {
			System.err.println("err in closeFile in OutputFile Class");
		}
	}
}
