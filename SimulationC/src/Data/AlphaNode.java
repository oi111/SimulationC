package Data;

import java.util.ArrayList;
import java.util.List;

public class AlphaNode {

	List<Double> ld = new ArrayList<Double>();
	AlphaParamECF ecf = new AlphaParamECF();
	AlphaParamFlom flom = new AlphaParamFlom();
	AlphaParamLog log = new AlphaParamLog();

	public AlphaNode() {
	}

	public void add(int type, double val) {
		switch (type) {
		case 1:
			ecf.add(val);
			break;
		case 2:
			flom.add(val);
			break;
		case 3:
			log.add(val);
			break;
		default:
			break;
		}
	}

	public void cal(int type, double val) {
		switch (type) {
		case 1:
			ecf.cal();
			break;
		case 2:
			flom.cal();
			break;
		case 3:
			log.cal();
			break;
		default:
			break;
		}
	}

	public void getParam(int type) {

	}

}
