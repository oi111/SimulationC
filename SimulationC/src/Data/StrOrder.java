package Data;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class StrOrder {
	public Long time;

	public ITK ask[] = new ITK[200];
	public ITK bid[] = new ITK[200];
	public int bidN[];
	public double mid;
	public int t1;

	public int t2;

	// public int min1;
	public String str1;

	public String str2;
	int bin;

	public StrOrder(String time, String ask, String bid) {
		this.time = Long.valueOf(time).longValue();
		init(ask, bid);
		str1 = bid;
		str2 = ask;
		this.bin = bin;
	}

	void init(String ask, String bid) {
		processOne(ask, this.ask);
		processOne(bid, this.bid);
		t1 = this.bid[0].t1;
		// min1 = this.ybid[this.ybid.length - 1].t1;
		t2 = this.ask[this.ask.length - 1].t1;
		mid = (t1 + t2) / 2.0;
	}

	void processOne(String str, ITK lt[]) {
		String pt[] = str.replace("[", "").replace("]", "").split(",");
		for (int i = 0; i < pt.length / 2; i++) {
			int index = (int) ((Double.valueOf(pt[2 * i]).doubleValue() + 1e-6) * 100);
			int volume = (int) ((Double.valueOf(pt[2 * i + 1]).doubleValue() + 1e-6) * 1000);
			lt[i] = new ITK(index, volume);
		}
	}

	public static int[] getBidDiff(StrOrder a, StrOrder b, int n, int bin) {
		double best_bid = Math.max(a.mid, b.mid);
		int pt[] = new int[n];
		// List<ITK> tmp = connectITK(a.bid, b.bid);
		for (int i = 0; i < pt.length; i++)
			pt[i] = 0;

		for (int i = 0; i < a.bid.length; i++) {
			int index = -(int) (Math.log(a.bid[i].t1 / best_bid) * 100000 / bin);
			if (index >= 0 && index < n)
				pt[index] += a.bid[i].t2;
			// if ((best_bid - a.bid[i].t1) >= 0 && best_bid - a.bid[i].t1 < n)
			// pt[best_bid - a.bid[i].t1] += a.bid[i].t2;
		}
		for (int i = 0; i < b.bid.length; i++) {
			int index = -(int) (Math.log(b.bid[i].t1 / best_bid) * 100000 / bin);
			if (index >= 0 && index < n)
				pt[index] -= b.bid[i].t2;
			// if (best_bid - b.bid[i].t1 >= 0 && best_bid - b.bid[i].t1 < n)
			// pt[best_bid - b.bid[i].t1] -= b.bid[i].t2;
		}
		return pt;
	}

	public static ITK[] getBidDiffITK(StrOrder a, StrOrder b, int n, int bin) {
		double best_bid = Math.max(a.mid, b.mid);
		ITK pt[] = new ITK[n];
		for (int i = 0; i < pt.length; i++)
			pt[i] = new ITK(0, 0);
		for (int i = 0; i < a.bid.length; i++) {
			int index = -(int) (Math.log(a.bid[i].t1 / best_bid) * 100000 / bin);
			if (index >= 0 && index < n)
				pt[index].t1 += a.bid[i].t2;
		}
		for (int i = 0; i < b.bid.length; i++) {
			int index = -(int) (Math.log(a.bid[i].t1 / best_bid) * 100000 / bin);
			if (index >= 0 && index < n) {
				pt[index].t1 -= b.bid[i].t2;
				pt[index].t2 += b.bid[i].t2;
			}
		}
		return pt;
	}

	public static List<ITK> connectITK(ITK a[], ITK b[]) {
		List<ITK> li = new ArrayList<ITK>();
		for (int i = 0; i < a.length; i++) {
			ITK tmp = new ITK(a[i].t1, a[i].t2);
			li.add(tmp);
		}
		for (int i = 0; i < b.length; i++) {
			ITK tmp = new ITK(b[i].t1, -b[i].t2);
			li.add(tmp);
		}
		Collections.sort(li, new Comparator<ITK>() {
			@Override
			public int compare(ITK arg0, ITK arg1) {
				// TODO Auto-generated method stub
				return arg0.t1 - arg1.t1;
			}
		});
		return li;
	}

}
