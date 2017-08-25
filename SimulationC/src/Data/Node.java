package Data;

import java.util.ArrayList;
import java.util.List;

public class Node {
	public long time;
	List<TK> data = new ArrayList<TK>();
	List<TK> ask = new ArrayList<TK>();
	List<TK> bid = new ArrayList<TK>();
	public int index;
	public double t1, t2, min1, max1, min2, max2;
	public TK trade = null;
	String gt[];

	Node() {
	}

	public Node(String[] pt, int index) {
		try {
			this.gt = pt;
			this.index = index;
			time = Long.valueOf(pt[0]).longValue();
			List<TK> tmp1 = new ArrayList<TK>();
			List<TK> tmp2 = new ArrayList<TK>();
			data.clear();
			for (int i = 2; i < pt.length; i += 4)
				tmp1.add(new TK(Double.valueOf(pt[i]), Double.valueOf(pt[i + 1])));

			for (int i = 4; i < pt.length; i += 4)
				tmp2.add(new TK(Double.valueOf(pt[i]), Double.valueOf(pt[i + 1])));
			for (int i = tmp2.size() - 1; i >= 0; i--)
				bid.add(tmp2.get(i));
			min1 = tmp2.get(tmp2.size() - 1).t1;
			max1 = tmp2.get(0).t1;
			min2 = tmp1.get(tmp2.size() - 1).t1;
			max2 = tmp1.get(0).t1;
			t1 = tmp2.get(0).t1;
			t2 = tmp1.get(tmp1.size() - 1).t1;
			for (int i = tmp1.size() - 1; i >= 0; i--)
				ask.add(tmp1.get(i));
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println(index + " " + pt.length);
		}
	}

	public List<TK> getLst(int index) {
		if (index == -1)
			return bid;
		return ask;
	}

	public double getNumber(double val, int type) {
		if (type == -1)
			return getNumber(bid, val);
		return getNumber(ask, val);
	}

	double getNumber(List<TK> data, double val) {
		int left = 0;
		int right = data.size() - 1;
		int mid;
		while (left <= right) {
			mid = (left + right) / 2;
			if (Math.abs(data.get(mid).t1 - val) < 0.0001)
				return data.get(mid).t2;
			else if (data.get(mid).t1 > val)
				right = mid - 1;
			else
				left = mid + 1;
		}
		return 0;
	}

	public String toString() {
		String str = "";
		str = "[" + t1 + ", " + t2 + "], ";
		for (int i = 0; i < data.size(); i++)
			str = str + "(" + data.get(i).t1 + ", " + data.get(i).t2 + "), ";
		return str;
	}
}
