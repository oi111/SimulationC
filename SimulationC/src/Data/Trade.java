package Data;

public class Trade {
	public long date, date_ms, tid;
	public int type;
	public double price, amount;
	public int lprice;

	public Trade(long date, long date_ms, long tid, int type, double price, double amount) {
		this.date = date;
		this.date_ms = date_ms;
		this.tid = tid;
		this.type = type;
		this.price = price;
		this.amount = amount;
		this.lprice = (int) ((price + 1e-6) * 100);
	}
}
