package particleRealDM;

public class FluidDM {
	//variables
	double viscosity;
	double pressure;
	double density;
	double temperature;
	double poten_energy;
	double airconst;
	double [] wind;
	
	public FluidDM(double v, double p, double d, double t, double pe, double a, double [] w){
		this.viscosity = v;
		this.pressure = p;
		this.density = d;
		this.temperature = t;
		this.poten_energy = pe;
		this.airconst = a;
		this.wind = w;
	}
	public FluidDM(double airconst){
		this.airconst = airconst;
		this.viscosity = 0;
		this.pressure = 0;
		this.density = 0;
		this.temperature = 0;
		this.poten_energy = 0;
		this.airconst = 0;
		//this.wind = new double [2];
	}
	
	public double getViscosity() {
		return viscosity;
	}
	public void setViscosity(double viscosity) {
		this.viscosity = viscosity;
	}
	public double getPressure() {
		return pressure;
	}
	public void setPressure(double pressure) {
		this.pressure = pressure;
	}
	public double getDensity() {
		return density;
	}
	public void setDensity(double density) {
		this.density = density;
	}
	public double getTemperature() {
		return temperature;
	}
	public void setTemperature(double temperature) {
		this.temperature = temperature;
	}
	public double getPoten_energy() {
		return poten_energy;
	}
	public void setPoten_energy(double poten_energy) {
		this.poten_energy = poten_energy;
	}
	public double getAirconst() {
		return airconst;
	}
	public void setAirconst(double airconst) {
		this.airconst = airconst;
	}
	
	
}
