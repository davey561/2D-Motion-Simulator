package particleRealDM;

import java.util.ArrayList;

import org.opensourcephysics.display.Circle;
import org.opensourcephysics.frames.PlotFrame;
/**
 * This is the main Particle Object Class. All its fields are characteristics of a ball or particle (at least macroscopically). The methods besides setters and getters, deal with calculating recursively the particle's next position, and doing some math related to the particle's angle of movement that was helpful for calculations elsewhere.
 * @author Davey Morse
 *
 *Bugs
 *1. With angles:
 *	initially presuming I could use all positive angles 0 thru 2pi, but this didn't seem to be allowed by the arctan restrictions
 *	was initially using Math.atan to calculate angles from cartesian components, but range is restricted from -pi/2 to pi/2
 *	then started using Math.atan2 to calculate angles from cartesian components, assuming this would have range of all positive angles in unit circle (0 thru 2pi)
 *	after looking at documentation, realized that range is from -pi to pi, which is extended, but not exactly the way I'd wanted. So now, I will change the way I generate the angles, as to fit the nature of this function
 *
 */
public class ParticleDM extends Circle{

	//global variables of particles
	String name = "bobo"; //particle's name

	double [] vi = new double [2]; //Initial Velocity: first item is angle, second is magnitude (not used)
	double [] v = new double [2]; //Velocity as a function of time: first item is angle, second is magnitude

	double ixpos= 0; //initial x position of particle
	double xpos = 0; //variable x position
	double iypos = 0; //initial y position
	double ypos = 0;//variable y position

	double yvi = 0; //initial y velocity's magnitude (not used)
	double yv = 0; //velocity variable
	double xvi = 0;// initial x velocity (not used)
	double xv = 0; //x velocity
	
	//the accelerations simulated in the test classes are all universal (apply to all particles--> gravity), so individual accelerations for each particle aren't really necessary. Nonetheless, here they are
	double Xai = 0;	//initial x acceleration's magnitude
	double Xa = 0; //variable x acceleration
	double Yai = 0;	//initial y acceleration's magnitude
	double Ya = 0; //variable y acceleration

	double ti = 0; //initial place in time
	double t = 0; //time variable

	double mass = 0; //mass of particle
	double radius = 0; // particle is spherical, proportional to cube root of volume

	//METHODS

	//Constructor
	public ParticleDM(double ixpos, double iypos, double yai, double xai, double angle_init, double speed_init, int radius){
		this.setXY(ixpos, iypos);
		this.v[0] = angle_init;
		this.v[1] = speed_init;
		this.xpos = ixpos;
		this.ypos = iypos;

		recalcCartesianComponents(); //sets yv and xv to proper values
		this.Xai = xai;
		this.Yai = yai;
		pixRadius = radius; //prolly does nothing, need to learn how to change acutal radius
	}
	public ParticleDM(){
		this.ixpos= 0; //initial x position of particle
		this.xpos = 0; //variable x position

		this.iypos = 0; //initial y position
		this.ypos = 0;//variable y position

		this.yvi = 0; //initial y velocity's magnitude
		this.yv = 0; //velocity variable

		this.xvi = 0;// initial x velocity
		this.xv = 0; //x velocty

		this.Xai = 0;	//initial acceleration's magnitude
		this.Xa = 0; //variable acceleration

		this.Yai = 0;	//initial acceleration's magnitude
		this.Ya = 0; //variable acceleration

		this.ti = 0; //initial place in time
		this.t = 0; //time variable
	}

	
	//Setters and Getters
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public double getIxpos() {
		return ixpos;
	}
	public void setIxpos(double ixpos) {
		this.ixpos = ixpos;
	}
	public double getXpos() {
		return xpos;
	}
	public void setXpos(double xpos) {
		this.xpos = xpos;
	}
	public double getIypos() {
		return iypos;
	}
	public void setIypos(double iypos) {
		this.iypos = iypos;
	}
	public double getYpos() {
		return ypos;
	}
	public void setYpos(double ypos) {
		this.ypos = ypos;
	}

	public double getYv() {
		return yv;
	}
	public void setYv(double yv) {
		this.yv = yv;
		
		//after a change in the y component, fix the particle's magnitude and angle so that they are now also consistent
		recalcPolarComponents();
	}
	public double getXv() {
		return xv;
	}
	public void setXv(double xv) {
		this.xv = xv;
		
		//after a change in the y component, fix the particle's magnitude and angle so that they are now also consistent
		recalcPolarComponents();
	}
	public double getXai() {
		return Xai;
	}
	public void setXai(double xai) {
		this.Xai = xai;
	}
	public double getYai() {
		return Yai;
	}
	public void setYai(double ai) {
		this.Yai = ai;
	}
	public double getXa() {
		return Xa;
	}
	public void setXa(double xai) {
		this.Xa = xai;
	}
	public double getYa() {
		return Ya;
	}
	public void setYa(double ai) {
		this.Ya = ai;
	}
	public double getTi() {
		return ti;
	}
	public void setTi(double ti) {
		this.ti = ti;
	}
	public double getT() {
		return t;
	}
	public void setT(double t) {
		this.t = t;
	}

	public double[] getVi() {
		return vi;
	}

	public void setVi(double[] vi) {
		this.vi = vi;
	}

	public double [] getV(){
		return v;
	}

	public void setV(double[] v) {
		this.v = v;
		
		//recalculates the X and Y components of velocity, which much change accordingly
		recalcCartesianComponents();
	}

	/**
	 * Allows user to set the angle within also addressing the velocity's magnitude- just makes code slightly more concise in some places
	 * @param angle
	 */
	public void setAngle(double angle){
		this.setV( new double[] {angle, this.getV()[1]});
	}
	/**
	 * Allows user to set the velocity magnitude within also addressing the velocity's direction- just makes code slightly more concise in some places
	 * @param angle
	 */
	public void setVelocityMag(double mag){
		this.setV( new double[] {this.getV()[0], mag});
	}

	/**
	 * calculates the particle's direction of movement based on its x and y components of velocity
	 * @return current direction, based of cartesian components
	 */
	public double calcAngle() {
		return (Math.atan2(yv, xv));
	}
	public double getMass() {
		return mass;
	}

	public void setMass(double mass) {
		this.mass = mass;
	}

	public double getRadius() {
		return radius;
	}

	public void setRadius(double radius) {
		this.radius = radius;
	}

	/**
	 * Resets the Particle's x and y velocity components. This method is commonly used after the function's Polar Velocity Components have been changed (direction and magnitude), for the sake of consistency.
	 */
	public void recalcCartesianComponents(){
		this.xv = v[1]*Math.cos(v[0]);
		this.yv = v[1]*Math.sin(v[0]);
	}
	
	/**
	 * Find's x and y velocity components for a given velocity v, based on its velocity's current magnitude and direction
	 * @param v the given velocity, for which the x and y components are found
	 * @return
	 */
	public double[] findCartesianComponents(double [] v){
		double [] components = new double [2]; //first term is magnititude of the x component, second term is mag of the y component
		components[0] = v[1]*Math.cos(v[0]); //the x component of v
		components[1] = v[1]*Math.sin(v[0]); //the y component of v
		return components;
	}
	
	/**
	 * Adds two vectors (written in terms of their polar components), and returns the resultant's Cartesian Components (X and Y)
	 * @param v1
	 * @param v2
	 * @return
	 */
	public double [] addVectorsCartesian(double [] v1, double [] v2){
		double [] resultant = new double [2]; //first term is x component of resultant, second term is y component of the resultant
		resultant[0] = findCartesianComponents(v1)[0] + findCartesianComponents(v2)[0]; //finds x components of each vector, and adds them
		resultant[1] = findCartesianComponents(v1)[1] + findCartesianComponents(v2)[1]; //find y components of each vector, and adds them
		return resultant;
	}
	public void recalcPolarComponents(){
		//if vertical
		this.v[0] = Math.atan2(yv,xv); //set angle based on x and y components
		//magnitude
		this.v[1]   = Math.sqrt(Math.pow(yv, 2) + Math.pow((xv), 2)); //find magnitude based on x and y components (pythagorean theorem)
	}
	
	/**
	 * finds the polar components of the resultant vector when the two given are added
	 * @param v1 velocity one, given in polar form
	 * @param v2 velocity two, given in polar form
	 * @return
	 */
	public double [] findPolarComponents(double [] v1, double [] v2){
		double [] resultant = new double [2]; //first term is angle, second term is magnitude
		double[] cartesianr = addVectorsCartesian(v1, v2); //first adds them and finds cartesian result
		
		//translates to polar
		resultant[0] = Math.atan2(cartesianr[1],cartesianr[0]); //set angle based on x and y components
		//magnitude
		resultant[1]  = Math.sqrt(Math.pow(cartesianr[0], 2) + Math.pow(cartesianr[1], 2)); //sets magnitude
		return resultant;
	}
	/**
	 * Just a mathematical difference quotient
	 * @param a double
	 * @param b double 
	 * @param c double
	 * @param d double
	 * @return
	 */
	
	public double difQuot (double a, double b, double c, double d){
		return (a-b)/(c-d); //difference between a and b divided by difference between c and d
	}
	//NO FLUID
	//average velocity
	public double averageXVelocity (){
		return difQuot(x, ixpos, t, ti);
	}
	//average acceleration
	public double averageXAccel (){
		return difQuot(xv, xvi, t, ti);
	}
	/**
	 * Recursively recalculates calculation particle's position every time step
	 * @param air the fluid in which the particle is traveling
	 * @param interval the time step value
	 * @param yaccel a Universal Y acceleration
	 * @param calc_accel_and_velocity telling the program whether or not to recalculate the velocity and acceleration. The function only avoids this step when the velocity has just been recalculated by the CollisionDM class. Otherwise, acceleration and velocity are always recalculated
	 */
	//components
	public void components(FluidDM air, double interval, double yaccel, boolean calc_accel_and_velocity){
		this.setYai(yaccel);	
		if(calc_accel_and_velocity == true){
			//calculating Y velocity and acceleration, dependent on each other
			this.setYa(this.getYai() - 0*air.getAirconst()*this.getYv());	 //changes acceleration moleculesy -alpha*p.getYv(), which represents air resistance as a function of velocity
			this.setYv(this.getYv() + this.getYa()*interval); //changes velocity due to acceleration in given time interval
			//calculating X velocity and acceleration, dependent on each other
			this.setXa(this.getXai() - 0*air.getAirconst()*this.getXv());	 //changes acceleration by -alpha*p.getXv(), which represents air resistance as a function of velocity
			this.setXv(this.getXv() + this.getXa()*interval); // changes velocity due to acceleration in given time interval
		}
		
		this.setYpos(this.getYpos() + this.getYv()*interval); //changes y position due to velocity in given interval
		this.setXpos(this.getXv()*interval + this.getXpos()); //changes position due to velocity over given time interval
	}
}
