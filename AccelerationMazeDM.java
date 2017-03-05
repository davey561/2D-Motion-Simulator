package particleRealDM;

import java.awt.Color;
import java.util.Random;

import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display.DrawableShape;
import org.opensourcephysics.frames.PlotFrame;

/*
 * Personal Project for 2D Motion Unit
 * 
 * Change to make: all speed, distance, etc. should depend on displacement and direction
 */
public class AccelerationMazeDM extends AbstractSimulation{
	int molecules_num; 						//number of balls launched for each given velocity
	ParticleDM [] molecules; 				//the air molecules	
	PlotFrame disp = new PlotFrame("x", "y", "Acceleration Maze Display");	//display frame
	//PlotFrame vt = new PlotFrame("time (seconds)", "velocity (m/s)","Velocity-Time Graph"); //velocity-time graph for single molecule
	double interval = .02;			//small time interval used to approximate the position of balls
	Random ran = new Random();
	double counter = 0;
	FluidDM air;
	double threshold; //the distance between point-particles for there to exist a collision 

	/**
	 * Initializes a whole batch of balls
	 * @param anglerange range of initial angles that each of these balls will possess
	 * @param anglemin minimum angle of this range
	 */
	public void initializeBalls (double anglerange, double anglemin){
		//for each ball in array
		double space = 1;
		double angle;
		for (int i = 0; i<molecules.length; i++){
			//			double angle = ran.nextDouble()*2*Math.PI;
			if (i == 0){
				angle = 0;
			}
			else{
				angle = Math.PI;
			}

			double speed = ran.nextDouble()*control.getDouble("Speed Max"); //generating random speed of this particle
			//double x = ran.nextDouble()*space;
			//double y = ran.nextDouble()*space;
			double x = -4 + 8*i; //xposition
			double y = 0; //yposition
			//initialize each ball: double ixpos, double iypos, double yai, double xai, double angle_init, double speed_init, double radius
			molecules[i] = new ParticleDM(x,	y,	 control.getDouble("Universal Y Acceleration"),		0, 	angle, speed, control.getInt("Radius")); 	//creates a molecule;	

		}
	}

	public void components(ParticleDM molecule){
		System.out.println("Components!");
		molecule.setYa(molecule.getYai() - air.getAirconst()*molecule.getYv());	 //changes acceleration moleculesy -alpha*p.getYv(), which represents air resistance as a function of velocity
		System.out.println("uno");
		molecule.setYv(molecule.getYv() + molecule.getYa()*interval); //changes velocity due to acceleration in given time interval
	System.out.println("dos");
		molecule.setYpos(molecule.getYpos() + molecule.getYv()*interval); //changes y position due to velocity in given interval
		
		System.out.println("Middle Components!");
		molecule.setXa(molecule.getXai() - air.getAirconst()*molecule.getXv());	 //changes acceleration by -alpha*p.getXv(), which represents air resistance as a function of velocity
		molecule.setXv(molecule.getXv() + molecule.getXa()*interval); // changes velocity due to acceleration in given time interval
		molecule.setXpos(molecule.getXv()*interval + molecule.getXpos()); //changes position due to velocity over given time interval
	}

	protected void doStep() {
		CollisionDM romance = new CollisionDM(molecules);
		//for each particle
		counter ++;
		for(int n = 0; n<molecules.length; n++)
		{
			romance.clash(threshold, n, 5, counter);
			components(molecules[n]); //establishes x and y components for the particle
			molecules[n].setXY(molecules[n].getXpos(), molecules[n].getYpos());
			disp.addDrawable(molecules[n]);
			//Velocity-Time Graph for Molecule 0
			if (n==0){
				//	vt.append(0, counter*interval, molecules[0].getV()[1]);
			}
		}

	}
	/**
	 * initializes necessary variables for simulation and formalities for the plotframe
	 */
	public void initialize() {
		air = new FluidDM(0);
		air.airconst = control.getDouble("Air Resistance Constant");
		threshold = control.getDouble("Threshold");
		molecules = new ParticleDM [control.getInt("Number of Particles")];
		initializeBalls(Math.PI/2 -.02, .01);
		interval = control.getDouble("Interval");//setting small time interval approximation
		disp.setVisible(true);
		disp.setDefaultCloseOperation(3);
		double x = 10;
		disp.setPreferredMinMax(-x, x, -x, x);

		disp.clearData();
		disp.clearDrawables();
		System.out.println("All initialized");
		for (int i = 0; i<molecules.length; i++){
			disp.addDrawable(molecules[i]);
		}

	}
	/**
	 * resets all control panel values
	 */
	public void reset() {
		control.setValue("Interval",interval); //default interval width for particles' trajectories
		control.setValue("Radius", 5); //default radius size
		control.setValue("Air Resistance Constant", 0);
		control.setValue("X Position", 0);
		control.setValue("Y Position", 0);
		control.setValue("Universal Y Acceleration", 0);
		control.setValue("Threshold", .5);
		control.setValue("Number of Particles", 2);
		control.setValue("Speed Max", 20);
	}
	/**
	 * launches the simulation
	 * @param args I don't know what this does
	 */
	public static void main(String[] args) {
		SimulationControl.createApp(new AccelerationMazeDM());
	}
}
