package particleRealDM;

import java.awt.Color;
import java.util.Random;

import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display.Circle;
import org.opensourcephysics.display.Trail;
import org.opensourcephysics.frames.DisplayFrame;
import org.opensourcephysics.frames.PlotFrame;
/**
 * UNIMPORTANT CLASS: JUST A DUPLICATE OF THE DIAGONALCLASHTESTDM CLASH: BUT THIS SHOWS THE DYNAMIC AMONG TWO BALLS INSTEAD (JUST SHIFTED A FEW VARIABLES (BALL NUMBER, SPEED, TOOK AWAY THE RANDOMNESS,ETC)
 * @author student
 *
 */
public class ClashRealTestDM extends AbstractSimulation{
	int molecules_num; 						//number of balls launched for each given velocity
	ParticleDM [] molecules; 				//the air molecules	
	PlotFrame disp = new PlotFrame("x", "y", "Diagonal Clash Test Display");	//display frame
	DisplayFrame dist = new DisplayFrame("Layer Number (in order of increasing height)", "Average Number of Particles", "Particle Distribution Cumulative Average per Layer");
	DisplayFrame instantaneous_dist = new DisplayFrame("Layer Number (in order of increasing height)", "Current Number of Particles", "Particle Distribution Instantaneous per Layer");

	//PlotFrame vt = new PlotFrame("time (seconds)", "velocity (m/s)","Velocity-Time Graph"); //velocity-time graph for single molecule
	double interval = .02;			//small time interval used to approximate the position of balls
	Circle yval = new Circle(); //mimicks y value of the first particle
	Circle xval = new Circle(); //mimicks x value of the single particle (sliding along x axis)
	double box_dims[] = new double [2]; //stores dimensions of box in which particle exist
	Random ran = new Random(); //random generator
	int ds_counter = 0; //do step counter. First do step is zero, each next one is labeled as the number one more than the previous
	FluidDM air;
	double threshold; //the distance between point-particles for there to exist a collision 


	CollisionDM romance;
	ParticleDistributionDM expanse;

	boolean [] righthits; //parallel array for molecules, keeps track of do step number of all collisions with left or right  hand wall (shouldn't have the same one again so soon). counts the right and left wall as same
	boolean [] lefthits; //parallel array for molecules, keeps track of do step number for all collisions with bottom  or top wall (counts top or bottom as same)
	boolean [] tophits;
	boolean [] downhits;
	boolean [][]	par_hits;  //unfortunately, a little less than half of these storage spaces will not be necessary due to symmetry ([m][n] would store same data as [n][m])

	double[] density_distribution;
	/**
	 * Initializes a whole batch of balls
	 * @param anglerange range of initial angles that each of these balls will possess
	 * @param anglemin minimum angle of this range
	 */
	public void initializeBalls (double anglerange, double anglemin){
		//for each ball in array
		double angle;	//angle
		boolean independence;
		//for every molecule
		for (int i = 0; i<molecules.length; i++){
			independence = false;
			double x;
			double y;
			if (i == 0){
				 x = 10;
				 y = -5;
				angle = 5*Math.PI/6;
			}
			else{
				 x = 10;
				 y = 5;
				angle = -5*Math.PI/6;
			}
			double speed = control.getDouble("Speed Max"); //generating random speed of this particle
			
			//initialize each ball: double ixpos, double iypos, double yai, double xai, double angle_init, double speed_init, double radius
			molecules[i] = new ParticleDM(x,	y,	 control.getDouble("Universal Y Acceleration"),		0, 	angle, speed, control.getInt("Radius")); 	//creates a molecule;	
		}
	}
	/**
	 * function meant as supplement to certain bug correction in ParticleInteraction.clash method. When a clash is detected, then the particle changes velocity and moves again, the same clash may still indeed be detected. This would occur because clash's are detected based on whether two entities are within a 'threshold' distance of each other (for which there is an double variable 'threshold'). Particles often stay within the threshold after the initial collision, and therefore have the clash function apply again, which causes errors with their velocities. In order to resolve this bug, there is are three arrays established to keep track of when (in which do step)  particles' collide with each other and the two walls. This function sets equal the returned arrays to the arrays in this class.
	 */
	public void initializeBugData(){
		righthits = romance.getRighthits(); //sets array of right and left wall interactions equal to the possibly modified version from the ParticleInteraction 'romance'
		lefthits = romance.getLefthits(); //'' 	'' 	'' 	top   and bottom '''' '' '' '' ' '' ' '' ''
		tophits = romance.getTophits();
		downhits = romance.getDownhits();
		par_hits = romance.getPar_hits(); // '''''''' particle-particle '' '' '' '' '' '' '' '' 
	}

	protected void doStep() {
		ds_counter ++; //keep track of do step number  ( consecutive positive integers)
		//update romance and expanse (doesn't seem to be a necessary step, I think I set the original as pointers to molecules, not as the value, although I am not sure)
		romance.molecules = molecules;
		expanse.molecules = molecules;

		//for each particle
		for(int n = 0; n<molecules.length; n++)
		{
			molecules = romance.clashReal(threshold, n, box_dims);
			if((romance.thereisaclash == true)){
				molecules[n].components(air, interval, control.getDouble("Universal Y Acceleration"), false); //establishes x and y components for the particle
			}
			else{
				molecules[n].components(air, interval, control.getDouble("Universal Y Acceleration"), true); //establishes x and y components for the particle
			}


			this.initializeBugData();
			molecules[n].setXY(molecules[n].getXpos(), molecules[n].getYpos());
			disp.addDrawable(molecules[n]);
			//Velocity-Time Graph for Molecule 0
			if (n==0){
				//	vt.append(0, counter*interval, molecules[0].getV()[1]);
			}

			yval.setXY(-box_dims[0]/2, molecules[n].getY());
			xval.setXY(molecules[n].getX(), -box_dims[1]/2);
			disp.addDrawable(yval);
			disp.addDrawable(xval);
		}
		expanse.averageDensityDistribution(ds_counter); //recalculates average density distribution
		//	expanse.printDensDistribution(); //prints density distribution
		expanse.plotDistribution(expanse.density_distribution, ds_counter, dist, 50); //plots distribution
		expanse.plotDistribution(expanse.densityDistribution(), ds_counter, instantaneous_dist, 4);
		
		
	}
	/**
	 * initializes necessary variables for simulation and formalities for the plotframe
	 */
	public void initialize() {
		air = new FluidDM(0); //establishing the resistance constant of the fluid in which the particles are moving (because the particles here compose the fluid itself, the "fluid" that they are in has zero resistance (there is no fluid in which they are in)
		//	air.airconst = control.getDouble("Air Resistance Constant"); //air resistance const.
		threshold = control.getDouble("Threshold"); //sets threshold variable equal to input from control panel
		molecules = new ParticleDM [control.getInt("Number of Particles")]; //initializing array of molecules
		interval = control.getDouble("Interval");//setting small time interval approximation
		box_dims[0] = control.getDouble("Box Width");
		box_dims[1] = control.getDouble("Box Height");
		density_distribution = new double [control.getInt("Density Distribution Layers")];

		initializeBalls(Math.PI/2 -.01, .01); //initialize all aspects of the molecules: their speeds, positions, etc.

		disp.setVisible(true); //set display frame as visible
		disp.setDefaultCloseOperation(3); //weird thing you just gotta do
		//		

		disp.clearData(); //wipe all data
		disp.clearDrawables(); //wipe all drawables
		dist.clearDrawables();
		disp.setPreferredMinMax(-.75*box_dims[1], .75*box_dims[1], -.75*box_dims[1], .75*box_dims[1]);

		//Draw box for molecules to bounce around in
		Trail walls = new Trail();
		walls.addPoint(-box_dims[0]/2, -box_dims[1]/2);
		walls.addPoint(box_dims[0]/2, -box_dims[1]/2);
		walls.addPoint(box_dims[0]/2, box_dims[1]/2);
		walls.addPoint(-box_dims[0]/2, box_dims[1]/2);
		walls.closeTrail();

		//for all molecules
		for (int i = 0; i<molecules.length; i++){
			disp.addDrawable(molecules[i]); //display each
		}
		disp.addDrawable(walls);

		//initialize these bug preventative record keeping of whether a different types of collisions have happened in the immediate past of any molecule that the program creates
		righthits = new boolean [molecules.length]; //parallel array for molecules, keeps track of do step number of all collisions with left or right  hand wall (shouldn't have the same one again so soon). counts the right and left wall as same
		lefthits = new boolean [molecules.length]; //parallel array for molecules, keeps track of do step number for all collisions with bottom  or top wall (counts top or bottom as same)
		tophits = new boolean [molecules.length];
		downhits = new boolean [molecules.length];
		par_hits = new boolean [molecules.length][molecules.length]; //each molecule has a 'romance' with each other molecule when measuring their collisions with each other

		romance = new CollisionDM(molecules, par_hits, righthits, lefthits, tophits, downhits); //intializes this temp particle interaction
		expanse = new ParticleDistributionDM(molecules, control.getInt("Density Distribution Layers"), box_dims[1]/control.getInt("Density Distribution Layers"), density_distribution);
		dist.setPreferredMinMax(0, expanse.layers + 1, 0, molecules.length*7/control.getInt("Density Distribution Layers"));
		dist.setVisible(true);
		dist.setDefaultCloseOperation(3);
		instantaneous_dist.setVisible(true);
		instantaneous_dist.setDefaultCloseOperation(3);
		instantaneous_dist.setPreferredMinMax(0, expanse.layers, 0, molecules.length*7/control.getInt("Density Distribution Layers"));
		//initialize density distribution
		expanse.density_distribution = expanse.densityDistribution();

		yval.color = Color.BLUE;
		xval.color = Color.BLUE;
	}
	/**
	 * resets all control panel values
	 */
	public void reset() {
		control.setValue("Interval",interval); //default interval width for particles' trajectories
		control.setValue("Radius", 20); //default radius size
		control.setValue("Air Resistance Constant", 0); //air resistance constant 0
		control.setValue("X Position", 0); // x position is zero
		control.setValue("Y Position", 0); //y position is zero
		control.setValue("Universal Y Acceleration", 0); //acceleration default is zero
		control.setValue("Threshold", 2); //threshold for collisions default is .5
		control.setValue("Number of Particles", 2); //2 particles
		control.setValue("Speed Max", 30); //maximum speed particles can have; default is 20m/s
		control.setValue("Box Width", 30);
		control.setValue("Box Height", 30);
		control.setValue("Density Distribution Layers", control.getInt("Box Height"));
	}
	/**
	 * launches the simulation
	 * @param args I don't know what this does
	 */
	public static void main(String[] args) {
		SimulationControl.createApp(new ClashRealTestDM());
	}

}


