package particleRealDM;

import java.awt.Color;
import java.util.Random;

import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display.Circle;
import org.opensourcephysics.display.DrawableShape;
import org.opensourcephysics.display.Trail;
import org.opensourcephysics.frames.DisplayFrame;
import org.opensourcephysics.frames.PlotFrame;

/**
 * * Personal Project for 2D Motion Unit
 *  @author student
 * 	Purpose: This is the mother class. It simulates particle collisions inside a box. Gravity is on a dial: it can be cranked up or shut off completely depending on the presets. The box's dimensions and number of particles are also variable, and can be shifted around how the user pleases. 
 * 
 * 
 * A note!: if you do shift around the size of the display window for the box, also shift either the threshold distance for collisions or the particle radius size so that they match each other. The threshold distance should be the radius of the particle, so that when two particles approach each other a collision is sense when the are overlapping at all.
 * 
 * One Bug: the timing is completely off with real world. It is programmed to have a due step every .02 seconds, but they do not occur with nearly that frequency. What is displayed is a slow-motion simulation\
 * Another Bug (FIXED): some balls would escape box every once and a while, meaning that the collision was being unnecessarily disallowed. This was because the condition that there was just a collision with that wall was fulfilled, thereby negating any more collision proceedings. However, when the was still in the threshold range after an initial collision with a wall, occasionally would hit another particle which would send it back towards wall. With the collision proceedings still disabled, the ball would go straight through the wall this time
 * 	The Fix: instead of having the collision proceedings only occur always when the particle had not just collided with that same wall, I added another condition where the collision proceedings could occur: whether or not the particle was still in the threshold after the collision, IF the particle's velocity was directing it back towards the wall for some reason (presumably because it hit another particle while still in the threshold range of the wall in question)
 * Bug3: When a clash is detected, then the particle changes velocity and moves again, the same clash may still indeed be detected. This would occur because clashs are detected based on whether two entities are within a 'threshold' distance of each other (for which there is an double variable 'threshold'). Particles often stay within the threshold after the initial collision, and therefore have the clash function apply again, which causes errors with their velocities. In 
 * 
 * 
 * Bonus to do: threshold length should correspond with radius of particles. displayframe.getBounds() will help to solve this issue
 */

public class DiagonalClashTestDM extends AbstractSimulation{
	int molecules_num; 								//the number of molecules
	ParticleDM [] molecules; 						//the air molecules	
	
	PlotFrame disp = new PlotFrame("x", "y", "Diagonal Clash Test Display");	//display frame for balls in box
	
	//average distribution display frame
	DisplayFrame dist = new DisplayFrame("Layer Number (in order of increasing height)", "Average Number of Particles", "Particle Distribution Cumulative Average per Layer");
	//displays the average distribution of
	DisplayFrame instantaneous_dist = new DisplayFrame("Layer Number (in order of increasing height)", "Current Number of Particles", "Particle Distribution Instantaneous per Layer");
	
	//PlotFrame vt = new PlotFrame("time (seconds)", "velocity (m/s)","Velocity-Time Graph"); //velocity-time graph for single molecule
	double interval = .01;			//small time interval used to approximate the position of balls
	Circle yval = new Circle(); //mimicks y value of the first particle
	Circle xval = new Circle(); //mimicks x value of the single particle (sliding along x axis)
	double box_dims[] = new double [2]; //stores dimensions of box in which particle exist
	Random ran = new Random(); //random generator
	int ds_counter = 0; //do step counter. First do step is zero, each next one is labeled as the number one more than the previous
	double threshold; //the distance between point-particles for there to exist a collision 
	FluidDM air;
	
	CollisionDM romance;				//the collision object that this class cycles through, gives access
	ParticleDistributionDM expanse;		//the particle distribution
	
	boolean [] righthits; //parallel array for molecules, keeps track of do step number of all collisions with  right  hand wall (shouldn't have the same one again so soon). counts the right and left wall as same
	boolean [] lefthits; //parallel array for molecules, keeps track of do step number for all collisions with left wall (counts top or bottom as same)
	boolean [] tophits; //parallel array for molecules, keeps track of do step number for all collisions with top wall (counts top or bottom as same)
	boolean [] downhits; //parallel array for molecules, keeps track of do step number for all collisions with down wall (counts top or bottom as same)
	boolean [][]	par_hits;  //unfortunately, a little less than half of these storage spaces will not be necessary due to symmetry ([m][n] would store same data as [n][m])
	double[] density_distribution; //stores the distribution of particles by layer, where the number of balls in the lowest slice of the box is (i.e. from -20 to -19) stored in the first index
	/**
	 * Initializes a whole batch of balls, and each of their characteristics (e.g., position, velocities, directions,) randomly
	 * @param anglerange range of initial angles that each of these balls will possess
	 * @param anglemin minimum angle of this range
	 */
	public void initializeBalls (){
		//for each ball in array
		double angle;	//angle temporary storage value;
		boolean independence; //whether or not this ball's position is independent from others: if there is overlap, then this is false, and the position is calculated again
		//for every molecule
		for (int i = 0; i<molecules.length; i++){
			independence = false; //initialize independence as false, so program goes through the while loop below at least once
			angle = ran.nextDouble()*2*Math.PI - Math.PI;
			double x = 0; //initial x coordinate of the given particle
			double y = 0; //initial y coordinate of the given particle
			double speed = ran.nextDouble()*control.getDouble("Speed Max"); //generating random speed of this particle
			//while this ball's position is overlapping with any other previously established one
			while(independence == false){
				x = ran.nextDouble()*box_dims[0]-box_dims[0]/2; //randomly generates x position within the box's constraints
				y = ran.nextDouble()*box_dims[1] - box_dims[1]/2; //randomly generates y position within the box's constraints
				independence = true; //set default as true, so that it won't go through while loop again and regenerate coordinates, unless...
				//for each previously established particle
				outer: for (int j = 0; j<i; j++){
					//if the x coords are within threshold
					if ((Math.abs(molecules[j].getX() - x)<threshold)){
						//if y coords are also within threshold
						if ((Math.abs(molecules[j].getY()- y)<threshold)){
							independence = false; //this particle's positioning is not truly independent
							break outer; //break out of this loop
						}
					}
				}
			}
			//initialize each ball: double ixpos, double iypos, double yai, double xai, double angle_init, double speed_init, double radius
			//initializes molecule i
			molecules[i] = new ParticleDM(x,	y,	 control.getDouble("Universal Y Acceleration"),		0, 	angle, speed, control.getInt("Radius")); 	//creates a molecule;	
		}
	}
	/**
	 * function meant as supplement to certain bug correction in ParticleInteraction.clash method. In order to resolve bug 3 (as outlined at top of program), there are three arrays established to keep track of when (in which do step)  particles' collide with each other and the two walls. This function sets equal the returned arrays to the arrays in this class.
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
		romance.molecules = molecules; //sets the molecule arrays in each the collision class and this mother class equal, so info we are working with is consistent
		expanse.molecules = molecules;//sets the molecule arrays in each the particleDistribution class and this mother class equal, so info we are working with is consistent
		
		//for each particle
		for(int n = 0; n<molecules.length; n++)
		{
			molecules = romance.clashReal(threshold, n, box_dims); //run the collision function, which essentially tests to see if there is a collision, and if there is, recalculates all involved particles' velocities
			//If there was just a clash
			if((romance.thereisaclash == true)){
				//do not recalculate velocity and acceleration. Doing so will mess with the recalculations done above in the collision method
				molecules[n].components(air, interval, control.getDouble("Universal Y Acceleration"), false); //establishes x and y components for the particle
			}
			else{
				//otherwise, if no collision, recalculate the positions, velocities, and accelerations recursively as normal
				molecules[n].components(air, interval, control.getDouble("Universal Y Acceleration"), true); //establishes x and y components for the particle
			}


			this.initializeBugData(); //set all the 'bug-catching' data in collision class equal to info in DiagonalClashTestDM have here, so information is consistent
			molecules[n].setXY(molecules[n].getXpos(), molecules[n].getYpos()); //reset the x and y coords to new positions
			disp.addDrawable(molecules[n]); //add the molecule to the DisplayFrame

			//Initialize two balls to stay on the x and y axis respectively to shadow molecules[0]. This will be an easy way to track one particle in particular visually
			yval.setXY(-box_dims[0]/2, molecules[n].getY()); //sets position this "yball" to left edge of the box, and the height of molecules[0] at this time
			xval.setXY(molecules[n].getX(), -box_dims[1]/2);//sets position this "xball" to bottom edge of the box, and the x location of molecules[0] at this time
			//add both to the display frame
			disp.addDrawable(yval);
			disp.addDrawable(xval);
		}
		expanse.averageDensityDistribution(ds_counter); //recalculates average density distribution per layer of the box
		expanse.plotDistribution(expanse.density_distribution, ds_counter, dist, 50); //plots distribution of the average number of particles per each layer (averaging from records of each doStep
		expanse.plotDistribution(expanse.densityDistribution(), ds_counter, instantaneous_dist, 4); //plots distribution of current number of particles per each layer of box
		
		
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
		box_dims[0] = control.getDouble("Box Width"); //box's width
		box_dims[1] = control.getDouble("Box Height"); //box's height
		density_distribution = new double [control.getInt("Density Distribution Layers")]; //the array has a spot for each layer of the box. What is stored in the index for each layer is the average number of balls that have been in that layer.

		initializeBalls(); //initialize all aspects of the molecules: their speeds, positions, etc.

		disp.setVisible(true); //set display frame as visible
		disp.setDefaultCloseOperation(3); //weird thing you just gotta do
		//		

		disp.clearData(); //wipe all data
		disp.clearDrawables(); //wipe all drawables
		dist.clearDrawables();
		disp.setPreferredMinMax(-.75*box_dims[1], .75*box_dims[1], -.75*box_dims[1], .75*box_dims[1]);

		//Draw boxs for molecules to bounce around in
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
		disp.addDrawable(walls);//display the walls

		//initialize these bug preventative record keeping of whether a different types of collisions have happened in the immediate past of any molecule that the program creates
		righthits = new boolean [molecules.length]; //parallel array for molecules, keeps track of do step number of all collisions with left or right  hand wall (shouldn't have the same one again so soon). counts the right and left wall as same
		lefthits = new boolean [molecules.length]; //parallel array for molecules, keeps track of do step number for all collisions with bottom  or top wall (counts top or bottom as same)
		tophits = new boolean [molecules.length];
		downhits = new boolean [molecules.length];
		par_hits = new boolean [molecules.length][molecules.length]; //each molecule has a 'romance' with each other molecule when measuring their collisions with each other
		
		romance = new CollisionDM(molecules, par_hits, righthits, lefthits, tophits, downhits); //intializes this temp particle interaction
		expanse = new ParticleDistributionDM(molecules, control.getInt("Density Distribution Layers"), box_dims[1]/control.getInt("Density Distribution Layers"), density_distribution);
	
		//unimportant settings for the plot distribution display tables
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
		control.setValue("Radius", 2); //default radius size
		control.setValue("Air Resistance Constant", 0); //air resistance constant 0
		control.setValue("X Position", 0); // x position is zero
		control.setValue("Y Position", 0); //y position is zero
		control.setValue("Universal Y Acceleration", -10); //acceleration default is zero
		control.setValue("Threshold", .2); //threshold for collisions default is .5
		control.setValue("Number of Particles", 200); //2 particles
		control.setValue("Speed Max", 30); //maximum speed particles can have; default is 20m/s
		control.setValue("Box Width", 2);
		control.setValue("Box Height", 20);
		control.setValue("Density Distribution Layers", control.getInt("Box Height"));
	}
	/**
	 * launches the simulation
	 * @param args I don't know what this does
	 */
	public static void main(String[] args) {
		SimulationControl.createApp(new DiagonalClashTestDM());
	}

}
