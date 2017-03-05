package particleRealDM;

import java.awt.Color;
import java.util.ArrayList;

import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display.DrawableShape;
import org.opensourcephysics.frames.PlotFrame;
/*
 * Green Monster Baseball Problem
 * @author Davey Morse
 * 
 * Provides both a physical and mathematical simulation to determine smallest velocity at which there is an initial angle for a ball to travel over the Green Monster
 * 
 * Define 'Batch' of Particles as a given set of ParticleDM's that has been initialized
 */
public class BaseballDM extends AbstractSimulation{
	int precision = 30; 													//number of balls launched for each given velocity
	ParticleDM [] b = new ParticleDM [precision]; 							//the particle that has the properties as ball on the screen
	PlotFrame homerun = new PlotFrame("x", "y", "Initial Problem");			//Baseball display frame
	double alpha = 0.02;													//air resistance constant - takes into account mass, fluid viscosity
	double interval = .1;													//small time interval used to approximate the position of balls
	double monsterlocation = 100;												//X coordinate of wall
	double monsterheight = 10;												//height of wall
	double initial_velocity = 0;											//initial velocity at which all balls in a given 'Batch' are launched
	double [] heights = new double [b.length]; 								//parallel array for particles: stores value of height of each when x = 100
	double [] initialangles = new double [b.length]; 						//parallel array for particles: stores value of initial angle at which each is launched
	double v_decrease = 3;
	
	//incrementing angle, not y vs. x components
	double anglemin = 0.01;													//minimum constraint of angle range for a given Batch of launched balls
	double anglerange = Math.PI/2-.01;										//size of range of angles for a given Batch
	double counter = 0;														//counter
	double balls_done = 0;													//number of balls that have landed so far in a given batch
	boolean [] ballsdone = new boolean [b.length];							//parallel array for particles that marks whether ball has landed or passed wall yet

	/**
	 * finds the index of the balls with smallest and largest angles to pass over the wall
	 *  based on the way that balls distribute; if angle is too high or low, balls will not reach x-position of wall. The purpose of this search is to find the range of angles that work for this velocity.
	 * @param height the height of the wall
	 * @param heights array of heights of particles at the wall (with negative values if they never reached the wall)
	 * @return int [] with index number of first and last particle
	 */
	public static int [] findFirstAndLast (double height, double [] heights){
		int first = 0; 	//index of first ball to make it over
		int last = 0; 	//index of last ball to make it over
		// for each ball, of least to greatest initial angle
		first: for (int i = 0; i< heights.length; i++){
			//if its height at the wall's location is greater than wall's height
			if (heights[i]> height){
				first = i; //record this ball's index
				break first; //stop the search for first ball's index
			}
		}
		//for each ball, of greatest to least initial angle
		last: for (int i = (heights.length) - 1; i>=0; i--){
			//if the ball has made it over
			if (heights[i]> height){
				last = i; //record this ball's index
				break last; //stop search for last ball
			}
		}
		return new int []{first, last}; //return [index of first ball], [index of last ball]
	}
	public double [] componentsFromAngleMag(double angle){
		double ycomp = Math.sqrt((Math.pow(initial_velocity, 2)*Math.pow((Math.tan(angle)), 2))/((Math.pow((Math.tan(angle)), 2)) + 1));
		return new double [] {ycomp, ycomp/Math.tan(angle)};
	}
	
	/**
	 * Initializes a whole batch of balls
	 * @param anglerange range of initial angles that each of these balls will possess
	 * @param anglemin minimum angle of this range
	 */
	public void initializeBalls (double anglerange, double anglemin){
		//for each particle
		for (int i = 0; i<b.length; i++){
			double angle = (((double) i)/((double)precision))*(anglerange) + anglemin;	//sets this ball's angle to a somewhere within range of angles
			initialangles [i] = angle;													//records this angle
			//initialize each ball: double ixpos, double iypos, double yai, double xai, double angle_init, double speed_init, int radius
			b[i] = new ParticleDM(control.getDouble("X Pos"),	control.getDouble("Y Pos"), 	control.getDouble("Initial Vertical Acceleration"),		control.getDouble("Initial Horizontal Acceleration"), 	angle, initial_velocity, 4); 								//creates a particle;												//initializes characteristics of the ball
			ballsdone[i] = false;	//reset this record of whether or not balls have finished their paths; these fresh babies have not
		}
	}
	/**
	 * Determines if a ball has hit the ground
	 * @param yposition y position of a given ball
	 * @return boolean, indicating whether ball is at ground level
	 */
	public boolean hitGround (double yposition){
		boolean hit_ground = false;
		if(yposition<=0){
			hit_ground = true;
		}
		return hit_ground;
	}
	/**
	 * Determines if a ball has arrived at x coordinate of the wall
	 * @param xposition the given ball's x position
	 * @return boolean, indicating whether ball's x coord is at the wall
	 */
	public boolean atWall (double xposition){
		boolean atWall = false;
		if (xposition>=100){
			atWall = true;
		}
		return atWall;
	}

	protected void doStep() {
		//narrowing pool of particles down to those that got over the wall
		//the point here is to find if 1) particle never reaches wall or 2) height it reaches at x value of wall

		//for each particle
		homerun.clearData();
		outer: for(int n = 0; n<b.length; n++)
		{
			if(ballsdone[n] ==true){
				continue outer;
			}
			//y component
			b[n].setYa(b[n].getYai() - alpha*b[n].getYv());	 //changes acceleration by -alpha*p.getYv(), which represents air resistance as a function of velocity
			b[n].setYv(b[n].getYv() + b[n].getYa()*interval); //changes velocity due to acceleration in given time interval
			b[n].setYpos(b[n].getYpos() + b[n].getYv()*interval); //changes y position due to velocity in given interval

			//x component
			b[n].setXa(b[n].getXai() - alpha*b[n].getXv());	 //changes acceleration by -alpha*p.getXv(), which represents air resistance as a function of velocity
			b[n].setXv(b[n].getXv() + b[n].getXa()*interval); // changes velocity due to acceleration in given time interval
			b[n].setXpos(b[n].getXv()*interval + b[n].getXpos()); //changes position due to velocity over given time interval


			//checks to see
			//if ball has hit ground
			if (hitGround(b[n].getYpos()) == true){
				heights[n] = -1; //disqualifies this particle and it's specs, for purpose of the calculation
				ballsdone[n] = true;
				balls_done++;
				homerun.setMarkerColor(0, Color.BLACK);
				continue outer;
			}

			//if it hasn't yet hit the ground, if x value is equal to monster location
			if (atWall(b[n].getXpos()) == true){
				heights[n] = b[n].getYpos();
				ballsdone[n] = true;
				balls_done++;
				homerun.setMarkerColor(0, Color.BLACK);
				continue outer;
			}
			
			//set marker color based on angle
			//add point on graph
			homerun.append(0, b[n].getXpos(), b[n].getYpos());
		}
	if (balls_done== b.length)
	{
		//calculate angle of first and last balls to make it over
		//findFirstAndLast(monsterheight, heights)[0]; index of first ball to make it over
		//System.out.println("findFirstAndLast(monsterheight, heights)[0]: " + findFirstAndLast(monsterheight, heights)[0]);

		anglemin = initialangles[findFirstAndLast(monsterheight, heights)[0]]; //sets minimum angle of new range to that of first ball to pass over wall (because they were launched in order of decreasing angle, first will have smallest angle)
		anglerange = initialangles[findFirstAndLast(monsterheight, heights)[1]] - anglemin; //sets angle range to last ball's angle minus anglemin
			if (anglerange == 0) //if the last and first ball are one and the same
			{ 
				System.out.println("BREAK");
				initial_velocity+= v_decrease;
				v_decrease/=10;
				anglerange = v_decrease;
				
			}
		initial_velocity -= v_decrease;

		counter++; //counter for each reset

		System.out.println(counter + ")  Minimum Angle (radians): " + anglemin);
		System.out.println("= "+ (anglemin/(2*Math.PI))*360 + " degrees");
		System.out.println("Initial Velocity: " + initial_velocity);

		//reinitializing balls
		initializeBalls(anglerange, anglemin);
		balls_done = 0;
		//clearing all drawables from display frame (the past paths)
		homerun.clearData();
		homerun.setMarkerColor(0,Color.RED);
		}
	else{
	}
	}
	/**
	 * initializes necessary variables for simulation and formalities for the plotframe
	 */
	public void initialize() {
		double width = 4;
		DrawableShape green_monster = DrawableShape.createRectangle(monsterlocation+ width/2, monsterheight/2, width, monsterheight);
		green_monster.setMarkerColor(Color.GREEN, Color.BLACK);
		homerun.addDrawable(green_monster);
		interval = control.getDouble("Interval");//setting small time interval approximation
		alpha = control.getDouble("Air Resistance Constant");
		initial_velocity = control.getDouble("Initial Velocity");
		//incrementing angle, not y vs. x components
		initializeBalls(anglerange, anglemin);
		
		//formalities for frames
		homerun.setVisible(true);
		homerun.setPreferredMinMax(0, 150, 0, 100);
		//for VT graph
		//				yvt.setDefaultCloseOperation(3);
		//				yvt.setVisible(true);
		//				//vt.setPreferredMinMax(xmin, xmax, ymin, ymax);
	}
	/**
	 * resets all control panel values
	 */
	public void reset() {
		control.setValue("X Pos", 0); //default x position of particle
		control.setValue("Y Pos", 0); //default y position of particle
		control.setValue("Interval",interval); //default interval width for particles' trajectories
		control.setValue("Initial Velocity", 60); //default velocity
		control.setValue("Initial Vertical Acceleration", -9.8); //default vertical acceleration (due to gravity)
		control.setValue("Initial Horizontal Acceleration", 0);
		control.setValue("Radius", .5); //default radius size
		control.setValue("Color", Color.BLUE);
		control.setValue("Air Resistance Constant", 0);

	}
	/**
	 * launches the simulation
	 * @param args I don't know what this does
	 */
	public static void main(String[] args) {
		SimulationControl.createApp(new BaseballDM());
	}
}
