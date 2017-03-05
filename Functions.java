package particleRealDM;

import ParticleDM.ParticleDM;

public class Functions {
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
			//initialize each ball
			b[i] = new ParticleDM(control.getDouble("X Pos"),	control.getDouble("Y Pos"),	 componentsFromAngleMag(angle)[0], 	control.getDouble("Initial Vertical Acceleration"),		control.getDouble("Initial Horizontal Acceleration"), 	angle); 								//creates a particle;												//initializes characteristics of the ball
			ballsdone[i] = false;	//reset this record of whether or not balls have finished their paths; these fresh babies have not
		}
	}
}
