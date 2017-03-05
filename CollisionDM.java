package particleRealDM;
/**
 * Very complicated (BONUS BONUS BONUS): CollsionDM simulates the dynamics of a elastic particle collisions- recalculates each particles' new speeds and direction after each interactions with a wall or other particle
 * @author Davey
 * **elastic means that momentum is conserved <--> energy in the system is also always conserved. Because all the particles in this system have the same mass, the total speeds of each of the particles when there is no acceleration is also conserved, although this is not true when acceleration is present.
 * 
 * Two Versions of the main function: clashSwap and clashReal. Both calculate the results of a collision of two objects. The place they differ is in their calculation of particle-particle collisions.
 * 	collisionSwap assumed that all collisions were head-on from the perspective of both particles. That is, it presumed that, in  a collision, if you create a reference frame either particle as stationary, the other particle will transfer fully it's velocity (according the aforementioned reference frame) to the other particle. 
 *			While this presumption may sound complex, what results from it is that in each collision, both particles SWAP velocities, hence the name.
 *			However, this system is too simplistic, as it does not take into account collisions that are not "head-on," i.e. collisions where particles do not completely transfer their velocities to each other. Because partial transferrences are almost infinitely more likely to occur than total swappings of velocity, I wrote the next version of the function, explained below.
 *	collisionReal includes all levels of transferences
 * 
 * Note: in comments, particle and molecule are used interchangeably to describes ParticleDM's. This hold true EXCEPT in comments in the collision function below, for which particle and molecule have different meanings
 * 
 * 
 * Bug: with a lot of balls, every once and a while a single one will escape. FIXED No longer! Introduced a new condition for collisions, will explain more below
 */		
public class CollisionDM{
	ParticleDM [] molecules;//molecules;
	FluidDM [] environ;

	//Used for both main functions clashSwap and clashReal
	//each of these arrays indicate that the given particle is still within threshold range of the left, right, top, bottom walls or another particle, respectively
	//in each array, there is an item for each separate particle
	boolean [] lefthits;
	boolean [] righthits;
	boolean [] tophits;
	boolean [] downhits;
	boolean [][] par_hits;

	boolean thereisaclash; //keeps track of whether there is a clash during this do step. If there was, the main class recalculates the positions of the involved particles slightly differently (although most of the calculations are taken care of below)

	//Constructor
	public CollisionDM(ParticleDM [] molecules, boolean [][] par_hits, boolean [] righthits, boolean [] lefthits,  boolean [] tophits,  boolean [] downhits){
		this.molecules = molecules;
		this.par_hits = par_hits;
		this.righthits = righthits;
		this.lefthits = lefthits;
		this.tophits = tophits;
		this.downhits = downhits;
	}
	//bug right now is that particle often stays within threshold of wall/particle it just hit, meaning the code is then applied again to that very particle even after it's left
	//current solution: making sure a couple do steps go by before its possible for those two entities to collide again
	//new and better solution (currently implemented instead): ensures that it exits clash range before the clash of same type (with vertical wall, horizontal wall, or other particle (three types)) happens again
	/**
	 * Produces information about the clash between a given particle and either a horizontal wall, a vertical wall, or another particle. 
	 * When particles hit, it swaps the velocities (direction and magnitudes. One particle becomes the other) (CURRENTLY IMPLEMENTED) This is too simplistic, as even though all momentum is conserved, different fractions of the x and y velocites are transferred depending on the angle of collision
	 * However, the actuality of elastic collisions more complicated than this (simulated in clashReal below the clashSwap method)
	 * Very Complicated, has taken a lot of refinement.
	 * @param thresh the max distance between two objects for there to be declared a collision
	 * @param particlenumber the number particle we are dealing with, as an index in the molecules array (the mother data set)
	 * @return molecules, the array of all molecules, parallel to the same array 'molecules' in the Simulation Class, but refined with particles' changes in velocity (direction and/or speed) when there is indeed a collision
	 */
	public ParticleDM[] clashSwap (double thresh, int particlenumber, double [] box_dims){
		//for each particle that has been initialized
		double [] tempvel = new double [2]; //storage value
		thereisaclash = false; //sets there is a clash to false, will change if there is one detected in proceeding code

		//for each molecule
		for(int n = 0; n<molecules.length; n++){
			//if we are comparing the given particle to a DIFFERENT one (because we don't want to be comparing a particle to itself)
			if (particlenumber != n){
				//if the distance between this two particle's is less than the threshold amount
				if (Math.sqrt(Math.pow(molecules[particlenumber].getY()-molecules[n].getY(), 2) + Math.pow(molecules[particlenumber].getX()-molecules[n].getX(), 2))<2*thresh){
					//if the last crash between these two particles happened at least 3 do steps ago
					if (par_hits[particlenumber][n] == false){
						//swap angles and velocities of the two particles that crash
						tempvel[0] = molecules[particlenumber].getV()[0]; //store particle velocity values: angle
						tempvel[1] = molecules[particlenumber].getV()[1]; //magnitude
						molecules[particlenumber].setV(new double [] {molecules[n].getV()[0], molecules[n].getV()[1]}); //set velocity of one particle to that of the other
						molecules[n].setV(new double [] {tempvel[0], tempvel[1]}); //set velocity of this molecule to other's past velocity
						par_hits[particlenumber][n] = true; //storing do step number of this particle collision
						par_hits[n][particlenumber] = true; //according to symmetry of the par_hits array, every relationship is repeated twice, [a][b] and [b][a]
						thereisaclash = true; //there was just a particle-particle clash
					}
				}
				else{
					par_hits[particlenumber][n] = false; //record if it is not within range
					par_hits[n][particlenumber] = false; //record if it is not within range (times two for the redudancy stated above)
				}
			}

		}
		//top wall
		//if particle is within threshold distance from top wall
		if((box_dims[1]/2 - molecules[particlenumber].getY()<thresh)){
			//if hasn't just previously hit the top wall or it is moving back towards the wall (presumably because another particle has hit it while still within threshold)
			if ((tophits[particlenumber] == false)|| (molecules[particlenumber].getYv()>0)){
				molecules[particlenumber].setAngle(molecules[particlenumber].getV()[0]*(-1)); //reflects direction of velocity about the wall
				tophits[particlenumber] = true; //recording that now particle has just hit the top wall
				thereisaclash = true; //there has just been a clash: true
			}
		}
		else{
			tophits[particlenumber] = false; //if not within range of the wall, record such
		}

		//bottom wall
		if((box_dims[1]/2 + molecules[particlenumber].getY()<thresh)){
			//if hasn't just hit the wall, or it is moving towards the wall (EITHER CONDITION)
			if ((downhits[particlenumber] == false) || (molecules[particlenumber].getYv()<0)){
				molecules[particlenumber].setAngle(molecules[particlenumber].getV()[0]*(-1)); //reflects the angle of the particle's velocity about the wall
				downhits[particlenumber] = true; //records that now particle has just hit the bottom wall
				thereisaclash = true; //there has just ben a clash: true
			}
		}
		//in every other case where it is not in threshold range
		else{
			downhits[particlenumber] = false; //record that it has not just hit wall
		}

		//right wall
		//if the particle is within the threshold distance of the right wall
		if((box_dims[0]/2-molecules[particlenumber].getX()<thresh)){
			//if the particle wasn't previously within range of the right wall, or it is moving back towards it
			if((righthits[particlenumber]==false)|| (molecules[particlenumber].getXv()>0)){
				molecules[particlenumber].setAngle(Math.PI - molecules[particlenumber].getV()[0]); //reflect the x component of the velocity back away from the wall
				righthits[particlenumber] = true; //stores do step number of this collision
				thereisaclash = true; //there has just been a clash: true
			}

		}
		//otherwise
		else{
			righthits[particlenumber] = false; //the particle is not within threshold range, record as so
		}

		//left wall
		//if within threshold range of left wall:
		if((box_dims[0]/2 + molecules[particlenumber].getX()<thresh)){
			//if it hasn't just hit or is moving back towards the wall:
			if((lefthits[particlenumber]==false)|| (molecules[particlenumber].getXv()<0)){
				molecules[particlenumber].setAngle(Math.PI - molecules[particlenumber].getV()[0]); //reflect the velocity (x component) of the particle back towards the right
				lefthits[particlenumber] = true; //stores do step number of this collision
				thereisaclash = true; //there has just been a clash: true
			}
		}
		//if not within threshold range of left wall
		else{
			lefthits[particlenumber] = false; //then record as such
		}

		return molecules; //return's the modified molecule class according to any collision that may have occurred
	}
	/**
	 * The more accurate particle clash simulating class: this function is now used in DiagonalClashTestDM, not the ClashSwap function
	 * This function also calculates elastic collisions, but does so by not just swapping the velocities, but taking into account the angle of collision and changing the respective affected and unaffected components of each particle's velocity accordingly
	 * The affected and unaffectd velocity components are perpendicular to each other
	 * @param thresh the threshold distance, within which a collision is said to have occurred
	 * @param particlenumber the number particle that we are dealing with (this is the index number of that particle within the molecules array)
	 * @param box_dims //the dimensions of the box within which these particles are interacting
	 * @return molecules the array of ParticleDM's that this function alters
	 */
	public ParticleDM[] clashReal (double thresh, int particlenumber, double [] box_dims){
		//'particle' is used as label for molecules[particlenumber]
		//'molecule' is used as shorthand for the molecule[n], which is rotated through in the for loop below
		//for each particle that has been initialized
		double [] tempvel = new double [2]; //storage value- unimportant
		thereisaclash = false; //sets there is a clash to false, will change if there is one detected in proceeding code
		double centerangle = 0; //initializes the central angle of the collision. This is the angle formed by the line r = theta and the ray with endpoint as 'particle's" center and in direction of the molecule's center

		//these velocity components are  for calculations, they represent the affected and unaffected vector components of the velocities of the particle and the molecule respectively. a1 u1 are for the particle, a2 u2 are for the molecule
		double [] a1 = new double [2]; //the component of the particle's velocity that is affected in collision
		double [] u1 = new double [2];// the component of the particle's velocity that goes unchanged by collision (perpendicular to the angle of crash)
		double [] a2 = new double [2]; //the component of the particle's velocity that is affected in collision
		double [] u2 = new double [2];//the component of the molecule's velocity that goes unaffected in collision

		double delta; //angle formed between initial velocity vector of the particle and the 'affected' velocity component
		double delta2; //angle formed between initial velocity vector of the molecule and the 'affected' velocity component


		//now, beginning the collision DETECTIoN PROCESS
		//for each molecule
		for(int n = 0; n<molecules.length; n++){
			//As a shorthand, in comments below I call the particle whose number was given a parameter for this method the Particle, I call the other molecule n the Molecule
			//if we are comparing the given particle to a DIFFERENT one
			if (particlenumber != n){
				//if the distance between this two particle's is less than the threshold amount (should eventually make threshold the sum of the particle's radii)
				if (Math.sqrt(Math.pow(molecules[particlenumber].getY()-molecules[n].getY(), 2) + Math.pow(molecules[particlenumber].getX()-molecules[n].getX(), 2))<2*thresh){
					if (par_hits[particlenumber][n] == false){ //if they haven't just hit
						//print that theres a crash, listing the number of the given particle, then the particle with which it is crashing
						System.out.println("CRASH: between Molecule #" + particlenumber + "and Molecule # " + n);

						//1) find angle of line between both particles' centers
						centerangle = Math.atan2((molecules[n].getY() - molecules[particlenumber].getY()), (molecules[n].getX()-molecules[particlenumber].getX()));
						//center angle is currently oriented in direction that the particle made impact when heading into the collision, it will be the direction of the force exerted onto the other particle

						//2) find the affected (A) and unaffected (U) components of each particle's velocity heading into the collision
						//particle's components

						//All the calculations here were done out on paper, and are easier to understand visually

						a1[0] = centerangle; //the affected component of particle's velocity will have the same direction as the center angle
						delta = findSmallestAngDiff(molecules[particlenumber].v[0], a1[0]); //the angles in the right triangle (made of velocity vector (hypotenuse), and two vector components a and u (legs)) formed by the affected vector and the total velocity vector
						//set it equal to the minimum of the differences between the particle velocity's angle and the centerangle OR the opposite of the particle's velocity angle and the center angle
						a1[1] = molecules[particlenumber].v[1]*Math.cos(delta); //the affected component's magnitude
						u1[0] = findOtherLegAngle(molecules[particlenumber].v[0], a1[0]); //uses a function to determine the angle of the other leg (i.e. the unaffected component). Is a fairly complex process
						u1[1] = molecules[particlenumber].v[1]*Math.sin(delta); //the unaffected component's magnitude

						a2[0] = findOppAngle(centerangle); //the affected component of the molecule's velocity will have the opposite direction of the center angle which points from particle to the molecule (rather than from molecule to the particle)
						delta2 =findSmallestAngDiff(molecules[n].v[0], a2[0]); //the angle formed by the molecule's initial velocity component  (v[0]) and its affected component a2
						a2[1] = molecules[n].v[1]*Math.cos(delta2); //the magnitude of the affected component
						u2[0] = findOtherLegAngle(molecules[n].v[0], a2[0]); //uses a method below to find the perpendicular unaffected component angle (  more complex than it might seem)
						u2[1] = molecules[n].v[1]*Math.sin(delta2); //the magnitude of the unaffected velocity component

						//swap affected components
						tempvel = a1; //stores the particle's affected component
						a1 = a2; //changes the particle's a to the molecule's a
						a2 = tempvel; // sets the molecule's the particle's stored a from line 203

						//recalculate new velocities
						molecules[particlenumber].setV(molecules[particlenumber].findPolarComponents(a1, u1));
						molecules[n].setV(molecules[n].findPolarComponents(a2, u2));

						par_hits[particlenumber][n] = true; //storing do step number of this particle collision
						par_hits[n][particlenumber] = true; //according to symmetry of area, every relationship is repeated twice
						thereisaclash = true; //there was just a particle-particle clash
					}
				}
				else{
					par_hits[particlenumber][n] = false; //if not within range, record such
					par_hits[n][particlenumber] = false; //redundancy resolved
				}
			}

		}
		//top wall
		//if particle is within threshold distance from top wall
		if((box_dims[1]/2 - molecules[particlenumber].getY()<thresh)){
			//if particle hasn't just previously hit the top wall or if it is moving towards the wall
			if ((tophits[particlenumber] == false)|| (molecules[particlenumber].getYv()>0)){
				molecules[particlenumber].setAngle(molecules[particlenumber].getV()[0]*(-1)); //reflects direction of velocity about the wall
				tophits[particlenumber] = true; //recording that now particle has just hit top wall
				thereisaclash = true; //this was indeed a clash
			}
		}
		else{
			tophits[particlenumber] = false; //records that particle was not within rnag eof the top wall
		}

		//bottom wall
		if((box_dims[1]/2 + molecules[particlenumber].getY()<thresh)){
			//if hasn't just hit
			if ((downhits[particlenumber] == false) || (molecules[particlenumber].getYv()<0)){
				molecules[particlenumber].setAngle(molecules[particlenumber].getV()[0]*(-1)); //reflects y velocity component, keeps x component unchanged
				downhits[particlenumber] = true; //records that now particle has just hit a horizontal wall
				thereisaclash = true; //this is indeed a clash
			}
		}
		//in every other case where it is not in threshold range
		else{
			downhits[particlenumber] = false; //record that it has not just hit wall
		}

		//right wall
		if((box_dims[0]/2-molecules[particlenumber].getX()<thresh)){
			if((righthits[particlenumber]==false)|| (molecules[particlenumber].getXv()>0)){
				molecules[particlenumber].setAngle(Math.PI - molecules[particlenumber].getV()[0]);
				righthits[particlenumber] = true; //stores do step number of this collision
				thereisaclash = true; //this is indeed a clash
			}

		}
		else{
			righthits[particlenumber] = false; //record that there was not a collision
		}

		//left wall
		//if within threshold range of left wall:
		if((box_dims[0]/2 + molecules[particlenumber].getX()<thresh)){
			//if it hasn't just hit or is moving in same direction:
			if((lefthits[particlenumber]==false)|| (molecules[particlenumber].getXv()<0)){
				molecules[particlenumber].setAngle(Math.PI - molecules[particlenumber].getV()[0]); //reflect the velocity about the wall
				lefthits[particlenumber] = true; //stores do step number of this collision
				thereisaclash = true; //this was indeed a clash, record as such
			}
		}
		//if not within threshold range of left wall
		else{
			lefthits[particlenumber] = false; //record that it is not within threshold range
		}

		return molecules;
	}

	public boolean isThereisaclash() {
		return thereisaclash;
	}
	public void setThereisaclash(boolean thereisaclash) {
		this.thereisaclash = thereisaclash;
	}
	public ParticleDM[] getMolecules() {
		return molecules;
	}

	public void setMolecules(ParticleDM[] molecules) {
		this.molecules = molecules;
	}

	public FluidDM[] getEnviron() {
		return environ;
	}

	public void setEnviron(FluidDM[] environ) {
		this.environ = environ;
	}

	public boolean[][] getPar_hits() {
		return par_hits;
	}

	public void setPar_hits(boolean[][] par_hits) {
		this.par_hits = par_hits;
	}
	public boolean[] getLefthits() {
		return lefthits;
	}
	public void setLefthits(boolean[] lefthits) {
		this.lefthits = lefthits;
	}
	public boolean[] getRighthits() {
		return righthits;
	}
	public void setRighthits(boolean[] righthits) {
		this.righthits = righthits;
	}
	public boolean[] getTophits() {
		return tophits;
	}
	public void setTophits(boolean[] tophits) {
		this.tophits = tophits;
	}
	public boolean[] getDownhits() {
		return downhits;
	}
	public void setDownhits(boolean[] downhits) {
		this.downhits = downhits;
	}

	/**
	 * Given an angle, finds the opposite angle within range (-pi, pi]
	 * @param ang
	 * @return
	 */
	public double findOppAngle(double ang){
		double opp;
		//if given angle is greater than zero
		if(ang>0){
			opp = ang-Math.PI; //then opposite angle is pi less
		}
		//else if angle given is less than zero
		else if(ang<0){
			opp = ang + Math.PI; //then opposite angle is pi more
		}
		//else if angle is 0
		else{
			opp = Math.PI;
		}
		return opp;
	}
	/**
	 * Finds the smallest angle difference between two lines (the smallest angle is acute or right always)
	 * @param ang1
	 * @param ang2
	 * @return
	 */
	public double findSmallestAngDiff(double ang1, double ang2)	{
		return Math.min(Math.min(Math.min(Math.abs(ang1-ang2), Math.abs(findOppAngle(ang1)-ang2)), Math.abs(findOppAngle(ang2) - ang1)), Math.abs(findOppAngle(ang1)- findOppAngle(ang2)));
	}
	/**
	 * Finds the angle of the other leg of a right triangle
	 * @param r a double [2] with magnitude and angle of hypotenuse
	 * @param a a double[2] with magnitude and angle of the given leg
	 * @return
	 */
	public double findOtherLegAngle(double r, double a){
		double u = 0;
		double pi = Math.PI;

		//all SPECIAL CASES
		//if in quad 2 and if a's also in quad 2
		if ((inRange(pi/2, pi, r))&&(inRange(pi/2, pi, a))){

			if(a>r){
				u = a-pi/2;
			}
			else if (a<r){
				u = a-3*pi/2;
			}
		}
		//else if r's in 2 and a's in three
		else if ((inRange(pi/2, pi, r))&&(inRange(-pi, -pi/2, a))){
			//a must be less than R
			u = a+3*pi/2;
		}
		//if r's in quadrant three and a's in 2
		else if((inRange(-pi, -pi/2, r))&& (inRange(pi/2, pi, a))){
			//if a's in quad 2
			//a must be greater than R
			u = a -3*pi/2;
		}
		//if both r and a are both in quad three
		else if ((inRange(-pi, -pi/2, r))&& (inRange(-pi, -pi/2, a))){
			if(a>r){
				u = a+3*pi/2;
			}
			else if (a<r){
				u = a + pi/2;
			}
		}
		//if r is pi and a is in 2
		else if((Math.abs(r) == pi)&&(inRange(pi/2, pi, a))){
			u = a-3*pi/2;
		}
		//else if r is pi and a is in 3
		else if((Math.abs(r) == pi)&&(inRange(-pi, -pi/2, a))){
			u = a+3*pi/2;
		}
		//for all other cases:
		else{
			if(a<r){
				u = a+pi/2;
			}
			else if(a>r){
				u = a-pi/2;
			}
		}
		return u;
	}
	/**
	 * Checks whether a given object is in the range given by the max and min
	 * @param min lower bound of the range
	 * @param max upper bound of the range
	 * @param q the double that this method is testing
	 * @return true or false whether or not q is in the range
	 */
	public boolean inRange (double min, double max, double q){
		if ((q>min) && (q<max)){
			return true;
		}
		else{
			return false;
		}
	}
}
