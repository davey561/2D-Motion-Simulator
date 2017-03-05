package particleRealDM;

import org.opensourcephysics.display.Trail;
import org.opensourcephysics.frames.DisplayFrame;
import org.opensourcephysics.frames.PlotFrame;
/**
 * This class deals with the analysis of particles' height distribution. In particle, the method of analysis generally goes as follows:
 * 	1. Box is split up into as many horizontal 'layers' as the user inputs. A layer is a horizontal slice of the box with a minimum and maximum y bound. Ex) layer 1 could be all the molecules currently between y = -20 and y = -19
 * 	2. after particles are set loose, for each layers in each do step, the number of particles in that layer are counted up.
 * 
 * This class answers the core question that was leading me with this assignment: how do particles vertically distribute when bouncing around against each other, when gravity is acting on them? The Display Frames that show up alongside the visual representation of box with particles inside, those graphs provide the answer to this question: that height seems to have a inversely linear correlation with the number particles in the layer
 * @author student
 *
 */
public class ParticleDistributionDM {
	ParticleDM [] molecules;
	int layers;
	double layerheight;
	double [] density_distribution;
	int modnumber = 10;
	
	//METHODS
	//Constructors
	public ParticleDistributionDM(ParticleDM[] m, int l, double lh, double [] dd){
		molecules = m;
		layers = l;
		layerheight = lh;
		density_distribution = dd;
	}
	
	/**
	 * Finds number of balls currently in each layer, at this point in time
	 * @return
	 */
		public double [] densityDistribution(){
			double [] densities = new double [layers]; //this is the array that stores all the necessary values
			//for each molecule
			for (int i = 0; i<molecules.length; i++){
				//for each layer
				layerloop: for (int l = 0; l<layers; l++){
					//if the given molecule is in that layer
					if ((molecules[i].getY()>= (l-layers/2)*layerheight)&&(molecules[i].getY()<(l-layers/2+1)*layerheight)){
						densities[l]++; //add one to the variable storing the number of balls per layer in the densities array
						break layerloop; //you won't find the same molecule in a different layer at this time, so don't continue to loop through the layers; move to the next molecule
					}
				}
			}
			return densities;
		}
		/**
		 * //a recursive function, to average, for each layer, the number of balls for all different do step times. For ex. ds1 --> 5 balls in layer x, ds2 --> 9 balls in layer x. This function should find that Layer x average is 7 balls 
			//density_distribution array will be changed. This array will always be set to the average number of balls per layer. The way the average is calculated is recursive: if the average was x by last dostep for layer l, and in this do step layer l has y balls: new average = (x*[past dostep number] + y)curent dostep number
		 * @param ds_number
		 */
		public void averageDensityDistribution(int ds_number){
			//for each layer
			for(int layer = 0; layer<layers; layer++){
				//do this recursive averaging step to calculate and record new average number of balls
				density_distribution[layer] = ((ds_number-1)*density_distribution[layer] + densityDistribution()[layer])/ds_number;
			}
			
		}
		/**
		 * Prints the density distribution array in an understandable way
		 */
		public void printDensDistribution(){
			for (int i = 0; i< density_distribution.length; i++){
				System.out.println("In layer " + i + ", " + (i-(layers/2)) + " to " + (i- layers/2 + 1) + "the average number of particles has been: " + density_distribution[i]);
			}
		}
		/**
		 * Plot the average distribution of balls for each layer. The x axis is the layer number, and the y axis is the number of balls
		 * Such a plot makes viewable the trends from layer to layer of the number of balls that are generally present in each. The relation between layer number and number of particles looks generally linear, although it may be asymptotic
		 * @param density_distribution the array that stores all the data about the number of balls
		 * @param ds_number do step number
		 * @param dist the Display frame to plot the density distribution
		 * @param modnumber //number of do steps before the display frame clears and resets before more data is added
		 */
		public void plotDistribution(double [] density_distribution, int ds_number, DisplayFrame dist, int modnumber){
			Trail temp = new Trail(); // a horizontal line for each layer, the height at each layer number (listed the x axis) marks the average number of balls in that layer
			//for each layer
			for (int layer = 0; layer<density_distribution.length; layer++){
				//plot this line at the height of the average number of balls there
				temp.addPoint(layer, density_distribution[layer]);
				temp.addPoint(layer+1, density_distribution[layer]);
			}
			dist.addDrawable(temp);
			
			// every modnumber dosteps, it clears the frame
			if(ds_number %modnumber == 0){
				dist.clearDrawables();
			}
		}
}
