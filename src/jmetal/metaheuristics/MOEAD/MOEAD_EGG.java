//  MOEAD.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.metaheuristics.MOEAD;

import jmetal.core.*;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.Ranking;
import jmetal.util.vector.TwoLevelWeightVectorGenerator;
import jmetal.util.vector.VectorGenerator;

import java.util.Vector;

import Jama.Matrix;
import jmetal.util.comparators.ObjectiveComparator;

public class MOEAD_EGG extends Algorithm {

  private int populationSize_;
  /**
   * Stores the population
   */
  private SolutionSet population_;
  /**
   * Z vector (ideal point)
   */
  double[] z_;
  
  double r_; 
  double t_;
  /**
   * Lambda vectors
   */
  //Vector<Vector<Double>> lambda_ ; 
  double[][] lambda_;
  /**
   * T: neighbour size
   */
  int T_;
  int H_ = 12;
  /**
   * Neighborhood
   */
//  int[][] neighborhood_;
  /**
   * delta: probability that parent solutions are selected from neighbourhood
   */
  double delta_;
  /**
   * nr: maximal number of solutions replaced by each child solution
   */
  int nr_;
  String functionType_;
  int evaluations_;
  
  /* if only one layer is adopted, div2_=0 */
	int div1_;  // divisions in the boundary layer
	int div2_;  // divisions in the inside layer
	
  /**
   * Operators
   */
  Operator crossover_;
  Operator mutation_;

  String dataDirectory_;

  /** 
   * Constructor
   * @param problem Problem to solve
   */
  public MOEAD_EGG(Problem problem) {
    super (problem) ;

    functionType_ = "_PBI";

  } // DMOEA

  public SolutionSet execute() throws JMException, ClassNotFoundException {
    int maxEvaluations;

    evaluations_ = 0;
    maxEvaluations = ((Integer) this.getInputParameter("maxEvaluations")).intValue();
    
    div1_ = ((Integer) this.getInputParameter("div1")).intValue();

	div2_ = ((Integer) this.getInputParameter("div2")).intValue();
	
	/* generate two-layer weight vectors */
	VectorGenerator vg = new TwoLevelWeightVectorGenerator(div1_, div2_, problem_.getNumberOfObjectives());
	lambda_ = vg.getVectors();
	
	/*the population size is the same with the number of weight vectors*/
	populationSize_ = vg.getVectors().length;
    
//    populationSize_ = ((Integer) this.getInputParameter("populationSize")).intValue();
//    dataDirectory_ = this.getInputParameter("dataDirectory").toString();
//    System.out.println("POPSIZE: "+ populationSize_) ;

    population_ = new SolutionSet(populationSize_);

    T_ = 20;
    delta_ = 0.9;
    nr_ = 2;
    
    r_ = 1.0;
    t_ = 0.0;
//    neighborhood_ = new int[populationSize_][T_];

    z_ = new double[problem_.getNumberOfObjectives()];
    //lambda_ = new Vector(problem_.getNumberOfObjectives()) ;
//    lambda_ = new double[populationSize_][problem_.getNumberOfObjectives()];

    crossover_ = operators_.get("crossover"); // default: DE crossover
    mutation_ = operators_.get("mutation");  // default: polynomial mutation

    // STEP 1. Initialization
    // STEP 1.1. Compute euclidean distances between weight vectors and find T
//    initUniformWeight();
    //for (int i = 0; i < 300; i++)
   // 	System.out.println(lambda_[i][0] + " " + lambda_[i][1]) ;
    
//    initNeighborhood();

    // STEP 1.2. Initialize population
    initPopulation();

    // STEP 1.3. Initialize z_
    initIdealPoint();
    // STEP 2. Update
    do {
      int[] permutation = new int[populationSize_];
      Utils.randomPermutation(permutation, populationSize_);
      Ranking ranking = new Ranking(population_);
      SolutionSet front = ranking.getSubfront(0);

      front.Suppress();
      SolutionSet KP = null;
      if(front.size()>problem_.getNumberOfObjectives())
    	  KP = findingKneePoint(front);
      if(KP==null){
    	  KP = population_;
      }

      for (int i = 0; i < populationSize_; i++) {
        int n = permutation[i];

        // STEP 2.1. Mating selection
		int r1,r2;
		r1 = PseudoRandom.randInt(0, KP.size() - 1);
		do {
			r2 = PseudoRandom.randInt(0, KP.size() - 1);
		} while (r2==r1);
        // STEP 2.2. Reproduction
        Solution child;
        Solution[] parents = new Solution[3];

        parents[1] = KP.get(r1);
        parents[2] = KP.get(r2);
        parents[0] = population_.get(n);
        child = (Solution) crossover_.execute(parents);
        // Apply mutation
        mutation_.execute(child);

        // Evaluation
        problem_.evaluate(child);
        
        evaluations_++;

        // STEP 2.3. Repair. Not necessary

        // STEP 2.4. Update z_
        updateReference(child);
        // STEP 2.5. Update of solutions
        updateProblem(child, n);
      } // for
    } while (evaluations_ < maxEvaluations);

    return population_;
  }

  /**
   * 
   */
  public void initPopulation() throws JMException, ClassNotFoundException {
    for (int i = 0; i < populationSize_; i++) {
      Solution newSolution = new Solution(problem_);

      problem_.evaluate(newSolution);
      evaluations_++;
      population_.add(newSolution) ;
    } // for
  } // initPopulation

  /**
	 * Initialize the ideal objective vector
	 * 
	 * @throws JMException
	 * @throws ClassNotFoundException
	 */
	void initIdealPoint() throws JMException, ClassNotFoundException {
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++)
			z_[i] = 1.0e+30;

		for (int i = 0; i < populationSize_; i++)
			updateReference(population_.get(i));
	} // initIdealPoint
	  
  /**
   * 
   * @param individual
   */
  void updateReference(Solution individual) {
    for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
      if (individual.getObjective(n) < z_[n]) {
        z_[n] = individual.getObjective(n);
      }
    }
  } // updateReference
  /**
   * @param individual
   * @param id
   * @param type
   */
  void updateProblem(Solution indiv, int id) {
	    // indiv: child solution
	    // id:   the id of current subproblem
	    // type: update solutions in - neighborhood (1) or whole population (otherwise)
	    int size=population_.size();
	    int time;

	    time = 0;
	    int i=0;

	    int[] perm = new int[size];

	    Utils.randomPermutation(perm, size);
	    for (i = 0; i < size; i++) {
	      int k = perm[i];      // calculate the values of objective function regarding the current subproblem
	      double f1, f2;

	      f1 = fitnessFunction(population_.get(k), lambda_[k]);
	      f2 = fitnessFunction(indiv, lambda_[k]);

	      if (f2 < f1) {
	    	population_.replace(k, new Solution(indiv));
	        //population[k].indiv = indiv;
	        time++;
	      }
	      // the maximal number of solutions updated is not allowed to exceed 'limit'
	      if (time >= nr_) {
	    	  return;
	      }
	    }
	  } // updateProblem

  /**
	 * Calculate the dot product of two vectors
	 * 
	 * @param vec1
	 * @param vec2
	 * @return
	 */
	public double innerproduct(double[] vec1, double[] vec2) {
		double sum = 0;

		for (int i = 0; i < vec1.length; i++)
			sum += vec1[i] * vec2[i];

		return sum;
	}

/**
	 * Calculate the norm of the vector
	 * 
	 * @param z
	 * @return
	 */
	public double norm_vector(double[] z) {
		double sum = 0;

		for (int i = 0; i < problem_.getNumberOfObjectives(); i++)
			sum += z[i] * z[i];

		return Math.sqrt(sum);
	}
	
	double fitnessFunction(Solution indiv, double[] lambda) {
		double fitness;
		fitness = 0.0;

		if (functionType_.equals("_TCHE1")) {
			double maxFun = -1.0e+30;

			for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
				double diff = Math.abs(indiv.getObjective(n) - z_[n]);

				double feval;
				if (lambda[n] == 0) {
					feval = 0.0001 * diff;
				} else {
					feval = diff * lambda[n];
				}
				if (feval > maxFun) {
					maxFun = feval;
				}
			} // for

			fitness = maxFun;
		} else if (functionType_.equals("_TCHE2")) {
			double maxFun = -1.0e+30;

			for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
				double diff = Math.abs(indiv.getObjective(i) - z_[i]);

				double feval;
				if (lambda[i] == 0) {
					feval = diff / 0.000001;
				} else {
					feval = diff / lambda[i];
				}
				if (feval > maxFun) {
					maxFun = feval;
				}
			} // for
			fitness = maxFun;
		} else if (functionType_.equals("_PBI")) {
			double theta; // penalty parameter
			theta = 5.0;

			// normalize the weight vector (line segment)
			double nd = norm_vector(lambda);
			for (int i = 0; i < problem_.getNumberOfObjectives(); i++)
				lambda[i] = lambda[i] / nd;

			double[] realA = new double[problem_.getNumberOfObjectives()];
			double[] realB = new double[problem_.getNumberOfObjectives()];

			// difference between current point and reference point
			for (int n = 0; n < problem_.getNumberOfObjectives(); n++)
				realA[n] = (indiv.getObjective(n) - z_[n]);

			// distance along the line segment
			double d1 = Math.abs(innerproduct(realA, lambda));

			// distance to the line segment
			for (int n = 0; n < problem_.getNumberOfObjectives(); n++)
				realB[n] = (indiv.getObjective(n) - (z_[n] + d1 * lambda[n]));
			double d2 = norm_vector(realB);

			fitness = d1 + theta * d2;
		} else {
			System.out.println("MOEAD.fitnessFunction: unknown type "
					+ functionType_);
			System.exit(-1);
		}
		return fitness;
	} // fitnessEvaluation
	/**
	 * It is note that the POP only the first front
	 * @param pop
	 * @param T
	 * @param r
	 * @param t
	 * @return
	 */
	SolutionSet findingKneePoint(SolutionSet pop) {
		Matrix Pop_M, ones_M, ExtremePoint_M, Hyperplane, Distance_M;
		SolutionSet KneePoint = new SolutionSet();
		SolutionSet ExtremePoint = new SolutionSet(problem_.getNumberOfObjectives());
		SolutionSet TempPop = new SolutionSet(pop.size());
		double [] f_max = new double [problem_.getNumberOfObjectives()];
		double [] f_min = new double [problem_.getNumberOfObjectives()];
		for (int i=0;i<problem_.getNumberOfObjectives();i++) {
			f_max[i]=pop.get(0).getObjective(i);f_min[i]=pop.get(0).getObjective(i);
		}
		for(int i=0;i<pop.size();i++) {
			TempPop.add(pop.get(i));
			for (int j=0;j<problem_.getNumberOfObjectives();j++) {
				if(pop.get(i).getObjective(j)>f_max[j]) f_max[j]=pop.get(i).getObjective(j);
				if(pop.get(i).getObjective(j)<f_min[j]) f_min[j]=pop.get(i).getObjective(j);
			}
		}
		ones_M = new Matrix(problem_.getNumberOfObjectives(),1,1.0);
		//1. find extrem solution
		for(int i=0;i<problem_.getNumberOfObjectives();i++){
			ObjectiveComparator oc = new ObjectiveComparator(i);
			int best = TempPop.indexBest(oc);
			ExtremePoint.add(TempPop.get(best));
			TempPop.remove(best);
		}
		for (int i=0;i<ExtremePoint.size();i++){
			KneePoint.add(ExtremePoint.get(i));
		}
		
		//2. calculate extreme hyperplane
		ExtremePoint_M = new Matrix(ExtremePoint.writeObjectivesToMatrix());
		if(!ExtremePoint_M.lu().isNonsingular()) return ExtremePoint;
		Hyperplane = ExtremePoint_M.inverse().times(ones_M);
		double normHyperplane = Hyperplane.norm2();
		//Calculate the distance between each solution in TempPop and L
		Pop_M = new Matrix(TempPop.writeObjectivesToMatrix());
		Distance_M = Pop_M.times(Hyperplane);
		
		//Calculate R point
		r_ = r_*Math.exp(-((1.0-t_/0.5)/(double)problem_.getNumberOfObjectives()));
		double[] R = new double[problem_.getNumberOfObjectives()];
		for(int i=0;i<problem_.getNumberOfObjectives();i++)
			R[i] = (f_max[i]-f_min[i])*r_;
		
		//Add KneePoint
		int remain = TempPop.size();
		boolean Choose[] = new boolean [Distance_M.getRowDimension()];
		for(int i=0;i<Distance_M.getRowDimension();i++) Choose[i]=false;
		while(remain>0){
			//find the largest distance in Distance_M that have contain in KneePoint
			double largest_dist = Double.NEGATIVE_INFINITY;
			int idx = -1;
			for(int i=0;i<Distance_M.getRowDimension();i++){
				if(!Choose[i]){// judge the ith individual selected/remove or not
					double dist = -(Distance_M.get(i, 0)-1.0)/normHyperplane;
					if(dist > largest_dist){
						largest_dist = dist;
						idx = i;
					}
				}
			}
			if(idx!=-1){
				KneePoint.add(TempPop.get(idx));
				Solution currentKneeP = TempPop.get(idx);
				//remove neighbors
				for(int i=0;i<TempPop.size();i++){
					if(!Choose[i]){
						Solution s = TempPop.get(i);
						boolean flag = true;
						for(int j=0;j<problem_.getNumberOfObjectives();j++){
							if(Math.abs(s.getObjective(j)-currentKneeP.getObjective(j))>R[j]) {flag = false;break;}
						}
						if(flag){
							Choose[i]=true;remain--;
						}
					}
				}
			}else{ // can not find a largest distance in the remain individual
				break;
			}
		}
		
		//update t_
		t_ = (double)KneePoint.size()/(double)pop.size();
		return KneePoint;
	}
} // MOEAD

