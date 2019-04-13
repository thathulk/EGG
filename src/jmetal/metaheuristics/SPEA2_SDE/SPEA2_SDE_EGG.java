//  SPEA2.java
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

package jmetal.metaheuristics.SPEA2_SDE;

import Jama.Matrix;
import jmetal.core.*;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.Spea2_SDEFitness;
import jmetal.util.comparators.ObjectiveComparator;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;

/** 
 * This class representing the SPEA2 algorithm
 */
public class SPEA2_SDE_EGG extends Algorithm{
          
  /**
   * Defines the number of tournaments for creating the mating pool
   */
  public static final int TOURNAMENTS_ROUNDS = 1;
  double r_;
  double t_;
  int t;
  double[] nz_;
  double[] zideal_;   // ideal point
  /**
  * Constructor.
  * Create a new SPEA2 instance
  * @param problem Problem to solve
  */
  public SPEA2_SDE_EGG(Problem problem,int run) {                
    super(problem) ;
    t=run;
  } // Spea2
   
  /**   
  * Runs of the Spea2 algorithm.
  * @return a <code>SolutionSet</code> that is a set of non dominated solutions
  * as a result of the algorithm execution  
   * @throws JMException 
  */  
  public SolutionSet execute() throws JMException, ClassNotFoundException {
    int populationSize, archiveSize, maxEvaluations, evaluations;
    Operator crossoverOperator, mutationOperator, selectionOperator;
    SolutionSet solutionSet, archive, offSpringSolutionSet;    
    
    //Read the params
    populationSize = ((Integer)getInputParameter("populationSize")).intValue();
    archiveSize    = ((Integer)getInputParameter("archiveSize")).intValue();
    maxEvaluations = ((Integer)getInputParameter("maxEvaluations")).intValue();
        
    //Read the operators
    crossoverOperator = operators_.get("crossover");
    mutationOperator  = operators_.get("mutation");
    selectionOperator = operators_.get("selection");        
        
    //Initialize the variables
    solutionSet  = new SolutionSet(populationSize);
    archive     = new SolutionSet(archiveSize);
    evaluations = 0;
        
    //-> Create the initial solutionSet
    Solution newSolution;
    for (int i = 0; i < populationSize; i++) {
      newSolution = new Solution(problem_);
      problem_.evaluate(newSolution);            
      problem_.evaluateConstraints(newSolution);
      evaluations++;
      solutionSet.add(newSolution);
    }                        
    initIdealPoint(solutionSet);  // initialize the ideal point
	initnzPoint(solutionSet);  
//	int generations_=1;
    while (evaluations < maxEvaluations){     
      SolutionSet union = ((SolutionSet)solutionSet).union(archive);
      Spea2_SDEFitness spea = new Spea2_SDEFitness(union);
      spea.fitnessAssign();
      archive = spea.environmentalSelection(archiveSize);     
//      archive.printObjectivesToFile("DRER_Analysis/DR0.5ER0.5/thetaDEA_"+problem_.getName()+"_"+problem_.getNumberOfObjectives()+"_T"+t+"_G"+(generations_+1));	
//      generations_++;
//      System.out.println(generations_);
      // Create a new offspringPopulation
      offSpringSolutionSet= new SolutionSet(populationSize);    
      Solution  [] parents = new Solution[3];
      for (int i = 0; i < archive.size(); i++){
      	Ranking ranking = new NondominatedRanking(archive);
  	    SolutionSet front = ranking.getSubfront(0);
  	    front.Suppress();
  	    updateIdealPoint(front);  // update the ideal point
  		updatenzPoint(front);
  	    SolutionSet KP = new SolutionSet();
  	    if(front.size()>problem_.getNumberOfObjectives())
  	    	KP = findingKneePoint(front);
  	    else//if(KP.size()<problem_.getNumberOfObjectives())
  	    	for (int ki=0;ki<archive.size();ki++)
  	    		KP.add(archive.get(ki));
  	    int r1,r2;
  		do {
  			r1 = PseudoRandom.randInt(0, KP.size() - 1);
  		} while (r1 == i);
  		do {
  			r2 = PseudoRandom.randInt(0, KP.size() - 1);
  		} while (r2 == i || r2==r1);
  		
  		parents[0] = archive.get(i);
  		parents[1] = KP.get(r1);
  		parents[2] = KP.get(r2);
          
        //make the crossover 
//        crossoverOperator.setParameter("P3", (double)evaluations/maxEvaluations);
        Solution offSpring = (Solution)crossoverOperator.execute(parents);           
        mutationOperator.execute(offSpring);            
        problem_.evaluate(offSpring);
        problem_.evaluateConstraints(offSpring);            
        offSpringSolutionSet.add(offSpring);
        evaluations++;
      } // while
      // End Create a offSpring solutionSet
      solutionSet = offSpringSolutionSet;                   
    } // while
        
    Ranking ranking = new NondominatedRanking(archive);
//    ranking.getSubfront(0).printObjectivesToFile("DRER_Analysis/DR0.2ER0.9/SPEA2_"+problem_.getName()+"_"+problem_.getNumberOfObjectives()+"_T"+t+"_G"+(generations_+1));
    return ranking.getSubfront(0);
  } // execute
  SolutionSet findingKneePoint(SolutionSet pop) {
		Matrix Pop_M, ones_M, ExtremePoint_M, Hyperplane, Distance_M;
		SolutionSet KneePoint = new SolutionSet();
		SolutionSet ExtremePoint = new SolutionSet(problem_.getNumberOfObjectives());
		SolutionSet TempPop = new SolutionSet(pop.size());
		for(int i=0;i<pop.size();i++) TempPop.add(pop.get(i));
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
			R[i] = (nz_[i]-zideal_[i])*r_;
		
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
void initIdealPoint(SolutionSet population_) {
		int obj = problem_.getNumberOfObjectives();
		zideal_ = new double[obj];
		for (int j = 0; j < obj; j++) {
			zideal_[j] = Double.MAX_VALUE;

			for (int i = 0; i < population_.size(); i++) {
				if (population_.get(i).getObjective(j) < zideal_[j])
					zideal_[j] = population_.get(i).getObjective(j);
			}
		}
	}
	
	
	void updateIdealPoint(SolutionSet pop){
		for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
			for (int i = 0; i < pop.size(); i++) {
				if (pop.get(i).getObjective(j) < zideal_[j])
					zideal_[j] = pop.get(i).getObjective(j);
			}
		}
	}
	void initnzPoint(SolutionSet population_) {
		int obj = problem_.getNumberOfObjectives();
		nz_ = new double[obj];
		for (int j = 0; j < obj; j++) {
			nz_[j] = Double.MIN_VALUE;

			for (int i = 0; i < population_.size(); i++) {
				if (population_.get(i).getObjective(j) > nz_[j])
					nz_[j] = population_.get(i).getObjective(j);
			}
		}
	}
	
	
	void updatenzPoint(SolutionSet pop){
		for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
			for (int i = 0; i < pop.size(); i++) {
				if (pop.get(i).getObjective(j) > nz_[j])
					nz_[j] = pop.get(i).getObjective(j);
			}
		}
	}
} // SPEA2
