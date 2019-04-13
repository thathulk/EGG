package jmetal.metaheuristics.NSGAIII;

import Jama.Matrix;
import jmetal.core.*;

import jmetal.util.*;
import jmetal.util.comparators.ObjectiveComparator;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;
import jmetal.util.vector.VectorGenerator;
import jmetal.util.vector.TwoLevelWeightVectorGenerator;

public class NSGAIII_EGG extends Algorithm {

	private int populationSize_;

	private int div1_;
	private int div2_;

	private SolutionSet population_;
	SolutionSet offspringPopulation_;
	SolutionSet union_;
	SolutionSet KP;
	
	int evaluations_;
	double r_;
	double t_;
	
	Operator crossover_;
	Operator mutation_;
	Operator selection_;

	double[][] lambda_; // reference points
	
	boolean normalize_; // do normalization or not


	public NSGAIII_EGG(Problem problem) {
		super(problem);
	} // NSGAII

	public SolutionSet execute() throws JMException, ClassNotFoundException {
		int maxEvaluations_;

		evaluations_ = 0;
		r_ = 1.0;
	    t_ = 0.0;
	    maxEvaluations_ = ((Integer) this.getInputParameter("maxEvaluations"))
				.intValue();

		div1_ = ((Integer) this.getInputParameter("div1")).intValue();
		div2_ = ((Integer) this.getInputParameter("div2")).intValue();
		
		
		normalize_ = ((Boolean) this.getInputParameter("normalize")).booleanValue();

		VectorGenerator vg = new TwoLevelWeightVectorGenerator(div1_, div2_,problem_.getNumberOfObjectives());
		lambda_ = vg.getVectors();

		populationSize_ = vg.getVectors().length;
		if (populationSize_ % 2 != 0)
			populationSize_ += 1;


		mutation_ = operators_.get("mutation");
		crossover_ = operators_.get("crossover");
		selection_ = operators_.get("selection");

		initPopulation();
		while (evaluations_ < maxEvaluations_) {
//			offspringPopulation_ = new SolutionSet(populationSize_);
//			Solution[] parents = new Solution[2];
//			for (int i = 0; i < (populationSize_ / 2); i++) {
//				if (generations_ < maxGenerations_) {
//					// obtain parents
//
//					parents = (Solution[]) selection_.execute(population_);
//
//					Solution[] offSpring = (Solution[]) crossover_
//							.execute(parents);
//
//					mutation_.execute(offSpring[0]);
//					mutation_.execute(offSpring[1]);
//
//					problem_.evaluate(offSpring[0]);
//					problem_.evaluateConstraints(offSpring[0]);
//					problem_.evaluate(offSpring[1]);
//					problem_.evaluateConstraints(offSpring[1]);
//
//					offspringPopulation_.add(offSpring[0]);
//					offspringPopulation_.add(offSpring[1]);
//
//				} // if
//			} // for
			Ranking ranking = new NondominatedRanking(population_);
		    SolutionSet front = ranking.getSubfront(0);
		    front.Suppress();
		    KP = new SolutionSet();
		    if(front.size()>problem_.getNumberOfObjectives())
		    	KP = findingKneePoint(front);
		    else//if(KP.size()<problem_.getNumberOfObjectives())
		    	for (int ki=0;ki<population_.size();ki++)
		    		KP.add(population_.get(ki));
			createOffSpringPopulation();
			union_ = ((SolutionSet) population_).union(offspringPopulation_);

			// Ranking the union
			ranking = new NondominatedRanking(union_);

			int remain = populationSize_;
			int index = 0;
			front = null;
			population_.clear();

			// Obtain the next front
			front = ranking.getSubfront(index);
			while ((remain > 0) && (remain >= front.size())) {

				for (int k = 0; k < front.size(); k++) {
					population_.add(front.get(k));
				} // for

				// Decrement remain
				remain = remain - front.size();

				// Obtain the next front
				index++;
				if (remain > 0) {
					front = ranking.getSubfront(index);
				} // if
			}

			if (remain > 0) { // front contains individuals to insert

				new Niching(population_, front, lambda_, remain, normalize_)
						.execute();
				remain = 0;
			}
			
			evaluations_=evaluations_+populationSize_;

		}

		Ranking ranking = new NondominatedRanking(population_);
		return ranking.getSubfront(0);

	}

	public void initPopulation() throws JMException, ClassNotFoundException {

		population_ = new SolutionSet(populationSize_);

		for (int i = 0; i < populationSize_; i++) {
			Solution newSolution = new Solution(problem_);

			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);

			population_.add(newSolution);
		} // for
	} // initPopulation
	
	void createOffSpringPopulation() throws JMException {
		offspringPopulation_ = new SolutionSet(populationSize_);

		for (int i = 0; i < populationSize_; i++) 
			doCrossover(i);
	}
	
	
	void doCrossover(int i) throws JMException{
		int r1,r2;
		do {
			r1 = PseudoRandom.randInt(0, KP.size() - 1);
//			r1 = PseudoRandom.randInt(0, population_.size() - 1);
		} while (r1 == i);
		do {
			r2 = PseudoRandom.randInt(0, KP.size() - 1);
//			r2 = PseudoRandom.randInt(0, population_.size() - 1);
		} while (r2 == i || r2==r1);
		Solution[] parents = new Solution[3];
		
		parents[0] = population_.get(i);
		parents[1] = KP.get(r1);
		parents[2] = KP.get(r2);
//		parents[1] = population_.get(r1);
//		parents[2] = population_.get(r2);
//		crossover_.setParameter("P3", (double)generations_/maxGenerations);
		Solution offSpring = (Solution) crossover_.execute(parents);
		
		mutation_.execute(offSpring);
		
		problem_.evaluate(offSpring);
		
		
		offspringPopulation_.add(offSpring);
	}
	
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
} // NSGA-III

