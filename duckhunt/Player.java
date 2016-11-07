import java.util.*;

class HitBird{
	public int direction;
	
	public int birdNumber;
	
	public double directionProb;

}

class Player {
   
	HashMap<Integer, ArrayList<HMM>> speciesHMM; 
	//ArrayList<HMM> hmm;
	private int timeStep;
    public Player() {
    	speciesHMM = new HashMap<Integer, ArrayList<HMM>>();
        for (int i = 0; i < Constants.COUNT_SPECIES; i++) {
            speciesHMM.put(i, new ArrayList<HMM>());
        }
        timeStep = 0;
    }
    ArrayList<HitBird> hitBirdList;


    /**
     * Shoot!
     *
     * This is the function where you start your work.
     *
     * You will receive a variable pState, which contains information about all
     * birds, both dead and alive. Each bird contains all past moves.
     *
     * The state also contains the scores for all players and the number of
     * time steps elapsed since the last time this function was called.
     *
     * @param pState the GameState object with observations etc
     * @param pDue time before which we must have returned
     * @return the prediction of a bird we want to shoot at, or cDontShoot to pass
     */
    public Action shoot(GameState pState, Deadline pDue) {
        
    	int numHMM = pState.getNumBirds();
    	
    		
        	if(timeStep == 0){
//        		hmm = new ArrayList<HMM>();
//        		//initialize each HMM model for each bird        		
//            	for(int i=0;i<numHMM;i++){
//            		hmm.add(new HMM());
//            	}
            	timeStep = timeStep + 1;
            	return cDontShoot;
            	
        	}
        	else{
                double max = 0;
                int hitDirection = -1;
                int hitBird = -1;
                
        		if (pState.getRound() == 0)
        			return cDontShoot;
        		hitBirdList= new ArrayList<HitBird>();;
        		//start for the last x timestaps
        		int x= 100-numHMM*5;
        		
        		if(timeStep >= x){
         		   //hit

                   for(int i=0;i<numHMM;i++){
        			   
                	  
        			   // get past sequence(seq) for bird i
        			   Bird bird = pState.getBird(i);
        			   if (!bird.isAlive())
        				   continue;
        			   
        		       ArrayList<Integer> obArrayList=new ArrayList<Integer>();
        		       
        		       for (int j = 0; j < bird.getSeqLength(); j++) {
        		           if (bird.wasAlive(j)) 
        		        	   obArrayList.add(bird.getObservation(j));
        		       }
        		       int[] seq = new int[obArrayList.size()];

        		       for(int j=0;j<obArrayList.size();j++){
        		    	   seq[j]=(Integer)obArrayList.get(j);
        		       }  		           		           	     		    	 
     		    	  //1.compute the species prob (for one bird): bestProb
 		 		      double bestSpeciesProb = 0.0;
 		 		      
 		 		      
 		 		      
 		              int species = -1;//unknown
 		              for (int j = 0; j < Constants.COUNT_SPECIES; j++) {
 		          	    ArrayList<HMM> hmmList = speciesHMM.get(j);
 		                  for (HMM hmm : hmmList) {
 		                      if (hmm != null) {
 		                         double prob = Math.abs(hmm.Probability_HMM_Obsequence(seq));		                        
 		                          if (prob < Integer.MAX_VALUE && prob > bestSpeciesProb) {
 		                        	 bestSpeciesProb = prob;
 		                             species = j;
 		                          }
 		                      }
 		                  }
 		              }
 		             if(species == 5)
 		            	 continue;
 		             
 		         // get bird i 's hmms box:list
 		            ArrayList<HMM> list = speciesHMM.get(species);

 		            int[] countdir = new int[9];
 		            double[] countprob = new double[9];
 		            for (int p = 0; p < 9; ++p)
 		            {
 		            	countdir[p] = 0;
 		            	countprob[p] = 0;
 		            }
 		            double seqProbsx[] = new double[9];
 		            
 		            for(HMM hmm:list){
 		            	seqProbsx = hmm.predict(seq);
 		            	double tempMax = 0.0;
 		            	int tempDirect = -1;
 		            	for(int y=0;y<9;y++){
 		            		if(seqProbsx[y]>tempMax){
 		            			tempMax=seqProbsx[y];
 		            			tempDirect = y;
 		            		}
 		            		
 		            	}
 		            	if (tempMax > 0)
 		            	{
 		            		countdir[tempDirect]++;
 		            		countprob[tempDirect] += tempMax;
 		            	}
 		            }
 		            
 		            int maxdir = -1;
 		            int maxcount = -1;
 		            for (int r = 0; r < 9; ++r)
 		            {
 		            	if (countdir[r] > maxcount)
 		            	{
 		            		maxdir = r;
 		            		maxcount = countdir[r];
 		            	}
 		            }
// 		            System.err.println(i + " " + maxdir + " " +
// 		            		(double)maxcount / list.size() + " " + countprob[maxdir] / maxcount);
 		            if (maxcount > list.size() * 0.5675)
 		            {
 		            	if (countprob[maxdir] > 0.7018 * maxcount)
 		            	{
 		            		HitBird hb = new HitBird();
 		            		hb.birdNumber = i;
 		            		hb.directionProb = countprob[maxdir] / maxcount;
 		            		hb.direction = maxdir;
 		            		hitBirdList.add(hb);
 		            	}
 		            }
 		            

      		        
             		                   
 		              
        		   }//--end for all birds

                   for(HitBird hb :hitBirdList){
                	   if(hb.directionProb > max){
                		   max = hb.directionProb;
                		   hitBird = hb.birdNumber;
                		   hitDirection= hb.direction;
                	   }
                   }
                  
                   
         		   
         		   
         		   
         		 
         	    }//--end if from x timestap(train and hit);
        		timeStep = timeStep + 1;
        		if (max != 0.0)
        		{
        			System.err.println("Shoot!~!~!" + hitBird + " " + hitDirection);
        			return new Action(hitBird, hitDirection);
        		}
        		return cDontShoot;
        		
        		
        	}
        	
        	//return cDontShoot;
		  
        		
                  
		    
		  
        	
        
    }

    /**
     * Guess the species!
     * This function will be called at the end of each round, to give you
     * a chance to identify the species of the birds for extra points.
     *
     * Fill the vector with guesses for the all birds.
     * Use SPECIES_UNKNOWN to avoid guessing.
     *
     * @param pState the GameState object with observations etc
     * @param pDue time before which we must have returned
     * @return a vector with guesses for all the birds
     */
    public int[] guess(GameState pState, Deadline pDue) {
    	//System.err.printf("Guess Time!"+'\n');
    	
    	timeStep= 0;
    	int[] lGuess = new int[pState.getNumBirds()];
    	//Guess from 2nd round.  	
        if(pState.getRound()>0){
        	//System.err.printf("this is round "+pState.getRound()+'\n');
            
        	for (int i = 0; i < pState.getNumBirds(); ++i){
            	Bird bird = pState.getBird(i);
            	ArrayList<Integer> seq2 = new ArrayList<Integer>();
    		    for (int j = 0; j < bird.getSeqLength(); j++) {
    		    	if (bird.wasAlive(j)) 
    		    		seq2.add(bird.getObservation(j));
    		    }
    		    int[] seq=new int[seq2.size()];//ArrarList !!!!!!
    		    for(int s=0;s<seq2.size();s++){
    		    	seq[s]=(int) seq2.get(s);
    		    }
  		      
  		        double bestProb = 0.0;
  		      
  		        int species = -1;//unknown
  		        for (int j = 0; j < Constants.COUNT_SPECIES; j++) {
  		        	ArrayList<HMM> hmmList = speciesHMM.get(j);
  		        	for (HMM hmm : hmmList) {
  		        		if (hmm != null) {
  		        			double prob = Math.abs(hmm.Probability_HMM_Obsequence(seq));
  		        			if (prob > bestProb) {
  		        				bestProb = prob;
  		        				species = j;
  		        			}
  		        		}
                      
  		        	}
                  
  		        }
  		        lGuess[i] = species;
              
            }
            
            return lGuess;
        }else{
        	for(int i=0;i<pState.getNumBirds();i++){
        		lGuess[i]=Constants.SPECIES_RAVEN;
        	}
        	return lGuess;
        }

    }

    /**
     * If you hit the bird you were trying to shoot, you will be notified
     * through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pBird the bird you hit
     * @param pDue time before which we must have returned
     */
    public void hit(GameState pState, int pBird, Deadline pDue) {
        System.err.println("HIT BIRD!!!");
    }

    /**
     * If you made any guesses, you will find out the true species of those
     * birds through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pSpecies the vector with species
     * @param pDue time before which we must have returned
     */
    /**
     * @param pState
     * @param pSpecies
     * @param pDue
     */
    public void reveal(GameState pState, int[] pSpecies, Deadline pDue) {
    	timeStep = 0;
    	for(int i=0;i<pSpecies.length;i++){
			   
    		// get past sequence(seq) for bird i
			Bird bird = pState.getBird(i);
			   
		    ArrayList<Integer> seq2 = new ArrayList<Integer>();
		    for (int j = 0; j < bird.getSeqLength(); j++) {
		    	if (bird.wasAlive(j)) 
		    		seq2.add(bird.getObservation(j));
		    }
		    int[] seq=new int[seq2.size()];//ArrarList !!!!!!
		    for(int s=0;s<seq2.size();s++){
		    	seq[s]=(int) seq2.get(s);
		    }
            // 2.train each model  train() method using observation sequence for parameter. 
   		    HMM hmm=new HMM();
   		    hmm.train(seq);
   		    speciesHMM.get(pSpecies[i]).add(hmm);
		}
    	
    	
    }

    public static final Action cDontShoot = new Action(-1, -1);
    
    
}
