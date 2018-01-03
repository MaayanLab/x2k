import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.Random;

public class ChEA {
	
	int[][] rank;
	private HashMap<String, HashSet<String>> allGenelists = new HashMap<String, HashSet<String>>();
	double pcut = 0.05;
	
	private HashSet<String> humanarchs4 = new HashSet<String>();
	private HashSet<String> mousearchs4 = new HashSet<String>();
	private HashSet<String> hg133 = new HashSet<String>();
	
	private LinkedList<TranscriptionFactor> tempRes = new LinkedList<TranscriptionFactor>();
	private LinkedList<TranscriptionFactor> transcriptionFactors = new LinkedList<TranscriptionFactor>();
	private LinkedList<TranscriptionFactor> transcriptionFactorsT = new LinkedList<TranscriptionFactor>();
	private LinkedList<TranscriptionFactor> transcriptionFactorsGB = new LinkedList<TranscriptionFactor>();
	
	private FisherExact fisherTest;
	private HashSet<String> geneBgSet = new HashSet<String>();
	private HashSet<String> geneBgSetMouse = new HashSet<String>();
	private HashSet<String> geneBgSetHuman = new HashSet<String>();
	private HashSet<String> geneBgSetBoth = new HashSet<String>();
	
	private HashSet<String> geneBgSetT = new HashSet<String>();
	private HashSet<String> geneBgSetMouseT = new HashSet<String>();
	private HashSet<String> geneBgSetHumanT = new HashSet<String>();
	private HashSet<String> geneBgSetBothT = new HashSet<String>();
	
	private HashSet<String> geneBgSetGB = new HashSet<String>();
	private HashSet<String> geneBgSetMouseGB = new HashSet<String>();
	private HashSet<String> geneBgSetHumanGB = new HashSet<String>();
	private HashSet<String> geneBgSetBothGB = new HashSet<String>();
	
	private String rankHumanChea = "data/ChEA2015/ChEA_ranks_human.tsv";
	private String rankMouseChea = "data/ChEA2015/ChEA_ranks_mouse.tsv";
	private String rankBothChea = "data/ChEA2015/ChEA_ranks.tsv";
	
	private final String rankTransfacMouse = "data/ChEA2015/mouse_TRANSFAC_ranks.txt";
	private final String rankTransfacHuman = "data/ChEA2015/human_TRANSFAC_ranks.txt";
	private final String rankTransfacBoth = "data/ChEA2015/combined_TRANSFAC_ranks.txt";
	
	private final String rankGBMouse = "data/ChEA2015/mouse_GB_ranks.txt";
	private final String rankGBHuman = "data/ChEA2015/human_GB_ranks.txt";
	private final String rankGBBoth = "data/ChEA2015/PWM-GB_ranks.txt";
	
	private HashSet<String> geneList;
	
	protected final static String CHEA_BACKGROUND = "data/ChEA2015/chea_background.csv";
	protected final static String TRANSFAC_BACKGROUND = "data/ChEA2015/transfac_background.csv";
	protected final static String PWM_GB_BACKGROUND = "data/ChEA2015/PWM-GB.csv";
	
	public static int port = 5000;
	
	// Output header
	protected final String HEADER = "TF,Target/Input,Targets/Database,Fraction/Input,Fraction/Database,Difference,P-value,Z-score,Combined Score,Genes";

	public static void main(String[] args) {
		ChEA chea = new ChEA();
		
		args = new String[]{"data/ChEA2015/chea_background.csv", "data/ChEA2015/combined_ChEA_ranks.txt", "data/testgmt", "data/output/output2.tsv", "both", "oddsratio"};
		
		String backgroundFile = args[0];
		String backgroundCorrection = args[1];
		String inputGeneList = args[2];
		String outputFile = args[3];		
		String species = "both";			// both, mouse, human
		String sorting = "oddsratio";	// oddsratio, combined-score
		
		if(args.length > 4){
			species = args[4];
		}
		if(args.length > 5){
			sorting = args[5];
		}
		
		chea.initialize(inputGeneList);
		
		System.out.println("CHEA server running at port: "+port);
		chea.startMultiServer();
		
	}
	
	// By default, load settings from file
	public ChEA() { }
	
	public void startMultiServer(){
		int portNumber = port;
		
		ServerSocket serverSocket = null;
        Socket socket = null;

        try {
            serverSocket = new ServerSocket(portNumber);
        } catch (IOException e) {
            e.printStackTrace();
        }
        while (true) {
            try {
                socket = serverSocket.accept();
            } catch (IOException e) {
                System.out.println("I/O error: " + e);
            }
            // new thread for a client
            new EchoThread(socket, this).start();
        }
	}
	
	public class EchoThread extends Thread {
	    protected Socket socket;
	    private ChEA cheaInstance;

	    public EchoThread(Socket clientSocket, ChEA _chea) {
	        this.socket = clientSocket;
	        cheaInstance = _chea;
	    }

	    public void run() {
	        try {
	            PrintWriter outs = new PrintWriter(socket.getOutputStream(), true);
				BufferedReader in = new BufferedReader(new InputStreamReader(socket.getInputStream()));

				String inputLine, outputLine;
				
				while ((inputLine = in.readLine()) != null) {
					outputLine = inputLine;
					//outs.println(outputLine);
					//System.out.println(outputLine);
					if (outputLine.equals("kill")){
						break;
					}
					else if(outputLine.startsWith("run")){
						String[] sp = outputLine.trim().split(",");
						//System.out.println("computing -------------");
						String tfres = cheaInstance.run(sp[1], sp[2], sp[3], sp[4], sp[5]);
						
						outs.println(tfres+"messageComplete");
					}
					else{
						
					}
				}   
	        }
	        catch(Exception e){
	        	e.printStackTrace();
	        }
	    }
	}
	
	public void readGMT(String _path){
		allGenelists = new HashMap<String, HashSet<String>>();
		
		File folder = new File(_path);
		File[] listOfFiles = folder.listFiles();
		
		for(File f : listOfFiles){
			try{
				BufferedReader br = new BufferedReader(new FileReader(f));
				String line = "";
				
				while((line = br.readLine()) != null){
					String[] sp = line.trim().toUpperCase().split("\t");
					String g1 = f.getName()+";"+sp[0];

					HashSet<String> genes = new HashSet<String>();
					for(int j=2; j<sp.length; j++){
						genes.add(sp[j]);
					}
					allGenelists.put(g1, genes);
				}
				br.close();
			}
			catch(Exception e){
				e.printStackTrace();
			}
		}
	}
	
	public void initialize(String _geneList){
		
		readGMT(_geneList);
		fisherTest = new FisherExact(60000);
		
		//geneList = readGenes(_geneList);
		//System.out.println(geneList.size());
		
		humanarchs4 = readGenes("data/input/archs4humangenes.txt");
		mousearchs4 = readGenes("data/input/archs4mousegenes.txt");
		hg133 = readGenes("data/input/hg133genes.txt");
		
		long time = System.currentTimeMillis();
		readBackground(CHEA_BACKGROUND, rankBothChea, transcriptionFactors, geneBgSet);
		readBackground(TRANSFAC_BACKGROUND, rankTransfacBoth, transcriptionFactorsT, geneBgSetT);
		readBackground(PWM_GB_BACKGROUND, rankGBBoth, transcriptionFactorsGB, geneBgSetGB);
		//System.out.println("Read files: "+(System.currentTimeMillis() - time));
	}
	
	// Run for file names
	public String run(String _sorting, String _species, String _type, String _geneBackground, String _tfLength) {
		
		int tfLength = Integer.parseInt(_tfLength);
		
		HashSet<String> chooseBackground;
		if(_geneBackground.equals("humanarchs4")){
			chooseBackground = humanarchs4;
		}
		else if(_geneBackground.equals("mousearchs4")){
			chooseBackground = mousearchs4;
		}
		else if(_geneBackground.equals("hg133")){
			chooseBackground = hg133;
		}
		else{
			chooseBackground = geneBgSet;
		}
		
		String output = "";
		if(_type.equals("chea")){
			output = computeEnrichment(allGenelists, _species, _sorting, transcriptionFactors, geneBgSet, chooseBackground, tfLength);
		}
		else if(_type.equals("transfac")){
			//System.out.println("Tranfac");
			output = computeEnrichment(allGenelists, _species, _sorting, transcriptionFactorsT, geneBgSetT, chooseBackground, tfLength);
		}
		else{
			output = computeEnrichment(allGenelists, _species, _sorting, transcriptionFactorsGB, geneBgSetGB, chooseBackground, tfLength);
		}
		
		return output;
	}
	
	public void writeFile(String _filename) {
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(_filename)));
			for(TranscriptionFactor tf : transcriptionFactors){
				if(tf.enrichedTargets.size() > 0){
					bw.write(tf.toString()+"\n");
				}
			}
			bw.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	// Used to pass information about transcription factors
	public Collection<TranscriptionFactor> getTopRanked(int ranks) {
		LinkedHashSet<TranscriptionFactor> topRanked = new LinkedHashSet<TranscriptionFactor>(ranks);
		
		Iterator<TranscriptionFactor> itr = transcriptionFactors.iterator();
		while (itr.hasNext() && topRanked.size() < ranks){
			topRanked.add(itr.next());
		}
		return topRanked;
	}
	
	// Used to print lists
	public Collection<String> getTopRankedList(int ranks) {
		LinkedHashSet<String> topRanked = new LinkedHashSet<String>(ranks);				
		
		Iterator<TranscriptionFactor> itr = transcriptionFactors.iterator();
		while (itr.hasNext() && topRanked.size() < ranks){
			topRanked.add(itr.next().getSimpleName());
		}
		
		return topRanked;
	}
	
	public Collection<TranscriptionFactor> getRankedList() {
		return transcriptionFactors;
	}
	
	private HashSet<String> readGenes(String _geneListFile){
		HashSet<String> genes = new HashSet<String>();
		try{
			BufferedReader br = new BufferedReader(new FileReader(new File(_geneListFile)));
			String line = "";
			
			while((line = br.readLine()) != null){
				genes.add(line.trim().toUpperCase());
			}
			
			br.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
		return genes;
	}
	
	protected void readBackground(String _background, String _rank, LinkedList<TranscriptionFactor> _transcriptionFactors, HashSet<String> _geneBgSet) {
		
		HashMap<String, TranscriptionFactor> tfMap = new HashMap<String, TranscriptionFactor>();
		
		try{
			BufferedReader br = new BufferedReader(new FileReader(new File(_background)));
			String line = "";
			
			while((line = br.readLine()) != null){
				String[] sp = line.toUpperCase().trim().split(",");
				String simpleName = sp[1];
				String name = sp[2];			
				String target = sp[3];
				String species = sp[7];
			
				_geneBgSet.add(target);
				
				if (tfMap.containsKey(name)) {
					TranscriptionFactor tf = tfMap.get(name);
					tf.addTarget(target);
					//tf.setSpecies(species);
				}
				else {
					tfMap.put(name, new TranscriptionFactor(name, simpleName, species, target));
				}
			}
			br.close();
			
			br = new BufferedReader(new FileReader(new File(_rank)));
			line = "";

			while((line = br.readLine()) != null){
				String[] sp = line.toUpperCase().split("\t");
				if(sp.length == 3 && tfMap.keySet().contains(sp[0])){
					tfMap.get(sp[0]).setRankStats(Double.parseDouble(sp[1]), Double.parseDouble(sp[2]));
				}
			}
			br.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
		
		_transcriptionFactors.addAll(tfMap.values());
	}
	
	private String computeEnrichment(HashMap<String, HashSet<String>> _genes, String _species, String _sort, LinkedList<TranscriptionFactor> _tfs, HashSet<String> _geneBgSet, HashSet<String> _geneBackground, int _tfLength) {
		
		// filter genes from input list that are not associated with an upstream transcription factor
		String output = "";
		int counterk = 0;
		
		HashSet<String> genelist = new HashSet<String>(_geneBackground);
		genelist.retainAll(_geneBgSet);
		
		HashMap<String, HashSet<String>> tftemp = new HashMap<String, HashSet<String>>();
		for(TranscriptionFactor tf : _tfs){
			HashSet<String> tfGenes = new HashSet<String>(tf.getTargets());
			tfGenes.retainAll(genelist);
			tftemp.put(tf.getSimpleName(), tfGenes);
		}
		
		//System.out.println("Background: "+genelist.size());
		
		for(String key : _genes.keySet()){
			counterk++;
			HashSet<String> geneInputSet = new HashSet<String>(_genes.get(key));
			geneInputSet.retainAll(genelist);
			
			long time = System.currentTimeMillis();
			double maxOver = 0;
			double maxOverlap = 0;
			double pval = 1;
			
			for(TranscriptionFactor tf : _tfs){
				// Target input genes is the intersection of target background genes and input genes
				HashSet<String> overlapGenes = new HashSet<String>(geneInputSet);
				overlapGenes.retainAll(tftemp.get(tf.getSimpleName()));
				
				int numTF = tftemp.get(tf.getSimpleName()).size();
				int totalBgGenes = genelist.size();
				int totalInputGenes = geneInputSet.size();
				int numOverlap = overlapGenes.size();
				double oddsRatio = (numOverlap*1.0/(totalInputGenes - numOverlap))/(numTF*1.0/(totalBgGenes - numTF));
				
				if (numOverlap > 0) {
					if(oddsRatio > maxOver) {
						maxOver = oddsRatio;
						maxOverlap = numOverlap;
					}
					//System.out.println(oddsRatio+" - "+numOverlap+" - "+(totalInputGenes - numOverlap)+" - "+numTF+" - "+(totalBgGenes - numTF));
					double pvalue = fisherTest.getRightTailedP(numOverlap,(totalInputGenes - numOverlap), numTF, (totalBgGenes - numTF));
					pval = Math.min(pvalue, pval);
					tf.setEnrichedTargets(overlapGenes);
					tf.setFractionOfTargetsInInput(numOverlap*1.0/totalInputGenes);
					tf.setFractionOfTargetsInBackground(numTF*1.0/totalBgGenes);
					tf.setPValue(pvalue);
					tf.setOddsRatio(oddsRatio);
					tf.setBackgroundIntersect(numTF);
				}
				else {
					tf.setPValue(1);
					tf.setEnrichedTargets(new HashSet<String>());
				}
			}
			
			// First, sort by p-value
			Collections.sort(_tfs);
			
			// Count current rank and compute z-score
			int counter = 1;
			
			for (TranscriptionFactor tf : _tfs) {
				if(_sort.equals("combined_score") || _sort.equals("rank")){
					tf.computeScore(counter, false);
				}
				else if(_sort.equals("oddsratio")){
					tf.computeScore(counter, true);
				}
				counter++;
			}
			
			if (_sort.equals("combined_score") || _sort.equals("oddsratio")) {
				// Sort by combined score
				Collections.sort(_tfs, new Comparator<TranscriptionFactor>() {
					@Override
					public int compare(TranscriptionFactor o1, TranscriptionFactor o2) {
						if (o1.getCombinedScore() < o2.getCombinedScore())				
							return 1;
						else if (o1.getCombinedScore() > o2.getCombinedScore())
							return -1;
						else
							return 0;
					}
				});
			}
			else if (_sort.equals("rank")) {
				// Sort by z-score
				Collections.sort(_tfs, new Comparator<TranscriptionFactor>() {
					@Override
					public int compare(TranscriptionFactor o1, TranscriptionFactor o2) {
						if (o1.getZScore() > o2.getZScore())				
							return 1;
						else if (o1.getZScore() < o2.getZScore())
							return -1;
						else
							return 0;
					}
				});
			}
			
			StringBuffer sb = new StringBuffer();
			sb.append(key).append(",");
			int buff = 0;
			for(int i=0; i<Math.min(_tfs.size(), _tfLength+buff); i++){
				String tfname = _tfs.get(i).getSimpleName();
				//String tfp = ""+_tfs.get(i).getPValue();
				//String tfodds = ""+_tfs.get(i).getOddsRatio();
				if(_tfs.get(i).getPValue() < 0.05){
					//sb.append(tfname).append(",").append(tfp).append(",").append(tfodds).append(";");
					sb.append(tfname).append(",");
				}
				else{
					buff++;
				}
			}
			
			System.out.println(_tfs.size()+"-"+buff);
			
			output += sb.toString()+"\n";
			System.out.println(counterk+" - "+key+" - "+geneInputSet.size()+" - "+pval+" - "+maxOver+" - "+maxOverlap+" - "+(System.currentTimeMillis() - time));
		}
		
		output = output.replace(",\n", "\n");
		//System.out.println(output);
		
		return output;
	}
	
	private void doRandom(){
	
		HashSet<String> genelist = readGenes("data/input/archs4humangenes.txt");
		genelist.retainAll(geneBgSet);
		String[] genelistarr = genelist.toArray(new String[0]);
		//System.out.println("Sized: "+genelist.size());
		
		Random rn = new Random();
		
		HashMap<String, int[]> rankMap = new HashMap<String, int[]>();
		for (int i=0; i<transcriptionFactors.size(); i++) {
			rankMap.put(transcriptionFactors.get(i).getName(), new int[10000]);
		}
		
		for(int i=0; i<10000; i++){
			HashSet<String> randomGenes = new HashSet<String>();
			while(randomGenes.size() < 1000){
				randomGenes.add(genelistarr[rn.nextInt(genelistarr.length)]);
			}
		
			//computeEnrichment(randomGenes, "both", "pvalue", transcriptionFactors, geneBgSet, humanarchs4);
			
			for (int j=0; j<transcriptionFactors.size(); j++) {
				rankMap.get(transcriptionFactors.get(j).getName())[i] = j;
			}
		}
		
		//System.out.println(rankMap.keySet().toString());
		//System.out.println(Arrays.toString(rankMap.get("THAP11-20581084")));
		
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File("data/output/rank_human.txt")));

			for(int j=0; j<transcriptionFactors.size(); j++){
				int[] tt = rankMap.get(transcriptionFactors.get(j).getName());
				bw.write(transcriptionFactors.get(j).getName());
				for(int i=0; i<10000; i++){
					bw.write("\t"+tt[i]);
				}
				bw.write("\n");
			}

			bw.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
}