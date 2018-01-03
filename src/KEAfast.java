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

public class KEAfast {
	
	private FisherExact fisherTest;
	int[][] rank;
	double pcut = 0.05;
	
	
	private HashMap<String, HashSet<String>> allGenelists = new HashMap<String, HashSet<String>>();
	
	private HashSet<String> humanarchs4 = new HashSet<String>();
	private HashSet<String> mousearchs4 = new HashSet<String>();
	private HashSet<String> hg133 = new HashSet<String>();
	
	private LinkedList<TranscriptionFactor> kinasesP = new LinkedList<TranscriptionFactor>();
	private LinkedList<TranscriptionFactor> kinasesKP = new LinkedList<TranscriptionFactor>();
	private LinkedList<TranscriptionFactor> tempRes = new LinkedList<TranscriptionFactor>();
	
	
	private HashSet<String> geneBgSetP = new HashSet<String>();
	private HashSet<String> geneBgSetKP = new HashSet<String>();

	private String ranks = "data/KEA/kea_ranks.txt";
	
	//private HashSet<String> geneList;
	
	protected final static String KP_BACKGROUND = "data/KEA/kinase-protein_interactions.csv";
	protected final static String PHOSPHO_BACKGROUND = "data/KEA/phosphorylation_reactions.csv";
	
	static int port = 5002; 
	
	public static void main(String[] args) {
		KEAfast kea = new KEAfast();
		//long time = System.currentTimeMillis();
		kea.initialize("data/testgmt");
		
		System.out.println("KEA server running at port: "+port);
		kea.startMultiServer();
		
		
		//System.out.println(System.currentTimeMillis() - time);
	}
	
	// By default, load settings from file
	public KEAfast() { }
	
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
	    private KEAfast keaInstance;

	    public EchoThread(Socket clientSocket, KEAfast _kea) {
	        this.socket = clientSocket;
	        keaInstance = _kea;
	    }

	    public void run() {
	        try {
	            PrintWriter outs = new PrintWriter(socket.getOutputStream(), true);
				BufferedReader in = new BufferedReader(new InputStreamReader(socket.getInputStream()));
				
				String inputLine, outputLine;
				
				while ((inputLine = in.readLine()) != null) {
					outputLine = inputLine;
					
					//System.out.println(outputLine);
					if (outputLine.equals("kill")){
						break;
					}
					else if(outputLine.startsWith("run")){
						String[] sp = outputLine.trim().split(",");
						//System.out.println("computing -------------");
						String tfres = keaInstance.run(sp[1], sp[2], sp[3], sp[4]);
						//System.out.println("KEA complete");
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
		
		humanarchs4 = readGenes("data/input/archs4humangenes.txt");
		mousearchs4 = readGenes("data/input/archs4mousegenes.txt");
		hg133 = readGenes("data/input/hg133genes.txt");
		
		long time = System.currentTimeMillis();
		readBackground(KP_BACKGROUND, ranks, kinasesKP, geneBgSetKP);
		readBackground(PHOSPHO_BACKGROUND, ranks, kinasesP, geneBgSetP);

		//System.out.println("Read files: "+(System.currentTimeMillis() - time));
	}
	
	// Run for file names
	public String run(String _sorting, String _type, String _geneBackground, String _tfLength) {
		
		int tfLength = Integer.parseInt(_tfLength);
		
		HashSet<String> chooseBackground;
		if(_geneBackground.equals("humanarchs4")){
			chooseBackground = humanarchs4;
		}
		else if(_geneBackground.equals("humanarchs4")){
			chooseBackground = mousearchs4;
		}
		else if(_geneBackground.equals("hg133")){
			chooseBackground = mousearchs4;
		}
		else{
			chooseBackground = geneBgSetKP;
		}
		
		String output = "";
		if(_type.equals("KP")){
			output = computeEnrichment(allGenelists, _sorting, kinasesKP, geneBgSetKP, chooseBackground, tfLength);
		}
		else if(_type.equals("P")){
			output = computeEnrichment(allGenelists, _sorting, kinasesP, geneBgSetP, chooseBackground, tfLength);
		}
		
		return output;
	}
	
	// Used to pass information about transcription factors
	public Collection<TranscriptionFactor> getTopRanked(int ranks) {
		LinkedHashSet<TranscriptionFactor> topRanked = new LinkedHashSet<TranscriptionFactor>(ranks);
		
		Iterator<TranscriptionFactor> itr = tempRes.iterator();
		while (itr.hasNext() && topRanked.size() < ranks){
			topRanked.add(itr.next());
		}
		return topRanked;
	}
	
	// Used to print lists
	public Collection<String> getTopRankedList(int ranks) {
		LinkedHashSet<String> topRanked = new LinkedHashSet<String>(ranks);				
		
		Iterator<TranscriptionFactor> itr = tempRes.iterator();
		while (itr.hasNext() && topRanked.size() < ranks){
			topRanked.add(itr.next().getSimpleName());
		}
		
		return topRanked;
	}
	
	public Collection<TranscriptionFactor> getRankedList() {
		return tempRes;
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
				String simpleName = sp[2];
				String name = sp[2];			
				String target = sp[3];
				String species = "egal";
			
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
	
	private String computeEnrichment(HashMap<String, HashSet<String>> _genes, String _sort, LinkedList<TranscriptionFactor> _tfs, HashSet<String> _geneBgSet, HashSet<String> _geneBackground, int _tfLength) {
		
		// filter genes from input list that are not associated with an upstream transcription factor
		String output = "";

		HashSet<String> genelist = new HashSet<String>(_geneBackground);
		genelist.retainAll(_geneBgSet);
		
		HashMap<String, HashSet<String>> tftemp = new HashMap<String, HashSet<String>>();
		for(TranscriptionFactor tf : _tfs){
			HashSet<String> tfGenes = new HashSet<String>(tf.getTargets());
			tfGenes.retainAll(genelist);
			tftemp.put(tf.getSimpleName(), tfGenes);
		}
		
		System.out.println("Background: "+genelist.size());
		
		for(String key : _genes.keySet()){

			HashSet<String> geneInputSet = new HashSet<String>(_genes.get(key));
			geneInputSet.retainAll(genelist);
			
			double maxOver = 0;
			
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
					}
					
					double pvalue = fisherTest.getRightTailedP(numOverlap,(totalInputGenes - numOverlap), numTF, (totalBgGenes - numTF));
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
				if(_sort.equals("oddsratio")){
					tf.computeScore(counter, true);
				}
				else if(_sort.equals("combined_score") || _sort.equals("rank")){
					tf.computeScore(counter, false);
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
				String tfp = ""+_tfs.get(i).getPValue();
				String tfodds = ""+_tfs.get(i).getOddsRatio();
				if(_tfs.get(i).getPValue() < pcut){
					//sb.append(tfname).append(",").append(tfp).append(",").append(tfodds).append(";");
					sb.append(tfname).append(",");
				}
				else{
					buff++;
				}
			}
			output += sb.toString()+"\n";
		}
		
		output = output.replace(",\n", "\n");
		System.out.println(output);
		
		return output;
	}
	
//	private void doRandom(){
//	
//		HashSet<String> genelist = readGenes("data/input/archs4humangenes.txt");
//		genelist.retainAll(geneBgSet);
//		String[] genelistarr = genelist.toArray(new String[0]);
//		System.out.println("Sized: "+genelist.size());
//		
//		Random rn = new Random();
//		
//		HashMap<String, int[]> rankMap = new HashMap<String, int[]>();
//		for (int i=0; i<transcriptionFactors.size(); i++) {
//			rankMap.put(transcriptionFactors.get(i).getName(), new int[10000]);
//		}
//		
//		for(int i=0; i<10000; i++){
//			HashSet<String> randomGenes = new HashSet<String>();
//			while(randomGenes.size() < 1000){
//				randomGenes.add(genelistarr[rn.nextInt(genelistarr.length)]);
//			}
//		
//			//computeEnrichment(randomGenes, "both", "pvalue", transcriptionFactors, geneBgSet, humanarchs4);
//			
//			for (int j=0; j<transcriptionFactors.size(); j++) {
//				rankMap.get(transcriptionFactors.get(j).getName())[i] = j;
//			}
//		}
//		
//		System.out.println(rankMap.keySet().toString());
//		System.out.println(Arrays.toString(rankMap.get("THAP11-20581084")));
//		
//		try{
//			BufferedWriter bw = new BufferedWriter(new FileWriter(new File("data/output/rank_human.txt")));
//
//			for(int j=0; j<transcriptionFactors.size(); j++){
//				int[] tt = rankMap.get(transcriptionFactors.get(j).getName());
//				bw.write(transcriptionFactors.get(j).getName());
//				for(int i=0; i<10000; i++){
//					bw.write("\t"+tt[i]);
//				}
//				bw.write("\n");
//			}
//
//			bw.close();
//		}
//		catch(Exception e){
//			e.printStackTrace();
//		}
//	}
}