import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;

public class FastGenes2Networks {

	HashMap<String, String> networkList = new HashMap<String, String>();
	HashMap<String, HashMap<String, HashSet<String>>> networks;
	HashMap<String, HashSet<String>> allGenelists = new HashMap<String, HashSet<String>>();
	String[] inputGenes = null;
	int pathLength = 1;
	int maxSize = 100;
	static int port = 5001;
	
	public static void main(String[] args) {
		FastGenes2Networks fg2n = new FastGenes2Networks();
		
		args = new String[1];
		args[0] = "TP53,SOX2,POU5F1,TEAD4,ESRRB,CREM,E2F1,STAT3";
		
		String[] geneList = args[0].split(",");
		
		String[] networkStrings = new String[]{"BIND","HPRD","BIOGRID","HUMAP"};
		
		fg2n.initialize(networkStrings);
		
		System.out.println("G2N server running at port: "+port);
		fg2n.startMultiServer();
		
		
		//HashSet<String> geneSet = fg2n.run(networkStrings, geneList, 1);
		//System.out.println(geneSet.size());
	}
	
	public String runG2N(String[] _networkString, HashMap<String, HashSet<String>> _genes, int _pathLength) {
		
		HashMap<String, HashSet<String>> combinedNetwork = new HashMap<String, HashSet<String>>();
		
		for(int i=0; i<_networkString.length; i++){
			HashMap<String, HashSet<String>> nn = networks.get(_networkString[i]);
			for(String key : nn.keySet()) {
				if(combinedNetwork.containsKey(key)) {
					combinedNetwork.get(key).addAll(nn.get(key));
				}
				else {
					HashSet<String> temp = new HashSet<String>(nn.get(key));
					combinedNetwork.put(key, temp);
				}
			}
		}
		
		String output = "";
		
		for(String key : _genes.keySet()) {
			
			HashSet<String> inGenes = new HashSet<String>(_genes.get(key));
			inGenes.retainAll(combinedNetwork.keySet());
			
			HashSet<String> geneNeighborhood = new HashSet<String>(inGenes);
			for(int j=0; j<_pathLength; j++){
				String[] genes = geneNeighborhood.toArray(new String[0]);
				for(String gene : genes){
					geneNeighborhood.addAll(combinedNetwork.get(gene));
				}
			}
			
			String[] genearr = geneNeighborhood.toArray(new String[0]);
			Arrays.sort(genearr);
			
			String templine = key;
			for(int i=0; i<genearr.length; i++) {
				templine += ","+genearr[i];
			}
			output += templine+"\n";
		}
		//System.out.println("g2n: "+output);
		return output;
	}
	
	public void  initialize(String[] _networkStrings){
		
		networkList.put("BIND", "data/ppi/BIND.sig");
		networkList.put("BIOCARTA", "data/ppi/Biocarta.sig");
		networkList.put("BIOGRID", "data/ppi/BioGRID.sig");	
		networkList.put("DIP", "data/ppi/DIP.sig");
		networkList.put("FIGEYS", "data/ppi/figeys.sig");
		networkList.put("HPRD", "data/ppi/HPRD.sig");
		networkList.put("INNATEDB", "data/ppi/InnateDB.sig");	
		networkList.put("INTACT", "data/ppi/IntAct.sig");
		networkList.put("KEA", "data/ppi/KEA.sig");
		networkList.put("KEGG", "data/ppi/KEGG.sig");
		networkList.put("MINT", "data/ppi/MINT.sig");
		networkList.put("MIPS", "data/ppi/MIPS.sig");
		networkList.put("MURPHY", "data/ppi/MURPHY.sig");
		networkList.put("PDZBASE", "data/ppi/PDZBASE.sig");
		networkList.put("PPID", "data/ppi/PPID.sig");
		networkList.put("PREDICTEDPPI", "data/ppi/PREDICTEDPPI.sig");
		networkList.put("SNAVI", "data/ppi/SNAVI.sig");
		networkList.put("STELZL", "data/ppi/STELZL.sig");
		networkList.put("VIDAL", "data/ppi/VIDAL.sig");
		networkList.put("HUMAP", "data/ppi/huMAP.sig");
		
		networks = new HashMap<String, HashMap<String, HashSet<String>>>();
		
		long time = System.currentTimeMillis();
		String[] keys = networkList.keySet().toArray(new String[0]);
		for(int i=0; i<keys.length; i++){
			HashMap<String, HashSet<String>> net = readNetwork(networkList.get(keys[i]));
			networks.put(keys[i], net);
		}
		//System.out.println(System.currentTimeMillis() - time);
	}
	
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
	
	public HashMap<String, HashSet<String>> readNetwork(String _file){
		HashMap<String, HashSet<String>> network = new HashMap<String, HashSet<String>>();

		try{
			BufferedReader br = new BufferedReader(new FileReader(new File(_file)));
			String line = "";
			
			while((line = br.readLine()) != null){
				String[] sp = line.trim().toUpperCase().split(" ");
				if(sp.length >= 5) {
					String g1 = sp[0];
					String g2 = sp[5];
					
					if(network.containsKey(g1)){
						network.get(g1).add(g2);
					}
					else{
						HashSet<String> temp = new HashSet<String>();
						temp.add(g2);
						network.put(g1, temp);
					}
					
					if(network.containsKey(g2)){
						network.get(g2).add(g1);
					}
					else{
						HashSet<String> temp = new HashSet<String>();
						temp.add(g1);
						network.put(g2, temp);
					}
				}
			}
			br.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
		
		return network;
	}
	
	public class EchoThread extends Thread {
	    protected Socket socket;
	    private FastGenes2Networks g2nInstance;
	    private HashMap<String, HashSet<String>> geneLists = new HashMap<String, HashSet<String>>();
	    public EchoThread(Socket clientSocket, FastGenes2Networks _chea) {
	        this.socket = clientSocket;
	        g2nInstance = _chea;
	    }
	    
	    public void run() {
	    		
	    		int pathLength = 1;
	    		String[] networkStrings = new String[0];
	    		
	        try {
	            PrintWriter outs = new PrintWriter(socket.getOutputStream(), true);
				BufferedReader in = new BufferedReader(new InputStreamReader(socket.getInputStream()));
				
				String inputLine, outputLine;
				
				while ((inputLine = in.readLine()) != null) {
					outputLine = inputLine;
					
					if (outputLine.equals("kill")){
						break;
					}
					else if(outputLine.startsWith("run")){
						String[] sp = outputLine.trim().split(";");
						networkStrings = sp[1].split(",");
						pathLength = Integer.parseInt(sp[2]);
					}
					else if(outputLine.startsWith("messageComplete")) {
						String out = runG2N(networkStrings, geneLists, pathLength);
						outs.println(out+"messageComplete");
					}
					else if(outputLine.length() > 4){
						String[] sp = outputLine.trim().split(";");
						String[] sp2 = sp[1].split(",");
						HashSet<String> genes = new HashSet<String>();
						for(int i=1; i<sp2.length; i++) {
							genes.add(sp2[i]);
						}
						geneLists.put(sp2[0], genes);
					}
				}
	        }
	        catch(Exception e){
	        		e.printStackTrace();
	        }
	    }
	}
}
