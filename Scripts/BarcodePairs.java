import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Alejandro
 */
public class BarcodePairs {

    ArrayList<String> samples = new ArrayList<>();
    //2d hash with sample, and hash map with sample and value, so reads assigned correctly
    //has barcodeHash[sample_a][sample_a] = #
    HashMap<String, HashMap<String, Integer>> barcodeHash = new HashMap<>();
    //Hash with k=bc, value=sample; used to assign the bc to samples
    HashMap<String, String> mapHash = new HashMap<>();

    public static void main(String args[]) {
        String unassigned = null;
        String log1 = null;
        String log2 = null;
        String bc_file = null;
        int bc_length = 12;
        String out_File = "std";
        if (args.length == 0 || args[0].equals("-h") || args[0].equals("--help")) {
            System.out.println(help());
            System.exit(0);
        }
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-u")) {
                try {
                    i++;
                    unassigned = args[i];
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -u  \nNeed help? use BarcodePairs -h | --help");
                    System.exit(1);

                }
            } else if (args[i].equals("-b")) {
                try {
                    i++;
                    bc_file = args[i];
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -b  \nNeed help? use BarcodePairs -h | --help");
                    System.exit(1);

                }
            } else if (args[i].equals("-l1")) {
                try {
                    i++;
                    log1 = args[i];
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -l1  \nNeed help? use BarcodePairs -h | --help");
                    System.exit(1);

                }
            } else if (args[i].equals("-l2")) {
                try {
                    i++;
                    log2 = args[i];
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -l2  \nNeed help? use BarcodePairs -h | --help");
                    System.exit(1);

                }
            } else if (args[i].equals("-l")) {
                try {
                    i++;
                    bc_length = Integer.parseInt(args[i]);

                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -l \nNeed help? use BarcodePairs -h | --help");
                    System.exit(1);
                } catch (NumberFormatException nfe) {
                    System.err.println("Numeric argument expected for -c (barcode column) \nNeed help? use BarcodeCorrector -h | --help");
                    System.exit(1);

                }
            } else if (args[i].equals("-o") || args[i].equals("--output")) {
                try {
                    i++;
                    out_File = args[i];
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -o (output file) option \nNeed help? use BarcodeCorrector -h | --help");
                    System.exit(1);

                }
        }
      }
        if(bc_file != null && log1 != null && unassigned != null){
        	BarcodePairs ih = new BarcodePairs();
        	ih.processMappingFile(bc_file);
        	ih.processLogFile(log1);
        	if (log2 != null) {
            		ih.processLogFile(log2);
        	}
        	ih.process_barcodes(unassigned, bc_length);
        	ih.flush_results(out_File);
	}else{
		System.out.println("Please provide required input. See `BarcodePairs.java --help`");
	}

    }


public void flush_results(String outfile) {
        boolean to_std = true;
        try {
            FileWriter writer = null;
            if (!outfile.equals("std")) {
                to_std = false;
            } else {
                writer = new FileWriter(outfile);
            }
          
            //write the header
            for (String sample : samples) {
                if (to_std) {
                    System.out.print("\t" + sample);
                } else {
                    writer.write("\t" + sample);
                }
            }
            if (to_std) {
                System.out.println();
            } else {
                writer.write("\n");
            }

            for (int i = 0; i < samples.size(); i++) {
                if (to_std) {
                    System.out.print(samples.get(i));
                } else {
                    writer.write(samples.get(i));
                }
                String fws = samples.get(i);
                //for (int j = i; j < samples.size(); j++) {
                for (int j = 0; j < samples.size(); j++) {
                    fws = samples.get(i);
                    //arreglar esta parte!!!
                    /*   for (int k = 0; k < j; k++) {
                        String rvs = samples.get(k);
                        if (fws.compareTo(rvs) > 0) {
                            String tmp = fws;
                            fws = rvs;
                            rvs = tmp;
                        }
                        if (barcodeHash.containsKey(fws)) {
                            if (barcodeHash.get(fws).containsKey(rvs)) {
                                writer.write("\t" + barcodeHash.get(fws).get(rvs));
                            } else {
                                writer.write("\t0");
                            }
                        } else {
                            writer.write("\t0");
                        }
                    }*/
                    String rvs = samples.get(j);
                    if (fws.compareTo(rvs) > 0) {
                        String tmp = fws;
                        fws = rvs;
                        rvs = tmp;
                    }
                    if (barcodeHash.containsKey(fws)) {
                        if (barcodeHash.get(fws).containsKey(rvs)) {
                            if (to_std) {
                                System.out.print("\t" + barcodeHash.get(fws).get(rvs));
                            } else {
                                writer.write("\t" + barcodeHash.get(fws).get(rvs));
                            }
                        } else {
                            if (to_std) {
                                System.out.print("\t0");
                            } else {
                                writer.write("\t0");
                            }
                        }
                    } else {
                        if (to_std) {
                            System.out.print("\t0");
                        } else {
                            writer.write("\t0");
                        }
                    }
                }
                if (to_std) {
                    System.out.println();
                } else {
                    writer.write("\n");
                }
                //Process No_BC 
                //Not needed anymore the tag No_BC is on the samples AL.
                /*if (samples.compareTo(rvs) > 0) {
                        String tmp = fws;
                        fws = rvs;
                        rvs = tmp;
                    }
                    if (barcodeHash.containsKey(fws)) {
                        if (barcodeHash.get(fws).containsKey(rvs)) {
                            writer.write("\t" + barcodeHash.get(fws).get(rvs));
                        } else {
                            writer.write("\t0");
                        }
                    } else {
                        writer.write("\t0");
                    }*/
            }
            writer.close();
        } catch (IOException ex) {
            Logger.getLogger(BarcodePairs.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     * Process the unassigned barcodes. The distinct barcodes are the result of
     * this command: less split_test/unassigned.fna | grep "^>" | sed
     * 's/.*new_bc=// ; s/ bc_diffs.*$//' | sort | uniq -c | sort -k1gr >
     * unassigned_barcodes.txt
     *
     * @param file
     */
    public void process_barcodes(String file, int bc_length) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(file));
            String line;
            while ((line = br.readLine()) != null) {
                String tokens[] = line.trim().split(" ");
                String fwbc = tokens[1].substring(0, bc_length);
                String rvbc = tokens[1].substring(bc_length);
                String fw_k;
                String rv_k;
                if (mapHash.containsKey(fwbc)) {
                    fw_k = mapHash.get(fwbc);
                } else {
                    fw_k = "No_BC";
                }
                if (mapHash.containsKey(rvbc)) {
                    rv_k = mapHash.get(rvbc);
                } else {
                    rv_k = "No_BC";
                }
                //always populate the hash in lexicographical order, so we avoid
                //different key combinatios. e.g h[a][b] and h[b][a]
                if (fw_k.compareTo(rv_k) > 0) {
                    String tmp = fw_k;
                    fw_k = rv_k;
                    rv_k = tmp;
                }
                if (barcodeHash.containsKey(fw_k)) {
                    if (barcodeHash.get(fw_k).containsKey(rv_k)) {
                        barcodeHash.get(fw_k).put(rv_k, barcodeHash.get(fw_k).get(rv_k) + Integer.parseInt(tokens[0]));
                    } else {
                        barcodeHash.get(fw_k).put(rv_k, Integer.parseInt(tokens[0]));
                    }
                    //map.put(key, map.get(key) + 1);
                    //barcodeHash.put(tokens[0],)
                } else {
                    HashMap<String, Integer> barcodeHash2 = new HashMap<>();
                    barcodeHash2.put(rv_k, Integer.parseInt(tokens[0]));
                    barcodeHash.put(fw_k, barcodeHash2);
                }
            }
            br.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(BarcodePairs.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(BarcodePairs.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     * Process split_libraries log message
     *
     * @param file
     */
    public void processLogFile(String file) {
        // Total number seqs written       11247527
        //Median sequence length: 365.00
        BufferedReader bf;
        try {
            bf = new BufferedReader(new FileReader(file));
            String line = "";
            //this will be a hash with bc, sample, so we will also put the rev_com of the barcodes so we only
            //iterate over one hashmap
            // HashMap<String, String> barcodeHash = new HashMap<>();
            int l = 0;
            while ((line = bf.readLine()) != null) {
                l++;
                if (l > 15) {

                    if (!line.startsWith("Unassigned") && line.trim().length() > 0 && !line.startsWith("Total number seqs")) {
                        String tokens[] = line.split("\t");
                        //create 2D hashmap
                        //we run this method twice, one for each log file, so we need to add the values
                        //during the second round
                        if (tokens.length > 1) {
                            if (barcodeHash.containsKey(tokens[0])) {
                                //map.put(key, map.get(key) + 1);
                                barcodeHash.get(tokens[0]).put(tokens[0], barcodeHash.get(tokens[0]).get(tokens[0]) + Integer.parseInt(tokens[1]));
                                //barcodeHash.put(tokens[0],)
                            } else {
                                HashMap<String, Integer> barcodeHash2 = new HashMap<>();
                                barcodeHash2.put(tokens[0], Integer.parseInt(tokens[1]));
                                barcodeHash.put(tokens[0], barcodeHash2);
                            }
                        }
                    }
                }
            }
            bf.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(BarcodePairs.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(BarcodePairs.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public void processMappingFile(String file) {
        try {
            BufferedReader bf = new BufferedReader(new FileReader(file));
            String line = "";
            //this will be a hash with bc, sample, so we will also put the rev_com of the barcodes so we only
            //iterate over one hashmap
            while ((line = bf.readLine()) != null) {
                if (!line.startsWith("#")) {
                    String tokens[] = line.split("\t");
                    String bc1 = tokens[1].substring(0, 12);
                    String bc2 = tokens[1].substring(12);
                    String bc1_rc = reverseComp(bc1);
                    String bc2_rc = reverseComp(bc2);
                    mapHash.put(bc1, tokens[0]);
                    mapHash.put(bc2, tokens[0]);
                    mapHash.put(bc1_rc, tokens[0]);
                    mapHash.put(bc2_rc, tokens[0]);
                    samples.add(tokens[0]);

                }
            }
            bf.close();
            Collections.sort(samples);
            //When there is no combination or the bc is unassigned we use this TAG
            samples.add("No_BC");

        } catch (FileNotFoundException ex) {
            Logger.getLogger(BarcodePairs.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(BarcodePairs.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public String reverseComp(String cadena) {
        String invertida = "";
        if (cadena != null) {
            for (int i = cadena.length() - 1; i >= 0; i--) {
                char base = cadena.charAt(i);
                if (base == 'A') {
                    base = 'T';
                } else if (base == 'T') {
                    base = 'A';
                } else if (base == 'C') {
                    base = 'G';
                } else if (base == 'G') {
                    base = 'C';
                } else if (base == 'N') {
                    base = 'N';
                } else {
                    System.err.println("Caracter No Esperado. SUtils.reversoComplementarioNuc: " + base + "\nSecuencia: " + cadena);
                }
                invertida += base;
            }
        }
        return invertida;
    }

    private static String help() {
        String help = "\n#######################################################\n"
                + "###                 Barcode pairs                   ###\n"
                + "###                    v 1.0                        ###\n"
                + "###                             @company       NIOZ ###\n"
                + "###                             @author   A. Abdala ###\n"
                + "#######################################################\n\n"
                + "Program implemented for Cascabel pipeline. It needs the split_libraries_fastq.py log file\n"
                + "and the unassignned fasta file. Then, it compute all the barcode pairs ocurrences between samples.\n\n"
                + "usage java BarcodePairs  -b Barcode_file -u unassigned_reads -l1 split_log1 [-l2 split_log2]\n\n"
                + "Mandatory arguments:\n"
                + "  -b\tTab separated file with sample name in the first column and barcode information in the second column.\n"
                + "  -u\tFile with the summary of the unassigned barcodes\n"
                + "      \tThis file is the result of running this command: less unassigned.fna | grep \"^>\" | sed\n"
                + "     's/.*new_bc=// ; s/ bc_diffs.*$//' | sort | uniq -c | sort -k1gr > unassigned_barcodes.txt\n"
                + "  -l1\tLog file from split_libraries_fastq.py.\n"
                + "\nOptional arguments:\n"
                + "  -l2\tLog file from split_libraries_fastq.py (second pass).\n"
                + "  -l\tBarcode length. (For paired-end data FW_BC+RV_BC). [Def. 12]\n"
                + "  -o \tOutputfile. [Def. std_out]\n";

        return help;
    }
}

