/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package clustermapper;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author aabdala
 */
public class FileProcessor {

    private ArrayList<String> samples = new ArrayList<>();
    private boolean verbose = true;
    private boolean strict = false;//abort on error - keys not found 

    public ArrayList<String> getSamples() {
        return samples;
    }

    public void setSamples(ArrayList<String> samples) {
        this.samples = samples;
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    public boolean isStrict() {
        return strict;
    }

    public void setStrict(boolean strict) {
        this.strict = strict;
    }

    /**
     * Method to process uc files. It creates a hash map with all the clusters
     * on the uc files passed
     *
     * @param uc_files. Comma separated lists of uc files to process
     * @param sequenceid_delim. Sequence delim, in general this can be a
     * semicolon (;) A uc file can have sequence identifieres like the
     * following: SEQID;size=### common while working with vsearch / usearch so
     * this parameter tell the method what character indicates the sequence ID
     * @param sample_delim. WHile parsing the file, in order to compute the
     * sequence frequencies, we also need to now the provenance of the
     * sequences, so in general we spect the SEQUENCEID, to be composed by the
     * SAMPLEID+some sequence number, so in general the SEQID looks like:
     * SAMPID_### so for this case the sample_delim should be "_"
     * @param process_all_uc the file processing stops as soon as the cluster
     * session is reached (default behaviour). We add this flag in case we need
     * to deal with some strange formating where the clustering section could
     * not be at the very end of the file.
     * @return a Hash table with all the observations. Examples of sequence
     * format: The sequences could be in this two formats:
     * ERR1755873_1144;size=16261312 or just: ERR1879724_311 The sequence id
     * should be: ERR1755873_1144 The sample id should be: ERR1755873 And for
     * the moment we really dont care about the abundance;
     */
    public HashMap<String, Observations> processUCMap(String uc_files, String sequenceid_delim, String sample_delim, boolean process_all_uc) {
        HashMap<String, Observations> hmap = new HashMap<>();
        if (verbose) {
            System.out.println("*************************");
            System.out.println("PROCESS UC FILE(S)");
            System.out.println("Checking files...");
        }
        long start = System.currentTimeMillis();
        try {
            for (String ucInput : uc_files.split(",")) {
                if (verbose) {
                    System.out.println("\tuc file: " + ucInput);
                    System.out.println("\tCreating map from uc file...");
                }
                BufferedReader reader = new BufferedReader(new FileReader(ucInput));
                int lnum = 0;
                String linea;
                while ((linea = reader.readLine()) != null) {
                    lnum++;
                    if (linea.startsWith("S")) {
                        String tokens[] = linea.split("\t");
                        String id = tokens[8].contains(sequenceid_delim) ? tokens[8].substring(0, tokens[8].indexOf(sequenceid_delim)) : tokens[8];
                        String sample_id = tokens[8].contains(sample_delim) ? tokens[8].substring(0, tokens[8].indexOf(sample_delim)) : tokens[8];
                        if (!samples.contains(sample_id)) {
                            samples.add(sample_id);
                        }
                        if (!hmap.containsKey(id)) {
                            Observations ob = new Observations();
                            ob.addSampleObservation(sample_id, 1);
                            ob.addSequence(id);
                            ob.setSeedSequence(id);
                            hmap.put(id, ob);

                        } else {
                            // hmap.get(id).addSampleObservation(sample_id, 1);
                            System.err.println("Duplicate Key Error: " + id + "\nAtl line:" + linea);
                            if (strict) {
                                System.err.println("Please verify that there are no repeated key clusters within your uc file.\n"
                                        + "Strict mode on. Aborting...");
                            }
                        }
                    } else if (linea.startsWith("H")) {
                        String tokens[] = linea.split("\t");
                        String id = tokens[9].contains(sequenceid_delim) ? tokens[9].substring(0, tokens[9].indexOf(sequenceid_delim)) : tokens[9];
                        String sample_id = tokens[8].contains(sample_delim) ? tokens[8].substring(0, tokens[8].indexOf(sample_delim)) : tokens[8];
                        String sequence_id = tokens[8].contains(sequenceid_delim) ? tokens[8].substring(0, tokens[8].indexOf(sequenceid_delim)) : tokens[8];
                        if (!samples.contains(sample_id)) {
                            samples.add(sample_id);
                        }
                        if (hmap.containsKey(id)) {
                            hmap.get(id).addSampleObservation(sample_id, 1);
                            hmap.get(id).addSequence(sequence_id);
                        } else {
                            /**
                             * This exception can be thrown when there is a hit
                             * H with a single sequence S but the S key is not
                             * found at the hashmap
                             */
                            System.err.println("Key Error: " + id
                                    + "\n -Reference cluster key not found on data. (H without S record on uc file)."
                                    + "\nAt line:" + linea);
                            if (strict) {
                                System.err.println("Clustering error within your uc file.\n"
                                        + "H record should be preceded by a complementary S record (seed for the cluster)\n"
                                        + "Strict mode on. Aborting...");
                            }
                        }
                    } else if (linea.startsWith("C") && !process_all_uc) {
                        if (verbose) {
                            System.out.println("Cluster section reached. Closing file...");
                        }
                        break;
                    }
                }
                long end = (System.currentTimeMillis() - start) / 1000;
                if (verbose) {
                    System.out.println("Map created on: " + end + "s.");
                    System.out.println("Lines processed: " + lnum);
                    System.out.println("Number of elements into the map: " + hmap.size());
                    System.out.println("********************");
                }
                reader.close();
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(FileProcessor.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(FileProcessor.class.getName()).log(Level.SEVERE, null, ex);
        }
        return hmap;
    }

    /**
     * This method recluster a hash map given a new uc file as reference. This
     * method is used when a second clusterization is performes. For instance if
     * first we do a dereplicartion and after that we want to pick OTUs, we need
     * to remap the abundances from the first dereplication using the previous
     * hashmap and the result of the new clusterization. This method can also be
     * useful when the dereplication is performed in more than one chunk,
     * therefore you first derep each chunk, then concat the output fasta of
     * each one and then dereplicate again. To create a hash map with the
     * remaped sequence abundance, we use this method.
     *
     * @param remap Hashmap with the first clustering results
     * @param uc_map reference for the new uc file
     * @param sequenceid_delim. Sequence delim, in general this can be a
     * semicolon (;) A uc file can have sequence identifieres like the
     * following: SEQID;size=### common while working with vsearch / usearch so
     * this parameter tell the method what character indicates the sequence ID
     * @param process_all_uc the file processing stops as soon as the cluster
     * session is reached (default behaviour). We add this flag in case we need
     * to deal with some strange formating where the clustering section could
     * not be at the very end of the file.
     * @param sample_delim. WHile parsing the file, in order to compute the
     * sequence frequencies, we also need to now the provenance of the
     * sequences, so in general we spect the SEQUENCEID, to be composed by the
     * SAMPLEID+some sequence number, so in general the SEQID looks like:
     * SAMPID_### so for this case the sample_delim should be "_"
     *
     * @return
     */
    public HashMap<String, Observations> reMapUcMap(HashMap<String, Observations> remap, String uc_map, String sequenceid_delim, String sample_delim, boolean process_all_uc) {
        HashMap<String, Observations> hmap = new HashMap<>();
        if (verbose) {
            System.out.println("*************************");
            System.out.println("REMAPPING CLUSTERS");
            System.out.println("\tUC map:" + uc_map);

        }
        long start = System.currentTimeMillis();
        try {
            BufferedReader reader = new BufferedReader(new FileReader(uc_map));
            int lnum = 0;
            String linea;
            while ((linea = reader.readLine()) != null) {
                lnum++;
                if (linea.startsWith("S")) {
                    String tokens[] = linea.split("\t");
                    String id = tokens[8].contains(sequenceid_delim) ? tokens[8].substring(0, tokens[8].indexOf(sequenceid_delim)) : tokens[8];
                    Observations obs = remap.get(id);
                    if (obs != null) {
                        obs.setSeedSequence(id);
                        hmap.put(id, obs);
                        remap.remove(id);
                    } else {
                        System.err.println("There is no entry for key: " + id
                                + "\nThis key wont be included!\n"
                                + "Line:" + lnum);
                        if (strict) {
                            System.err.println("Inconsistent data.\n"
                                    + "S record should be referenced on previous HashMap\n"
                                    + "Strict mode on. Aborting...");
                        }
                    }
                } else if (linea.startsWith("H")) {
                    String tokens[] = linea.split("\t");
                    String id = tokens[9].contains(sequenceid_delim) ? tokens[9].substring(0, tokens[9].indexOf(sequenceid_delim)) : tokens[9];
                    String old_id = tokens[8].contains(sequenceid_delim) ? tokens[8].substring(0, tokens[8].indexOf(sequenceid_delim)) : tokens[8];
                    Observations obs = remap.get(old_id);
                    if (obs != null) {
                        hmap.get(id).mergeObservations(obs);
                        remap.remove(old_id);
                    } else {
                        System.err.println("There is no entry for key: " + id
                                + "\nThis key wont be included!\n"
                                + "Line:" + lnum);
                        if (strict) {
                            System.err.println("Inconsistent data.\n"
                                    + "H record should be referenced on previous HashMap\n"
                                    + "Strict mode on. Aborting...");
                        }
                    }

                } else if (linea.startsWith("C") && !process_all_uc) {
                    if (verbose) {
                        System.out.println("Cluster section reached. Closing file...");
                    }
                    break;
                }
            }
            long end = (System.currentTimeMillis() - start) / 1000;
            if (verbose) {
                System.out.println("Re Map finished on: " + end + "s.");
                System.out.println("Number of elements into the new map: " + hmap.size());
                System.out.println("********************");
            }
            reader.close();

        } catch (FileNotFoundException ex) {
            Logger.getLogger(FileProcessor.class.getName()).log(Level.SEVERE, null, ex);
            System.err.println("There is no file: " + uc_map);
            System.err.println("Aborting...");
            System.exit(21);
        } catch (IOException ex) {
            Logger.getLogger(FileProcessor.class.getName()).log(Level.SEVERE, null, ex);
            System.err.println("Error reading file: " + uc_map);
            System.err.println("Aborting...");
            System.exit(21);
        }
        return hmap;
    }

    /**
     * This method recluster a hash map given a OTU map file as reference. This
     * method is very similar to reMapUcMap but it use used when a second
     * clusterization is performes. For instance if first we do a dereplicartion
     * and after that we want to pick OTUs, we need to remap the abundances from
     * the first dereplication using the previous hashmap and the result of the
     * new clusterization. This method can also be useful when the dereplication
     * is performed in more than one chunk, therefore you first derep each
     * chunk, then concat the output fasta of each one and then dereplicate
     * again. To create a hash map with the remaped sequence abundance, we use
     * this method. a OTU map record looks like the following:
     * OTUID\tSEQ1\tSEQ2...SEQN \nOTUID2\tSEQ1tSEQ2\SEQN \n...SEQM
     *
     * @param map The path to the OTU txt map
     * @param remap Hashmap with previous clustering results
     * @param sequenceid_delim. Sequence delim, in general this can be a
     * semicolon (;) A OTU file can have sequence identifieres like the
     * following: SEQID;size=### common while working with vsearch / usearch so
     * this parameter tell the method what character indicates the sequence ID
     * @param record_delim the delimiter between records, default tab
     * @param has_label if the first element is a sequence ID by itself (false)
     * @return
     */
    public HashMap<String, Observations> reMapOtuMap(HashMap<String, Observations> remap, String map, String sequenceid_delim, String record_delim, boolean has_label) {
        HashMap<String, Observations> hmap = new HashMap<>();
        if (verbose) {
            System.out.println("*************************");
            System.out.println("REMAPPING CLUSTERS");
            System.out.println("\tOTU map:" + map);
        }

        long start = System.currentTimeMillis();
        int collapsed = 0;
        int initial_size = remap.size();
        try {
            BufferedReader reader = new BufferedReader(new FileReader(map));
            int lnum = 0;
            String linea;
            while ((linea = reader.readLine()) != null) {
                lnum++;
                String tokens[] = linea.split(record_delim);
                //String id = tokens[0];                
                int idx = 0;
                String id = tokens[idx].contains(sequenceid_delim) ? tokens[idx].substring(0, tokens[idx].indexOf(sequenceid_delim)) : tokens[idx];
                //el primer elemento solo es un nombre / label, como OTU 1, de lo contrario el id es el id de una secuencia en el hash
                if (has_label) {
                    idx = 1;
                }
                // String id = tokens[idx].contains(sequenceid_delim) ? tokens[idx].substring(0, tokens[idx].indexOf(sequenceid_delim)) : tokens[idx];
                boolean isFirst = true;
                for (int i = idx; i < tokens.length; i++) {
                    String old_id = tokens[i].contains(sequenceid_delim) ? tokens[i].substring(0, tokens[i].indexOf(sequenceid_delim)) : tokens[i];
                    Observations obs = remap.get(old_id);
                    remap.remove(old_id);
                    if (obs != null) {
                        if (isFirst) {
                            isFirst = false;
                            hmap.put(id, obs);
                        } else {
                            hmap.get(id).mergeObservations(obs);
                            collapsed++;
                        }
                    } else {
                        System.err.println("There is no entry for key: " + id
                                + "\nThis key wont be included!\n"
                                + "Line:" + lnum);
                        if (strict) {
                            System.err.println("Inconsistent data.\n"
                                    + "Sequence ID should be referenced on previous HashMap\n"
                                    + "Strict mode on. Aborting...");
                        }
                    }
                }
            }
            long end = (System.currentTimeMillis() - start) / 1000;
            if (verbose) {
                System.out.println("Collapsed observations: " + collapsed + " from " + initial_size);
                System.out.println("Observations collapsed on: " + end + "s.");
                System.out.println("Number of elements into the new map: " + hmap.size());
                System.out.println("********************");
            }
            reader.close();

        } catch (FileNotFoundException ex) {
            Logger.getLogger(FileProcessor.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(FileProcessor.class.getName()).log(Level.SEVERE, null, ex);
        }
        return hmap;
    }

    /**
     * Method to print the hashtable. This method creates a OTU map table ready
     * to be process by either QIIME scripts or biom utilities.
     *
     * @param outfile the target outputfile where to write the map.
     * @param map The hash map to be flushed into the OTU Map
     * @param otu_str This kind of tables contains on the first column the
     * cluster ID. The default name of this column is given by this string
     * (OTUID by default). The other names cames from the name of the samples
     * @param relabel true if use a different name for the OTU ids, by default
     * is the name of the sequence seed id
     * @param label if relabel true, the name like
     * @param labelidx incremental number for new label id The format for this
     * file is clusterid\tseed_id\tseqs_ids\t as required by QIIME
     */
    public void printOTUMap(String outfile, HashMap<String, Observations> map, String otu_str, boolean relabel, String label, int labelidx) {
        String nl = System.lineSeparator();

        try {
            FileWriter writer = null;
            if (outfile.trim().length() > 1) {
                writer = new FileWriter(outfile);
            }
            long start = System.currentTimeMillis();
            if (verbose) {
                System.out.println("*************************");
                System.out.println("PRINTING OTU MAP");
                System.out.println("\ttarget outputfile: " + outfile);
            }

            for (Map.Entry<String, Observations> entry : map.entrySet()) {
                Observations obs = entry.getValue();
                if (relabel) {
                    writer.write(label + labelidx);
                    ++labelidx;
                } else {
                    writer.write(entry.getKey());
                }
                for (String sequence : obs.getSequences()) {
                    writer.write("\t" + sequence);
                }
                writer.write(nl);

            }
            if (writer != null) {
                writer.close();
            }
            long end = (System.currentTimeMillis() - start) / 1000;
            System.out.println("finished on: " + end + "s.");
            System.out.println("********************");
        } catch (IOException ex) {
            Logger.getLogger(FileProcessor.class.getName()).log(Level.SEVERE, null, ex);
            System.err.println("Error writing file: " + outfile);
            System.err.println("Aborting...");
            System.exit(21);
        }

    }

    /**
     * Method to print the hashtable. This method creates a OTU map table ready
     * to be process by either QIIME scripts or biom utilities.
     *
     * @param outfile the target outputfile where to write the map.
     * @param map The hash map to be flushed into the OTU Map
     * @param otu_str This kind of tables contains on the first column the
     * cluster ID. The default name of this column is given by this string
     * (OTUID by default). The other names cames from the name of the samples
     * @param relabel true if use a different name for the OTU ids, by default
     * is the name of the sequence seed id
     * @param label if relabel true, the name like
     * @param labelidx incremental number for new label id The output from this
     * file can be used directly with biom utilities to create a biom table
     */
    public void printOTUTable(String outfile, HashMap<String, Observations> map, String otu_str, boolean relabel, String label, int labelidx) {
        String nl = System.lineSeparator();

        try {
            FileWriter writer = null;
            if (outfile.trim().length() > 1) {
                writer = new FileWriter(outfile);
            }
            long start = System.currentTimeMillis();
            if (verbose) {
                System.out.println("*************************");
                System.out.println("PRINTING OTU Table");
                System.out.println("\ttarget outputfile: " + outfile);
            }
            writer.write("#" + otu_str);
            for (String sample : samples) {
                if (writer != null) {
                    writer.write("\t" + sample);
                }
            }
            if (writer != null) {
                writer.write(nl);
            }
            for (Map.Entry<String, Observations> entry : map.entrySet()) {
                Observations obs = entry.getValue();

                if (writer != null) {
                    if (relabel) {
                        writer.write(label + labelidx);
                        ++labelidx;
                    } else {
                        writer.write(entry.getKey());
                    }
                } else {
                    System.out.print(entry.getKey());
                }
                for (String sample : samples) {
                    if (writer != null) {
                        writer.write("\t" + obs.getAbundanceBySample(sample));
                    } else {
                        System.out.println("\t" + obs.getAbundanceBySample(sample));
                    }

                }
                if (writer != null) {
                    writer.write(nl);
                } else {
                    System.out.println();
                }
            }
            if (writer != null) {
                writer.close();
            }
            long end = (System.currentTimeMillis() - start) / 1000;
            System.out.println("finished on: " + end + "s.");
            System.out.println("********************");
        } catch (IOException ex) {
            Logger.getLogger(FileProcessor.class.getName()).log(Level.SEVERE, null, ex);
            System.err.println("Error writing file: " + outfile);
            System.err.println("Aborting...");
            System.exit(21);
        }

    }
}
