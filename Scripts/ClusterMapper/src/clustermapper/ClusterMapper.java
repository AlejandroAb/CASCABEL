/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package clustermapper;

import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author aabdala
 */
public class ClusterMapper {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        String mode = "";
        if (args.length < 1) {
            printHelp();
            System.exit(0);
        }
        long start = System.currentTimeMillis();
        String sample_delim = "_";
        String abundance_delim = "=";
        String sequenceid_delim = ";";
        String ucFile = "";
        String ucFile2 = "";
        String outFile = "table.out";
        String otuFile = "";
        String otuHeader = "OTUID";
        String txtFieldSeparator = "\t";
        boolean verbose = true;
        boolean relabel = false;
        boolean print_otu_tbl = false;
        String label = "otu";
        int labelidx = 0;
        /**
         * for the mapping if the uc file is well formated, we only need to
         * process it until we reach the cluster section. In cases where this
         * file is not ordered, maybe this flag could solve later issues... set
         * the flag to true in order to continue processing the file until it
         * ends...
         */
        boolean process_all_uc = false;
        boolean strict = false;
        boolean has_label = true;
        /*  String uc_chunks = args[0];//"/home/NIOZ.NL/aabdala/projects/ucmap/test_files/results_100_vsearch.uc,/home/NIOZ.NL/aabdala/projects/ucmap/test_files/results_200_vsearch.uc,/home/NIOZ.NL/aabdala/projects/ucmap/test_files/results_300_vsearch.uc";//args[0]
        String uc_result = args[1];//"/home/NIOZ.NL/aabdala/projects/ucmap/test_files/results_all.uc";//args[1]
        String fasta_relabel = args[2];//"/home/NIOZ.NL/aabdala/projects/ucmap/test_files/results_all_relabel.fasta";
        String swarms = args[3];//"/home/NIOZ.NL/aabdala/projects/ucmap/test_files/swarm.out";
         */
        ArrayList<String> modes = new ArrayList<>();
        modes.add("uc2uc");
        modes.add("uc2otu");
        for (int i = 0; i < args.length; i++) {
            if (i == 0 && (!args[i].equals("-h") && !args[i].equals("--help"))) {
                mode = args[i].toLowerCase();
                if (!modes.contains(mode)) {
                    System.out.println("Invalid mode!");
                    System.out.println("Valid values are: uc2uc, uc2otu");
                    System.out.println("Need help? use java ClusterMapper -h | --help");
                    printHelp();
                    System.exit(1);
                }
            } else if (args[i].equals("-h") || args[i].equals("--help")) {
                printHelp();
                System.exit(0);
            } else if (args[i].equals("-uc")) {
                try {
                    ucFile = args[i + 1];
                    i++;
                } catch (ArrayIndexOutOfBoundsException aiobe) {
                    System.out.println("Option -uc. Argument expected!");
                    printHelp();
                    System.exit(1);
                }
            } else if (args[i].equals("-uc2")) {
                try {
                    ucFile2 = args[i + 1];
                    i++;
                } catch (ArrayIndexOutOfBoundsException aiobe) {
                    System.out.println("Option -uc2. Argument expected!");
                    printHelp();
                    System.exit(1);
                }
            } else if (args[i].equals("-otu")) {
                try {
                    otuFile = args[i + 1];
                    i++;
                } catch (ArrayIndexOutOfBoundsException aiobe) {
                    System.out.println("Option -otu. Argument expected!");
                    printHelp();
                    System.exit(1);
                }
            } else if (args[i].equals("-o")) {
                try {
                    outFile = args[i + 1];
                    i++;
                } catch (ArrayIndexOutOfBoundsException aiobe) {
                    System.out.println("Option -o. Argument expected!");
                    printHelp();
                    System.exit(1);
                }
            } else if (args[i].equals("-otuh")) {
                try {
                    otuHeader = args[i + 1];
                    i++;
                } catch (ArrayIndexOutOfBoundsException aiobe) {
                    System.out.println("Option -otuh. Argument expected!");
                    printHelp();
                    System.exit(1);
                }
            } else if (args[i].equals("-l")) {
                try {
                    label = args[i + 1];
                    i++;
                } catch (ArrayIndexOutOfBoundsException aiobe) {
                    System.out.println("Option -l. Argument expected!");
                    printHelp();
                    System.exit(1);
                }
            } else if (args[i].equals("-lidx")) {
                try {
                    labelidx = Integer.parseInt(args[i + 1]);
                    i++;
                } catch (ArrayIndexOutOfBoundsException aiobe) {
                    System.out.println("Option -lidx. Argument expected!");
                    printHelp();
                    System.exit(1);
                } catch (NumberFormatException nfe) {
                    System.out.println("Option -lidx. Numerical argument expected!\nFound:" + args[i + 1]);
                    printHelp();
                    System.exit(1);

                }
            } else if (args[i].equals("--relabel")) {
                relabel = true;
            } else if (args[i].equals("--no-otuid")) {
                has_label = false;
            } else if (args[i].equals("--silent")) {
                verbose = false;
            } else if (args[i].equals("--full-uc")) {
                process_all_uc = true;
            } else if (args[i].equals("--strict")) {
                strict = true;
            } else if (args[i].equals("--otu-tbl")) {
                print_otu_tbl = true;
            } else {
                System.err.println("Unknown option: " + args[i] + "\nNeed help? use java ClusterMapper -h | --help");
                System.exit(10);
            }
        }
        FileProcessor fp = new FileProcessor();
        fp.setStrict(strict);
        fp.setVerbose(verbose);

        if (mode.equals("uc2uc")) {
            if (ucFile.length() > 1 && ucFile2.length() > 1) {
                HashMap<String, Observations> firstUC = fp.processUCMap(ucFile, sequenceid_delim, sample_delim, process_all_uc);
                //processUCMap(uc_chunks, sequenceid_delim, sample_delim);
                HashMap<String, Observations> finalMap = fp.reMapUcMap(firstUC, ucFile2, sequenceid_delim, sample_delim, process_all_uc);
                //reMapMap(firstUC, uc_result, sequenceid_delim, sample_delim);
                fp.printOTUMap(outFile, finalMap, otuHeader, relabel, label, labelidx);
                if (print_otu_tbl) {
                    fp.printOTUTable(outFile + ".biom.txt", finalMap, otuHeader, relabel, label, labelidx);
                }
            } else {
                System.err.println("For mode uc2uc please supply both uf files (-uc and -uc2)\nNeed help? use java ClusterMapper -h | --help");
                printHelp();
            }
        } else if (mode.equals("uc2otu")) {
            if (ucFile.length() > 1 && otuFile.length() > 1) {
                HashMap<String, Observations> firstUC = fp.processUCMap(ucFile, sequenceid_delim, sample_delim, process_all_uc);
                //processUCMap(uc_chunks, sequenceid_delim, sample_delim);
                HashMap<String, Observations> finalMap = fp.reMapOtuMap(firstUC, otuFile, sequenceid_delim, txtFieldSeparator, has_label);
                //reMapMap(firstUC, uc_result, sequenceid_delim, sample_delim);
                fp.printOTUMap(outFile, finalMap, otuHeader, relabel, label, labelidx);
                if (print_otu_tbl) {
                    fp.printOTUTable(outFile + ".biom.txt", finalMap, otuHeader, relabel, label, labelidx);
                }
            } else {
                System.err.println("For mode uc2otu please supply both uf files (-uc and -otu)\nNeed help? use java ClusterMapper -h | --help");
                printHelp();
            }
        } else {
            System.err.println("Invalid mode!\nNeed help? use java ClusterMapper -h | --help");
            printHelp();
            System.exit(10);
        }

        //first we map the sequences
        //  HashMap<String, Observations> firstUC = processUCMap(uc_chunks, sequenceid_delim, sample_delim);
        // HashMap<String, Observations> finalMap = reMapMap(firstUC, uc_result, sequenceid_delim, sample_delim);
        // relabelTable(finalMap, sequenceid_delim, fasta_relabel);
        //  HashMap<String, Observations> collapsedMap = processTxtMap(finalMap, swarms, sequenceid_delim, false);
        // printTable("table.out", collapsedMap);
        //printRelabeledTable("table.out", finalMap, sequenceid_delim, fasta_relabel);
        long end = (System.currentTimeMillis() - start) / 1000;
        System.out.println("****TOTAL TIME: " + end + "s.");
        System.exit(0);
    }

    private static void printHelp() {
        System.out.println("************************************************");
        System.out.println("*              Cluster Mapping tool            *");
        System.out.println("*                                              *");
        System.out.println("**@company: NIOZ                               *");
        System.out.println("**@author:  Alejandro Abdala                  **");
        System.out.println("**@version: 2.0                               **");
        System.out.println("************************************************");
        System.out.println("Usage: java ClusterMapper <mode> [options]");
        System.out.println("Modes:");
        System.out.println("  uc2uc      Map two different uc files. i.e. First uc file from dereplication");
        System.out.println("             and the second one from the OTU picking. Usefull for methods like");
        System.out.println("             uclust, vsearch, etc.");
        System.out.println("             Mandatory arguments for this mode: -uc and -uc2");
        System.out.println("  uc2otu     Map a uc file and the second clustering should be a OTU map on txt");
        System.out.println("             format. i.e QIIME output from pick_otus.py script");
        System.out.println("             Mandatory arguments for this mode: -uc and -otu");
        System.out.println("Input options:");
        System.out.println(" -uc [file]  UC file generated by usearch, vsearch or any alike clustering method");
        System.out.println("             This file should be the result of the first clusterization / dereplication.");
        System.out.println("             This argument accepts comma separated list of uc files");
        System.out.println(" -uc2 [file] UC file generated by usearch, vsearch or any alike clustering method");
        System.out.println("             This file should be the result of the second clusterization / dereplication.");
        System.out.println(" -otu [file] OTU txt table. Frequency table with the name of the OTU ids (first column) ");
        System.out.println("             and the abundance per sample.");
        System.out.println("             This file should be the result of the second clusterization / dereplication.");
        System.out.println("Output options:");
        System.out.println(" -o [file]   Output file to write new \"re-mapped\" file. This output is like the QIIME's");
        System.out.println("             default otu map after pick_otus.py script:clusterid seed_seq_id seq_id_2");
        System.out.println("             Default table.out");
        System.out.println(" --otu-tbl   If you want a biom like OTU table (frequency table) also pass this flag.");
        System.out.println("             The output name will be the same as the default output file (-o) plus \".biom.txt\"");
        System.out.println(" -otuh [str] Header for the OTU column on the otu table file (first column). Default: OTUID");
        System.out.println(" --relabel   Default output have the ID of the seed sequence as OTU name. Use this flag");
        System.out.println("             to relabel the out table with -l label. Apply for both -o and --otu-tbl");
        System.out.println(" -l [str]    Label string for naming OTUs in the outputfile.");
        System.out.println("             Only apply when relabeling is enable (--relabel). default \"otu\"");
        System.out.println(" -lidx       If --relabel otu table, initial number to be appended to the otu label (-l).");
        System.out.println("             The name of the new label will be \"label+lidx\"");
        System.out.println(" --silent    Do not output progress / actions of the mapping tool");
        System.out.println("File parsing options:");
        System.out.println(" -samd [str] Sample delimiter. Character or string used to delimiter the sample within");
        System.out.println("             a given sequence accession. For instance NIOZ90.1_328 for this accession");
        System.out.println("             NIOZ90.1 is the sample name and 328 the sequence number, therefore the");
        System.out.println("             sample delimiter will be an underscore \"_\". Default value:\"_\" ");
        System.out.println(" -seqd [str] Sequence delimiter. Character or string used to delimiter the sequence ID");
        System.out.println("             within a given sequence accession. For instance some clustering tools");
        System.out.println("             append abundance info to the sequence accession. i. e: ");
        System.out.println("             NIOZ90.1_328;size=47 for this case we need to exclude the abundance information");
        System.out.println("             and treat only the first part (NIOZ90.1_328) as the sequence accession.");
        System.out.println("             Default value: \";\"");
        System.out.println(" -otud [str] Record delimiter used for parsing the otu map (-otu). Most of the clustering    ");
        System.out.println("             tools split this file as tsv (tab separated values) but programs like \"swarm\"");
        System.out.println("             creates a similar output but with records separated with the space character");
        System.out.println("             Default value: \"\\t\" (tab character)");
        System.out.println(" --no-otuid  While processing the otu maps (-otu) is very likly that the first column is the");
        System.out.println("             name of the new otuid, if the first column is the seed sequence, pass this flag");
        System.out.println("             for the correct processing of the otu file.");
        System.out.println(" --full-uc   Process uc files till the end. Normal behaviour is to stop processing the");
        System.out.println("             the file as soon as the cluster section has been reached.");
    }

}
