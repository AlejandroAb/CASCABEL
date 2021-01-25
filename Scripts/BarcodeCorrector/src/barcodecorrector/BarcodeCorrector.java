/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package barcodecorrector;

/**
 *
 * @author Alegandro
 */
public class BarcodeCorrector {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        /*Un comment this for single tests
        DistanceTools dt = new DistanceTools(-1, -1, 1);
          float distance = dt.compute("CGTATAAATGCGGTATGCGAGGTC", "CGTATAAATGCGTATGCGAGGTCG", false);
          int distance2 = dt.levenshteinDistance("CGTATAAATGCGGTATGCGAGGTC", "CGTATAAATGCGTATGCGAGGTCG");
          int distance3 = dt.levenshteinDistanceReplacement("CGTATAAATGCGGTATGCGAGGTC", "CGTATAAATGCGTATGCGAGGTCG");
           System.out.println("motif: " + distance);
           System.out.println("levenshteinDistance " + distance2);
           System.out.println("levenshtein Replacement " + distance3);
           */
        if (args.length == 0 || args[0].equals("-h") || args[0].equals("--help")) {
            System.out.println(help());
            System.exit(0);
        }
        // System.out.println("|_.Comparition|_.Levenshtein|_.Levenshtein Substitution|_.Motif implementation|");
        long start = System.currentTimeMillis();
        // String fastq_file = "C:\\Users\\Alegandro\\Documents\\Projects_\\NIOZ\\AmpliPipeline\\BarcodeCorrector\\barcodes.10k.with_wrongs.fastq";
        String fastq_file = null;
        String barcode_file = null;
        String output_file = null;
        int barcode_column = 2;
        boolean debug = false;
        int misMatch = 2;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-f") || args[i].equals("--fastq")) {
                try {
                    i++;
                    fastq_file = args[i];
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -f (fastq file) \nNeed help? use BarcodeCorrector -h | --help");
                    System.exit(1);

                }
            } else if (args[i].equals("-b") || args[i].equals("--barcodes")) {
                try {
                    i++;
                    barcode_file = args[i];
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -b (barcode file) option \nNeed help? use BarcodeCorrector -h | --help");
                    System.exit(1);

                }
            } else if (args[i].equals("-o") || args[i].equals("--output")) {
                try {
                    i++;
                    output_file = args[i];
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -o (output file) option \nNeed help? use BarcodeCorrector -h | --help");
                    System.exit(1);

                }
            } else if (args[i].equals("-d") || args[i].equals("--debug")) {
                debug = true;

            } else if (args[i].equals("-m") || args[i].equals("--max-mismatch")) {
                try {
                    i++;
                    misMatch = Integer.parseInt(args[i]);

                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -m (max mismatch) option \nNeed help? use BarcodeCorrector -h | --help");
                    System.exit(1);
                } catch (NumberFormatException nfe) {
                    System.err.println("Numeric argument expected for -m (max mismatch) \nNeed help? use BarcodeCorrector -h | --help");
                    System.exit(1);

                }
            } else if (args[i].equals("-c") || args[i].equals("--barcode-column")) {
                try {
                    i++;
                    barcode_column = Integer.parseInt(args[i]);

                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -c (barcode column) option \nNeed help? use BarcodeCorrector -h | --help");
                    System.exit(1);
                } catch (NumberFormatException nfe) {
                    System.err.println("Numeric argument expected for -c (barcode column) \nNeed help? use BarcodeCorrector -h | --help");
                    System.exit(1);

                }
            }
        }
        // String fastq_file = "C:\\Users\\Alegandro\\Documents\\Projects_\\NIOZ\\AmpliPipeline\\BarcodeCorrector\\barcodes.10k.fastq";
        // String barcode_file = "C:\\Users\\Alegandro\\Documents\\Projects_\\NIOZ\\AmpliPipeline\\BarcodeCorrector\\mapNIOZ140_AA.txt";
        if (fastq_file == null) {
            System.err.println("Please supply a fastq file. Option -f");
            System.exit(0);
        }
        if (barcode_file == null) {
            System.err.println("Please supply a barcode file. Option -b");
            System.exit(0);
        }
        if (output_file == null) {
            output_file = fastq_file + "_corrected";
        }
        FastqBCComparator comparator = new FastqBCComparator(fastq_file, barcode_file);
        comparator.setDebug(debug);
        comparator.setBc_column(barcode_column);
        comparator.correctBarcodes(misMatch, output_file);

        /*
         //String pattern = "3214";
        //String toCompare[] = {"3214", "214", "321", "4214", "314", "312", "3124", "31423", "34", "32141", "13214", "4214", "321412", "3321412", "32321414"};

        //DistanceTools dt = new DistanceTools(-1, -1, 1);
           String pattern = "ATGCTTGATTAAATGCTTGATTAA";
        String toCompare[] = {"TTGCTTGAATGCTTGATTAATTAA", "ATCCATGCTTGATTAATTGGTTAA", "ATGCTTGATTAAATGCTGGCTTAC",
            "GGGCTTATGCTTGATTAAGAAAAA", "GGGCAAATGCTTGATTAAGAGTAA", "TCCGTTGAATACATGCTTGATTAA",
            "ATGCTTGATTATATGCTTGTTAA", "CACCATGCATGCTTGATTAAACGA", "TACGAACATGCTTGATTAAATTTT",
            "ATGCTTGATTAAGAGCAACTCCTC", "TACGAACGAATAATGCTTGATTAA", "ATGCTTGATTAATACGAACTAATT",
            "ATGCTTGATTAAATGCTTGATTAA", "AtGCaTGCTTaaATGCTTGATTAA", "ATGCTTGATTAANNNNTTCATTAa"};
        for (int j = 0; j <= 2000000; j++) {
            for (int i = 0; i < toCompare.length; i++) {

                //float distance = dt.compute(pattern.toUpperCase(), toCompare[i].toUpperCase(), false);
                //int distance2 = dt.levenshteinDistance(pattern.toUpperCase(), toCompare[i].toUpperCase());
                int distance3 = dt.levenshteinDistanceReplacement(pattern.toUpperCase(), toCompare[i].toUpperCase());
                String color;
                if (i % 2 <= 0) {
                    color = "green";
                } else {
                    color = "blue";
                }
                //   System.out.println("|%{color:"+color+"}"+pattern+"%|/2. "+distance2+"|/2. "+distance3+"|/2. "+distance+"|");
                //  System.out.println("|%{color:"+color+"}"+toCompare[i]+"%|");
                //System.out.println(pattern + " vs " + toCompare[i] + " = "+ distance);
                //System.out.println("levenshteinDistance "+ pattern + " vs " + toCompare[i] + " = "+ distance2);
                //System.out.println("levenshtein Replacement "+ pattern + " vs " + toCompare[i] + " = "+ distance3);
                //System.out.println("**************************************");
                // System.exit(0);
            }
        }*/
        double end = (double) (System.currentTimeMillis() - start) / 1000;
       
        System.out.println("Time: " + end);
        System.out.println("**************************************");
    }

    private static String help() {
        String help = "\n#######################################################\n"
                + "###               Barcode corrector                 ###\n"
                + "###                    v 1.0                        ###\n"
                + "###                             @company       NIOZ ###\n"
                + "###                             @author   A. Abdala ###\n"
                + "#######################################################\n\n"
                + "Program implemented for Cascabel pipeline. It correct barcodes on a fastq file\n"
                + "using the  Levenshtein (edit) distance.\n"
                + "usage java BarcodeCorrector -f Fastq file -b Barcode file [options]\n\n"
                + "Mandatory arguments:\n"
                + "  -f\t--fastq           \tFastq file with the barcodes to correct (output from extract_barcodes.py from Qiime)\n"
                + "  -b\t--barcodes        \tTab separated file with barcode information in one of its columns [Def. 2]\n"
                + "\nOptional arguments:\n"
                + "  -c\t--barcode-column  \tThe column on the barcode file where the barcodes are found [Def. 2]\n"
                + "  -m\t--max-mismatch    \tMaximum number of mismatrches between a barcode in the fastq file and the mapping file to be found. [Def. 2]\n"
                + "  -o\t--output          \tName of the output file. If not supplied is equal to input file + the suffix \"_corrected\".\n"
                + "  -d\t--debug           \tOn debug mode 3 extra files are generated:'\n"
                + "                        \t  - File with all the corrected barcodes. \n"
                + "                        \t     * Format: original header (new line) original barcode (new line) new barcode.\n"
                + "                        \t     * File name equal to output +  \"_replacements\"\n"
                + "                        \t  - File with all the barcodes with ambiguity (two or more possible barcodes to be assigned).\n"
                + "                        \t     * Format: original header (new line) original barcode (new line) possible barcodes separated with new lines.\n"
                + "                        \t     * File name equal to output +  \"__ambiguos\"\n"
                + "                        \t  - File with all the barcodes that couldn't be assigned.\n"
                + "                        \t     * Format: fastq file with not_assigned barcodes.\n"
                + "                        \t     * File name equal to output +  \"__noMatch\"\n";

        return help;
    }
}
