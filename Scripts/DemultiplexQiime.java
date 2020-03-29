import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author aabdala
 * @version 1.0
 */
public class DemultiplexQiime {

    private String fastqFile;
    private String outDir = "";
    private String suffix = "";
    private String preffix = "";
    private boolean debug = false;
    private String forward = "";
    private String reverse = "";
    private String paired = "";
    private boolean fastq = true;
    private boolean raw_fastq = true;
    private String orientation = "fw_rv";
    int bc_length = 0;
    private boolean removeHeader = false;
    private boolean overWrite = false;
    private boolean writeSummaryFiles = false;

    public static void main(String args[]) {
        DemultiplexQiime dmx = null;
        if (args.length < 1) {
            printHelp();
            System.exit(0);
        }
        dmx = new DemultiplexQiime();
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-h") || args[i].equals("--help")) {
                printHelp();
                System.exit(0);
            } else if (args[i].equals("-d") || args[i].equals("--dmx-file")) {
                try {
                    i++;
                    dmx.setFastqFile(args[i]);
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -d -dmx-file \nNeed help? use DemultiplexQiime -h | --help");
                    System.exit(1);

                }

            } else if (args[i].equals("-o") || args[i].equals("--out")) {
                try {
                    i++;
                    dmx.setOutDir(args[i]);
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for --out option \nNeed help? use DemultiplexQiime -h | --help");
                    System.exit(1);

                }
            } else if (args[i].equals("-s") || args[i].equals("--sufix")) {
                try {
                    i++;
                    dmx.setSuffix(args[i]);
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -s --sufix \nNeed help? use DemultiplexQiime -h | --help");
                    System.exit(1);

                }

            } else if (args[i].equals("-p") || args[i].equals("--prefix")) {
                try {
                    i++;
                    dmx.setPreffix(args[i]);
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -p --prefix \nNeed help? use DemultiplexQiime -h | --help");
                    System.exit(1);

                }

            } else if (args[i].equals("-r1") || args[i].equals("--forward")) {
                try {
                    i++;
                    dmx.setForward(args[i]);
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -r1 --forward \nNeed help? use DemultiplexQiime -h | --help");
                    System.exit(1);

                }

            } else if (args[i].equals("-r2") || args[i].equals("--reverse")) {
                try {
                    i++;
                    dmx.setReverse(args[i]);
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -r2 --reverse \nNeed help? use DemultiplexQiime -h | --help");
                    System.exit(1);

                }

            } else if (args[i].equals("-r") || args[i].equals("--single")) {
                try {
                    i++;
                    dmx.setPaired(args[i]);
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -r --single \nNeed help? use DemultiplexQiime -h | --help");

                }

            } else if (args[i].equals("-a") || args[i].equals("--assignation")) {
                try {
                    i++;
                    dmx.setOrientation(args[i]);
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -s --strand \nNeed help? use DemultiplexQiime -h | --help");

                }

            } else if (args[i].equals("-b") || args[i].equals("--barcodes")) {
                try {
                    i++;
                    dmx.setBc_length(Integer.parseInt(args[i]));
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -b --barcodes \nNeed help? use DemultiplexQiime -h | --help");

                } catch (NumberFormatException nfe) {
                    System.err.println("Numerical argument expected for -b --barcodes. Found: " + args[i] + "\nNeed help? use DemultiplexQiime -h | --help");

                }

            } else if (args[i].equals("--fasta")) {
                dmx.setFastq(false);
            } else if (args[i].equals("--out-fasta")) {
                dmx.setRaw_fastq(false);
            } else if (args[i].equals("-v") || args[i].equals("--verbose")) {
                dmx.setDebug(true);
            } else if (args[i].equals("--remove-header")) {
                dmx.setRemoveHeader(true);
            } else if (args[i].equals("--over-write")) {
                dmx.setOverWrite(true);
            } else if (args[i].equals("--write-summary")) {
                dmx.setWriteSummaryFiles(true);
            } else {
                System.err.println("Unknown option: " + args[i] + "\nNeed help? use DemultiplexQiime -h | --help");
                System.exit(1);
            }
        }

        if (dmx.getFastqFile().length() < 1) {
            System.err.println("Please supply a fasta / fastq file to process\n"
                    + "Need help? use DemultiplexQiime -h | --help");
            System.exit(1);

        }
        if (dmx.getOutDir().length() > 0 && !dmx.getOutDir().endsWith("/")) {
            dmx.setOutDir(dmx.getOutDir() + "/");
        }
        if (dmx.getForward().length() > 0 && dmx.getReverse().length() > 0) {
            File fw = new File(dmx.getForward());
            if (!fw.exists()) {
                System.err.println("Please supply a valid forward read file\n"
                        + "File: " + dmx.getForward() + " does not exist\n"
                        + "Need help? use DemultiplexQiime -h | --help");
                System.exit(1);
            }
            File rv = new File(dmx.getReverse());
            if (!rv.exists()) {
                System.err.println("Please supply a valid reverse read file\n"
                        + "File: " + dmx.getReverse() + " does not exist\n"
                        + "Need help? use DemultiplexQiime -h | --help");
                System.exit(1);
            }
            dmx.processFastQMap2Raw();
        } else if (dmx.getPaired().length() > 0) {
            File single = new File(dmx.getPaired());
            if (!single.exists()) {
                System.err.println("Please supply a valid single read file\n"
                        + "File: " + dmx.getPaired() + " does not exist\n"
                        + "Need help? use DemultiplexQiime -h | --help");
                System.exit(1);
            }
            dmx.processFastQMap2Raw();
        } else {
            dmx.processFastQ();
        }

    }

    /**
     * Demultiplex fastq file according to acc on sequence level. In order to
     * demultiplex the fastq file, the header or sequence tittle should be in a
     * qiime format: @<SAMPLE_NAME>_# old raw header This will create as many
     * fastq files as different <SAMPLE_NAME>s are found on the fastq file
     */
    public void processFastQ() {
        HashMap<String, FileWriter> writers = new HashMap<>();
        HashMap<String, Integer> counts = new HashMap<>();

        try {
            BufferedReader reader = new BufferedReader(new FileReader(fastqFile));
            String linea;
            int numL = 0;
            String nl = System.getProperty("line.separator");
            FileWriter currentWriter = null;
            int factor = fastq ? 4 : 2;
            String file_extenssion = fastq ? "fastq" : "fasta";
            while ((linea = reader.readLine()) != null) {
                String sample;
                if (numL % factor == 0 && linea.length() > 0) {
                    sample = linea.substring(1, linea.indexOf("_"));
                    currentWriter = writers.get(sample);
                    if (currentWriter == null) {
                        // currentWriter = new FileWriter(outDir + preffix + sample + suffix + "." + file_extenssion);
                        currentWriter = new FileWriter(outDir + preffix + sample + suffix + "_" + orientation + ".txt");

                        writers.put(sample, currentWriter);
                        counts.put(sample, 1);
                    } else {
                        Integer tmp = counts.get(sample);
                        tmp++;
                        counts.put(sample, tmp);
                    }
                    currentWriter.write(linea + nl);
                } else {
                    currentWriter.write(linea + nl);
                }
                numL++;
            }
            reader.close();
            for (String key : writers.keySet()) {
                writers.get(key).close();
            }
            FileWriter writer = new FileWriter(outDir + "summary.txt");
            writer.write("#Sample\tCounts\tFile(s)\n");
            for (String key : counts.keySet()) {
                writer.write(key + "\t" + counts.get(key) + "\t" + outDir + preffix + key + suffix + "." + file_extenssion + "\n");
            }
            writer.close();

        } catch (FileNotFoundException ex) {
            Logger.getLogger(DemultiplexQiime.class.getName()).log(Level.SEVERE, null, ex);
            System.err.println("There is no file: " + fastqFile);
            System.err.println("Program aborting");
            System.exit(1);
        } catch (IOException ioe) {
            Logger.getLogger(DemultiplexQiime.class.getName()).log(Level.SEVERE, null, ioe);
            System.err.println("Error reading file: " + fastqFile);
            System.err.println("Program aborting");
        }
    }

    public void processFastQMap2Raw() {
        HashMap<String, FileWriter> writers = new HashMap<>();
        HashMap<String, Integer> counts = new HashMap<>();
        ArrayList<String> listFiles = new ArrayList<>();

        try {
            BufferedReader reader = new BufferedReader(new FileReader(fastqFile));
            String linea;
            int numL = 0;
            String nl = System.getProperty("line.separator");
            FileWriter currentWriter = null;
            int factor = fastq ? 4 : 2;
            while ((linea = reader.readLine()) != null) {
                try {
                    String sample;
                    String raw_header;
                    if (numL % factor == 0 && linea.length() > 0) {
                        sample = linea.substring(1, linea.indexOf("_"));
                        raw_header = linea.substring(linea.indexOf(" ") + 1);//first removes the QIIME acc
                        raw_header = raw_header.substring(0, raw_header.indexOf(" ")); //the only keep the raw acc
                        currentWriter = writers.get(sample);
                        if (currentWriter == null) {
                            currentWriter = new FileWriter(outDir + preffix + sample + suffix + "_" + orientation + ".txt");
                            writers.put(sample, currentWriter);
                            counts.put(sample, 1);
                            listFiles.add(outDir + preffix + sample + suffix + "_" + orientation + ".txt");
                        } else {
                            Integer tmp = counts.get(sample);
                            tmp++;
                            counts.put(sample, tmp);
                        }
                        currentWriter.write(raw_header + nl);
                    }
                    /*else {
                    currentWriter.write(linea + nl);
                }*/
                    numL++;
                } catch (Exception ex) {
                    System.err.println("Error processing line: " + numL++);
                    ex.printStackTrace();
                    System.exit(1);
                    reader.close();
                }
            }
            reader.close();
            for (String key : writers.keySet()) {
                writers.get(key).close();
                //String file = outDir + preffix + key + suffix + ".txt";
                String file = outDir + preffix + key + suffix + "_" + orientation + ".txt";
                if (forward.length() > 0) {
                    createFastq(key, file, forward, "_1");
                }
                if (reverse.length() > 0) {
                    createFastq(key, file, reverse, "_2");
                }
                if (paired.length() > 0) {
                    createFastq(key, file, paired, "");
                }
            }
            //Delete files unless --write-summary
            if (!writeSummaryFiles) {
                for (String fname : listFiles) {
                    File f = new File(fname);
                    f.delete();
                }

            }
            FileWriter writer = new FileWriter(outDir + "summary_" + orientation + ".txt");
            writer.write("#Samples assigned in " + orientation + " orientation \n");
            writer.write("#Sample\tCounts\tFile(s)\n");
            String file_extenssion = raw_fastq ? "fastq" : "fasta";
            for (String key : counts.keySet()) {
                writer.write(key + "\t" + counts.get(key) + "\t");
                if (forward.length() > 0) {
                    writer.write(outDir + preffix + key + suffix + "_1" + "." + file_extenssion + ".gz");
                }
                if (reverse.length() > 0) {
                    writer.write(";" + outDir + preffix + key + suffix + "_2" + "." + file_extenssion + ".gz\n");
                } else if (paired.length() > 0) {
                    writer.write(outDir + preffix + key + suffix + "" + "." + file_extenssion + ".gz\n");
                }
            }
            writer.close();

        } catch (FileNotFoundException ex) {
            Logger.getLogger(DemultiplexQiime.class.getName()).log(Level.SEVERE, null, ex);
            System.err.println("There is no file: " + fastqFile);
            System.err.println("Program aborting");
            System.exit(1);
        } catch (IOException ioe) {
            Logger.getLogger(DemultiplexQiime.class.getName()).log(Level.SEVERE, null, ioe);
            System.err.println("Error reading file: " + fastqFile);
            System.err.println("Program aborting");
        }
    }

    /**
     * This method creates a fasta | fastq files by fetching sequences from
     * original raw sequences
     *
     * @param sample the id of any sample on the library
     * @param ids the file containing the accessions belonging to such samples
     * @param raw_file the raw file with the original sequences before any
     * quality trimming
     * @param extra extra string to put to the names of the resultant files, for
     * example for forward reads "_1" or "_2" for reverse reads or "" for single
     * libraries
     * @return
     */
    public String createFastq(String sample, String ids, String raw_file, String extra) {
        String command = "";
        String file_extenssion = raw_fastq ? "fastq" : "fasta";
        int start_idx = bc_length + 1;
        String header = removeHeader ? "$1" : "$0";//all the header
        String append = overWrite ? ">" : ">>";
        try {
            if (start_idx > 1 && raw_fastq) {
                command = "seqtk subseq  " + raw_file + " " + ids
                        + " | awk '{if(NR%4==2 || NR%4==0){print substr($0," + start_idx + ")}else if(NR%4==1){print " + header + "}else{print $0}}' "
                        + " | gzip " + append + " " + outDir + preffix + sample + suffix + extra + "." + file_extenssion + ".gz";
            } else if (start_idx > 1) { //is fasta file
                command = "seqtk subseq  " + raw_file + " " + ids
                        + " | awk '{if(NR%2==0){print substr($0," + start_idx + ")}else if(NR%4==1){print " + header + "}else{print $0}}' "
                        + " | gzip " + append + " " + outDir + preffix + sample + suffix + extra + "." + file_extenssion + ".gz";
            } else { //no bc removal 
                command = "seqtk subseq  " + raw_file + " " + ids
                        + " | gzip " + append + " " + outDir + preffix + sample + suffix + extra + "." + file_extenssion + ".gz";
            }

            if (debug) {
                System.out.println("command: " + command);
            }
            String cmdArr[] = {"/bin/sh", "-c", command};
            Process proc = Runtime.getRuntime().exec(cmdArr);
            //Process proc = Runtime.getRuntime().exec("c:/Program Files/R/R-3.0.3/bin/Rscript " + scriptDir + "script_abundancia.R " + workingDir +" "+ scriptDir + " " +file + " " + imageName);
            proc.waitFor();
            InputStream inputstream = proc.getInputStream();
            InputStreamReader inputstreamreader = new InputStreamReader(inputstream);
            BufferedReader bufferedreader = new BufferedReader(inputstreamreader);
            String line = "";

            while ((line = bufferedreader.readLine()) != null) {
                if (debug) {
                    System.out.println(line);
                }
            }
            //
            bufferedreader.close();
            if (proc.waitFor() != 0) {
                System.err.println("exit value = " + proc.exitValue());
                System.err.println("command = " + command);
                writeErrStream(proc.getErrorStream());
                proc.destroy();
                return "Error";
            } else {
                proc.destroy();
                return "command: " + command;
            }
        } catch (Exception e) {
            System.out.println("" + e.getLocalizedMessage());
            e.printStackTrace();
            return "Error runing: " + command;
        }
    }

    /**
     * Print help menu
     */
    private static void printHelp() {
        String helpHeader = "\n#################################################################\n"
                + "###                      Sequence Demultiplexer               ###\n"
                + "###                            v 1.1                          ###\n"
                + "###                                       @company       NIOZ ###\n"
                + "###                                       @author   A. Abdala ###\n"
                + "#################################################################\n";
        String description = "This program is created to demultiplex sample labeled Qiime like fastq|fasta files.\n"
                + "The expected format for the sequence headers should be:\n"
                + "@|><SAMPLE_NAME>_# <original fastq raw reads header> example:\n"
                + "@NIOZ80.39_5 M01102:251:000000000-B8LBW:1:1101:17625:1850\n\n"
                + "Usage:\n\tjava DemultiplexQiime [options] -d <demultiplexed_file>\n"
                + "The program will write individual demultiplexed files into -o <output_dir>\n"
                + "And a tab delimited file at <output_dir>/summary_<orientation>.txt  with the following columns:\n"
                + "#Sample<tab>Counts<tab>File(s)\n\n"
                + "Options:\n";

        String help = ""
                + "  -d  --dmx-file <file>    Fastq/Fasta file with the the sequences already demultiplexed on sequences header. \n"
                + "      --fasta              If the input file (-d) is on fasta format pass this flag (defult input is fastq)\n"
                + "  -o  --out <dir>          Base directory to write the demultiplexed files\n"
                + "      --out-fasta          If the output will be on fasta format pass this flag to use propper extension\n"
                + "      --over-write         Overwrite generated files if exists at output directory. Default False\n"
                + "  -s  --sufix <str>        Sufix to output file names. Default <blank> no sufix\n"
                + "  -p  --prefix <str>       Prefix to output file names. Default <blank> no prefix\n"
                + "                           *The individual demultiplexed files will be called:\n"
                + "                           <base_directory>/<prefix><sample><sufix>.<fastq|fasta>.gz\n"
                + "  -r1 --forward <file>     Forward raw reads to extract the sequences according to the demultiplexed file (-d)\n"
                + "  -r2 --reverse <file>     Reverse raw reads to extract the sequences according to the demultiplexed file (-d)\n"
                + "  -r  --single <file>      Single raw reads to extract the sequences according to the demultiplexed file (-d)\n"
                + "                           ----------------------\n"
                + "                           *Specify -r1 <fw.reads> -r2 <rv.reads> for paired end data.\n"
                + "                           *Specify -r <single.reads> for single end data.\n"
                + "                           *Dont specify any reads and the demultiplex will be performed only with the input file (-d)\n"
                + "                           ----------------------\n"
                + "  -b  --barcodes <num>     Remove barcodes from raw reads. Default 0 -removes 0 bases-\n"
                + "  -a  --assignation <str>  This option refers to the orientation of the reads at the moment of its assignation.\n"
                + "                           This value would be used only on the summary file and it is used to distinguish between strand\n"
                + "                           specific configurations. Default fw_rev\n"
                + "      --remove-header      Pass this flag to only keep the first part of the header.\n"
                + "                            i.e., @M01102:251:000000000-B8LBW:1:2110:14503:27454 1:N:0:ATTCCT\n"
                + "                            to    @M01102:251:000000000-B8LBW:1:2110:14503:27454\n"
                + "                            This option is useful when the reads FW and RV are interchanged between raw and demultiplexed\n"
                + "                            files, in order to have all of them in the correct strand.\n"
                + "     --write-summary       This option will generate one file per sample, per strand with the headers of the original fastq\n"
                + "                           that belogns to such sample.\n"
                + "  -v  --verbose            Output details during processing\n";

        System.out.println(helpHeader + description + help);
    }

    public void writeErrStream(InputStream stream) throws IOException {
        InputStreamReader inputstreamreader = new InputStreamReader(stream);
        BufferedReader bufferedreader = new BufferedReader(inputstreamreader);
        String line = "";
        while ((line = bufferedreader.readLine()) != null) {
            System.out.println(line);

        }
    }

    public boolean isDebug() {
        return debug;
    }

    public void setDebug(boolean debug) {
        this.debug = debug;
    }

    public String getForward() {
        return forward;
    }

    public void setForward(String forward) {
        this.forward = forward;
    }

    public String getReverse() {
        return reverse;
    }

    public boolean isFastq() {
        return fastq;
    }

    public void setFastq(boolean fastq) {
        this.fastq = fastq;
    }

    public void setReverse(String reverse) {
        this.reverse = reverse;
    }

    public String getPaired() {
        return paired;
    }

    public void setPaired(String paired) {
        this.paired = paired;
    }

    public String getSuffix() {
        return suffix;
    }

    public void setSuffix(String suffix) {
        this.suffix = suffix;
    }

    public String getPreffix() {
        return preffix;
    }

    public void setPreffix(String preffix) {
        this.preffix = preffix;
    }

    public String getFastqFile() {
        return fastqFile;
    }

    public void setFastqFile(String fastqFile) {
        this.fastqFile = fastqFile;
    }

    public String getOutDir() {
        return outDir;
    }

    public void setOutDir(String outDir) {
        this.outDir = outDir;
    }

    public DemultiplexQiime(String fastqFile) {
        this.fastqFile = fastqFile;

    }

    public DemultiplexQiime() {

    }

    public boolean isRaw_fastq() {
        return raw_fastq;
    }

    public void setRaw_fastq(boolean raw_fastq) {
        this.raw_fastq = raw_fastq;
    }

    public String getOrientation() {
        return orientation;
    }

    public void setOrientation(String orientation) {
        this.orientation = orientation;
    }

    public int getBc_length() {
        return bc_length;
    }

    public void setBc_length(int bc_length) {
        this.bc_length = bc_length;
    }

    public boolean isRemoveHeader() {
        return removeHeader;
    }

    public void setRemoveHeader(boolean removeHeader) {
        this.removeHeader = removeHeader;
    }

    public boolean isOverWrite() {
        return overWrite;
    }

    public void setOverWrite(boolean overWrite) {
        this.overWrite = overWrite;
    }

    public boolean isWriteSummaryFiles() {
        return writeSummaryFiles;
    }

    public void setWriteSummaryFiles(boolean writeSummaryFiles) {
        this.writeSummaryFiles = writeSummaryFiles;
    }

}

