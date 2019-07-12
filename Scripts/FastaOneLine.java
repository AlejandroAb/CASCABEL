/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author aabdala
 */
public class FastaOneLine {

    private String fasta = "/home/NIOZ.NL/aabdala/fasta.tmp";
    private int min_length = 1;
    private int max_length = 0;
    private String outfile = "";
    private boolean overwrite = true;
    private boolean writeDiscarded = false;

    public static void main(String args[]) {
        FastaOneLine fasta = null;
        if (args.length < 0) {
            printHelp();
            System.exit(0);
        }
        fasta = new FastaOneLine();
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-h") || args[i].equals("--help")) {
                printHelp();
                System.exit(0);
            } else if (args[i].equals("-f") || args[i].equals("--fasta")) {
                try {
                    i++;
                    fasta.setFasta(args[i]);
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -f --fasta \nNeed help? use FastaOneLine -h | --help");
                    System.exit(1);

                }

            } else if (args[i].equals("-o") || args[i].equals("--out")) {
                try {
                    i++;
                    fasta.setOutfile(args[i]);
                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -o --out option \nNeed help? use FastaOneLine -h | --help");
                    System.exit(1);

                }
            } else if (args[i].equals("-m") || args[i].equals("--minimum")) {
                try {
                    i++;
                    try {
                        fasta.setMin_length(Integer.parseInt(args[i]));
                    } catch (NumberFormatException nfe) {
                        System.err.println("NUmeric value expected for --minimum value. Found: " + args[i]);
                    }

                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -m --minimum option \nNeed help? use FastaOneLine -h | --help");
                    System.exit(1);

                }
            } else if (args[i].equals("-M") || args[i].equals("--maximum")) {
                try {
                    i++;
                    try {
                        fasta.setMax_length(Integer.parseInt(args[i]));
                    } catch (NumberFormatException nfe) {
                        System.err.println("NUmeric value expected for --maximum value. Found: " + args[i]);
                    }

                } catch (ArrayIndexOutOfBoundsException aoie) {
                    System.err.println("Argument expected for -m --maximum option \nNeed help? use FastaOneLine -h | --help");
                    System.exit(1);

                }
            } else if (args[i].equals("--write-discarded")) {
                fasta.setWriteDiscarded(true);
            }
        }
        if (fasta.getOutfile().length() > 1) {
            fasta.setOverwrite(false);
        }
        if (fasta.getFasta().length() < 2) {
            System.err.println("Please supply a fasta file to process\n"
                    + "Need help? use FastaOneLine -h | --help");
            System.exit(1);
        } else {
            File fasta_file = new File(fasta.getFasta());
            if (!fasta_file.exists()) {
                System.err.println("Please supply a valid fasta file\n"
                        + "File: " + fasta.getFasta() + " does not exist\n"
                        + "Need help? use FastaOneLine -h | --help");
                System.exit(1);
            } else {
                fasta.writeOneLiner();
            }
        }
        

    }

    public int getMax_length() {
        return max_length;
    }

    public void setMax_length(int max_length) {
        this.max_length = max_length;
    }

    public boolean isOverwrite() {
        return overwrite;
    }

    public void setOverwrite(boolean overwrite) {
        this.overwrite = overwrite;
    }

    public boolean isWriteDiscarded() {
        return writeDiscarded;
    }

    public void setWriteDiscarded(boolean writeDiscarded) {
        this.writeDiscarded = writeDiscarded;
    }

    /**
     * Print help menu
     */
    private static void printHelp() {
        String helpHeader = "\n#################################################################\n"
                + "###                     Fasta 2 One Line                      ###\n"
                + "###                                                    v 1.0  ###\n"
                + "###                                       @company       NIOZ ###\n"
                + "###                                       @author   A. Abdala ###\n"
                + "#################################################################\n";
        String description = "This program is created to convert fasta files with multiple lines\n"
                + "to single line fasta files. This program can also filter the sequences\n"
                + "based on the user selected minimum and maximum lengths.\n\n"
                + "Usage:\n\tjava FastaOneLine [options] -f <fasta_file>\n"
                + "The program will overwrite the multi line sequence fasta file into a single\n"
                + "sequence line fasta file.\n"
                + "In order to avoid input file overwriting, please supply any output file (-o)\n"
                + "Options:\n";

        String help = ""
                + "\t-f  --fasta    [file]   Fasta file to convert from to multiline to single line sequences. \n"
                + "\t-o  --out      [file]   Name of the output file, if this is not supplied, the input file will be\n"
                + "\t                        overwritten.\n"
                + "\t-m  --minimum  [number] Minimum length to retain any sequence. default = 0\n"
                + "\t-M  --maximum  [number] Maximum length to retain any sequence. default = any\n"
                + "\t--write-discarded       Write discarded reads into new file.\n"
                + "\t                        Discarded file will be called like the input + .discarded. default = false\n\n";

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

    public void writeOneLiner() {
        try {
            BufferedReader reader = new BufferedReader(new FileReader(fasta));
            FileWriter writer = null;
            if (overwrite) {
                writer = new FileWriter(fasta + ".tmp");
            } else {
                writer = new FileWriter(outfile);
            }
            FileWriter writerDiscard = null;
            if (writeDiscarded) {
                writerDiscard = new FileWriter(fasta + ".discarded");
            }
            String line;
            String header = "";
            StringBuilder sequence = new StringBuilder();
            boolean wasHeader = false; //if previous line wasHeader
            int lineas = 0;
            int sequences = 0;
            int headers = 0;
            int others = 0;
            int ok = 0;
            int nok = 0;
            int min = 0, max = 0;
            long start = System.currentTimeMillis();
            while ((line = reader.readLine()) != null) {
                lineas++;
                if (line.trim().startsWith(">")) {
                    headers++;
                    //si la anterior linea era un header
                    if (wasHeader) {
                        if (min_length == 0) {
                            writer.write(header + "\n");
                            ok++;
                        } else {
                            if (writeDiscarded) {
                                writerDiscard.write(header + "\n");
                            }
                            nok++;
                            min++;
                        }
                    }

                    if (!wasHeader && sequence.length() >= min_length && lineas > 1) {
                        if ((max_length > 0 && sequence.length() <= max_length) || max_length == 0) {
                            writer.write(header + "\n");
                            writer.write(sequence.toString() + "\n");
                            ok++;
                        } else {
                            nok++;
                            max++;
                            if (writeDiscarded) {
                                writerDiscard.write(header + "\n");
                                writerDiscard.write(sequence.toString() + "\n");
                            }
                        }
                        sequence = new StringBuilder();
                    } else if (!wasHeader && sequence.length() < min_length && lineas > 1) {
                        if (writeDiscarded) {
                            writerDiscard.write(header + "\n");
                        }
                        nok++;
                        min++;
			sequence = new StringBuilder();
                    }
                    header = line;
                    wasHeader = true;
                } else if (line.trim().length() > 0) {
                    sequence.append(line.trim());
                    wasHeader = false;
                    sequences++;
                } else {
                    others++;
                }
            }
            //saliendo del while flush the remaining sequence /  header
            if (wasHeader && min_length == 0) {
                writer.write(header + "\n");
                ok++;
            } else if (wasHeader) {
                if (writeDiscarded) {
                    writerDiscard.write(header + "\n");
                }
                nok++;
                min++;
            }
            if (!wasHeader && sequence.length() >= min_length) {
                if ((max_length > 0 && sequence.length() <= max_length) || max_length == 0) {
                    writer.write(header + "\n");
                    writer.write(sequence.toString() + "\n");
                    ok++;
                } else {
                    nok++;
                    max++;
                    if (writeDiscarded) {
                        writerDiscard.write(header + "\n");
                        writerDiscard.write(sequence.toString() + "\n");
                    }
                }
            } else if (!wasHeader && sequence.length() < min_length) {
                if (writeDiscarded) {
                    writerDiscard.write(header + "\n");
                }
                nok++;
                min++;
            }
            writer.close();
            if (writeDiscarded) {
                writerDiscard.close();
            }
            reader.close();
            if (overwrite) {
                File f = new File(fasta);
                File ftmp = new File(fasta + ".tmp");
                ftmp.renameTo(f);
            }
            long end = (System.currentTimeMillis() - start) / 1000;
            System.out.println("**************************");
            System.out.println("Time: " + end + "s.");
            System.out.println("Lines processed: " + lineas);
            System.out.println("Sequence headers found: " + headers);
            System.out.println("Lines containing sequences: " + sequences);
            System.out.println("Sequences writted: " + ok);
            System.out.println("Sequences smaller than " + min_length + ": " + min);
            if (max_length > 0) {
                System.out.println("Sequences larger than " + max_length + ": " + max);
            }
            System.out.println("Blank / Other lines:" + others);
            System.out.println("***********END************");

        } catch (FileNotFoundException ex) {
            Logger.getLogger(FastaOneLine.class.getName()).log(Level.SEVERE, null, ex);
            System.err.println("Error reading file: " + fasta);
            System.exit(1);
        } catch (IOException ex) {
            Logger.getLogger(FastaOneLine.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public String getFasta() {
        return fasta;
    }

    public void setFasta(String fasta) {
        this.fasta = fasta;
    }

    public int getMin_length() {
        return min_length;
    }

    public void setMin_length(int min_length) {
        this.min_length = min_length;
    }

    public String getOutfile() {
        return outfile;
    }

    public void setOutfile(String outfile) {
        this.outfile = outfile;
    }

}
