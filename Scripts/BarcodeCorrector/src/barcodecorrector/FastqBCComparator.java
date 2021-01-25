/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package barcodecorrector;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Alegandro
 */
public class FastqBCComparator {

    private String fastqFile;
    private String barcodes;
    private boolean debug = false;
    private int bc_column = 2;

    public boolean isDebug() {
        return debug;
    }

    public int getBc_column() {
        return bc_column;
    }

    public void setBc_column(int bc_column) {
        this.bc_column = bc_column;
    }

    public void setDebug(boolean debug) {
        this.debug = debug;
    }

    public String getFastqFile() {
        return fastqFile;
    }

    public void setFastqFile(String fastqFile) {
        this.fastqFile = fastqFile;
    }

    public String getBarcodes() {
        return barcodes;
    }

    public void setBarcodes(String barcodes) {
        this.barcodes = barcodes;
    }

    public FastqBCComparator(String fastqFile, String barcodes) {
        this.fastqFile = fastqFile;
        this.barcodes = barcodes;
    }

    public void correctBarcodes(int maxMismatch, String outputFile) {
        try {
            BufferedReader fqreader = new BufferedReader(new FileReader(fastqFile));
            BufferedReader barcodeReader = new BufferedReader(new FileReader(barcodes));
            FileWriter outputWriter = new FileWriter(outputFile);
            FileWriter replacements = null;
            FileWriter ambiguos = null;
            FileWriter notMatch = null;
            if (debug) {
                replacements = new FileWriter(outputFile + "_replacements");
                ambiguos = new FileWriter(outputFile + "_ambiguos");
                notMatch = new FileWriter(outputFile + "_noMatch");
            }
            String bcLine = "";
            ArrayList<String> barcodes = new ArrayList();
            DistanceTools dt = new DistanceTools();
            while ((bcLine = barcodeReader.readLine()) != null) {
                String bcRow[] = bcLine.split("\t");
                if (bcLine.trim().length() > 1 && !bcLine.startsWith("#") && bcRow.length > (bc_column - 1)) {

                    //fixed to second column of mapping file, later we can improve this
                    barcodes.add(bcRow[(bc_column - 1)]);
                }
            }
            barcodeReader.close();
            String fqLine = "";
            double lidx = 0;
            ArrayList<String> optBarcodes = new ArrayList();
            int perfect_match = 0;
            int replaced = 0;
            int ambiguity = 0;
            int not_match = 0;
            int reads = 0;
            String header = "";
            String quality = "";
            boolean match = true;
            while ((fqLine = fqreader.readLine()) != null) {
                lidx++;
                if (lidx % 4 == 2) {
                    match = true;
                    int distance = maxMismatch + 1;
                    boolean writen = false;
                    reads++;
                    for (int i = 0; i < barcodes.size(); i++) {
                        int newd = dt.levenshteinDistance(fqLine, barcodes.get(i));
                        if (newd == 0) {
                            outputWriter.write(fqLine + "\n");
                            writen = true;
                            perfect_match++;
                            break;
                        } else if (newd <= maxMismatch) {
                            if (newd < distance) {
                                distance = newd;
                                optBarcodes = new ArrayList();
                                optBarcodes.add(barcodes.get(i));
                            } else if (newd == distance) {
                                optBarcodes.add(barcodes.get(i));
                            }
                        }
                    }
                    if (!writen) {
                        //we only write the barcode if there is one and only one option
                        if (optBarcodes.size() == 1) {
                            outputWriter.write(optBarcodes.get(0) + "\n");
                            if (debug) {
                                replacements.write(header + "\n");
                                replacements.write(fqLine + "\n");
                                replacements.write(optBarcodes.get(0) + "\n");
                            }
                            replaced++;
                        } else if (optBarcodes.size() > 1) {
                            ambiguity++;
                            if (debug) {
                                ambiguos.write(header + "\n");
                                ambiguos.write(fqLine + "\n");
                                for (int j = 0; j < optBarcodes.size(); j++) {
                                    ambiguos.write(optBarcodes.get(j) + "\n");
                                }
                            }
                            outputWriter.write(fqLine + "\n");
                        } else {
                            not_match++;
                            if (debug) {
                                notMatch.write(header + "\n");
                                notMatch.write(fqLine + "\n");
                                notMatch.write("+\n");
                                //notMatch.write(quality + "\n");
                                match = false;
                            }
                            outputWriter.write(fqLine + "\n");
                        }
                    }

                } else {
                    if (debug && lidx % 4 == 1) {
                        header = fqLine;
                    }
                    if (debug && lidx % 4 == 0 && !match) {
                        // quality = fqLine;
                        if (debug) {
                           // notMatch.write(header + "\n");
                            //notMatch.write(fqLine + "\n");
                            //notMatch.write("+\n");
                            notMatch.write(fqLine + "\n");
                            //match = false;
                        }
                    }
                    outputWriter.write(fqLine + "\n");
                }
                optBarcodes = new ArrayList();

            }
            outputWriter.close();
            fqreader.close();
            if (debug) {
                notMatch.close();
                ambiguos.close();
                replacements.close();
            }
            System.out.println("Reads processed: " + reads);
            System.out.println("Reads with perfect match: " + perfect_match);
            System.out.println("Reads with barcode corrected: " + replaced);
            System.out.println("Reads with barcode ambiguity: " + ambiguity);
            System.out.println("Reads without a match: " + not_match);

        } catch (FileNotFoundException ex) {
            Logger.getLogger(FastqBCComparator.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(FastqBCComparator.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
