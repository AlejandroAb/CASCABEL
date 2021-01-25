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
public class DistanceTools {

    private float gapValue = -.5f;
    private float notEqualsValue = -.5f;
    private float equalValue = 1;

    public DistanceTools() {
    }

    public DistanceTools(float gapValue, float notEqualsValue, float equalValue) {
        this.gapValue = gapValue;
        this.notEqualsValue = notEqualsValue;
        this.equalValue = equalValue;
    }

    //defs:-.5f,-.5f,1
    public float getGapValue() {
        return gapValue;
    }

    public void setGapValue(float gapValue) {
        this.gapValue = gapValue;
    }

    public float getNotEqualsValue() {
        return notEqualsValue;
    }

    public void setNotEqualsValue(float notEqualsValue) {
        this.notEqualsValue = notEqualsValue;
    }

    public float getEqualValue() {
        return equalValue;
    }

    public void setEqualValue(float equalValue) {
        this.equalValue = equalValue;
    }

    public float maximum(float a, float b, float c) {
        return Math.max(Math.max(a, b), c);
    }

    public String getPath(float u, float l, float d) {
        if (d >= u) {
            if (d >= l) {
                return "D";
            } else {
                return "L";
            }
        } else if (l >= u) {
            return "L";
        } else {
            return "U";
        }
    }

    public float compute(CharSequence str1,
            CharSequence str2, boolean debug) {

        float[][] distance = new float[str1.length() + 1][str2.length() + 1];
        String[][] camino = new String[str1.length() + 1][str2.length() + 1];
        for (int i = 0; i <= str1.length(); i++) {
            if (i == 0) {
                distance[i][0] = 0;
                camino[i][0] = "E";
            } else if (i == 1) {
                distance[i][0] = 0;
                camino[i][0] = "T";
            } else if (i == str1.length()) {
                distance[i][0] = distance[i - 1][0] - gapValue;
                camino[i][0] = "T";
            } else {
                distance[i][0] = distance[i - 1][0] + gapValue;
                camino[i][0] = "T";
            }
        }

        for (int j = 0; j <= str2.length(); j++) {

            if (j == 0) {
                camino[0][j] = "E";
                distance[0][j] = 0;
            } else if (j == 1) {
                camino[0][j] = "L";
                distance[0][j] = 0;
            } else if (j == str2.length()) {
                camino[0][j] = "L";
                distance[0][j] = distance[0][j - 1] - gapValue;;
            } else {
                distance[0][j] = distance[0][j - 1] + gapValue;
                camino[0][j] = "L";
            }

        }

        for (int i = 1; i <= str1.length(); i++) {
            for (int j = 1; j <= str2.length(); j++) {
                //para que gaps o permutaciones en los extremos cuesten o no
                if ((i == 1 && j == 1) || (i == str1.length() && j == str2.length())) {
                    distance[i][j] = maximum(
                            distance[i - 1][j], //ya no se suma gap para quitar costo
                            distance[i][j - 1],
                            distance[i - 1][j - 1]
                            + ((str1.charAt(i - 1) == str2.charAt(j - 1)) ? equalValue
                            : notEqualsValue));//0 para quitar costo a cambio
                    camino[i][j] = getPath(
                            distance[i - 1][j],
                            distance[i][j - 1],
                            distance[i - 1][j - 1]
                            + ((str1.charAt(i - 1) == str2.charAt(j - 1)) ? equalValue
                            : notEqualsValue));
                } else {
                    distance[i][j] = maximum(
                            distance[i - 1][j] + gapValue,
                            distance[i][j - 1] + gapValue,
                            distance[i - 1][j - 1]
                            + ((str1.charAt(i - 1) == str2.charAt(j - 1)) ? equalValue
                            : notEqualsValue));
                    camino[i][j] = getPath(
                            distance[i - 1][j] + gapValue,
                            distance[i][j - 1] + gapValue,
                            distance[i - 1][j - 1]
                            + ((str1.charAt(i - 1) == str2.charAt(j - 1)) ? equalValue
                            : notEqualsValue));
                }
            }
        }
        if (debug) {
            System.out.append(str1 + " vs " + str2 + "\n");
            for (int i = 0; i < distance.length; i++) {
                for (int j = 0; j < distance[0].length; j++) {
                    System.out.print(distance[i][j] + "\t");
                }
                System.out.println();
            }
            for (int i = 0; i < camino.length; i++) {
                for (int j = 0; j < camino[0].length; j++) {
                    System.out.print(camino[i][j] + "\t");
                }
                System.out.println();
            }

        }
        return distance[str1.length()][str2.length()];
    }

    /**
     *
     * @param lhs
     * @param rhs
     * @return
     */
    public int levenshteinDistance(CharSequence lhs, CharSequence rhs) {
        int len0 = lhs.length() + 1;
        int len1 = rhs.length() + 1;

        // the array of distances                                                       
        int[] cost = new int[len0];
        int[] newcost = new int[len0];

        // initial cost of skipping prefix in String s0                                 
        for (int i = 0; i < len0; i++) {
            cost[i] = i;
        }

        // dynamically computing the array of distances                                  
        // transformation cost for each letter in s1                                    
        for (int j = 1; j < len1; j++) {
            // initial cost of skipping prefix in String s1                             
            newcost[0] = j;

            // transformation cost for each letter in s0                                
            for (int i = 1; i < len0; i++) {
                // matching current letters in both strings                             
                int match = (lhs.charAt(i - 1) == rhs.charAt(j - 1)) ? 0 : 1;

                // computing cost for each transformation                               
                int cost_replace = cost[i - 1] + match;
                int cost_insert = cost[i] + 1;
                int cost_delete = newcost[i - 1] + 1;

                // keep minimum cost                                                    
                newcost[i] = Math.min(Math.min(cost_insert, cost_delete), cost_replace);
            }

            // swap cost/newcost arrays                                                 
            int[] swap = cost;
            cost = newcost;
            newcost = swap;
        }

        // the distance is the cost for transforming all letters in both strings        
        return cost[len0 - 1];
    }
/**
 * Like levenshteinDistance implementation but only allowing replacement. 
 * @param lhs
 * @param rhs
 * @return 
 */
    public int levenshteinDistanceReplacement(CharSequence lhs, CharSequence rhs) {
        int len0 = lhs.length() + 1;
        int len1 = rhs.length() + 1;

        // the array of distances                                                       
        int[] cost = new int[len0];
        int[] newcost = new int[len0];

        // initial cost of skipping prefix in String s0                                 
        for (int i = 0; i < len0; i++) {
            cost[i] = i;
        }

        // dynamically computing the array of distances                                  
        // transformation cost for each letter in s1                                    
        for (int j = 1; j < len1; j++) {
            // initial cost of skipping prefix in String s1                             
            newcost[0] = j;

            // transformation cost for each letter in s0                                
            for (int i = 1; i < len0; i++) {
                // matching current letters in both strings                             
                int match = (lhs.charAt(i - 1) == rhs.charAt(j - 1)) ? 0 : 1;

                // computing cost for each transformation                               
                int cost_replace = cost[i - 1] + match;
                int cost_insert = cost[i] + 1;
                int cost_delete = newcost[i - 1] + 1;

                // keep minimum cost                                                    
               // newcost[i] = Math.min(Math.min(cost_insert, cost_delete), cost_replace);
               newcost[i] = cost_replace;
            }

            // swap cost/newcost arrays                                                 
            int[] swap = cost;
            cost = newcost;
            newcost = swap;
        }

        // the distance is the cost for transforming all letters in both strings        
        return cost[len0 - 1];
    }
}
