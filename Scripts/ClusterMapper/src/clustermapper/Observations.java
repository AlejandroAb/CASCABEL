/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package clustermapper;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author aabdala
 */
public class Observations {

    private HashMap<String, Integer> map = new HashMap<>();
    private ArrayList<String> sequences = new ArrayList<>();
    private String seed_sequence;
    public Observations() {
    }

    public String getSeedSequence() {
        return seed_sequence;
    }

    public void setSeedSequence(String seed) {
        this.seed_sequence = seed;
    }

    public void addSequence(String sequence) {
        sequences.add(sequence);
    }

    public ArrayList<String> getSequences() {
        return sequences;
    }

    public void setSequences(ArrayList<String> sequences) {
        this.sequences = sequences;
    }
    

    /**
     * En teoria el merge busca el elemento si existe aplica la funcion, en este
     * caso sum y le suma la abundancia. Si la llave no existe, la creo.
     *
     * @param idSample
     * @param abundance
     */
    public void addSampleObservation(String idSample, int abundance) {
        map.merge(idSample, abundance, Integer::sum);
    }

    public void mergeObservations(Observations toMerge) {
        for (Map.Entry<String, Integer> entry : toMerge.getMap().entrySet()) {
            map.merge(entry.getKey(), entry.getValue(), Integer::sum);
        }
        this.sequences.addAll(toMerge.getSequences());
    }

    /**
     * Return the abundance given a sample.
     *
     * @return the abundance of the sample, 0 if is not present
     */
    public int getAbundanceBySample(String sample) {
        //Integer t = map.get(sample);
        return map.getOrDefault(sample, 0);
    }

    /*public void addSampleObservation(String idSample, String abundance) {
        for (SampleObs ob : obs) {
            if (ob.getId().equals(idSample)) {
                ob.incrementStrAbundance(abundance);
            }
        }
        this.obs.add(o);
    }

    public SampleObs getObservation(String id) {

        for (SampleObs o : obs) {
            if (o.getId().equals(id)) {
                return o;
            }
        }
        return null;
    }*/
    public HashMap<String, Integer> getMap() {
        return map;
    }

    public void setMap(HashMap<String, Integer> map) {
        this.map = map;
    }
}
