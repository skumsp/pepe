/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepe;

import ErrorCorrection.DataSet;
import ErrorCorrection.Read;
import ErrorCorrection.ReadFreqComparator;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;

/**
 *
 * @author pavel
 */
public class test1 {
    public static void main(String[] args) throws IOException
    {
        int thr = 1;
        String type = "power";
        String folder = "test_res";
        String trueSeq = type + ".fas";
        String algRes = folder + File.separator + type + "-reads.fas_corrected.fas";
        DataSet ds = new DataSet(algRes);
        ds.delGaps();
        ds.PrintReads(algRes + "_noGaps.fas");
        ds = new DataSet(algRes + "_noGaps.fas");
        ReadFreqComparator rfc = new ReadFreqComparator(); 
        Collections.sort(ds.reads, rfc);
        ds.PrintUniqueReadsThr(algRes + "_unique.fas",thr);
        ds = new DataSet(algRes + "_unique.fas","ET");
        System.out.println("All reads: " + ds.getTotalNReads());
        System.out.println("Unique reads: " + ds.reads.size());
        DataSet dsTrue = new DataSet(trueSeq);
        int nTrueFound = 0;
        for (Read r : dsTrue.reads)
            if (ds.containRead(r))
                nTrueFound++;
        System.out.println("nTrueFound: " + nTrueFound);
        int freqTrueFound = 0;
        for (Read r : dsTrue.reads)
            freqTrueFound += ds.getFrequency(r);
        System.out.println("freqTrueFound: " + freqTrueFound);
            
        
        double maxgap = 0;
        int igap = 0;
        for (int i = 1; i < ds.reads.size(); i++)
        {
            double d = ((double) ds.reads.get(i-1).getFreq()) / ds.reads.get(i).getFreq();
            if (d > maxgap)
            {
                maxgap = d;
                igap = i-1;
            }
        }
        System.out.println(maxgap + " " + igap);
        
    }
    
}
