/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepe;

import java.io.File;
import java.io.IOException;

/**
 *
 * @author pavel
 */
public class GenerateReads {
    public static void main(String[] args) throws IOException
    {
        double errorRate = 0.12;
        int nReads = 100000;
        
        double prob = 1;
        
        String refsFold = "popRefs1";
        String outFold = refsFold + "_simul";
        File folder = new File(refsFold);
        
        File[] listFiles = folder.listFiles();
        Simulator sm = new Simulator();
        for (int i = 0; i < listFiles.length; i++)
        {
            double p = Math.random();
            if (p > prob)
                continue;
            String refName = listFiles[i].getPath();
            System.out.println(refName);
            sm.generateReads(refName, errorRate, nReads, outFold);
        }
    }
    
}
