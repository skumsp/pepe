/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepe;

import ErrorCorrection.DataSet;
import ErrorCorrection.Read;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.StringTokenizer;

/**
 *
 * @author pavel
 */
public class Simulator 
{
    void generateReads(String refFile, double errorRate, int nReads, String outDir) throws IOException
    {
        DataSet ds = new DataSet(refFile,"ET");
        int n = ds.reads.size();
        String freqFileName = "freq.txt";
        int len = ds.reads.get(0).getLength();
        
        int totCount = 0;
        for (int i = 0; i < n; i++)
            totCount+= ds.reads.get(i).getFreq();
        
        FileWriter fw = new FileWriter(freqFileName);
        for (Read r : ds.reads)
        {
            double freq = 100*((double) r.getFreq())/totCount;
            fw.write(r.name + " " + freq + "\n");            
        }
        fw.close();
        
        StringTokenizer st = new StringTokenizer(refFile,"_" + File.separator);
        st.nextToken();
        String prefix = st.nextToken();
        
        System.out.println("Running Grinder \n");
        Runtime run=Runtime.getRuntime();
        Process p=null;
        String s = "grinder -reference_file " + refFile + " -tr " + nReads + " -rd " + len + " -id " + len;
        s +=  " -md uniform " + errorRate + " -mr 100 0 " + "-af " + freqFileName + " -bn " + prefix + " -od " + outDir;

        p = run.exec(s);
        BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));

        System.out.println("Grinder output:\n");
        while ((s = stdInput.readLine()) != null) 
            System.out.println(s);
//        if (p.exitValue()!= 0)
//            System.out.println("Error in Grinder");
        System.out.println("Finished!");
                
    }
    
}
