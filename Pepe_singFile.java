/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepe;

import ErrorCorrection.DataSet;
import ErrorCorrection.Read;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;
import java.util.concurrent.ExecutionException;

/**
 *
 * @author pavel
 */
public class Pepe_singFile 
{
    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException 
    {
        int gapop = 15;
        int gapext =6;
        double percBad = 0.5;
        int k = 25;
        int minlen = 306;
        int maxIter = 5;
        
        String folder_name = "popRefs_simul";
        File folder = new File(folder_name);
                
        File[] list_files = folder.listFiles();
       
        
        FileWriter fw = new FileWriter("stats.txt");
                
        for (int i = 0; i < list_files.length; i++)
        {
            String dset_file = list_files[i].getPath();
            
            long startTime = System.currentTimeMillis();
            
            PairedDataSet ds = new PairedDataSet(dset_file,"grinder");
            
            
            fw.write(list_files[i].getName() + " ");
            fw.write(ds.getNPairedReads() + " ");
            
            ds.comparePairedEndPrintStat(dset_file + "_stats.txt", gapop, gapext);
            ds.delShortReads(minlen);
            DataSet perf = ds.getPerfectReads();
            perf.PrintUniqueReadsNoFreqTag(dset_file + "_precorrected.fas");
            fw.write(perf.reads.size()+ " ");
//            ds.printReadsOneFile(dset_file1 + "_clipped_together.fas");
            ds.createJointDS();
            //            ds.joint.PrintReadsStat(folder_name);            
            int nErrCorr = ds.correctErrors(k);
            int nIter = 1;
            while ((nErrCorr > 0) && (nIter <= maxIter))
            {
                ds.delJointDS();
                ds.createJointDS();
                nErrCorr = ds.correctErrors(k);
                nIter++;
            }
            perf = ds.getPerfectReads();
            perf.PrintUniqueReadsNoFreqTag(dset_file + "_corrected.fas");
            fw.write(perf.reads.size()+ " ");
            long tm = System.currentTimeMillis() - startTime;
            fw.write("" + ((double) tm) / 1000);
            
            fw.write("\n");
            
/*            DataSet check = new DataSet("one-reads.fas_corrected.fas_unique.fas");
            for (Read r : perf.reads)
                if (r.nucl.equalsIgnoreCase(check.reads.get(0).nucl))
                    System.out.println(r.name);*/
            
        }
        fw.close();
        System.exit(0);
    }
    
}
