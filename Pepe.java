/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepe;

import ErrorCorrection.DataSet;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;

/**
 *
 * @author kki8
 */
public class Pepe {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException 
    {
        int gapop = 15;
        int gapext =6;
        double percBad = 0.5;
        int k = 25;
        int minlen = 180;
        
        String folder_name = "NY_fas_all";
        File folder = new File(folder_name);
                
        File[] list_files = folder.listFiles();
        
        String refFile_name = "HCV_HVR1_264_single_subgen.fas";
        DataSet refs = new DataSet(refFile_name,'c');
        
        FileWriter fw = new FileWriter("stats.txt");
                
        for (int i = 0; i < list_files.length; i+=2)
        {
            String dset_file1 = list_files[i].getPath();
            DataSet ds1 = new DataSet(dset_file1,'c');
            String dset_file2 = list_files[i+1].getPath();
            DataSet ds2 = new DataSet(dset_file2,'c');
            
            long startTime = System.currentTimeMillis();
            
            PairedDataSet ds = new PairedDataSet(ds1,ds2);
            fw.write(list_files[i].getName() + " ");
            fw.write(ds.getNPairedReads() + " ");
            
/*            ds.removeBadReads(percBad);
            
            ds.forward.fixDirectionGenotypingRefParallel(refs, gapop, gapext);
            ds.reverse.fixDirectionGenotypingRefParallel(refs, gapop, gapext);
            ds.forward.PrintUniqueReadsNoFreqTag(dset_file1 + "_reversed.fas");
            ds.reverse.PrintUniqueReadsNoFreqTag(dset_file2 + "_reversed.fas");*/
            
            ds.delAllNs();
//            ds.printReadsOneFile(dset_file1 + "_together.fas");
            ds.comparePairedEndPrintStat(dset_file1 + "_stats.txt", gapop, gapext);
            ds.delShortReads(minlen);
            DataSet perf = ds.getPerfectReads();
            fw.write(perf.reads.size()+ " ");
//            ds.printReadsOneFile(dset_file1 + "_clipped_together.fas");
            ds.createJointDS();
            //            ds.joint.PrintReadsStat(folder_name);            
            ds.correctErrors(k);
            perf = ds.getPerfectReads();
            fw.write(perf.reads.size()+ " ");
            long tm = System.currentTimeMillis() - startTime;
            fw.write("" + ((double) tm) / 1000);
            
//            ds.comparePairedEndPrintStat(dset_file1 + "_stats_corrected.txt", gapop, gapext);
//            ds.printReadsOneFile(dset_file1 + "_corrected_together.fas");
            fw.write("\n");
            
        }
        fw.close();
        System.exit(0);
    }
    
}
