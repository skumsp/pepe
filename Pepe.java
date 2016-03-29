/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepe;

import ErrorCorrection.DataSet;
import ErrorCorrection.DataSetXL;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;
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
        int minlen = 220;
        int maxIter = 3;
        
        String folder_name = "Illumina_Raw_demultiplexed_tag_MID_fas";
        File folder = new File(folder_name);
                
        File[] list_files = folder.listFiles();
        
        String refFile_name = "HCV_HVR1_264.fas";
        DataSet refs = new DataSet(refFile_name,'c');
        
        
        ArrayList<String> forwFiles = new ArrayList();
        ArrayList<String> revFiles = new ArrayList();
        
        for (int i = 0; i < list_files.length; i++)
        {
            String dset_file1 = list_files[i].getName();
            if (dset_file1.contains("R1"))
            {
                StringTokenizer st = new StringTokenizer(dset_file1,"_.");
                String ds_name = "";
                String s = st.nextToken();
                while (!s.equalsIgnoreCase("R1"))
                {
                    ds_name = ds_name + s + "_";
                    s = st.nextToken();
                }
                for (int j = 0; j < list_files.length; j++)
                {
                    String dset_file2 = list_files[j].getName();
                    if (dset_file2.startsWith(ds_name) && !dset_file2.equalsIgnoreCase(dset_file1))
                    {
                        forwFiles.add(list_files[i].getPath());
                        revFiles.add(list_files[j].getPath());
                    }
                        
                }
            }
            
            
        }
        
        FileWriter fw = new FileWriter("stats.txt");
                
    for (int i = 0; i < forwFiles.size(); i++)
        {
            String dset_file1 = forwFiles.get(i);
            DataSetXL ds1 = new DataSetXL(dset_file1,'c');
            String dset_file2 = revFiles.get(i);
            DataSetXL ds2 = new DataSetXL(dset_file2,'c');
            
            long startTime = System.currentTimeMillis();
            
            PairedDataSetXL pds = new PairedDataSetXL(ds1,ds2);
            fw.write(dset_file1 + " ");
            fw.write(pds.getNPairedReads() + " ");
            
            pds.delAllNs();
//            ds.printReadsOneFile(dset_file1 + "_together.fas");
            pds.comparePairedEndPrintStat(dset_file1 + "_stats.txt", gapop, gapext);
            pds.delShortReads(minlen);
            DataSet perf = pds.getPerfectReads();
            fw.write(perf.reads.size()+ " ");
//            ds.printReadsOneFile(dset_file1 + "_clipped_together.fas");
            pds.createJointDS();
            //            ds.joint.PrintReadsStat(folder_name);            
            int nErrCorr = pds.correctErrors(k);
            int nIter = 1;
            while ((nErrCorr > 0) && (nIter <= maxIter))
            {
                pds.delJointDS();
                pds.createJointDS();
                nErrCorr = pds.correctErrors(k);
                nIter++;
            }
            perf = pds.getPerfectReads();
            perf.PrintUniqueReadsNoFreqTag(dset_file1 + "_perfect.fas");
            fw.write(perf.reads.size()+ " ");
            long tm = System.currentTimeMillis() - startTime;
            fw.write("" + ((double) tm) / 1000);
            
            pds.comparePairedEndPrintStat(dset_file1 + "_stats_corrected.txt", gapop, gapext);            
            pds.printReadsOneFile(dset_file1 + "_corrected_together.fas"); 
            pds.forward.PrintReads(dset_file1 + "_corrected_forward.fas");
            fw.write("\n");
            
        }
        fw.close();
        System.exit(0);
    }
    
}
