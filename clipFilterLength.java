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
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutionException;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;

/**
 *
 * @author kki8
 */
public class clipFilterLength {
    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException, CompoundNotFoundException 
    {
		int gapop = 15;
                int gapext = 6; 
                int chimThr = 30;
                int refThr = 3;
                int minlen = 250;
                
                String folder_name  = "bootstrapping_pepe";
                File fl_ref = new File("clones_cut.fas");                
                String refFile_name = fl_ref.getPath();
                
                                
                DataSet refs = new DataSet(refFile_name,'c');
                
                File folder = new File(folder_name);
                File[] list_files = folder.listFiles();
                
                for (int i = 0; i < list_files.length; i++)
                {
                    String dset_file = list_files[i].getPath();
                    DataSet ds = new DataSet(dset_file,"ET");
                    
                    ds.setAvProc(Runtime.getRuntime().availableProcessors());
                    ds.clipToRefParallel(refs, gapop, gapext);
                    Collections.sort(ds.reads, new ReadFreqComparator());
                    ds.PrintUniqueReads(dset_file + "_clipped.fas");
                }
                System.exit(0);
    }
}
