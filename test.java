/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepe;

import ErrorCorrection.DataSet;
import ErrorCorrection.ReadFreqComparator;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;

/**
 *
 * @author pavel
 */
public class test {
    public static void main(String[] args) throws IOException
    {
        
        String folder_name = "test_res";
        File folder = new File(folder_name);
                
        File[] list_files = folder.listFiles();
        int thr = 1;

                
        for (int i = 0; i < list_files.length; i++)
        {
            String dset_file = list_files[i].getPath();
            DataSet ds = new DataSet(dset_file);
            DataSet ds1 = new DataSet(dset_file,'c');
            ReadFreqComparator rfc = new ReadFreqComparator(); 
            Collections.sort(ds.reads, rfc);
            ds.PrintUniqueReadsThr(dset_file+ "_unique.fas",thr);
        }
    }
    
}
