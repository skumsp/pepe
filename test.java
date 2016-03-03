/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepe;

import ErrorCorrection.DataSet;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;

/**
 *
 * @author kki8
 */
public class test {
    public static void main(String[] args) throws IOException, CompoundNotFoundException
    {
        String folder_name = "test2";
        String refAddr = "refs.fas";
        DataSet refs = new DataSet(refAddr,'c');
        File folder = new File(folder_name);
                
        File[] list_files = folder.listFiles();
        
        for (int i = 0; i < list_files.length; i++)
        {
            String dset_file = list_files[i].getPath();
            DataSet ds = new DataSet(dset_file);
            ds.clipToRef(refs, 15, 6);
            ds.findHaplotypes("KEC","Frequency");
            ds.PrintHaplotypes(dset_file + "_hapl_clipped.fas");
            
        }
        
    }
    
}
