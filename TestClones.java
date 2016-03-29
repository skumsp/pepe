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
import java.util.Map;
import java.util.concurrent.ExecutionException;

/**
 *
 * @author kki8
 */
public class TestClones {
    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException 
    {
		int gapop = 15;
                int gapext = 6; 
                int chimThr = 30;
                int refThr = 3;
                int minlen = 250;
                
                String folder_name  = "Filtered_by_tag_both_MID_primer_casper_05_fas";
                File fl_ref = new File("clones_cut.fas");                
                String refFile_name = fl_ref.getPath();
                
                                
                DataSet refs = new DataSet(refFile_name,'c');
                
                File folder = new File(folder_name);
                File[] list_files = folder.listFiles();
                
                for (int i = 0; i < list_files.length; i++)
                {
                    String dset_file = list_files[i].getPath();
                    DataSet ds = new DataSet(dset_file);
                    
                    ds.delShortReads(minlen);
                    Collections.sort(ds.reads, new ReadFreqComparator());
                    ds.PrintUniqueReads(dset_file + "_unique.fas");
                    
                    ds.setAvProc(Runtime.getRuntime().availableProcessors());
                    ds.fixDirectionGenotypingRefParallel(refs, gapop, gapext);
                    
                    HashMap<String,Integer> totFreqRef0 = new HashMap();
                    HashMap<String,Integer> maxFreqRef0 = new HashMap();
                    HashMap<String,Integer> totFreqRef3 = new HashMap();
                    HashMap<String,Integer> maxFreqRef3 = new HashMap();
                    HashSet<String> allRefs = new HashSet();
                    int freqchim = 0;
                    int maxchim = 0;
                    
                    FileWriter fw = new FileWriter(dset_file + "_stats_reads.txt");
                    for (Read r : ds.reads)
                    {
                        allRefs.add(r.genotype);
                        fw.write(r.name + " " + r.getFreq() + " " + r.getLength() + " " + r.genotype + " " + r.distgen + "\n");
                        
                        if (r.distgen >= chimThr)
                        {
                            freqchim = freqchim + r.getFreq();
                            if (r.getFreq() > maxchim)
                                maxchim = r.getFreq();
                        }
                        if (r.distgen == 0)
                        {
                            if (totFreqRef0.containsKey(r.genotype))
                            {
                                int f = totFreqRef0.get(r.genotype);
                                totFreqRef0.put(r.genotype, f + r.getFreq());
                                int m = maxFreqRef0.get(r.genotype);
                                if (r.getFreq() > m)
                                    maxFreqRef0.put(r.genotype, r.getFreq());
                            }
                            else
                            {
                                totFreqRef0.put(r.genotype, r.getFreq());
                                maxFreqRef0.put(r.genotype, r.getFreq());
                            }
                        }
                        if (r.distgen <= refThr)
                        {
                            if (totFreqRef3.containsKey(r.genotype))
                            {
                                int f = totFreqRef3.get(r.genotype);
                                totFreqRef3.put(r.genotype, f + r.getFreq());
                                int m = maxFreqRef3.get(r.genotype);
                                if (r.getFreq() > m)
                                    maxFreqRef3.put(r.genotype, r.getFreq());
                            }
                            else
                            {
                                totFreqRef3.put(r.genotype, r.getFreq());
                                maxFreqRef3.put(r.genotype, r.getFreq());
                            }
                        }
                        
                        if (r.distgen > refThr)
                        {
                            if (!totFreqRef0.containsKey(r.genotype))
                            {
                                totFreqRef0.put(r.genotype, 0);
                                maxFreqRef0.put(r.genotype, 0);
                            }
                            if (!totFreqRef3.containsKey(r.genotype))
                            {
                                totFreqRef3.put(r.genotype, 0);
                                maxFreqRef3.put(r.genotype, 0);
                            }
                        }
                    }
                    fw.close();
                                        
                    int totalReads = ds.getNreads();
                    fw = new FileWriter(dset_file + "_stats_total.txt");
                    
                    double fc = 100*((double) freqchim)/totalReads;
                    double mc = 100*((double) maxchim)/totalReads;
                    fw.write("chim" + " " + freqchim + " " + maxchim + " " + fc + " " + mc + " " + "\n");
                    
                    for (String s : allRefs)
                    {
                        int f0 = 0;
                        int m0 = 0;
                        if (totFreqRef0.containsKey(s))
                        {
                            f0 = totFreqRef0.get(s);
                            m0 = maxFreqRef0.get(s);
                        }
                        double ff0 = 100*((double) f0)/totalReads;
                        double mf0 = 100*((double) m0)/totalReads;
                        
                        
                        int f3 = 0;
                        int m3 = 0;
                        if (totFreqRef3.containsKey(s))
                        {
                            f3 = totFreqRef3.get(s);
                            m3 = maxFreqRef3.get(s);
                        }
                        double ff3 = 100*((double) f3)/totalReads;
                        double mf3 = 100*((double) m3)/totalReads;
                        
                        if (f3 > 0)
                        {
                            DecimalFormat formatter = new DecimalFormat("#0.000000");                        
                            fw.write(s + " " + f0 + " " + m0 + " " + formatter.format(ff0) + " " + formatter.format(mf0) + " " + "\n");
                            fw.write(s + " " + f3 + " " + m3 + " " + formatter.format(ff3) + " " + formatter.format(mf3) + " " + "\n");
                        }
                    }
                    
                    fw.close();
                }
                System.exit(0);
    }
}
