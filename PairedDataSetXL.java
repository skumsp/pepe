/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepe;

import ErrorCorrection.DataSet;
import ErrorCorrection.DataSetXL;
import ErrorCorrection.Read;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 *
 * @author kki8
 */
public class PairedDataSetXL {
    DataSetXL forward;
    DataSetXL reverse;
    DataSetXL joint;
    HashMap<Integer,Integer> mapForwJ;
    HashMap<Integer,Integer> mapRevJ;
    
    public PairedDataSetXL(DataSetXL ds1, DataSetXL ds2)
    {
        forward = ds1;
        reverse = ds2;
    }
    
    public void comparePairedEndPrintStat(String addr, int gapop, int gapext) throws InterruptedException, ExecutionException, IOException
    {
        int n = forward.reads.size();
        ArrayList<Read> badForw = new ArrayList();
        ArrayList<Read> badRev = new ArrayList();
        ArrayList<ArrayList<Integer>> allErr = new ArrayList();
        List< Future > futuresList = new ArrayList< Future >();
        int nrOfProcessors = Runtime.getRuntime().availableProcessors();
//        int nrOfProcessors = 1;
        ExecutorService eservice = Executors.newFixedThreadPool(nrOfProcessors);
        
        Iterator irForw = this.forward.reads.iterator();
        Iterator irRev = this.reverse.reads.iterator();
        while (irForw.hasNext())
        {
            Read forr = (Read) irForw.next();
            Read revr = (Read) irRev.next();
            processPairedEndReadsTask rt = new processPairedEndReadsTask(forr,revr,gapop,gapext);
            futuresList.add(eservice.submit(rt));
			 
        }
        irForw = this.forward.reads.iterator();
        irRev = this.reverse.reads.iterator();
        for(int i = 0; i < futuresList.size(); i++) 
        {
            Read forr = (Read) irForw.next();
            Read revr = (Read) irRev.next();
            Future future = futuresList.get(i);
            ArrayList<Integer> ar = (ArrayList<Integer>) future.get();
            System.out.println("Local alignment: read "+i + "/" + n);
            if (ar.size() == 1)
            {
                badForw.add(forr);
                badRev.add(revr);
            }
            else
                allErr.add(ar);
        }
        forward.reads.removeAll(badForw);
        reverse.reads.removeAll(badRev);
        FileWriter fw = new FileWriter(addr);
        for (int i = 0; i < allErr.size(); i++)
        // indel subs totalerr length
            fw.write((i+1) + " " + allErr.get(i).get(0) + " " + allErr.get(i).get(1) + " " + 
                    (allErr.get(i).get(0)+allErr.get(i).get(1)) + " " + allErr.get(i).get(2) + "\n");
        fw.close();
        
    }
    public void createJointDS()
    {
        joint = new DataSetXL();
        HashMap<String,Integer> seq = new HashMap<String, Integer>();
        Iterator irForw = this.forward.reads.iterator();
        Iterator irRev = this.reverse.reads.iterator();
        int i = 0;
        while (irForw.hasNext())
        {
            System.out.println("Joint DS: " + (i++) + "\\" + this.getNPairedReads());
            Read rForw = (Read) irForw.next();
            Read rRev = (Read) irRev.next();
            String s = rForw.getSeqNoGaps();
            if (seq.containsKey(s))
            {
                int f = seq.get(s);
                seq.put(s, f + rForw.getFreq());
            }
            else
            {
                seq.put(s, rForw.getFreq());
            }
            
            s = rRev.getSeqNoGaps();
            if (seq.containsKey(s))
            {
                int f = seq.get(s);
                seq.put(s, f + rRev.getFreq());
            }
            else
            {
                seq.put(s, rRev.getFreq());
            }
        }
        int nUnique = seq.size();
        i = 0;
        Iterator ir = seq.entrySet().iterator();
        while (ir.hasNext())
        {
            System.out.println("Create joint DS reads: " + i + "\\" + nUnique);
            Map.Entry me = (Map.Entry) ir.next();
            String nucl = (String) me.getKey();
            String name = "read" + i;
            int f = (Integer) me.getValue();
            joint.reads.add(new Read(nucl,name,f));
            i++;
        }
        seq = null;
        
    }
    public void delJointDS()
    {
        joint.reads = null;
        joint = null;
        System.gc();
    }
    int getNPairedReads()
    {
        return forward.reads.size();
    }
    int correctErrors(int k) throws InterruptedException, ExecutionException
    {
        joint.setK(k);
        joint.calculateKMersAndKCounts();
        joint.reads = null;
        
        List< Future > futuresList = new ArrayList< Future >();
//        int nrOfProcessors = Runtime.getRuntime().availableProcessors();
        int nrOfProcessors = 32;
        ExecutorService eservice = Executors.newFixedThreadPool(nrOfProcessors);
            
        Iterator irForw = this.forward.reads.iterator();
        Iterator irRev = this.reverse.reads.iterator();
        
        while (irForw.hasNext())
        {
            Read rFor = (Read) irForw.next();
            Read rRev = (Read) irRev.next();
            ErrorCorrectionTaskXL et = new ErrorCorrectionTaskXL(rFor,rRev,joint.kmers,k);
            futuresList.add(eservice.submit(et));
        }
        
        int i = 0;
        int nErrCorr = 0;
        for(Future future:futuresList) 
        {
            int ncorr = (Integer) future.get();
            nErrCorr+=ncorr;
            System.out.println("Error correction: read "+i + "/" + this.getNPairedReads());
            i++;
        }
        
        System.out.println();
        return nErrCorr;
    }
    void printReadsOneFile(String addr) throws IOException
    {
        FileWriter fw = new FileWriter(addr);
        Iterator irForw = this.forward.reads.iterator();
        Iterator irRev = this.reverse.reads.iterator();
        while (irForw.hasNext())
        {
            Read rForw = (Read) irForw.next();
            Read rRev = (Read) irRev.next();
            fw.write(">" + rForw.name + "\n" + rForw.nucl + "\n");
            fw.write(">" + rRev.name + "\n" + rRev.nucl + "\n");
        }
        fw.close();
    }
    void delShortReads(int cutoff)
   {
        Iterator irForw = this.forward.reads.iterator();
        Iterator irRev = this.reverse.reads.iterator();
        int i = 1;
        while (irForw.hasNext())
        {
            System.out.println("Del short: read " + (i++) + "//" + this.getNPairedReads());
            Read rForw = (Read) irForw.next();
            Read rRev = (Read) irRev.next();
            if ((rForw.getLength() < cutoff) || (rRev.getLength() < cutoff))
            {
                irForw.remove();
                irRev.remove();
            }
        }
    }
    DataSet getPerfectReads()
    {
        DataSet ds = new DataSet();
        Iterator irForw = this.forward.reads.iterator();
        Iterator irRev = this.reverse.reads.iterator();
        int i = 1;
        while (irForw.hasNext())
        {
            Read rForw = (Read) irForw.next();
            Read rRev = (Read) irRev.next();
            System.out.println("Getting perfect: " + (i++) + "//" + this.getNPairedReads());
            if (rForw.nucl.equalsIgnoreCase(rRev.nucl))
                ds.reads.add(rForw);
        }
        return ds;
    }
    void delAllNs()
    {

        Iterator irForw = this.forward.reads.iterator();
        Iterator irRev = this.reverse.reads.iterator();
        int i = 1;
        while (irForw.hasNext())
        {
            Read rForw = (Read) irForw.next();
            Read rRev = (Read) irRev.next();
            System.out.println("Removing Ns: read " + (i++) + "//" + this.getNPairedReads());
            int nN = 0;
            for (int j = 0; j < rForw.nucl.length(); j++)
                if (rForw.nucl.charAt(j) == 'N')
                    nN++;
            if (nN == rForw.nucl.length())
            {
                irForw.remove();
                irRev.remove();
                continue;
            }
            nN = 0;
            for (int j = 0; j <rRev.nucl.length(); j++)
                if (rRev.nucl.charAt(j) == 'N')
                    nN++;
            if (nN == rRev.nucl.length())
            {
                irForw.remove();
                irRev.remove();
                continue;
            }
        }
    }
}
