/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepe;

import ErrorCorrection.DataSet;
import ErrorCorrection.Read;
import ErrorCorrection.RevCompGenTask;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
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
public class PairedDataSet {
    DataSet forward;
    DataSet reverse;
    DataSet joint;
    HashMap<Integer,Integer> mapForwJ;
    HashMap<Integer,Integer> mapRevJ;
    
    public PairedDataSet(DataSet ds1, DataSet ds2)
    {
        forward = ds1;
        reverse = ds2;
    }
    public PairedDataSet(String addr,String type) throws FileNotFoundException, IOException
    {
        BufferedReader br = new BufferedReader(new FileReader(addr));
        forward = new DataSet();
        reverse = new DataSet();
        if (type.equalsIgnoreCase("grinder"))
        {
            int i = 0;
            int j = 0;
            String name = br.readLine();
            while (name != null)
            {
                if (i==0)
                    j++;


                String s = br.readLine();
                String seq = "";
                while ((s != null) && (!s.startsWith(">")))
                {
                    seq+=s;
                    s = br.readLine();
                }
                if (i==0)
                    System.out.println("Reading read " + j + " forward");
                else
                    System.out.println("Reading read " + j + " reverse");
            
                Read r = new Read(seq,name);
                if (r.name.contains("complement"))
                    r.RevComp();
                if (i==0)
                {
                    forward.reads.add(r);
                    i=1;
                }
                else
                {
                    reverse.reads.add(r);
                    i=0;
                }
                name = s;  

            }
        }
        br.close();
    }
    public void removeBadReads(double perc)
    {
        ArrayList<Read> removeForw = new ArrayList();
        ArrayList<Read> removeRev = new ArrayList();
        for (int i = 0; i < forward.reads.size(); i++)
        {
            Read r = forward.reads.get(i);
            int nBad = 0;
            for (int j = 0; j < r.getLength(); j++)
                if (r.getNucl(j) == 'N')
                    nBad++;
            double percN = ((double) nBad)/r.getLength();
            if (percN >= perc)
            {
                removeForw.add(forward.reads.get(i));
                removeRev.add(reverse.reads.get(i));
            }
        }
        for (int i = 0; i < reverse.reads.size(); i++)
        {
            Read r = reverse.reads.get(i);
            int nBad = 0;
            for (int j = 0; j < r.getLength(); j++)
                if (r.getNucl(j) == 'N')
                    nBad++;
            double percN = ((double) nBad)/r.getLength();
            if (percN >= perc)
            {
                removeForw.add(forward.reads.get(i));
                removeRev.add(reverse.reads.get(i));
            }
                          
        }
        forward.reads.removeAll(removeForw);
        reverse.reads.removeAll(removeRev);

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
                 
        for (int i = 0; i < n; i++)
        {
            processPairedEndReadsTask rt = new processPairedEndReadsTask(forward.reads.get(i),reverse.reads.get(i),gapop,gapext);
            futuresList.add(eservice.submit(rt));
			 
        }
        for(int i = 0; i < futuresList.size(); i++) 
        {
            Future future = futuresList.get(i);
            ArrayList<Integer> ar = (ArrayList<Integer>) future.get();
            System.out.println("Local alignment: read "+i + "/" + n);
            if (ar.size() == 1)
            {
                badForw.add(forward.reads.get(i));
                badRev.add(reverse.reads.get(i));
            }
            else
                allErr.add(ar);
        }
        forward.reads.removeAll(badForw);
        reverse.reads.removeAll(badRev);
        FileWriter fw = new FileWriter(addr);
        for (int i = 0; i < allErr.size(); i++)
            fw.write((i+1) + " " + forward.reads.get(i).name + " " + allErr.get(i).get(0) + " " + allErr.get(i).get(1) + " " + 
                    (allErr.get(i).get(0)+allErr.get(i).get(1)) + " " + allErr.get(i).get(2) + "\n");
        fw.close();
        
    }
    public void createJointDS()
    {
        joint = new DataSet();
        mapForwJ = new HashMap<Integer,Integer>();
        mapRevJ = new HashMap<Integer,Integer>();
        HashMap<String,ArrayList<String>> seq = new HashMap<String, ArrayList<String>>();
        for (int i = 0; i < this.getNPairedReads(); i++)
        {
            System.out.println("Joint DS: " + i + "\\" + this.getNPairedReads());
            String s = forward.reads.get(i).getSeqNoGaps();
            if (seq.containsKey(s))
            {
                ArrayList<String> ar = seq.get(s);
                ar.add(("f"+i));
                seq.put(s, ar);
            }
            else
            {
                ArrayList<String> ar = new ArrayList();
                ar.add(("f"+i));
                seq.put(s, ar);
            }
            
            s = reverse.reads.get(i).getSeqNoGaps();
            if (seq.containsKey(s))
            {
                ArrayList<String> ar = seq.get(s);
                ar.add(("r"+i));
                seq.put(s, ar);
            }
            else
            {
                ArrayList<String> ar = new ArrayList();
                ar.add(("r"+i));
                seq.put(s, ar);
            }
        }
        int nUnique = seq.size();
        int i = 0;
        Iterator ir = seq.entrySet().iterator();
        while (ir.hasNext())
        {
            System.out.println("Create joint DS reads: " + i + "\\" + nUnique);
            Map.Entry me = (Map.Entry) ir.next();
            String nucl = (String) me.getKey();
            String name = "read" + i;
            ArrayList<String> ar = (ArrayList<String>) me.getValue();
            joint.reads.add(new Read(nucl,name,ar.size()));
            for (String s : ar)
            {
                int j = Integer.parseInt(s.substring(1));
                if (s.charAt(0) == 'f')
                    this.mapForwJ.put(j, i);
                else
                    this.mapRevJ.put(j, i);
            }
            i++;
        }
        seq = null;
        
    }
    public void delJointDS()
    {
        for (Read r : joint.reads)
            r.kmers = null;
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
        
        List< Future > futuresList = new ArrayList< Future >();
//        int nrOfProcessors = Runtime.getRuntime().availableProcessors();
        int nrOfProcessors = 32;
        ExecutorService eservice = Executors.newFixedThreadPool(nrOfProcessors);
                         
        for (int i = 0; i < this.getNPairedReads(); i++)
        {
            Read rFor = this.forward.reads.get(i);
            Read rRev = this.reverse.reads.get(i);
            Read rForJoint = joint.reads.get(this.mapForwJ.get(i));
            Read rRevJoint = joint.reads.get(this.mapRevJ.get(i));
            ErrorCorrectionTask et = new ErrorCorrectionTask(rFor,rRev,rForJoint,rRevJoint,k);
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
        for (int i = 0; i < this.getNPairedReads(); i++)
        {
            fw.write(">" + forward.reads.get(i).name + "\n" + forward.reads.get(i).nucl + "\n");
            fw.write(">" + reverse.reads.get(i).name + "\n" + reverse.reads.get(i).nucl + "\n");
        }
        fw.close();
    }
    void delShortReads(int cutoff)
    {
        ArrayList<Read> toDelFor = new ArrayList();
        ArrayList<Read> toDelRev = new ArrayList();
        for (int i = 0; i < this.getNPairedReads(); i++)
            if ((this.forward.reads.get(i).getLength() < cutoff) || (this.reverse.reads.get(i).getLength() < cutoff))
            {
                toDelFor.add(this.forward.reads.get(i));
                toDelRev.add(this.reverse.reads.get(i));
            }
        this.forward.reads.removeAll(toDelFor);
        this.reverse.reads.removeAll(toDelRev);
    }
    DataSet getPerfectReads()
    {
        DataSet ds = new DataSet();
        for (int i = 0; i < this.getNPairedReads(); i++)
        {
            System.out.println("Getting perfect: " + i + "//" + this.getNPairedReads());
            if (this.forward.reads.get(i).nucl.equalsIgnoreCase(this.reverse.reads.get(i).nucl))
                ds.reads.add(this.forward.reads.get(i));
        }
        return ds;
    }
    void delAllNs()
    {
        ArrayList<Read> toDelFor = new ArrayList();
        ArrayList<Read> toDelRev = new ArrayList();
        for (int i = 0; i < this.getNPairedReads(); i++)
        {
            System.out.println("Removing Ns: read " + i + "//" + this.getNPairedReads());
            int nN = 0;
            for (int j = 0; j < this.forward.reads.get(i).nucl.length(); j++)
                if (this.forward.reads.get(i).nucl.charAt(j) == 'N')
                    nN++;
            if (nN == this.forward.reads.get(i).nucl.length())
            {
                toDelFor.add(this.forward.reads.get(i));
                toDelRev.add(this.reverse.reads.get(i));
            }
            nN = 0;
            for (int j = 0; j < this.reverse.reads.get(i).nucl.length(); j++)
                if (this.reverse.reads.get(i).nucl.charAt(j) == 'N')
                    nN++;
            if (nN == this.reverse.reads.get(i).nucl.length())
            {
                toDelFor.add(this.forward.reads.get(i));
                toDelRev.add(this.reverse.reads.get(i));
            }
        }
        this.forward.reads.removeAll(toDelFor);
        this.reverse.reads.removeAll(toDelRev);
    }
}
