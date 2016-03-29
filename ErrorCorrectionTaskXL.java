/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepe;

import ErrorCorrection.Kmer;
import ErrorCorrection.Read;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.Callable;

/**
 *
 * @author kki8
 */
public class ErrorCorrectionTaskXL implements Callable
{
    Read rFor;
    Read rRev;
    HashMap<String,Integer> kmers;
    int k;

    ErrorCorrectionTaskXL(Read r1, Read r2, HashMap hm, int i)
    {
        rFor = r1;
        rRev = r2;
        kmers = hm;
        k = i;
    }
    public Object call()
    {
//        if (rFor.name.equalsIgnoreCase(">40410/1 reference=one_1_1 position=1..306 errors=292%A"))
//                System.out.println();
        Read rForJoint = new Read(rFor.getSeqNoGaps());
        Read rRevJoint = new Read(rRev.getSeqNoGaps()); 
        ArrayList<Integer> errors = new ArrayList();
        ArrayList<Integer> errorsForJoint = new ArrayList();
        ArrayList<Integer> errorsRevJoint = new ArrayList();
        errors.add(-1);
        errorsForJoint.add(-1);
        errorsRevJoint.add(-1);
        int iForJoint = 0;
        int iRevJoint = 0;
        for (int j = 0; j < rFor.getLength(); j++)
        {
            if (rFor.nucl.charAt(j)!=rRev.nucl.charAt(j))
            {
                errorsForJoint.add(iForJoint);
                errorsRevJoint.add(iRevJoint);
                errors.add(j);
            }
            if (rFor.nucl.charAt(j) != '-')
                iForJoint++;
            if (rRev.nucl.charAt(j) != '-')
                iRevJoint++;
        }
        errors.add(rFor.getLength()+1);
        errorsForJoint.add(rForJoint.getLength());
        errorsRevJoint.add(rRevJoint.getLength());
        
        if (errors.size() == 2)
            return 0;
        
        
        StringBuilder sbFor = new StringBuilder(rFor.nucl);
        StringBuilder sbRev = new StringBuilder(rRev.nucl);
        
        int nErrCorr = 0;
        for (int i = 1; i < errors.size() -1; i++)
        {
            try
            {
                int bFor;
                int eFor;
                if (rFor.nucl.charAt(errors.get(i)) != '-')
                {
                    bFor = Math.max(errorsForJoint.get(i)-k+1, errorsForJoint.get(i-1)+1);
                    eFor = Math.min(errorsForJoint.get(i), errorsForJoint.get(i+1)-k);
                }
                else
                {
                    bFor = Math.max(errorsForJoint.get(i)-k+1, errorsForJoint.get(i-1)+1);
                    eFor = Math.min(errorsForJoint.get(i)-1, errorsForJoint.get(i+1)-k);
                }
                if (eFor < bFor)
                    continue;

                int bRev;
                int eRev;
                if (rRev.nucl.charAt(errors.get(i)) != '-')
                {
                    bRev = Math.max(errorsRevJoint.get(i)-k+1, errorsRevJoint.get(i-1)+1);
                    eRev = Math.min(errorsRevJoint.get(i), errorsRevJoint.get(i+1)-k);
                }
                else
                {
                    bRev = Math.max(errorsRevJoint.get(i)-k+1, errorsRevJoint.get(i-1)+1);
                    eRev = Math.min(errorsRevJoint.get(i)-1, errorsRevJoint.get(i+1)-k);
                }
                if (eRev < bRev)
                    continue;

                double avKcountFor = 0;
                double avKcountRev = 0;
                for (int j = bFor; j <= eFor; j++)
                {
                    String kmer = rForJoint.nucl.substring(j, j+k);
                    int kcount = kmers.get(kmer);
                    avKcountFor += kcount;
                }
                avKcountFor /= (eFor - bFor + 1);

                for (int j = bRev; j <= eRev; j++)
                {
                    String kmer = rRevJoint.nucl.substring(j, j+k);
                    int kcount = kmers.get(kmer);
                    avKcountRev += kcount;
                }
                avKcountRev /= (eRev - bRev + 1);

                if (avKcountFor > avKcountRev)
                    sbRev.setCharAt(errors.get(i), sbFor.charAt(errors.get(i)));
                if (avKcountRev > avKcountFor)
                    sbFor.setCharAt(errors.get(i), sbRev.charAt(errors.get(i)));
                
                nErrCorr++;

            }
            catch (Exception e)
            {
                System.out.println(i);
            }
                
        }
        rFor.nucl = new String(sbFor);
        rRev.nucl = new String(sbRev);
        return nErrCorr;
    }
}
