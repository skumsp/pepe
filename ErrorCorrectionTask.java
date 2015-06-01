/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepe;

import ErrorCorrection.Kmer;
import ErrorCorrection.Read;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.Callable;

/**
 *
 * @author kki8
 */
public class ErrorCorrectionTask implements Callable
{
    Read rFor;
    Read rRev;
    Read rForJoint;
    Read rRevJoint;
    int k;
    ErrorCorrectionTask(Read r1, Read r2, Read r3, Read r4, int i)
    {
        rFor = r1;
        rRev = r2;
        rForJoint = r3;
        rRevJoint = r4;
        k = i;
    }
    public Object call()
    {
//        if (rFor.name.equalsIgnoreCase(">40410/1 reference=one_1_1 position=1..306 errors=292%A"))
//                System.out.println();
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
        
        HashMap<Integer,Integer> kmersFor = new HashMap();
        HashMap<Integer,Integer> kmersRev = new HashMap();
        for (Kmer km : rForJoint.kmers)
            kmersFor.put(km.getStartPos(), km.getKCount());
        for (Kmer km : rRevJoint.kmers)
            kmersRev.put(km.getStartPos(), km.getKCount());
        
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
                    avKcountFor += kmersFor.get(j);
                avKcountFor /= (eFor - bFor + 1);

                for (int j = bRev; j <= eRev; j++)
                    avKcountRev += kmersRev.get(j);
                avKcountRev /= (eRev - bRev + 1);

                if (avKcountFor > avKcountRev)
                    sbRev.setCharAt(errors.get(i), sbFor.charAt(errors.get(i)));
                if (avKcountRev > avKcountFor)
                    sbFor.setCharAt(errors.get(i), sbRev.charAt(errors.get(i)));
                
                nErrCorr++;
                
                
                
/*                FileWriter fw = new FileWriter("read_corr" + i + ".fas");
                fw.write(">uncor1" + "\n" + rFor.nucl + "\n" + ">uncor2" + "\n" + rRev.nucl + "\n" 
                        + ">cor1" + "\n" + sbFor + "\n" + ">cor2" + "\n" + sbRev + "\n");
                fw.close(); */
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
