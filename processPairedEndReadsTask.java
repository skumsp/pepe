/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepe;

import ErrorCorrection.DataSet;
import ErrorCorrection.Read;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;
import java.util.concurrent.Callable;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;

/**
 *
 * @author kki8
 */
public class processPairedEndReadsTask implements Callable 
{
    Read rForw;
    Read rRev;
    int gapop;
    int gapext;
    public  processPairedEndReadsTask() 
    {
        
    }
   public  processPairedEndReadsTask(Read r1, Read r2,int gapop1,int gapext1) 
   { 
       rForw = r1;
       rRev = r2;
       gapop = gapop1;
       gapext = gapext1;
   }

   public Object call() throws FileNotFoundException, IOException 
   {
        int nSub = 0;
        int nIndel = 0;
        DNASequence target = new DNASequence(rForw.nucl, AmbiguityDNACompoundSet.getDNACompoundSet());
        DNASequence query = new DNASequence(rRev.nucl, AmbiguityDNACompoundSet.getDNACompoundSet());
 
        SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();
        
 
	SimpleGapPenalty gapP = new SimpleGapPenalty();
	gapP.setOpenPenalty((short)gapop);
	gapP.setExtensionPenalty((short)gapext);
 
        SequencePair<DNASequence, NucleotideCompound> psa = null;
        
        try
        {
            psa =
				Alignments.getPairwiseAlignment(query, target,
						Alignments.PairwiseSequenceAlignerType.LOCAL, gapP, matrix);
        }
        catch(Exception e)      
        {
            ArrayList<Integer> errors = new ArrayList();
            errors.add(-1);
            return errors;
        }
        

//        System.out.println(psa.getAlignedSequence(1));
//        System.out.println(psa.getAlignedSequence(2));
        
            int st = 1;
            while (psa.hasGap(st))
                st++;
            int end = psa.getLength();
            while (psa.hasGap(end))
                end--;

            String seq1 = "";
            String seq2 = "";
            for (int i = st; i <= end; i++)
            {
                seq1+=psa.getCompoundAt(1, i);
                seq2+=psa.getCompoundAt(2, i);
                if (psa.getCompoundAt(1, i) != psa.getCompoundAt(2, i))
                {
                    if (psa.hasGap(i))
                        nIndel++;
                    else
                        nSub++;
                }
            }
            rForw.nucl = seq1;
            rRev.nucl = seq2;
            ArrayList<Integer> errors = new ArrayList();
            errors.add(nIndel);
            errors.add(nSub);
            errors.add(end-st+1);
            return errors;
   }
}
