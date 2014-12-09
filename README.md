MIGMA
=====

Written by: Maria Chatzou & Pantelis Bagos


Description
-----------
Briefly, MIGMA is based on a Mixture Transitiom Distribution model
where each lag contributes to the prediction of the current letter in 
a separate and additive way.

For a more detailed description, as well as literature pointers please
see the application note under the doc/ directory, which accompanies the code.


Compile
---------
Compile MIGMA using "make"
Delete MIGMA executable using "make clean"
Delete MIGMA executable and object files using "make wipe"


Notes
-----
Both training and test sequence files must be in FASTA format (see http://tigrblast.tigr.org/web-hmm/fasta.html ).



Example
-------
In the example/ directory you will find the following files:

inputs:
7tm_3		 - a family from Pfam, which we will be used for prediction
7tm_3_TRAIN   - 4/5 of the family sequences, which will be used for training

(you can copy these to a new directory, say example2, and compare results)


_________

Training		
_________

Gerneral:

    ./migma -i_train <input_training_file> -o_mtd <mtd_output_file>  -a <alphabet (1:protein, 2:DNA)> -alg <algorithm (default: EM, 1: EM, 2:Viterbi)> -or <order> 

(in the example/ directory, type the following command line)

In Linux:

    ../migma -i_train 7tm_3_TRAIN.fa -o_mtd 7tm_3.mtd -a 1 -or 20 -alg 1

In Windows:

    ../migma.exe -i_train 7tm_3_TRAIN.fa -o_mtd 7tm_3.mtd -a 1 -or 20 -alg 1


yields the following output files:

7tm_3.mtd    - the MTD file that results from training from 7tm_3_TRAIN.txt



___________

Prediction	
___________

General:

    ./migma -m <mtd_file> -i <input_file> -o <output_file>

(in the example/ directory, type the following command line)

In Linux:

    ../migma  -m 7tm_3.mtd  -i 7tm_3.fa  -o 7tm_3.pred
    
In Windows:

    ../migma.exe  -m 7tm_3.mtd  -i 7tm_3.fa  -o 7tm_3.pred
 
  
yields the following output file:

7tm_3.pred  - the results of the prediction

option:

  -p 	when used it also yields a second output file which contains the predisctions of the amino acid probabilities 



______________________

Training & Prediction:
______________________

General:	

     ./migma -i_train <input_training_file> -o_mtd mtd_output_file> -a <alphabet (1:protein, 2:DNA)> -alg <algorithm (default: EM, 1: EM, 2:Viterbi)> -or <mtd_order> -i <input_file> -o <output_file>

(in the example/ directory, type the following command line)

In Linux:

    ../migma -i_train 7tm_3_TRAIN.fa -o_mtd 7tm_3.mtd -a 1 -or 20 -alg 1 -i 7tm_3.fa -o 7tm_3.pred

In Windows:

    ../migma.exe -i_train 7tm_3_TRAIN.fa -o_mtd 7tm_3.mtd -a 1 -or 20 -alg 1 -i 7tm_3.fa -o 7tm_3.pred


yields the following output files:

7tm_3.mtd    - the MTD file that results from training from 7tm_3_TRAIN.txt

7tm_3.pred   - the results of the prediction using 7tm_3.mtd as the model



Contact information
-------------------
Maria Chatzou
e-mail: maria.chatzou@crg.eu

