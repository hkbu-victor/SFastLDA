NOTE:

    A. The program should be compiled with GCC 4.8.2 or above.

    B. To create the executable file:
        1. Change to the directory containing the source file and Makefile.
        2. Type "make" at the command line.

    C. The programs only accept the UCI bag-of-words format (please refer to https://archive.ics.uci.edu/ml/datasets/Bag+of+Words). 

    D. The program super_fast_final1_d only accept doc-by-doc UCI bag-of-words format.

    E. The program super_fast_final1_w only accept word-by-word UCI bag-of-words format.

    F. A conversion tool "doc_uci_to_w_uci" is provided here to convert the doc-by-doc UCI dataset to the word-by-word format data file.
   

USAGE:

    1. Command for modeling doc-by-doc datasets:
          super_fast_final1_d -alpha <double> -beta <double> -ntopics <int> -niters <int> -dfile <string> -seed <int> -ofile <string>

    2. Command for modeling word-by-word datasets:
          super_fast_final1_w -alpha <double> -beta <double> -ntopics <int> -niters <int> -dfile <string> -seed <int> -ofile <string>

    3. Command for converting doc-by-doc UCI data file to word-by-word UCI data file:
          doc_uci_to_w_uci -inputFile <string> -outputFile <string>
          e.g. doc_uci_to_w_uci -inputFile abc.txt -outputFile abc_w.txt

    Remark: i) If the program arguement -ofile is not given, no output file will be created.
           ii) If -ofile is specified (e.g., -ofile results), one output file will be created after going through the specified Gibbs isterations:
                 The file with ".Z" extension gives the topic index values of the token w_{d,n} (i.e. the n-th word of document d). Each line
                    has 3 entries corresponding to
                             doc_id word_id topic_id

LICENSE:

  HKBU

CITATION:

  Victor C.W. Cheng and William K. Cheung, SFastLDA: Super Fast Latent Dirichlet Allocation [Software], 2016, available from
http://jklsdjkldjkl
 

ENQUIRIES:

    Please email to victor@comp.hkbu.edu.hk        
