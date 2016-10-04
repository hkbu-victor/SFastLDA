INTRODUCTION:

    SFastLDA is a super fast implementation of the standard collapsed Gibbs sampling (CGS) for latent Dirichlet allocation (LDA). In contrast to some approximate implementations such as LightLDA or WarpLDA, it is an exact implementation of the CGS. Hence, the models inferred have the same quality as that from CGS.
    
    We have compared SFastLDA with WarpLDA (on ENRON, NYTimes, and mPubMed datasets) and find that 400 iterations are sufficient for SFastLDA to infer models better than that learnt with WarpLDA after running 5,000 iterations.
    
    SFastLDA has two main programs: super_fast_final1_d.cpp accepts data files with document-by-document UCI bag-of-words format and super_fast_final1_w.cpp deals with word-by-word UCI bag-of-words format. In our experiments, super_fast_final1_w is 2-5 times faster than superfast_final1_d. A tool is provided here to convert the document-by-document format data file to word-by-word format data file.
    
    The current version of SFastLDA only uses a single thread of the CPU, next version will support multi-threads.
    
NOTE:

    A. The program should be compiled with GCC 4.8.2 or above and assumed to be executed on 64-bit Linux platforms.

    B. To create the executable file:
        1. Change to the directory containing the source files and Makefile.
        2. Type "make" at the command line.

    C. The programs only accept the UCI bag-of-words format (please refer to https://archive.ics.uci.edu/ml/datasets/Bag+of+Words). 

    D. The program super_fast_final1_d only accept doc-by-doc UCI bag-of-words format.

    E. The program super_fast_final1_w only accept word-by-word UCI bag-of-words format.

    F. A conversion tool "doc_uci_to_word_uci" is provided here to convert the doc-by-doc UCI dataset to the word-by-word format data file.
    
    g. Some document-by-document UCI bag-of-words datasets can be obtained from https://archive.ics.uci.edu/ml/machine-learning-databases/bag-of-words/ .

USAGE:

    1. Command for modeling doc-by-doc datasets:
        super_fast_final1_d -alpha <double> -beta <double> -ntopics <int> -niters <int> -dfile <string> -seed <int> -ofile <string>
        e.g. ./super_fast_final1_d -alpha 0.125 -beta 0.01 -ntopics 400 -niters 400 -seed 1 -dfile enron_orig.txt -ofile m.txt            
    2. Command for modeling word-by-word datasets:
        super_fast_final1_w -alpha <double> -beta <double> -ntopics <int> -niters <int> -dfile <string> -seed <int> -ofile <string>
        e.g. ./super_fast_final1_w -alpha 0.125 -beta 0.01 -ntopics 400 -niters 400 -seed 1 -dfile enron_word.txt -ofile m.txt            
    3. Command for converting doc-by-doc UCI bag-of-words data file to word-by-word UCI bag-of-words data file:
        doc_uci_to_w_uci -inputFile <string> -outputFile <string>
        e.g. ./doc_uci_to_word_uci -inputFile enron_orig.txt -outputFile enron_word.txt

    Remark: i) If the program option -ofile is not specified, no output file will be created.
            ii) If -ofile is specified (e.g., -ofile results.txt), the output file will be created after going through the specified Gibbs isterations:
            iii) The output file has a number of lines and each of them corresponds to one word token appearing in the input data file. The format of the output file is as follows:
            
                             doc_id word_id topic_id
                             doc_id word_id topic_id
                             ......
                             ......
                             ......
                             

LICENSE:

  HKBU

CITATION:

  Victor C.W. Cheng and William K. Cheung, SFastLDA: Super Fast Latent Dirichlet Allocation (SFastLDA) [Software], 2016, available from
https://github.com/hkbu-victor/SFastLDA .
 

ENQUIRIES:

    Please email to victor@comp.hkbu.edu.hk        
