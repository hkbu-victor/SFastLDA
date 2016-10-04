#include <stdio.h>
#include <algorithm>
#include <string.h>

struct SparseData{
	int doc_id;
	int word_id;
	int count;
};

struct compare_DoublePair{
	bool operator()(const SparseData & lhs, const SparseData & rhs)
	{
		if (lhs.word_id == rhs.word_id)
			return (lhs.doc_id < rhs.doc_id);
		return (lhs.word_id < rhs.word_id);
	}
};

char *inputFile;
char *outputFile;

void show_help(){
	printf("\nUCI data set conversion tool:\n");
	printf("\n doc-by-doc UCI data set ---> word-by-word data set\n");
	printf("=====================================\n");
	printf("Command line usage:\n");
	printf("doc_uci_to_w_uci -inputFile <string> -outputFile <string>\n");
	printf("=====================================\n");
}

int parse_args(int argc, char **argv){
	char *arg;

	inputFile = NULL;
	outputFile = NULL;

	if (argc < 5){
		show_help();
		return -1;
	}

	int i = 0;
	while (i < argc){
		arg = argv[i];

		if (strcmp(arg,"-inputFile")==0){
			inputFile = argv[++i];
		} else if (strcmp(arg,"-outputFile")==0){
			outputFile = argv[++i];
		} else {
            // any more
		}
		i++;
	}
	if (strcmp(inputFile, outputFile) == 0){
		printf("\n***Input file and output file have the same path and name !!!\n");
		return -1;
	}
	return 0;
}


int main(int argc, char **argv){
  	char line[ 80 ];
    int i, docId, wordId, count;
    int nlines, D, W;

    SparseData *data_lines;
    FILE *fr;

    i = parse_args(argc, argv);

    if (i < 0)
    	return -1;

    nlines = 0;
    printf("\n*** Conversion Tool: doc-by-doc UCI --->  word-by-word UCI");
    printf("\n*** Opening and reading dataset file %s...\n", inputFile);

    fr = fopen(inputFile,"rt");
    if (fr == NULL){
          	printf("Cannot find the file %s !!!\n", inputFile);
          	exit(1);
     }
    fgets(line, 80, fr);
    sscanf(line, "%d", &D);
    fgets(line, 80, fr);
    sscanf(line, "%d", &W);
    fgets(line, 80, fr);
    sscanf(line, "%d", &nlines);

    data_lines = new SparseData[nlines];

    int lineNo = 0;
    while (fgets(line, 80, fr) != NULL){
      	sscanf( line, "%d %d %d", &docId, &wordId, &count);

      	data_lines[lineNo].doc_id = docId;
      	data_lines[lineNo].word_id = wordId;
      	data_lines[lineNo].count = count;
      	lineNo++;

      	if (((lineNo+3) % 2000000) == 0)
			printf("processing line %d \n", lineNo);
    }
    printf("Total data lines read : %d\n\n", lineNo + 3);
    fclose(fr);

    std::sort(data_lines, data_lines + nlines, compare_DoublePair());


    FILE *fh;
    fh = fopen(outputFile, "w");
        if (fh == NULL){
          printf("\nFile writting error : %s \n!!!", outputFile);
          exit (-1);
        }
        //write 3 lines of meta-information
        printf("*** Write word-by-word UCI data to : %s.\n\n", outputFile);
        fprintf(fh, "%d\n", D);
        fprintf(fh, "%d\n", W);
        fprintf(fh, "%d\n", nlines);
        for (i=0; i<nlines; i++){
        	fprintf(fh, "%d %d %d\n", data_lines[i].doc_id, data_lines[i].word_id, data_lines[i].count);
        	if ((i % 2000000) == 0)
        		printf("writing line %d \n", i);
        }

        fclose(fh);
        printf("Total data lines written : %d\n", i + 3);
    printf("\nFinished!!!\n");
}




