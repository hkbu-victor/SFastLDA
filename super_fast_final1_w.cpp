#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <ctime>

//#define UINT32_MAX 4294967295.0

char *dFile, *oFile;
long *square_increment, *square_decrement, *square_of;
int T, NN, SEED, D, W, d_pivotal_topic_count, pivotal_topic_count, common_topic_count, nTokens, square_limit;
int *ztot, *pivotal_topic_indexes, *d_pivotal_topic_indexes, *common_topic_pos, *d, *w, *z, *doc_lens, *word_lens;
int *D_pair_count, *W_pair_count;
double argALPHA, BETA, WBETA, ALPHA_BETA_Z, ONE_OVER_Z_2_norm;
double *ALPHA, *pivotal_topic_values, *d_pivotal_topic_values, *ONE_OVER_Z, *BETA_D_2_norm;
int *wi_pair_val;
int **wpVal;
int *wi_array;
int *wi_index_to_pos;

//******************************************************************************************
time_t t_begin_setup, t_end_setup, t_begin_iteration, t_end_iteration;

struct W_DoublePair{
	int index;
	float w_over_z_value;
};
W_DoublePair *wi_pair_idx, **wpIdx;

struct D_DoublePair{
	int index;
	int value;
};
D_DoublePair *di_pair, **dpPair;

void show_help(){
	printf("\nHKBU SuperFast Gibb's sampling for LDA\n");
	printf("\n word-by-word LDA modeling\n");
	printf("=====================================\n");
	printf("Command line usage:\n");
	printf("super_fast_final1_w -alpha <double> -beta <double> -ntopics <int> -niters <int> -dfile <string> -seed <int> -ofile <string>\n");
	printf("=====================================\n");
	printf("Parameters \t Type \t\t Descriptions\n");
	printf("-alpha\t\t<double> Hyper-parameter for documents-topic distributions\n");
    printf("-beta\t\t<double> Hyper-parameter for word-topic distributions \n");
    printf("-ntopics\t\t<int> Number of topics \n");
    printf("-niters\t\t<int> Number of Gibbs sampling iterations \n");
    printf("-seed\t\t<int> Random seed \n");
    printf("-dfile\t\t<file> Input training data file \n");
	printf("-ofile\t\t<file> The output WP file and Z file (in text format)\n");
}

int parse_args(int argc, char **argv){
	char *arg;

	dFile = NULL; oFile = NULL;
	T = 0; NN = 0; SEED = 0;
	argALPHA = 0.0, BETA = 0.01;

	if (argc < 10){
		show_help();
		return -1;
	}

	int i = 0;
	while (i < argc){
		arg = argv[i];

		if (strcmp(arg,"-dfile")==0){
			dFile = argv[++i];
		} else if (strcmp(arg,"-ofile")==0){
			oFile = argv[++i];
		} else if (strcmp(arg,"-alpha")==0){
			argALPHA = atof(argv[++i]);
		} else if (strcmp(arg,"-beta")==0){
			BETA = atof(argv[++i]);
		} else if (strcmp(arg,"-ntopics")==0){
			T = atoi(argv[++i]);
		} else if (strcmp(arg,"-niters")==0){
			NN = atoi(argv[++i]);
		} else if (strcmp(arg,"-seed")==0){
			SEED = atoi(argv[++i]);
		} else {
            // any more
		}

		i++;
	}
	return 0;
}

/*
int export_wp(char *filename)
{
    int wi, j;

    FILE *fh;
    fh = fopen(filename, "w");
    if (fh == NULL){
      printf("\nFile writting error !!!");
      exit (-1);
    }
    for (wi=0; wi<W; wi++){
    	wi_pair_idx = wpIdx[wi];
    	wi_pair_val = wpVal[wi];
    	for (j=0; j<W_pair_count[wi]; j++)
          fprintf(fh,"%d %d %d\n", wi + 1, wi_pair_idx[j].index + 1, wi_pair_val[j]);
    }
    fclose(fh);
    return 0;
}
*/

int export_results(char *filename)
{
	int i;

	FILE *fh;
	fh = fopen(filename, "w");
    if (fh == NULL){
      printf("\nFile writting error !!!");
      exit (-1);
    }

    for (i=0; i<nTokens; i++){
    	fprintf(fh, "%d %d %d\n", d[ i ] + 1, w[ i ] + 1, z[ i ] + 1);
    }

	fclose(fh);
	return 0;
}

void eval_BETA_D_2_norm(){
	double dtemp;
	for (int di=0; di<D; di++){
		BETA_D_2_norm[ di ] = 0.0;
		di_pair = dpPair[di];
		for (int j=0; j<D_pair_count[di]; j++){
			dtemp = (double) di_pair[j].value * BETA;
			BETA_D_2_norm[di] += dtemp * dtemp;
		}
	}
}

void eval_ALPHA_BETA_Z(){
     ALPHA_BETA_Z = 0.0;
     for (int topic=0; topic<T; topic++){
         ALPHA_BETA_Z += BETA*ALPHA[ topic ] / (WBETA + (double) ztot[ topic ]);
     }
}

void eval_ONE_OVER_Z_2_norm(){
     double dtemp;
     ONE_OVER_Z_2_norm = 0.0;
     for (int topic=0; topic<T; topic++){
         dtemp = 1.0/(WBETA + (double) ztot[topic]);
         ONE_OVER_Z_2_norm += dtemp*dtemp;
     }
}

void eval_ONE_OVER_Z(){
    for (int topic=0; topic<T; topic++){
        ONE_OVER_Z[topic] = 1.0/(WBETA + (double) ztot[topic]);
    }
}

void eval_square_increment(){
    for (long i=0; i<square_limit; i++){
        square_increment[i] = (i+1L)*(i+1L) - i*i;
    }
}

void eval_square(){
	for (long i=0; i<square_limit; i++){
		square_of[i] = i*i;
	}
}

void eval_square_decrement(){
    for (long i=0; i<square_limit; i++){
        square_decrement[i] = (i-1L) * (i-1L) - i * i;
    }
}

struct compare_W_DoublePair{
	bool operator()(const W_DoublePair & lhs, const W_DoublePair & rhs)
	{
		return (lhs.w_over_z_value > rhs.w_over_z_value);
	}
};


struct compare_D_DoublePair{
	bool operator()(const D_DoublePair & lhs, const D_DoublePair & rhs)
	{
		return (lhs.value > rhs.value);
	}
};

void sort_dpPairIdx(){
	int di, pair_count;

	for (di=0; di<D; di++){
		di_pair = dpPair[di];
		pair_count = D_pair_count[di];
		std::sort(di_pair, di_pair + pair_count, compare_D_DoublePair());
	}
}

void eval_likelihood(){

	//assume symmetric alpha and beta
  double log_likelihood = 0.0;
  double ALPHA_BAR = ALPHA[0] * T;
  double BETA_BAR = WBETA;
  int wi, di, j;

  double lgamma_beta = lgamma(BETA);
  double lgamma_alpha = lgamma(ALPHA[0]);

  for (wi=0; wi<W; wi++){
     wi_pair_val = wpVal[wi];
     for (j=0; j<W_pair_count[wi]; j++){
    	log_likelihood += lgamma( BETA + wi_pair_val[ j ] ) - lgamma_beta;
     }
  }

  for (di=0; di<D; di++){
	  di_pair = dpPair[di];
	  for (j=0; j<D_pair_count[di]; j++){
		  log_likelihood += lgamma(ALPHA[0] + di_pair[j].value) - lgamma_alpha;
	  }
  }

  for (j=0; j<T; j++)
	  log_likelihood += lgamma(BETA_BAR) - lgamma(BETA_BAR + ztot[j] );

  for (di=0; di<D; di++){
	  log_likelihood += lgamma(ALPHA_BAR) - lgamma(ALPHA_BAR + doc_lens[di] );
  }

  log_likelihood = log_likelihood / nTokens;

  printf("LOG LIKELIHOOD (per token) : %lf\n", log_likelihood);

}


void GibbsSamplerLDA(){
	int i, j, n, wi, di, topic, itemp, itemp1, vtemp, di_pair_count, topic_temp, pos;
	int i_start, i_mid, i_end, sampled_topic, di_old, di_value_o, di_value_n, wi_value_o, wi_value_n;
	int pivotal_topic_count_o, check_threshold, di_nonzero_index, wi_pair_count;
	long di_2_norm_old, di_2_norm_temp;
	double dtemp1, dtemp2, r, one_over_z_o, one_over_z_n, ALPHA_W_Z, wi_OVER_Z_2_norm;
	double real_BETA_D_Z, Z_r_total, Z_e = 0.0, r_Z_e = 0.0, res1, res2, wi_2_norm_temp;
	double Z_r_total_o, one_over_Z_e_o, Z_r_total_n;
	bool success_sampled;
	int count1, count2, count3, count4, count5;
	//DoublePair *p_src, *p_dest, *p_end;

    printf("Starting Random initialization\n" );
    srand48(2*SEED + 1);
    for (i=0; i<nTokens; i++){
    	z[ i ] = (int) (T * drand48());
    }
    for (i=0; i<nTokens; i++){
        topic = z[ i ];
    	wi = w[ i ];
        wi_pair_idx = wpIdx[wi];
        wi_pair_val = wpVal[wi];
        for (j=0; j<T; j++){
        	if (wi_pair_idx[j].index==-1){
        		wi_pair_idx[j].index = topic;
        		wi_pair_val[j] = 1;
        		W_pair_count[wi]++;
        		break;
        	}else{
        		if (wi_pair_idx[j].index==topic){
        			wi_pair_val[j]++;
        			break;
        		}
        	}
        }

        di = d[i];
        di_pair = dpPair[di];
        for (j=0; j<T; j++){
        	if (di_pair[j].index==-1){
        		di_pair[j].index= topic;
        		di_pair[j].value = 1;
        		D_pair_count[di]++;
        		break;
        	}else{
        		if (di_pair[j].index==topic){
        			di_pair[j].value++;
        			break;
        		}
        	}
        }

        ztot[ topic ]++; // increment ztot matrix
    }
    printf("Random Initialization Finished !!!\n");

    //do initialization and evaluation
    WBETA = W*BETA;
	eval_BETA_D_2_norm();
	eval_ALPHA_BETA_Z();
	eval_ONE_OVER_Z_2_norm();
	eval_ONE_OVER_Z();
	eval_square_increment();
	eval_square_decrement();
	eval_square();

	time(&t_begin_iteration);

	std::fill_n(wi_array, T, 0);
    for (int iter=0; iter<NN; iter++){
    	printf("\nIteration : %d of %d...\n", (iter+1), NN );
    	count1 = 0; count2 = 0; count3 = 0; count4 = 0; count5 = 0;
    	//in this program, Alpha is symmetric and not changed. If you need to revise it, do it here.
        int token_i = 0;
        while (token_i < nTokens){

        	//for each word
        	wi = w[token_i];
        	wi_pair_idx = wpIdx[wi];
        	wi_pair_val = wpVal[wi];
        	ALPHA_W_Z = 0.0;
        	wi_OVER_Z_2_norm = 0.0;
        	wi_pair_count = W_pair_count[wi];
        	for (j=0; j<wi_pair_count; j++){
        		itemp = wi_pair_idx[j].index;
        		vtemp = wi_pair_val[j];
        		wi_array[itemp] = vtemp;
        		dtemp2 = (double) vtemp * ONE_OVER_Z[itemp];
        		ALPHA_W_Z += ALPHA[itemp] * dtemp2;
        		wi_OVER_Z_2_norm += dtemp2 * dtemp2;
        		wi_pair_idx[j].w_over_z_value = (float) dtemp2;
        	}

        	if (iter > 0)
            	std::sort(wi_pair_idx, wi_pair_idx + wi_pair_count, compare_W_DoublePair());

            for (j=0; j<wi_pair_count; j++)
            	wi_index_to_pos[ wi_pair_idx[j].index ] = j;

            sampled_topic = -1;
	        di_old = -1;
	        di_2_norm_old = -1L;

            for (n=0; n<word_lens[wi]; n++){
            	di = d[token_i];
	            di_pair = dpPair[di];
	            di_pair_count = D_pair_count[di];

	            topic = z[ token_i ];

	            wi_value_o = wi_array[topic];
	            wi_array[topic]--;
	            wi_value_n = wi_value_o - 1;

	            pos = 0;
	            while (di_pair[pos++].index != topic); //pos + 1 is the target!!!
	            di_value_o = di_pair[ --pos ].value;
	            di_pair[ pos ].value--;
	            di_value_n = di_value_o - 1;

	            if (di_value_n==0){
	            	//while (di_pair[j++].index != topic); can be identified by "pos"
	            	if (di_pair_count > (pos+1))
	            		memmove(di_pair + pos, di_pair + pos + 1, (di_pair_count - pos - 1)*sizeof(D_DoublePair));
	            	di_pair_count--;

	            	if (di == di_old){
	            		j = pos;
	            		for (i=0; i<common_topic_count; i++){
	            			if (common_topic_pos[i] > j){
	            				common_topic_pos[i]--;
	            			}else{
	            				if (common_topic_pos[i] == j)
	            					pos = i;
	            			}
	            		}
	            		if (common_topic_count > (pos+1))
	            			memmove(common_topic_pos + pos, common_topic_pos + pos + 1, (common_topic_count - pos - 1)*sizeof(int));
	            		common_topic_count--;
	            		di_2_norm_old--;
	            	}
	            }

	            one_over_z_o = ONE_OVER_Z[ topic ];
	            ztot[topic]--;
	            one_over_z_n = 1.0/(WBETA + (double) ztot[topic]);

	            if (di == di_old){
	            	if (wi_value_n==0){
	            		if (di_value_n > 0){
	            			//w=0, d>0
	            			j=0;
	            			while (di_pair[ common_topic_pos[j++] ].index != topic);
	            			if (common_topic_count > j)
	            				memmove(common_topic_pos + j - 1, common_topic_pos + j, (common_topic_count - j)*sizeof(int));
	            			common_topic_count--;
	            			di_2_norm_old -= square_of[ di_value_o ];
	            		}
	            	}else{
	            		if (di_value_n > 0)
	            			di_2_norm_old += square_decrement[ di_value_o ];
	            	}
                }

                dtemp1 = wi_value_o*one_over_z_o;
	            dtemp2 = wi_value_n*one_over_z_n;
	            wi_OVER_Z_2_norm += dtemp2*dtemp2 - dtemp1*dtemp1;
                BETA_D_2_norm[di] += BETA * BETA * (double) square_decrement[ di_value_o ];
	            ONE_OVER_Z[ topic ] = one_over_z_n;
	            ONE_OVER_Z_2_norm += one_over_z_n*one_over_z_n - one_over_z_o*one_over_z_o;
	            ALPHA_W_Z += ALPHA[topic]*(wi_value_n*one_over_z_n - wi_value_o*one_over_z_o);
	            ALPHA_BETA_Z += ALPHA[topic]*BETA*(one_over_z_n - one_over_z_o);

	            r = drand48();
	            Z_r_total = 0.0;
	            res2 = ALPHA_W_Z + ALPHA_BETA_Z;
	            success_sampled = false;
	            pivotal_topic_count = 0;
	            di_2_norm_temp = 0L;
	            wi_2_norm_temp = wi_OVER_Z_2_norm + 1.0e-15;
                Z_r_total_o = 0.0;
                one_over_Z_e_o = 0.0;
                pivotal_topic_count_o = 0;

                res1 = sqrt(BETA_D_2_norm[di]*ONE_OVER_Z_2_norm);

                if (di == di_old){
               		di_2_norm_temp = di_2_norm_old;  //usewi_2_norm_old if doc and word are the same as the previous token
               	}else{
               		common_topic_count = 0;
               		for (di_nonzero_index=0; (di_nonzero_index < di_pair_count); di_nonzero_index++){
               			itemp = di_pair[ di_nonzero_index ].index;
           				if (wi_array[ itemp ] > 0){
           					//common topic
           					di_2_norm_temp += square_of[ di_pair[ di_nonzero_index ].value ];
          					common_topic_pos[ common_topic_count ] = di_nonzero_index;
           					common_topic_count++;
           				}
               		}
               		di_2_norm_old = di_2_norm_temp;
               	}

               	check_threshold = common_topic_count / 10;
               	int common_topic_count_minus_one = common_topic_count - 1;
               	int check_threshold_increment = 5 + common_topic_count / 10;
               	for (di_nonzero_index=0; di_nonzero_index < common_topic_count; di_nonzero_index++){
               		topic_temp = di_pair[common_topic_pos[ di_nonzero_index ]].index;
               		dtemp1 = wi_array[topic_temp] * ONE_OVER_Z[ topic_temp ];
               		wi_2_norm_temp -= dtemp1*dtemp1;
               		itemp1 = di_pair[ common_topic_pos[di_nonzero_index] ].value;
               		di_2_norm_temp -= square_of[ itemp1 ];

       				Z_r_total += dtemp1* (double) itemp1;
       				pivotal_topic_indexes[ pivotal_topic_count ] = topic_temp;
       				pivotal_topic_values[ pivotal_topic_count ] = Z_r_total + Z_r_total_o;
       				pivotal_topic_count++;

               		if (di_nonzero_index >= check_threshold){
               			//check point reached
                		dtemp1 = wi_2_norm_temp * (double) di_2_norm_temp;

                		if (dtemp1 > 1.0e-18)
                  			Z_e = sqrt(dtemp1) + Z_r_total + Z_r_total_o + res1 + res2;
                  		else
                  			Z_e = Z_r_total + Z_r_total_o + res1 + res2;

                  		Z_r_total_n = Z_r_total + (1.0 - Z_e * one_over_Z_e_o) * Z_r_total_o;
                   		//r_n = (r - one_over_Z_e_o * Z_r_total_o);
                  		r_Z_e = (r - one_over_Z_e_o * Z_r_total_o) * Z_e;

                  		if (r_Z_e <= Z_r_total_n){
                  			success_sampled = true;
                  			count1++;
                  			dtemp1 = (1.0 - Z_e * one_over_Z_e_o);
                  			dtemp2 = dtemp1 * Z_r_total_o;
                  			if (r_Z_e <= dtemp2){
                  				r_Z_e = r_Z_e / dtemp1;

                  				i_start = 0;
	                            i_end = pivotal_topic_count_o - 1;
	                            i_mid = (i_start + i_end) / 2;
	                            if (r_Z_e < pivotal_topic_values[0]){
	                               	sampled_topic = pivotal_topic_indexes[0];
	                            }else{
	                               	while (i_start < i_mid){
	                               		if (pivotal_topic_values[i_mid] <= r_Z_e)
	                               			i_start = i_mid;
	                               		else
	                               			i_end = i_mid;
	                               		i_mid = (i_start + i_end) / 2;
	                               	}
	                               	sampled_topic = pivotal_topic_indexes[i_end];
	                            }
	                  		}else{
	                  			r_Z_e -= dtemp2 - Z_r_total_o;

	                  			i_start = pivotal_topic_count_o;
	                           	i_end = pivotal_topic_count - 1;
	                           	i_mid = (i_start + i_end) / 2;
	                          	if (r_Z_e < pivotal_topic_values[pivotal_topic_count_o]){
	                           		sampled_topic = pivotal_topic_indexes[pivotal_topic_count_o];
	                           	}else{
	                           		while (i_start < i_mid){
	                           			if (pivotal_topic_values[i_mid] <= r_Z_e)
	                           				i_start = i_mid;
	                           			else
	                           				i_end = i_mid;
	                           			i_mid = (i_start + i_end) / 2;
	                           		}
	                           		sampled_topic = pivotal_topic_indexes[i_end];
	                           	}
	                  		}
	                  		break;
	                  	}else{
	                  		Z_r_total_o += Z_r_total;
	                  		Z_r_total = 0.0;
	                  		one_over_Z_e_o = 1.0/Z_e;
	                  		pivotal_topic_count_o = pivotal_topic_count;

	               			check_threshold = di_nonzero_index + check_threshold_increment;
	               			if (check_threshold > common_topic_count_minus_one)
	               				check_threshold = common_topic_count_minus_one;
	               		}
	              	}
	            }

	            if (!success_sampled){
	               	d_pivotal_topic_count = 0;
	               	real_BETA_D_Z = 0.0;
	               	for (i=0; i<di_pair_count; i++){
	               		itemp = di_pair[i].index;
	               		real_BETA_D_Z += BETA * di_pair[i].value * ONE_OVER_Z[itemp];
	               		d_pivotal_topic_indexes[ d_pivotal_topic_count ] = itemp;
	               		d_pivotal_topic_values[ d_pivotal_topic_count ] = real_BETA_D_Z;
	               		d_pivotal_topic_count++;
	               	}

	               	Z_e = (Z_r_total + Z_r_total_o + real_BETA_D_Z + res2);
	              	r_Z_e = r * Z_e;
	               	if (r_Z_e < (Z_r_total + Z_r_total_o)){
	               		success_sampled = true;
	               		count2++;
	               		//r_n = (r - one_over_Z_e_o * Z_r_total_o);
	               		r_Z_e = (r - one_over_Z_e_o * Z_r_total_o) * Z_e / (1.0 - Z_e * one_over_Z_e_o);
	               		i_start = 0;
	                   	i_end = pivotal_topic_count_o - 1;
	                   	i_mid = (i_start + i_end) / 2;
	                   	if (r_Z_e < pivotal_topic_values[0]){
	                   		sampled_topic = pivotal_topic_indexes[0];
	                   	}else{
	                   		while (i_start < i_mid){
	                   			if (pivotal_topic_values[i_mid] <= r_Z_e)
	                   				i_start = i_mid;
	                   			else
	                   				i_end = i_mid;
	                   			i_mid = (i_start + i_end) / 2;
	                   		}
	                   		sampled_topic = pivotal_topic_indexes[i_end];
	                   	}
	               	}else{
	                   		//r_z_e is the actural Z value
	               		r_Z_e -= (Z_r_total + Z_r_total_o);
   						if (r_Z_e <= real_BETA_D_Z){
   							success_sampled = true;
   							count3++;
   	  						i_start = 0;
   	  						i_end = d_pivotal_topic_count - 1;
   	  						i_mid = (i_start + i_end) / 2;
   	  						if (r_Z_e < d_pivotal_topic_values[0]){
   	  							sampled_topic = d_pivotal_topic_indexes[0];
   	  						}else{
   	  							while (i_start < i_mid){
   	  								if (d_pivotal_topic_values[i_mid] <= r_Z_e)
   	  									i_start = i_mid;
   	  								else
   	  									i_end = i_mid;
   	  								i_mid = (i_start + i_end)/2;
   	  							}
   	  							sampled_topic = d_pivotal_topic_indexes[i_end];
   	  						}
   						}else
   							r_Z_e -= real_BETA_D_Z + 1.0e-15; //1.0e-10 is for rounding error
	                }
                }

	            if (!success_sampled){
	            	if (r_Z_e <= ALPHA_W_Z){
	            		count4++;
	   					for (i=0; i<wi_pair_count; i++){
	   						sampled_topic = wi_pair_idx[i].index;
   							r_Z_e -= ALPHA[sampled_topic] * wi_array[sampled_topic]*ONE_OVER_Z[sampled_topic];
	   						if (r_Z_e <= 0.0 ){
	   							break;
	   						}
	   					}
	   				}else{
	   					count5++;
	   					r_Z_e -= ALPHA_W_Z;
	   					r_Z_e /= BETA;
	   					for (sampled_topic=0; sampled_topic<T; sampled_topic++){
	   						r_Z_e -= ALPHA[sampled_topic]*ONE_OVER_Z[sampled_topic];
	   						if (r_Z_e <= 0.0){
	   							break;
	   						}
	   					}
	   				}
	   			}

	            z[ token_i ] = sampled_topic;
	            bool common_topic_added = false;
	            pos = 0;
	            while ((di_pair[pos].index != sampled_topic) && (pos<di_pair_count))
	            	pos++;
	            if (pos == di_pair_count){
	            	//new topic
	               	di_pair[di_pair_count].index = sampled_topic;
	               	di_pair[di_pair_count].value = 0;
	               	common_topic_pos[ common_topic_count ] = di_pair_count;

	               	di_pair_count++;
	               	common_topic_count++;
	               	common_topic_added = true;
	            }
	            di_value_o = di_pair[ pos ].value;
	            di_pair[ pos ].value++;
	            di_value_n = di_value_o + 1;

	            if (wi_index_to_pos[sampled_topic]==-1){
	                //new topic
	               	wi_pair_idx[ wi_pair_count ].index = sampled_topic;
	               	wi_index_to_pos[sampled_topic] = wi_pair_count;
	               	wi_pair_count++;

	               	if (!common_topic_added){
	               		common_topic_pos[ common_topic_count ] = pos;
	               		common_topic_count++;
	               		common_topic_added = true;
	               	}
	            }

	            if (common_topic_added)
	            	di_2_norm_old += square_of[di_value_n];
	            else
	            	di_2_norm_old += square_increment[ di_value_o ];

	            wi_value_o = wi_array[sampled_topic];
	            wi_array[sampled_topic]++;
	            wi_value_n = wi_value_o + 1;

	            one_over_z_o = ONE_OVER_Z[ sampled_topic ];
	            ztot[sampled_topic]++;
	            one_over_z_n = 1.0/(WBETA + (double) ztot[sampled_topic]);

	            dtemp1 = wi_value_o*one_over_z_o;
	            dtemp2 = wi_value_n*one_over_z_n;
	            wi_OVER_Z_2_norm += dtemp2*dtemp2 - dtemp1*dtemp1;
	            BETA_D_2_norm[di] += BETA * BETA * (double) square_increment[ di_value_o ];
	            ONE_OVER_Z[ sampled_topic ] = one_over_z_n;
	            ONE_OVER_Z_2_norm += one_over_z_n*one_over_z_n - one_over_z_o*one_over_z_o;
	            ALPHA_W_Z += ALPHA[sampled_topic]*(wi_value_n*one_over_z_n - wi_value_o*one_over_z_o);
	            ALPHA_BETA_Z += ALPHA[sampled_topic]*BETA*(one_over_z_n - one_over_z_o);
	            token_i++;
	            di_old = di;
	            D_pair_count[di] = di_pair_count;
	        }
	        //end of a word
/*
	        for (j=0; j<wi_pair_count; j++){
	        	wi_pair_val[j] = wi_array[wi_pair_idx[j].index];
	        	wi_array[wi_pair_idx[j].index] = 0;
	        }
	        W_pair_count[wi] = wi_pair_count;
*/
            i = 0;
        	for (j=0; j<wi_pair_count; j++){
        		itemp = wi_pair_idx[j].index;
        		if (wi_array[itemp] > 0){
        			wi_pair_idx[i].index = wi_pair_idx[j].index;
        			wi_pair_val[i] = wi_array[itemp];
        			wi_array[itemp] = 0;
        			i++;
        		}
        		wi_index_to_pos[ itemp ] = -1;
        	}
        	W_pair_count[wi] = i;
	    }//end of an iter
        if (iter > 0)
        	sort_dpPairIdx();
	    time(&t_end_iteration);
	    printf("Total Gibbs iteration time used : %f\n", difftime(t_end_iteration, t_begin_iteration));
	    //printf("%f, %f, %f, %f, %f\n", (double)count1/nTokens, (double)count2/nTokens, (double)count3/nTokens, (double)count4/nTokens, (double)count5/nTokens);
	    eval_likelihood();
    }
	       //end of all iterations
}
//************************************************************************************************



int main(int argc, char **argv){
  	char line[ 80 ];
    int n, i, j, docId, wordId, count;
    FILE *fr;

    parse_args(argc, argv);

    square_limit = -1;
    printf("\nFinal Sparse word-by-word: Opening and reading dataset file...\n");

    D = 0; W = 0; nTokens = 0;
    fr = fopen(dFile,"rt");
    if (fr == NULL){
       	printf("\n\nCannot find the UCI bag-of-words format input corpus file %s !!!\n\n\n", dFile);
       	exit(1);
    }
    fgets(line, 80, fr);
    sscanf(line, "%d", &D);
    fgets(line, 80, fr);
    sscanf(line, "%d", &W);
    fgets(line, 80, fr);

    int lineNo = 3;

    while (fgets(line, 80, fr) != NULL){
      	sscanf( line, "%d %d %d", &docId, &wordId, &count);
      	lineNo++;
      	nTokens += count;
      	if ((lineNo % 2000000) == 0)
			printf("first parse : processing line %d \n", lineNo);
    }
    printf("\nTotal documents: %d, Vocab Size : %d, Num. Tokens : %d\n", D, W, nTokens);

    //allocate memory
    d = (int *) calloc(nTokens, sizeof(int));
    w = (int *) calloc(nTokens, sizeof(int));
    z = (int *) calloc(nTokens, sizeof(int));
    doc_lens = (int *) calloc(D, sizeof(int));
    word_lens = (int *) calloc(W, sizeof(int));

    rewind(fr);

    fgets(line, 80, fr); fgets(line, 80, fr); fgets(line, 80, fr);
    lineNo = 3;
    n = 0;
    while (fgets(line, 80, fr) != NULL){
    	sscanf( line, "%d %d %d", &docId, &wordId, &count);
    	lineNo++;
    	doc_lens[ docId-1 ] += count;
    	word_lens[ wordId-1 ] += count;
    	for (i=0; i<count; i++){
    		d[n] = docId - 1;
    		w[n] = wordId - 1;
    		n++;
    	}
    	if ((lineNo % 1000000) == 0)
    		printf("second parse : processing line %d \n", lineNo);
    }
    fclose(fr);

    for (i=0; i<W; i++){
    	if (word_lens[i] > square_limit)
    		square_limit = word_lens[i];
    }

    printf( "Arguments:\n" );
    printf( "\tNumber of words      W = %d\n"    , W );
    printf( "\tNumber of docs       D = %d\n"    , D );
    printf( "\tNumber of topics     T = %d\n"    , T );
    printf( "\tNumber of iterations N = %d\n"    , NN );
    printf( "\tHyperparameter   ALPHA = %4.4f\n" , argALPHA );
    printf( "\tHyperparameter    BETA = %4.4f\n" , BETA );
    printf( "\tSeed number            = %d\n"    , SEED );
    printf( "\tNumber of tokens       = %d\n"    , nTokens );
    printf( "\tSquare Limit           = %d\n"    , square_limit);

    if (oFile != NULL){
    	printf("\n\tOutput file: %s.\n\n", oFile);
    }else{
    	printf("\n\tNo output file specified.\n\n");
    }


    ztot = new int[T]; std::fill_n(ztot, T, 0);
    ALPHA = new double[T]; std::fill_n(ALPHA, T, argALPHA);
    BETA_D_2_norm = new double[D];	//summation over topics of ALPHA*wi
    D_pair_count = new int[D]; std::fill_n(D_pair_count, D, 0);
    W_pair_count = new int[W]; std::fill_n(W_pair_count, W, 0);
    pivotal_topic_indexes = new int[T];
    pivotal_topic_values = new double[T];
    d_pivotal_topic_indexes = new int[T];
    d_pivotal_topic_values = new double[T];
    common_topic_pos = new int[T];
    ONE_OVER_Z = new double[T];
    square_increment = new long[square_limit];
    square_decrement = new long[square_limit];
    square_of = new long[square_limit];
    wi_array = new int[T];
    wi_index_to_pos = new int[T];

    wpIdx = new W_DoublePair*[W];
    wpVal = new int*[W];
    for (i=0; i<W; i++){
   		wpIdx[i] = new W_DoublePair[T];
    	wpVal[i] = new int[T];
    	for (j=0; j<T; j++){
    		wpIdx[i][j].index = -1;
    	}
		std::fill_n(wpVal[i], T, 0);
    }

    dpPair = new D_DoublePair*[D];
    for (i=0; i<D; i++){
    	if (doc_lens[i] < T){
        	dpPair[i] = new D_DoublePair[doc_lens[i]];
        	di_pair = dpPair[i];
    		for (j=0; j<doc_lens[i]; j++){
    			di_pair[j].index = -1;
    			di_pair[j].value = 0;
    		}
    	}else{
        	dpPair[i] = new D_DoublePair[T];
        	di_pair = dpPair[i];
    		for (j=0; j<T; j++){
    			di_pair[j].index = -1;
    			di_pair[j].value = 0;
    		}
    	}
    }

    std::fill_n(wi_index_to_pos, T, -1);

    GibbsSamplerLDA();

    if (oFile != NULL){
    	export_results(oFile);
    	printf("\nOutput file is written to %s", oFile);
    }

    printf("\nFinished !!!\n");
}
