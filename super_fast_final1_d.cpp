#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <assert.h>

//#define UINT32_MAX 4294967295.0

char *dFile, *oFile;
long *square_increment, *square_decrement, *square_of;
int T, NN, SEED, D, W, w_pivotal_topic_count, pivotal_topic_count, common_topic_count, nTokens, square_limit;
int *ztot, *pivotal_topic_indexes, *w_pivotal_topic_indexes, *common_topic_indexes, *d, *w, *doc_lens, *word_lens;
//victor int *z;
int *z;
int *D_pair_count, *W_pair_count;
double argALPHA, BETA, WBETA, ALPHA_BETA_Z, ONE_OVER_Z_2_norm;
double *ALPHA, *ALPHA_square, *pivotal_topic_values, *w_pivotal_topic_values, *ONE_OVER_Z, *ALPHA_W_2_norm;
//victor int *dp_pair_val;
int *dp_pair_val;
int *wi_array;
//victor int *wi_idx;
int *wi_idx;
int **wpVal;
//victor int **wpIdx;
int **wpIdx;
//victor int **dpVal;
int **dpVal;
int *di_array;
int success_count1, success_count2, success_count3;


time_t t_begin_setup, t_end_setup, t_begin_iteration, t_end_iteration;

struct DoublePair{
	//victor int index;
	int index;
	//victor double d_over_z_value;
	float d_over_z_value;
};
DoublePair *dp_pair_idx, **dpIdx;

void show_help(){
	printf("\nHKBU SuperFast Gibb's sampling for LDA\n");
	printf("\n doc-by-doc LDA modeling\n");
	printf("=====================================\n");
	printf("Command line usage:\n");
	printf("super_fast_final1_d -alpha <double> -beta <double> -ntopics <int> -niters <int> -dfile <string> -seed <int> -ofile <string>\n");
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
    	wi_idx = wpIdx[wi];
    	wi_array = wpVal[wi];
    	for (j=0; j<W_pair_count[wi]; j++)
          fprintf(fh,"%d %d %d\n", wi + 1, wi_idx[j] + 1, wi_array[wi_idx[j]]);
    }
    fclose(fh);
    return 0;
}
*/

void evaluate_ALPHA_square(){
	for (int i=0; i<T; i++)
		ALPHA_square[i] = ALPHA[i]*ALPHA[i];
}

void eval_ALPHA_W_2_norm(){
    double dtemp;
    for (int wi=0; wi<W; wi++){
        ALPHA_W_2_norm[ wi ] = 0.0;
        wi_idx = wpIdx[wi];
        wi_array = wpVal[wi];
        for (int j=0; j<W_pair_count[wi]; j++){
         	dtemp = (double) wi_array[ wi_idx[j] ] * ALPHA[ wi_idx[j] ];
            ALPHA_W_2_norm[ wi ] += dtemp*dtemp;
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

struct compare_DoublePair{
	bool operator()(const DoublePair & lhs, const DoublePair & rhs)
	{
		return (lhs.d_over_z_value > rhs.d_over_z_value);
	}
};

void eval_likelihood(){

	//assume symmetric alpha and beta
  double log_likelihood = 0.0;
  double ALPHA_BAR = ALPHA[0] * T;
  double BETA_BAR = WBETA;
  int wi, di, j;

  double lgamma_beta = lgamma(BETA);
  double lgamma_alpha = lgamma(ALPHA[0]);

  for (wi=0; wi<W; wi++){
     wi_idx = wpIdx[wi];
     wi_array = wpVal[wi];
     for (j=0; j<W_pair_count[wi]; j++){
     	log_likelihood += lgamma(BETA + wi_array[ wi_idx[j] ]) - lgamma_beta;
     }
  }

  for (di=0; di<D; di++){
	  dp_pair_idx = dpIdx[di];
	  dp_pair_val = dpVal[di];

	  for (j=0; j<D_pair_count[di]; j++){
		  log_likelihood += lgamma(ALPHA[0] + dp_pair_val[j]) - lgamma_alpha;
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
	int i, j, n, wi, di, topic, itemp, itemp1, vtemp, di_pair_count, topic_temp;
	int i_start, i_mid, i_end, sampled_topic, wi_old, di_value_o, di_value_n, wi_value_o;
	int pivotal_topic_count_o, check_threshold, di_nonzero_index, wi_pair_count;
	long wi_2_norm_old, wi_2_norm_temp;
	double dtemp1, dtemp2, r, one_over_z_o, one_over_z_n, BETA_D_Z, di_OVER_Z_2_norm;
	double real_ALPHA_W_Z, Z_r_total, Z_e = 0.0, r_Z_e = 0.0, res1, res2, di_2_norm_temp;
	double Z_r_total_o, one_over_Z_e_o, Z_r_total_n;
	bool success_sampled;
	//DoublePair *p_src, *p_dest, *p_end;

    printf("Starting Random initialization\n" );
    srand48(2*SEED + 1);
    for (i=0; i<nTokens; i++){
    	z[ i ] = (int) (T * drand48());
    }
    for (i=0; i<nTokens; i++){
        topic = z[ i ];
    	wi = w[ i ];
        wi_idx = wpIdx[wi];
        wi_array = wpVal[wi];
        for (j=0; j<T; j++){
        	if (wi_idx[j] == -1){
        		wi_idx[j] = topic;
        		wi_array[topic] = 1;
        		W_pair_count[wi]++;
        		break;
        	}else{
        		if (wi_idx[j] == topic){
        			wi_array[topic]++;
        			break;
        		}
        	}
        }

        di = d[i];
        dp_pair_idx = dpIdx[di];
        dp_pair_val = dpVal[di];
        for (j=0; j< T; j++ ){
        	if (dp_pair_idx[j].index == -1){
        		dp_pair_idx[j].index = topic;
                dp_pair_val[j] = 1;
                D_pair_count[di]++;
                break;
            }else{
                if (dp_pair_idx[j].index == topic){
                   dp_pair_val[j]++;
                   break;
                }
            }
        }
        ztot[ topic ]++; // increment ztot matrix
    }
    printf("Random Initialization Finished !!!\n");

	    //do initialization and evaluation
    WBETA = W*BETA;
	        //initialize di_index_to_pos to -1 once only as it will be restored to -1 when a document is processed.
	eval_ALPHA_W_2_norm();
	eval_ALPHA_BETA_Z();
	eval_ONE_OVER_Z_2_norm();
	eval_ONE_OVER_Z();
	eval_square_increment();
	eval_square_decrement();
	eval_square();

	time(&t_begin_iteration);

	std::fill_n(di_array, T, 0);
    for (int iter=0; iter<NN; iter++){
    	success_count1 = 0;
    	success_count2 = 0;
    	success_count3 = 0;
    	printf("\nIteration : %d of %d...\n", (iter+1), NN );

    	//in this program, Alpha is symmetric and not changed. If you need to revise it, do it here.
	    evaluate_ALPHA_square();
        int token_i = 0;
        while (token_i < nTokens){
        	//for each document
        	di = d[token_i];
	        BETA_D_Z = 0.0;
	        dp_pair_idx = dpIdx[di];
	        dp_pair_val = dpVal[di];
	        di_OVER_Z_2_norm = 0.0;
	        di_pair_count = D_pair_count[di];
            for (j=0; j<di_pair_count; j++){
	            itemp = dp_pair_idx[j].index;
	            vtemp = dp_pair_val[j];
	            di_array[itemp] = vtemp;
	            dtemp2 = (double)vtemp * ONE_OVER_Z[itemp];
	            BETA_D_Z += dtemp2;
	            di_OVER_Z_2_norm += dtemp2*dtemp2;
	            dp_pair_idx[j].d_over_z_value = dtemp2;
	        }
            BETA_D_Z *= BETA;

            //if ((iter > 3) && (iter % 2)==0)
            if (iter > 3)
            	std::sort(dp_pair_idx, dp_pair_idx + di_pair_count, compare_DoublePair());

            sampled_topic = -1;
	        wi_old = -1;
	        wi_2_norm_old = -1L;
            for (n=0; n<doc_lens[di]; n++){
            	wi = w[token_i];
            	wi_pair_count = W_pair_count[wi];
	            wi_idx = wpIdx[wi];
	            wi_array = wpVal[wi];
	            topic = z[ token_i ];

	            di_value_o = di_array[topic];
	            di_array[topic]--;
	            di_value_n = di_value_o - 1;

	            if (di_value_n==0){
	            	j=0;
	            	while (dp_pair_idx[j++].index != topic);
                   	memmove(dp_pair_idx + j - 1, dp_pair_idx + j, (di_pair_count - j)*sizeof(DoublePair));
                   	di_pair_count--;
                }

	            wi_value_o = wi_array[topic];
	            wi_array[topic]--;
	            if (wi_array[topic] == 0){
	            	j=0;
	            	while (wi_idx[j++] != topic);
	            	//victor memmove(wi_idx + j - 1, wi_idx + j, (wi_pair_count - j)*sizeof(int));
	            	memmove(wi_idx + j - 1, wi_idx + j, (wi_pair_count - j)*sizeof(int));
	            	wi_pair_count--;
	            }

	            one_over_z_o = ONE_OVER_Z[ topic ];
	            ztot[topic]--;
	            one_over_z_n = 1.0/(WBETA + (double) ztot[topic]);

	            if (wi == wi_old){
	            	if (di_value_n==0){
	            		j=0;
	            		while (common_topic_indexes[j++] != topic);
	            		memmove(common_topic_indexes + j - 1, common_topic_indexes + j, (common_topic_count - j)*sizeof(int));
	                   	common_topic_count--;
	            		wi_2_norm_old -= square_of[ wi_value_o ];
	            	}else
	            		wi_2_norm_old += square_decrement[ wi_value_o ];
                }

                dtemp1 = di_value_o*one_over_z_o;
	            dtemp2 = di_value_n*one_over_z_n;
	            di_OVER_Z_2_norm += dtemp2*dtemp2 - dtemp1*dtemp1;
                ALPHA_W_2_norm[wi] += ALPHA_square[topic] * (double) square_decrement[ wi_value_o ];
	            ONE_OVER_Z[ topic ] = one_over_z_n;
	            ONE_OVER_Z_2_norm += one_over_z_n*one_over_z_n - one_over_z_o*one_over_z_o;
	            BETA_D_Z += BETA*(di_value_n*one_over_z_n - di_value_o*one_over_z_o);
	            ALPHA_BETA_Z += ALPHA[topic]*BETA*(one_over_z_n - one_over_z_o);

	            r = drand48();
	            Z_r_total = 0.0;
	            res2 = BETA_D_Z + ALPHA_BETA_Z;
	            success_sampled = false;
	            pivotal_topic_count = 0;
	            di_2_norm_temp = di_OVER_Z_2_norm + 1.0e-18;
                Z_r_total_o = 0.0;
                one_over_Z_e_o = 0.0;
                pivotal_topic_count_o = 0;
                wi_2_norm_temp = 0L;

               	res1 = sqrt(ALPHA_W_2_norm[wi]*ONE_OVER_Z_2_norm);

               	if (wi == wi_old){
               		wi_2_norm_temp = wi_2_norm_old;  //usewi_2_norm_old if doc and word are the same as the previous token
               	}else{
               		common_topic_count = 0;
               		for (di_nonzero_index=0; di_nonzero_index < di_pair_count; di_nonzero_index++){
          				itemp = wi_array[ dp_pair_idx[di_nonzero_index].index ];
           				if (itemp > 0){
           					wi_2_norm_temp += square_of[ itemp ];
          					common_topic_indexes[ common_topic_count ] = dp_pair_idx[di_nonzero_index].index;
           					common_topic_count++;
           				}
               		}
               		wi_2_norm_old = wi_2_norm_temp;
               	}
//10,5,10
               	check_threshold = common_topic_count / 10;
               	int common_topic_count_minus_one = common_topic_count - 1;
               	int check_threshold_increment = 5 + common_topic_count / 10;
               	for (di_nonzero_index=0; di_nonzero_index < common_topic_count; di_nonzero_index++){
               		topic_temp = common_topic_indexes[ di_nonzero_index ];
               		dtemp1 = di_array[topic_temp] * ONE_OVER_Z[ topic_temp ];
               		di_2_norm_temp -= dtemp1*dtemp1;
               		itemp1 = wi_array[topic_temp];
               		wi_2_norm_temp -= square_of[ itemp1 ];
       				Z_r_total += dtemp1* (double) itemp1;
       				pivotal_topic_indexes[ pivotal_topic_count ] = topic_temp;
       				pivotal_topic_values[ pivotal_topic_count ] = Z_r_total + Z_r_total_o;
       				pivotal_topic_count++;
               		if (di_nonzero_index >= check_threshold){
               			//check point reached
                		dtemp1 = di_2_norm_temp* (double) wi_2_norm_temp;

                		if (dtemp1 > 1.0e-18)
                  			Z_e = sqrt(dtemp1) + Z_r_total + Z_r_total_o + res1 + res2;
                  		else
                  			Z_e = Z_r_total + Z_r_total_o + res1 + res2;

                  		Z_r_total_n = Z_r_total + (1.0 - Z_e * one_over_Z_e_o) * Z_r_total_o;
                   		//r_n = (r - one_over_Z_e_o * Z_r_total_o);
                  		r_Z_e = (r - one_over_Z_e_o * Z_r_total_o) * Z_e;

                  		if (r_Z_e <= Z_r_total_n){
                  			success_sampled = true;
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

               	if (success_sampled) success_count1++;

	            if (!success_sampled){
	               	//if (Z_bound_count <=0){
	               	//w_pivotal_topic_count = 0;
	               	real_ALPHA_W_Z = 0.0;
	               	for (i=0; i<wi_pair_count; i++){
	               		itemp = wi_idx[i];
	               		real_ALPHA_W_Z += ALPHA[itemp] * wi_array[itemp] * ONE_OVER_Z[itemp];
	               		w_pivotal_topic_indexes[ i ] = itemp;
	               		w_pivotal_topic_values[ i ] = real_ALPHA_W_Z;
	               		//w_pivotal_topic_count++;
	               	}
	               	w_pivotal_topic_count = wi_pair_count;

	               	Z_e = (Z_r_total + Z_r_total_o + real_ALPHA_W_Z + res2);
	              	r_Z_e = r * Z_e;
	               	if (r_Z_e < (Z_r_total + Z_r_total_o)){
	               		//use the real_ALPHA_W_Z to check whether sampling is successful with the first term
	               		success_sampled = true;
	               		success_count2++;
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
	               		//check sampling with the 2nd term
	                   		//r_z_e is the actural Z value
	               		r_Z_e -= (Z_r_total + Z_r_total_o);
   						if (r_Z_e <= real_ALPHA_W_Z){
   							success_sampled = true;
   							success_count2++;
   	  						i_start = 0;
   	  						i_end = w_pivotal_topic_count - 1;
   	  						i_mid = (i_start + i_end) / 2;
   	  						if (r_Z_e < w_pivotal_topic_values[0]){
   	  							sampled_topic = w_pivotal_topic_indexes[0];
   	  						}else{
   	  							while (i_start < i_mid){
   	  								if (w_pivotal_topic_values[i_mid] <= r_Z_e)
   	  									i_start = i_mid;
   	  								else
   	  									i_end = i_mid;
   	  								i_mid = (i_start + i_end)/2;
   	  							}
   	  							sampled_topic = w_pivotal_topic_indexes[i_end];
   	  						}
   						}else
   							r_Z_e -= real_ALPHA_W_Z + 1.0e-18; //1.0e-18 is for rounding error
	                }
                }

	            if (!success_sampled){
	            	//check sampling with the 3rd and 4th term
	            	success_count3++;
	            	if (r_Z_e <= BETA_D_Z){
	   					r_Z_e /= BETA;
	   					for (i=0; i<di_pair_count; i++){
	   						sampled_topic = dp_pair_idx[i].index;
	   						r_Z_e -= di_array[sampled_topic]*ONE_OVER_Z[sampled_topic];
	   						if (r_Z_e <= 0.0 ){
	   							break;
	   						}
	   					}
	   				}else{
	   					r_Z_e -= BETA_D_Z;
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
	            if (wi_array[sampled_topic] == 0){
	               	wi_idx[ wi_pair_count ] = sampled_topic;
	                wi_pair_count++;

	               	common_topic_indexes[ common_topic_count ] = sampled_topic;
	               	common_topic_count++;
	               	common_topic_added = true;
	            }
	            wi_value_o = wi_array[sampled_topic];
	            wi_array[sampled_topic]++;

	            if (di_array[sampled_topic]==0){
	                //new topic
	               	dp_pair_idx[ di_pair_count ].index = sampled_topic;
	               	di_pair_count++;

	               	if (!common_topic_added){
	               		common_topic_indexes[ common_topic_count ] = sampled_topic;
	               		common_topic_count++;
	               		common_topic_added = true;
	               	}
	            }

	            if (common_topic_added)
	            	wi_2_norm_old += square_of[wi_value_o + 1];
	            else
	            	wi_2_norm_old += square_increment[ wi_value_o ];

	            di_value_o = di_array[sampled_topic];
	            di_array[sampled_topic]++;
	            di_value_n = di_value_o + 1;

	            one_over_z_o = ONE_OVER_Z[ sampled_topic ];
	            ztot[sampled_topic]++;
	            one_over_z_n = 1.0/(WBETA + (double) ztot[sampled_topic]);


	            dtemp1 = di_value_o*one_over_z_o;
	            dtemp2 = di_value_n*one_over_z_n;
	            di_OVER_Z_2_norm += dtemp2*dtemp2 - dtemp1*dtemp1;
	            ALPHA_W_2_norm[wi] += ALPHA_square[sampled_topic] * (double) square_increment[ wi_value_o ];
	            ONE_OVER_Z[ sampled_topic ] = one_over_z_n;
	            ONE_OVER_Z_2_norm += one_over_z_n*one_over_z_n - one_over_z_o*one_over_z_o;
	            BETA_D_Z += BETA*(di_value_n*one_over_z_n - di_value_o*one_over_z_o);
	            ALPHA_BETA_Z += ALPHA[sampled_topic]*BETA*(one_over_z_n - one_over_z_o);
	            W_pair_count[wi] = wi_pair_count;
	            wi_old = wi;
	            token_i++;
	        }
	        //end of a doc
	        for (j=0; j<di_pair_count; j++){
	        	dp_pair_val[j] = di_array[dp_pair_idx[j].index];
	        	di_array[dp_pair_idx[j].index] = 0;
	        }
	        D_pair_count[di] = di_pair_count;
	    }//end of an iter

	    time(&t_end_iteration);
	    printf("Total Gibbs iteration time used : %f\n", difftime(t_end_iteration, t_begin_iteration));
	    //printf("count1 = %d, count2 = %d, count3 = %d \n", success_count1, success_count2, success_count3);
	    eval_likelihood();
    }
	       //end of all iterations
}
//************************************************************************************************

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



int main(int argc, char **argv){
  	char line[ 80 ];
    int n, i, j, docId, wordId, count, max_doc, max_word;
    FILE *fr;

    parse_args(argc, argv);

    //square_limit = 10000000;
    square_limit = 0;
    printf("\nSFastLDA doc-by-doc: Opening and reading dataset file...\n");
    printf("*******************************************************\n");

    D = 0; W = 0; nTokens = 0; max_doc = 0; max_word = 0;
    fr = fopen(dFile,"rt");
    if (fr == NULL){
    	printf("\n\nCannot find the UCI Bag-of-words format input corpus file %s !!!\n\n\n", dFile);
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
      	if (docId > max_doc) max_doc = docId;
      	if (wordId > max_word) max_word = wordId;
      	if ((docId <= 0)||(wordId<=0)||(count<=0)){
      		printf("***********************************");
      		return -1;
      	}
      	if ((lineNo % 2000000) == 0)
			printf("reading input data (first pass) : processing line %d \n", lineNo);
    }
    printf("\nTotal documents: %d, Vocab Size : %d, Num. Tokens : %d\n", D, W, nTokens);
    printf("\nTotal documents: %d, Vocab Size : %d, Num. Tokens : %d\n", max_doc, max_word, nTokens);

    //allocate memory
    d = (int *) calloc(nTokens, sizeof(int));
    w = (int *) calloc(nTokens, sizeof(int));
    //victor z = (int *) calloc(nTokens, sizeof(int));
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
    	if ((lineNo % 2000000) == 0)
    		printf("reading input data (second pass)  : processing line %d \n", lineNo);
    }
    fclose(fr);

    for (i=0; i<W; i++){
    	if (word_lens[i] > square_limit)
    		square_limit = word_lens[i];
    }

    printf("Square limit: %d\n", square_limit);
    printf( "Arguments:\n" );
    printf( "\tNumber of words      W = %d\n"    , W );
    printf( "\tNumber of docs       D = %d\n"    , D );
    printf( "\tNumber of topics     T = %d\n"    , T );
    printf( "\tNumber of iterations N = %d\n"    , NN );
    printf( "\tHyperparameter   ALPHA = %4.4f\n" , argALPHA );
    printf( "\tHyperparameter    BETA = %4.4f\n" , BETA );
    printf( "\tSeed number            = %d\n"    , SEED );
    printf( "\tNumber of tokens       = %d\n"    , nTokens );
    if (oFile != NULL){
    	printf("\n\tOutput file: %s.\n\n", oFile);
    }else{
    	printf("\n\tNo output file specified.\n\n");
    }


    ztot = new int[T]; std::fill_n(ztot, T, 0);
    ALPHA = new double[T]; std::fill_n(ALPHA, T, argALPHA);
    ALPHA_W_2_norm = new double[W];	//summation over topics of ALPHA*wi
    D_pair_count = new int[D]; std::fill_n(D_pair_count, D, 0);
    W_pair_count = new int[W]; std::fill_n(W_pair_count, W, 0);
    pivotal_topic_indexes = new int[T];
    pivotal_topic_values = new double[T];
    w_pivotal_topic_indexes = new int[T];
    w_pivotal_topic_values = new double[T];
    common_topic_indexes = new int[T];
    ONE_OVER_Z = new double[T];
    square_increment = new long[square_limit];
    square_decrement = new long[square_limit];
    square_of = new long[square_limit];
    ALPHA_square = new double[T];
    di_array = new int[T];

    dpIdx = new DoublePair*[D];
    //victor dpVal = new int*[D];
    dpVal = new int*[D];
    for (i=0; i<D; i++){
    	if (doc_lens[i] < T){
    		dpIdx[i] = new DoublePair[doc_lens[i]];
    		//victor dpVal[i] = new int[doc_lens[i]];
    		dpVal[i] = new int[doc_lens[i]];
    		for (j=0; j<doc_lens[i]; j++){
    			dpIdx[i][j].index = -1;
    		}
    		std::fill_n(dpVal[i], doc_lens[i], 0);
    	}else{
    		dpIdx[i] = new DoublePair[T];
    		//victor dpVal[i] = new int[T];
    		dpVal[i] = new int[T];
    		for (j=0; j<T; j++){
    			dpIdx[i][j].index = -1;
    		}
			std::fill_n(dpVal[i], T, 0);
    	}
    }

    //victor wpIdx = new int*[W];
    wpIdx = new int*[W];
    wpVal = new int*[W];
    for (i=0; i<W; i++){
    	if (word_lens[i] < T){
    		//victor
    		wpIdx[i] = new int[word_lens[i]];
    		std::fill_n(wpIdx[i], word_lens[i], -1);
    	}else{
    		//victor
    		wpIdx[i] = new int[T];
    		std::fill_n(wpIdx[i], T, -1);
    	}
    	wpVal[i] = new int[T];
    	std::fill_n(wpVal[i], T, 0);
    }


    GibbsSamplerLDA();

    if (oFile != NULL){
    	export_results(oFile);
    	printf("\nOutput file is written to %s", oFile);
    }

    printf("\nFinished !!!\n");
}
