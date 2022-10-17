functions {
  vector game_logl(vector params, vector mu_game, real[] time_arr, int[] marks_arr) {
    int L = size(time_arr);
    int M = 30;
    real alpha;
    matrix[M, 3] delta_mat;
    matrix[M, M] beta[3];
    matrix[M, M] gamma[3];
    matrix[M, M] beta_mat;
    matrix[M, M] gamma_mat;
    matrix[M, M-1] gamma_tr;
    matrix[M, M-1] gamma_inc;
    vector[M] gamma_col;
    int   i = 2;
    real ll = 0;
    int hist_vec[L] = marks_arr[1:L];
    int loc_vec[L] = marks_arr[(L+1):(2*L)];
    int mark_vec[L] = marks_arr[(2*L+1):(3*L)];
    vector[L] time_vec = to_vector(time_arr);
    vector[M] loglik;
      
    for(loc in 1:3){
      
      beta[loc] = to_matrix(params[(1+(loc-1)*M*M):(loc*M*M)], M, M);
      gamma_tr  = to_matrix(params[(3*M*M+1+(loc-1)*M*(M-1)):(3*M*M+loc*M*(M-1))], M, M-1) + 
        rep_matrix(to_row_vector(mu_game), M);
      delta_mat[, loc] = params[(3*M*M+3*M*(M-1)+(loc-1)*M+1):(3*M*M+3*M*(M-1)+(loc*M))];
  
      for(j in 1:M){
        gamma_inc[j] = exp(gamma_tr[j])/(1 + sum(exp(gamma_tr[j])));
        gamma_col[j] = 1 - sum(gamma_inc[j]);
      }
      
      gamma[loc] = append_col(gamma_inc, gamma_col);
    }
    alpha = params[(3*M*M+3*M*(M-1)+(3*M)+1)];
    
    while( (i <= L) && (mark_vec[i] <= M) ){ // likelihood for marks
        
        beta_mat = beta[loc_vec[i]];
        gamma_mat = gamma[loc_vec[i]];
        
        for(k in 1:M){
          loglik[k] = log_sum_exp( append_row( log(delta_mat[k, loc_vec[i]]), 
              alpha - beta_mat[mark_vec[(i-hist_vec[i]):(i-1)], k] .* (time_vec[i] - time_vec[(i-hist_vec[i]):(i-1)] ) + 
              log(gamma_mat[mark_vec[(i-hist_vec[i]):(i-1)], k]) ) );
        }
        ll += loglik[mark_vec[i]] - log_sum_exp(loglik);
      
        i += 1;
    }
    
    return [ll]';
  }
}

data{
  int<lower=1> G; // No. of games
  int<lower=1> L; // Max game length
  int<lower=1> M; // Number of marks
  int<lower=1> T; // Number of teams
  int<lower=1> B1; // Number of beta pars
  int<lower=1> indx_def[B1, 2]; // Beta indices
  int<lower=1> indx_mid[B1, 2]; // Beta indices
  int<lower=1> indx_atk[B1, 2]; // Beta indices
  int<lower=1> meta[G, 3]; // Team info
  real<lower=0> times[G, L]; // observed times 
  int<lower=0> marks[G, 3*L]; // observed marks and states
}

parameters{
  vector<lower=0>[B1] beta_def_raw; // unknown beta
  vector<lower=0>[B1] beta_mid_raw; // unknown beta
  vector<lower=0>[B1] beta_atk_raw; // unknown beta
  vector[2*B1] theta_def_raw; // unknown theta
  vector[2*B1] theta_mid_raw; // unknown theta
  vector[2*B1] theta_atk_raw; // unknown theta
  simplex[M] delta[3]; // unknown delta
  real alpha; // unknown alpha
  matrix[T-1, M-1] mu_atk; // unknown mu_atk
}

model{
  matrix[M, M] beta_zone_1 = rep_matrix(100, M, M);
  matrix[M, M] beta_zone_2 = rep_matrix(100, M, M);
  matrix[M, M] beta_zone_3 = rep_matrix(100, M, M);
  matrix[M, M-1] theta_zone_1 = rep_matrix(-100, M, M-1);
  matrix[M, M-1] theta_zone_2 = rep_matrix(-100, M, M-1);
  matrix[M, M-1] theta_zone_3 = rep_matrix(-100, M, M-1);
  vector[3*M*M + 3*M*(M-1) + 3*M + 1] params;
  matrix[M, 3] delta_mat;
  matrix[T, M-1] mu_atk_full = append_row(mu_atk, rep_row_vector(0, M-1));
  vector[M-1] mu_game[G];
  
  for (i in 1:B1){
    beta_zone_1[indx_def[i, 1], indx_def[i, 2]] = beta_def_raw[i];
    theta_zone_1[indx_def[i, 1], indx_def[i, 2]] = theta_def_raw[i];
    if(indx_def[i, 1] < (M/2 + 1)){
      beta_zone_3[indx_def[i, 1]+(M/2), indx_def[i, 2]+(M/2)] = beta_def_raw[i];
      theta_zone_3[indx_def[i, 1]+(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i];
    }else{
      beta_zone_3[indx_def[i, 1]-(M/2), indx_def[i, 2]+(M/2)] = beta_def_raw[i];
      theta_zone_3[indx_def[i, 1]-(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i];
    }
  
    beta_zone_2[indx_mid[i, 1], indx_mid[i, 2]] = beta_mid_raw[i];
    theta_zone_2[indx_mid[i, 1], indx_mid[i, 2]] = theta_mid_raw[i];
    if(indx_mid[i, 1] < (M/2 + 1)){
      beta_zone_2[indx_mid[i, 1]+(M/2), indx_mid[i, 2]+(M/2)] = beta_mid_raw[i];
      theta_zone_2[indx_mid[i, 1]+(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i];
    }else{
      beta_zone_2[indx_mid[i, 1]-(M/2), indx_mid[i, 2]+(M/2)] = beta_mid_raw[i];
      theta_zone_2[indx_mid[i, 1]-(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i];
    }

    beta_zone_3[indx_atk[i, 1], indx_atk[i, 2]] = beta_atk_raw[i];
    theta_zone_3[indx_atk[i, 1], indx_atk[i, 2]] = theta_atk_raw[i];
    if(indx_atk[i, 1] < (M/2 + 1)){
      beta_zone_1[indx_atk[i, 1]+(M/2), indx_atk[i, 2]+(M/2)] = beta_atk_raw[i];
      theta_zone_1[indx_atk[i, 1]+(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i];
    }else{
      beta_zone_1[indx_atk[i, 1]-(M/2), indx_atk[i, 2]+(M/2)] = beta_atk_raw[i];
      theta_zone_1[indx_atk[i, 1]-(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i];
    }
  }
  
  // Priors
  beta_def_raw ~ exponential(2);
  beta_mid_raw ~ exponential(2);
  beta_atk_raw ~ exponential(2);
  theta_def_raw ~ std_normal();
  theta_mid_raw ~ std_normal();
  theta_atk_raw ~ std_normal();
  alpha ~ std_normal();
  to_vector(mu_atk) ~ std_normal();
  
  for (j in 1:3) {
    delta[j] ~ dirichlet(rep_vector(1.0, M));
    delta_mat[, j] = delta[j];
  }
  
  params =  append_row( append_row( append_row( append_row( append_row( append_row( append_row( 
  to_vector(beta_zone_1) ,
  to_vector(beta_zone_2) ),
  to_vector(beta_zone_3) ),
  to_vector(theta_zone_1) ),
  to_vector(theta_zone_2) ),
  to_vector(theta_zone_3) ),
  to_vector(delta_mat) ),
  alpha );
  
  for(i in 1:G){
   mu_game[i] = to_vector(append_col( mu_atk_full[meta[i,2],1:(M/2)],  mu_atk_full[meta[i,3],(M/2 + 1):(M-1)]));
  }
  
  target += sum(map_rect(game_logl, params, mu_game, times, marks));
}
