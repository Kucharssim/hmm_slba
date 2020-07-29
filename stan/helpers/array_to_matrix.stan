  matrix array_to_matrix(vector[] tran_prob, int N_sts){
    matrix[N_sts, N_sts] tran_prob_mat;
    
    for(k in 1:N_sts){
      tran_prob_mat[k,] = to_row_vector(tran_prob[k]);
    }
    
    return tran_prob_mat;
  }