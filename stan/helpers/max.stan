  int max_int(int[,] x){
    int N_sts = size(x);
    int out[N_sts];
    // int cont = 0;
    
    for(k in 1:N_sts) out[k] = max(x[k]);
    
    // while(cont < max(out)) cont += 1;
    return max(out);
  }