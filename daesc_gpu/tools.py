# Define a beta-binomial functios
def beta_binomial_log_likelihood(params, y, n):
    alpha = np.exp(params[0])
    beta = np.exp(params[1])
    epsilon = 1e-10 
    log_lik = betaln(y + alpha + epsilon, n - y + beta + epsilon) - betaln(alpha + epsilon, beta + epsilon)

    return -np.sum(log_lik)

# Optimize the beta-binomial functions
def optimize_multiple_starts(y, n, num_starts=10):
    best_phi = None
    best_result = None

    for _ in range(num_starts):
      params_init = [np.log(np.random.uniform(0.5, 5.0)), np.log(np.random.uniform(0.5, 5.0))]

      result = minimize(beta_binomial_log_likelihood, 
                        params_init,args=(y, n),
                        method='Nelder-Mead', options={'maxiter': 1000, 'xatol': 1e-6, 'fatol': 1e-6} )

      if best_result is None or result.fun < best_result.fun:
        best_result = result
        alpha_hat = np.exp(result.x[0])
        beta_hat = np.exp(result.x[1])
        best_phi = 1 / (1 + alpha_hat + beta_hat)
    if best_phi is not None:
        phi = 1 / (1 / best_phi - 1)  # 根据 R 的公式
        phi = max(phi, 0.05)         # 确保 phi >= 0.05
    else:
        phi = 0.05
    return phi

# Some useful functions
# Clear used memroy in GPU
def clear_memory():
  memory_pool = cp.get_default_memory_pool()
  total_bytes_before = memory_pool.total_bytes()
  used_bytes_before = memory_pool.used_bytes()

  total_gb_before = total_bytes_before / (1024**3)
  used_gb_before = used_bytes_before / (1024**3)

  print()
  print(f'Before resetting the memory pool:')
  print(f'Total bytes currently allocated in the memory pool: {total_gb_before:.6f} GB')
  print(f'Total bytes currently used in the memory pool: {used_gb_before:.6f} GB')
  
  for obj in gc.get_objects():
      if isinstance(obj, cp.ndarray):
          del obj
  gc.collect()
  cp.get_default_memory_pool().free_all_blocks()
  cp.get_default_pinned_memory_pool().free_all_blocks()
  total_bytes_after = memory_pool.total_bytes()
  used_bytes_after = memory_pool.used_bytes()

  total_gb_after = total_bytes_after / (1024**3)
  used_gb_after = used_bytes_after / (1024**3)

  print(f'After resetting the memory pool:')
  print(f'Total bytes currently allocated in the memory pool: {total_gb_after:.6f} GB')
  print(f'Total bytes currently used in the memory pool: {used_gb_after:.6f} GB \n \n ')

# split data into n parts with equal number
def split_array(array, n):
    avg = len(array) / float(n)
    out = []
    last = 0.0
    while last < len(array):
        out.append(array[int(last):int(last + avg)])
        last += avg
    return out
