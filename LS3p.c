/******************************************************************************** 
C implementation of three-parameter Least Squares sine fitting method

Supplementary file for LS4p_parab.c

 Input parameters: 
              - x: Data samples                                    
              - w: relative angular frequency estimate                          
              - N: number of samples    
              - pM_vect: parameter vector containing A, B, C
	      - CF: cost function
			 
 Output parameters:
               -  pM_vect will contain the optimal parameter set after the 
                      evaluation
	       -  CF will contain the cost function after the evaluation
		      
                                                                   
 Written by: Vilmos Pálfi                                        
                                                                   
 Last modified: April 8, 2021                                      
********************************************************************************/

void LS3p(float32_t* x, float32_t w, uint32_t N, float32_t* pM_vect, float32_t* CF)
{
	float32_t D_matrix_vect[3*N], DT_vect[3*N], DTD_vect[9], DTD_inv_vect[9], DTx_vect[3];
	float32_t yfit_vect[N], e_vect[N];
	arm_matrix_instance_f32 D, DT, DTD, DTD_inv, xM, DTx, pM, yfit;
	arm_status status;
	uint32_t n;
	float32_t c;

	/* First, columns of matrix D is created, that is, columns corresponding to parameter A, B and C */
	
	for (n = 0; n < N; n++)
	{
		
		D_matrix_vect[3*n] = arm_cos_f32((n+1)*w);
		D_matrix_vect[3*n + 1] = arm_sin_f32((n+1)*w);
		D_matrix_vect[3*n + 2] = 1.0;
	}
	
	/* Matrix initialisation in order to be able to perform matrix operation in the ARM 
	   environment   */
	
	arm_mat_init_f32(&D, N, 3, D_matrix_vect);  // From vector D_matrix_vect we get matrix D with reshaping
	arm_mat_init_f32(&DT, 3, N, DT_vect);
	arm_mat_init_f32(&DTD, 3, 3, DTD_vect);
	arm_mat_init_f32(&DTD_inv, 3, 3, DTD_inv_vect);
	arm_mat_init_f32(&xM, N, 1, x);
	arm_mat_init_f32(&DTx, 3, 1, DTx_vect);
	arm_mat_init_f32(&pM, 3, 1, pM_vect);
	arm_mat_init_f32(&yfit, N, 1, yfit_vect);
	status = arm_mat_trans_f32(&D, &DT);     // In matrix DT we have the transpose of matrix D
	if ( status != ARM_MATH_SUCCESS)
	{
		print_error_func("Error: LS3p, status = arm_mat_trans_f32(&D, &DT);"); //We check whether the operation was performed without errors
	}
	status = arm_mat_mult_f32(&DT, &D, &DTD);  // In DTD we have D'*D (' denotes the transpose operator)
	if ( status != ARM_MATH_SUCCESS)
	{
		print_error_func("Error: LS3p, status = arm_mat_mult_f32(&DT, &D, &DTD);");
	}
	status = arm_mat_inverse_f32(&DTD, &DTD_inv); // In DTD_inv we have (D'*D)^-1
	if ( status != ARM_MATH_SUCCESS)
	{
		print_error_func("Error: LS3p, status = arm_mat_inverse_f32(&DTD, &DTD_inv);");
	}
	status = arm_mat_mult_f32(&DT, &xM, &DTx); // In DTx we have D'*x
	if ( status != ARM_MATH_SUCCESS)
	{
		print_error_func("Error: LS3p, status = arm_mat_mult_f32(&DT, &xM, &DTx);");
	}
	status = arm_mat_mult_f32(&DTD_inv, &DTx, &pM); // We calculate the parameter vector: pM = (D'*D)^-1 *(D'*x), see (9) in the article
	if ( status != ARM_MATH_SUCCESS)
	{
		print_error_func("Error: LS3p, status = arm_mat_mult_f32(&DTD_inv, &DTx, &pM);");
	}

	if (CF)   // If CF is not NULL, we evaluate the cost function
	{
		status = arm_mat_mult_f32(&D, &pM, &yfit);  // The fitted vector yfit is given by D*pM, see (24) in the article
		if ( status != ARM_MATH_SUCCESS)
 		{
			print_error_func("Error: LS3p_mod, status = arm_mat_mult_f32(&DTD_inv, &DTx, &pM);");
		}

		arm_sub_f32(x, yfit_vect, e_vect, N);  //e_vect contains the difference between the fitted sine and the original samples (error vector)

		arm_dot_prod_f32(e_vect, e_vect, N, CF); // Dot product is calculated (the elements are multiplied one-by-one and they are summed). This is the 
		                                         // Cost function, see (25) in the article
	}
}
