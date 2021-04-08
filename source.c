/******************************************************************************** 
C implementation of four-parameter Least Squares sine fitting method proposed in 
"A Computationally Efficient Non-iterative Four-parameter Sine Fitting Method, 
IET Signal Processing, 2021

 Input parameters: 
              - x: Data samples                                    
              - w_init: initial relative angular frequency estimate    
              - dw: deviation from w0, where the three-parameter   
                    LS cost function is evaluated                  
              - N: number of samples    
              - pM_vect: parameter vector containing A, B, C and relative 
	                 angular frequency
			 
 Output parameters: pM_vect will contain the optimal parameter set after the 
                    evaluation
                                                                   
 Written by: Vilmos Pálfi                                        
                                                                   
 Last modified: April 8, 2021                                      

********************************************************************************/

void LS4p_parab(float32_t* x, float32_t w_init, float32_t dw, uint32_t N, float32_t* pM_vect)
{
	float32_t pM_vect_0[9], pM_vect_opt[3], CF_opt;
	
	float32_t dw_opt, w_opt;
	
	/* Matrix M contains the elements of the matrix in (20) of the articla
	   which is needed to determine optimal parameters                */
	
	float32_t M_vect[9] = { 0.5, -1.0, 0.5, -0.5, 0.0, 0.5, 0.0, 1.0, 0.0 };
	float32_t parab_vect[3], CF_vect[3]; // CF_vect contains 3 evaluations at 3 frequencies
	arm_matrix_instance_f32 M, parab, CF;
	arm_status status;
	
	/* Three LS3p fitting is evaluated at three different frequencies: w_init, w_init-dw, 
	   w_init+dw  */ 
	
	LS3p(x, w_init - dw, N, pM_vect_0, CF_vect);
	LS3p(x, w_init, N, pM_vect_0 + 3, CF_vect + 1);
	LS3p(x, w_init + dw, N, pM_vect_0 + 6, CF_vect + 2);
	
	/* Matrix initialisation in order to be able to perform matrix operation in the ARM 
	   environment   */
	
	arm_mat_init_f32(&M, 3, 3, M_vect);
	arm_mat_init_f32(&CF, 3, 1, CF_vect);
	arm_mat_init_f32(&parab, 3, 1, parab_vect);
	
	/* We perform matrix operation in (20) in in the article */
	
	status = arm_mat_mult_f32(&M, &CF, &parab);    // Check whether the operation was performed succesfully
	if ( status != ARM_MATH_SUCCESS)
	{
		print_error_func("Error: LS4p_parab, status = arm_mat_mult_f32(&M, &CF, &parab);");
	}
	dw_opt = -parab_vect[1]/2/parab_vect[0]*dw;  // We perform (21) in the article
	w_opt = w + dw_opt;       // Optimal relative angular frequency is calculated
	
	LS3p(x, w_opt, N, pM_vect_opt, NULL);  // Another LS3p method is performed to obtain optimal parameters at the optimal frequency
	pM_vect[0] = pM_vect_opt[0];
	pM_vect[1] = pM_vect_opt[1];
	pM_vect[2] = pM_vect_opt[2];
	pM_vect[3] = w_opt;
}
