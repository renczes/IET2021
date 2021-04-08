// Three-parameters least squares fit

void LS3p(float32_t* x, float32_t w, uint32_t N, float32_t* pM_vect, float32_t* CF)
{
	float32_t X_matrix_vect[3*SAMPLES/2], DT_vect[3*SAMPLES/2], DTD_vect[9], DTD_inv_vect[9], DTx_vect[3];
	float32_t yfit_vect[SAMPLES/2], e_vect[SAMPLES/2];
	arm_matrix_instance_f32 D, DT, DTD, DTD_inv, xM, DTx, pM, yfit;
	arm_status status;
	uint32_t n;
	float32_t c;

	for (n = 0; n < SAMPLES/2; n++)
	{
		
		X_matrix_vect[3*n] = arm_cos_f32((n+1)*w);
		X_matrix_vect[3*n + 1] = arm_sin_f32((n+1)*w);
		X_matrix_vect[3*n + 2] = 1.0;
	}
	arm_mat_init_f32(&D, SAMPLES/2, 3, X_matrix_vect);
	arm_mat_init_f32(&DT, 3, SAMPLES/2, DT_vect);
	arm_mat_init_f32(&DTD, 3, 3, DTD_vect);
	arm_mat_init_f32(&DTD_inv, 3, 3, DTD_inv_vect);
	arm_mat_init_f32(&xM, SAMPLES/2, 1, x);
	arm_mat_init_f32(&DTx, 3, 1, DTx_vect);
	arm_mat_init_f32(&pM, 3, 1, pM_vect);
	arm_mat_init_f32(&yfit, SAMPLES/2, 1, yfit_vect);
	status = arm_mat_trans_f32(&D, &DT);
	if ( status != ARM_MATH_SUCCESS)
	{
		print_error_func("Error: LS3p, status = arm_mat_trans_f32(&D, &DT);");
	}
	status = arm_mat_mult_f32(&DT, &D, &DTD);
	if ( status != ARM_MATH_SUCCESS)
	{
		print_error_func("Error: LS3p, status = arm_mat_mult_f32(&DT, &D, &DTD);");
	}
	status = arm_mat_inverse_f32(&DTD, &DTD_inv);
	if ( status != ARM_MATH_SUCCESS)
	{
		print_error_func("Error: LS3p, status = arm_mat_inverse_f32(&DTD, &DTD_inv);");
	}
	status = arm_mat_mult_f32(&DT, &xM, &DTx);
	if ( status != ARM_MATH_SUCCESS)
	{
		print_error_func("Error: LS3p, status = arm_mat_mult_f32(&DT, &xM, &DTx);");
	}
	status = arm_mat_mult_f32(&DTD_inv, &DTx, &pM);
	if ( status != ARM_MATH_SUCCESS)
	{
		print_error_func("Error: LS3p, status = arm_mat_mult_f32(&DTD_inv, &DTx, &pM);");
	}

	if (CF)
	{
		status = arm_mat_mult_f32(&D, &pM, &yfit);
		if ( status != ARM_MATH_SUCCESS)
		{
			print_error_func("Error: LS3p_mod, status = arm_mat_mult_f32(&DTD_inv, &DTx, &pM);");
		}

		arm_sub_f32(x, yfit_vect, e_vect, SAMPLES/2);

		arm_dot_prod_f32(e_vect, e_vect, SAMPLES/2, CF);
	}
}
