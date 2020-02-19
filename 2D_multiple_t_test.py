#!/usr/bin/python

'''
This script computes the statistical significance (Welch's t-test) between each point in a 2D matrix.
'''

import sys, os
from scipy import stats
import numpy as np
import string
#from statsmodels.stats import multitest
from statsmodels.sandbox.stats.multicomp import multipletests
#from multipy.fwer import sidak

def read_density_file(density_file):
    '''Read the density file'''

    FL_list = []
    TMD_list = []

    print "Reading density file..."
    for line in density_file:
        line = line.strip()
        line_list = line.split()
        
        #print line_list[7:10]

        FL_list.append(np.array(line_list[2:7], dtype=float))
        TMD_list.append(np.array(line_list[9:15], dtype=float))

    return FL_list, TMD_list

def t_test(FL, TMD, my_alpha):
    '''
    Perform t-test (Student's or Welch's), comparing the statistical significance of the difference
    between each point of 2 2D matrices. Also computes the corrected p-values with the Holm-Sidak method.
    '''

    ttest_p_value_list = []

    #print "Calculating p-values..."
    
    for lipids in zip(FL, TMD):
		# Student's t-test
		#t_value, p_value = stats.ttest_ind(lipids[0], lipids[1])

		# Welch's t-test
		t_value, p_value = stats.ttest_ind(lipids[0], lipids[1], equal_var = False)      
		
		if np.isnan(p_value) == True:
			p_value = 1
            
		ttest_p_value_list.append(p_value)

    print "Calculating corrected p-values..."
    #pvals_corr_hs, alphaSidak_hs, alphaBonf_hs = multipletests(ttest_p_value_list, alpha = my_alpha, method='hs')[1:] # FWER, Holm-Sidak correction
    pvals_corr_fdr_bh = multipletests(ttest_p_value_list, alpha = my_alpha, method='fdr_bh')[1]                       # FDR, Benjamini-Hochberg correction
    #pvals_corr_fdr_by = multipletests(ttest_p_value_list, alpha=my_alpha, method='fdr_by')[1]                         # FDR, Benjamini-Yekutieli correction
    #pvals_corr_fdr_tsbh = multipletests(ttest_p_value_list, alpha=my_alpha, method='fdr_tsbh')[1]                     # FDR, two stage fdr correction (non-negative)
    #pvals_corr_fdr_tsbky = multipletests(ttest_p_value_list, alpha=my_alpha, method='fdr_tsbky')[1]                   # FDR, two stage fdr correction (non-negative)
    
    ## Multipy implementation
    #significant_pvals = sidak(ttest_p_value_list, alpha = 0.05)
    
    #return ttest_p_value_list, pvals_corr_hs, pvals_corr_fdr_bh, pvals_corr_fdr_by, alphaSidak_hs, alphaBonf_hs
    return ttest_p_value_list, pvals_corr_fdr_bh

def write_values(p_values, pvals_corr_fdr_bh, my_alpha, out_name):
    '''Write the p-values for each pixel'''

    output_filename = 'p_values_' + out_name + '_' + str(my_alpha) + '_BH.dat'
    p_values_file = open(output_filename, 'w')

    p_values_file.write('#p-value' + '\t' + '#Benjamini-Hochberg' + '\n')

    print "Writing p-values for alpha = %s cut-off...\n" % str(my_alpha)

    for p in zip(p_values, pvals_corr_fdr_bh):
        p_values_file.write(str(p[0]) + '\t' + str(p[1]) + '\n')

    p_values_file.close()

# Main
density_file = open(sys.argv[1])
my_alpha = float(sys.argv[2])

out_name = os.path.splitext(os.path.basename(sys.argv[1]))[0][8:]

FL_list, TMD_list = read_density_file(density_file)

# Compute the p-values with the t-test, and the corrected p-values with FWER or FDR correction methods
#p_values, pvals_corr_hs, pvals_corr_fdr_bh, pvals_corr_fdr_by, alphaSidak_hs, alphaBonf_hs = t_test(FL_list, TMD_list, my_alpha)
p_values, pvals_corr_fdr_bh = t_test(FL_list, TMD_list, my_alpha)

# Write p-values
write_values(p_values, pvals_corr_fdr_bh, my_alpha, out_name)
#write_values(p_values, sig_pvalues)

print "Min p-value in the dataset \t \t \t \t \t:", np.amin(p_values)
#print "Min corrected p-value (Holm-Sidak correction, FWER) \t \t:", np.amin(pvals_corr_hs)
print "Min corrected p-value (Benjamini-Hochberg correction, FDR) \t:", np.amin(pvals_corr_fdr_bh)
#print "Min corrected p-value (Benjamini-Yekutieli correction, FDR) \t:", np.amin(pvals_corr_fdr_by)
#print "Min corrected p-value (two stage FDR correction 1) \t \t:", np.amin(pvals_corr_fdr_tsbh)
#print "Min corrected p-value (two stage FDR correction 2) \t \t:", np.amin(pvals_corr_fdr_tsbky)

#print "\nHolm-Sidak corrected alpha\t \t \t \t \t:", alphaSidak_hs
#print "Bonferroni corrected alpha\t \t \t \t \t:", alphaBonf_hs,"\n"

#if np.amin(pvals_corr_hs) < my_alpha:
    #print "*** At least one corrected p-value is significant after the Sidak's correction (FWER) *** \n"
if np.amin(pvals_corr_fdr_bh) < my_alpha:
    print "\n*** At least one corrected p-value is significant after the Benjamini-Hochberg correction (FDR) *** \n"
else:
    print "\n*** None of the corrected p-values is significant *** \n"
	
	
