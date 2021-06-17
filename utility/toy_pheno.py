# Generate a toy phenotype.txt ATOMM file (16 samples for a 16x16 factor design)
import numpy as np

host = np.tile(np.arange(1,17), 16)
patho = np.tile(np.arange(1,17), (16,1)).T.flatten()
y_intercept = np.ones(len(host), dtype=int)
pheno = np.random.randint(2, size=len(host))
pheno_mat = np.column_stack((host, patho, y_intercept, pheno))
np.savetxt("/work_bgfs/d/dgallinson/outputs/intermediates/8_joint-variants_DEP/phenotype.txt", pheno_mat, fmt="%i", delimiter="\t")