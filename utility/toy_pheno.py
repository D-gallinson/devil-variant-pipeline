# Generate a toy phenotype.txt ATOMM file (16 samples for a 16x16 factor design)
import numpy as np

host = np.tile(np.arange(1,17), 16)
patho = np.tile(np.arange(1,17), (16,1)).T.flatten()
y_intercept = np.ones(len(host), dtype=int)
co1 = np.random.randint(2, size=len(host))
co2 = np.random.normal(0, 5, [len(host), 1])
pheno = np.random.normal(500, 10, [len(host), 1])
pheno_mat = np.column_stack((host, patho, y_intercept, co1, co2, pheno))
np.savetxt("phenotype.txt", pheno_mat, fmt="%f", delimiter="\t")