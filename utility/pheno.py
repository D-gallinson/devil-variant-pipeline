from samples import Samples, TumorSamples, Compare
import pandas as pd


# samples = TumorSamples(tissue="both", id_paths=Samples.prelim)
# samples.survival_proxy()
# samples.subset_non_nan("devil_survival_days")
# samples.normalize("devil_survival_days", method="log")
# samples.normalize("devil_survival_days")
# samples.to_factor("Sex")
# samples.ATOMM_pheno("../../results/ATOMM/input/master_pheno.txt", ["Sex", "devil_survival_days"])


samples = Samples(id_paths=Samples.batch_ids, assume_YOB=False)

# samples.estimate_age()
# samples.to_factor("Sex")
# samples.subset_non_nan("infection_age")
# samples.normalize("infection_age", method="log")
# samples.normalize("infection_age")
# samples.ATOMM_pheno("../../master_phenotype.txt", ["Sex", "infection_age"])


# samples.survival_proxy()
# samples.subset_non_nan("devil_survival_days")
# print(samples.count("Sex", "Male"))
# print(samples.count("Sex", "Female"))

# comp = Compare(samples.sample_df, ["devil_survival_days", "complementary_surv_estimate"])
# comp.compare_strict()
# comp.show_diff()
# comp.stats()