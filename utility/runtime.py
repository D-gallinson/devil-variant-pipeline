#Calculate basic execution time summary stats for a .err file formatted as: \nExecution Time: 1234s

import os
import math


dir = "scripts/align/align-logs"
files = os.listdir(dir)
files = [f"{dir}/{file}" for file in files if file.endswith(".err")]

times = []
for file in files:
	with open(file) as handle:
		handle.readline()	#take this out if the file doesn't start with \n
		time = handle.readline().rstrip().split()[-1][:-1]
		times.append(int(time))

total = sum(times)
mean = total / len(times)
std = math.sqrt(sum([(time - mean)**2 for time in times]) / total)

print("-----Mean Execution Times-----")
print(f"Total runs: {len(times)}")
print(f"{mean:.2f}s +/- {std:.2f}s")
print(f"{mean/60:.2f} min +/- {std/60:.2f} min")
print(f"{mean/60**2:.2f} hours +/- {std/60**2:.2f} hours")
print(f"{mean/(60**2 * 24):.2f} days +/- {std/(60**2 * 24):.2f} days")