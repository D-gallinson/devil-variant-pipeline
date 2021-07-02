# Read variables from a .env file. "vars" should be a list
class Env:
	env_path = "/work_bgfs/d/dgallinson/scripts/master"

	def __init__(self, env):
		self.vars = {}
		with open(f"{Env.env_path}/{env}") as h:
			for line in h:
				if line[0] == "#" or line == "\n":
					continue
				line = line.rstrip().split("=")
				self.vars[line[0]] = line[1]
		
	def get_var(self, var):
		return self.vars[var]