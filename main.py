class CSP:
    def __init__(self, n_var, domain_size,density,tightness):
        """
        Initialize the CSP problem.
        """
        self.variables = generate_variables(n_var)
        self.domains = generate_domains(self.variables, domain_size)
        self.density = density
        self.tightness = tightness
        self.constrains = generate_constrains(self.domains,self.density,self.tightness)

    def generate_variables(self,n_var):
        """
        generate list of variables
        """
        variables = []
        return variables
    def generate_domains(self,variables,domain_size):
        """
        Generate domains:
        """
        domains = {} # {0:[1,2,3], 1:[1,2,3]}
        return domains