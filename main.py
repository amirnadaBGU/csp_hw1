import copy
import random
from itertools import combinations
class CSP:
    def __init__(self, n_var, domain_size,density,tightness):
        """
        Initialize the CSP problem.
        """
        self.density = density
        self.tightness = tightness
        self.variables = self.generate_variables(n_var)
        self.init_domains = self.generate_domains(self.variables, domain_size) # Initial domains - to not be changed
        self.domains = copy.deepcopy(self.init_domains) # dynamic domains
        self.constrains = self.generate_constraints(self.variables, self.domains, self.density, self.tightness)
        self.constrains_checks = 0

    def generate_variables(self,n_var):
        """
        generate list of variables.
        input:
            integer meaning number of total variables (int)
        output:
            list of variables: each variable is an integer from 0 to n-1 (list)
        """
        variables = list(range(n_var))
        return variables

    def generate_domains(self,variables,domain_size):
        """
        Generate domains:
        input:
            [list of variables (int), domain size (int)]
        output:
            domains dictionary: key - variable (int), value - list of possible values (list of int)
        """
        domains = {var: list(range(domain_size)) for var in variables} # {var_i:[val_a,val_b,val_c], var_j:[val_a,val_b,val_c]}
        return domains

    def generate_constraints(self,variables,domains,density,tightness):
        """
        Generate constraints based on density (p1) and tightness (p2).
        Function saves two-sided constraints:
        input:
            [variables (list of int), domain size (dict as described above), density (float), tightness (float)]
        output:
            { (var_i, var_j): set of forbidden pairs } (dict with key - tuple of var (int), value - set of forbidden pairs)
        """
        constraints = {}
        n = len(variables)
        # --- Part 1 Decide about constrained variables---
        # Max number of constraints possible:
        max_constraints = (n * (n - 1)) // 2
        # Define number of constraints:
        num_constraints = int(max_constraints * density)
        # Choose random pairs to be constrained:
        all_possible_pairs = list(combinations(variables, 2))
        chosen_pairs = random.sample(all_possible_pairs, num_constraints)

        # --- Part 2 Decide about constrained values---
        # Get domain size
        domain_values = domains[variables[0]]
        domain_size = len(domain_values)
        # All value pairs possible - cartesian multiplication:
        all_value_pairs = [(v1, v2) for v1 in domain_values for v2 in domain_values]
        # Get num of conflicts:
        max_conflicts = domain_size * domain_size
        num_conflicts = int(max_conflicts * tightness)

        for (var1, var2) in chosen_pairs:
            # Get Random pairs of values
            forbidden_pairs = set(random.sample(all_value_pairs, num_conflicts))

            # Save the original constrain
            constraints[(var1, var2)] = forbidden_pairs

            # Save the opposite constrain as well
            reversed_forbidden_pairs = {(v2, v1) for (v1, v2) in forbidden_pairs}
            constraints[(var2, var1)] = reversed_forbidden_pairs

        return constraints

    def check_conflict(self, var1, val1, var2, val2):
        """
        Check if two variables and their assigned values violate a constraint.
        input:
            var1, var2: The variable identifiers (integers).
            val1, val2: The values assigned to them.
        output:
            True if there is a conflict (violation), False otherwise.
        """
        # Add constraints check in 1:
        self.constraint_checks += 1
        # Check if a constraint exists between these two variables
        # (We only need to check one direction because we stored both directions)
        if (var1, var2) in self.constraints:
            forbidden_pairs = self.constraints[(var1, var2)]
            # Check if the specific value pair is in the forbidden set
            if (val1, val2) in forbidden_pairs:
                return True  # Conflict detected!
        return False  # No conflict



if __name__ == "__main__":
    mycsp = CSP(n_var=3,domain_size=3,density=0.34,tightness=0.33)