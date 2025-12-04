import copy
import random
from itertools import combinations
from typing import Optional

from matplotlib import pyplot as plt
from zope.interface.common import optional


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
        self.constraints = self.generate_constraints(self.variables, self.domains, self.density, self.tightness)
        self.constraint_checks = 0

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

    def check_assigment_validity(self, assignment):
        """
        Checks if the entire assignment is valid by verifying that no two
        assigned variables violate a constraint.
        input:
            assignment: dict {var: val}
        output:
            True if assignment is consistent, False otherwise.
        """
        # Check all combinations in the assigment
        for var1, var2 in combinations(assignment.keys(), 2):
            val1 = assignment[var1]
            val2 = assignment[var2]

            if self.check_conflict(var1, val1, var2, val2):
                return False

        return True

    def reset_metrics(self):

        self.constraint_checks = 0

# Algorithms:
def backtracking_search(csp):
    """
    Main entry point for Backtracking Search.
    input:
        csp: CSP instance
    output:
        dict {var: val} if solution found, else None.
    """
    return backtrack_recursive_call({}, csp)

def backtrack_recursive_call(assignment, csp):
    """
    Recursive function for Backtracking - Classic Implementation.
    input:
        assignment: Dictionary containing current assignments {0: 1, 1: 0...}
        csp: The CSP problem object
    output:
        dictionary containing current solution {0: 1, 1: 0...} or none
    """

    # 1. Base Case: Have we assigned values to all variables?
    if len(assignment) == len(csp.variables):
        return assignment

    # 2. Select Unassigned Variable
    # Pick the next variable in order.
    # Since variables are indices (0, 1, 2...), the next variable is the current length of assignment.
    var = csp.variables[len(assignment)]

    # 3. Domain Iteration: Try every value in the variable's domain
    for value in csp.domains[var]:

        # 4. Consistency Check
        is_consistent = True

        # Iterate over variables that have already been assigned (the past)
        for prev_var, prev_val in assignment.items():
            # Check for conflict between the new variable (var) and the old one (prev_var)
            # Note: The check_conflict method automatically increments csp.constraint_checks
            if csp.check_conflict(var, value, prev_var, prev_val):
                is_consistent = False
                break  # One conflict is enough to invalidate this value

        # 5. If the value is consistent, proceed
        if is_consistent:
            assignment[var] = value  # Assign the value

            # Recursive call to the next step
            result = backtrack_recursive_call(assignment, csp)

            # If the recursion returned a result (not None), a solution was found!
            if result is not None:
                return result

            # 6. Backtrack: If we reached here, the path failed.
            # Remove the assignment and try the next value in the loop.
            del assignment[var]

    # 7. If we tried all values and none worked - return failure to the previous level
    return None

class FC_Base_Solver:
    def __init__(self, csp):
        self.csp = csp
        self.assignment = {}
        self.domains = copy.deepcopy(csp.init_domains)
        self.backtracks = 0
        self.conflict_set = {v: set() for v in csp.variables}

    def _get_removed_values(self, assigned_var: int, assigned_val: int, var_j: int, domain_j: list[int]) -> list[int]:
        removed_vals = []
        if (assigned_var, var_j) in self.csp.constraints:
            for val_j in domain_j:
                if self.csp.check_conflict(assigned_var, assigned_val, var_j, val_j):
                    removed_vals.append(val_j)
        return removed_vals

    def _check_domain_wipeout(self, var_j: int, domain_j: list[int], assigned_var: int) -> tuple[Optional[int], set[int]]:
        if not domain_j:
            return var_j, {assigned_var}
        return None, set()

    def select_unassigned_variable(self):
        """
        Returns next unassigned variable.
        """
        for var in self.csp.variables:
            if var not in self.assignment:
                return var
        return None

    def order_domain_values(self, var):
        """
        Returns the domain ov a variable
        """
        return self.domains[var]

    def forward_check(self, assigned_var, assigned_val, current_domains, assignment):
        """
        Do FC
        """
        # Copy current domains
        new_domains = {v: d.copy() for v, d in current_domains.items()}

        # Going via all future variables to be assigned
        for var_j in self.csp.variables:
            if var_j not in assignment and var_j != assigned_var:

                # 1. קבלת ערכים להסרה מהמשתנה var_j
                removed_vals = self._get_removed_values(assigned_var, assigned_val, var_j, new_domains[var_j])

                if removed_vals:
                    # 2. הסרה
                    for val in removed_vals:
                        new_domains[var_j].remove(val)

                    # 3. בדיקת ריקון
                    failed_var, conflict_fc = self._check_domain_wipeout(var_j, new_domains[var_j], assigned_var)

                    if failed_var is not None:
                        return new_domains, failed_var, conflict_fc

        return new_domains, None, set()

    def solve(self):
        return self._search()

    def _search(self):
        if len(self.assignment) == len(self.csp.variables):
            return self.assignment

        var_i = self.select_unassigned_variable()
        if var_i is None: return self.assignment

        old_domains = copy.deepcopy(self.domains)

        for val_i in self.order_domain_values(var_i):
            self.assignment[var_i] = val_i

            new_domains, failed_var_j, _ = self.forward_check(var_i, val_i, self.domains, self.assignment)

            if failed_var_j is None:
                self.domains = new_domains
                result = self._search()
                if result is not None:
                    return result

            # Backtrack
            del self.assignment[var_i]
            self.backtracks += 1
            self.domains = old_domains

        return None

class FC_CBJ_Solver(FC_Base_Solver):
    def __init__(self, csp, verbose=False):
        super().__init__(csp)
        self.verbose = verbose
        self.active_jump_target = None
        self.conflicting_ancestors: dict[int, set[int]] = {v: set() for v in csp.variables}

    def log(self, msg):
        if getattr(self, 'verbose', False):
            indent = "  " * len(self.assignment)
            print(f"{indent}{msg}")

    def find_backjump_var(self, conflict_set):
        if not conflict_set: return -1
        max_var = -1
        for var in conflict_set:
            # מוודאים ש:
            # 1. המשתנה הוא חלק מההשמה הנוכחית
            # 2. הוא הגדול ביותר שמצאנו
            # 3. הוא לא המשתנה שאנחנו נמצאים בו כרגע (למרות שטיפלנו בזה בהסרה)
            if var in self.assignment and var > max_var:
                max_var = var
        return max_var

    def forward_check(self, assigned_var, assigned_val, current_domains, assignment):
        new_domains = {v: d.copy() for v, d in current_domains.items()}

        for var_j in self.csp.variables:
            if var_j not in assignment and var_j != assigned_var:

                # שמירת היסטוריה
                if (assigned_var, var_j) in self.csp.constraints:
                    self.conflicting_ancestors[var_j].add(assigned_var)

                removed_vals = self._get_removed_values(assigned_var, assigned_val, var_j, new_domains[var_j])

                if removed_vals:
                    for val in removed_vals:
                        new_domains[var_j].remove(val)

                    # DWO Detected
                    if not new_domains[var_j]:
                        # מחזירים את האבות של המשתנה שהתרוקן + המשתנה הנוכחי
                        culprits = self.conflicting_ancestors[var_j].copy()
                        culprits.add(assigned_var)
                        return new_domains, var_j, culprits

        return new_domains, None, set()

    def _cleanup_ancestors(self, variable_to_remove):
        for v in self.csp.variables:
            self.conflicting_ancestors[v].discard(variable_to_remove)

    def _search(self):
        if len(self.assignment) == len(self.csp.variables):
            self.log(f"✅ Solution Found: {self.assignment}")
            return self.assignment

        var_i = self.select_unassigned_variable()
        if var_i is None: return self.assignment

        old_domains = copy.deepcopy(self.domains)

        # 1. אתחול עם ה-Ancestors (מהעבר)
        self.conflict_set[var_i] = self.conflicting_ancestors[var_i].copy()

        for val_i in self.order_domain_values(var_i):
            self.assignment[var_i] = val_i

            new_domains, wiped_var, culprits = self.forward_check(var_i, val_i, self.domains, self.assignment)

            if wiped_var is None:
                # FC עבר בהצלחה
                self.domains = new_domains
                result = self._search()

                if result is not None:
                    return result

                # חזרה מקפיצה
                if self.active_jump_target is not None:
                    if self.active_jump_target < var_i:
                        self._cleanup_ancestors(var_i)
                        del self.assignment[var_i]
                        self.domains = old_domains
                        return None
                    else:
                        self.active_jump_target = None
            else:
                # FC נכשל (DWO)
                self.backtracks += 1
                culprits.discard(var_i)
                self.conflict_set[var_i].update(culprits)

            # Backtrack
            self._cleanup_ancestors(var_i)
            del self.assignment[var_i]
            self.domains = old_domains


        # Dead End
        var_k = self.find_backjump_var(self.conflict_set[var_i])

        if var_k == -1:
            return None

        # מורישים את הקונפליקטים למשתנה היעד
        if var_k in self.conflict_set: # שורת הגנה
            # לפני שמעבירים, מוודאים שוב ש-var_i לא נמצא שם (ליתר ביטחון)
            self.conflict_set[var_i].discard(var_i)
            self.conflict_set[var_k].update(self.conflict_set[var_i])

        self.active_jump_target = var_k
        return None


# Helper functions
def generate_graph(algorithm,p1,n_var,domain_size,iterations = 10):
    """
    Function that generates graphs according to requirements
    input:
        algorithm: algorithm name (string) # backtrack, forward_check, forward_check_cbj
        p1: number between 0 and 1 (float)
        n_var: number of variables (int)
        domain_size: number of variables (int)
    output:
        void - function generates graphs according to requirements
    """
    # Generate runs
    constraint_checks = {}
    for i in range(1, 10):
        p2 = i / 10
        constraint_checks[p2]=[]
        for k in range(iterations):
            csp = CSP(n_var=n_var,domain_size=domain_size,density=p1,tightness=p2)
            if algorithm == 'backtrack':
                _ = backtracking_search(csp)
            else:
                print(f"Algorithm: {algorithm} is still not implemented")
                return
            constraint_checks[p2].append(csp.constraint_checks)

    # Generate data for graphs
    constraint_checks_averages = {}
    for key in constraint_checks:
        constraint_checks_averages[key] = sum(constraint_checks[key])/len(constraint_checks[key])

    # Generate matplotlib graph

    x_values = constraint_checks_averages.keys()
    y_values = [constraint_checks_averages[x] for x in x_values]

    # Graph creation
    plt.figure(figsize=(10, 6))  # canvas size
    plt.plot(x_values, y_values, marker='o', linestyle='-', color='b', label=algorithm)

    # Titles
    plt.title(f'CSP Performance: {algorithm}\n(N={n_var}, D={domain_size}, Density={p1})')
    plt.xlabel('Tightness (p2)')
    plt.ylabel('Average Constraint Checks')

    # Design
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.show()

def compare_algorithms(p1, n_var, domain_size, iterations=5):
    tightness_values = [i / 10.0 for i in range(1, 10)]

    avg_checks_fc = []
    avg_checks_cbj = []

    print(f"\n📈 Starting Graph Generation (N={n_var}, D={domain_size}, Density={p1}, Iterations={iterations})...")
    print("   Please wait, this might take a moment...")

    for p2 in tightness_values:
        total_checks_fc = 0
        total_checks_cbj = 0

        for k in range(iterations):
            csp = CSP(n_var=n_var, domain_size=domain_size, density=p1, tightness=p2)

            # FC Run
            csp.reset_metrics()
            solver_fc = FC_Base_Solver(csp)
            solver_fc.solve()
            total_checks_fc += csp.constraint_checks

            # CBJ Run
            csp.reset_metrics()
            solver_cbj = FC_CBJ_Solver(csp)
            solver_cbj.solve()
            total_checks_cbj += csp.constraint_checks

        avg_fc = total_checks_fc / iterations
        avg_cbj = total_checks_cbj / iterations
        avg_checks_fc.append(avg_fc)
        avg_checks_cbj.append(avg_cbj)
        print(f"   > Tightness {p2}: FC={avg_fc:.1f}, CBJ={avg_cbj:.1f}")

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(tightness_values, avg_checks_fc, marker='o', linestyle='-', color='blue', label='Forward Checking (FC)')
    plt.plot(tightness_values, avg_checks_cbj, marker='s', linestyle='--', color='red', label='FC + CBJ')

    plt.title(f'Performance Comparison: FC vs FC+CBJ\n(N={n_var}, D={domain_size}, Density={p1})')
    plt.xlabel('Tightness (p2)')
    plt.ylabel('Average Constraint Checks')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)

    print("📊 Graph generated successfully.")
    plt.show() #

if __name__ == "__main__":
    # mycsp = CSP(n_var=2,domain_size=2,density=0.5,tightness=0.5)
    # solution = backtracking_search(mycsp)
    # print(f"Solution: {solution} Correctness is {mycsp.check_assigment_validity(solution)}")
    # print(f"Number of constraints checks: {mycsp.constraint_checks}")
    # generate_graph('backtrack',0.4,10,10)
    # print("--- Sanity Check ---")
    # small_csp = CSP(n_var=10, domain_size=10, density=0.7, tightness=0.5)
    # solver = FC_CBJ_Solver(small_csp)
    # sol = solver.solve()
    compare_algorithms(0.3, 10, 10, iterations=100)