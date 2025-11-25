import copy
import random
from itertools import combinations
from typing import Optional, Set, List, Tuple


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

    def reset_metrics(self):
        """
        פונקציית עזר לאיפוס המדדים לפני הרצה של אלגוריתם חדש על אותו מופע
        """
        self.constraint_checks = 0


class FC_Base_Solver:
    def __init__(self, csp):
        self.csp = csp
        self.assignment = {}
        # הדומיינים הדינמיים (משתנים בהתאם ל-FC)
        self.domains = copy.deepcopy(csp.init_domains)
        self.backtracks = 0
        # ה-Conflict Set קיים בבסיס כמידע, אך לא משומש ב-FC הרגיל
        self.conflict_set = {v: set() for v in csp.variables}

        # --- פונקציות ליבה משותפות ---

    def _get_removed_values(self, assigned_var: int, assigned_val: int, var_j: int, domain_j: List[int]) -> List[int]:
        """
        מחזירה את רשימת הערכים שיש להסיר מ-var_j עקב ההקצאה של assigned_var.
        """
        removed_vals = []
        # בדיקה אם קיים אילוץ כלשהו בין המשתנים
        if (assigned_var, var_j) in self.csp.constrains:
            for val_j in domain_j:
                # שימוש בפונקציה של המחלקה CSP לבדיקת הקונפליקט
                if self.csp.check_conflict(assigned_var, assigned_val, var_j, val_j):
                    removed_vals.append(val_j)
        return removed_vals

    def _check_domain_wipeout(self, var_j: int, domain_j: List[int], assigned_var: int) -> Tuple[
        Optional[int], Set[int]]:
        """
        בודקת אם הדומיין התרוקן.
        מחזירה: (המשתנה שנכשל או None, סט הקונפליקטים)
        """
        if not domain_j:
            # הדומיין התרוקן -> החזר את המשתנה שנכשל ואת הגורם (assigned_var)
            return var_j, {assigned_var}
        return None, set()

    def select_unassigned_variable(self):
        # ... (כפי שמומש קודם: MRV או סדר אינדקסים) ...
        pass

    def order_domain_values(self, var):
        # ... (כפי שמומש קודם: LCV או סדר רגיל) ...
        return self.domains[var]

    def forward_check(self, assigned_var, assigned_val, current_domains, assignment):
        # ... (הפונקציה הממומשת שחילקנו אותה קודם, ללא שינוי) ...
        # מחזירה: new_domains, failed_var, conflict_fc
        pass

    # -----------------------------------------------
    # ⭐️ פונקציה מרכזית לפתרון - Backtracking בסיסי
    # -----------------------------------------------

    def solve(self):
        return self._search()

    def _search(self):
        if len(self.assignment) == len(self.csp.variables):
            return self.assignment

        var_i = self.select_unassigned_variable()
        if var_i is None: return self.assignment

        # שמירת הדומיינים לפני ניסיון ההקצאה
        old_domains = copy.deepcopy(self.domains)

        for val_i in self.order_domain_values(var_i):
            self.assignment[var_i] = val_i

            # 1. Forward Check
            new_domains, failed_var_j, _ = self.forward_check(var_i, val_i, self.domains, self.assignment)

            if failed_var_j is None:
                # 2. אין כשל -> המשך רקורסיה
                self.domains = new_domains
                result = self._search()

                if result is not None:
                    return result

            # 3. טיפול בכשל / חזרה אחורה (Backtrack)

            # הסרת ההקצאה
            del self.assignment[var_i]
            self.backtracks += 1

            # שחזור הדומיינים מהעותק שנשמר
            self.domains = old_domains

        # 4. נכשלו כל הערכים -> חזרה אחורה ברקורסיה
        return None


class FC_CBJ_Solver(FC_Base_Solver):
    # אין צורך ב-__init__ חדש, יורש מ-FC_Base_Solver

    def find_backjump_var(self, conflict_set):
        """
        פונקציה ייחודית ל-CBJ.
        מוצאת את המשתנה המוקצה המאוחר ביותר ב-ConflictSet.
        """
        # ... (המימוש הקודם שהצעתי) ...
        pass

        # -----------------------------------------------

    # ⭐️ פונקציה מרכזית לפתרון - CBJ דורס את Base Search
    # -----------------------------------------------

    def _search(self):
        if len(self.assignment) == len(self.csp.variables):
            return self.assignment

        var_i = self.select_unassigned_variable()
        if var_i is None: return self.assignment

        old_domains = copy.deepcopy(self.domains)

        # ⚠️ אתחול Conflict Set לפני ניסיון הקצאה
        self.conflict_set[var_i] = set()

        for val_i in self.order_domain_values(var_i):
            self.assignment[var_i] = val_i

            # 1. Forward Check
            # conflict_fc מכיל את {var_i} אם היה כשל
            new_domains, failed_var_j, conflict_fc = self.forward_check(var_i, val_i, self.domains, self.assignment)

            if failed_var_j is None:
                # 2. אין כשל -> המשך רקורסיה
                self.domains = new_domains
                result = self._search()

                if result is not None:
                    return result

                # 3. חזרה רגילה לאחר כשל רקורסיבי (אם לא בוצע Backjump)
                # אם הפונקציה החזירה None, אנחנו לא יודעים אם זה כישלון או קפיצה.
                # במודל CBJ, אם הוחזר None, זה אומר שיש להשתמש ב-ConflictSet

                # **הערה חשובה:** ב-Python, אנחנו לא מעבירים "סיבת כשל" מפורשת.
                # אנחנו מסתמכים על כך שה-ConflictSet של var_i התעדכן
                # בקריאה הרקורסיבית שכשלה.

            else:
                # 4. ⭐️ טיפול בכשל FC (Domain Wipeout)
                self.backtracks += 1

                # 4.1. איחוד קבוצת הקונפליקטים (var_i + הגורמים לכישלון)
                self.conflict_set[var_i].update(conflict_fc)

                # 5. הסרת ההקצאה (תמיד) ושחזור דומיינים
            del self.assignment[var_i]
            self.domains = old_domains

        # 6. ⭐️ נכשלו כל הערכים -> Backjump!

        var_k = self.find_backjump_var(self.conflict_set[var_i])

        if var_k == -1:
            return None  # כשל גלובלי

        # 6.1. איחוד קבוצת הקונפליקטים של var_i לתוך var_k
        self.conflict_set[var_k].update(self.conflict_set[var_i])

        # 6.2. ביצוע Backjump - מחזירים None כדי שהרקורסיה הקודמת תמשיך
        return None


if __name__ == "__main__":
    mycsp = CSP(n_var=3,domain_size=3,density=0.34,tightness=0.33)