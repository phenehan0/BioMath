
class Allele:
    def __init__(self, index=0, symbol="", dominance=1, frequency=0):
        self.index = index
        self.symbol = symbol
        self.dominance = dominance
        self.frequency = frequency

    def to_dict(self):
        return {
            "index": self.index,
            "name": self.symbol,
            "dominance": self.dominance, 
            "frequency": self.frequency
        }

class Genotype:
    def __init__(self, alleles, genotype_id=0, exp_frequency=0, obs_frequency=0):
        self.alleles = alleles
        self.genotype_id = genotype_id
        self.is_homozygote = True
        self.exp_frequency = exp_frequency
        self.obs_frequency = obs_frequency
    
    @property
    def is_homozygote(self):
        return self._is_homozygote

    @is_homozygote.setter
    def is_homozygote(self, alleles):
        self._is_homozygote = True
        init_idx = self.alleles[0].index
        for i in range(1, len(self.alleles)):
            if self.alleles[i].index != init_idx:
                self._is_homozygote = False
                break        
    
    @property
    def genotype_id(self):
        return self._genotype_id

    @genotype_id.setter
    def genotype_id(self, alleles):
        self._genotype_id = 0
        for i in self.alleles:
            self._genotype_id += i.index

class Gene:
    def __init__(self, alleles, popsize, genotypes=[]):
        self.alleles = alleles
        self.genotypes = genotypes
        self.popsize = popsize

        p, q = self.alleles[0], self.alleles[1]
        if (p.frequency + q.frequency)**2 != 1.0:
            raise ValueError("[ERROR] The sum of all genotype frequencies must equal 1.")
        
    @property
    def genotypes(self):
        return self._genotypes

    @genotypes.setter
    def genotypes(self, alleles):
        self._genotypes = []
        allele_frequencies = [a.frequency for a in self.alleles]
        # print(allele_frequencies)
        p = _polynomial_expansion(allele_frequencies)
        for (a1_idx, a2_idx) in p:
            freq = p[(a1_idx, a2_idx)]
            g = Genotype((self.alleles[a1_idx], self.alleles[a2_idx]), exp_frequency = freq)
            self._genotypes.append(g)

def _polynomial_expansion(values):
    #TODO this is only for binomials, should be generalized for polynomials of an arbitrary degree
    result = 0
    combos = []
    result = {}
    result_str = ""
    vals = enumerate(values)
    vals = list(vals)
    vals2 = vals.copy()
    for idx1, v1 in vals:
        for idx2, v2 in vals2:
            if (idx1, idx2) not in combos and idx1 == idx2:
                result[(idx1, idx2)] = v1*v2
                # print(str(v1)+"^2")
            elif (idx1, idx2) not in combos:
                result[(idx1, idx2)] = 2*v1*v2
                # print("2(" + str(v1) + ")("+str(v2)+")")
            combos.append((idx1, idx2))
            combos.append((idx2, idx1))
    return result

def chi_square_critical_value(dof, p=0.05):
    rows = []
    with open("chi_square_critical_values.txt", "r") as f:

        table = f.read()
        for i in table.split("\n"):
            rows.append(i.split())
        f.close()
    # print(rows)
    col_num = None
    for col_idx, col in enumerate(rows[0]):
        if float(col) == 1-p:
            col_num = col_idx
    if col_num and dof in range(len(rows)):
        return float(rows[dof][col_num])
    return

def chi_squared_test(expected, observed):
    result = 0
    if len(observed) != len(expected):
        raise ValueError("[ERROR] The observed and expected data sets are not of the same size.")
    n = len(expected)
    for e,o in zip(expected, observed):
        result += (1 / e) * ((e-o)**2)
    return result

def reject_null_hypothesis(expected, observed, dof, p=0.05):
    chi_squared = chi_squared_test(expected, observed)
    critical_value = chi_square_critical_value(dof, p=p)
    print("chi squared: "+str(chi_squared))
    print("critical value: "+str(critical_value))
    if chi_squared > critical_value:
        return True
    return False


if __name__ == "__main__":
    # example from Wiki page
    A = Allele(index = 0, symbol="A", dominance=1, frequency=0.954)
    a = Allele(index =1, symbol="a", dominance=0, frequency=0.046)
    g = Gene((A, a), 1612)
    
    # print(g.alleles)
    genotypes = g.genotypes
    popsize = g.popsize
    expected = []
    observed = [1469, 138, 5]
    for g in genotypes:
        expected.append(g.exp_frequency * popsize)

    print("expected: "+str(expected))
    print("observed: "+str(observed))
    print(reject_null_hypothesis(expected, observed, 1))