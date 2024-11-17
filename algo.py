import numpy as np

def is_complementary(base1, base2):
    """Check if two bases are complementary."""
    pairs = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
    return base2 == pairs.get(base1, '')

def nussinov(sequence):
    """Nussinov algorithm to predict RNA secondary structure."""
    n = len(sequence)
    dp = np.zeros((n, n), dtype=int)

    # Fill the DP table
    for l in range(1, n):  # Length of the subsequence
        for i in range(n - l):
            j = i + l
            pair_score = 1 if is_complementary(sequence[i], sequence[j]) else 0
            dp[i, j] = max(
                dp[i+1, j],               # Skip i
                dp[i, j-1],               # Skip j
                dp[i+1, j-1] + pair_score # Pair (i, j)
            )

    # Traceback to build the structure
    structure = ['.'] * n
    def traceback(i, j):
        if i >= j:
            return
        if dp[i, j] == dp[i+1, j]:
            traceback(i+1, j)
        elif dp[i, j] == dp[i, j-1]:
            traceback(i, j-1)
        elif dp[i, j] == dp[i+1, j-1] + (1 if is_complementary(sequence[i], sequence[j]) else 0):
            structure[i] = '('
            structure[j] = ')'
            traceback(i+1, j-1)

    traceback(0, n-1)
    return ''.join(structure)

def run_nussinov(sequences):
    """Run Nussinov algorithm on a list of sequences."""
    results = []
    for seq in sequences:
        structure = nussinov(seq["sequence"])
        results.append({
            "id": seq["id"],
            "algorithm": "Nussinov",
            "structure": structure,
            "sequence": seq["sequence"]
        })
    return results

import RNA

def zuker(sequence):
    """Zuker algorithm to predict RNA secondary structure."""
    fc = RNA.fold_compound(sequence)
    mfe_structure, mfe = fc.mfe()
    return mfe_structure, mfe

def run_zuker(sequences):
    """Run Zuker algorithm on a list of sequences."""
    results = []
    for seq in sequences:
        structure, mfe = zuker(seq["sequence"])
        results.append({
            "id": seq["id"],
            "algorithm": "Zuker",
            "structure": structure,
            "mfe": mfe,
            "sequence": seq["sequence"]
        })
    return results

def extract_base_pairs(dot_bracket):
    """Extract base pair indices from dot-bracket notation."""
    stack = []
    pairs = set()
    for i, char in enumerate(dot_bracket):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                j = stack.pop()
                pairs.add((j, i))
    return pairs

def compare_structures(struct1, struct2):
    """Compare two RNA structures using base pair overlap."""
    pairs1 = extract_base_pairs(struct1)
    pairs2 = extract_base_pairs(struct2)
    common = pairs1.intersection(pairs2)
    return len(common) / max(len(pairs1), len(pairs2), 1)  # Fraction of shared pairs