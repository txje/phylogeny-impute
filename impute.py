# --------------------------------------------------------------
# Jeremy Wang
# Dec 11, 2010
# 
# License: MIT (see LICENSE)
# 
# Impute sequence based on local phylogenetic structure
# (see https://github.com/txje/compatible-intervals)
#
# Note: intentionally does NOT impute a strain with itself
# --------------------------------------------------------------

import sys

def main(seqFile, cliqueFile, outputFile, target):

    print "  Reading cliques..."
    cliques, cliqueStrains = readCliques(cliqueFile)
    # our target strain must appear in the cliques
    if not target in cliqueStrains:
        print "%s not found in %s." % (target, cliqueFile)
        sys.exit(1)
    target_clique_index = cliqueStrains.index(target)

    print "  Reading genotype data..."
    genotypes, genotypeStrains = readGenotypes(seqFile)
    # map genotyped strains to clique indices
    genotyped_indices = [] # list of clique indices for which there exist genotypes
    genotype_map = {} # map clique indices to genotype indices
    for g in genotypeStrains:
        if g in cliqueStrains and g != target: # do NOT impute a strain with itself
            clique_index = cliqueStrains.index(g)
            genotyped_indices.append(clique_index)
            genotype_map[clique_index] = genotypeStrains.index(g)
    genotyped_indices.sort()
    print "    Using these genotyped strains to impute: " + str([cliqueStrains[g] for g in genotyped_indices])
    
    # -------------------------------------------------------
    # Imputation Method
    # -------------------------------------------------------

    print "  Initializing imputed arrays..."
    imputed = [None for i in genotypes]

    snps = len(genotypes)
    
    print "  Imputing %i variants..." % snps
    # for each SNP
    inv = 0 # carry the interval index along also
    for snp_index in xrange(len(genotypes)):
        while inv < len(cliques) and cliques[inv][1] < genotypes[snp_index][0]:
            inv += 1
        imputed[snp_index] = impute(genotypes, target, snp_index, genotyped_indices, genotype_map, target_clique_index, cliques, inv)
    
    print "  Writing imputed data..."
    output(outputFile, target, genotypes, imputed)
        
#-------------------------------------
# impute
# 
# - impute a given snp_index for the target strain
# - use overlapping intervals' cliques (HC),
#   expanding to use others as necessary (MC)
# - returns (allele, confidence {'h','m','l'})
# -------------------------------------
def impute(genotypes, target, snp_index, genotyped_indices, genotype_map, target_clique_index, cliques, inv):
    
    group = genotyped_indices + [target_clique_index] # default: everybody
    
    alleles = [genotypes[snp_index][1][genotype_map[g]] for g in genotyped_indices]
    
    # find source genotypes with private alleles or N/H
    private = []
    for i in xrange(len(alleles)):
        if alleles[i] in ['N', 'H'] or alleles.count(alleles[i]) == 1:
            private.append(i)
    # remove them in reverse order (to not screw up the list)
    for p in sorted(private, key = lambda a: a*-1):
        group.pop(p)
    
    counts = []
    
    # interval index offset should go 0, 1, -1, 2, -2, etc.
    # anything past 1 will be medium confidence
    inv_offset = 0
    
    while inv_offset > -10:
    
        # ran into the beginning or end
        if inv + inv_offset < 0 or inv + inv_offset >= len(cliques):
            if counts != []:
                return (counts[-1][0], 'm')
            else:
                return ('N', 'l')
        
        other = [g for g in cliques[inv+inv_offset][2] if target_clique_index in g][0]
        group = intersect(group, other, target_clique_index)
    
        alleles = [genotypes[snp_index][1][genotype_map[g]] for g in group if g != target_clique_index]
        
        oldcounts = counts # keep a backup in case we reduce too far
        counts = {}
        for a in alleles:
            counts[a] = counts.get(a, 0) + 1
        counts = sorted([(k,v) for k,v in counts.iteritems()], key = lambda a: a[1])
        
        if len(counts) == 1:
            if inv_offset >= 0 and inv_offset <= 1: # HC
                return (alleles[0], 'h')
            else: # MC (expanded)
                return (alleles[0], 'm')
        elif len(counts) == 0:
            if len(oldcounts) == 0: # LC, no source (ever)
                return ('N', 'l')
            else: # MC
                return (oldcounts[-1][0], 'm') # get most common before last reduction
    
        inv_offset *= -1
        if inv_offset >= 0:
            inv_offset += 1
    
    # so far unresolved, use majority
    return (counts[-1][0], 'm') # MC

# -------------------------------------
# intersect
# 
# - intersect group and other as sets
# - return the intersection with the target in it
# - reduces imputation sources for MC
#   based on neighboring intervals
# -------------------------------------

def intersect(group, other, target):
    both = set(group) & set(other)
    if target in both:
        return list(both)
    else:
        return []

# -------------------------------------
# output
# 
# - writes the given imputed genotypes
#   for a target strain to the output
#   file specified (outputFile) in the
#   form:
#
#   position,target_strn,confidence
#   pos,target_strn_allele,{0,1,2}
#   pos,target_strn_allele,{0,1,2}
#   ...
# 
# -------------------------------------
def output(outputFile, target, genotypes, imputed):
    fout = open(outputFile, 'w')
    fout.write("position," + target.replace('/', '_') + ",confidence")
    for i in xrange(len(imputed)):
        fout.write("\n" + str(genotypes[i][0]) + "," + str(imputed[i][0]) + "," + ("2" if imputed[i][1] == 'h' else ("1" if imputed[i][1] == 'm' else "0")))
    fout.close()

# -------------------------------------
# readCliques
#
# - given a filename, reads in the format:
#   start_pos,end_pos,strn0,strn1,strn2 ...
#   start position,end position,strn0 clique,strn1 clique,strn2 clique ...
#
# - returns a list of the form:
#   [
#       [start, end, [[clique0 strain indices], [clique1 strain indices], [clique2 strain indices] ...]]
#       [start, end, [[clique0 strain indices], [clique1 strain indices], [clique2 strain indices] ...]]
#       ...
#   ]
# -------------------------------------
def readCliques(fname):
    cliques = []

    data = [line.strip().split(',') for line in open(fname, 'r').read().strip().split('\n')]
    strains = data[0][2:]
    
    for line in data[1:]:
        d = [int(i) for i in line]

        groups = {}
        for i in xrange(2, len(d)):
            if not groups.has_key(d[i]):
                groups[d[i]] = []
            groups[d[i]].append(i - 2)
        groups = [v for k,v in groups.iteritems()] # in no particular order
        
        cliques.append([d[0], d[1], groups])
    
    return cliques, strains

# -------------------------------------
# readGenotypes
#
# - given a filename, reads in the compressed format:
#   pos,pos,pos,pos,strn0,strn1,strn2 ... '\n'
#   <pos - 4 bytes><allele - 1 byte><allele - 1 byte><allele - 1 byte> ...
#
# - the header line is ended with a newline
# - subsequent data lines have a fixed length and are not terminated by a newline
#
# - returns a list of the form:
#   [
#       [pos, [strn0 allele, strn1 allele, strn2 allele ...]]
#       [pos, [strn0 allele, strn1 allele, strn2 allele ...]]
#       ...
#   ]
# -------------------------------------
def readGenotypes(fname):
    data = open(fname, 'rb').read()
    
    header = data[:data.index('\n')].strip().split(',')
    strains = header[4:]

    genotype = []
    data = data[data.index('\n')+1:] # skip the header
    bytes_per = len(header)
    index = 0
    while index < len(data):
        line = data[index:index+bytes_per]
        genotype.append([ord(line[0])*(2**24) + ord(line[1])*(2**16) + ord(line[2])*(2**8) + ord(line[3]), line[4:]])
        index += bytes_per
    
    return genotype, strains

# -------------------------------------
# main: see usage
# -------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 5:
        print "usage: python impute.py <seqFile> <cliqueFile> <outputFile> <targetStrain>"
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
