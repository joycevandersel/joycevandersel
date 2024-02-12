"""
https://rosalind.info/problems/iev/

Six nonnegative integers, each of which does not exceed 20,000.
The integers correspond to the number of couples in a population possessing
each genotype pairing for a given factor. In order, the six given integers
represent the number of couples having the following genotypes:

1. AA-AA
2. AA-Aa
3. AA-aa
4. Aa-Aa
5. Aa-aa
6. aa-aa

Return: The expected number of offspring displaying the dominant phenotype in
the next generation, under the assumption that every couple has exactly
two offspring.
"""

def calc_exp_dom_offspring(genotype):
    """ Calculating offspring with dominant phenotype

    :param genotype: string of 6 numbers
    :return: number of offspring with a dominant phenotype based on 2 offspring

    Here there are two groups based on mendel in category 1,2 and 3 all
    individuals will display the dominant feature. in category 4 75% and in 5
    50%. Category 6 will not result in a dominant phenotype and is excluded
    """
    genotype=genotype.split(' ')
    genotype=[int(x) for x in genotype]
    dom_off= (genotype[0]+genotype[1]+genotype[2])*2+genotype[3]*1.5+genotype[4]
    return dom_off

def main():
    """this is the main function of the script"""
    genotype= "16928 18810 18240 16049 19869 19834"
    exp_off=calc_exp_dom_offspring(genotype)
    print(exp_off)


if __name__ == "__main__":
    main()