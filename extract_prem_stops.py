#!/usr/bin/env python

from __future__ import print_function
import argparse
import gzip
import subprocess
import pysam
from vcf2raw_sfs import get_minor_freq
from degen_to_bed import start_stop_ok, cds_coord_to_codon_coords


def reverse_strand(codon_coords):

    """
    identifies if codon coordinates go backwards
    :param codon_coords: list
    :return: bool
    """

    if codon_coords[0] > codon_coords[1] > codon_coords[2]:
        return True

    elif codon_coords[0] < codon_coords[1] < codon_coords[2]:
        return False

    else:
        raise ValueError


def multi_nucleotide_poly(snp_position, snp_set, codon_pos):

    """
    identifies if adjacent bases in the codon are polymorphic
    :param snp_position: int
    :param snp_set: set
    :param codon_pos: int
    :return: bool
    """

    if codon_pos == 0 and snp_position + 1 in snp_set:
        return True

    elif codon_pos == 1:
        if snp_position + 1 in snp_set:
            return True
        if snp_position - 1 in snp_set:
            return True

    elif codon_pos == 2 and snp_position - 1 in snp_set:
        return True

    else:
        return False


def snp_makes_stop(pysam_snp, codon_pos, codon, codon_coords):

    """
    takes a biallelic snp and checks to see if it
    makes a premature stop codon
    :param pysam_snp: pysam.VariantRecord
    :param codon_pos: int
    :param codon: str
    :param codon_coords: list
    :return: bool
    """

    ref = pysam_snp.ref.upper()
    alt = pysam_snp.alts[0].upper()

    # reverse comp if neccessary
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    if reverse_strand(codon_coords):
        ref = comp[ref]
        alt = comp[alt]

    stop = prem_stop(codon, codon_pos, bases=(ref, alt))

    return stop


def prem_stop(codon, pos, bases=('A', 'T', 'G', 'C')):

    """
    goes through a list of bases and sees if any make a stop codon
    :param codon: str
    :param pos: int
    :param bases: set
    :return: bool
    """

    stops = {'TAG', 'TGA', 'TAA'}

    for base in bases:

        # construct possible codons
        if pos == 0:
            alt_codon = base + codon[1:]
        elif pos == 1:
            alt_codon = codon[0] + base + codon[2]
        elif pos == 2:
            alt_codon = codon[0:2] + base
        else:
            raise IndexError

        if alt_codon in stops:
            return True

    return False


def cds_snp_coords(vcf):

    """
    greps CDS snp coords from vcf and makes a dict of sets
    :param vcf: str
    :return: dict
    """

    grep = 'zgrep ANNO=CDS {} | cut -f 1,2'.format(vcf)

    cds_snps = subprocess.Popen(grep, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
    cds_snps_dict = {}

    for snp in cds_snps:
        chromo, pos = snp.split()[0], int(snp.split()[1])
        if chromo not in cds_snps_dict.keys():
            cds_snps_dict[chromo] = set()

        cds_snps_dict[chromo] |= {pos}

    return cds_snps_dict


def process_transcript(sequence, coords, trans_name, nonsense_data, call_seq, snps, chromo, vcf_records, n):

    """
    calculates length, call sites and nonsense snps for a transcript
    :param sequence: str
    :param coords: list
    :param trans_name: str
    :param nonsense_data: dict
    :param call_seq: str
    :param snps: set
    :param chromo: str
    :param vcf_records: pysam.VariantFile
    :param n: int
    :return: None
    """

    if sequence != '' and len(sequence) % 3 == 0 and start_stop_ok(sequence):
        for i in range(0, len(sequence) - 3, 3):

            length = len(sequence)
            codon = sequence[i:i + 3]

            # skip codons with Ns
            if 'N' in codon:
                continue

            base_pos = coords[i:i + 3]

            for pos in range(0, 3):
                site_pos = base_pos[pos]
                if trans_name not in nonsense_data.keys():
                    nonsense_data[trans_name] = {'length': length, 'call': 0, 'freqs': []}

                # is a premature stop possible?
                if prem_stop(codon, pos):

                    # skip if site is not callable
                    call = call_seq[site_pos - 1]

                    if call.upper() != 'K':
                        continue

                    nonsense_data[trans_name]['call'] += 1

                    # if snp at site
                    if site_pos in snps[chromo]:

                        # check not an MNP
                        if multi_nucleotide_poly(site_pos, snps[chromo], pos):
                            continue

                        # does it introduce a stop?
                        snp_data = list(vcf_records.fetch(chromo, site_pos - 1, site_pos))
                        snp_record = snp_data[0]

                        # skip snps that don't make stop
                        if not snp_makes_stop(snp_record, pos, codon, base_pos):
                            continue

                        frequency = get_minor_freq(snp_record, 'snp', n)
                        nonsense_data[trans_name]['freqs'].append(frequency)

                    # move on if no snp at position
                    else:
                        continue

                # if site cannot give a premature stop
                else:
                    continue


def main():
    # arguments
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-cds_fa', help='Fasta file with CDS sequences in', required=True)
    parser.add_argument('-chr', help='chromosome', required=True)
    parser.add_argument('-vcf', help='SNP vcf path', required=True)
    parser.add_argument('-call_fa', help='Callable sites fasta file', required=True)
    parser.add_argument('-n', help='Sample size', required=True, type=int)
    parser.add_argument('-out', help='output file', required=True)
    args = parser.parse_args()

    # variables
    fa = args.cds_fa
    nonsense_data = {}
    snps = cds_snp_coords(args.vcf)
    vcf_records = pysam.VariantFile(args.vcf)
    call_fa = pysam.FastaFile(args.call_fa)
    call_seq = call_fa.fetch(args.chr)
    n = args.n

    # loop through fasta
    sequence, chromo, coords, skip, trans_name = '', '', [], False, ''

    for line in gzip.open(fa):

        # skip sequence line following skipped header
        if skip is True:
            skip = False
            continue

        # if a header line
        if line.startswith('>'):

            # skips records with both and forward and reverse strand in coords
            if 'complement' in line and 'join' in line:
                skip = True
                continue

            # process prev sequence
            process_transcript(sequence, coords, trans_name, nonsense_data, call_seq, snps, chromo, vcf_records, n)

            # reset holders and store details of next sequence
            sequence = ''
            header_info = line.split(';')[1]
            chromo = header_info.split(':')[0].replace('loc=', '').strip()
            trans_name = line.split(' ')[0].strip('>')
            if '(' in header_info:
                method = header_info.split('(')[0].split(':')[1]
                coords = cds_coord_to_codon_coords(header_info.split('(')[1].rstrip(')'), method)
            else:
                method = 'join'
                coords = cds_coord_to_codon_coords(header_info.split(':')[1], method)

        # if a sequence line add to seq string
        else:
            sequence += line.rstrip()

    # process final sequence
    process_transcript(sequence, coords, trans_name, nonsense_data, call_seq, snps, chromo, vcf_records, n)

    # output sites
    with open(args.out, 'w') as out:
        for trans in sorted(nonsense_data.keys()):

            gene = '-'.join(trans.split('-')[:-1])
            transcript = trans.split('-')[-1]
            trans_len = nonsense_data[trans]['length']
            n_call = nonsense_data[trans]['call']
            allele_freqs = nonsense_data[trans]['freqs']
            if len(allele_freqs) == 0:
                allele_freqs = 'nan'
            else:
                allele_freqs = ','.join([str(x) for x in allele_freqs])

            print(gene, transcript, args.chr, trans_len, n_call, allele_freqs, sep='\t', file=out)


if __name__ == '__main__':
    main()
