def get_fasta_from_label(labelfile,fastafile):
    fasta = open(fastafile, 'w')
    with open(labelfile, 'r') as label:
        chain0 = ''
        chain=[]; sequence=[]
        for aaterm in label:
            aaterm = aaterm.split()
            chain1 = '>' + aaterm[0]
            residue = aaterm[1]
            if chain1 != chain0 and chain0:
                chain.append(chain0)
                fasta.write(''.join(chain) + '\n')
                outfasta = ''.join(sequence) + '\n'
                fasta.write(outfasta)
                chain = [];sequence=[]
                sequence.append(residue)
            else:
                sequence.append(residue)
            chain0 = chain1
        chain.append(chain0)
        fasta.write(''.join(chain) + '\n')
        outfasta = ''.join(sequence) + '\n'
        fasta.write(outfasta)
    fasta.close()
