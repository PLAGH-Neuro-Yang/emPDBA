dna_codes = {
    'DA': 'A', 'DT': 'T', 'DG': 'G', 'DC': 'C',
}

def massDNA(filename):
    import struct
    in_file = open(filename)
    pdb_format = '6s5s1s4s1s3s1s1s4s1s3s8s8s8s6s6s6s4s2s3s'
    total = 0
    aminomass = []
    for line in in_file:
        if line[0:4] == "ATOM":
            col = struct.unpack(pdb_format, line.encode())
            amino = col[5].strip().decode("utf-8")
            numer = col[8].strip().decode("utf-8")
            if amino in dna_codes:
                if [amino,numer] not in aminomass:
                    aminomass.append([amino,numer])
    total = (len(aminomass)*303.7)+79
    return total