import mysql.connector
from Bio import Entrez
import re
import time


def main():
    """ Calls all other functions."""
    dic_seq_hits = results_and_input()
    tax_organisms_erd = determining_taxonomy(dic_seq_hits)
    proteins_erd = determining_proteins(dic_seq_hits)
    too_database(dic_seq_hits, tax_organisms_erd, proteins_erd)


def results_and_input():
    """ Reads a file containing BLAST results and puts them in a dictionary."""
    file = open("BLAST_results_goed.txt", "r")
    dic_seq_hits = {}
    parameters = []
    hits = []
    header = ""

    for line in file:
        if line.startswith("~"):
            dic_seq_hits[header] = hits
            hits = []
            header = line.strip("~").strip('\n')
        elif line.startswith("$"):
            parameters.append(line.strip('\n').strip('$'))
            for _ in range(6):
                parameter = file.readline().strip('\n')
                parameter = parameter.replace("'", '"')
                parameters.append(parameter)
            if "hypothetical" not in parameters[1].lower():
                hits.append(parameters)
            parameters = []
    dic_seq_hits[header] = hits
    dic_seq_hits.pop('')

    return dic_seq_hits


def determining_taxonomy(dic_seq_hits):
    """ Detemines the taxonomy of accession codes from a dictionary containing BLAST results."""
    tax_organisms_erd = []
    organism_id = 1
    total = 1

    for key, hits in dic_seq_hits.items():
        for hit in hits:
            acc_code = hit[6]
            print(acc_code)
            Entrez.email = "larsmaasmill@gmail.com"  # Always tell NCBI who you are
            handle = Entrez.efetch(db='protein', id=acc_code, retmode='xml')
            record = Entrez.read(handle)
            taxonomy = record[0]["GBSeq_taxonomy"]
            taxonomy = taxonomy.split(';')
            tax_organisms_erd, organism_id = save_taxonomy(taxonomy, tax_organisms_erd, organism_id)
            found = False
            for tax in tax_organisms_erd:
                if tax[1] == taxonomy[-1]:
                    hit.append(tax[0])
                    found = True
                    break
            if found is False:
                hit.append('Null')
            print(str(total) + ' gedaan')
            total += 1
            time.sleep(0.4)

    return tax_organisms_erd


def save_taxonomy(taxonomy, tax_organisms_erd, organism_id):
    """ Saves the taxonomy obtained from an accession code."""
    organism_found = False
    organism = []
    for taxon in taxonomy:
        for taxons_erd in tax_organisms_erd:
            if taxon == taxons_erd[1]:
                organism_found = True
                break
        if organism_found is False:
            organism.append(organism_id)
            organism.append(taxon)
            reference = taxonomy.index(taxon) - 1
            if reference >= 0:
                for tax in tax_organisms_erd:
                    if tax[1] == taxonomy[reference]:
                        organism.append(tax[0])
                        break
            else:
                organism.append('Null')
            tax_organisms_erd.append(organism)
            organism_id += 1
        organism = []
        organism_found = False

    return tax_organisms_erd, organism_id


def determining_proteins(dic_seq_hits):
    """ Detemines the proteins from BLAST resultts saved in a dictionary."""
    proteins_erd = []
    protein_id = 1

    for hits in dic_seq_hits.values():
        for hit in hits:
            protein = re.split(r'^.*?\|.*?\|', hit[1])
            protein = protein[1]
            protein = re.split(r' \[.*', protein)
            protein = protein[0].strip(" ")
            protein = protein.replace("'", '"')
            proteins_erd, protein_id = save_proteins(proteins_erd, protein, protein_id)
            found = False
            for protein_entry in proteins_erd:
                if protein_entry[1] == protein:
                    hit.append(protein_entry[0])
                    found = True
                    break
            if found is False:
                hit.append('Null')

    return proteins_erd


def save_proteins(proteins_erd, protein, protein_id):
    """ Saves a protein name obtained from a BLAST result description."""
    protein_function = []
    found = False
    for protein_erd in proteins_erd:
        if protein.lower() == protein_erd[1].lower():
            found = True
            break
    if found is False:
        protein_function.append(protein_id)
        protein_function.append(protein)
        protein_id += 1
        proteins_erd.append(protein_function)

    return proteins_erd, protein_id


def too_database(dic_seq_hits, tax_organisms_erd, proteins_erd):
    """ Exports data obtained via BLAST results to a database."""
    connection = mysql.connector.connect(host="hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com",
                                         user="ophia@hannl-hlo-bioinformatica-mysqlsrv",
                                         db="ophia",
                                         password='594990',
                                         port="3306")
    cursor = connection.cursor()
    for protein in proteins_erd:
        cursor.execute("""insert into eiwitten (`Eiwit_id`, `Eiwit_naam`)
                        values ({}, '{}')""".format(protein[0], protein[1]))

    for taxon in tax_organisms_erd:
        cursor.execute("""insert into organismen (`Organisme_id`, `Taxonomische_Rang`, `Relatie_Tax`)
                                        value ({}, '{}', {})""".format(taxon[0], taxon[1], taxon[2]))

    for param, resultaten in dic_seq_hits.items():
        param = param.split('\t')
        cursor.execute("""insert into input (`Header`, `Sequentie`, `Sequentie_id`)
                           value ('{}', '{}', {})""".format(param[0], param[1], param[2]))
        for hit in resultaten:
            cursor.execute("""insert into resultaten (`Resultaat_id`, `Description`, `Bit_score`, `E_Value`, 
            `Perc_identity`, `Query_coverage`, `Accessiecode`, `Sequentie_id`, `Organisme_id`, `Eiwit_id`)
                            value ('{}', '{}', {}, '{}', '{}', '{}', '{}', {}, {}, {})""".format(
                            hit[0], hit[1], hit[2], hit[3], hit[4], hit[5], hit[6], param[2], hit[7], hit[8]))
    connection.commit()
    print('done')


main()
