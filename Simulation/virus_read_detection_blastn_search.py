
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

E_VALUE_THRESH = 1e-10

"""
    note:
    Not using num_threads since is its known about blast that it may not use all the threads 
    allocated by the user. 
    Refernce: http://voorloopnul.com/blog/how-to-correctly-speed-up-blast-using-num_threads/
"""
def BLASTn_CHO_DB_search(temp_fasta_file=''):

    BLAST_xml_file= temp_fasta_file.split('_Threading.fa')[0]+'_BlastResult.xml'
    blastn_cline = NcbiblastnCommandline(query=temp_fasta_file, 
                                # db="/home/gridsan/raeuf/Raeuf_Notebook/Blast/db/RefSeq/CHO/CHO-K1", 
                                db="/home/gridsan/alexsli/Alexander_Notebook/Blast/db/CHO_K1/CHO_K1", 
                                         evalue=E_VALUE_THRESH,
                                         outfmt=5,
                                         num_descriptions= 10,
                                         num_alignments= 10,
                                         out=BLAST_xml_file, 
                                         num_threads=1)

    # stdout, stderr = blastn_cline()
    blastn_cline()

    result_handle= open(BLAST_xml_file, 'r')
    blast_records= NCBIXML.parse(result_handle)
    blast_record= next(blast_records)

    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                alignmentTitle= alignment.title
                if 'cricetulus griseus' in alignmentTitle.lower():
                    return True





def BLASTn_RVDB_DB_search(temp_fasta_file=''):
    
    BLAST_xml_file= temp_fasta_file.split('_Threading.fa')[0]+'_BlastResult.xml'
    blastn_cline = NcbiblastnCommandline(query=temp_fasta_file, 
                                # db="/home/gridsan/raeuf/Raeuf_Notebook/Blast/db/RVDB/RVDBv18.0", 
                                db="/home/gridsan/alexsli/Alexander_Notebook/Blast/db/RVDB/RVDBv18.0", 
                                         evalue=E_VALUE_THRESH,
                                         outfmt=5,
                                         num_descriptions= 5,
                                         num_alignments= 5,
                                         out=BLAST_xml_file,
                                         num_threads=1)

    # stdout, stderr = blastn_cline()
    blastn_cline()

    result_handle= open(BLAST_xml_file, 'r')
    blast_records= NCBIXML.parse(result_handle)
    blast_record= next(blast_records)

    result_list=[]
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                # if ('feline' in alignment.hit_def.lower() or 
                    # 'leukemia' in alignment.hit_def.lower()):
                # if ('minute' in alignment.hit_def.lower() or 
                    # 'mvm' in alignment.hit_def.lower()):
                # if ('porcine' in alignment.hit_def.lower() or 
                    # 'circovirus' in alignment.hit_def.lower()):
                    result=[alignment.hit_def, alignment.hit_id, 
                             alignment.accession, alignment.length,
                             blast_record.query, blast_record.query_length,
                             hsp.score, hsp.expect, 
                             hsp.identities, hsp.positives, 
                             hsp.gaps, hsp.align_length,
                             hsp.strand, hsp.query_start, 
                             hsp.query_end, hsp.query, 
                             hsp.match, hsp.sbjct, 
                             hsp.sbjct_start, hsp.sbjct_end]
                    result_list.append(result)
    return result_list






def BLASTn_RVDB_DB_search_result_parser(resultFile_txt, RVDB_result, alignment_num):

    print()
    print(f'\n****Alignment:{alignment_num}****')
    resultFile_txt.write(f'\n****Alignment:{alignment_num}****\n')

    print('sbjct_title:', RVDB_result[0])
    resultFile_txt.write(f'sbjct_title: {RVDB_result[0]}\n')

    print('sbjct_id:', RVDB_result[1])
    resultFile_txt.write(f'sbjct_id: {RVDB_result[1]}\n')

    print('sbjct_accession:',RVDB_result[2])
    resultFile_txt.write(f'sbjct_accession: {RVDB_result[2]}\n')

    print('sbjct_length:', RVDB_result[3])
    resultFile_txt.write(f'sbjct_length: {RVDB_result[3]}\n')

    print(f"query_id: {RVDB_result[4]}")
    resultFile_txt.write(f'query_id: {RVDB_result[4]}\n')

    print('query_length:',RVDB_result[5])
    resultFile_txt.write(f'query_length: {RVDB_result[5]}\n')

    print('score:',RVDB_result[6])
    resultFile_txt.write(f'score: {RVDB_result[6]}\n')

    print('evalue:',RVDB_result[7])
    resultFile_txt.write(f'evalue: {RVDB_result[7]}\n')

    print('identities:',RVDB_result[8])
    resultFile_txt.write(f'identities: {RVDB_result[8]}\n')

    print('positives:',RVDB_result[9])
    resultFile_txt.write(f'positives: {RVDB_result[9]}\n')

    print('gaps:',RVDB_result[10])
    resultFile_txt.write(f'gaps: {RVDB_result[10]}\n')

    print('align_length:',RVDB_result[11])
    resultFile_txt.write(f'align_length: {RVDB_result[11]}\n')

    print('strand:',RVDB_result[12])
    resultFile_txt.write(f'strand: {RVDB_result[12]}\n')

    print('query_start:',RVDB_result[13])
    resultFile_txt.write(f'query_start: {RVDB_result[13]}\n')

    print('query_end:',RVDB_result[14])
    resultFile_txt.write(f'query_end: {RVDB_result[14]}\n')

    print('query:',RVDB_result[15])
    resultFile_txt.write(f'query: {RVDB_result[15]}\n')

    print('match:',RVDB_result[16])
    resultFile_txt.write(f'match: {RVDB_result[16]}\n')

    print('sbjct:',RVDB_result[17])
    resultFile_txt.write(f'sbjct: {RVDB_result[17]}\n')

    print('sbjct_start:',RVDB_result[18])
    resultFile_txt.write(f'sbjct_start: {RVDB_result[18]}\n')

    print('sbjct_end:',RVDB_result[19])
    print('*****************\n')
    resultFile_txt.write(f'sbjct_end: {RVDB_result[19]}\n')
    resultFile_txt.write('*****************\n\n')
