import uuid
import os
from unittest import TestCase
from cartselection import helper
from tgmi.transcripts import Transcript

class TestHelper(TestCase):


    def test_read_excluded_transcripts(self):
        expected = {
            'NM_000247': 'no_mapping',
            'NM_000348': 'incorrect_cds_length',
            'NM_000451': 'multiple_mapping'
        }
        assert helper.read_excluded_transcripts('test/unit/data/excluded_test.txt') == expected


    def test_read_appris_file(self):
        result, result_principal = helper.read_appris_file('test/unit/data/appris_test.txt')

        assert result['84645'] == ['NM_032561']
        assert result_principal['NM_032561'] == 'PRINCIPAL:1'

        assert result['1759'] == ['NM_004408']
        assert result_principal['NM_004408'] == 'PRINCIPAL:3'

        assert result['79603'] == ['XM_011528293', 'NM_024552', 'XM_011528294']
        assert result_principal['XM_011528293'] == 'PRINCIPAL:1'
        assert result_principal['NM_024552'] == 'PRINCIPAL:1'
        assert result_principal['XM_011528294'] == 'PRINCIPAL:1'


    def test_read_refseqscan_output(self):
        expected = {
            'NM_001005484': '.',
            'NM_152486': 'c.1027T>C',
            'NM_015658': 'c.898A>G,c.1182T>C'
        }
        assert helper.read_refseqscan_output('test/unit/data/refseqscan_output_test.txt') == expected

    def test_read_gene_dict(self):
        expected = {
            'HGNC:100': '41',
            'HGNC:1000': '603',
            'HGNC:1101': '675'
        }
        assert helper.read_gene_dict('test/unit/data/genenames_test.txt') == expected


    def test_translate_gene_id(self):
        gene_dict = helper.read_gene_dict('test/unit/data/genenames_test.txt')

        prefix = str(uuid.uuid4())
        missing = open(prefix + '_missing.txt', 'w')
        log = open(prefix + '_log.txt', 'w')

        assert helper.translate_gene_id(gene_dict, 'HGNC:1000', missing, log) == '603'
        assert helper.translate_gene_id(gene_dict, 'HGNC:1234', missing, log) is None

        missing.close()
        log.close()

        with open(prefix + '_missing.txt') as f:
            s = f.read().strip()
        assert s == 'HGNC:1234\tno_NCBI_GeneID'

        with open(prefix + '_log.txt') as f:
            s = f.read().strip()
        assert s == '-> NCBI Gene ID: 603\n\n\nThe gene has no NCBI Gene ID\n\nGene added to Missing List'

        os.remove(prefix + '_missing.txt')
        os.remove(prefix + '_log.txt')


    def test_get_appris_principal_isoforms_for_gene(self):
        prefix = str(uuid.uuid4())
        missing = open(prefix + '_missing.txt', 'w')
        log = open(prefix + '_log.txt', 'w')

        appris, _ = helper.read_appris_file('test/unit/data/appris_test.txt')

        assert helper.get_appris_principal_isoforms_for_gene(appris, '79603', missing, log, 'xyz') == ['XM_011528293', 'NM_024552', 'XM_011528294']
        assert helper.get_appris_principal_isoforms_for_gene(appris, '12345', missing, log, 'xyz') is None

        missing.close()
        log.close()

        with open(prefix + '_missing.txt') as f:
            s = f.read().strip()
        assert s == 'xyz\tnot_in_APPRIS'

        with open(prefix + '_log.txt') as f:
            s = f.read().strip()
        assert s == 'The gene has the following APPRIS principal isoforms: ' + str(appris['79603']) + '\n\n\nThe gene has no APPRIS principal isoforms \n\nGene added to Missing List'

        os.remove(prefix + '_missing.txt')
        os.remove(prefix + '_log.txt')


    def test_get_nms(self):
        prefix = str(uuid.uuid4())
        missing = open(prefix + '_missing.txt', 'w')
        log = open(prefix + '_log.txt', 'w')

        appris_principal = ['XM_123456', 'NM_123456', 'XM_1234567', 'NM_1234567', 'NM_12345678', 'XM_12345678']
        assert helper.get_nms(appris_principal, missing, log, 'xyz') == ['NM_123456', 'NM_1234567', 'NM_12345678']

        appris_principal = ['XM_123456', 'XM_1234567', 'XM_12345678']
        assert helper.get_nms(appris_principal, missing, log, 'xyz') is None

        missing.close()
        log.close()

        with open(prefix + '_missing.txt') as f:
            s = f.read().strip()
        assert s == 'xyz\tonly_XM'

        with open(prefix + '_log.txt') as f:
            s = f.read().strip()
        assert s == 'APPRIS contains only XM principal isoforms for this gene; gene added to Missing List'

        os.remove(prefix + '_missing.txt')
        os.remove(prefix + '_log.txt')


    def test_same_cds(self):
        t1 = Transcript()
        record = {'strand': '+', 'exons': '10000-20000,30000-40000,50000-60000', 'coding_start': '12000', 'coding_end': '53000'}
        t1.read_from_database_record(record)

        t2 = Transcript()
        record = {'strand': '+', 'exons': '9000-20000,30000-40000,50000-68000', 'coding_start': '12000', 'coding_end': '53000'}
        t2.read_from_database_record(record)

        t3 = Transcript()
        record = {'strand': '+', 'exons': '3000-20000,30000-40000,50000-68000', 'coding_start': '12000', 'coding_end': '53000'}
        t3.read_from_database_record(record)

        t4 = Transcript()
        record = {'strand': '+', 'exons': '10000-20000,30000-40000,50000-60000', 'coding_start': '12500', 'coding_end': '53000'}
        t4.read_from_database_record(record)

        assert helper.same_cds(t1, t2)
        assert helper.same_cds(t1, t3)
        assert not helper.same_cds(t2, t4)


    def test_all_have_same_cds(self):
        t1 = Transcript()
        record = {'strand': '+', 'exons': '10000-20000,30000-40000,50000-60000', 'coding_start': '12000', 'coding_end': '53000'}
        t1.read_from_database_record(record)

        t2 = Transcript()
        record = {'strand': '+', 'exons': '9000-20000,30000-40000,50000-68000', 'coding_start': '12000', 'coding_end': '53000'}
        t2.read_from_database_record(record)

        t3 = Transcript()
        record = {'strand': '+', 'exons': '3000-20000,30000-40000,50000-68000', 'coding_start': '12000', 'coding_end': '53000'}
        t3.read_from_database_record(record)

        t4 = Transcript()
        record = {'strand': '+', 'exons': '10000-20000,30000-40000,50000-60000', 'coding_start': '12500', 'coding_end': '53000'}
        t4.read_from_database_record(record)

        assert not helper.all_have_same_cds([t1, t2, t3, t4])
        assert helper.all_have_same_cds([t1, t2, t3])


    def test_check_nms_in_refseq_db(self):
        pass


    def test_create_transcript_objects(self):
        pass