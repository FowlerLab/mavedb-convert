import os
from unittest import TestCase

from ..main import parse_args


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(BASE_DIR + '/data/')
SRC = 'wt.fasta'


class TestArgParsing(TestCase):
    
    @staticmethod
    def mock_args(program, src, dst=None, wt_seq='AAA', offset=0,
                  one_based=False, score_column=None, hgvs_column=None,
                  input_type='scores', skip_header='0', skip_footer='0',
                  sheet_name=None, is_coding=True,
                  ):
        return {
            'enrich': True if program == 'enrich' else False,
            'enrich2': True if program == 'enrich2' else False,
            'empiric': True if program == 'empiric' else False,
            '<src>': os.path.join(DATA_DIR, src),
            '--dst': os.path.join(DATA_DIR, dst) if dst else dst,
            '--score_column': score_column,
            '--hgvs_column': hgvs_column,
            '--skip_header': skip_header,
            '--skip_footer': skip_footer,
            '--sheet_name': sheet_name,
            '--wtseq': wt_seq,
            '--offset': offset,
            '--input_type': input_type,
            '--one_based': one_based,
            '--is_coding': is_coding,
        }
        
    def test_io_error_invalid_file(self):
        with self.assertRaises(SystemExit):
            parse_args(
                self.mock_args(program='enrich2', src='missing_file.csv')
            )

    def test_error_input_type_not_scores_or_counts(self):
        with self.assertRaises(SystemExit):
            parse_args(
                self.mock_args(
                    program='enrich2', src=SRC, input_type='unknown')
            )

    def test_error_score_column_is_none_for_scores_input_empiric_or_enrich2(self):
        with self.assertRaises(SystemExit):
            parse_args(
                self.mock_args(
                    program='enrich', src=SRC, input_type='scores',
                    wt_seq='AAA', offset=0, score_column=None)
            )
        with self.assertRaises(SystemExit):
            parse_args(
                self.mock_args(
                    program='empiric', src=SRC, input_type='scores',
                    wt_seq='AAA', offset=0, score_column=None)
            )

    def test_ignore_score_column_for_counts_input_empiric_or_enrich2(self):
        parse_args(
            self.mock_args(
                program='enrich', src=SRC, input_type='counts',
                wt_seq='AAA', offset=0, score_column=None)
        )
        parse_args(
            self.mock_args(
                program='empiric', src=SRC, input_type='counts',
                wt_seq='AAA', offset=0, score_column=None)
        )

    def test_src_path_is_normalised(self):
        _, kwargs = parse_args(
            self.mock_args(program='enrich2', src=SRC)
        )
        self.assertEqual(
            kwargs['src'],
            os.path.normpath(os.path.expanduser(os.path.join(DATA_DIR, SRC)))
        )
        
    def test_dst_path_is_normalised(self):
        _, kwargs = parse_args(
            self.mock_args(program='enrich2', src=SRC, dst=BASE_DIR + '//data')
        )
        self.assertEqual(kwargs['dst'], DATA_DIR)

    def test_makes_dst_directory_tree(self):
        _, kwargs = parse_args(
            self.mock_args(program='enrich2',
                           src=SRC, dst=BASE_DIR + '/data/subdir/')
        )
        self.assertTrue(os.path.isdir(kwargs['dst']))
        os.rmdir(kwargs['dst'])
        
    def test_program_enrich_correctly_deduced(self):
        expected_p, _ = parse_args(
            self.mock_args(program='enrich', score_column='score',
                           src=SRC, wt_seq='AAA', offset=0)
        )
        self.assertEqual(expected_p, 'enrich')
        
    def test_program_enrich2_correctly_deduced(self):
        expected_p, _ = parse_args(
            self.mock_args(program='enrich2', src=SRC)
        )
        self.assertEqual(expected_p, 'enrich2')
        
    def test_program_empiric_correctly_deduced(self):
        expected_p, _ = parse_args(
            self.mock_args(program='empiric', src=SRC,
                           wt_seq='AAA', offset=0, score_column='score')
        )
        self.assertEqual(expected_p, 'empiric')
        
    def test_exit_non_int_offset(self):
        with self.assertRaises(SystemExit):
            parse_args(
                self.mock_args(
                    program='enrich2', src=SRC, wt_seq='AAA', offset='a')
            )
    
    def test_exit_is_coding_and_offset_is_not_mult_of_three(self):
        with self.assertRaises(SystemExit):
            parse_args(
                self.mock_args(
                    program='enrich2', src=SRC,
                    wt_seq='AAA', is_coding=True, offset=4)
            )

    def test_error_wtseq_not_multiple_of_three(self):
        for program in ('enrich', 'empiric', 'enrich2'):
            with self.assertRaises(SystemExit):
                parse_args(
                    self.mock_args(
                        program=program, src=SRC, wt_seq='AAAA')
                )
            
    def test_exit_none_wt_seq(self):
        for program in ('enrich', 'empiric'):
            with self.assertRaises(SystemExit):
                parse_args(
                    self.mock_args(
                        program=program, src=SRC, 
                        wt_seq=None, offset=0)
                )
            
    def test_exit_invalid_dna_chars(self):
        for program in ('enrich', 'empiric'):
            with self.assertRaises(SystemExit):
                parse_args(
                    self.mock_args(
                        program=program, src=SRC, 
                        wt_seq='atxkjsf', offset=0)
                )
                
    def test_can_load_sequence_from_fasta(self):
        for program in ('enrich', 'empiric'):
            fasta = os.path.join(DATA_DIR, 'lower.fa')
            _, kwargs = parse_args(
                self.mock_args(
                    program=program, src=SRC, wt_seq=fasta, offset=0,
                    score_column='score'
                )
            )
            expected = (
                "ACAGTTGGATATAGTAGTTTGTACGAGTTGCTTGTGGCTT"
                "CGCCAGCGCATACCAGCATAGTAAAGGCAACGGCCTCTGA"
                "GAGGCTACGATCGTGCCTTGTGGCAAGTCTTCGCTCGCAC"
                "GCCCTTCCTACCGTGCTATGAGAGGAAATCTCGGGCGTA"
            )
            self.assertEqual(kwargs['wtseq'], expected)

    def test_parses_skip_footer_and_header(self):
        _, args = parse_args(
            self.mock_args(program='enrich2', src=SRC,
                           skip_header='1', skip_footer='0')
        )
        self.assertEqual(args['skip_header'], 1)
        self.assertEqual(args['skip_footer'], 0)
        
    def test_parses_one_based_as_false(self):
        _, args = parse_args(
            self.mock_args(program='enrich2', src=SRC, one_based=False)
        )
        self.assertEqual(args['one_based'], False)

    def test_parses_sheet_name(self):
        _, args = parse_args(
            self.mock_args(program='enrich2', src=SRC, sheet_name='TEF')
        )
        self.assertEqual(args['sheet_name'], "TEF")
